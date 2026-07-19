/*
 * ucmuon_transport_pumas.c
 * UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
 *
 * PUMAS-based muon transport for UCMuon.
 *
 * Modes:
 *   0 = forward  — reads surface muon file, transports forward, writes 18-col output
 *   1 = backward — backward MC flux integration, writes per-event flux file
 *
 * Stdin protocol (one token per line):
 *   mode          0=forward | 1=backward
 *   mdf_path      absolute path to external/pumas-master/examples/data/materials.xml
 *   dump_path     path to .pumas physics dump cache (auto-created on first run)
 *   mat_name      StandardRock | Water     (Water covers Ice/Seawater via rho)
 *   rho_gcm3      density [g/cm³] (0 = material default)
 *   depth_m       overburden depth [m]
 *   energy_loss   0=CSDA | 1=MIXED | 2=STRAGGLED
 *   scattering    0=disabled | 1=mixed
 *
 *   [if mode==0 forward]:
 *     infile         path to muons_surface.dat
 *     outfile        path for muons_underground.dat (18-col)
 *     transport_all  0=hit_flag==1 only | 1=all muons
 *     depth_axis     0=depth is X | 1=depth is Y | 2=depth is Z
 *     seed           RNG seed (0 = time-based)
 *
 *   [if mode==1 backward]:
 *     outfile        path for pumas_bwd_events.dat (per-event flux file)
 *     E_min_GeV      minimum detector kinetic energy [GeV]
 *     E_max_GeV      maximum detector kinetic energy [GeV]
 *     cos_theta_min  minimum cos(zenith)  (1.0 = vertical only)
 *     cos_theta_max  maximum cos(zenith)  (0.0 = horizontal)
 *     n_events       number of MC events
 *     spectrum_id    0=GCCLY (Guan et al.) | 1=Gaisser
 *     seed           RNG seed (0 = time-based)
 */

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "pumas.h"

/* ── Constants ─────────────────────────────────────────────────────────── */
#define MMUON_GEV   0.10566          /* muon rest mass [GeV] */
#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#define PROGRESS_EVERY  500

/* ── Global PUMAS handles ───────────────────────────────────────────────── */
static struct pumas_physics *g_physics = NULL;
static struct pumas_context *g_context = NULL;

/* ── Medium: uniform slab with density override ─────────────────────────── */
static double g_rho_kgm3 = 2650.0;  /* [kg/m³] */

static double locals_uniform(struct pumas_medium *medium,
    struct pumas_state *state, struct pumas_locals *locals)
{
    locals->density = g_rho_kgm3;
    return 0.;  /* uniform → no stepping hint */
}

static struct pumas_medium g_medium = { 0, &locals_uniform };

/* Infinite medium (geometry handled by distance limit per muon) */
static enum pumas_step medium_infinite(struct pumas_context *ctx,
    struct pumas_state *state, struct pumas_medium **med_ptr,
    double *step_ptr)
{
    if (med_ptr)  *med_ptr  = &g_medium;
    if (step_ptr) *step_ptr = 0.;
    return PUMAS_STEP_CHECK;
}

/* ── xoshiro256++ RNG ───────────────────────────────────────────────────── */
static uint64_t rng_s[4];

static inline uint64_t rng_rotl(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

static uint64_t rng_next(void)
{
    const uint64_t r = rng_rotl(rng_s[0] + rng_s[3], 23) + rng_s[0];
    const uint64_t t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0];
    rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2];
    rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t;
    rng_s[3]  = rng_rotl(rng_s[3], 45);
    return r;
}

static double ucmuon_random(struct pumas_context *ctx)
{
    return (rng_next() >> 11) * (1.0 / (double)(1LL << 53));
}

static void rng_seed(uint64_t seed)
{
    for (int i = 0; i < 4; i++) {
        seed += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng_s[i] = z ^ (z >> 31);
    }
}

/* ── Atmospheric flux models ────────────────────────────────────────────── */
static double charge_fraction(double charge)
{
    const double cr = 1.2766;  /* CMS charge ratio */
    if      (charge < 0.) return 1. / (1. + cr);
    else if (charge > 0.) return cr / (1. + cr);
    else                  return 1.;
}

static double cos_theta_star(double cos_theta)
{
    /* Volkova's parameterization */
    const double p[] = { 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285 };
    const double cs2 =
        (cos_theta * cos_theta + p[0]*p[0] + p[1]*pow(cos_theta, p[2]) +
         p[3]*pow(cos_theta, p[4])) / (1. + p[0]*p[0] + p[1] + p[3]);
    return cs2 > 0. ? sqrt(cs2) : 0.;
}

/* GCCLY: Guan et al. 2015 (https://arxiv.org/abs/1509.06176) */
static double flux_gccly(double cos_theta, double ke_GeV, double charge)
{
    const double Emu = ke_GeV + MMUON_GEV;
    const double cs  = cos_theta_star(cos_theta);
    const double ec  = 1.1 * Emu * cs;
    const double rpi = 1. + ec / 115.;
    const double rK  = 1. + ec / 850.;
    const double f_gaisser = 1.4e3 * pow(Emu, -2.7) * (1./rpi + 0.054/rK)
                             * charge_fraction(charge);
    return pow(1. + 3.64 / (Emu * pow(cs, 1.29)), -2.7) * f_gaisser;
}

/* Gaisser (PDG ch.30) */
static double flux_gaisser(double cos_theta, double ke_GeV, double charge)
{
    const double Emu = ke_GeV + MMUON_GEV;
    const double ec  = 1.1 * Emu * cos_theta;
    const double rpi = 1. + ec / 115.;
    const double rK  = 1. + ec / 850.;
    return 1.4e3 * pow(Emu, -2.7) * (1./rpi + 0.054/rK)
           * charge_fraction(charge);
}

typedef double (*flux_fn_t)(double, double, double);

/* ── Physics init ───────────────────────────────────────────────────────── */
static void physics_init(const char *mdf_path, const char *dump_path,
    const char *mat_name, int energy_loss_id, int scattering_id)
{
    /* Try dump first (fast, ~100 ms vs. ~10 s for MDF parse) */
    FILE *f = fopen(dump_path, "rb");
    if (f) {
        fprintf(stdout, "  Loading physics from dump: %s\n", dump_path);
        fflush(stdout);
        enum pumas_return rc = pumas_physics_load(&g_physics, f);
        fclose(f);
        if (rc != PUMAS_RETURN_SUCCESS) g_physics = NULL;
    }

    if (!g_physics) {
        fprintf(stdout, "  Building physics from MDF (first run, ~10 s): %s\n",
                mdf_path);
        fflush(stdout);
        if (pumas_physics_create(&g_physics, PUMAS_PARTICLE_MUON,
                mdf_path, NULL, NULL) != PUMAS_RETURN_SUCCESS) {
            fprintf(stderr, "ERROR: pumas_physics_create failed for '%s'\n",
                    mdf_path);
            exit(1);
        }
        f = fopen(dump_path, "wb");
        if (f) {
            pumas_physics_dump(g_physics, f);
            fclose(f);
            fprintf(stdout, "  Physics dump saved: %s\n", dump_path);
            fflush(stdout);
        }
    }

    if (pumas_physics_material_index(g_physics, mat_name, &g_medium.material)
            != PUMAS_RETURN_SUCCESS) {
        fprintf(stderr, "ERROR: material '%s' not in MDF\n", mat_name);
        exit(1);
    }

    pumas_context_create(&g_context, g_physics, 0);
    g_context->medium = &medium_infinite;
    g_context->random = &ucmuon_random;

    /* Energy loss mode: 0=CSDA, 1=MIXED, 2=STRAGGLED */
    if      (energy_loss_id == 0) g_context->mode.energy_loss = PUMAS_MODE_CSDA;
    else if (energy_loss_id == 1) g_context->mode.energy_loss = PUMAS_MODE_MIXED;
    else                          g_context->mode.energy_loss = PUMAS_MODE_STRAGGLED;

    /* Scattering: 0=disabled, 1=mixed */
    g_context->mode.scattering =
        (scattering_id == 0) ? PUMAS_MODE_DISABLED : PUMAS_MODE_MIXED;

    fprintf(stdout, "  PUMAS ready: mat=%s  energy_loss=%d  scattering=%d\n",
            mat_name, energy_loss_id, scattering_id);
    fflush(stdout);
}

/* ── String helper ──────────────────────────────────────────────────────── */
static int read_line(char *buf, int sz, FILE *src)
{
    buf[0] = '\0';
    while (fgets(buf, sz, src)) {
        /* strip trailing newline/space */
        int n = (int)strlen(buf);
        while (n > 0 && (buf[n-1] == '\n' || buf[n-1] == '\r' ||
                         buf[n-1] == ' '  || buf[n-1] == '\t'))
            buf[--n] = '\0';
        if (n > 0) return 1;
    }
    return 0;
}

/* ── FORWARD MODE ───────────────────────────────────────────────────────── */
static void run_forward(const char *infile, const char *outfile,
    int transport_all, int depth_axis, double depth_m)
{
    FILE *fin = fopen(infile, "r");
    if (!fin) {
        fprintf(stderr, "ERROR: cannot open '%s': %s\n", infile, strerror(errno));
        exit(1);
    }

    /* Count lines to know total */
    int total = 0;
    {
        char tmp[512];
        while (fgets(tmp, sizeof(tmp), fin)) {
            if (tmp[0] == '#' || tmp[0] == '\n') continue;
            total++;
        }
        rewind(fin);
    }

    FILE *fout = fopen(outfile, "w");
    if (!fout) {
        fprintf(stderr, "ERROR: cannot open '%s': %s\n", outfile, strerror(errno));
        fclose(fin);
        exit(1);
    }

    fprintf(fout,
        "# UCMuon PUMAS forward transport\n"
        "# depth=%.2f m  rho=%.3f g/cm3  energy_loss=%d  scattering=%d\n"
        "# EventID xs ys zs Es theta_s phi_s charge alive "
        "x y z E cx cy cz theta phi\n",
        depth_m, g_rho_kgm3 / 1e3,
        (int)g_context->mode.energy_loss,
        (g_context->mode.scattering == PUMAS_MODE_DISABLED) ? 0 : 1);

    /* Configure forward mode */
    g_context->mode.direction = PUMAS_MODE_FORWARD;
    g_context->event |= PUMAS_EVENT_LIMIT_DISTANCE;

    int transported = 0, survived = 0;
    char line[512];
    while (fgets(line, sizeof(line), fin)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        int   evid, charge_i;
        double xs, ys, zs, p_srf, px, py, pz, theta_s, phi_s, E_srf;
        int   hit_flag = 1;

        int nc = sscanf(line,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
            &evid, &xs, &ys, &zs, &p_srf, &px, &py, &pz,
            &theta_s, &phi_s, &E_srf, &charge_i, &hit_flag);
        if (nc < 12) continue;
        if (!transport_all && hit_flag != 1) continue;

        double KE_srf = E_srf - MMUON_GEV;
        if (KE_srf <= 0.) continue;

        /* Direction cosines from momentum or theta/phi */
        double cx, cy, cz;
        if (p_srf > 0.) {
            cx = px / p_srf;
            cy = py / p_srf;
            cz = pz / p_srf;
        } else {
            cx =  sin(theta_s) * cos(phi_s);
            cy =  sin(theta_s) * sin(phi_s);
            cz = -cos(theta_s);
        }

        /* Slant path based on depth axis */
        double cos_depth;
        if      (depth_axis == 0) cos_depth = fabs(cx);
        else if (depth_axis == 1) cos_depth = fabs(cy);
        else                      cos_depth = fabs(cz);
        if (cos_depth < 0.01) cos_depth = 0.01;
        double slant_m = depth_m / cos_depth;

        /* PUMAS state: position in metres, KE in GeV */
        struct pumas_state state = {
            .charge    = (double)charge_i,
            .energy    = KE_srf,
            .weight    = 1.,
            .position  = { xs * 1e-2, ys * 1e-2, zs * 1e-2 },
            .direction = { cx, cy, cz }
        };

        g_context->limit.distance = slant_m;
        enum pumas_event event;
        pumas_context_transport(g_context, &state, &event, NULL);

        /* alive: muon reached the target depth */
        int alive = (event & PUMAS_EVENT_LIMIT_DISTANCE) ? 1 : 0;

        double Ef, xf, yf, zf, cxf, cyf, czf, theta_f, phi_f;
        if (alive) {
            Ef    = state.energy + MMUON_GEV;
            xf    = state.position[0] * 100.;
            yf    = state.position[1] * 100.;
            zf    = state.position[2] * 100.;
            cxf   = state.direction[0];
            cyf   = state.direction[1];
            czf   = state.direction[2];
            theta_f = acos(fmax(-1., fmin(1., -czf)));
            phi_f   = atan2(cyf, cxf);
        } else {
            Ef  = 0.;
            xf  = xs; yf = ys; zf = zs;
            cxf = 0.; cyf = 0.; czf = -1.;
            theta_f = 0.; phi_f = 0.;
        }

        fprintf(fout,
            "%10d %13.4f %13.4f %13.4f %13.6f %13.6f %13.6f %4d %2d "
            "%13.4f %13.4f %13.4f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f\n",
            evid, xs, ys, zs, E_srf, theta_s, phi_s, charge_i, alive,
            xf, yf, zf, Ef, cxf, cyf, czf, theta_f, phi_f);

        transported++;
        survived += alive;
        if (transported % PROGRESS_EVERY == 0) {
            fprintf(stdout, "  Transported: %d  Survived: %d  Total: %d\n",
                    transported, survived, total);
            fflush(stdout);
        }
    }

    fclose(fin);
    fclose(fout);

    fprintf(stdout, "  Transported: %d  Survived: %d  Total: %d\n",
            transported, survived, transported);
    fprintf(stdout, "  Muons transported: %d\n", transported);
    fprintf(stdout, "  Survived:          %d\n", survived);
    fprintf(stdout, "  Survival rate:     %.4f %%\n",
            transported > 0 ? 100. * survived / transported : 0.);
    fprintf(stdout, "  Output:            %s\n", outfile);
    fflush(stdout);
}

/* ── BACKWARD MODE ──────────────────────────────────────────────────────── */
static void run_backward(const char *outfile,
    double E_min, double E_max,
    double cos_theta_min, double cos_theta_max,
    int n_events, int spectrum_id, double depth_m)
{
    flux_fn_t flux_fn = (spectrum_id == 0) ? flux_gccly : flux_gaisser;

    /* Configure backward mode: stop when slant path distance is reached */
    g_context->mode.direction = PUMAS_MODE_BACKWARD;
    g_context->event |= PUMAS_EVENT_LIMIT_DISTANCE;

    const double rk         = log(E_max / E_min);   /* log-energy range */
    const double w_angle    = 2. * M_PI * (cos_theta_min - cos_theta_max);
    /* w_angle = ∫ dΩ over [theta_min, theta_max] = 2π*(cos_min - cos_max) */

    FILE *fout = fopen(outfile, "w");
    if (!fout) {
        fprintf(stderr, "ERROR: cannot open '%s': %s\n", outfile, strerror(errno));
        exit(1);
    }

    fprintf(fout,
        "# PUMAS backward MC  depth=%.2f m  rho=%.3f g/cm3\n"
        "# E_range=[%.3f, %.3f] GeV (KE)  cos_theta=[%.3f, %.3f]\n"
        "# spectrum=%s  N=%d\n"
        "# ev  E_det_GeV  cos_theta  charge  E_surf_GeV  "
        "flux_contribution[m-2 s-1 GeV-1 sr-1]\n",
        depth_m, g_rho_kgm3 / 1e3,
        E_min, E_max, cos_theta_max, cos_theta_min,
        spectrum_id == 0 ? "GCCLY" : "Gaisser",
        n_events);

    double sum_w = 0., sum_w2 = 0.;
    int ok = 0;

    for (int i = 0; i < n_events; i++) {
        /* Sample detector kinetic energy: log-uniform [E_min, E_max] */
        double E_det = E_min * exp(rk * ucmuon_random(g_context));
        double wf_E  = E_det * rk;  /* 1/PDF for log-uniform sampling */

        /* Sample cos(theta): uniform [cos_theta_max, cos_theta_min]
         * (cos_theta_max < cos_theta_min for zenith range [0, theta_max]) */
        double cos_t = cos_theta_max
                       + (cos_theta_min - cos_theta_max) * ucmuon_random(g_context);
        double sin_t = sqrt(fmax(0., 1. - cos_t * cos_t));

        /* Charge: randomise equally */
        double cf    = (ucmuon_random(g_context) > 0.5) ? 1. : -1.;
        double wf    = wf_E * 2.;  /* energy × charge sampling weights */

        /* Slant path for this zenith */
        double slant_m = depth_m / cos_t;

        /* PUMAS backward: start at "detector" (position irrelevant for flux)
         * direction = muon direction at detection (downward into detector)
         * In backward mode PUMAS moves opposite → upward through rock
         */
        struct pumas_state state = {
            .charge    = cf,
            .energy    = E_det,
            .weight    = wf,
            .direction = { sin_t, 0., -cos_t },  /* pointing downward */
            .position  = { 0., 0., 0. }
        };

        g_context->limit.distance = slant_m;
        enum pumas_event event;
        pumas_context_transport(g_context, &state, &event, NULL);

        /* Muon reached surface if distance limit hit (not stopped/decayed) */
        if (!(event & PUMAS_EVENT_LIMIT_DISTANCE)) continue;

        /* Surface cos(theta): backward mode flips direction, so surface
         * direction = -state.direction (muon was going downward, now going up)
         * cos_theta_surface = -state.direction[2] */
        double cos_t_surf = -state.direction[2];
        if (cos_t_surf < 0.) cos_t_surf = 0.;

        /* Kinetic energy at surface */
        double KE_surf = state.energy;

        /* Flux at surface × PUMAS weight = flux contribution at detector */
        double f_surf  = flux_fn(cos_t_surf, KE_surf, state.charge);
        double contrib = state.weight * f_surf;

        /* Per-sr per-GeV (dΦ/dE_det at fixed direction): */
        double flux_val = contrib;

        sum_w  += flux_val;
        sum_w2 += flux_val * flux_val;
        ok++;

        fprintf(fout,
            "%8d %13.6f %13.6f %6.0f %13.6f %15.6e\n",
            i, E_det, cos_t, cf, KE_surf, flux_val);

        if ((i + 1) % PROGRESS_EVERY == 0) {
            fprintf(stdout, "  Transported: %d  Survived: %d  Total: %d\n",
                    i + 1, ok, n_events);
            fflush(stdout);
        }
    }

    fclose(fout);

    /* Total integrated flux rate [m-2 s-1]:
     * Rate = w_angle × (1/N) × Σ flux_val
     * (w_angle accounts for solid angle integration)
     */
    double rate     = w_angle * sum_w  / n_events;
    double rate_err = w_angle * sqrt(fmax(0., (sum_w2 / n_events)
                                         - (sum_w / n_events) * (sum_w / n_events))
                                     / n_events);

    fprintf(stdout, "  Transported: %d  Survived: %d  Total: %d\n",
            n_events, ok, n_events);
    fprintf(stdout, "  Muons transported: %d\n", n_events);
    fprintf(stdout, "  Survived:          %d\n", ok);
    fprintf(stdout, "  Rate:    %.4e +/- %.4e m-2 s-1\n", rate, rate_err);
    fprintf(stdout, "  Output:  %s\n", outfile);
    fflush(stdout);
}

/* ── main ───────────────────────────────────────────────────────────────── */
int main(void)
{
    char buf[1024];

#define NEXT(def) (read_line(buf, sizeof(buf), stdin) ? buf : (def))

    int mode = atoi(NEXT("1"));

    char mdf_path[1024], dump_path[1024], mat_name[64];
    strncpy(mdf_path,  NEXT("materials.xml"), sizeof(mdf_path)  - 1);
    strncpy(dump_path, NEXT("pumas.pumas"),   sizeof(dump_path) - 1);
    strncpy(mat_name,  NEXT("StandardRock"),  sizeof(mat_name)  - 1);

    double rho_gcm3    = atof(NEXT("2.65"));
    double depth_m     = atof(NEXT("100.0"));
    int    energy_loss = atoi(NEXT("0"));
    int    scattering  = atoi(NEXT("0"));

    g_rho_kgm3 = (rho_gcm3 > 0.) ? rho_gcm3 * 1e3 : 2650.;

    fprintf(stdout, "UCMuon PUMAS driver  mode=%s  mat=%s  rho=%.3f g/cm3"
            "  depth=%.1f m\n",
            mode == 0 ? "forward" : "backward",
            mat_name, g_rho_kgm3 / 1e3, depth_m);
    fflush(stdout);

    physics_init(mdf_path, dump_path, mat_name, energy_loss, scattering);

    if (mode == 0) {
        /* ── forward ─────────────────────────────────────────────────── */
        char infile[1024], outfile[1024];
        strncpy(infile,  NEXT("muons_surface.dat"),    sizeof(infile)  - 1);
        strncpy(outfile, NEXT("muons_underground.dat"), sizeof(outfile) - 1);
        int transport_all = atoi(NEXT("0"));
        int depth_axis    = atoi(NEXT("2"));

        /* Seed the RNG: without this the xoshiro state is all-zero and
         * random() returns 0 forever — NaN energies in straggled mode and
         * infinite rejection loops with scattering enabled. */
        uint64_t seed = (uint64_t)atol(NEXT("0"));
        if (seed == 0) seed = (uint64_t)time(NULL);
        rng_seed(seed);

        fprintf(stdout, "  Forward: infile=%s  outfile=%s  axis=%d  seed=%llu\n",
                infile, outfile, depth_axis, (unsigned long long)seed);
        fflush(stdout);

        run_forward(infile, outfile, transport_all, depth_axis, depth_m);

    } else {
        /* ── backward ────────────────────────────────────────────────── */
        char outfile[1024];
        strncpy(outfile, NEXT("pumas_bwd_events.dat"), sizeof(outfile) - 1);

        double E_min         = atof(NEXT("1.0"));
        double E_max         = atof(NEXT("1000.0"));
        double cos_theta_min = atof(NEXT("1.0"));   /* = cos(0°) = vertical */
        double cos_theta_max = atof(NEXT("0.0"));   /* = cos(90°) = horizontal */
        int    n_events      = atoi(NEXT("50000"));
        int    spectrum_id   = atoi(NEXT("0"));
        uint64_t seed        = (uint64_t)atol(NEXT("0"));

        if (seed == 0) seed = (uint64_t)time(NULL);
        rng_seed(seed);

        fprintf(stdout, "  Backward: E=[%.2f, %.2f] GeV  cos_theta=[%.3f, %.3f]"
                "  N=%d  spectrum=%s  seed=%llu\n",
                E_min, E_max, cos_theta_max, cos_theta_min, n_events,
                spectrum_id == 0 ? "GCCLY" : "Gaisser",
                (unsigned long long)seed);
        fflush(stdout);

        run_backward(outfile, E_min, E_max,
                     cos_theta_min, cos_theta_max,
                     n_events, spectrum_id, depth_m);
    }

    pumas_context_destroy(&g_context);
    pumas_physics_destroy(&g_physics);

    return 0;
}
