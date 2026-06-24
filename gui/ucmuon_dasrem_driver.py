# ucmuon_dasrem_driver.py — guaranteed-hit surface muon generator
# UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
# MIT License 2026
#
# Reverse-sampling strategy (100 % detector hits) following the approach of:
#   Yao et al., J. Appl. Phys. 138, 144901 (2025)
#
# Every generated muon is guaranteed to hit the detector (hit_flag = 1).
# Pure Python / NumPy — no Fortran or subprocess needed.
#
# Algorithm (horizontal XY source at z = z_src, detector top at z = z_top < z_src):
#   1. Sample (θ, φ) from user-selected angular distribution.
#   2. Sample momentum p from selected spectrum model.
#   3. Sample hit point (x_d, y_d) uniformly on detector face.
#   4. Back-project to source surface via straight-line geometry (D = z_src − z_top):
#        x_src = x_d − sinθ cosφ · D / cosθ
#        y_src = y_d − sinθ sinφ · D / cosθ
#        z_src = source_z_cm  (user-configured; default 0)
#
# Output: standard 14-column UCMuon format (hit_flag = 1, det_mask = 1 for all rows).

import math, sys, time
import numpy as np
from pathlib import Path

M_MU_GEV    = 0.10566        # muon rest mass [GeV]
MU_PLUS_FRAC = 1.27 / 2.27   # P(μ+) from μ+/μ− = 1.27 charge ratio


# ── Spectrum samplers ─────────────────────────────────────────────────────────

def _cosmo_sample(p_min, p_max, N, rng):
    """CosmoALEPH: dN/dp ∝ p^-3.1952 — exact inverse-CDF."""
    ap1 = -3.1952 + 1.0   # = -2.1952
    lo  = p_min**ap1;  hi = p_max**ap1
    return (lo + rng.random(N) * (hi - lo)) ** (1.0 / ap1)


def _power_sample(p_min, p_max, N, rng, alpha=-2.7):
    """Power-law dN/dp ∝ p^alpha — inverse-CDF (log-uniform when alpha = -1)."""
    ap1 = alpha + 1.0
    if abs(ap1) < 1e-10:
        return np.exp(rng.uniform(math.log(p_min), math.log(p_max), N))
    lo = p_min**ap1;  hi = p_max**ap1
    return (lo + rng.random(N) * (hi - lo)) ** (1.0 / ap1)


def _rejection_sample(flux_fn, theta_mean_deg, p_min, p_max, N, rng):
    """
    Rejection sampling for theta-dependent spectra (Guan, Frosin, Bugaev, Reyna).
    Uses the batch mean zenith angle for the envelope fit, which is consistent
    with how the Fortran generator samples p and theta independently.

    Envelope: g(p) = A · p^-2.7  (conservative upper bound for all supported models).
    """
    alpha = -2.7
    ap1   = alpha + 1.0   # = -1.7
    lo    = p_min**ap1;  hi = p_max**ap1

    # Compute envelope scale from a p-grid at the batch mean theta
    p_grid = np.geomspace(p_min, p_max, 300)
    T_grid = np.sqrt(p_grid**2 + M_MU_GEV**2) - M_MU_GEV
    f_grid = flux_fn(T_grid, theta_mean_deg)
    env_grid = p_grid**alpha
    valid = env_grid > 1e-300
    scale = np.nanmax(f_grid[valid] / env_grid[valid]) * 1.15   # 15 % headroom

    result  = np.empty(N)
    n_done  = 0
    while n_done < N:
        n_need = (N - n_done) * 4
        p_try  = (lo + rng.random(n_need) * (hi - lo)) ** (1.0 / ap1)
        T_try  = np.sqrt(p_try**2 + M_MU_GEV**2) - M_MU_GEV
        f_try  = flux_fn(T_try, theta_mean_deg)
        env    = scale * p_try**alpha
        accept = rng.random(n_need) < f_try / np.where(env > 1e-300, env, 1e-300)
        p_acc  = p_try[accept]
        take   = min(len(p_acc), N - n_done)
        result[n_done:n_done + take] = p_acc[:take]
        n_done += take
    return result


# ── Angular distribution samplers ─────────────────────────────────────────────

def _sample_theta(angular_mode, theta_max_deg, N, rng):
    """Sample zenith angles from the selected angular distribution."""
    if angular_mode == 1:        # vertical pencil beam (θ = 0)
        return np.zeros(N)
    cos_max = math.cos(math.radians(theta_max_deg))
    u = rng.random(N)
    if angular_mode == 2:        # cos²θ  (recommended for muons)
        cos_th = (u * (1.0 - cos_max**3) + cos_max**3) ** (1.0 / 3.0)
    elif angular_mode == 5:      # cos³θ  (Reyna-Bugaev)
        cos_th = (u * (1.0 - cos_max**4) + cos_max**4) ** 0.25
    else:                        # uniform in cosθ (isotropic hemisphere)
        cos_th = cos_max + u * (1.0 - cos_max)
    return np.arccos(np.clip(cos_th, 0.0, 1.0))


# ── Detector hit-point samplers ───────────────────────────────────────────────

def _sample_cylinder(det, N, rng):
    """Uniform disk sample on the top face of a cylinder."""
    cx = (det["ax"] + det["bx"]) / 2.0
    cy = (det["ay"] + det["by"]) / 2.0
    r  = det["r"]
    r_s   = r * np.sqrt(rng.random(N))
    phi_s = rng.random(N) * 2.0 * math.pi
    return cx + r_s * np.cos(phi_s), cy + r_s * np.sin(phi_s)


def _sample_box(det, N, rng):
    """Uniform sample on the top face of an axis-aligned box."""
    x = rng.uniform(det["xmin"], det["xmax"], N)
    y = rng.uniform(det["ymin"], det["ymax"], N)
    return x, y


def _det_depth_cm(det, source_z_cm=0.0):
    """
    Distance from the source surface (z = source_z_cm) to the TOP face of the
    detector [cm].  Returns a positive value when the source is above the detector.
    """
    if det["shape"] == 1:   # Cylinder: top face = max(az, bz)
        z_top = max(det.get("az", 0.0), det.get("bz", 0.0))
    else:                   # Box AABB: top face = max(zmin, zmax)
        z_top = max(det.get("zmin", 0.0), det.get("zmax", 0.0))
    return source_z_cm - z_top


# ── Main generator ────────────────────────────────────────────────────────────

def generate_dasrem(N, spectrum_mode, e_min, e_max, angular_mode, theta_max,
                    detectors, output_file,
                    source_z_cm=0.0, seed=None, progress_fn=None):
    """
    Generate N surface muons with guaranteed detector hits (DAS-REM).

    Parameters
    ----------
    N             : int   — number of muons to write
    spectrum_mode : int   — 1=CosmoALEPH, 2=Power-law E^-3.7, 4=Guan, 5=Frosin,
                           6=Bugaev, 7=Reyna-Bugaev (3/8 unsupported → fall back to 1)
    e_min, e_max  : float — kinetic energy range [GeV]
    angular_mode  : int   — 1=vertical (θ=0), 2=cos²θ, 3=uniform cone, 5=cos³θ
                           (4 unsupported → falls back to uniform cone)
    theta_max     : float — maximum zenith angle [degrees]
    detectors     : list  — detector dicts from GUI session state (uses detectors[0])
    output_file   : str   — path to output .dat file (14-column UCMuon format)
    source_z_cm   : float — z coordinate of the source surface [cm] (default 0)
    seed          : int or None — RNG seed for reproducibility
    progress_fn   : callable(n_done, N, line) or None — called after each batch

    Returns
    -------
    (n_written, elapsed_s)
    """
    _gui_dir = Path(__file__).parent
    if str(_gui_dir) not in sys.path:
        sys.path.insert(0, str(_gui_dir))
    from fast_flux_estimator import _reyna_bugaev, _bugaev, _guan_2015, _frosin_2025

    rng   = np.random.default_rng(seed)
    p_min = math.sqrt(max(e_min, M_MU_GEV + 1e-6)**2 - M_MU_GEV**2)
    p_max = math.sqrt(max(e_max, M_MU_GEV + 1e-6)**2 - M_MU_GEV**2)

    det      = detectors[0]
    depth_cm = _det_depth_cm(det, source_z_cm)
    if depth_cm <= 0.0:
        if det["shape"] == 1:
            z_top = max(det.get("az", 0.0), det.get("bz", 0.0))
        else:
            z_top = max(det.get("zmin", 0.0), det.get("zmax", 0.0))
        raise ValueError(
            f"DAS-REM requires the source surface to be above the detector top face. "
            f"Source z = {source_z_cm:.1f} cm, detector top = {z_top:.1f} cm. "
            f"Either lower the detector or raise the source surface z."
        )

    _flux_fns = {
        4: _guan_2015, 5: _frosin_2025, 6: _bugaev, 7: _reyna_bugaev,
    }

    BATCH  = min(50_000, N)
    n_done = 0
    t0     = time.time()

    with open(output_file, "w") as fh:
        fh.write("# DAS-REM generator — UCLouvain Muography Group\n")
        fh.write("# Yao et al., J. Appl. Phys. 138, 144901 (2025)\n")
        fh.write(f"# Source surface z = {source_z_cm:.1f} cm  "
                 f"Detector depth (top face): {depth_cm:.1f} cm  "
                 f"shape={'Cylinder' if det['shape']==1 else 'Box'}\n")
        fh.write(f"# Spectrum mode {spectrum_mode}  "
                 f"E=[{e_min},{e_max}] GeV  angular_mode={angular_mode}  "
                 f"theta_max={theta_max} deg\n")
        fh.write("# EventID x[cm] y[cm] z[cm] p[GeV/c] px py pz "
                 "theta[rad] phi[rad] E[GeV] charge hit_flag det_mask\n")

        while n_done < N:
            n_batch = min(BATCH, N - n_done)

            # ── Angles ───────────────────────────────────────────────────────
            theta  = _sample_theta(angular_mode, theta_max, n_batch, rng)
            phi    = rng.random(n_batch) * 2.0 * math.pi
            cos_th = np.cos(theta)
            sin_th = np.sin(theta)
            cx     = sin_th * np.cos(phi)
            cy     = sin_th * np.sin(phi)
            cz     = -cos_th     # downward (negative z)

            # ── Momenta ───────────────────────────────────────────────────────
            if spectrum_mode == 1:
                p = _cosmo_sample(p_min, p_max, n_batch, rng)
            elif spectrum_mode == 2:
                p = _power_sample(p_min, p_max, n_batch, rng, alpha=-2.7)
            elif spectrum_mode in _flux_fns:
                th_mean = math.degrees(float(np.mean(theta)))
                p = _rejection_sample(_flux_fns[spectrum_mode], th_mean,
                                      p_min, p_max, n_batch, rng)
            else:
                p = _cosmo_sample(p_min, p_max, n_batch, rng)

            # ── Detector hit point → back-project to surface ─────────────────
            if det["shape"] == 1:
                x_d, y_d = _sample_cylinder(det, n_batch, rng)
            else:
                x_d, y_d = _sample_box(det, n_batch, rng)

            t_path = depth_cm / cos_th   # path length from source surface to detector [cm]
            x_src  = x_d - cx * t_path
            y_src  = y_d - cy * t_path
            z_src  = np.full(n_batch, source_z_cm)

            # ── Momentum components & energy ──────────────────────────────────
            px_v = p * cx
            py_v = p * cy
            pz_v = p * cz
            E    = np.sqrt(p**2 + M_MU_GEV**2)
            chg  = np.where(rng.random(n_batch) < MU_PLUS_FRAC, 1, -1)

            # ── Write rows ────────────────────────────────────────────────────
            for i in range(n_batch):
                fh.write(
                    f"{n_done + i + 1:10d}"
                    f" {x_src[i]:14.4f} {y_src[i]:14.4f} {z_src[i]:12.4f}"
                    f" {p[i]:12.6f} {px_v[i]:12.6f} {py_v[i]:12.6f} {pz_v[i]:12.6f}"
                    f" {theta[i]:10.6f} {phi[i]:10.6f}"
                    f" {E[i]:12.6f} {chg[i]:3d}"
                    f" 1 1\n"
                )

            n_done  += n_batch
            elapsed  = time.time() - t0
            rate     = n_done / max(elapsed, 0.001)
            # No thousands-separator in Saved/tried: parse_progress uses \d+ which
            # stops at commas, so "50,000" would be read as "50".
            line     = (f"  Saved {n_done} ... tried {n_done}"
                        f"  ({100.0 * n_done / N:.1f}%)"
                        f"  {rate:,.0f} muons/s")
            if progress_fn:
                progress_fn(n_done, N, line)

    return n_done, time.time() - t0
