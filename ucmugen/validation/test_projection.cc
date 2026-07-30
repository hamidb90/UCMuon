// Analytic validation of the projection-correct generator.
//
// Phase 1 validates the physics (spectra, angular models, charge ratio) against
// the Fortran bit-for-bit. That says nothing about the surface projection,
// which the Fortran does not implement. So the projection is validated instead
// against closed-form results and against geometric identities that must hold
// regardless of the spectrum.
//
// Build:
//   c++ -std=c++17 -O2 -I../include test_projection.cc -o test_projection

#include "UCMuGen.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using namespace ucmugen;

static int g_failures = 0;

static void check(const char* what, double got, double expect, double rtol,
                  const char* note = "") {
  const double rel = std::fabs(got - expect) / std::fabs(expect);
  const bool ok = rel <= rtol;
  if (!ok) ++g_failures;
  std::printf("  [%s] %-46s got=%.6g  expect=%.6g  rel=%.2e (tol %.0e) %s\n",
              ok ? "ok  " : "FAIL", what, got, expect, rel, rtol, note);
}

static Generator make(Spectrum s, std::shared_ptr<Surface> surf,
                      double theta_max_deg) {
  Generator g;
  g.setSpectrum(s)
      .setSurface(std::move(surf))
      .setEnergyRange(1.0, 1000.0)
      .setThetaRange(0.0, theta_max_deg * kPi / 180.0);
  return g;
}

int main() {
  // ---------------------------------------------------------------------
  // 1. Closed form for an isotropic spectrum on a horizontal plane.
  //
  // CosmoALEPH has no zenith dependence, so the rate factorises exactly:
  //     R = A * I_p * int_{theta<theta_max} cos(theta) dOmega
  //       = A * I_p * pi * sin^2(theta_max)
  // This is the tightest available test of the whole normalisation chain:
  // surface area, solid angle, the momentum integral, and the projection.
  // ---------------------------------------------------------------------
  std::printf("1. isotropic spectrum on a horizontal plane (closed form)\n");
  for (double tmax : {30.0, 60.0, 70.0, 89.0}) {
    auto plane = std::make_shared<Plane>(100.0, 100.0, Vec3{0, 0, 0});
    Generator g = make(Spectrum::CosmoALEPH, plane, tmax);
    const double I_p = g.cdf().integrated_flux();
    const double s = std::sin(tmax * kPi / 180.0);
    const double expect = g.surfaceArea() * I_p * kPi * s * s;
    char label[64];
    std::snprintf(label, sizeof(label), "theta_max = %.0f deg", tmax);
    check(label, g.rate(400000), expect, 5e-3);
  }

  // ---------------------------------------------------------------------
  // 2. Shape independence of a flat horizontal surface.
  //
  // The rate through a flat horizontal surface depends on its area and the
  // sky, not on its shape. Disk and Plane use completely different position
  // samplers, so agreement here exercises both.
  // ---------------------------------------------------------------------
  std::printf("\n2. flat surface: rate depends on area, not shape\n");
  {
    const double R = 100.0;
    const double side = std::sqrt(kPi * R * R);      // equal-area square
    auto disk = std::make_shared<Disk>(R, Vec3{0, 0, 0});
    auto plane = std::make_shared<Plane>(side / 2.0, side / 2.0, Vec3{0, 0, 0});
    Generator gd = make(Spectrum::Guan, disk, 70.0);
    Generator gp = make(Spectrum::Guan, plane, 70.0);
    check("disk area == plane area", gd.surfaceArea(), gp.surfaceArea(), 1e-12);
    check("disk rate == plane rate", gd.rate(400000), gp.rate(400000), 5e-3);
  }

  // ---------------------------------------------------------------------
  // 3. Tilted plane, closed form for an isotropic sky.
  //
  // With an isotropic sky spanning the downward hemisphere, a plane whose
  // normal is tilted by alpha from vertical collects
  //     R(alpha) / R(0) = (1 + cos(alpha)) / 2
  //
  // Note this is NOT the naive projected-area law cos(alpha). That law holds
  // for a parallel beam, or for a field isotropic over the full sphere. Here
  // the field exists only over the downward hemisphere, so tilting the plane
  // brings previously-excluded near-horizon directions into view at the same
  // time as it loses others, and the rate falls off far more slowly. At 45 deg
  // the difference is 21%: cos gives 0.707, the truth is 0.854. Verified
  // against direct numerical integration of
  //     int max(0, sin a sin t cos f + cos a cos t) sin t dt df
  // over t in [0, pi/2].
  //
  // This whole weighting is what the Fortran omits for tilted surfaces.
  // ---------------------------------------------------------------------
  std::printf("\n3. tilted plane vs isotropic sky (the weighting Fortran omits)\n");
  {
    auto flat = std::make_shared<Plane>(100.0, 100.0, Vec3{0, 0, 0});
    Generator g0 = make(Spectrum::CosmoALEPH, flat, 89.999);
    const double r0 = g0.rate(1000000);
    for (double alpha : {15.0, 30.0, 45.0, 60.0}) {
      const double a = alpha * kPi / 180.0;
      auto tilted = std::make_shared<Plane>(
          100.0, 100.0, Vec3{0, 0, 0}, Vec3{std::sin(a), 0.0, std::cos(a)});
      Generator g = make(Spectrum::CosmoALEPH, tilted, 89.999);
      char label[64];
      std::snprintf(label, sizeof(label),
                    "R(%.0f deg)/R(0) = (1+cos a)/2", alpha);
      char note[64];
      std::snprintf(note, sizeof(note), "[naive cos a would be %.4f]", std::cos(a));
      check(label, g.rate(1000000) / r0, 0.5 * (1.0 + std::cos(a)), 5e-3, note);
    }
  }

  // ---------------------------------------------------------------------
  // 4. Cylinder caps.
  //
  // EcoMug's cylinder is lateral-surface-only, so its J' carries sin^2(theta)
  // and a vertical muon can never be generated. Including the caps is not a
  // refinement: without them the vertical flux onto a squat cylinder is
  // missing entirely. For a short, wide cylinder the cap contribution
  // dominates.
  // ---------------------------------------------------------------------
  std::printf("\n4. cylinder caps (EcoMug's cylinder omits these)\n");
  {
    struct Case { double r, h; const char* shape; };
    for (Case c : {Case{100.0, 400.0, "tall  (h/r=4)"},
                   Case{200.0,  50.0, "squat (h/r=0.25)"}}) {
      auto with = std::make_shared<Cylinder>(c.r, c.h, Vec3{0, 0, 0}, true);
      auto without = std::make_shared<Cylinder>(c.r, c.h, Vec3{0, 0, 0}, false);
      Generator gw = make(Spectrum::Guan, with, 89.0);
      Generator gn = make(Spectrum::Guan, without, 89.0);
      const double rw = gw.rate(400000), rn = gn.rate(400000);
      std::printf("  [info] %-46s with caps=%.4g Hz  lateral only=%.4g Hz  "
                  "missing %.1f%%\n", c.shape, rw, rn, 100.0 * (1.0 - rn / rw));
    }
  }

  // ---------------------------------------------------------------------
  // 5. Generated sample reproduces the integrated rate.
  //
  // rate() integrates the density; generate() samples it by rejection. They
  // are separate code paths sharing only the weight function, so the mean
  // projection over generated muons must match the ratio the integrator
  // implies. This catches an envelope that silently truncates the tail.
  // ---------------------------------------------------------------------
  std::printf("\n5. sampler and integrator agree on <cos(theta)>\n");
  {
    auto plane = std::make_shared<Plane>(100.0, 100.0, Vec3{0, 0, 0});
    Generator g = make(Spectrum::CosmoALEPH, plane, 70.0);
    // Isotropic sky weighted by cos(theta): <cos> = int c^2 dc / int c dc
    // over [cos(theta_max), 1] = (2/3)(1-c0^3)/(1-c0^2).
    const double c0 = std::cos(70.0 * kPi / 180.0);
    const double expect = (2.0 / 3.0) * (1.0 - c0 * c0 * c0) / (1.0 - c0 * c0);
    double sum = 0.0;
    const int n = 400000;
    for (int i = 0; i < n; ++i) sum += -g.generate().direction.z;
    check("<cos theta> of generated muons", sum / n, expect, 5e-3);
    std::printf("  [info] acceptance = %.1f%%, envelope = %.3f\n",
                100.0 * g.acceptance(), g.envelope());
  }

  std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "ALL PASSED",
              g_failures, g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
