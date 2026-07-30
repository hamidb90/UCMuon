// Validation of detector-directed sampling.
//
// Directed sampling is a variance-reduction trick, and the failure mode of
// every such trick is that it stays plausible while quietly changing the
// normalisation. So nothing here checks that the generator "runs": each test
// pins the absolute rate to something known independently of the sampler.
//
//   1. far field    the rate into a distant detector is fixed by its projected
//                   area alone, which gives a closed form for four shapes and
//                   exercises all three intersection routines
//   2. equivalence  directed sampling and the plain acceptance cut must agree
//                   on the rate; they share no sampling code path
//   3. containment  a detector that intercepts every ray must reproduce the
//                   undetected rate, including the degenerate full-sphere cone
//   4. speed        what the whole exercise was for
//
// Build:
//   c++ -std=c++17 -O2 -I../include test_detector.cc -o test_detector

#include "UCMuGen.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>

using namespace ucmugen;

static int g_failures = 0;

static void check(const char* what, double got, double expect, double rtol,
                  const char* note = "") {
  const double rel = std::fabs(got - expect) / std::fabs(expect);
  const bool ok = rel <= rtol;
  if (!ok) ++g_failures;
  std::printf("  [%s] %-44s got=%.6g  expect=%.6g  rel=%.2e (tol %.0e) %s\n",
              ok ? "ok  " : "FAIL", what, got, expect, rel, rtol, note);
}

/// For counts, where a relative tolerance would divide by zero.
static void checkExact(const char* what, long long got, long long expect) {
  const bool ok = got == expect;
  if (!ok) ++g_failures;
  std::printf("  [%s] %-44s got=%lld  expect=%lld\n",
              ok ? "ok  " : "FAIL", what, got, expect);
}

// Geometry shared by the far-field tests. The detector sits directly below a
// small source disk, far enough that it subtends a narrow cone: then the solid
// angle is its projected area over D^2 and the obliquity factor is 1.
//
// D is deliberately enormous relative to the detectors. A detector has depth,
// so the face that actually intercepts the rays is nearer than its centre by
// its own half-thickness, and the rate goes as 1/distance^2: at D = 200 m a
// 6 m cylinder is off by (20000/19700)^2 = 3%, which would swamp everything
// being tested here. At 20 km the same effect is 3e-4. The sphere is the one
// shape immune to this, since 2*pi*(1-cos asin(R/D)) is exact at any distance.
static constexpr double kD = 2.0e6;        // detector depth, cm
static constexpr double kSrcR = 50.0;      // source disk radius, cm

static Generator makeGen(std::shared_ptr<Detector> det) {
  Generator g;
  // CosmoALEPH has no zenith dependence, so J(p,c)/J(p,1) = 1 and the rate is
  // pure geometry. That is what makes a closed form available at all.
  g.setSpectrum(Spectrum::CosmoALEPH)
      .setSurface(std::make_shared<Disk>(kSrcR, Vec3{0, 0, 0}))
      .setEnergyRange(1.0, 1000.0)
      .setDetector(std::move(det));
  return g;
}

static void farField(const char* label, std::shared_ptr<Detector> det,
                     double projected_area_cm2, int npoints = 4000000) {
  Generator g = makeGen(std::move(det));
  const double expect =
      g.surfaceArea() * g.cdf().integrated_flux() * projected_area_cm2 / (kD * kD);
  double r, err;
  g.rateAndError(r, err, npoints);
  char note[64];
  std::snprintf(note, sizeof(note), "[MC err %.1e]", err / r);
  check(label, r, expect, 8e-3, note);
}

int main() {
  std::printf("1. far field: rate is set by projected area (closed form)\n");
  // A sphere of radius R projects to pi*R^2 from every direction.
  farField("sphere, A_proj = pi R^2",
           std::make_shared<SphereDetector>(Vec3{0, 0, -kD}, 200.0),
           kPi * 200.0 * 200.0);

  // A box seen from directly above projects to its top face.
  farField("box, A_proj = 4 hx hy",
           BoxDetector::centred(Vec3{0, 0, -kD}, 150.0, 250.0, 100.0),
           4.0 * 150.0 * 250.0);

  // A vertical cylinder seen from above projects to its cap, pi r^2. Its
  // length must not matter, which is the point of the cap tests.
  farField("cylinder upright, A_proj = pi r^2",
           std::make_shared<CylinderDetector>(Vec3{0, 0, -kD - 300.0},
                                              Vec3{0, 0, -kD + 300.0}, 200.0),
           kPi * 200.0 * 200.0);

  // Laid on its side the same cylinder projects to a rectangle 2r by L, so
  // this is the test that the arbitrary-orientation path is right.
  farField("cylinder on its side, A_proj = 2 r L",
           std::make_shared<CylinderDetector>(Vec3{-300.0, 0, -kD},
                                              Vec3{300.0, 0, -kD}, 200.0),
           2.0 * 200.0 * 600.0);

  std::printf("\n2. directed sampling == plain acceptance cut\n");
  {
    auto det = [] {
      return std::make_shared<CylinderDetector>(Vec3{120.0, -80.0, -6000.0},
                                                Vec3{120.0, -80.0, -5400.0},
                                                250.0);
    };
    Generator fast = makeGen(det());
    Generator slow = makeGen(det());
    slow.setDirectedSampling(false);

    double rf, ef, rs, es;
    fast.rateAndError(rf, ef, 4000000);
    slow.rateAndError(rs, es, 40000000);   // blind needs far more samples
    // Compare against the combined statistical error of the two estimates.
    const double sigma = std::sqrt(ef * ef + es * es);
    const double tol = 4.0 * sigma / rs;
    char note[80];
    std::snprintf(note, sizeof(note), "[4 sigma = %.1e, off-axis detector]", tol);
    check("rate agrees within 4 sigma", rf, rs, tol, note);

    // Same question asked of the generated events rather than the integrator.
    Generator g = makeGen(det());
    long long n = 0;
    for (int i = 0; i < 20000; ++i) {
      const Muon m = g.generate();
      if (det()->intersects(m.position, m.direction)) ++n;
    }
    checkExact("every generated muon hits the detector", n, 20000);
    checkExact("no envelope violations", g.envelopeViolations(), 0);
  }

  std::printf("\n3. a detector that intercepts everything == no detector\n");
  {
    // A slab 1 cm below the source, wide enough that no ray inside the
    // 60 degree window can escape it. Its bounding sphere swallows the source
    // disk, so this also exercises the degenerate full-sphere cone.
    auto build = [](std::shared_ptr<Detector> det) {
      Generator g;
      g.setSpectrum(Spectrum::Guan)
          .setSurface(std::make_shared<Disk>(kSrcR, Vec3{0, 0, 0}))
          .setEnergyRange(1.0, 1000.0)
          .setThetaRange(0.0, 60.0 * kPi / 180.0)
          .setDetector(std::move(det));
      return g;
    };
    Generator plain = build(nullptr);
    Generator caught =
        build(BoxDetector::centred(Vec3{0, 0, -50.0}, 1.0e7, 1.0e7, 49.0));

    double r0, e0, r1, e1;
    plain.rateAndError(r0, e0, 2000000);
    caught.rateAndError(r1, e1, 2000000);
    const double tol = 4.0 * std::sqrt(e0 * e0 + e1 * e1) / r0;
    char note[80];
    std::snprintf(note, sizeof(note), "[4 sigma = %.1e]", tol);
    check("intercepting slab reproduces open rate", r1, r0, tol, note);
    std::printf("  [info] %-44s %.3f sr (full sphere = %.3f)\n",
                "cone solid angle used", caught.coneSolidAngleMax(),
                4.0 * kPi);
  }

  std::printf("\n4. what it buys\n");
  {
    // A realistic muography setup rather than the 20 km used above: a 2 m
    // detector under 200 m of rock, generated over a 40 m disk. The blind mode
    // is quadratically worse in D, so it is given fewer muons and the
    // comparison is made per muon.
    auto build = [](bool directed) {
      Generator g;
      g.setSpectrum(Spectrum::Guan)
          .setSurface(std::make_shared<Disk>(4000.0, Vec3{0, 0, 0}))
          .setEnergyRange(1.0, 1000.0)
          .setDetector(std::make_shared<SphereDetector>(Vec3{0, 0, -20000.0},
                                                        200.0))
          .setDirectedSampling(directed);
      return g;
    };
    double per_muon[2] = {0.0, 0.0};
    int idx = 0;
    for (bool directed : {false, true}) {
      Generator g = build(directed);
      const int n = directed ? 200000 : 2000;
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < n; ++i) g.generate();
      const double ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - t0).count();
      per_muon[idx++] = ms / n;
      std::printf("  [info] %-16s %9.4f ms/muon   acceptance = %.3g%%\n",
                  directed ? "directed" : "acceptance cut", ms / n,
                  100.0 * g.acceptance());
    }
    std::printf("  [info] %-16s %9.0fx\n", "speed-up",
                per_muon[0] / per_muon[1]);
  }

  std::printf("\n%s (%d failure%s)\n",
              g_failures ? "FAILURES PRESENT" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
