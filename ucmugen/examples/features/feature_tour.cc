// A tour of every UCMuGen feature, in one runnable file.
//
// The Geant4 example next door is deliberately minimal: one spectrum, one
// surface, one detector, because its job is to show the integration and be
// copied. This file is the opposite. It exercises the whole API and prints what
// each option does to the numbers, so you can see the effect of a choice before
// making it in your own application.
//
// It needs no Geant4 and no external data:
//
//   c++ -std=c++17 -O2 -I../../include feature_tour.cc -o feature_tour
//   ./feature_tour
//
// With the PARMA header generated (see ucmugen/tools/make_parma_header.py):
//
//   c++ -std=c++17 -O2 -I../../include -DUCMUGEN_TOUR_PARMA \
//       feature_tour.cc -o feature_tour
//
// Everything is in UCMuGen's units throughout: centimetres, GeV, seconds.

#include "UCMuGen.h"

#ifdef UCMUGEN_TOUR_PARMA
#include "UCMuGen_PARMA.h"
#endif

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace ucmugen;

namespace {

// Rates are Monte Carlo integrals. Fewer points here than the default, because
// the tour computes a lot of them and 1% is plenty to see the differences.
constexpr int kRatePoints = 60000;

void heading(const char* title) {
  std::printf("\n=============================================================="
              "==\n  %s\n===================================================="
              "============\n", title);
}

/// A generator aimed at a 2 m x 2 m x 20 cm plate at the origin from a 6 m
/// sky plane 3 m up. The same setup as the Geant4 example, so the numbers
/// printed here are directly comparable with the numbers it prints.
Generator baseline() {
  Generator g;
  g.setSpectrum(Spectrum::Guan)
      .setSurface(std::make_shared<Plane>(300.0, 300.0, Vec3{0.0, 0.0, 300.0}))
      .setEnergyRange(1.0, 1000.0)
      .setThetaRange(0.0, 70.0 * kPi / 180.0)
      .setDetector(BoxDetector::centred(Vec3{0.0, 0.0, 0.0}, 100.0, 100.0, 10.0));
  return g;
}

}  // namespace

// ---------------------------------------------------------------------------
// 1. The spectra
// ---------------------------------------------------------------------------
// Eight of them. Six carry an absolute normalisation and so give a rate in Hz;
// PowerLaw and Electron are shapes only, and rate() returns a negative number
// for those rather than a meaningless one. That distinction is the reason to
// check the sign instead of printing whatever comes back.
void tour_spectra() {
  heading("1. Spectra");

  struct Entry { Spectrum s; const char* name; const char* note; };
  const std::vector<Entry> entries = {
      {Spectrum::CosmoALEPH,    "CosmoALEPH",    "cosmoALEPH fit"},
      {Spectrum::PowerLaw,      "PowerLaw",      "shape only, no normalisation"},
      {Spectrum::Guan,          "Guan",          "modified Gaisser, 2015"},
      {Spectrum::Frosin,        "Frosin",        "Guan form, refitted 2025"},
      {Spectrum::GaisserBugaev, "GaisserBugaev", "no atmospheric correction"},
      {Spectrum::ReynaBugaev,   "ReynaBugaev",   "Reyna 2006 parametrisation"},
      {Spectrum::Electron,      "Electron",      "cosmic e+/e-, shape only"},
  };

  std::printf("\n  %-16s %14s   %s\n", "spectrum", "rate [Hz]", "note");
  std::printf("  %-16s %14s   %s\n", "----------------", "--------------",
              "----");
  for (const auto& e : entries) {
    Generator g = baseline();
    g.setSpectrum(e.s);
    const double r = g.rate(kRatePoints);
    if (r < 0.0)
      std::printf("  %-16s %14s   %s\n", e.name, "--", e.note);
    else
      std::printf("  %-16s %14.3f   %s\n", e.name, r, e.note);
  }

#ifdef UCMUGEN_TOUR_PARMA
  // PARMA is the eighth. It is the only spectrum that knows where and when you
  // are, so it is the only one whose rate moves when you change altitude.
  std::printf("\n  Parma (site-aware), same geometry:\n");
  for (const double alt_km : {0.0, 1.0, 3.0}) {
    parma::Site site;
    site.cutoff_GV = 3.25;
    site.depth_gcm2 = parma::depth_from_altitude(alt_km, 50.67);
    parma::install(site);
    Generator g = baseline();
    g.setSpectrum(Spectrum::Parma);
    std::printf("    altitude %4.1f km  depth %7.1f g/cm2   rate %10.3f Hz\n",
                alt_km, site.depth_gcm2, g.rate(kRatePoints));
  }
  parma::uninstall();
#else
  std::printf("\n  Parma: build with -DUCMUGEN_TOUR_PARMA to include it.\n");
#endif
}

// ---------------------------------------------------------------------------
// 2. The generation surfaces
// ---------------------------------------------------------------------------
// This is the part with no equivalent in a plain spectrum sampler. The
// projection cos factor lives in the Surface, so a tilted plane, a dome and a
// cylinder are all first-class: none of them needs a hand-derived integrator.
//
// The rates below are NOT expected to match each other. Each surface has its
// own area and its own view of the sky, and each is aimed at the same plate, so
// what should match is any pair that fully encloses the detector.
void tour_surfaces() {
  heading("2. Generation surfaces");

  struct Entry { std::shared_ptr<Surface> s; const char* name; const char* note; };
  const std::vector<Entry> entries = {
      {std::make_shared<Plane>(300.0, 300.0, Vec3{0, 0, 300}),
       "Plane (flat)", "6 m x 6 m sky plane, normal +Z"},
      {std::make_shared<Plane>(300.0, 300.0, Vec3{0, 0, 300},
                               Vec3{0.0, std::sin(0.5), std::cos(0.5)}),
       "Plane (tilted)", "same plane tilted 28.6 deg"},
      {std::make_shared<Disk>(300.0, Vec3{0, 0, 300}),
       "Disk", "3 m radius, normal +Z"},
      {std::make_shared<HSphere>(500.0, Vec3{0, 0, 0}),
       "HSphere", "5 m dome over the detector"},
      {std::make_shared<Cylinder>(400.0, 600.0, Vec3{0, 0, 0}, true),
       "Cylinder", "4 m radius, 6 m tall, with caps"},
  };

  std::printf("\n  %-16s %12s %14s   %s\n", "surface", "area [cm2]",
              "rate [Hz]", "note");
  std::printf("  %-16s %12s %14s   %s\n", "----------------", "------------",
              "--------------", "----");
  for (const auto& e : entries) {
    Generator g = baseline();
    g.setSurface(e.s);
    // A dome and a cylinder see the whole sky, including muons arriving from
    // below the horizon of a flat plane, so open the zenith range for them.
    g.setThetaRange(0.0, kPi / 2.0);
    std::printf("  %-16s %12.3g %14.3f   %s\n", e.name, e.s->area(),
                g.rate(kRatePoints), e.note);
  }

  std::printf("\n  A tilted plane is the case worth noticing: the Fortran\n"
              "  generator applies no projection weighting to a tilted or\n"
              "  vertical surface, so its sky intensity is not a\n"
              "  through-surface flux there. Here the cos factor comes from\n"
              "  the surface itself, so the tilt is handled by construction.\n");
}

// ---------------------------------------------------------------------------
// 3. The detector shapes, and directed sampling
// ---------------------------------------------------------------------------
// A detector does two things: it restricts what is generated to muons whose ray
// reaches it, and it turns rate() into the rate *into the detector*.
//
// Directed sampling is the speed-up. Instead of drawing over the whole sky and
// discarding almost everything, directions come from the cone the detector's
// bounding sphere subtends from each sampled surface point. It is not an
// approximation: the two modes must agree on rate() to within Monte Carlo
// error, which is what the second table checks.
void tour_detectors() {
  heading("3. Detector shapes, and directed sampling");

  struct Entry { std::shared_ptr<Detector> d; const char* name; };
  const std::vector<Entry> entries = {
      {BoxDetector::centred(Vec3{0, 0, 0}, 100.0, 100.0, 10.0), "BoxDetector"},
      {std::make_shared<SphereDetector>(Vec3{0, 0, 0}, 100.0), "SphereDetector"},
      {std::make_shared<CylinderDetector>(Vec3{0, 0, -50}, Vec3{0, 0, 50}, 100.0),
       "CylinderDetector"},
  };

  std::printf("\n  %-18s %14s %14s %12s\n", "detector", "directed [Hz]",
              "blind [Hz]", "agree?");
  std::printf("  %-18s %14s %14s %12s\n", "------------------", "--------------",
              "--------------", "------------");
  for (const auto& e : entries) {
    Generator gd = baseline();
    gd.setDetector(e.d).setDirectedSampling(true);
    double rd = 0.0, ed = 0.0;
    gd.rateAndError(rd, ed, kRatePoints);

    Generator gb = baseline();
    gb.setDetector(e.d).setDirectedSampling(false);
    double rb = 0.0, eb = 0.0;
    gb.rateAndError(rb, eb, kRatePoints);

    // Compare against the combined 1-sigma error of the two estimates. Both
    // are Monte Carlo, so the test is statistical, not exact equality.
    const double sigma = std::sqrt(ed * ed + eb * eb);
    const bool ok = std::fabs(rd - rb) < 4.0 * sigma;
    std::printf("  %-18s %8.3f+-%-4.2f %8.3f+-%-4.2f %12s\n", e.name, rd, ed,
                rb, eb, ok ? "yes" : "NO");
  }

  // The reason to bother. Acceptance is the fraction of proposals that
  // survive, so the ratio of the two is roughly the saving in wall time.
  //
  // The size of the win is not a property of the generator, it is a property
  // of how much of the sky the detector covers. A plate that fills the view
  // from the source plane gains little, because blind sampling was already
  // hitting it most of the time. Shrink the detector and the gain grows
  // without limit, which is exactly the regime real muography lives in.
  std::printf("\n  %-22s %12s %12s %10s\n", "detector half-size",
              "directed", "blind", "speed-up");
  std::printf("  %-22s %12s %12s %10s\n", "----------------------",
              "------------", "------------", "----------");
  for (const double h : {100.0, 30.0, 10.0, 3.0, 1.0}) {
    auto det = BoxDetector::centred(Vec3{0, 0, 0}, h, h, 10.0);

    Generator gd = baseline();
    gd.setDetector(det).setDirectedSampling(true);
    for (int i = 0; i < 20000; ++i) gd.generate();

    // Blind sampling on a small detector is slow by construction, which is the
    // point being made, so ask it for only a few hundred muons. acceptance()
    // is accepted/tried, so a few hundred is already a stable estimate.
    Generator gb = baseline();
    gb.setDetector(det).setDirectedSampling(false);
    for (int i = 0; i < 300; ++i) gb.generate();
    const double acc_blind = gb.acceptance();

    char buf[32];
    std::snprintf(buf, sizeof buf, "%.0f x %.0f cm", 2 * h, 2 * h);
    std::printf("  %-22s %12.2e %12.2e %9.0fx\n", buf, gd.acceptance(),
                acc_blind, gd.acceptance() / acc_blind);
  }
}

// ---------------------------------------------------------------------------
// 4. Angular ranges
// ---------------------------------------------------------------------------
// setThetaRange restricts zenith, setPhiRange restricts azimuth. Both are
// honoured by the rate integral as well as by generate(), so a restricted
// generator still reports the rate for the region it actually covers, and not
// some fraction of a full-sky number.
void tour_angles() {
  heading("4. Angular ranges");

  std::printf("\n  %-34s %14s\n", "range", "rate [Hz]");
  std::printf("  %-34s %14s\n", "----------------------------------",
              "--------------");

  struct Entry { double t0, t1, p0, p1; const char* name; };
  const std::vector<Entry> entries = {
      {0.0, 70.0 * kPi / 180.0, 0.0, 2 * kPi, "theta 0-70 deg, all azimuth"},
      {0.0, 30.0 * kPi / 180.0, 0.0, 2 * kPi, "theta 0-30 deg, all azimuth"},
      {0.0, 10.0 * kPi / 180.0, 0.0, 2 * kPi, "theta 0-10 deg (near vertical)"},
      {30.0 * kPi / 180.0, 70.0 * kPi / 180.0, 0.0, 2 * kPi,
       "theta 30-70 deg (annulus)"},
      {0.0, 70.0 * kPi / 180.0, 0.0, kPi, "theta 0-70, azimuth half-sky"},
      {0.0, 70.0 * kPi / 180.0, 0.0, kPi / 2.0, "theta 0-70, azimuth quadrant"},
  };
  for (const auto& e : entries) {
    Generator g = baseline();
    g.setThetaRange(e.t0, e.t1).setPhiRange(e.p0, e.p1);
    std::printf("  %-34s %14.3f\n", e.name, g.rate(kRatePoints));
  }

  std::printf("\n  Halving the azimuth range roughly halves the rate, as it\n"
              "  must for a symmetric setup. Restricting zenith cuts much\n"
              "  harder than solid angle alone, because the flux falls with\n"
              "  zenith and the projection factor falls with it too.\n");
}

// ---------------------------------------------------------------------------
// 5. Random numbers and reproducibility
// ---------------------------------------------------------------------------
// Two ways to control the stream. setSeed uses UCMuGen's own LCG; setRng takes
// any Rng, which is how the Geant4 example hands control to G4UniformRand so a
// run is reproducible from the Geant4 seed.
void tour_rng() {
  heading("5. Random numbers and reproducibility");

  auto first_five = [](Generator& g) {
    std::string s;
    char buf[32];
    for (int i = 0; i < 5; ++i) {
      std::snprintf(buf, sizeof buf, "%.4f ", g.generate().momentum);
      s += buf;
    }
    return s;
  };

  Generator a = baseline(); a.setSeed(12345);
  Generator b = baseline(); b.setSeed(12345);
  Generator c = baseline(); c.setSeed(999);
  const std::string sa = first_five(a), sb = first_five(b), sc = first_five(c);

  std::printf("\n  seed 12345      : %s\n", sa.c_str());
  std::printf("  seed 12345 again: %s  -> %s\n", sb.c_str(),
              sa == sb ? "identical, as it must be" : "MISMATCH");
  std::printf("  seed 999        : %s  -> %s\n", sc.c_str(),
              sa != sc ? "different, as it must be" : "MISMATCH");

  // An external engine. Any callable returning [0,1) works; in Geant4 this is
  // G4UniformRand, which is what makes /random/setSeeds control the primaries.
  UCMuonLcg external(4242);
  Generator d = baseline();
  d.setRng(std::unique_ptr<Rng>(new CallbackRng([&external] {
    return external.flat();
  })));
  std::printf("  external engine : %s\n", first_five(d).c_str());
  std::printf("\n  provenance: %s\n", d.provenance().c_str());
}

// ---------------------------------------------------------------------------
// 6. Normalisation, exposure, and what a run is worth
// ---------------------------------------------------------------------------
// The practical payoff. rate() converts a count of simulated muons into a real
// exposure, which is what lets a simulated histogram be compared with data.
void tour_normalisation() {
  heading("6. Normalisation and exposure");

  Generator g = baseline();
  double r = 0.0, err = 0.0;
  g.rateAndError(r, err, 200000);

  std::printf("\n  rate into the detector : %.3f +- %.3f Hz  (%.2f %%)\n", r,
              err, 100.0 * err / r);
  std::printf("\n  %14s %16s %16s\n", "muons", "live time [s]", "live time");
  std::printf("  %14s %16s %16s\n", "--------------", "----------------",
              "----------------");
  for (const long long n : {1000LL, 100000LL, 1000000LL, 10000000LL}) {
    const double t = n / r;
    char human[64];
    if (t < 120.0)        std::snprintf(human, sizeof human, "%.1f s", t);
    else if (t < 7200.0)  std::snprintf(human, sizeof human, "%.1f min", t / 60.0);
    else                  std::snprintf(human, sizeof human, "%.2f h", t / 3600.0);
    std::printf("  %14lld %16.3f %16s\n", n, t, human);
  }

  std::printf("\n  A generator without this number can tell you the shape of\n"
              "  a distribution but not how long you waited for it.\n");
}

// ---------------------------------------------------------------------------
// 7. The legacy driver
// ---------------------------------------------------------------------------
// A second, separate API whose contract is bit-exact reproduction of the
// Fortran generator, including its conventions and its quirks. Use it to
// cross-check against archived Fortran output; use the projection-correct
// Generator above for new work.
void tour_legacy() {
  heading("7. The legacy Fortran-compatible driver");

  legacy::Config cfg;
  cfg.spectrum = Spectrum::Guan;
  cfg.angular = Angular::Cos2;
  cfg.e_min = 1.0;
  cfg.e_max = 1000.0;
  cfg.theta_max_deg = 70.0;
  cfg.source = legacy::Source::Disk;
  cfg.radius_cm = 200.0;
  cfg.seed = 20260729;

  legacy::Generator g(cfg);
  std::printf("\n  %10s %12s %12s %10s %8s\n", "i", "p [GeV/c]", "E [GeV]",
              "theta[deg]", "pdg");
  std::printf("  %10s %12s %12s %10s %8s\n", "----------", "------------",
              "------------", "----------", "--------");
  for (int i = 0; i < 5; ++i) {
    const legacy::Muon m = g.generate();
    std::printf("  %10d %12.4f %12.4f %10.3f %8d\n", i, m.p, m.energy,
                m.theta * 180.0 / kPi, m.pdg);
  }

  std::printf("\n  The five angular models (Vertical, Cos2, UniformCone,\n"
              "  GuanSelfConsistent, Cos3) belong to this driver. The\n"
              "  projection-correct Generator does not need them: its angular\n"
              "  distribution follows from the flux and the surface.\n");
}

int main() {
  std::printf("UCMuGen feature tour\n");
  std::printf("Units throughout: cm, GeV, s.\n");

  tour_spectra();
  tour_surfaces();
  tour_detectors();
  tour_angles();
  tour_rng();
  tour_normalisation();
  tour_legacy();

  std::printf("\nDone.\n");
  return 0;
}
