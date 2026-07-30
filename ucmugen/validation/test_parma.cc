// Validation of the PARMA spectrum.
//
// UCMuGen_PARMA.h is generated from JAEA's own subroutines.cpp with three
// mechanical edits (see ucmugen/tools/make_parma_header.py). Two of those edits
// restrict species loops to muons so that only the muon tables need embedding,
// and that is the one place where the generated header could silently diverge
// from the stock build. It does not: at generation time the two are compared
// over 129,600 points spanning every W index, cutoff rigidity, atmospheric
// depth, local-geometry mode, energy and zenith angle in the model's range, and
// they agree bit-for-bit.
//
// Run that comparison under BOTH compilers, not one. Doing it only under clang
// hid a real bug: the stock source calls unqualified abs() on a double, and
// libc++ hoists the double overload into the global namespace while libstdc++
// does not, so under gcc it silently bound to C's integer abs and truncated
// -135.4 to -135. That moved getFFPfromW by a factor of nine and the flux by
// 2e-4, quietly, on one compiler only. The generated header now names the std
// overloads explicitly. Current state: bit-identical to stock under clang, and
// within 1 ULP (6e-16) under gcc, against 7.3e-13 between the two stock builds
// themselves.
//
// That comparison needs an unpacked parma_cpp, so it cannot run here. What runs
// here instead:
//
//   1. reference values lifted from the stock JAEA build, so a regeneration
//      that changed the physics would be caught without parma_cpp present
//   2. the site dependence that is the entire reason for wanting PARMA:
//      altitude, geomagnetic cutoff, solar modulation
//   3. the wiring into UCMuGen: unit conversion, Spectrum::Parma dispatch,
//      and the absolute rate that follows from it
//
// Build:
//   c++ -std=c++17 -O2 -I../include test_parma.cc -o test_parma

#include "UCMuGen.h"
#include "UCMuGen_PARMA.h"

#include <cmath>
#include <cstdio>
#include <memory>

using namespace ucmugen;

static int g_failures = 0;

static void check(const char* what, double got, double expect, double rtol) {
  const double rel = std::fabs(got - expect) / std::fabs(expect);
  const bool ok = rel <= rtol;
  if (!ok) ++g_failures;
  std::printf("  [%s] %-42s got=%.10g  expect=%.10g  rel=%.1e\n",
              ok ? "ok  " : "FAIL", what, got, expect, rel);
}

static void expect_true(const char* what, bool ok, const char* detail = "") {
  if (!ok) ++g_failures;
  std::printf("  [%s] %-42s %s\n", ok ? "ok  " : "FAIL", what, detail);
}

int main() {
  // -------------------------------------------------------------------------
  // 1. Values from the stock JAEA build.
  //
  // Produced by linking the official subroutines.cpp directly and printing
  // getMuonSpec(mu+) + getMuonSpec(mu-) times getSpecAngFinal, at 17 digits.
  // The tolerance is 1e-12 rather than 0 only to absorb the difference between
  // an -O2 fused multiply-add here and there; in practice these match exactly.
  // -------------------------------------------------------------------------
  std::printf("1. against the stock JAEA build (flux in /cm2/s/MeV/sr)\n");
  struct Ref { double w, r, d, e, g, c, expect; const char* tag; };
  static const Ref refs[] = {
    {0, 3.25, 1033.23, 1000, 0.15, 1, 2.2825743560210692e-06, "sea level, 1 GeV, vertical"},
    {0, 3.25, 1033.23, 10000, 0.15, 1, 1.1983011755760992e-07, "sea level, 10 GeV, vertical"},
    {0, 3.25, 1033.23, 100000, 0.15, 1, 2.8643635719767066e-10, "sea level, 100 GeV, vertical"},
    {0, 3.25, 1033.23, 10000, 0.15, 0.5, 4.1079780824970483e-08, "sea level, 10 GeV, 60 deg"},
    {0, 3.25, 700, 10000, 0.15, 1, 1.5681269591765259e-07, "3 km altitude, 10 GeV, vertical"},
    {150, 3.25, 1033.23, 10000, 0.15, 1, 1.11877633527772e-07, "solar max, 10 GeV, vertical"},
    {0, 15, 1033.23, 10000, 0.15, 1, 1.2029477447928929e-07, "high cutoff, 10 GeV, vertical"},
  };
  for (const Ref& r : refs) {
    namespace d = parma::detail;
    const double sp = d::getMuonSpecCpp(1, r.w, r.r, r.d, r.e)
                    + d::getMuonSpecCpp(2, r.w, r.r, r.d, r.e);
    const double an = d::getSpecAngFinalCpp(4, r.w, r.r, r.d, r.e, r.g, r.c);
    check(r.tag, sp * an, r.expect, 1e-12);
  }

  // The embedded altitude table, likewise from the stock build.
  std::printf("\n2. altitude to atmospheric depth (embedded AtomDepth table)\n");
  check("0 km", parma::depth_from_altitude(0.0, 50.67), 1033.227453, 1e-12);
  check("1 km", parma::depth_from_altitude(1.0, 50.67), 916.48014360000002, 1e-12);
  check("3 km", parma::depth_from_altitude(3.0, 50.67), 715.03520570000001, 1e-12);
  check("5 km", parma::depth_from_altitude(5.0, 50.67), 551.13621880000005, 1e-12);

  // -------------------------------------------------------------------------
  // 3. The site dependence, which is the whole point: a sea-level
  //    parametrisation cannot express any of these.
  // -------------------------------------------------------------------------
  std::printf("\n3. site dependence (10 GeV/c vertical)\n");
  auto site_at = [](double depth, double rc, double w) {
    parma::Site s;
    s.depth_gcm2 = depth; s.cutoff_GV = rc; s.w_index = w;
    return s;
  };
  const double p = 10.0;
  const double f_sea = parma::intensity(site_at(1033.23, 3.25, 0.0), p, 1.0);
  const double f_3km = parma::intensity(site_at(715.04, 3.25, 0.0), p, 1.0);
  const double f_5km = parma::intensity(site_at(551.14, 3.25, 0.0), p, 1.0);
  std::printf("  [info] sea level %.4g, 3 km %.4g, 5 km %.4g /cm2/s/sr/(GeV/c)\n",
              f_sea, f_3km, f_5km);
  expect_true("flux rises with altitude", f_3km > f_sea && f_5km > f_3km);
  std::printf("  [info] 3 km is %.2fx sea level, 5 km is %.2fx\n",
              f_3km / f_sea, f_5km / f_sea);

  // The geomagnetic cutoff is a low-energy effect, and asserting otherwise is
  // a trap: at 10 GeV/c the flux is flat in Rc to about 1% and not even
  // monotonic, because removing low-rigidity primaries hardens the surviving
  // spectrum. So test it where it acts, at 1 GeV/c, and on the integral.
  auto integrated = [&](double rc) {
    const parma::Site s = site_at(1033.23, rc, 0.0);
    double sum = 0.0;
    constexpr int kN = 2000;
    for (int i = 0; i < kN; ++i) {
      const double lo = 0.1 * std::pow(10.0, 4.0 * i / kN);
      const double hi = 0.1 * std::pow(10.0, 4.0 * (i + 1) / kN);
      sum += 0.5 * (parma::intensity(s, lo, 1.0) + parma::intensity(s, hi, 1.0))
             * (hi - lo);
    }
    return sum;
  };
  const double i_low = integrated(1.0), i_mid = integrated(8.0),
               i_high = integrated(20.0);
  std::printf("  [info] vertical integral 0.1-1000 GeV/c: Rc=1 %.5g, "
              "Rc=8 %.5g, Rc=20 %.5g\n", i_low, i_mid, i_high);
  expect_true("integrated flux falls with cutoff rigidity",
              i_high < i_mid && i_mid < i_low);
  expect_true("1 GeV/c flux falls with cutoff rigidity",
              parma::intensity(site_at(1033.23, 20.0, 0.0), 1.0, 1.0) <
              parma::intensity(site_at(1033.23, 1.0, 0.0), 1.0, 1.0));
  // PARMA clamps Rc below 1 GV: muon fluxes are identical there by
  // construction (getMuonSpec does r = max(1.0, r1)).
  expect_true("cutoff below 1 GV is clamped, as PARMA documents",
              parma::intensity(site_at(1033.23, 0.1, 0.0), 1.0, 1.0) ==
              parma::intensity(site_at(1033.23, 1.0, 0.0), 1.0, 1.0));

  const double f_smin = parma::intensity(site_at(1033.23, 3.25, 0.0), p, 1.0);
  const double f_smax = parma::intensity(site_at(1033.23, 3.25, 150.0), p, 1.0);
  expect_true("solar maximum suppresses the flux", f_smax < f_smin);

  const double f_horiz = parma::intensity(site_at(1033.23, 3.25, 0.0), p, 0.1);
  expect_true("flux falls towards the horizon at 10 GeV", f_horiz < f_sea);

  // -------------------------------------------------------------------------
  // 4. Wiring into the generator.
  // -------------------------------------------------------------------------
  std::printf("\n4. Spectrum::Parma through the generator\n");
  expect_true("throws until a provider is installed", [] {
    parma::uninstall();
    try {
      flux::intensity(Spectrum::Parma, 10.0, 1.0);
      return false;
    } catch (const std::invalid_argument&) { return true; }
  }());

  parma::Site site;
  site.cutoff_GV = 3.25;
  site.depth_gcm2 = parma::depth_from_altitude(0.0, 50.67);
  parma::install(site);

  check("flux::intensity matches parma::intensity",
        flux::intensity(Spectrum::Parma, 10.0, 1.0),
        parma::intensity(site, 10.0, 1.0), 1e-15);

  auto make = [](const parma::Site& s) {
    parma::install(s);
    Generator g;
    g.setSpectrum(Spectrum::Parma)
        .setSurface(std::make_shared<Plane>(100.0, 100.0, Vec3{0, 0, 0}))
        .setEnergyRange(1.0, 1000.0)
        .setThetaRange(0.0, 70.0 * kPi / 180.0);
    return g;
  };

  Generator gen = make(site);
  const double rate_sea = gen.rate(200000);
  std::printf("  [info] rate through 2x2 m at sea level: %.4g Hz\n", rate_sea);
  expect_true("rate is finite and positive",
              std::isfinite(rate_sea) && rate_sea > 0.0);

  // A 1 m^2 horizontal detector sees very roughly 100 muons/s at sea level.
  // This is an order-of-magnitude guard against a unit slip, not a physics
  // measurement: the surface here is 4 m^2 and the acceptance is capped at 70
  // degrees, so anything within a factor of a few of 400 Hz is reasonable.
  expect_true("rate is the right order of magnitude",
              rate_sea > 40.0 && rate_sea < 4000.0, "[expect O(100) Hz/m^2]");

  parma::Site alt = site;
  alt.depth_gcm2 = parma::depth_from_altitude(3.0, 50.67);
  Generator gen_alt = make(alt);
  const double rate_3km = gen_alt.rate(200000);
  std::printf("  [info] same detector at 3 km: %.4g Hz (%.2fx)\n",
              rate_3km, rate_3km / rate_sea);
  expect_true("altitude raises the generated rate", rate_3km > rate_sea);

  // Generated muons must be sane, not merely numerous.
  parma::install(site);
  Generator g2 = make(site);
  long long bad = 0, mu_plus = 0;
  constexpr int kN = 20000;
  for (int i = 0; i < kN; ++i) {
    const Muon m = g2.generate();
    if (!(m.momentum > 0.0) || !std::isfinite(m.momentum)) ++bad;
    if (m.direction.z >= 0.0) ++bad;              // must be downward
    if (std::abs(m.pdg) != 13) ++bad;
    if (m.pdg == -13) ++mu_plus;
  }
  expect_true("20k generated muons are all well formed", bad == 0);
  const double ratio = double(mu_plus) / double(kN - mu_plus);
  std::printf("  [info] mu+/mu- from UCMuGen's ratio = %.3f; PARMA's own at "
              "10 GeV/c = %.3f\n", ratio, parma::charge_ratio(site, 10.0));

  std::printf("\n%s (%d failure%s)\n",
              g_failures ? "FAILURES PRESENT" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
