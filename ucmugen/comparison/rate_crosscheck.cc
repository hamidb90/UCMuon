// Item 2: absolute normalisation, UCMuGen vs EcoMug.
//
// The two tools use different built-in flux models, so comparing their default
// rates would compare parametrisations, not integrators. Instead both are given
// the IDENTICAL differential flux (UCMuGen's Guan model, handed to EcoMug
// through SetDifferentialFlux) over the same momentum and zenith ranges. Any
// remaining difference is then a difference of integration, which is the thing
// worth cross-checking: two independently written estimators agreeing on the
// same input validates both.
//
// EcoMug's custom-J sky integrator computes
//     rate = V * <J(p,theta) cos(theta) sin(theta)>,  V = dp dtheta dphi
// which is the flux through unit horizontal area, in whatever area unit J uses.
// UCMuGen's rate() is the rate through the whole surface, so it is divided by
// the area to compare. J is in cm^-2 s^-1 sr^-1 (GeV/c)^-1 throughout.

#include "UCMuGen.h"
#include "EcoMug.h"

#include <cmath>
#include <cstdio>
#include <memory>

using namespace ucmugen;

int main() {
  const double p_min = 1.0, p_max = 1000.0;          // GeV/c
  const double th_max_deg = 70.0;
  const double th_max = th_max_deg * kPi / 180.0;

  // The shared physics: UCMuGen's Guan parametrisation.
  auto J = [](double p_GeV, double theta_rad) {
    return flux::intensity(Spectrum::Guan, p_GeV, std::cos(theta_rad));
  };

  // ---- UCMuGen -------------------------------------------------------------
  // setEnergyRange takes TOTAL energy, so convert the momentum limits.
  const double m = kMuonMass;
  const double e_min = std::sqrt(p_min * p_min + m * m);
  const double e_max = std::sqrt(p_max * p_max + m * m);
  const double half = 100.0;                          // cm, 2 m x 2 m plane
  Generator g;
  g.setSpectrum(Spectrum::Guan)
      .setSurface(std::make_shared<Plane>(half, half, Vec3{0, 0, 0}))
      .setEnergyRange(e_min, e_max)
      .setThetaRange(0.0, th_max);
  double r_uc, e_uc;
  g.rateAndError(r_uc, e_uc, 20000000);
  const double area_cm2 = g.surfaceArea();
  const double uc_per_cm2 = r_uc / area_cm2;
  const double uc_err = e_uc / area_cm2;

  // ---- EcoMug --------------------------------------------------------------
  EcoMug eco;
  eco.SetUseSky();
  eco.SetSkySize({{2.0, 2.0}});                       // metres, area cancels
  eco.SetSkyCenterPosition({{0.0, 0.0, 0.0}});
  eco.SetDifferentialFlux(J);
  eco.SetMinimumMomentum(p_min);
  eco.SetMaximumMomentum(p_max);
  eco.SetMinimumTheta(0.0);
  eco.SetMaximumTheta(th_max);
  // EcoMug samples uniformly in momentum over a steeply falling spectrum, so
  // its estimator has much higher variance than its quoted error suggests at
  // modest statistics. Repeat over independent seeds and use the spread of the
  // means, which is the honest uncertainty on its answer.
  double r_eco = 0.0, e_eco = 0.0;
  const int kSeeds = 8;
  double vals[kSeeds], sum = 0.0;
  for (int i = 0; i < kSeeds; ++i) {
    eco.SetSeed(20260730 + 7919 * i);
    double r, e;
    eco.GetAverageGenRateAndError(r, e, 50000000);
    vals[i] = r; sum += r;
  }
  r_eco = sum / kSeeds;
  double var = 0.0;
  for (int i = 0; i < kSeeds; ++i) var += (vals[i] - r_eco) * (vals[i] - r_eco);
  e_eco = std::sqrt(var / (kSeeds - 1) / kSeeds);   // error on the mean

  std::printf("Same differential flux (Guan), p in [%.0f, %.0f] GeV/c, "
              "theta in [0, %.0f] deg\n\n", p_min, p_max, th_max_deg);
  std::printf("  UCMuGen  rate/area = %.6g +- %.1g  cm^-2 s^-1\n",
              uc_per_cm2, uc_err);
  std::printf("  EcoMug   rate/area = %.6g +- %.1g  cm^-2 s^-1\n",
              r_eco, e_eco);

  const double rel = std::fabs(uc_per_cm2 - r_eco) / r_eco;
  const double sigma = std::sqrt(uc_err * uc_err + e_eco * e_eco);
  std::printf("\n  difference = %.2e relative, %.2f combined sigma\n",
              rel, std::fabs(uc_per_cm2 - r_eco) / sigma);

  // For scale: the classic sea-level figure for a horizontal surface.
  std::printf("  for scale, UCMuGen = %.1f muons m^-2 s^-1 over this window\n",
              uc_per_cm2 * 1.0e4);

  // ---- The live-time gap ---------------------------------------------------
  std::printf("\nLive time for 1e6 muons through the same 2 m x 2 m plane:\n");
  std::printf("  UCMuGen  liveTime(1e6) = %.4g s\n", g.liveTime(1000000));
  std::printf("  EcoMug   GetEstimatedTime(1e6) = %.4g s",
              eco.GetEstimatedTime(1000000));
  std::printf("   <-- returns 0 for a user-supplied J\n");
  return 0;
}
