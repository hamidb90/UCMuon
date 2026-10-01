// intensity_dump.cc -- UCMuGen side of the cross-implementation flux test.
//
// Prints, for every spectrum with a built-in absolute normalisation, the
// differential intensity I(p, cos theta) on a fixed grid, the vertical
// integrated flux the momentum CDF reports, and the rate through a horizontal
// surface both by quadrature (the same algorithm as intensity_dump.f90 and
// test_flux_consistency.py) and from Generator::rateAndError, the Monte Carlo
// estimator users actually call. test_flux_consistency.py compares the lot.
//
//   c++ -std=c++17 -O2 -I ucmugen/include tests/flux/intensity_dump.cc
#include "UCMuGen.h"

#include <cstdio>

using namespace ucmugen;

namespace {

const int kModes[] = {1, 4, 5, 6, 7};
const double kP[] = {1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 2000.0};
const double kCos[] = {1.0, 0.8660254037844387, 0.5, 0.25881904510252074, 0.1};
const double kEmin[] = {1.0, 10.0, 100.0};
const double kEmax = 1500.0;
const double kThetaMaxDeg = 85.0;

double pOfE(double e) { return std::sqrt(e * e - kMuonMass * kMuonMass); }

// R = 2 pi * int_{cos theta_max}^1 c dc int_{p_min}^{p_max} I(p, c) dp, per cm^2.
// Trapezoid on a 2000-point log grid in p, midpoint rule on 400 cells in c.
double horizontalRate(Spectrum s, double p_min, double p_max) {
  const int np = 2000, nc = 400;
  const double c_lo = std::cos(kThetaMaxDeg * kPi / 180.0);
  const double lr = std::log(p_max / p_min);
  double total = 0.0;
  for (int k = 0; k < nc; ++k) {
    const double c = c_lo + (k + 0.5) * (1.0 - c_lo) / nc;
    double inner = 0.0, prev_p = p_min, prev_f = flux::intensity(s, p_min, c);
    for (int j = 1; j < np; ++j) {
      const double p = p_min * std::exp(lr * j / (np - 1));
      const double f = flux::intensity(s, p, c);
      inner += 0.5 * (f + prev_f) * (p - prev_p);
      prev_p = p;
      prev_f = f;
    }
    total += c * inner * (1.0 - c_lo) / nc;
  }
  return 2.0 * kPi * total;
}

}  // namespace

int main() {
  for (int m : kModes)
    for (double p : kP)
      for (double c : kCos)
        std::printf("I %d %.17g %.17g %.17g\n", m, p, c,
                    flux::intensity(Spectrum(m), p, c));

  for (int m : kModes)
    for (double e : kEmin) {
      const double p_min = pOfE(e), p_max = pOfE(kEmax);
      MomentumCdf cdf(p_min, p_max, Spectrum(m));
      std::printf("V %d %.17g %.17g\n", m, e, cdf.integrated_flux());
      std::printf("Q %d %.17g %.17g\n", m, e,
                  horizontalRate(Spectrum(m), p_min, p_max));

      Generator g;
      g.setSpectrum(Spectrum(m)).setPrintWarnings(false)
       .setSurface(std::make_shared<Disk>(1.0 / std::sqrt(kPi), Vec3{0, 0, 0}))
       .setEnergyRange(e, kEmax)
       .setThetaRange(0.0, kThetaMaxDeg * kPi / 180.0);
      double r, err;
      g.rateAndError(r, err, 400000);            // 1 cm^2 disk: rate per cm^2
      std::printf("M %d %.17g %.17g %.17g\n", m, e, r, err);
      std::printf("W %d %.17g %zu\n", m, e, g.warnings().size());
    }
  return 0;
}
