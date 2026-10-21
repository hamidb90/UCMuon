// parma_altitude.cc -- PARMA muon intensity at the MURAVES site (40.8271 N,
// 14.4006 E) at sea level and at the detector's 608 m, for inversion_tests.py.
// Needs PARMA (python3 tools/install_parma.py): UCMuGen_PARMA.h and data/EXPACS.
//
// Prints "alt_km zenith_deg p_GeV I" with I in cm^-2 s^-1 sr^-1 (GeV/c)^-1.
#include "UCMuGen.h"
#include "UCMuGen_PARMA.h"

#include <cmath>
#include <cstdio>

int main(int argc, char** argv) {
  ucmugen::parma::set_data_dir(argc > 1 ? argv[1] : "data/EXPACS/parma");
  const double lat = 40.8271, lon = 14.4006;
  for (double alt : {0.0, 0.608}) {
    ucmugen::parma::Site s;
    s.cutoff_GV = ucmugen::parma::cutoff_from_location(lat, lon);
    s.depth_gcm2 = ucmugen::parma::depth_from_altitude(alt, lat);
    s.w_index = 58.5;                         // mean W 1951-2026
    for (double th : {75.0, 80.0, 85.0}) {
      const double c = std::cos(th * M_PI / 180.0);
      for (int i = 0; i <= 400; ++i) {
        const double p = 0.5 * std::pow(10.0, i * std::log10(4.0e4) / 400.0);
        std::printf("%.3f %.1f %.8g %.10g\n", alt, th, p,
                    ucmugen::parma::intensity(s, p, c));
      }
    }
  }
  return 0;
}
