// parma_dump.cc -- PARMA (spectrum 3) values for test_flux_reference.py.
//
// Separate from intensity_dump.cc because UCMuGen_PARMA.h is JAEA's and
// non-commercial-only: the reference test skips PARMA when it is absent.
// Site: ucmugen::parma::Site defaults (sea level, 3 GV cutoff, W = 0).
//
// Prints "I p cos I" on the grid the reference test asks for on stdin.
#include "UCMuGen.h"
#include "UCMuGen_PARMA.h"

#include <cstdio>

int main() {
  ucmugen::parma::install(ucmugen::parma::Site{});
  double p, c;
  while (std::scanf("%lf %lf", &p, &c) == 2)
    std::printf("%.17g %.17g %.17g\n", p, c,
                ucmugen::flux::intensity(ucmugen::Spectrum::Parma, p, c));
  return 0;
}
