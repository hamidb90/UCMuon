// Dump UCMuGen legacy-driver events in the exact column format of the Fortran
// generator's muons_surface.dat, so the Python validation harness can compare
// the two implementations without any format adaptation.
//
// Build:
//   c++ -std=c++17 -O2 -I../include dump_ucmugen.cc -o dump_ucmugen
//
// Usage (all arguments required, in this order):
//   dump_ucmugen emin emax spectrum angular theta_max_deg source \
//                radius_cm half_lx_cm half_ly_cm source_z_cm nmuons seed
//
// spectrum: 1 CosmoALEPH, 2 PowerLaw, 4 Guan, 5 Frosin, 6 GaisserBugaev,
//           7 ReynaBugaev, 8 Electron          (3 = PARMA is not supported)
// angular : 1 vertical, 2 cos^2, 3 uniform cone, 4 Guan self-consistent,
//           5 cos^3
// source  : 1 disk, 2 rectangle, 3 hemisphere

#include "UCMuGen.h"

#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv) {
  if (argc != 13) {
    std::fprintf(stderr,
                 "usage: %s emin emax spectrum angular theta_max_deg source "
                 "radius_cm half_lx_cm half_ly_cm source_z_cm nmuons seed\n",
                 argv[0]);
    return 2;
  }

  int a = 1;
  ucmugen::legacy::Config cfg;
  cfg.e_min          = std::atof(argv[a++]);
  cfg.e_max          = std::atof(argv[a++]);
  cfg.spectrum       = static_cast<ucmugen::Spectrum>(std::atoi(argv[a++]));
  cfg.angular        = static_cast<ucmugen::Angular>(std::atoi(argv[a++]));
  cfg.theta_max_deg  = std::atof(argv[a++]);
  cfg.source         = static_cast<ucmugen::legacy::Source>(std::atoi(argv[a++]));
  cfg.radius_cm      = std::atof(argv[a++]);
  cfg.half_lx_cm     = std::atof(argv[a++]);
  cfg.half_ly_cm     = std::atof(argv[a++]);
  cfg.source_z_cm    = std::atof(argv[a++]);
  const long n       = std::atol(argv[a++]);
  cfg.seed           = static_cast<std::int32_t>(std::atol(argv[a++]));

  if (cfg.spectrum == ucmugen::Spectrum::Parma) {
    std::fprintf(stderr, "PARMA is not available in the permissive build\n");
    return 3;
  }

  ucmugen::legacy::Generator gen(cfg);

  std::printf("# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV"
              "  theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask\n");
  for (long i = 1; i <= n; ++i) {
    const ucmugen::legacy::Muon m = gen.generate();
    // Matches the Fortran write format at ucmuon_gen_omp.f90:944.
    std::printf("%10ld %13.4f %13.4f %13.4f %13.6f %13.6f %13.6f %13.6f "
                "%13.9f %13.9f %13.6f %4d %2d %8d\n",
                i, m.x, m.y, m.z, m.p, m.px, m.py, m.pz,
                m.theta, m.phi, m.energy, m.charge, 1, 0);
  }
  return 0;
}
