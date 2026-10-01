// ucmugen_sample.cc -- events and rate from ucmugen::Generator, the reference
// sampler of J(p, theta) * max(0, -n.d), for the joint-sampling and live-time
// tests. Geometry is given in the Fortran generator's conventions (cm).
//
//   ucmugen_sample SPEC EMIN EMAX THMAX_DEG SURFACE ... [det ...] N SEED
//     SURFACE = disk   R CX CY CZ NX NY NZ
//             | plane  HX HY CX CY CZ NX NY NZ
//             | hsphere R CZ
//     det     = cyl AX AY AZ BX BY BZ R MARGIN   (Fortran's inflated cylinder)
//             | box XMIN XMAX YMIN YMAX ZMIN ZMAX MARGIN
//             | none
//     N = 0 prints only the rate.
//
// Output: "# rate R ERR" (s^-1; into the detector when one is given), then
// one line per event: p cx cy cz x y z.
//
// PARMA (SPEC 3) needs -DUCMUGEN_TEST_PARMA and four more leading arguments:
// W CUTOFF_GV DEPTH_GCM2 LOCAL_G, before SPEC.
#include "UCMuGen.h"
#ifdef UCMUGEN_TEST_PARMA
#include "UCMuGen_PARMA.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace ucmugen;

int main(int argc, char** argv) {
  int a = 1;
  auto next = [&]() -> double {
    if (a >= argc) { std::fprintf(stderr, "missing argument\n"); std::exit(2); }
    return std::atof(argv[a++]);
  };
#ifdef UCMUGEN_TEST_PARMA
  parma::Site site;
  site.w_index = next();
  site.cutoff_GV = next();
  site.depth_gcm2 = next();
  site.local_g = next();
  parma::install(site);
#endif
  const int spec = int(next());
  const double emin = next(), emax = next(), thmax = next() * kPi / 180.0;

  std::shared_ptr<Surface> surf;
  const char* kind = argv[a++];
  if (!std::strcmp(kind, "disk")) {
    const double r = next(), cx = next(), cy = next(), cz = next();
    const double nx = next(), ny = next(), nz = next();
    surf = std::make_shared<Disk>(r, Vec3{cx, cy, cz}, Vec3{nx, ny, nz});
  } else if (!std::strcmp(kind, "plane")) {
    const double hx = next(), hy = next(), cx = next(), cy = next(), cz = next();
    const double nx = next(), ny = next(), nz = next();
    surf = std::make_shared<Plane>(hx, hy, Vec3{cx, cy, cz}, Vec3{nx, ny, nz});
  } else {
    const double r = next(), cz = next();
    surf = std::make_shared<HSphere>(r, Vec3{0.0, 0.0, cz});
  }

  std::shared_ptr<Detector> det;
  const char* dk = argv[a++];
  if (!std::strcmp(dk, "cyl")) {
    Vec3 A{next(), next(), next()}, B{next(), next(), next()};
    const double r = next(), m = next();
    const Vec3 u = (B - A).unit();
    det = std::make_shared<CylinderDetector>(A - u * m, B + u * m, r + m);
  } else if (!std::strcmp(dk, "box")) {
    const double x0 = next(), x1 = next(), y0 = next(), y1 = next();
    const double z0 = next(), z1 = next(), m = next();
    det = std::make_shared<BoxDetector>(Vec3{x0 - m, y0 - m, z0 - m},
                                        Vec3{x1 + m, y1 + m, z1 + m});
  }
  const long long n = (long long)next();
  const int seed = int(next());

  Generator g;
  g.setSpectrum(Spectrum(spec)).setPrintWarnings(false).setSurface(surf)
   .setEnergyRange(emin, emax).setThetaRange(0.0, thmax).setSeed(seed);
  if (det) g.setDetector(det);
  double r, err;
  g.rateAndError(r, err, 2000000);
  std::printf("# rate %.10g %.10g\n", r, err);
  for (long long i = 0; i < n; ++i) {
    const Muon m = g.generate();
    std::printf("%.9g %.9g %.9g %.9g %.6f %.6f %.6f\n", m.momentum,
                m.direction.x, m.direction.y, m.direction.z,
                m.position.x, m.position.y, m.position.z);
  }
  return 0;
}
