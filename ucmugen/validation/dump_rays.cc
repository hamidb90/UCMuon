// Dump generated primary rays for the blind/directed figure.
//
// Writes one CSV row per generated muon: which sampling mode produced it, where
// it started on the generation surface, where it is going, and whether its ray
// reaches the detector. That is exactly the generator's output, so the figure
// drawn from it shows what the sampler does rather than a re-enactment.
//
// Build:
//   c++ -std=c++17 -O2 -I../include dump_rays.cc -o dump_rays
//   ./dump_rays > rays.csv

#include "UCMuGen.h"

#include <cstdio>
#include <memory>

using namespace ucmugen;

// A muography-scale geometry rather than the Geant4 example's: a 1 x 1 m plate
// under a 10 x 10 m sky plane 5 m above. The example uses a large plate close to
// the surface, where blind sampling is only ~10x worse; the contrast the figure
// is about grows as the detector shrinks and recedes, which is the real case.
static const double kDetHX = 50.0, kDetHY = 50.0, kDetHZ = 10.0;
static const double kSkyHalf = 500.0, kSkyZ = 500.0;

int main() {
  auto detector = BoxDetector::centred(Vec3{0, 0, 0}, kDetHX, kDetHY, kDetHZ);

  std::printf("mode,x,y,z,dx,dy,dz,hit\n");
  for (int directed = 0; directed <= 1; ++directed) {
    Generator g;
    g.setSpectrum(Spectrum::Guan)
        .setSurface(std::make_shared<Plane>(kSkyHalf, kSkyHalf,
                                            Vec3{0, 0, kSkyZ}))
        .setEnergyRange(1.0, 1000.0)
        .setThetaRange(0.0, 70.0 * kPi / 180.0)
        .setSeed(20260730);
    if (directed) g.setDetector(detector);

    // Equal counts per panel so the two are visually comparable; the point of
    // the figure is where the rays go, not how many were drawn.
    const int kN = 220;
    for (int i = 0; i < kN; ++i) {
      const Muon m = g.generate();
      const bool hit = detector->intersects(m.position, m.direction);
      std::printf("%s,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%d\n",
                  directed ? "directed" : "blind",
                  m.position.x, m.position.y, m.position.z,
                  m.direction.x, m.direction.y, m.direction.z, hit ? 1 : 0);
    }
    // The acceptance is the quantity the figure is really about, so record it
    // for the caption rather than leaving it to be quoted from elsewhere.
    std::fprintf(stderr, "%s: acceptance %.5f%%\n",
                 directed ? "directed" : "blind", 100.0 * g.acceptance());
  }
  return 0;
}
