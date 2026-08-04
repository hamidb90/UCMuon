// The whole UCMuGen integration, and the only file you need to copy into your
// own application.
//
// Compare with the EcoMug path in ucmugen/comparison/cylinder_g4.cc, which
// needs the caller to unpack (theta, phi, |p|), know that theta has already
// been flipped by pi, rebuild the momentum vector by hand, and convert total
// energy to kinetic energy.

// Enables FireG4 below. The header is arranged so that defining this after a
// first plain include is not silently fatal, but defining it first is tidier.
#define UCMUGEN_WITH_GEANT4

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#ifdef UCMUGEN_EXAMPLE_PARMA
#include "UCMuGen_PARMA.h"
#endif

#include <memory>

using namespace ucmugen;

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  fGun = new G4ParticleGun(1);

  // Everything below is in UCMuGen's units: centimetres and GeV.
  // The sensitive volume is a 2 m x 2 m x 20 cm plate at the origin, matching
  // the geometry in ucmugen_example.cc.
  auto detector = BoxDetector::centred(Vec3{0.0, 0.0, 0.0},
                                       100.0, 100.0, 10.0);

  fGen.setSpectrum(Spectrum::Guan)
      // A 6 m x 6 m horizontal sky plane, 3 m above the detector.
      .setSurface(std::make_shared<Plane>(300.0, 300.0, Vec3{0.0, 0.0, 300.0}))
      .setEnergyRange(1.0, 1000.0)            // total energy, GeV
      .setThetaRange(0.0, 70.0 * kPi / 180.0)
      // Aim at the detector: only muons whose ray reaches it are generated,
      // and rate() below becomes the rate into it rather than through the sky
      // plane. Comment this out and the example still works, just far slower.
      .setDetector(detector)
      // Draw from Geant4's engine, so a run is reproducible from the Geant4
      // seed and stays consistent under /random/setSeeds.
      .setRng(std::make_unique<CallbackRng>([] { return G4UniformRand(); }));

#ifdef UCMUGEN_EXAMPLE_PARMA
  // Location- and date-aware flux. Sea level here; pass a real altitude to
  // depth_from_altitude() for a mountain or volcano site, which is the whole
  // reason this option exists.
  parma::Site site;
  site.cutoff_GV = 3.25;                                   // geomagnetic cutoff
  site.depth_gcm2 = parma::depth_from_altitude(0.0, 50.67);
  parma::install(site);
  fGen.setSpectrum(Spectrum::Parma);
  G4cout << "[UCMuGen] PARMA flux: cutoff " << site.cutoff_GV
         << " GV, atmospheric depth " << site.depth_gcm2 << " g/cm2" << G4endl;
#endif

  // Absolute normalisation. This is a Monte Carlo integral, so compute it once
  // rather than per event. Every muography analysis needs it, and it is
  // available for any spectrum and surface combination.
  fRate = fGen.rate();

  G4cout << "[UCMuGen] " << fGen.provenance() << G4endl;
  G4cout << "[UCMuGen] rate into the detector : " << fRate << " Hz" << G4endl;
  G4cout << "[UCMuGen] 1e6 muons correspond to: " << LiveTime(1000000)
         << " s of live time" << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fGun; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // This one line is the integration. It converts cm to mm, GeV to MeV, total
  // energy to kinetic energy, and the charge to a PDG code whose sign is the
  // opposite of the naive guess.
  FireG4(fGen, fGun, event);
}
