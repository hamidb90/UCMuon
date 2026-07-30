// Minimal Geant4 primary generator using UCMuGen.
//
// This is the whole integration. Compare with the EcoMug path in
// benchmark/geant4_muon_rock_v5/src/PrimaryGeneratorAction.cc, which needs the
// caller to unpack (theta, phi, |p|), know that theta has already been flipped
// by pi, rebuild the momentum vector by hand, and convert total energy to
// kinetic energy.

#include "PrimaryGeneratorAction.hh"

#define UCMUGEN_WITH_GEANT4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include "UCMuGen.h"

using namespace ucmugen;

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  fGun = new G4ParticleGun(1);

  // A 2 m x 2 m horizontal sky plane, 5 m above the origin, Guan spectrum,
  // 1 GeV to 1 TeV, out to 70 degrees from vertical.
  fGen.setSpectrum(Spectrum::Guan)
      .setSurface(std::make_shared<Plane>(100.0 * units::cm,
                                          100.0 * units::cm,
                                          Vec3{0.0, 0.0, 500.0 * units::cm}))
      .setEnergyRange(1.0, 1000.0)
      .setThetaRange(0.0, 70.0 * kPi / 180.0)
      // Aim at the sensitive volume: a 1 m x 1 m x 20 cm box at the origin.
      // Only muons whose ray reaches it are generated, and rate() below
      // becomes the rate into the box rather than through the sky plane.
      // Comment this line out and the run still works, just far slower.
      .setDetector(BoxDetector::centred(Vec3{0.0, 0.0, 0.0},
                                        50.0 * units::cm,
                                        50.0 * units::cm,
                                        10.0 * units::cm))
      // Draw from Geant4's engine so the run is reproducible from the G4 seed
      // and stays consistent under /random/setSeeds.
      .setRng(std::make_unique<CallbackRng>([] { return G4UniformRand(); }));

  // Absolute normalisation: how much real time this run represents. Every
  // muography analysis needs this, and it is available for any spectrum and
  // surface combination, including a user-supplied flux.
  G4cout << fGen.provenance() << G4endl;
  G4cout << "  rate      : " << fGen.rate() << " Hz" << G4endl;
  G4cout << "  1e6 muons = " << fGen.liveTime(1000000) << " s of live time"
         << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fGun; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  FireG4(fGen, fGun, event);
}
