// The only file you need to copy into your own Geant4 application.
//
// The generator is held by value: it owns its RNG and its cached momentum CDF,
// and one instance per thread is exactly what Geant4's worker model wants.

#ifndef UCMUGEN_EXAMPLE_PRIMARYGENERATORACTION_HH
#define UCMUGEN_EXAMPLE_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "UCMuGen.h"

class G4Event;
class G4ParticleGun;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  /// Rate through the generation surface, s^-1. With a detector attached this
  /// is the rate *into the detector*, which is what a run normalises to.
  double Rate() const { return fRate; }

  /// Seconds of real time that `n` generated muons correspond to.
  double LiveTime(long long n) const { return (fRate > 0.0) ? n / fRate : 0.0; }

  /// Exposed so the run action can report on the sampling.
  const ucmugen::Generator& generator() const { return fGen; }

 private:
  G4ParticleGun* fGun = nullptr;
  ucmugen::Generator fGen;
  double fRate = 0.0;
};

#endif  // UCMUGEN_EXAMPLE_PRIMARYGENERATORACTION_HH
