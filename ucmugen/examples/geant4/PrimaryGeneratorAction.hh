// Header for the minimal UCMuGen primary generator.
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

  /// Exposed so the run action can report live time for the run.
  const ucmugen::Generator& generator() const { return fGen; }

 private:
  G4ParticleGun* fGun = nullptr;
  ucmugen::Generator fGen;
};

#endif  // UCMUGEN_EXAMPLE_PRIMARYGENERATORACTION_HH
