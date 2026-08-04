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
  ///
  /// Static and computed once. Every worker configures an identical generator,
  /// so the rate is a property of the setup rather than of a thread, and the
  /// Monte Carlo integral behind it is not worth repeating per worker. Being
  /// static is also what lets the master thread report it: in MT the master
  /// has no PrimaryGeneratorAction of its own.
  static double Rate();

  /// Seconds of real time that `n` generated muons correspond to.
  static double LiveTime(long long n) {
    const double r = Rate();
    return (r > 0.0) ? n / r : 0.0;
  }

  /// Exposed so the run action can report on the sampling.
  const ucmugen::Generator& generator() const { return fGen; }

  /// Build the configured generator. One place, so the generator the workers
  /// run and the generator the rate is computed from cannot drift apart.
  static ucmugen::Generator MakeGenerator();

  /// Install the PARMA flux, once, before any worker thread starts.
  ///
  /// PARMA is installed into a process-wide provider rather than into a
  /// generator, so this must NOT go in the constructor: Geant4 builds one
  /// PrimaryGeneratorAction per worker, and that would have every thread
  /// writing the same provider while the others read it. Call it from main()
  /// before the run manager exists. Without PARMA compiled in it does nothing,
  /// so it is always safe to call.
  static void InstallFlux();

 private:
  G4ParticleGun* fGun = nullptr;
  ucmugen::Generator fGen;
};

#endif  // UCMUGEN_EXAMPLE_PRIMARYGENERATORACTION_HH
