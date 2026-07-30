// End-to-end check of the Geant4 hand-off.
//
// The physics is validated elsewhere. What is checked here is the boundary:
// FireG4 converts cm to mm, GeV to MeV, total energy to kinetic energy, and a
// charge to a PDG code whose sign is the opposite of the naive guess. Those
// four conversions are exactly the ones hand-rolled integrations get wrong, and
// none of them is exercised by a test that does not link Geant4.
//
// The trick that makes it checkable: run two generators with the same seed, one
// through FireG4 and one calling generate() directly. They produce the same
// muon, so the primary vertex Geant4 ends up holding can be compared against
// the muon that produced it.
//
// A run manager, a world box and a two-particle physics list are set up because
// Geant4 refuses access to the particle table until a physics list has been
// registered (G4ParticleTable::CheckReadiness, fatal since version 8.0). That
// is the minimum; there is no tracking here.
//
// Local variables avoid the names `g` and `m`: CLHEP defines them as gram and
// metre, and -Wshadow (which geant4-config turns on) flags every one.
//
// Build (Geant4 must be on the path):
//   c++ -std=c++17 -O2 $(geant4-config --cflags) -I../include \
//       test_geant4.cc $(geant4-config --libs) -o test_geant4

#define UCMUGEN_WITH_GEANT4

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PVPlacement.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "Randomize.hh"

#include "UCMuGen.h"

#include <cmath>
#include <cstdio>
#include <memory>

using namespace ucmugen;

static int g_failures = 0;

static void fail(const char* what, double got, double expect) {
  ++g_failures;
  std::printf("  [FAIL] %-42s got=%.10g  expect=%.10g\n", what, got, expect);
}

static void nearEq(const char* what, double got, double expect, double atol) {
  if (std::fabs(got - expect) > atol) fail(what, got, expect);
}

/// Just enough geometry for the run manager to accept an initialisation.
class World : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* vac = new G4Material("vac", 1.0, 1.008 * g / mole, 1.0e-25 * g / cm3);
    auto* box = new G4Box("world", 10.0 * m, 10.0 * m, 10.0 * m);
    auto* log = new G4LogicalVolume(box, vac, "world");
    return new G4PVPlacement(nullptr, {}, log, "world", nullptr, false, 0);
  }
};

/// Muons only. The gun looks them up by PDG code and nothing is tracked.
class TwoParticles : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override {
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
  }
  void ConstructProcess() override { AddTransportation(); }
};

static Generator build() {
  Generator gen;
  gen.setSpectrum(Spectrum::Guan)
      .setSurface(std::make_shared<Plane>(100.0, 100.0, Vec3{0.0, 0.0, 500.0}))
      .setEnergyRange(1.0, 1000.0)
      .setThetaRange(0.0, 70.0 * kPi / 180.0)
      .setDetector(BoxDetector::centred(Vec3{0.0, 0.0, 0.0}, 50.0, 50.0, 10.0))
      .setSeed(424242);
  return gen;
}

int main() {
  auto* run = new G4RunManager;
  run->SetUserInitialization(new World);
  run->SetUserInitialization(new TwoParticles);
  run->Initialize();

  Generator fired = build();     // drives FireG4
  Generator shadow = build();    // same seed, so the same muon sequence
  auto* gun = new G4ParticleGun(1);
  auto detector = BoxDetector::centred(Vec3{0.0, 0.0, 0.0}, 50.0, 50.0, 10.0);

  constexpr int kN = 5000;
  int mu_plus = 0;
  std::printf("1. unit and convention conversions across FireG4 (%d events)\n",
              kN);

  for (int i = 0; i < kN; ++i) {
    G4Event event;
    FireG4(fired, gun, &event);
    const Muon mu = shadow.generate();

    const G4PrimaryVertex* v = event.GetPrimaryVertex();
    if (!v) { fail("primary vertex exists", 0, 1); break; }
    const G4PrimaryParticle* p = v->GetPrimary();
    if (!p) { fail("primary particle exists", 0, 1); break; }

    // Position: UCMuGen works in cm, Geant4 stores mm.
    nearEq("position x [mm] = 10 * cm", v->GetX0() / mm, mu.position.x * 10.0,
          1e-6);
    nearEq("position y [mm] = 10 * cm", v->GetY0() / mm, mu.position.y * 10.0,
          1e-6);
    nearEq("position z [mm] = 10 * cm", v->GetZ0() / mm, mu.position.z * 10.0,
          1e-6);

    // Energy: the gun is set with *kinetic* energy in GeV, stored in MeV.
    // Handing it the total energy instead is the classic slip; at 1 GeV that
    // is an 11% error in kinetic energy and it would still look plausible.
    nearEq("kinetic [MeV] = 1000 * GeV", p->GetKineticEnergy() / MeV,
          mu.kinetic * 1000.0, 1e-3);

    // PDG: mu+ is -13. Getting this backwards silently swaps the charge ratio.
    if (p->GetPDGcode() != mu.pdg) fail("pdg code", p->GetPDGcode(), mu.pdg);
    if (mu.pdg == -13) ++mu_plus;

    const G4ThreeVector d = p->GetMomentumDirection();
    nearEq("direction x", d.x(), mu.direction.x, 1e-9);
    nearEq("direction y", d.y(), mu.direction.y, 1e-9);
    nearEq("direction z", d.z(), mu.direction.z, 1e-9);
    nearEq("direction is a unit vector", d.mag(), 1.0, 1e-9);

    // The muon Geant4 is about to track must still be aimed at the detector:
    // directed sampling has to survive the hand-off, in Geant4's own units.
    const Vec3 pos_from_g4{v->GetX0() / cm, v->GetY0() / cm, v->GetZ0() / cm};
    const Vec3 dir_from_g4{d.x(), d.y(), d.z()};
    if (!detector->intersects(pos_from_g4, dir_from_g4))
      fail("G4 primary still aimed at the detector", 0, 1);

    if (g_failures) break;      // one report is enough
  }
  if (!g_failures)
    std::printf("  [ok  ] all conversions exact over %d events\n", kN);

  std::printf("\n2. charge ratio survives the hand-off\n");
  {
    const double ratio = double(mu_plus) / double(kN - mu_plus);
    const bool ok = ratio > 1.1 && ratio < 1.45;
    if (!ok) ++g_failures;
    std::printf("  [%s] mu+/mu- = %.3f (expect ~1.2 to 1.4)\n",
                ok ? "ok  " : "FAIL", ratio);
  }

  std::printf("\n3. driving the generator from Geant4's own engine\n");
  {
    // The documented multithreading pattern: one stream per worker, owned by
    // Geant4, so /random/setSeeds reproduces a run.
    auto make = [] {
      Generator gen;
      gen.setSpectrum(Spectrum::Guan)
          .setSurface(std::make_shared<Disk>(200.0, Vec3{0.0, 0.0, 500.0}))
          .setEnergyRange(1.0, 1000.0)
          .setRng(std::make_unique<CallbackRng>([] { return G4UniformRand(); }));
      return gen;
    };
    G4Random::setTheSeed(987654321);
    Generator a = make();
    double first[8];
    for (int i = 0; i < 8; ++i) first[i] = a.generate().momentum;

    G4Random::setTheSeed(987654321);
    Generator b = make();
    bool same = true;
    for (int i = 0; i < 8; ++i)
      if (b.generate().momentum != first[i]) same = false;

    if (!same) ++g_failures;
    std::printf("  [%s] same G4 seed reproduces the muon sequence\n",
                same ? "ok  " : "FAIL");
  }

  delete gun;
  std::printf("\n%s (%d failure%s)\n",
              g_failures ? "FAILURES PRESENT" : "ALL PASSED", g_failures,
              g_failures == 1 ? "" : "s");
  return g_failures ? 1 : 0;
}
