// A complete, runnable Geant4 application driven by UCMuGen.
//
// Everything specific to UCMuGen lives in PrimaryGeneratorAction.{hh,cc}, which
// is the only file worth copying into your own application. This file is the
// scaffolding that makes the example run on its own: a geometry, a physics
// list, a counter and a main().
//
//   ./ucmugen_example run.mac      batch, prints a rate and a live time
//   ./ucmugen_example tracks.mac   writes tracks for a figure
//   ./ucmugen_example              interactive with visualisation

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SteppingVerbose.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4Run.hh"
#include "FTFP_BERT.hh"

#include <atomic>

namespace {
std::atomic<long long> gMuonsInDetector{0};
}

// ---------------------------------------------------------------------------
// Geometry: air, with a 2 m x 2 m x 20 cm plastic-scintillator plate at the
// origin. Deliberately plain: this example is about the primaries.
// ---------------------------------------------------------------------------
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* air = nist->FindOrBuildMaterial("G4_AIR");
    auto* scint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    auto* worldBox = new G4Box("World", 6 * m, 6 * m, 6 * m);
    auto* worldLog = new G4LogicalVolume(worldBox, air, "World");
    auto* worldPhys = new G4PVPlacement(nullptr, {}, worldLog, "World",
                                        nullptr, false, 0);

    // Must match the BoxDetector in PrimaryGeneratorAction.cc.
    auto* detBox = new G4Box("Detector", 1 * m, 1 * m, 10 * cm);
    auto* detLog = new G4LogicalVolume(detBox, scint, "Detector");
    new G4PVPlacement(nullptr, {}, detLog, "Detector", worldLog, false, 0);

    return worldPhys;
  }
};

// ---------------------------------------------------------------------------
// Count muons entering the plate, once per track.
// ---------------------------------------------------------------------------
class SteppingAction : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* step) override {
    const auto* track = step->GetTrack();
    if (std::abs(track->GetDefinition()->GetPDGEncoding()) != 13) return;
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    if (!pre || !post || !post->GetTouchableHandle()->GetVolume()) return;
    const G4String from = pre->GetTouchableHandle()->GetVolume()->GetName();
    const G4String to = post->GetTouchableHandle()->GetVolume()->GetName();
    if (from != "Detector" && to == "Detector") ++gMuonsInDetector;
  }
};

// ---------------------------------------------------------------------------
// Turn the count into a rate, which is the point of the exercise.
// ---------------------------------------------------------------------------
class RunAction : public G4UserRunAction {
 public:
  explicit RunAction(const PrimaryGeneratorAction* gen) : fGen(gen) {}

  void BeginOfRunAction(const G4Run*) override { gMuonsInDetector = 0; }

  void EndOfRunAction(const G4Run* run) override {
    if (!IsMaster()) return;
    const long long n = run->GetNumberOfEvent();
    if (n == 0) return;
    const long long hits = gMuonsInDetector.load();
    const double live = fGen ? fGen->LiveTime(n) : 0.0;
    G4cout << "\n---------------- UCMuGen example ----------------\n"
           << "  generated muons        : " << n << "\n"
           << "  muons entering the plate: " << hits << "  ("
           << (100.0 * hits / n) << " %)\n";
    if (live > 0.0) {
      G4cout << "  rate into the detector : " << fGen->Rate() << " Hz\n"
             << "  live time of this run  : " << live << " s\n"
             << "  measured hit rate      : " << (hits / live) << " Hz\n";
    }
    G4cout << "-------------------------------------------------" << G4endl;
  }

 private:
  const PrimaryGeneratorAction* fGen;
};

class ActionInitialization : public G4VUserActionInitialization {
 public:
  void BuildForMaster() const override { SetUserAction(new RunAction(nullptr)); }
  void Build() const override {
    auto* gen = new PrimaryGeneratorAction;
    SetUserAction(gen);
    SetUserAction(new RunAction(gen));
    SetUserAction(new SteppingAction);
  }
};

int main(int argc, char** argv) {
  G4UIExecutive* ui = (argc == 1) ? new G4UIExecutive(argc, argv) : nullptr;

  // Serial by default: the example prints a single rate, and one generator per
  // worker with independent streams is a separate discussion. UCMuGen is fine
  // in MT, one Generator per thread.
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new FTFP_BERT);
  runManager->SetUserInitialization(new ActionInitialization);

  auto* vis = new G4VisExecutive;
  vis->Initialize();
  auto* uiManager = G4UImanager::GetUIpointer();

  if (ui) {
    uiManager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete vis;
  delete runManager;
  return 0;
}
