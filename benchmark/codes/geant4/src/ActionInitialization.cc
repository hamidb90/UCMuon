#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* det)
    : G4VUserActionInitialization(), fDetector(det) {}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    auto* run   = new RunAction();
    auto* event = new EventAction(run);
    auto* step  = new SteppingAction(event);
    auto* gen   = new PrimaryGeneratorAction();

    SetUserAction(gen);
    SetUserAction(run);
    SetUserAction(event);
    SetUserAction(step);
}