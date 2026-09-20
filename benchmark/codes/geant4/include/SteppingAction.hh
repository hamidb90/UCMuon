#ifndef SteppingAction_hh
#define SteppingAction_hh

#include "G4UserSteppingAction.hh"
#include <set>
#include <utility>

class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
    explicit SteppingAction(EventAction* ea);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

private:
    EventAction* fEventAction;

    // Prevents recording the same (eventID, depthCm) plane twice
    std::set<std::pair<int,int>> fRecordedPlanes;

    // Accumulated per-process energy loss for the current primary muon
    double fELossIon  = 0.0;
    double fELossBrem = 0.0;
    double fELossPair = 0.0;
    double fELossNucl = 0.0;
};

#endif