#ifndef RunAction_hh
#define RunAction_hh

#include "G4UserRunAction.hh"
#include <string>

// "muon" → muon CSV,  "sec" → secondaries CSV
void RunAction_WriteRow(const std::string& stream, const std::string& row);
// Track min/max initial KE (GeV) across all recorded muon events
void RunAction_UpdateEnergy(double initKE_GeV);

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction() override;
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction  (const G4Run*) override;
};

#endif