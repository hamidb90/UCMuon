#include "EventAction.hh"
#include "RunAction.hh"
#include "SimConfig.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>
#include <iomanip>

EventAction::EventAction(RunAction* ra)
    : G4UserEventAction(), fRunAction(ra) {}
EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    fHits.clear();
    fSecs.clear();
    fInitKE_GeV = fInitTheta_deg = fInitPhi_deg = 0.0;
    fInitPDG       = 0;
    fStopDepth_cm  = 0.0;
    fDidStop       = false;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
    int evtID = evt->GetEventID();

    if (!fHits.empty())
        RunAction_UpdateEnergy(fInitKE_GeV);

    for (const auto& h : fHits) {
        std::ostringstream row;
        row << std::fixed << std::setprecision(6)
            << evtID             << "," << h.depthCm       << "," << h.mwe          << ","
            << h.pdg             << "," << fInitKE_GeV     << "," << fInitTheta_deg << ","
            << fInitPhi_deg      << "," << h.exitKE_GeV    << "," << h.exitPx       << ","
            << h.exitPy          << "," << h.exitPz        << "," << h.exitX_cm     << ","
            << h.exitY_cm        << "," << h.eLossIon      << "," << h.eLossBrem    << ","
            << h.eLossPair       << "," << h.eLossNucl     << "," << h.eLossTotal   << ","
            << h.angleScat_deg   << "," << h.latDisp_cm    << "\n";
        RunAction_WriteRow("muon", row.str());
    }

    for (const auto& s : fSecs) {
        std::ostringstream row;
        row << std::fixed << std::setprecision(6)
            << evtID      << "," << s.pdg    << "," << s.depthCm << ","
            << s.ke_MeV   << "," << s.px_MeV << "," << s.py_MeV  << ","
            << s.pz_MeV   << "," << s.x_cm   << "," << s.y_cm    << "\n";
        RunAction_WriteRow("sec", row.str());
    }

    if (fDidStop && fInitKE_GeV > 0.0) {
        std::ostringstream row;
        row << std::fixed << std::setprecision(6)
            << evtID          << "," << fInitPDG       << ","
            << fInitKE_GeV    << "," << fInitTheta_deg << ","
            << fInitPhi_deg   << "," << fStopDepth_cm  << "\n";
        RunAction_WriteRow("stop", row.str());
    }
}