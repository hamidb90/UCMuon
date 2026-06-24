#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>
#include <string>

SteppingAction::SteppingAction(EventAction* ea)
    : G4UserSteppingAction(), fEventAction(ea) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    int  pdg    = track->GetDefinition()->GetPDGEncoding();
    bool isMuon = (std::abs(pdg) == 13);
    bool isPrimary = (track->GetTrackID() == 1);

    // ── 1. Cache primary kinematics on first step ─────────────────────────────
    if (isMuon && isPrimary && track->GetCurrentStepNumber() == 1) {
        double ke_MeV = track->GetKineticEnergy();
        G4ThreeVector mom = track->GetMomentum();
        // theta = angle from downward axis (-Z)
        double px = mom.x(), py = mom.y(), pz = mom.z();
        double p  = mom.mag();
        double theta_deg = (p > 0) ? std::acos(-pz / p) * 180.0 / M_PI : 0.0;
        double phi_deg   = std::atan2(py, px) * 180.0 / M_PI;

        fEventAction->SetInitKE(ke_MeV);          // MeV → stored as GeV inside
        fEventAction->SetInitTheta(theta_deg);
        fEventAction->SetInitPhi(phi_deg);
        fEventAction->SetInitPDG(pdg);

        // Reset accumulated energy loss for this new primary
        fELossIon = fELossBrem = fELossPair = fELossNucl = 0.0;
        fRecordedPlanes.clear();
    }

    // ── 2. Accumulate per-process energy loss for primary muon ────────────────
    if (isMuon && isPrimary) {
        double deTotal = step->GetTotalEnergyDeposit() / GeV;
        if (deTotal > 0.0) {
            auto* postProc = step->GetPostStepPoint()->GetProcessDefinedStep();
            if (postProc) {
                const std::string& proc = postProc->GetProcessName();
                if      (proc.find("Ioni")    != std::string::npos) fELossIon  += deTotal;
                else if (proc.find("Brem")    != std::string::npos) fELossBrem += deTotal;
                else if (proc.find("PairProd")!= std::string::npos) fELossPair += deTotal;
                else if (proc.find("Nuclear") != std::string::npos ||
                         proc.find("nuclear") != std::string::npos) fELossNucl += deTotal;
                else                                                 fELossIon  += deTotal;
            }
        }
    }

    // ── 2b. Detect primary muon stopping inside rock ─────────────────────────
    // Only count range-out / decay inside the rock or scoring planes.
    // fStopAndKill also fires when Geant4 kills a track at the world boundary,
    // so we guard against recording side-exit muons as "stopped".
    if (isMuon && isPrimary && track->GetTrackStatus() == fStopAndKill) {
        auto* preVol2 = step->GetPreStepPoint()->GetPhysicalVolume();
        if (preVol2) {
            G4String vn = preVol2->GetLogicalVolume()->GetName();
            bool inRock = (vn == "Rock") ||
                          (vn.size() >= 6 && vn.substr(0, 6) == "Score_");
            if (inRock) {
                G4ThreeVector stopPos = track->GetPosition();
                double stopDepth_cm = -stopPos.z() / 10.0;  // G4 mm → cm
                if (stopDepth_cm >= 0.0)
                    fEventAction->SetStopInfo(stopDepth_cm);
            }
        }
    }

    // ── 3. Volume boundary check ──────────────────────────────────────────────
    auto* preVol  = step->GetPreStepPoint() ->GetPhysicalVolume();
    auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    if (!preVol || !postVol) return;

    G4String postVolName = postVol->GetLogicalVolume()->GetName();
    G4String preVolName  = preVol ->GetLogicalVolume()->GetName();

    bool postIsScore = (postVolName.size() >= 6 && postVolName.substr(0,6) == "Score_");
    bool preIsScore  = (preVolName .size() >= 6 && preVolName .substr(0,6) == "Score_");
    bool entering    = postIsScore && !preIsScore;
    if (!entering) return;

    // Parse depth: "Score_50cm" → 50
    std::string vn = std::string(postVolName);
    int depthCm = 0;
    try {
        size_t end = vn.rfind("cm");
        if (end == std::string::npos) return;
        depthCm = std::stoi(vn.substr(6, end - 6));
    } catch (...) { return; }

    double mwe  = depthCm * 2.65 / 100.0;
    int    evID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // ── 3a. Primary muon scoring ──────────────────────────────────────────────
    if (isMuon && isPrimary) {
        auto key = std::make_pair(evID, depthCm);
        if (fRecordedPlanes.count(key)) return;
        fRecordedPlanes.insert(key);

        G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
        G4ThreeVector mom = step->GetPostStepPoint()->GetMomentum();
        double exitKE_GeV = step->GetPostStepPoint()->GetKineticEnergy() / GeV;
        double px = mom.x(), py = mom.y(), pz = mom.z();
        double p  = mom.mag();
        double angleScat_deg = (p > 0)
            ? std::acos(std::max(-1.0, std::min(1.0, -pz/p))) * 180.0 / M_PI : 0.0;
        double latDisp_cm = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y()) / 10.0;

        ScoreHit h;
        h.depthCm       = depthCm;
        h.mwe           = mwe;
        h.pdg           = pdg;
        h.initKE_GeV    = fEventAction->GetInitKE_GeV();
        h.initTheta_deg = fEventAction->GetInitTheta_deg();
        h.initPhi_deg   = fEventAction->GetInitPhi_deg();
        h.exitKE_GeV    = exitKE_GeV;
        h.exitPx        = px / GeV;
        h.exitPy        = py / GeV;
        h.exitPz        = pz / GeV;
        h.exitX_cm      = pos.x() / 10.0;
        h.exitY_cm      = pos.y() / 10.0;
        h.eLossIon      = fELossIon;
        h.eLossBrem     = fELossBrem;
        h.eLossPair     = fELossPair;
        h.eLossNucl     = fELossNucl;
        h.eLossTotal    = fELossIon + fELossBrem + fELossPair + fELossNucl;
        h.angleScat_deg = angleScat_deg;
        h.latDisp_cm    = latDisp_cm;
        fEventAction->AddHit(h);
    }

    // ── 3b. Secondary particles (from primary muon parent) ────────────────────
    if (!isMuon && track->GetParentID() == 1) {
        double ke_MeV = step->GetPostStepPoint()->GetKineticEnergy() / MeV;
        if (ke_MeV < 0.1) return;
        G4ThreeVector sp = step->GetPostStepPoint()->GetPosition();
        G4ThreeVector sm = step->GetPostStepPoint()->GetMomentum();
        SecHit s;
        s.pdg     = pdg;
        s.depthCm = depthCm;
        s.ke_MeV  = ke_MeV;
        s.px_MeV  = sm.x() / MeV;
        s.py_MeV  = sm.y() / MeV;
        s.pz_MeV  = sm.z() / MeV;
        s.x_cm    = sp.x() / 10.0;
        s.y_cm    = sp.y() / 10.0;
        fEventAction->AddSecondary(s);
    }
}