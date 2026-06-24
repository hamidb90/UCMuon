#ifndef EventAction_hh
#define EventAction_hh

#include "G4UserEventAction.hh"
#include "G4Types.hh"
#include <vector>
#include <string>

class RunAction;

struct ScoreHit {
    int    depthIdx;
    double depthCm, mwe;
    int    pdg;
    double initKE_GeV, initTheta_deg, initPhi_deg;
    double exitKE_GeV;
    double exitPx, exitPy, exitPz;
    double exitX_cm, exitY_cm;
    double eLossIon, eLossBrem, eLossPair, eLossNucl, eLossTotal;
    double angleScat_deg, latDisp_cm;
};

struct SecHit {
    int    pdg;
    double depthCm;
    double ke_MeV;
    double px_MeV, py_MeV, pz_MeV;
    double x_cm, y_cm;
};

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* ra);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction  (const G4Event*) override;

    void AddHit      (const ScoreHit& h) { fHits.push_back(h); }
    void AddSecondary(const SecHit&   s) { fSecs.push_back(s); }

    void SetInitKE   (double ke)  { fInitKE_GeV    = ke / 1000.0; }
    void SetInitTheta(double th)  { fInitTheta_deg = th; }
    void SetInitPhi  (double ph)  { fInitPhi_deg   = ph; }
    void SetInitPDG  (int    pdg) { fInitPDG       = pdg; }

    double GetInitKE_GeV()    const { return fInitKE_GeV; }
    double GetInitTheta_deg() const { return fInitTheta_deg; }
    double GetInitPhi_deg()   const { return fInitPhi_deg; }
    int    GetInitPDG()       const { return fInitPDG; }

    void   SetStopInfo(double depth_cm) { fStopDepth_cm = depth_cm; fDidStop = true; }
    bool   DidStop()            const   { return fDidStop; }
    double GetStopDepth_cm()    const   { return fStopDepth_cm; }

private:
    RunAction*           fRunAction;
    std::vector<ScoreHit> fHits;
    std::vector<SecHit>   fSecs;

    double fInitKE_GeV    = 0.0;
    double fInitTheta_deg = 0.0;
    double fInitPhi_deg   = 0.0;
    int    fInitPDG       = 0;

    double fStopDepth_cm  = 0.0;
    bool   fDidStop       = false;
};

#endif