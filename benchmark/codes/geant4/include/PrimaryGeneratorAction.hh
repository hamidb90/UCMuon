#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include <vector>

// EcoMug.h is NOT included here — only in PrimaryGeneratorAction.cc.
// fEcoMugPtr is void* so no EcoMug type leaks into this header.

class G4Event;

struct MuonRecord {
    int    pdg;
    double x_mm, y_mm, z_mm;
    double px_MeV, py_MeV, pz_MeV;
    double ke_MeV;
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event* event) override;

    double GetInitKE_MeV() const { return fInitKE_MeV; }
    double GetInitTheta()  const { return fInitTheta; }
    double GetInitPhi()    const { return fInitPhi; }
    int    GetInitPDG()    const { return fInitPDG; }

private:
    void   LoadSourceFile();
    void   FireFromFile(G4Event*);
    void   FirePowerLaw(G4Event*);
    void   FireEcoMug(G4Event*);
    double SamplePowerLaw();
    void   SetGunFromRecord(double x, double y, double z,
                            double px, double py, double pz,
                            double ke, int pdg);
    void   RecordInitialState(double ke, double px, double py, double pz, int pdg);

    G4ParticleGun*          fGun      = nullptr;
    std::vector<MuonRecord> fRecords;
    int                     fIndex    = 0;
    void*                   fEcoMugPtr = nullptr;   // EcoMug* cast in .cc

    double fInitKE_MeV = 0.0;
    double fInitTheta  = 0.0;
    double fInitPhi    = 0.0;
    int    fInitPDG    = 13;
};
#endif
