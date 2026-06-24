#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(DetectorConstruction* detector);
    virtual ~ActionInitialization() = default;

    virtual void BuildForMaster() const override;
    virtual void Build() const override;

private:
    DetectorConstruction* fDetector;
};

#endif