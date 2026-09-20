#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include <vector>

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume* Construct() override;
    const std::vector<G4LogicalVolume*>& GetScoringVolumes() const;
private:
    std::vector<G4LogicalVolume*> fScoringVolumes;
};
#endif
