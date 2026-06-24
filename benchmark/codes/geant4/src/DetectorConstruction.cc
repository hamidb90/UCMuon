#include "DetectorConstruction.hh"
#include "SimConfig.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

DetectorConstruction::DetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();

    G4Material* rock = new G4Material("StandardRock",
                                      SimConfig::ROCK_Z_EFF,
                                      SimConfig::ROCK_A_EFF * (g/mole),
                                      SimConfig::ROCK_DENSITY * (g/cm3));
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

    // -----------------------------------------------------------------------
    // Geometry layout
    //
    //   z =  +5 cm  ← world top (5 cm above rock entry)
    //   z =   0     ← rock entry face  (muons enter here)
    //   z =  -D     ← scoring plane at depth D
    //   z =  -rockDepth_cm ← rock bottom (deepest scoring depth + margin)
    //   z = -(rockDepth_cm + 5 cm) ← world bottom
    //
    // Scoring planes are children of the ROCK logical volume.
    // Their local z inside rock = world_z - rockCentreZ_world
    //   = -D - (-rockHalfZ) = rockHalfZ - D
    //   (positive = upper part of rock, negative = lower)
    // -----------------------------------------------------------------------

    const auto& depths   = SimConfig::SCORING_DEPTHS_CM;
    double maxDepth_cm   = depths.back();
    double rockDepth_cm  = maxDepth_cm + 5.0;    // 5 cm margin below deepest plane
    double rockHalfZ_cm  = rockDepth_cm / 2.0;
    double rockHalfZ     = rockHalfZ_cm * cm;

    // Rock centre in world: entry face at z=0, so centre at z = -rockHalfZ
    double rockCentreZ   = -rockHalfZ;

    double halfXY        = SimConfig::SLAB_HALF_XY_CM * cm;

    // World: 5 cm headroom above rock, 5 cm below rock bottom
    double worldHalfZ = rockDepth_cm * cm + 10.0*cm;

    // -----------------------------------------------------------------------
    // World
    // -----------------------------------------------------------------------
    G4Box* worldBox = new G4Box("World", halfXY + 1*cm, halfXY + 1*cm, worldHalfZ);
    G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, air, "World");
    worldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VPhysicalVolume* worldPhys =
        new G4PVPlacement(nullptr, G4ThreeVector(), worldLog, "World",
                          nullptr, false, 0, true);

    // -----------------------------------------------------------------------
    // Rock slab — placed inside World, centre at z = -rockHalfZ (world)
    // -----------------------------------------------------------------------
    G4Box* rockBox = new G4Box("Rock", halfXY, halfXY, rockHalfZ);
    G4LogicalVolume* rockLog = new G4LogicalVolume(rockBox, rock, "Rock");
    G4VisAttributes* rockVis = new G4VisAttributes(G4Colour(0.6, 0.5, 0.4, 0.3));
    rockVis->SetForceSolid(false);
    rockLog->SetVisAttributes(rockVis);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, rockCentreZ),
                      rockLog, "Rock", worldLog, false, 0, true);

    // -----------------------------------------------------------------------
    // Geometry summary
    // -----------------------------------------------------------------------
    G4cout << "\n[Detector] Geometry\n"
           << "  World half-Z     : " << worldHalfZ/cm   << " cm\n"
           << "  Rock entry face  : z =    0.0 cm (world)\n"
           << "  Rock bottom      : z = " << -rockDepth_cm << " cm (world)\n"
           << "  Rock half-Z      : " << rockHalfZ_cm    << " cm\n"
           << "  Rock half-XY     : " << SimConfig::SLAB_HALF_XY_CM << " cm\n"
           << "  Rock density     : " << SimConfig::ROCK_DENSITY << " g/cm3\n\n";

    // -----------------------------------------------------------------------
    // Scoring planes — placed INSIDE the rock logical volume
    //
    // Rock local frame: centre of rock is origin.
    //   Top of rock  (world z=0)           → rock local z = +rockHalfZ_cm
    //   Depth D from entry face (world z=-D)→ rock local z = +rockHalfZ_cm - D
    // -----------------------------------------------------------------------
    const double planeHalfZ = 0.5 * mm;
    fScoringVolumes.clear();

    std::vector<G4Colour> colours = {
        G4Colour(1,0,0),     G4Colour(1,0.5,0),   G4Colour(1,1,0),
        G4Colour(0,1,0),     G4Colour(0,1,1),      G4Colour(0,0,1),
        G4Colour(0.5,0,1),   G4Colour(1,0,1),      G4Colour(0.7,0.7,0),
        G4Colour(0,0.7,0.7), G4Colour(0.7,0,0.7),  G4Colour(0.5,0.5,0.5)
    };

    G4cout << "[Detector] Scoring planes (children of Rock volume):\n"
           << "  Idx | Depth_cm |    m.w.e. | z_world_cm | z_rock_cm\n"
           << "  ----+----------+-----------+------------+----------\n";

    for (int i = 0; i < (int)depths.size(); ++i) {
        double d_cm      = depths[i];
        double mwe       = d_cm * SimConfig::ROCK_DENSITY / 100.0;   // g/cm2 -> m.w.e.
        double z_world   = -d_cm;                                      // cm, world frame
        double z_rock_cm = rockHalfZ_cm - d_cm;                       // cm, rock frame

        std::ostringstream nm;
        nm << "Score_" << (int)std::round(d_cm) << "cm";

        G4Box* box = new G4Box(nm.str(), halfXY, halfXY, planeHalfZ);
        G4LogicalVolume* lv = new G4LogicalVolume(box, air, nm.str());

        G4VisAttributes* va =
            new G4VisAttributes(colours[i % (int)colours.size()]);
        va->SetForceWireframe(true);
        lv->SetVisAttributes(va);

        // Place inside rockLog — position in ROCK's local coordinates
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, z_rock_cm * cm),
                          lv, nm.str(), rockLog, false, i, true);

        fScoringVolumes.push_back(lv);

        G4cout << "  " << std::setw(3) << i
               << " | " << std::setw(8) << d_cm
               << " | " << std::setw(9) << std::fixed << std::setprecision(3) << mwe
               << " | " << std::setw(10) << z_world
               << " | " << std::setw(9) << z_rock_cm << "\n";
    }
    G4cout << "\n";

    // -----------------------------------------------------------------------
    // Production cuts region for rock
    // -----------------------------------------------------------------------
    G4Region* rgn = new G4Region("RockRegion");
    rockLog->SetRegion(rgn);
    rgn->AddRootLogicalVolume(rockLog);
    G4ProductionCuts* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(SimConfig::PRODUCTION_CUT_MM * mm);
    rgn->SetProductionCuts(cuts);

    return worldPhys;
}

const std::vector<G4LogicalVolume*>& DetectorConstruction::GetScoringVolumes() const
{
    return fScoringVolumes;
}
