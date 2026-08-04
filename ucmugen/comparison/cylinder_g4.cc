// Item 1: cylindrical generation surfaces compared inside Geant4.
//
// EcoMug's cylindrical surface is the lateral surface only: its area is
// dphi*r*h and its J' carries sin^2(theta), so a muon arriving from directly
// overhead is never generated. For a tall thin detector that is harmless. For a
// squat one, the vertical component is exactly the part that matters.
//
// This is a three-way comparison rather than a two-way one, so that it shows
// which answer is right rather than only that two differ:
//
//   (a) UCMuGen on a large flat plane above the detector  -- the reference,
//       since a horizontal plane has no cap ambiguity at all
//   (b) UCMuGen on a cylinder enclosing the detector, caps included
//       -- must agree with (a); the rate into a detector cannot depend on the
//          surface chosen to generate through
//   (c) EcoMug on the same cylinder, lateral surface only
//
// Both generators are given the IDENTICAL differential flux, so nothing here
// depends on their default parametrisations differing.
//
// Muons are tracked by Geant4 through a vacuum world with transportation only:
// what is being compared is geometric acceptance, so adding physics processes
// would only add noise to the thing under test.

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#define UCMUGEN_WITH_GEANT4
#include "UCMuGen.h"
#include "EcoMug.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

using namespace ucmugen;

// ---- geometry, in centimetres to match UCMuGen ------------------------------
static const double kDetHalfX = 50.0;    // 1 m x 1 m plate
static const double kDetHalfY = 50.0;
static const double kDetHalfZ = 5.0;     // 10 cm thick
static const double kCylR = 150.0;       // enclosing cylinder, radius 1.5 m
static const double kCylH = 100.0;       // height 1 m  -> h/r = 0.67, squat
static const double kPlaneHalf = 350.0;  // reference plane, 7 m x 7 m
static const double kPlaneZ = 100.0;     // 1 m above the detector
// Sanity: a 70 deg muon from z=100 travels 100*tan(70)=275 cm laterally,
// plus the 50 cm detector half-width = 325 cm < 350 cm, so the plane
// cannot truncate the acceptance. Getting this wrong silently lowers the
// reference and makes every comparison against it meaningless.

static const double kPlaneHalf2 = 250.0;  // second reference, 5 m x 5 m
static const double kPlaneZ2 = 50.0;      // 0.5 m above the detector
static const double kPmin = 1.0, kPmax = 1000.0;      // GeV/c
static const double kThMax = 70.0 * kPi / 180.0;

// The shared physics.
static double sharedJ(double p_GeV, double theta_rad) {
  return flux::intensity(Spectrum::Guan, p_GeV, std::cos(theta_rad));
}

// ---- counting ---------------------------------------------------------------
static long long g_hits = 0;
static bool g_hit_this_event = false;
static double g_sum_costheta_hit = 0.0;

class Stepping : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* step) override {
    if (g_hit_this_event) return;
    const auto* vol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (vol && vol->GetName() == "Det") {
      g_hit_this_event = true;
      ++g_hits;
    }
  }
};

class Geometry : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vac = nist->FindOrBuildMaterial("G4_Galactic");
    auto* world = new G4Box("World", 10.0 * m, 10.0 * m, 10.0 * m);
    auto* wlog = new G4LogicalVolume(world, vac, "World");
    auto* wphys = new G4PVPlacement(nullptr, {}, wlog, "World", nullptr,
                                    false, 0);
    auto* det = new G4Box("Det", kDetHalfX * cm, kDetHalfY * cm,
                          kDetHalfZ * cm);
    auto* dlog = new G4LogicalVolume(det, vac, "Det");
    new G4PVPlacement(nullptr, {}, dlog, "Det", wlog, false, 0);
    return wphys;
  }
};

class Physics : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override {
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
  }
  void ConstructProcess() override { AddTransportation(); }
};

// ---- the three generators ---------------------------------------------------
enum class Mode { UcPlane, UcCylinder, UcCylinderNoCaps, EcoCylinder,
                  UcPlane2 };

class Primaries : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Primaries(Mode m) : mode_(m) {
    gun_ = new G4ParticleGun(1);
    const double e_min = std::sqrt(kPmin * kPmin + kMuonMass * kMuonMass);
    const double e_max = std::sqrt(kPmax * kPmax + kMuonMass * kMuonMass);
    ucm_.setSpectrum(Spectrum::Guan)
        .setEnergyRange(e_min, e_max)
        .setThetaRange(0.0, kThMax)
        .setSeed(20260730);
    if (mode_ == Mode::UcPlane) {
      ucm_.setSurface(std::make_shared<Plane>(kPlaneHalf, kPlaneHalf,
                                              Vec3{0, 0, kPlaneZ}));
    } else if (mode_ == Mode::UcPlane2) {
      ucm_.setSurface(std::make_shared<Plane>(kPlaneHalf2, kPlaneHalf2,
                                              Vec3{0, 0, kPlaneZ2}));
    } else if (mode_ == Mode::UcCylinder) {
      ucm_.setSurface(std::make_shared<Cylinder>(kCylR, kCylH, Vec3{0, 0, 0},
                                                 /*caps=*/true));
    } else if (mode_ == Mode::UcCylinderNoCaps) {
      ucm_.setSurface(std::make_shared<Cylinder>(kCylR, kCylH, Vec3{0, 0, 0},
                                                 /*caps=*/false));
    } else {
      eco_.SetUseCylinder();
      eco_.SetCylinderRadius(kCylR * 1.0e-2);      // EcoMug works in metres
      eco_.SetCylinderHeight(kCylH * 1.0e-2);
      eco_.SetCylinderCenterPosition({{0.0, 0.0, 0.0}});
      eco_.SetDifferentialFlux(&sharedJ);
      eco_.SetMinimumMomentum(kPmin);
      eco_.SetMaximumMomentum(kPmax);
      eco_.SetMinimumTheta(0.0);
      eco_.SetMaximumTheta(kThMax);
      eco_.SetSeed(20260730);
    }
  }
  ~Primaries() override { delete gun_; }

  void GeneratePrimaries(G4Event* evt) override {
    g_hit_this_event = false;
    if (mode_ == Mode::EcoCylinder) {
      eco_.Generate();
      const auto pos = eco_.GetGenerationPosition();
      std::array<double, 3> p3;
      eco_.GetGenerationMomentum(p3);
      const double pmag = eco_.GetGenerationMomentum();
      const double ke = std::sqrt(pmag * pmag + kMuonMass * kMuonMass)
                      - kMuonMass;
      gun_->SetParticleDefinition(
          G4ParticleTable::GetParticleTable()->FindParticle(
              eco_.GetCharge() > 0 ? -13 : 13));
      gun_->SetParticlePosition(G4ThreeVector(pos[0] * CLHEP::m,
                                              pos[1] * CLHEP::m,
                                              pos[2] * CLHEP::m));
      gun_->SetParticleMomentumDirection(
          G4ThreeVector(p3[0], p3[1], p3[2]).unit());
      gun_->SetParticleEnergy(ke * CLHEP::GeV);
      gun_->GeneratePrimaryVertex(evt);
      last_cos_ = -G4ThreeVector(p3[0], p3[1], p3[2]).unit().z();
    } else {
      const Muon mu = ucm_.generate();
      last_cos_ = -mu.direction.z;
      gun_->SetParticleDefinition(
          G4ParticleTable::GetParticleTable()->FindParticle(mu.pdg));
      gun_->SetParticlePosition(G4ThreeVector(mu.position.x * CLHEP::cm,
                                              mu.position.y * CLHEP::cm,
                                              mu.position.z * CLHEP::cm));
      gun_->SetParticleMomentumDirection(
          G4ThreeVector(mu.direction.x, mu.direction.y, mu.direction.z));
      gun_->SetParticleEnergy(mu.kinetic * CLHEP::GeV);
      gun_->GeneratePrimaryVertex(evt);
    }
  }

  double lastCos() const { return last_cos_; }
  Generator& ucm() { return ucm_; }
  EcoMug& eco() { return eco_; }

 private:
  Mode mode_;
  G4ParticleGun* gun_;
  Generator ucm_;
  EcoMug eco_;
  double last_cos_ = 0.0;
};

class Actions : public G4VUserActionInitialization {
 public:
  explicit Actions(Mode m) : mode_(m) {}
  void Build() const override {
    SetUserAction(new Primaries(mode_));
    SetUserAction(new Stepping);
  }
 private:
  Mode mode_;
};

// Surface rate from each tool's own estimator. Deliberately computed with
// plain objects that never touch Geant4, so no particle gun is constructed
// before the run manager has a physics list.
static double surfaceRate(Mode mode, double& err) {
  if (mode == Mode::EcoCylinder) {
    EcoMug e;
    e.SetUseCylinder();
    e.SetCylinderRadius(kCylR * 1.0e-2);
    e.SetCylinderHeight(kCylH * 1.0e-2);
    e.SetCylinderCenterPosition({{0.0, 0.0, 0.0}});
    e.SetDifferentialFlux(&sharedJ);
    e.SetMinimumMomentum(kPmin);
    e.SetMaximumMomentum(kPmax);
    e.SetMinimumTheta(0.0);
    e.SetMaximumTheta(kThMax);
    // Average over seeds: the custom-J estimator samples uniformly in momentum
    // over a steeply falling spectrum and is noisy at modest statistics.
    double sum = 0.0;
    const int kSeeds = 8;
    double vals[kSeeds];
    for (int i = 0; i < kSeeds; ++i) {
      e.SetSeed(20260730 + 7919 * i);
      double r, e1;
      e.GetAverageGenRateAndError(r, e1, 20000000);
      vals[i] = r;
      sum += r;
    }
    // Custom-J integration returns a rate per unit area in the units of J
    // (cm^-2 s^-1); GetGenSurfaceArea is in m^2.
    double var = 0.0;
    for (int i = 0; i < kSeeds; ++i)
      var += (vals[i] - sum / kSeeds) * (vals[i] - sum / kSeeds);
    const double scale = e.GetGenSurfaceArea() * 1.0e4;
    err = std::sqrt(var / (kSeeds - 1) / kSeeds) * scale;
    return (sum / kSeeds) * scale;
  }
  const double e_min = std::sqrt(kPmin * kPmin + kMuonMass * kMuonMass);
  const double e_max = std::sqrt(kPmax * kPmax + kMuonMass * kMuonMass);
  Generator g;
  g.setSpectrum(Spectrum::Guan)
      .setEnergyRange(e_min, e_max)
      .setThetaRange(0.0, kThMax)
      .setSeed(20260730);
  if (mode == Mode::UcPlane)
    g.setSurface(std::make_shared<Plane>(kPlaneHalf, kPlaneHalf,
                                         Vec3{0, 0, kPlaneZ}));
  else if (mode == Mode::UcPlane2)
    g.setSurface(std::make_shared<Plane>(kPlaneHalf2, kPlaneHalf2,
                                         Vec3{0, 0, kPlaneZ2}));
  else
    g.setSurface(std::make_shared<Cylinder>(
        kCylR, kCylH, Vec3{0, 0, 0}, mode == Mode::UcCylinder));
  double r;
  g.rateAndError(r, err, 8000000);
  return r;
}

int main(int argc, char** argv) {
  // One case per process: Geant4 does not take kindly to several run managers
  // in one run, and this keeps each measurement independent.
  const int which = (argc > 1) ? std::atoi(argv[1]) : 0;
  const int nev = (argc > 2) ? std::atoi(argv[2]) : 200000;
  const Mode modes[] = {Mode::UcPlane, Mode::UcCylinder,
                        Mode::UcCylinderNoCaps, Mode::EcoCylinder,
                        Mode::UcPlane2};
  const char* names[] = {"UCMuGen, flat plane (reference)",
                         "UCMuGen, cylinder with caps",
                         "UCMuGen, cylinder lateral only",
                         "EcoMug,  cylinder (lateral only)",
                         "UCMuGen, plane #2 (indep. reference)"};
  const Mode mode = modes[which];

  auto* rm = new G4RunManager;
  rm->SetUserInitialization(new Geometry);
  rm->SetUserInitialization(new Physics);
  rm->SetUserInitialization(new Actions(mode));
  rm->Initialize();
  g_hits = 0;
  rm->BeamOn(nev);

  double surf_err = 0.0;
  const double surf = surfaceRate(mode, surf_err);
  const double frac = double(g_hits) / double(nev);
  // Binomial error on the hit fraction, combined with the Monte Carlo error on
  // the surface rate. Omitting the second term understates the uncertainty and
  // would make a 1% agreement look like a 2-sigma discrepancy.
  const double frac_rel = (g_hits > 0) ? 1.0 / std::sqrt(double(g_hits)) : 0.0;
  const double surf_rel = (surf > 0.0) ? surf_err / surf : 0.0;
  const double det = surf * frac;
  const double det_rel = std::sqrt(frac_rel * frac_rel + surf_rel * surf_rel);
  std::printf("RESULT %d | %-34s | hits %8lld / %d | frac %10.6f%% | "
              "surface %12.6g +- %.3g Hz | detector %12.6g +- %.3g Hz "
              "(%.2f%%)\n",
              which, names[which], g_hits, nev, 100.0 * frac, surf, surf_err,
              det, det * det_rel, 100.0 * det_rel);
  delete rm;
  return 0;
}
