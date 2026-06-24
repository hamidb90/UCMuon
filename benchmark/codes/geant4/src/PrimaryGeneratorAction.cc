#include "PrimaryGeneratorAction.hh"
#include "SimConfig.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iomanip>

#if defined(USE_ECOMUG)
#include "EcoMug.h"
#endif

// ─────────────────────────────────────────────────────────────────────────────
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fGun = new G4ParticleGun(1);

    if (SimConfig::SOURCE_MODE == SimConfig::SourceMode::FILE)
        LoadSourceFile();

#if defined(USE_ECOMUG)
    if (SimConfig::SOURCE_MODE == SimConfig::SourceMode::ECOMUG) {
        EcoMug* em = new EcoMug();
        em->SetUseSky();
        em->SetSkySize({SimConfig::ECOMUG_PLANE_HALFX_MM * 2.0,
                        SimConfig::ECOMUG_PLANE_HALFY_MM * 2.0});
        em->SetSkyCenterPosition({0.0, 0.0,
                                  SimConfig::ECOMUG_PLANE_Z_MM});
        em->SetMinimumMomentum(SimConfig::ECOMUG_PMIN_MEV / 1000.0);
        em->SetMaximumMomentum(SimConfig::ECOMUG_PMAX_MEV / 1000.0);
        em->SetMaximumTheta(SimConfig::ECOMUG_MAX_THETA_RAD);
        fEcoMugPtr = static_cast<void*>(em);
        G4cout << "[PrimaryGenerator] EcoMug SKY (flat rect plane):\n"
               << "  Size     : " << SimConfig::ECOMUG_PLANE_HALFX_MM*2.0/10.0
               << " x " << SimConfig::ECOMUG_PLANE_HALFY_MM*2.0/10.0 << " cm\n"
               << "  z        : " << SimConfig::ECOMUG_PLANE_Z_MM << " mm\n"
               << "  KE range : " << SimConfig::SPEC_EMIN_GEV << " - "
               << SimConfig::SPEC_EMAX_GEV << " GeV\n"
               << "  thetaMax : " << SimConfig::ECOMUG_MAX_THETA_RAD*180.0/M_PI << " deg\n"
               << "  Angular  : cos^(N+1)(theta), N=max(0.1, 2.856-0.655*ln(p/GeV)) [EcoMug]\n";
    }
#endif
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
#if defined(USE_ECOMUG)
    if (fEcoMugPtr) delete static_cast<EcoMug*>(fEcoMugPtr);
#endif
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    switch (SimConfig::SOURCE_MODE) {
        case SimConfig::SourceMode::FILE:     FireFromFile(event);  break;
        case SimConfig::SourceMode::POWERLAW: FirePowerLaw(event);  break;
        case SimConfig::SourceMode::ECOMUG:   FireEcoMug(event);    break;
    }
}

// ── FILE ─────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::LoadSourceFile()
{
    const auto& cfg = SimConfig::CSV_COLS;
    std::ifstream f(SimConfig::SOURCE_FILE);
    if (!f.is_open())
        throw std::runtime_error("Cannot open source file: " + SimConfig::SOURCE_FILE);

    if (cfg.has_header) { std::string d; std::getline(f, d); }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::vector<std::string> cols;
        std::string tok;
        while (ss >> tok) cols.push_back(tok);

        int need = std::max({cfg.col_pdg, cfg.col_x, cfg.col_y, cfg.col_z,
                             cfg.col_px,  cfg.col_py, cfg.col_pz, cfg.col_ke}) + 1;
        if ((int)cols.size() < need) continue;

        MuonRecord r;
        r.pdg    = cfg.pdg_is_string
                   ? (cols[cfg.col_pdg] == "mu-" ? 13 : -13)
                   : std::stoi(cols[cfg.col_pdg]);
        r.x_mm   = std::stod(cols[cfg.col_x])  * cfg.pos_to_mm;
        r.y_mm   = std::stod(cols[cfg.col_y])  * cfg.pos_to_mm;
        r.z_mm   = std::stod(cols[cfg.col_z])  * cfg.pos_to_mm;
        r.px_MeV = std::stod(cols[cfg.col_px]);
        r.py_MeV = std::stod(cols[cfg.col_py]);
        r.pz_MeV = std::stod(cols[cfg.col_pz]);
        r.ke_MeV = std::stod(cols[cfg.col_ke]) * cfg.ke_to_MeV;
        fRecords.push_back(r);
    }
    if (fRecords.empty())
        throw std::runtime_error("No valid records in: " + SimConfig::SOURCE_FILE);

    // Clamp XY to inside slab face; override z to rock entry face
    double margin_mm = SimConfig::SLAB_HALF_XY_CM * 10.0 - 10.0;
    int nClamped = 0;
    for (auto& r : fRecords) {
        double ox = r.x_mm, oy = r.y_mm;
        r.x_mm = std::max(-margin_mm, std::min(margin_mm, r.x_mm));
        r.y_mm = std::max(-margin_mm, std::min(margin_mm, r.y_mm));
        r.z_mm = 1.0;  // 1 mm above z=0 rock entry face
        if (ox != r.x_mm || oy != r.y_mm) ++nClamped;
    }

    G4cout << "[PrimaryGenerator] FILE mode: " << fRecords.size()
           << " muons from " << SimConfig::SOURCE_FILE << "\n"
           << "  First : PDG=" << fRecords[0].pdg
           << "  KE=" << fRecords[0].ke_MeV << " MeV\n"
           << "  XY clamped to +-" << margin_mm/10.0 << " cm: "
           << nClamped << "/" << fRecords.size() << " records adjusted\n"
           << "  All z forced to 1 mm (rock entry face)\n";
}

void PrimaryGeneratorAction::FireFromFile(G4Event* event)
{
    const MuonRecord& r = fRecords[fIndex++ % (int)fRecords.size()];
    SetGunFromRecord(r.x_mm, r.y_mm, r.z_mm,
                     r.px_MeV, r.py_MeV, r.pz_MeV, r.ke_MeV, r.pdg);
    RecordInitialState(r.ke_MeV, r.px_MeV, r.py_MeV, r.pz_MeV, r.pdg);
    fGun->GeneratePrimaryVertex(event);
}

// ── POWERLAW ──────────────────────────────────────────────────────────────────
double PrimaryGeneratorAction::SamplePowerLaw()
{
    double Emin = SimConfig::SPEC_EMIN_GEV * 1000.0;
    double Emax = SimConfig::SPEC_EMAX_GEV * 1000.0;
    double n    = 1.0 - SimConfig::SPEC_GAMMA;
    double u    = G4UniformRand();
    if (std::fabs(n) < 1e-9)
        return Emin * std::exp(u * std::log(Emax / Emin));
    return std::pow(u*(std::pow(Emax,n) - std::pow(Emin,n)) + std::pow(Emin,n), 1.0/n);
}

void PrimaryGeneratorAction::FirePowerLaw(G4Event* event)
{
    double ke  = SamplePowerLaw();
    int    pdg = (G4UniformRand() < SimConfig::POWERLAW_MU_MINUS_FRACTION) ? 13 : -13;
    double hXY = SimConfig::SLAB_HALF_XY_CM * cm / mm;
    double x   = (G4UniformRand()*2.0 - 1.0) * hXY;
    double y   = (G4UniformRand()*2.0 - 1.0) * hXY;
    double z   = 1.0;  // 1 mm above rock entry
    double m_mu = 105.6583755;
    double E    = ke + m_mu;
    double p_MeV = std::sqrt(E*E - m_mu*m_mu);

    // Sample theta from cos^(N+1)(theta)*sin(theta) via exact inverse CDF:
    //   cos(theta) = (1 - u*(1 - cosMax^(N+2)))^(1/(N+2))
    double p_GeV = p_MeV / 1000.0;
    double N     = std::max(0.1, 2.856 - 0.655 * std::log(p_GeV));
    double cosMax = std::cos(SimConfig::POWERLAW_MAX_THETA_RAD);
    double u      = G4UniformRand();
    double cosT   = std::pow(1.0 - u * (1.0 - std::pow(cosMax, N + 2.0)), 1.0 / (N + 2.0));
    double sinT   = std::sqrt(1.0 - cosT * cosT);
    double phi    = G4UniformRand() * 2.0 * M_PI;

    double px =  p_MeV * sinT * std::cos(phi);
    double py =  p_MeV * sinT * std::sin(phi);
    double pz = -p_MeV * cosT;   // downward

    SetGunFromRecord(x, y, z, px, py, pz, ke, pdg);
    RecordInitialState(ke, px, py, pz, pdg);
    fGun->GeneratePrimaryVertex(event);
}

// ── ECOMUG ────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::FireEcoMug(G4Event* event)
{
#if defined(USE_ECOMUG)
    EcoMug* em = static_cast<EcoMug*>(fEcoMugPtr);
    em->Generate();

    const auto& pos = em->GetGenerationPosition();  // mm
    double theta    = em->GetGenerationTheta();      // rad from vertical (0=straight down)
    double phi      = em->GetGenerationPhi();        // rad azimuth
    double p_GeV    = em->GetGenerationMomentum();   // GeV/c scalar
    int    charge   = em->GetCharge();               // +1=mu+, -1=mu-
    int    pdg      = (charge > 0) ? -13 : 13;      // mu+=PDG -13, mu-=PDG 13

    double p_MeV = p_GeV * 1000.0;
    double px    =  p_MeV * std::sin(theta) * std::cos(phi);
    double py    =  p_MeV * std::sin(theta) * std::sin(phi);
    double pz    =  p_MeV * std::cos(theta);  // EcoMug stores theta in [π/2,π] after M_PI flip, so cos(theta)<0 = downward

    double m_mu = 105.6583755;  // MeV/c2
    double ke   = std::sqrt(p_MeV*p_MeV + m_mu*m_mu) - m_mu;

    SetGunFromRecord(pos[0], pos[1], pos[2], px, py, pz, ke, pdg);
    RecordInitialState(ke, px, py, pz, pdg);
    fGun->GeneratePrimaryVertex(event);
#else
    G4Exception("PrimaryGeneratorAction::FireEcoMug", "SOURCE_MODE", FatalException,
                "ECOMUG mode requires -DUSE_ECOMUG=ON at cmake time.");
#endif
}

// ── Helpers ───────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::SetGunFromRecord(
    double x, double y, double z,
    double px, double py, double pz,
    double ke, int pdg)
{
    auto* pt = G4ParticleTable::GetParticleTable();
    fGun->SetParticleDefinition(pt->FindParticle(pdg == 13 ? "mu-" : "mu+"));
    fGun->SetParticlePosition(G4ThreeVector(x*mm, y*mm, z*mm));
    fGun->SetParticleEnergy(ke * MeV);
    fGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz).unit());
}

void PrimaryGeneratorAction::RecordInitialState(
    double ke, double px, double py, double pz, int pdg)
{
    fInitKE_MeV = ke;
    fInitPDG    = pdg;
    double p    = std::sqrt(px*px + py*py + pz*pz);
    fInitTheta  = (p > 0) ? std::acos(std::fabs(pz) / p) : 0.0;
    fInitPhi    = std::atan2(py, px);
}
