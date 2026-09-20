#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "SimConfig.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "FTFP_BERT.hh"
#include "G4ios.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

// ── Helpers ──────────────────────────────────────────────────────────────────
static void PrintUsage(const char* prog)
{
    std::cerr
        << "\nUsage:\n"
        << "  " << prog << " [options] <nEvents>\n\n"
        << "Options:\n"
        << "  -f <file>     Source file path (FILE mode only)\n"
        << "  -m <mode>     Source mode: file | powerlaw | ecomug\n"
        << "  -emin <GeV>   Min energy/momentum in GeV (powerlaw: KE, ecomug: p)\n"
        << "  -emax <GeV>   Max energy/momentum in GeV (powerlaw: KE, ecomug: p)\n"
        << "  -theta <deg>  Max zenith angle in degrees for powerlaw/ecomug  [default: 30]\n"
        << "  -o <dir>      Output directory  [default: outputs]\n"
        << "  -g            GUI mode — open Qt window (vis.mac must be in build dir)\n"
        << "  -q            Quiet mode — redirect Geant4 log to outputs/geant4.log\n"
        << "  -n <events>   Number of events  (same as positional arg)\n\n"
        << "Positional shorthand (bare numbers = [emin_GeV] [emax_GeV] nEvents):\n"
        << "  " << prog << " -m ecomug 1000           → nEvents=1000, defaults\n"
        << "  " << prog << " -m ecomug 1 100 1000     → emin=1 GeV, emax=100 GeV, nEvents=1000\n"
        << "  " << prog << " -f muons.txt 50000\n\n";
}

static std::string Timestamp()
{
    std::time_t t = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", std::localtime(&t));
    return buf;
}

static void MkdirP(const std::string& dir)
{
    // Simple single-level mkdir; for nested paths each segment is created.
    std::string cur;
    for (char c : dir + "/") {
        if (c == '/') {
            if (!cur.empty()) ::mkdir(cur.c_str(), 0755);
        }
        cur += c;
    }
}

// KE (GeV) → momentum (MeV/c) for muons
static double KEtoP_MeV(double ke_GeV)
{
    constexpr double m = 0.10566;  // muon mass GeV/c²
    return std::sqrt(ke_GeV * ke_GeV + 2.0 * ke_GeV * m) * 1000.0;
}

static void SetEnergyRange(double emin_GeV, double emax_GeV)
{
    SimConfig::SPEC_EMIN_GEV   = emin_GeV;
    SimConfig::SPEC_EMAX_GEV   = emax_GeV;
    SimConfig::ECOMUG_PMIN_MEV = KEtoP_MeV(emin_GeV);
    SimConfig::ECOMUG_PMAX_MEV = KEtoP_MeV(emax_GeV);
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char** argv)
{
    G4int       nEvents   = 1000;
    bool        quietMode = false;
    bool        guiMode   = false;
    std::string outDir    = "../outputs";
    std::vector<double> posArgs;  // bare numbers: [emin] [emax] nEvents

    // Initialise ECOMUG momentum limits from the default KE range
    SetEnergyRange(SimConfig::SPEC_EMIN_GEV, SimConfig::SPEC_EMAX_GEV);

    for (int i = 1; i < argc; ++i) {
        if      (std::strcmp(argv[i],"-f")==0 && i+1<argc) {
            SimConfig::SOURCE_FILE = argv[++i];
            SimConfig::SOURCE_MODE = SimConfig::SourceMode::FILE;
        }
        else if (std::strcmp(argv[i],"-m")==0 && i+1<argc) {
            std::string m = argv[++i];
            if      (m=="file")     SimConfig::SOURCE_MODE = SimConfig::SourceMode::FILE;
            else if (m=="powerlaw") SimConfig::SOURCE_MODE = SimConfig::SourceMode::POWERLAW;
            else if (m=="ecomug")   SimConfig::SOURCE_MODE = SimConfig::SourceMode::ECOMUG;
            else { std::cerr<<"[Error] Unknown mode: "<<m<<"\n"; PrintUsage(argv[0]); return 1; }
        }
        else if (std::strcmp(argv[i],"-theta")==0 && i+1<argc) {
            double rad = std::stod(argv[++i]) * M_PI / 180.0;
            SimConfig::POWERLAW_MAX_THETA_RAD = rad;
            SimConfig::ECOMUG_MAX_THETA_RAD   = rad;
        }
        else if (std::strcmp(argv[i],"-emin")==0 && i+1<argc) {
            SimConfig::SPEC_EMIN_GEV   = std::stod(argv[++i]);
            SimConfig::ECOMUG_PMIN_MEV = KEtoP_MeV(SimConfig::SPEC_EMIN_GEV);
        }
        else if (std::strcmp(argv[i],"-emax")==0 && i+1<argc) {
            SimConfig::SPEC_EMAX_GEV   = std::stod(argv[++i]);
            SimConfig::ECOMUG_PMAX_MEV = KEtoP_MeV(SimConfig::SPEC_EMAX_GEV);
        }
        else if (std::strcmp(argv[i],"-o")==0 && i+1<argc) { outDir = argv[++i]; }
        else if (std::strcmp(argv[i],"-g")==0)              { guiMode  = true; }
        else if (std::strcmp(argv[i],"-q")==0)              { quietMode = true; }
        else if (std::strcmp(argv[i],"-n")==0 && i+1<argc) { nEvents = std::stoi(argv[++i]); }
        else if (std::strcmp(argv[i],"-h")==0) { PrintUsage(argv[0]); return 0; }
        else { try { posArgs.push_back(std::stod(argv[i])); } catch(...) {
            std::cerr<<"[Error] Unrecognised argument: "<<argv[i]<<"\n";
            PrintUsage(argv[0]); return 1; } }
    }

    // Positional numbers: nEvents  |  emin nEvents  |  emin emax nEvents
    if (posArgs.size() == 1) {
        nEvents = (int)posArgs[0];
    } else if (posArgs.size() == 2) {
        SetEnergyRange(posArgs[0], SimConfig::SPEC_EMAX_GEV);
        nEvents = (int)posArgs[1];
    } else if (posArgs.size() >= 3) {
        SetEnergyRange(posArgs[0], posArgs[1]);
        nEvents = (int)posArgs[2];
    }

    if (SimConfig::SOURCE_MODE == SimConfig::SourceMode::FILE
        && SimConfig::SOURCE_FILE.empty()) {
        std::cerr<<"[Error] FILE mode selected but no source file given.\n"
                 <<"  Use: "<<argv[0]<<" -f /path/to/muons.txt <nEvents>\n\n";
        return 1;
    }

    // ── Create output directory ───────────────────────────────────────────────
    MkdirP(outDir);
    SimConfig::OUTPUT_DIR = outDir;

    std::string ts = Timestamp();

    // ── Redirect Geant4 output to log file if -q ──────────────────────────────
    std::ofstream logFile;
    std::streambuf* coutBuf = nullptr;
    std::string logPath = outDir + "/geant4_" + ts + ".log";

    if (quietMode) {
        logFile.open(logPath);
        if (!logFile.is_open()) {
            std::cerr << "[Warning] Cannot open log file " << logPath
                      << " — running verbose.\n";
            quietMode = false;
        } else {
            coutBuf = std::cout.rdbuf();
            std::cout.rdbuf(logFile.rdbuf());
        }
    }

    // ── Header (always to real stdout even in quiet mode) ─────────────────────
    auto& out = quietMode ? std::cerr : std::cout;
    out << "\n=== MuonRock v6 ===\n";
    switch(SimConfig::SOURCE_MODE){
        case SimConfig::SourceMode::FILE:
            out << "  Source : FILE -> " << SimConfig::SOURCE_FILE << "\n"; break;
        case SimConfig::SourceMode::POWERLAW:
            out << "  Source : POWERLAW  E^-" << SimConfig::SPEC_GAMMA
                << "  [" << SimConfig::SPEC_EMIN_GEV << "-"
                << SimConfig::SPEC_EMAX_GEV << " GeV]"
                << "  thetaMax=" << SimConfig::POWERLAW_MAX_THETA_RAD*180.0/M_PI << " deg\n"; break;
        case SimConfig::SourceMode::ECOMUG:
            out << "  Source : ECOMUG  KE=["
                << SimConfig::SPEC_EMIN_GEV << "-"
                << SimConfig::SPEC_EMAX_GEV << " GeV]"
                << "  thetaMax=" << SimConfig::ECOMUG_MAX_THETA_RAD*180.0/3.14159 << " deg\n"; break;
    }
    out << "  Events : " << nEvents << "\n"
        << "  OutDir : " << outDir  << "\n";
    if (quietMode)
        out << "  Log    : " << logPath << "\n";
    out << "==================\n\n";

    // ── Geant4 setup ──────────────────────────────────────────────────────────
    G4UIExecutive* ui = guiMode ? new G4UIExecutive(argc, argv, "Qt") : nullptr;

    auto* det    = new DetectorConstruction();
    auto* runMgr = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
    runMgr->SetUserInitialization(det);
    runMgr->SetUserInitialization(new FTFP_BERT());
    runMgr->SetUserInitialization(new ActionInitialization(det));
    runMgr->Initialize();

    G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/setCut "
                     + std::to_string(SimConfig::PRODUCTION_CUT_MM) + " mm");

    if (guiMode) {
        UI->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    } else {
        UI->ApplyCommand("/run/beamOn " + std::to_string(nEvents));
    }

    delete visManager;
    delete runMgr;

    // ── Restore stdout ────────────────────────────────────────────────────────
    if (quietMode && coutBuf) {
        std::cout.rdbuf(coutBuf);
        logFile.close();
        std::cout << "[Done] Geant4 log saved to: " << logPath << "\n";
    }

    return 0;
}
