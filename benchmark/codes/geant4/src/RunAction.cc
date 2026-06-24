#include "RunAction.hh"
#include "SimConfig.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <sys/stat.h>
#include <sys/time.h>

// ── Shared state ─────────────────────────────────────────────────────────────
static std::ofstream gMuonCSV, gSecCSV, gStopCSV;
static std::string   gBaseName;
static int           gMuonHits  = 0;
static int           gSecHits   = 0;
static int           gStopCount = 0;
static double        gEmin      =  1e18;
static double        gEmax      = -1e18;
static struct timeval gT0;

void RunAction_WriteRow(const std::string& stream, const std::string& row)
{
    if (stream == "muon" && gMuonCSV.is_open()) { gMuonCSV << row; ++gMuonHits;  }
    if (stream == "sec"  && gSecCSV .is_open()) { gSecCSV  << row; ++gSecHits;   }
    if (stream == "stop" && gStopCSV.is_open()) { gStopCSV << row; ++gStopCount; }
}

void RunAction_UpdateEnergy(double initKE_GeV)
{
    if (initKE_GeV < gEmin) gEmin = initKE_GeV;
    if (initKE_GeV > gEmax) gEmax = initKE_GeV;
}

static void MkdirP(const std::string& dir)
{
    std::string cur;
    for (char c : dir + "/") {
        if (c == '/') { if (!cur.empty()) ::mkdir(cur.c_str(), 0755); }
        cur += c;
    }
}

RunAction::RunAction()  {}
RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
    gMuonHits = gSecHits = gStopCount = 0;
    gEmin =  1e18;
    gEmax = -1e18;
    gettimeofday(&gT0, nullptr);

    // Timestamp
    std::time_t now = std::time(nullptr);
    char tsbuf[32];
    std::strftime(tsbuf, sizeof(tsbuf), "%Y%m%d_%H%M%S", std::localtime(&now));

    // Mode tag
    std::string modeTag;
    switch (SimConfig::SOURCE_MODE) {
        case SimConfig::SourceMode::FILE:     modeTag = "file";     break;
        case SimConfig::SourceMode::POWERLAW: modeTag = "powerlaw"; break;
        case SimConfig::SourceMode::ECOMUG:   modeTag = "ecomug";   break;
    }

    MkdirP(SimConfig::OUTPUT_DIR);

    gBaseName = "run_" + modeTag + "_" + tsbuf;
    std::string muonPath = SimConfig::OUTPUT_DIR + "/" + gBaseName + "_muons.csv";
    std::string secPath  = SimConfig::OUTPUT_DIR + "/" + gBaseName + "_secondaries.csv";

    std::string stopPath  = SimConfig::OUTPUT_DIR + "/" + gBaseName + "_stopped.csv";

    gMuonCSV.open(muonPath);
    gSecCSV .open(secPath);
    gStopCSV.open(stopPath);

    gMuonCSV << "EventID,DepthCm,MWE,PDG,InitKEGeV,InitThetaDeg,InitPhiDeg,"
                "ExitKEGeV,ExitPxGeVc,ExitPyGeVc,ExitPzGeVc,ExitXcm,ExitYcm,"
                "ELossIonGeV,ELossBremGeV,ELossPairGeV,ELossNuclGeV,ELossTotalGeV,"
                "AngleScatDeg,LatDispCm\n";

    gSecCSV  << "EventID,PDG,DepthCm,KE_MeV,Px_MeV,Py_MeV,Pz_MeV,Xcm,Ycm\n";
    gStopCSV << "EventID,PDG,InitKEGeV,InitThetaDeg,InitPhiDeg,StopDepthCm\n";

    G4cout << "[RunAction] Output base: " << SimConfig::OUTPUT_DIR
           << "/" << gBaseName << "\n";
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    struct timeval t1; gettimeofday(&t1, nullptr);
    double elapsed = (t1.tv_sec - gT0.tv_sec)
                   + (t1.tv_usec - gT0.tv_usec) * 1e-6;

    int nEvents = run->GetNumberOfEvent();

    if (gMuonCSV.is_open()) gMuonCSV.close();
    if (gSecCSV .is_open()) gSecCSV .close();
    if (gStopCSV.is_open()) gStopCSV.close();

    // ── Write run summary TXT (same base name) ────────────────────────────────
    std::string summaryPath = SimConfig::OUTPUT_DIR + "/" + gBaseName + "_summary.txt";
    std::ofstream sf(summaryPath);

    std::string modeStr;
    switch (SimConfig::SOURCE_MODE) {
        case SimConfig::SourceMode::FILE:
            modeStr = "FILE  -> " + SimConfig::SOURCE_FILE; break;
        case SimConfig::SourceMode::POWERLAW:
            modeStr = "POWERLAW  E^-" + std::to_string(SimConfig::SPEC_GAMMA)
                    + "  [" + std::to_string((int)SimConfig::SPEC_EMIN_GEV)
                    + "-"   + std::to_string((int)SimConfig::SPEC_EMAX_GEV) + " GeV]"; break;
        case SimConfig::SourceMode::ECOMUG:
            modeStr = "ECOMUG  pMin=" + std::to_string((int)(SimConfig::ECOMUG_PMIN_MEV/1000))
                    + " GeV/c  thetaMax="
                    + std::to_string((int)(SimConfig::ECOMUG_MAX_THETA_RAD*180.0/3.14159)) + " deg"; break;
    }

    bool hasHits = (gEmin < 1e17);
    sf << "=================================================\n"
       << "  MuonRock Run Summary\n"
       << "=================================================\n"
       << "  Run name   : " << gBaseName        << "\n"
       << "  Source     : " << modeStr           << "\n"
       << "  Events     : " << nEvents           << "\n"
       << "  Muon plane crossings : " << gMuonHits         << "\n"
       << "  Sec hits             : " << gSecHits           << "\n"
       << "  Stopped muons        : " << gStopCount         << "\n";
    if (hasHits)
        sf << "  KE min     : " << std::fixed << std::setprecision(1) << gEmin << " GeV\n"
           << "  KE max     : " << std::fixed << std::setprecision(1) << gEmax << " GeV\n";
    sf << "  Elapsed    : " << std::fixed << std::setprecision(2)
                            << elapsed << " s\n"
       << "  Rate       : " << std::fixed << std::setprecision(1)
                            << (elapsed > 0 ? nEvents/elapsed : 0) << " evt/s\n"
       << "  CPU/evt    : " << std::fixed << std::setprecision(3)
                            << (nEvents > 0 ? elapsed/nEvents*1000 : 0) << " ms\n"
       << "-------------------------------------------------\n"
       << "  Files:\n"
       << "    " << gBaseName << "_muons.csv\n"
       << "    " << gBaseName << "_secondaries.csv\n"
       << "    " << gBaseName << "_stopped.csv\n"
       << "    " << gBaseName << "_summary.txt\n"
       << "=================================================\n";
    sf.close();

    // ── Terminal summary ──────────────────────────────────────────────────────
    std::cerr
        << "\n╔══════════════════════════════════════════╗\n"
        << "║          MuonRock Run Summary            ║\n"
        << "╠══════════════════════════════════════════╣\n"
        << "║  Events    : " << std::setw(8) << nEvents  << "                  ║\n"
        << "║  Muon crossings : " << std::setw(8) << gMuonHits<< "             ║\n"
        << "║  Sec hits       : " << std::setw(8) << gSecHits  << "             ║\n"
        << "║  Stopped muons  : " << std::setw(8) << gStopCount<< "             ║\n";
    if (hasHits)
        std::cerr
        << "║  KE min    : " << std::setw(7) << std::fixed << std::setprecision(1)
                             << gEmin << " GeV"           << "                 ║\n"
        << "║  KE max    : " << std::setw(7) << std::fixed << std::setprecision(1)
                             << gEmax << " GeV"           << "                 ║\n";
    std::cerr
        << "║  Time      : " << std::setw(7) << std::fixed << std::setprecision(1)
                             << elapsed << " s"            << "                  ║\n"
        << "║  Rate      : " << std::setw(7) << std::fixed << std::setprecision(1)
                             << (elapsed>0 ? nEvents/elapsed : 0) << " evt/s"
                             << "              ║\n"
        << "╠══════════════════════════════════════════╣\n"
        << "║  " << SimConfig::OUTPUT_DIR << "/" << gBaseName << "\n"
        << "╚══════════════════════════════════════════╝\n\n";
}