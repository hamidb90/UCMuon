#!/usr/bin/env python3
"""
benchmark_analysis.py — MuonRock v6 / UCLouvain CP3
Multi-code muon propagation benchmark: Geant4, PHITS, MUSIC, PROPOSAL

Usage (Geant4 only):
    python3 benchmark_analysis.py --geant4 outputs/run_powerlaw_20260512_120000
    python3 benchmark_analysis.py --geant4 outputs/        # auto-picks latest run_*

Usage (multi-code comparison):
    python3 benchmark_analysis.py \\
        --geant4  outputs/run_powerlaw_... \\
        --phits   phits_results/ \\
        --music   music_results/ \\
        --proposal proposal_results/ \\
        --outdir  figures/

Backward-compatible positional form (Geant4 only):
    python3 benchmark_analysis.py outputs/run_powerlaw_...

Outputs:
    benchmark_summary.csv    per-depth statistics for every loaded code
    fig00_summary.png        2x3 overview panel (multi-code; only if >1 code)
    fig01_transmission.png   survival curves
    fig02_energy_loss.png    mean dE + process fractions + distributions
    fig03_angular.png        mean angle + Highland formula + distributions
    fig04_lateral.png        lateral displacement RMS vs depth
    fig05_exit_spectrum.png  exit KE histograms at each depth
    fig06_dedx.png           dE/dx vs initial energy + Bethe-Bloch reference
    fig07_charge_ratio.png   mu-/mu+ ratio vs depth
    fig08_secondaries.png    secondary types + KE spectrum  [Geant4 only]
    fig09_stopped.png        stopped muon KE + stop depth   [Geant4 only]
    fig10_timing.png         code speed comparison
    fig11_exit_position.png  exit (x,y) positions at each depth  [codes with position data]

Geant4 column names (written by RunAction.cc):
    muons CSV    : EventID,DepthCm,MWE,PDG,InitKEGeV,InitThetaDeg,InitPhiDeg,
                   ExitKEGeV,ExitPxGeVc,ExitPyGeVc,ExitPzGeVc,ExitXcm,ExitYcm,
                   ELossIonGeV,ELossBremGeV,ELossPairGeV,ELossNuclGeV,ELossTotalGeV,
                   AngleScatDeg,LatDispCm
    secondaries  : EventID,PDG,DepthCm,KE_MeV,Px_MeV,Py_MeV,Pz_MeV,Xcm,Ycm
    stopped      : EventID,PDG,InitKEGeV,InitThetaDeg,InitPhiDeg,StopDepthCm

PHITS / MUSIC / PROPOSAL common CSV format (normalize your output to this):
    muons CSV    : EventID,DepthCm,PDG,InitKE_GeV,ExitKE_GeV[,AngleScat_Deg,LatDisp_cm]
    summary CSV  : DepthCm,MWE,N_transmitted,Transmission_%,MeanExitKE_GeV
    timing file  : plain text containing "Elapsed : <seconds>"  (same as G4 summary)

PHITS known limitations vs Geant4 (documented in BENCHMARK_FEEDBACK.md §3):
    1. No per-event output — PHITS produces aggregate tallies only (T-Cross histograms).
       Per-muon InitKE, ExitKE, scatter angle, and position are not available.
       MCS angles and lateral displacement are marked — in all comparison tables.
    2. MCS model difference — PHITS uses Lynch-Molière (nspred=-2); Geant4 uses
       Urban MscModel (G4UrbanMscModel via FTFP_BERT). Both derive from Molière
       theory but differ in step-limited accumulation. Not a configuration error.
    3. Material I value — Geant4 rock: Z=11 single-element → I(Na)=149 eV (NIST).
       PHITS rock: Mg(75%)+O(25%) mixture → composite I≈138 eV. Effect on
       ionisation dE/dx < 0.3% at 300 GeV. Not the cause of the −12% KE discrepancy.
    4. PHITS dE/dx discrepancy — mean exit KE is −12.2% vs Geant4 at 530 MWE,
       growing monotonically with depth. Likely cause: different cross-section
       parametrisations for muon bremsstrahlung (imubrm) and pair production (imuppd)
       at 100–300 GeV. Requires dedicated single-layer test to confirm.
"""

import sys, os, re, glob, argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# ── Publication-quality global style ──────────────────────────────────────────
plt.rcParams.update({
    "font.family":          "serif",
    "font.serif":           ["Times New Roman", "Palatino Linotype", "DejaVu Serif"],
    "font.size":            11,
    "axes.labelsize":       13,
    "axes.titlesize":       14,
    "axes.titlepad":        8,
    "axes.labelpad":        6,
    "axes.linewidth":       1.2,
    "xtick.labelsize":      10,
    "ytick.labelsize":      10,
    "xtick.direction":      "in",
    "ytick.direction":      "in",
    "xtick.top":            True,
    "ytick.right":          True,
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,
    "xtick.major.size":     5,
    "xtick.minor.size":     3,
    "ytick.major.size":     5,
    "ytick.minor.size":     3,
    "legend.fontsize":      9,
    "legend.framealpha":    0.85,
    "legend.edgecolor":     "0.7",
    "lines.linewidth":      2.0,
    "grid.alpha":           0.3,
    "grid.linewidth":       0.6,
    "figure.dpi":           150,
    "savefig.dpi":          300,
    "savefig.bbox":         "tight",
})

# ── Constants ─────────────────────────────────────────────────────────────────
ROCK_DENSITY  = 2.65     # g/cm³ standard rock
M_MU          = 0.10566  # GeV/c²
X0_ROCK_GCM2  = 27.74    # radiation length of standard rock [g/cm²]

# Depths excluded from all comparison figures and tables.
# 1000 cm (26.5 MWE): Geant4's 5 GeV muon range ≈ 996 cm, so ~54% of 5 GeV muons
# pass the 10 m plane while all other codes stop them before 10 m (range ~860 cm).
# This Bragg-peak boundary effect creates a 7–9 pp transmission gap that propagates
# into exit KE and MCS angle comparisons at that depth, making codes incomparable.
EXCLUDE_DEPTH_CM = {1000}

PDG_NAMES = {
    22: r"$\gamma$", 11: r"$e^-$", -11: r"$e^+$", 13: r"$\mu^-$", -13: r"$\mu^+$",
    2212: "p", 2112: "n", 211: r"$\pi^+$", -211: r"$\pi^-$", 111: r"$\pi^0$",
    321: r"$K^+$", -321: r"$K^-$", 130: r"$K^0_L$", 310: r"$K^0_S$",
}

CODE_COLORS  = {"Geant4": "#2196F3", "PHITS": "#E91E63", "PHITS*": "#E91E63",
                "MUSIC":  "#FF9800", "PROPOSAL": "#4CAF50",
                "BB":     "#9C27B0", "UCMuon": "#00BCD4"}
LINE_STYLES  = ["-", "--", "-.", ":"]
MARKERS      = ["o", "s", "^", "D"]

def line_kw(idx, color, **extra):
    return dict(color=color, linewidth=2,
                linestyle=LINE_STYLES[idx % 4],
                marker=MARKERS[idx % 4], markersize=6, **extra)

# ── Theory curves ─────────────────────────────────────────────────────────────
def bethe_bloch_mwe(E_GeV):
    """dE/dx [GeV / m.w.e.] for muons in rock (Groom 2001: a + b*E)."""
    return 0.259 + 3.65e-4 * np.asarray(E_GeV)

def highland_deg(depth_mwe, p_GeV, beta=1.0):
    """
    Highland / Lynch-Dahl MCS formula: theta_0 in degrees.
    depth_mwe : rock thickness in m.w.e.  (= depth_cm * 2.65 / 100)
    p_GeV     : muon momentum in GeV/c
    """
    x_X0  = np.asarray(depth_mwe) * 100.0 / X0_ROCK_GCM2
    p_MeV = np.asarray(p_GeV) * 1000.0
    t_rad = (13.6 / (beta * p_MeV)) * np.sqrt(x_X0) * (1.0 + 0.038 * np.log(x_X0))
    return np.degrees(t_rad)

# ── Column alias resolver ─────────────────────────────────────────────────────
ALIASES = {
    "EventID":        ["eventid", "event", "evt"],
    "DepthCm":        ["depthcm", "depth_cm", "depth"],
    "MWE":            ["mwe"],
    "PDG":            ["pdg", "pdgid"],
    "InitKE_GeV":     ["initke_gev", "initkegev", "initke", "ke_in"],
    "InitThetaDeg":   ["initthetadeg", "inittheta", "theta_in"],
    "InitPhiDeg":     ["initphideg",  "initphi",   "phi_in"],
    "ExitKE_GeV":     ["exitke_gev",  "exitkegev", "exitke", "ke_out"],
    "ExitPx_GeVc":    ["exitpx_gevc", "exitpxgevc", "exitpx"],
    "ExitPy_GeVc":    ["exitpy_gevc", "exitpygevc", "exitpy"],
    "ExitPz_GeVc":    ["exitpz_gevc", "exitpzgevc", "exitpz"],
    "ExitX_cm":       ["exitx_cm",    "exitxcm",   "exitx"],
    "ExitY_cm":       ["exity_cm",    "exitycm",   "exity"],
    "ELoss_Ion_GeV":  ["eloss_ion_gev",  "elossiongev"],
    "ELoss_Brem_GeV": ["eloss_brem_gev", "elossbremgev"],
    "ELoss_Pair_GeV": ["eloss_pair_gev", "elosspairgev"],
    "ELoss_Nucl_GeV": ["eloss_nucl_gev", "elossnuclgev"],
    "ELoss_Total_GeV":["eloss_total_gev","elosstotalgev","elosstotal","eloss_gev"],
    "AngleScat_Deg":  ["anglescat_deg",  "anglescatdeg", "angle_deg"],
    "LatDisp_cm":     ["latdisp_cm",     "latdispcm"],
}

def resolve_cols(df):
    col_map = {c.lower().replace("_","").replace(" ",""): c for c in df.columns}
    rename  = {}
    for canon, variants in ALIASES.items():
        if canon in df.columns:
            continue
        for v in variants:
            k = v.lower().replace("_","").replace(" ","")
            if k in col_map:
                rename[col_map[k]] = canon
                break
    return df.rename(columns=rename) if rename else df

def derive_cols(df):
    """Add derived columns if missing."""
    if "InitKE_GeV" in df.columns and "InitP_GeVc" not in df.columns:
        df = df.copy()
        df["InitP_GeVc"] = np.sqrt((df["InitKE_GeV"] + M_MU)**2 - M_MU**2)
    if "ExitKE_GeV" in df.columns and "ExitP_GeVc" not in df.columns:
        df["ExitP_GeVc"] = np.sqrt(np.maximum((df["ExitKE_GeV"] + M_MU)**2 - M_MU**2, 0))
    if "ELoss_Total_GeV" not in df.columns:
        procs = [c for c in ["ELoss_Ion_GeV","ELoss_Brem_GeV",
                              "ELoss_Pair_GeV","ELoss_Nucl_GeV"] if c in df.columns]
        if procs:
            df["ELoss_Total_GeV"] = df[procs].sum(axis=1)
        elif "InitKE_GeV" in df.columns and "ExitKE_GeV" in df.columns:
            df["ELoss_Total_GeV"] = df["InitKE_GeV"] - df["ExitKE_GeV"]
    if "MWE" not in df.columns and "DepthCm" in df.columns:
        df["MWE"] = df["DepthCm"] * ROCK_DENSITY / 100.0
    return df

# ── CSV loader ─────────────────────────────────────────────────────────────────
def load_csv(path, tag=""):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.readlines()
    hdr = next((i for i, l in enumerate(lines)
                if l.strip() and not l.startswith("#")), 0)
    df = pd.read_csv(path, skiprows=hdr, sep=r"[\s,]+", engine="python")
    df.columns = df.columns.str.strip()
    if tag:
        print(f"  {tag}: {len(df):,} rows  cols={list(df.columns)}")
    return resolve_cols(df)

def parse_summary(path):
    """Extract n_input, timing_s, rate_evts from a RunAction *_summary.txt."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            m = re.search(r'Events\s*:\s*(\d+)', line)
            if m: out["n_input"] = int(m.group(1))
            m = re.search(r'Elapsed\s*:\s*([\d.]+)', line)
            if m: out["timing_s"] = float(m.group(1))
            m = re.search(r'Rate\s*:\s*([\d.]+)\s*evt', line)
            if m: out["rate"] = float(m.group(1))
    return out

# ── Per-depth statistics ──────────────────────────────────────────────────────
def compute_stats(muons, n_input):
    rows = []
    for d in sorted(muons["DepthCm"].unique()):
        if d in EXCLUDE_DEPTH_CM:
            continue
        sub = muons[muons["DepthCm"] == d]
        n   = len(sub)
        mwe = d * ROCK_DENSITY / 100.0
        row = {"DepthCm": int(d), "MWE": round(mwe, 3),
               "N_transmitted": n,
               "Transmission_%": round(100.0 * n / max(n_input, 1), 3)}
        for col, key in [("InitKE_GeV",     "MeanInitKE_GeV"),
                          ("ExitKE_GeV",     "MeanExitKE_GeV"),
                          ("ELoss_Total_GeV","MeanELoss_GeV"),
                          ("AngleScat_Deg",  "MeanAngle_deg"),
                          ("LatDisp_cm",     "MeanLatDisp_cm")]:
            if col in sub.columns:
                row[key] = round(sub[col].mean(), 4)
        for col, key in [("ELoss_Total_GeV","StdELoss_GeV"),
                          ("AngleScat_Deg",  "StdAngle_deg")]:
            if col in sub.columns:
                row[key] = round(sub[col].std(), 4)
        if "AngleScat_Deg" in sub.columns:
            row["P95Angle_deg"] = round(sub["AngleScat_Deg"].quantile(0.95), 4)
            if n > 1 and "StdAngle_deg" in row:
                row["SEM_Angle_deg"] = round(row["StdAngle_deg"] / np.sqrt(n), 6)
        procs = ["ELoss_Ion_GeV","ELoss_Brem_GeV","ELoss_Pair_GeV","ELoss_Nucl_GeV"]
        if all(c in sub.columns for c in procs):
            # Use sum of per-process columns as denominator so fractions sum to 100%
            # regardless of whether ELoss_Total_GeV was overridden to the kinematic value.
            proc_sum = max(sub[procs].sum().sum(), 1e-12)
            row["FracIon_%"]  = round(100*sub["ELoss_Ion_GeV"].sum()  / proc_sum, 2)
            row["FracBrem_%"] = round(100*sub["ELoss_Brem_GeV"].sum() / proc_sum, 2)
            row["FracPair_%"] = round(100*sub["ELoss_Pair_GeV"].sum() / proc_sum, 2)
            row["FracNucl_%"] = round(100*sub["ELoss_Nucl_GeV"].sum() / proc_sum, 2)
        rows.append(row)
    return pd.DataFrame(rows)

# ═════════════════════════════════════════════════════════════════════════════
#  CODE LOADERS
# ═════════════════════════════════════════════════════════════════════════════

def load_geant4(run_dir, n_override=None):
    """
    Load Geant4 MuonRock output.  Accepts three forms:
      1. Flat outputs/ directory  ->  picks the latest *_muons.csv inside it
      2. File prefix  (outputs/run_file_20260512_180248)
         ->  resolves to outputs/run_file_20260512_180248_muons.csv
      3. A subdirectory containing *_muons.csv files (legacy layout)
    """
    print(f"\n[Geant4] Loading: {run_dir}")

    # Resolve to the common file base (path without _muons.csv suffix)
    if os.path.isdir(run_dir):
        # Directory: find latest *_muons.csv directly inside it
        candidates = glob.glob(os.path.join(run_dir, "*_muons.csv"))
        if not candidates:
            # Legacy: look in run_* subdirectories
            for sub in glob.glob(os.path.join(run_dir, "run_*")):
                if os.path.isdir(sub):
                    candidates.extend(glob.glob(os.path.join(sub, "*_muons.csv")))
        if not candidates:
            print(f"  [ERROR] No *_muons.csv found in {run_dir}")
            return None
        # Pick most recently modified (not alphabetically last)
        muon_path = max(candidates, key=os.path.getmtime)
        base = muon_path[:-len("_muons.csv")]
        print(f"  Auto-selected  : {os.path.basename(base)}")
    else:
        # File prefix: try appending _muons.csv directly
        muon_path = run_dir + "_muons.csv"
        if not os.path.exists(muon_path):
            parent  = os.path.dirname(run_dir)
            matches = sorted(glob.glob(
                os.path.join(parent, os.path.basename(run_dir) + "*_muons.csv")))
            if not matches:
                print(f"  [ERROR] No muon CSV found for prefix: {run_dir}")
                return None
            muon_path = matches[-1]
        base = muon_path[:-len("_muons.csv")]

    sec_path  = base + "_secondaries.csv"
    stop_path = base + "_stopped.csv"
    summ_path = base + "_summary.txt"

    muons = load_csv(muon_path, os.path.basename(muon_path))
    if muons is None or muons.empty:
        print("  [ERROR] muon CSV is empty")
        return None
    muons = derive_cols(muons)

    # Recompute AngleScat_Deg = true 3-D MCS deflection from initial direction.
    # RunAction.cc writes the absolute exit zenith angle (angle from vertical), which
    # is dominated by the source zenith spread and is NOT the MCS scatter angle.
    # Here we recompute from exit momentum vs. initial direction so that Geant4 is
    # directly comparable to MUSIC (which also reports MCS deflection).
    if all(c in muons.columns for c in
           ["InitThetaDeg", "InitPhiDeg", "ExitPx_GeVc", "ExitPy_GeVc", "ExitPz_GeVc"]):
        th_i = np.radians(muons["InitThetaDeg"])
        ph_i = np.radians(muons["InitPhiDeg"])
        cx_i =  np.sin(th_i) * np.cos(ph_i)
        cy_i =  np.sin(th_i) * np.sin(ph_i)
        cz_i = -np.cos(th_i)                           # downward = −z
        ptot = np.sqrt(muons["ExitPx_GeVc"]**2
                     + muons["ExitPy_GeVc"]**2
                     + muons["ExitPz_GeVc"]**2).clip(lower=1e-12)
        dot  = (cx_i * muons["ExitPx_GeVc"] / ptot
              + cy_i * muons["ExitPy_GeVc"] / ptot
              + cz_i * muons["ExitPz_GeVc"] / ptot).clip(-1.0, 1.0)
        muons["AngleScat_Deg"] = np.degrees(np.arccos(dot))
        print("  AngleScat_Deg  : recomputed as MCS deflection from initial direction")

    # LatDisp_cm stored by RunAction includes the source-position spread (muons are
    # injected at various (xs, ys) across the surface sampling area, giving values
    # ~1800 cm even at 1 m depth).  Without surface entry positions in the CSV we
    # cannot separate MCS transverse spread from geometric drift, so we set LatDisp_cm
    # to NaN.  Exit (x, y) columns are preserved for the position map (fig11).
    if "LatDisp_cm" in muons.columns:
        muons["LatDisp_cm"] = float("nan")
        print("  LatDisp_cm     : set to NaN (source position not in CSV; see fig11)")

    # ELoss_Total_GeV from the CSV = sum of GetTotalEnergyDeposit() over muon steps.
    # This captures only local muon-step deposits; energy carried away by secondary
    # photons / e+e- pairs is NOT included.  At shallow depths ~67% of the true
    # kinematic loss is captured; at 200 m only ~53% (radiative processes dominate).
    # Replace with the kinematic definition (InitKE − ExitKE) so that Geant4 is
    # directly comparable to the grid codes (MUSIC/BB/UCMuon/PROPOSAL) which compute
    # ELoss from the energy difference.  Keep the original per-process columns
    # (ELoss_Ion_GeV etc.) unchanged — the process fractions in fig02 use their sum
    # as their own denominator and remain internally consistent.
    if "InitKE_GeV" in muons.columns and "ExitKE_GeV" in muons.columns:
        muons["ELoss_Total_GeV"] = muons["InitKE_GeV"] - muons["ExitKE_GeV"]
        print("  ELoss_Total_GeV: replaced with InitKE−ExitKE (kinematic; comparable to grid codes)")

    sec  = load_csv(sec_path  if os.path.exists(sec_path)  else None, os.path.basename(sec_path))
    stop = load_csv(stop_path if os.path.exists(stop_path) else None, os.path.basename(stop_path))
    meta = parse_summary(summ_path)

    n_input  = n_override or meta.get("n_input", muons["EventID"].nunique())
    timing_s = meta.get("timing_s", 0.0)

    depths = sorted(muons["DepthCm"].unique())
    print(f"  Injected events: {n_input:,}")
    print(f"  Muon rows      : {len(muons):,}")
    print(f"  Depths (cm)    : {[int(d) for d in depths]}")
    if timing_s:
        print(f"  Timing         : {timing_s:.1f} s  ({meta.get('rate',0):.0f} evt/s)")

    return {"label": "Geant4", "color": CODE_COLORS["Geant4"],
            "muons": muons, "secondaries": sec, "stopped": stop,
            "n_input": n_input, "timing_s": timing_s,
            "stats": compute_stats(muons, n_input),
            "base": base}   # e.g. ../outputs/run_file_20260512_180248


def load_phits(run_dir, n_override=None):
    """
    Load PHITS output.

    Provide a CSV normalized to the common schema:
      File: phits_muons.csv
      Required columns: EventID, DepthCm, PDG, InitKE_GeV, ExitKE_GeV
      Optional columns: AngleScat_Deg, LatDisp_cm

    OR a summary-only file:
      File: phits_summary.csv
      Columns: DepthCm, MWE, N_transmitted, Transmission_%, MeanExitKE_GeV

    Include timing (same format as Geant4 summary.txt) in phits_timing.txt.

    PHITS geometry / physics notes (muon_rock.inp reviewed 2026-05-23):
      - Rock: Mg(75%)+O(25%) mass fraction, density 2.65 g/cm³ → Z_eff=11, A_eff=22.23
        Geant4 uses Z=11, A=22 single-element → I value differs by ~8% (149 vs 138 eV)
        Effect on ionisation dE/dx < 0.3% — not the cause of the −12% KE discrepancy.
      - Scoring depths match Geant4: 100, 1000, 2500, 5000, 10000, 20000 cm.
      - emin(6)=1 MeV (µ+), emin(7)=1 MeV (µ−) — fixed 2026-05-23 (emin(5) is π−).
      - MCS: Lynch-Molière (nspred=-2) vs Geant4 Urban MscModel — model difference,
        not a config error; explains part of the ~8% MCS angle offset vs MUSIC.
      - No per-event output: all PHITS tallies are aggregate histograms.
        AngleScat_Deg and LatDisp_cm will always be NaN for PHITS results.
      - dE/dx discrepancy (−12.2% exit KE at 530 MWE): likely different cross-section
        parametrisations for muon brem/pair at 100–300 GeV. Run single-layer test to confirm.
    """
    print(f"\n[PHITS] Loading: {run_dir}")
    muon_path = os.path.join(run_dir, "phits_muons.csv")
    summ_path = os.path.join(run_dir, "phits_summary.csv")
    time_path = os.path.join(run_dir, "phits_timing.txt")

    if os.path.exists(muon_path):
        muons = load_csv(muon_path, "phits_muons.csv")
        if muons is None: return None
        muons    = derive_cols(muons)
        n_input  = n_override or muons["EventID"].nunique()
        timing_s = parse_summary(time_path).get("timing_s", 0.0)
        stats    = compute_stats(muons, n_input)
        return {"label": "PHITS", "color": CODE_COLORS["PHITS"],
                "muons": muons, "secondaries": None, "stopped": None,
                "n_input": n_input, "timing_s": timing_s, "stats": stats}

    if os.path.exists(summ_path):
        stats    = pd.read_csv(summ_path)
        timing_s = parse_summary(time_path).get("timing_s", 0.0)
        n_input  = n_override or 0
        if not n_input and "Transmission_%" in stats.columns and "N_transmitted" in stats.columns:
            tr = stats.iloc[0]["Transmission_%"]
            n  = stats.iloc[0]["N_transmitted"]
            if tr > 0: n_input = int(round(n / (tr / 100.0)))
        print(f"  Summary-only mode ({len(stats)} depths); no per-event data")
        # Mark as old/mismatched run only if fewer than 50K events were used
        old_run = (n_input > 0 and n_input < 50000)
        label   = "PHITS*" if old_run else "PHITS"
        result  = {"label": label, "color": CODE_COLORS["PHITS"],
                   "muons": None, "secondaries": None, "stopped": None,
                   "n_input": n_input, "timing_s": timing_s, "stats": stats}
        if old_run:
            result["note"] = "old run / mismatched source — re-run for apples-to-apples"
        return result

    print(f"  [SKIP] No PHITS output found — need phits_muons.csv or phits_summary.csv")
    return None


# ── Per-plane .dat loader (shared by MUSIC and PROPOSAL) ─────────────────────
_DEPTH_FROM_FNAME = {
    "1m": 100, "10m": 1000, "25m": 2500,
    "50m": 5000, "100m": 10000, "200m": 20000,
}

_DAT_COLS = ["EventID","xs","ys","zs","Es","ts","ps","charge","alive",
             "x","y","z","E","cx","cy","cz","theta_ug","phi_ug"]

def _load_dat_files(run_dir, prefix, n_override):
    """
    Read per-plane ASCII .dat files written by MUSIC or PROPOSAL.

    File naming: <PREFIX>_<depth>.dat  (e.g. MUSIC_10m.dat, PROPOSAL_25m.dat)
    Column layout (18 columns, same for both codes):
        EventID xs ys zs Es theta_srf phi_srf charge alive
        x y z E cx cy cz theta_ug phi_ug

    Returns (muons_df, n_input) where muons_df has standard benchmark columns:
        EventID, DepthCm, PDG, InitKE_GeV, ExitKE_GeV,
        ExitX_cm, ExitY_cm, AngleScat_Deg, LatDisp_cm
    Only alive==1 rows are kept; n_input is the source particle count.
    """
    files = sorted(glob.glob(os.path.join(run_dir, f"{prefix}_*.dat")))
    if not files:
        return None, 0

    frames  = []
    n_input = 0
    for fpath in files:
        tag = os.path.basename(fpath).replace(f"{prefix}_", "").replace(".dat", "")
        tag = tag.replace("bench_", "")   # handle <PREFIX>_bench_<depth>.dat naming
        depth_cm = _DEPTH_FROM_FNAME.get(tag)
        if depth_cm is None:
            continue

        df = pd.read_csv(fpath, sep=r"\s+", comment="#", header=None,
                         names=_DAT_COLS, engine="python")
        if n_input == 0:
            n_input = len(df)

        alive = df[df["alive"] == 1].copy()
        if alive.empty:
            continue

        alive["DepthCm"]    = depth_cm
        alive["PDG"]        = np.where(alive["charge"] == 1, -13, 13)
        # Es/E are total energies in the .dat format; subtract muon mass to get KE
        alive["InitKE_GeV"] = alive["Es"] - M_MU
        alive["ExitKE_GeV"] = np.maximum(alive["E"] - M_MU, 0.0)
        alive["EntryX_cm"]  = alive["xs"]
        alive["EntryY_cm"]  = alive["ys"]
        alive["ExitX_cm"]   = alive["x"]
        alive["ExitY_cm"]   = alive["y"]

        # Direction cosines at surface from (theta_srf=ts, phi_srf=ps)
        cx_s = np.sin(alive["ts"]) * np.cos(alive["ps"])
        cy_s = np.sin(alive["ts"]) * np.sin(alive["ps"])
        cz_s = -np.cos(alive["ts"])   # downward = −z direction

        # Detect no-MCS mode: if the exit direction cosines (cx, cy) have zero
        # variance, every muon exited with the same direction → scattering disabled.
        # This is geometry-independent and works for any source direction, unlike
        # comparing theta_ug to theta_srf (fragile for non-vertical sources).
        cx_std = alive["cx"].std() if len(alive) > 1 else 0.0
        cy_std = alive["cy"].std() if len(alive) > 1 else 0.0
        no_mcs = (cx_std < 1e-6) and (cy_std < 1e-6)
        if no_mcs:
            print(f"    WARNING: {os.path.basename(fpath)} — cx=cy=const for all muons"
                  f" (transport ran without MCS; direction/scatter data unavailable)")
            alive["AngleScat_Deg"] = float("nan")
            alive["LatDisp_cm"]    = float("nan")
        else:
            dot = (cx_s * alive["cx"] + cy_s * alive["cy"]
                   + cz_s * alive["cz"]).clip(-1.0, 1.0)
            alive["AngleScat_Deg"] = np.degrees(np.arccos(dot))
            path  = np.abs(alive["z"]) / np.abs(cz_s).clip(lower=1e-6)
            exp_x = alive["xs"] + cx_s * path
            exp_y = alive["ys"] + cy_s * path
            alive["LatDisp_cm"] = np.sqrt((alive["x"] - exp_x)**2
                                          + (alive["y"] - exp_y)**2)

        keep = ["EventID","DepthCm","PDG","InitKE_GeV","ExitKE_GeV",
                "EntryX_cm","EntryY_cm",
                "ExitX_cm","ExitY_cm","AngleScat_Deg","LatDisp_cm"]
        frames.append(alive[keep].reset_index(drop=True))
        print(f"    {os.path.basename(fpath)}: {len(alive):,} surviving / {len(df):,} total")

    if not frames:
        return None, n_input
    return pd.concat(frames, ignore_index=True), n_input


def load_music(run_dir, n_override=None):
    """
    Load MUSIC output from per-plane MUSIC_<depth>.dat files.

    Expected files in run_dir: MUSIC_1m.dat, MUSIC_10m.dat, MUSIC_25m.dat,
    MUSIC_50m.dat, MUSIC_100m.dat, MUSIC_200m.dat.

    Columns per file (18, space-separated):
        EventID xs ys zs Es theta_srf phi_srf charge alive
        x y z E cx cy cz theta_ug phi_ug
    where (xs,ys,zs) is surface entry, (x,y,z) is underground position,
    alive=1 means the muon traversed the full depth.

    Provides: InitKE_GeV, ExitKE_GeV, ExitX_cm, ExitY_cm,
              AngleScat_Deg (3-D scatter from initial direction),
              LatDisp_cm (transverse displacement from straight-line path).
    """
    print(f"\n[MUSIC] Loading: {run_dir}")
    time_path = os.path.join(run_dir, "music_timing.txt")

    muons, n_raw = _load_dat_files(run_dir, "MUSIC", n_override)
    if muons is None:
        print(f"  [SKIP] No MUSIC_*.dat files found in {run_dir}")
        return None

    muons    = derive_cols(muons)
    n_input  = n_override or n_raw
    timing_s = parse_summary(time_path).get("timing_s", 0.0)
    stats    = compute_stats(muons, n_input)
    print(f"  n_input={n_input:,}  surviving rows={len(muons):,}")
    return {"label": "MUSIC", "color": CODE_COLORS["MUSIC"],
            "muons": muons, "secondaries": None, "stopped": None,
            "n_input": n_input, "timing_s": timing_s, "stats": stats}


def load_proposal(run_dir, n_override=None):
    """
    Load PROPOSAL output from per-plane PROPOSAL_<depth>.dat files.

    Expected files in run_dir: PROPOSAL_1m.dat, PROPOSAL_10m.dat, etc.
    Column layout is identical to MUSIC (18 columns).

    Provides: InitKE_GeV, ExitKE_GeV, ExitX_cm, ExitY_cm,
              AngleScat_Deg, LatDisp_cm.
    """
    print(f"\n[PROPOSAL] Loading: {run_dir}")
    time_path = os.path.join(run_dir, "proposal_timing.txt")

    muons, n_raw = _load_dat_files(run_dir, "PROPOSAL", n_override)
    if muons is None:
        print(f"  [SKIP] No PROPOSAL_*.dat files found in {run_dir}")
        return None

    muons    = derive_cols(muons)
    n_input  = n_override or n_raw
    timing_s = parse_summary(time_path).get("timing_s", 0.0)
    stats    = compute_stats(muons, n_input)
    print(f"  n_input={n_input:,}  surviving rows={len(muons):,}")
    return {"label": "PROPOSAL", "color": CODE_COLORS["PROPOSAL"],
            "muons": muons, "secondaries": None, "stopped": None,
            "n_input": n_input, "timing_s": timing_s, "stats": stats}


def load_bb(run_dir, n_override=None):
    """
    Load BB engine output from per-plane BB_bench_<depth>.dat files.
    Column layout is identical to MUSIC/PROPOSAL (18 columns).
    """
    print(f"\n[BB] Loading: {run_dir}")
    time_path = os.path.join(run_dir, "bb_timing.txt")

    muons, n_raw = _load_dat_files(run_dir, "BB", n_override)
    if muons is None:
        print(f"  [SKIP] No BB_*.dat files found in {run_dir}")
        return None

    muons    = derive_cols(muons)
    n_input  = n_override or n_raw
    timing_s = parse_summary(time_path).get("timing_s", 0.0)
    stats    = compute_stats(muons, n_input)
    print(f"  n_input={n_input:,}  surviving rows={len(muons):,}")
    return {"label": "BB", "color": CODE_COLORS["BB"],
            "muons": muons, "secondaries": None, "stopped": None,
            "n_input": n_input, "timing_s": timing_s, "stats": stats}


def load_ucmuon(run_dir, n_override=None):
    """
    Load UCMuon engine output from per-plane UCMuon_bench_<depth>.dat files.
    Column layout is identical to MUSIC/PROPOSAL (18 columns).
    """
    print(f"\n[UCMuon] Loading: {run_dir}")
    time_path = os.path.join(run_dir, "ucmuon_timing.txt")

    muons, n_raw = _load_dat_files(run_dir, "UCMuon", n_override)
    if muons is None:
        print(f"  [SKIP] No UCMuon_*.dat files found in {run_dir}")
        return None

    muons    = derive_cols(muons)
    n_input  = n_override or n_raw
    timing_s = parse_summary(time_path).get("timing_s", 0.0)
    stats    = compute_stats(muons, n_input)
    print(f"  n_input={n_input:,}  surviving rows={len(muons):,}")
    return {"label": "UCMuon", "color": CODE_COLORS["UCMuon"],
            "muons": muons, "secondaries": None, "stopped": None,
            "n_input": n_input, "timing_s": timing_s, "stats": stats}

# ═════════════════════════════════════════════════════════════════════════════
#  FIGURE FUNCTIONS
# ═════════════════════════════════════════════════════════════════════════════

def _cmap_depths(n):
    return plt.cm.viridis(np.linspace(0.1, 0.9, max(n, 1)))


# ── Fig 01: Transmission ──────────────────────────────────────────────────────
def fig_transmission(results, outdir):
    # Two-panel layout: full linear range (top) + zoomed view of deep depths (bottom)
    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(9, 9),
                                            gridspec_kw={"hspace": 0.45})
    has_old_src = False
    n_col = "N_transmitted"

    for i, res in enumerate(results):
        s = res["stats"]
        if "Transmission_%" not in s.columns: continue

        s_plot = s[s["Transmission_%"] > 0].copy()
        if s_plot.empty: continue

        kw = line_kw(i, res["color"])
        label = res["label"]

        for ax in (ax_full, ax_zoom):
            ax.plot(s_plot["MWE"], s_plot["Transmission_%"],
                    label=label, **kw)
            label = "_nolegend_"  # only add legend entry once

            if n_col in s_plot.columns:
                s_low = s_plot[s_plot[n_col] < 5]
                if not s_low.empty:
                    ax.plot(s_low["MWE"], s_low["Transmission_%"],
                            color=res["color"], marker=MARKERS[i % 4],
                            markersize=11, linestyle="none",
                            fillstyle="none", markeredgewidth=1.8,
                            alpha=0.8, zorder=5)

        # N annotations for Geant4 only on the full panel
        if i == 0 and n_col in s_plot.columns:
            for _, row in s_plot.iterrows():
                ax_full.annotate(f"N={int(row[n_col]):,}",
                                 (row["MWE"], row["Transmission_%"]),
                                 textcoords="offset points", xytext=(0, 8),
                                 ha="center", fontsize=7, color="#555")

        if "note" in res:
            has_old_src = True

    # Full-range panel: linear 0–105 %
    ax_full.set_xlabel("Depth (m.w.e.)", fontsize=13)
    ax_full.set_ylabel("Muon transmission (%)", fontsize=13)
    ax_full.set_title("Muon Survival vs Rock Depth  (linear)", fontsize=14)
    ax_full.set_ylim(0, 108)
    ax_full.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax_full.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax_full.legend(fontsize=10)
    ax_full.grid(True, alpha=0.3)

    # Zoomed panel: deep depths only (≥ 66.25 mwe), spread of codes is visible
    ax_zoom.set_xlabel("Depth (m.w.e.)", fontsize=13)
    ax_zoom.set_ylabel("Muon transmission (%)", fontsize=13)
    ax_zoom.set_title("Zoomed: deep depths  (≥ 66.25 m.w.e.)", fontsize=13)
    ax_zoom.set_xlim(55, None)
    ax_zoom.set_ylim(0, 108)
    ax_zoom.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax_zoom.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax_zoom.grid(True, alpha=0.3)

    note = "PHITS*: old 10K run on previous source  |  " if has_old_src else ""
    ax_zoom.text(0.98, 0.02,
                 note + "open rings = N < 5 (low statistics, treat with caution)",
                 transform=ax_zoom.transAxes, ha="right", va="bottom",
                 fontsize=8, color="#777", style="italic")
    _save(fig, outdir, "fig01_transmission.png")


# ── Fig 02: Energy loss ────────────────────────────────────────────────────────
def fig_energy_loss(results, outdir):
    fig = plt.figure(figsize=(16, 5))
    gs  = GridSpec(1, 3, figure=fig, wspace=0.40)

    ax1 = fig.add_subplot(gs[0])
    for i, res in enumerate(results):
        s = res["stats"]
        if "MeanELoss_GeV" not in s.columns: continue
        ax1.errorbar(s["MWE"], s["MeanELoss_GeV"],
                     yerr=s.get("StdELoss_GeV"),
                     fmt="o", color=res["color"], linewidth=2, capsize=3,
                     linestyle=LINE_STYLES[i % 4], label=res["label"])
    ax1.set_xlabel("Depth (m.w.e.)", fontsize=12)
    ax1.set_ylabel("Mean energy loss (GeV)", fontsize=12)
    ax1.set_title("Mean ΔE vs Depth", fontsize=13)
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

    # Process fractions — Geant4 only (index-based x so bars are always visible)
    ax2 = fig.add_subplot(gs[1])
    g4 = _get(results, "Geant4")
    if g4 is not None:
        s = g4["stats"]
        fracs = [("Ionisation",    "FracIon_%",  "#2196F3"),
                 ("Bremsstrahlung","FracBrem_%",  "#E91E63"),
                 ("Pair prod.",    "FracPair_%",  "#FF9800"),
                 ("Nuclear",       "FracNucl_%",  "#9C27B0")]
        avail = [(p, c, col) for p, c, col in fracs if c in s.columns]
        if avail:
            x   = np.arange(len(s))
            bot = np.zeros(len(s))
            for proc, col, color in avail:
                ax2.bar(x, s[col].values, bottom=bot,
                        color=color, label=proc, width=0.72, edgecolor="white", linewidth=0.4)
                bot += s[col].values
            ax2.set_xticks(x)
            ax2.set_xticklabels([f"{mwe:.1f}" for mwe in s["MWE"]],
                                rotation=40, ha="right", fontsize=9)
            ax2.set_xlabel("Depth (m.w.e.)")
            ax2.set_ylabel("Fraction of ΔE (%)")
            ax2.set_title("ΔE by Process (Geant4)")
            ax2.legend(fontsize=8, loc="upper right"); ax2.set_ylim(0, 115)
            ax2.grid(True, alpha=0.3, axis="y")

    # ΔE distributions — log y, step histograms, skip sparse depths
    ax3 = fig.add_subplot(gs[2])
    first = _first_with(results, "ELoss_Total_GeV")
    if first:
        m, depths = first["muons"], sorted(first["muons"]["DepthCm"].unique())
        cmap = _cmap_depths(len(depths))
        plotted = 0
        for i, d in enumerate(depths):
            sub = m[m["DepthCm"] == d]["ELoss_Total_GeV"].dropna()
            if len(sub) < 10:
                continue
            mwe_d = d * ROCK_DENSITY / 100.0
            ax3.hist(sub, bins=50, color=cmap[i], linewidth=1.8,
                     label=f"{mwe_d:.1f} m.w.e.", density=True, histtype="step")
            plotted += 1
        if plotted > 0:
            ax3.set_yscale("log")
            ax3.set_xlabel("Energy loss (GeV)")
            ax3.set_ylabel("Probability density")
            ax3.set_title(f"ΔE Distribution ({first['label']})")
            ax3.legend(fontsize=7, ncol=2)
            ax3.grid(True, alpha=0.3, which="both")

    _save(fig, outdir, "fig02_energy_loss.png")


# ── Fig 03: Angular scattering + Highland ─────────────────────────────────────
def fig_angular(results, outdir):
    fig = plt.figure(figsize=(13, 5))
    gs  = GridSpec(1, 2, figure=fig, wspace=0.35)

    ax1 = fig.add_subplot(gs[0])
    for i, res in enumerate(results):
        s = res["stats"]
        if "MeanAngle_deg" not in s.columns: continue
        # Skip codes where all angle values are NaN (PROPOSAL, PHITS)
        valid = s[s["MeanAngle_deg"].notna()]
        if valid.empty: continue
        sem_col = "SEM_Angle_deg" if "SEM_Angle_deg" in s.columns else None
        valid_yerr = s.loc[s["MeanAngle_deg"].notna(), sem_col] if sem_col else None
        ax1.errorbar(valid["MWE"], valid["MeanAngle_deg"],
                     yerr=valid_yerr,
                     fmt="o", color=res["color"], linewidth=2, capsize=3,
                     linestyle=LINE_STYLES[i % 4], label=res["label"])

    # Highland theory reference (uses mean momentum of surviving Geant4 muons)
    g4 = _get(results, "Geant4")
    if g4 is not None and g4["muons"] is not None and "InitP_GeVc" in g4["muons"].columns:
        mean_p = g4["muons"]["InitP_GeVc"].mean()
        s = g4["stats"]
        mwe_fine = np.linspace(max(s["MWE"].min() * 0.05, 0.01), s["MWE"].max() * 1.02, 300)
        ax1.plot(mwe_fine, highland_deg(mwe_fine, mean_p),
                 "k--", linewidth=1.5, alpha=0.7,
                 label=fr"Highland ($\langle p\rangle_{{G4}}$={mean_p:.0f} GeV/c)")

    ax1.set_xlabel("Depth (m.w.e.)", fontsize=12)
    ax1.set_ylabel("Mean MCS deflection angle (°)", fontsize=12)
    ax1.set_title("MCS Scattering Angle vs Depth\n"
                  r"(3-D deflection from each muon's initial direction)", fontsize=12)
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)
    ax1.text(0.98, 0.97,
             "Error bars: ±SEM  |  26.5 MWE excluded (5 GeV Bragg boundary)",
             transform=ax1.transAxes, ha="right", va="top",
             fontsize=8, color="#777", style="italic")

    # Distribution panel: show first code that has real (non-NaN) angle data
    first = next((r for r in results
                  if r["muons"] is not None
                  and "AngleScat_Deg" in r["muons"].columns
                  and not r["muons"]["AngleScat_Deg"].isna().all()), None)
    ax2 = fig.add_subplot(gs[1])
    if first:
        m, depths = first["muons"], sorted(first["muons"]["DepthCm"].unique())
        cmap = _cmap_depths(len(depths))
        plotted = 0
        for i, d in enumerate(depths):
            sub = m[m["DepthCm"] == d]["AngleScat_Deg"].dropna()
            if len(sub) < 10: continue
            mwe_d = d * ROCK_DENSITY / 100.0
            ax2.hist(sub, bins=50, color=cmap[i], linewidth=1.8,
                     label=f"{mwe_d:.1f} m.w.e.", density=True, histtype="step")
            plotted += 1
        if plotted > 0:
            ax2.set_yscale("log")
            ax2.set_xlabel("MCS deflection angle (°)")
            ax2.set_ylabel("Probability density")
            ax2.set_title(f"MCS Angle Distribution ({first['label']})")
            ax2.legend(fontsize=7, ncol=2)
            ax2.grid(True, alpha=0.3, which="both")

    _save(fig, outdir, "fig03_angular.png")


# ── Fig 04: Lateral displacement ──────────────────────────────────────────────
def fig_lateral(results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.40)

    # Left: MCS transverse displacement — only codes where LatDisp_cm is valid (non-NaN).
    # Geant4 LatDisp_cm is set to NaN because the source-position spread (~1000s cm)
    # dominates the exit position; MUSIC computes the true MCS-only transverse spread.
    plotted_l = 0
    for i, res in enumerate(results):
        s = res["stats"]
        if "MeanLatDisp_cm" not in s.columns: continue
        valid = s[s["MeanLatDisp_cm"].notna()]
        if valid.empty: continue
        axes[0].plot(valid["MWE"], valid["MeanLatDisp_cm"],
                     label=res["label"], **line_kw(i, res["color"]))
        plotted_l += 1
    axes[0].set_xlabel("Depth (m.w.e.)", fontsize=12)
    axes[0].set_ylabel("Mean MCS transverse displacement (cm)", fontsize=12)
    axes[0].set_title("MCS Lateral Displacement vs Depth\n"
                      "(deviation from straight-line trajectory)", fontsize=12)
    if plotted_l == 0:
        axes[0].text(0.5, 0.5, "No MCS lateral data available",
                     ha="center", va="center", transform=axes[0].transAxes,
                     fontsize=12, color="#777")
    axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)
    axes[0].text(0.98, 0.02,
                 "Geant4 & PROPOSAL excluded\n(source-position spread not in CSV)",
                 transform=axes[0].transAxes, ha="right", va="bottom",
                 fontsize=8, color="#777", style="italic")

    # Right: lateral displacement from entry position — codes that record both
    # entry (xs,ys) and exit (x,y) positions so comparisons are on equal footing.
    disp_codes = [r for r in results
                  if r["muons"] is not None
                  and "EntryX_cm" in r["muons"].columns
                  and "ExitX_cm"  in r["muons"].columns]
    if disp_codes:
        for res in disp_codes:
            m    = res["muons"]
            cidx = results.index(res)
            rows = []
            for d in sorted(m["DepthCm"].unique()):
                if d in EXCLUDE_DEPTH_CM:
                    continue
                sub = m[m["DepthCm"] == d]
                if len(sub) < 2: continue
                disp = np.sqrt((sub["ExitX_cm"] - sub["EntryX_cm"])**2
                               + (sub["ExitY_cm"] - sub["EntryY_cm"])**2).mean()
                rows.append({"MWE": d * ROCK_DENSITY / 100.0, "R": disp})
            if rows:
                df_r = pd.DataFrame(rows)
                axes[1].plot(df_r["MWE"], df_r["R"],
                             label=res["label"], **line_kw(cidx, res["color"]))
        axes[1].set_xlabel("Depth (m.w.e.)", fontsize=12)
        axes[1].set_ylabel(r"Mean $\sqrt{(\Delta x)^2+(\Delta y)^2}$ (cm)", fontsize=12)
        axes[1].set_title("Lateral Drift from Entry Position\n"
                          "(MCS + slant; comparable across codes)", fontsize=12)
        axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)
        axes[1].text(0.98, 0.02,
                     "Geant4/PHITS excluded\n(entry position not in output CSV)",
                     transform=axes[1].transAxes, ha="right", va="bottom",
                     fontsize=8, color="#777", style="italic")
    else:
        axes[1].text(0.5, 0.5,
                     "Entry position not available\n(dat-format codes only)",
                     ha="center", va="center", transform=axes[1].transAxes,
                     fontsize=11, color="#777")

    _save(fig, outdir, "fig04_lateral.png")


# ── Fig 05: Exit KE spectra at each depth ────────────────────────────────────
def fig_exit_spectrum(results, outdir):
    all_depths = set()
    for r in results:
        if r["muons"] is not None and "DepthCm" in r["muons"].columns:
            all_depths.update(r["muons"]["DepthCm"].unique())
    depths = sorted(d for d in all_depths if d not in EXCLUDE_DEPTH_CM)
    if not depths: return

    ncols = min(3, len(depths))
    nrows = (len(depths) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes_flat = np.array(axes).flatten()

    for k, d in enumerate(depths):
        ax  = axes_flat[k]
        mwe = d * ROCK_DENSITY / 100.0
        plotted = 0
        for res in results:
            m = res["muons"]
            if m is None or "ExitKE_GeV" not in m.columns: continue
            sub = m[m["DepthCm"] == d]["ExitKE_GeV"]
            if len(sub) < 2: continue
            ax.hist(sub, bins=60, alpha=0.7, color=res["color"],
                    label=res["label"], density=True,
                    histtype="step", linewidth=2)
            plotted += 1
        ax.set_title(f"Depth {int(d)} cm  ({mwe:.1f} m.w.e.)", fontsize=11)
        ax.set_xlabel("Exit KE (GeV)", fontsize=10)
        ax.set_ylabel("Probability density", fontsize=10)
        if plotted > 0:
            ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    for j in range(len(depths), len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle("Exit Kinetic Energy Spectra by Depth", fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, outdir, "fig05_exit_spectrum.png")


# ── Fig 06: dE/dx vs initial energy + Bethe-Bloch ────────────────────────────
def fig_dedx(results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.subplots_adjust(wspace=0.35)

    ax1 = axes[0]
    for i, res in enumerate(results):
        m = res["muons"]
        if m is None or "InitKE_GeV" not in m.columns: continue
        if "ELoss_Total_GeV" not in m.columns: continue
        depths = sorted(m["DepthCm"].unique())
        # Use shallowest depth so all energy groups are present
        d_use = depths[0]
        sub = m[m["DepthCm"] == d_use].copy()
        depth_mwe = d_use * ROCK_DENSITY / 100.0
        sub["dEdx"] = sub["ELoss_Total_GeV"] / depth_mwe  # GeV / m.w.e.
        try:
            lo = sub["InitKE_GeV"].quantile(0.02)
            hi = sub["InitKE_GeV"].quantile(0.98)
            if lo <= 0: lo = sub["InitKE_GeV"].min()
            if hi <= lo: hi = sub["InitKE_GeV"].max()
            if hi <= lo or np.isnan(lo) or np.isnan(hi):
                continue
            e_bins = np.logspace(np.log10(lo), np.log10(hi), 25)
            if len(np.unique(e_bins)) < 2:
                continue
            sub["E_bin"] = pd.cut(sub["InitKE_GeV"], bins=e_bins)
            grp = sub.groupby("E_bin", observed=True)["dEdx"].agg(["mean","std"]).dropna()
            e_ctr = [iv.mid for iv in grp.index]
            ax1.errorbar(e_ctr, grp["mean"], yerr=grp["std"],
                         fmt="o", color=res["color"], linewidth=2, capsize=3,
                         linestyle=LINE_STYLES[i % 4], label=res["label"])
        except Exception:
            pass

    e_th = np.logspace(0.5, np.log10(400), 200)
    ax1.plot(e_th, bethe_bloch_mwe(e_th), "k--", linewidth=1.5, alpha=0.7,
             label="Bethe-Bloch + radiative\n(Groom 2001)")
    ax1.set_xscale("log")
    ax1.set_xlabel("Initial KE (GeV)", fontsize=12)
    ax1.set_ylabel("Mean dE/dx  (GeV / m.w.e.)", fontsize=12)
    ax1.set_title("Energy Loss Rate vs Initial Energy\n(avg. over shallowest scoring plane)", fontsize=12)
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

    # Initial vs exit KE scatter (Geant4, sampled for speed)
    ax2 = axes[1]
    g4 = _get(results, "Geant4")
    if g4 is not None and g4["muons"] is not None:
        m = g4["muons"]
        if "InitKE_GeV" in m.columns and "ExitKE_GeV" in m.columns:
            depths = sorted(m["DepthCm"].unique())
            cmap   = _cmap_depths(len(depths))
            for k, d in enumerate(depths):
                s2 = m[m["DepthCm"] == d]
                s2 = s2.sample(min(len(s2), 3000), random_state=42)
                ax2.scatter(s2["InitKE_GeV"], s2["ExitKE_GeV"],
                            alpha=0.2, s=3, color=cmap[k], label=f"{int(d)} cm")
            emax = m["InitKE_GeV"].quantile(0.99)
            ax2.plot([0, emax], [0, emax], "k--", linewidth=1, alpha=0.4, label="no loss")
            ax2.set_xlabel("Initial KE (GeV)", fontsize=12)
            ax2.set_ylabel("Exit KE (GeV)", fontsize=12)
            ax2.set_title("Initial vs Exit KE (Geant4)", fontsize=13)
            ax2.legend(fontsize=7, ncol=2, markerscale=4)
            ax2.grid(True, alpha=0.3)

    _save(fig, outdir, "fig06_dedx.png")


# ── Fig 07: Muon charge ratio ─────────────────────────────────────────────────
def fig_charge_ratio(results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.35)

    ax1 = axes[0]
    for i, res in enumerate(results):
        m = res["muons"]
        if m is None or "PDG" not in m.columns: continue
        rows = []
        for d in sorted(m["DepthCm"].unique()):
            if d in EXCLUDE_DEPTH_CM:
                continue
            sub = m[m["DepthCm"] == d]
            nm  = (sub["PDG"] ==  13).sum()
            np_ = (sub["PDG"] == -13).sum()
            tot = nm + np_
            if tot > 0:
                rows.append({"MWE": d * ROCK_DENSITY / 100.0,
                             "ratio": nm / tot})
        if not rows: continue
        cr = pd.DataFrame(rows)
        ax1.plot(cr["MWE"], cr["ratio"], label=res["label"],
                 **line_kw(i, res["color"]))

    ax1.axhline(0.43, color="k", linestyle=":", linewidth=1.5, alpha=0.7,
                label=r"Surface $\mu^-$ fraction ($\approx$0.43)")
    ax1.set_ylim(0.2, 0.8)
    ax1.set_xlabel("Depth (m.w.e.)")
    ax1.set_ylabel(r"$\mu^-$ fraction  [$N(\mu^-)\,/\,N(\mu^-\!+\!\mu^+)$]")
    ax1.set_title("Muon Charge Ratio vs Depth")
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

    # µ⁻ vs µ⁺ initial KE at shallowest plane
    ax2 = axes[1]
    first = next((r for r in results if r["muons"] is not None
                  and "PDG" in r["muons"].columns
                  and "InitKE_GeV" in r["muons"].columns), None)
    if first:
        m  = first["muons"]
        d0 = sorted(m["DepthCm"].unique())[0]
        sub = m[m["DepthCm"] == d0]
        lo, hi = sub["InitKE_GeV"].min(), sub["InitKE_GeV"].quantile(0.99)
        bins = np.linspace(lo, hi, 60)
        ax2.hist(sub[sub["PDG"]==  13]["InitKE_GeV"], bins=bins,
                 alpha=0.6, color="#2196F3", label=r"$\mu^-$", density=True)
        ax2.hist(sub[sub["PDG"]== -13]["InitKE_GeV"], bins=bins,
                 alpha=0.6, color="#E91E63", label=r"$\mu^+$", density=True)
        ax2.set_xlabel("Initial KE (GeV)")
        ax2.set_ylabel("Probability density")
        ax2.set_title(fr"$\mu^-$ vs $\mu^+$ Spectrum ({first['label']}, {int(d0)} cm)")
        ax2.legend(fontsize=10); ax2.grid(True, alpha=0.3)

    _save(fig, outdir, "fig07_charge_ratio.png")


# ── Fig 08: Secondaries (Geant4 only) ─────────────────────────────────────────
def fig_secondaries(results, outdir):
    g4 = _get(results, "Geant4")
    if g4 is None or g4.get("secondaries") is None:
        print("Skipping fig08 — no secondary data")
        return
    sec = g4["secondaries"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.35)

    if "PDG" in sec.columns:
        counts = sec["PDG"].value_counts().head(10)
        labels = [PDG_NAMES.get(int(p), f"PDG {int(p)}") for p in counts.index]
        axes[0].bar(range(len(counts)), counts.values, color="#9C27B0", alpha=0.85)
        axes[0].set_xticks(range(len(counts)))
        axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
        axes[0].set_ylabel("Count", fontsize=12)
        axes[0].set_title(f"Secondary Particle Types  (total {len(sec):,})", fontsize=13)
        axes[0].grid(True, alpha=0.3, axis="y")
        for k, v in enumerate(counts.values):
            axes[0].text(k, v + counts.max() * 0.01, str(v), ha="center", fontsize=8)

    ke_col = next((c for c in ["KE_MeV","ke_mev","KE"] if c in sec.columns), None)
    if ke_col:
        clip = sec[ke_col].quantile(0.99)
        axes[1].hist(sec[ke_col].clip(upper=clip), bins=80,
                     color="#9C27B0", alpha=0.7)
        axes[1].set_xlabel("Kinetic energy (MeV)", fontsize=12)
        axes[1].set_ylabel("Count", fontsize=12)
        axes[1].set_title("Secondary KE Spectrum (99th pct clip)", fontsize=13)
        axes[1].set_yscale("log"); axes[1].grid(True, alpha=0.3)

    _save(fig, outdir, "fig08_secondaries.png")


# ── Fig 09: Stopped muons (Geant4 only) ───────────────────────────────────────
def fig_stopped(results, outdir):
    g4 = _get(results, "Geant4")
    if g4 is None or g4.get("stopped") is None:
        print("Skipping fig09 — no stopped-muon data (*_stopped.csv)")
        return
    stop     = g4["stopped"]
    n_input  = g4["n_input"]

    # Resolve columns
    stop = resolve_cols(stop)
    for can, variants in [("InitKE_GeV",  ["initke_gev","initkegev","initke"]),
                           ("StopDepth_cm",["stopdepthcm","stopdepth_cm","stopdepth"])]:
        if can not in stop.columns:
            sl = {c.lower().replace("_",""): c for c in stop.columns}
            for v in variants:
                k = v.lower().replace("_","")
                if k in sl:
                    stop = stop.rename(columns={sl[k]: can}); break

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.35)

    if "InitKE_GeV" in stop.columns:
        pct = 100.0 * len(stop) / max(n_input, 1)
        axes[0].hist(stop["InitKE_GeV"], bins=60, color="#F44336", alpha=0.85)
        axes[0].set_xlabel("Initial KE (GeV)", fontsize=12)
        axes[0].set_ylabel("Count", fontsize=12)
        axes[0].set_title(f"Stopped Muons — Initial KE\n"
                           f"({len(stop):,} / {n_input:,} = {pct:.1f}%)", fontsize=12)
        axes[0].grid(True, alpha=0.3)

    if "StopDepth_cm" in stop.columns:
        mwe_stop = stop["StopDepth_cm"] * ROCK_DENSITY / 100.0
        axes[1].hist(mwe_stop, bins=60, color="#F44336", alpha=0.85)
        axes[1].set_xlabel("Stopping depth (m.w.e.)", fontsize=12)
        axes[1].set_ylabel("Count", fontsize=12)
        axes[1].set_title("Stopping Depth Distribution", fontsize=13)
        axes[1].grid(True, alpha=0.3)

    _save(fig, outdir, "fig09_stopped.png")


# ── Fig 10: Timing comparison ─────────────────────────────────────────────────
def fig_timing(results, outdir):
    data = [(r["label"], r["n_input"], r["timing_s"])
            for r in results if r["timing_s"] > 0]
    if not data:
        print("Skipping fig10 — no timing information")
        return

    labels  = [d[0] for d in data]
    times   = [d[2] for d in data]
    rates   = [d[1] / d[2] for d in data]
    colors  = [CODE_COLORS.get(l, "#78909C") for l in labels]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.subplots_adjust(wspace=0.45)

    for ax, vals, ylabel, title in [
        (axes[0], times, "Wall time (s)",        "Total Run Time"),
        (axes[1], rates, "Throughput (events/s)", "Simulation Speed"),
    ]:
        bars = ax.bar(labels, vals, color=colors, alpha=0.85, edgecolor="white")
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.grid(True, alpha=0.3, axis="y")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2,
                    b.get_height() + max(vals) * 0.01,
                    f"{v:.1f}" if v < 1e4 else f"{v:.0f}",
                    ha="center", fontsize=9)

    _save(fig, outdir, "fig10_timing.png")


# ── Fig 00: Summary overview (multi-code only) ────────────────────────────────
def fig_summary(results, outdir):
    fig = plt.figure(figsize=(18, 10))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.38)
    # (title, stat_key, ylim_or_"log", y_label)
    panels = [
        ("Transmission (%)",           "Transmission_%",  "log",       "Transmission (%)"),
        ("Mean ΔE (GeV)",              "MeanELoss_GeV",   None,        "Mean ΔE (GeV)"),
        ("Mean MCS Angle (°)",         "MeanAngle_deg",   None,        "Mean MCS deflection (°)"),
        ("MCS Lateral Disp. (cm)",     "MeanLatDisp_cm",  None,        "Mean MCS LatDisp (cm)"),
        ("Mean Exit KE (GeV)",         "MeanExitKE_GeV",  None,        "Mean exit KE (GeV)"),
        ("µ⁻ Fraction",                None,              (0.2, 0.8),  r"$\mu^-$ fraction"),
    ]
    for k, (title, key, ylim, ylabel) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        for i, res in enumerate(results):
            if key:
                s = res["stats"]
                if key not in s.columns: continue
                # Only plot non-NaN rows; for transmission also require > 0
                mask = s[key].notna()
                if key == "Transmission_%":
                    mask = mask & (s[key] > 0)
                sv = s[mask]
                if sv.empty: continue
                ax.plot(sv["MWE"], sv[key], label=res["label"],
                        **line_kw(i, res["color"]))
            else:
                # charge ratio
                m = res["muons"]
                if m is None or "PDG" not in m.columns: continue
                rows = [{"MWE": d * ROCK_DENSITY / 100.0,
                         "r": (m[m["DepthCm"]==d]["PDG"]==13).sum() /
                              max(len(m[m["DepthCm"]==d]), 1)}
                        for d in sorted(m["DepthCm"].unique())
                        if d not in EXCLUDE_DEPTH_CM]
                if rows:
                    cr = pd.DataFrame(rows)
                    ax.plot(cr["MWE"], cr["r"], label=res["label"],
                            **line_kw(i, res["color"]))
        if key is None:
            ax.axhline(0.43, color="k", linestyle=":", alpha=0.5, linewidth=1)
        if ylim == "log":
            ax.set_yscale("log")
            ax.set_ylim(1e-3, 200)
            ax.grid(True, alpha=0.25, which="both")
        else:
            if ylim:
                ax.set_ylim(*ylim)
            ax.grid(True, alpha=0.25)
        ax.set_xlabel("Depth (m.w.e.)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)

    fig.suptitle("Benchmark Overview — All Codes", fontsize=15, fontweight="bold")
    _save(fig, outdir, "fig00_summary.png")


# ── Fig 11: Exit position (x, y) at each depth ───────────────────────────────
def fig_exit_position(results, outdir):
    codes = [r for r in results
             if r["muons"] is not None
             and "ExitX_cm" in r["muons"].columns
             and "ExitY_cm" in r["muons"].columns]
    if not codes:
        print("Skipping fig11 — no exit position data (ExitX_cm / ExitY_cm)")
        return

    all_depths = set()
    for r in codes:
        all_depths.update(r["muons"]["DepthCm"].unique())
    depths = sorted(d for d in all_depths if d not in EXCLUDE_DEPTH_CM)

    ncols  = len(codes)
    nrows  = len(depths)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 4.0 * nrows),
                             squeeze=False)

    for row, d in enumerate(depths):
        mwe = d * ROCK_DENSITY / 100.0
        for col, res in enumerate(codes):
            ax  = axes[row][col]
            m   = res["muons"]
            sub = m[m["DepthCm"] == d]
            if len(sub) < 2:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="#aaa")
                ax.set_title(f"{res['label']}  {int(d)} cm ({mwe:.1f} mwe)", fontsize=9)
                continue
            lim = sub[["ExitX_cm","ExitY_cm"]].abs().quantile(0.99).max()
            lim = max(lim * 1.05, 100.0)
            ax.hexbin(sub["ExitX_cm"], sub["ExitY_cm"],
                      gridsize=40, cmap="YlOrRd", mincnt=1,
                      extent=(-lim, lim, -lim, lim))
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_xlabel("Exit x (cm)", fontsize=9)
            ax.set_ylabel("Exit y (cm)", fontsize=9)
            ax.set_title(f"{res['label']}  {int(d)} cm  ({mwe:.1f} mwe)\n"
                         f"N={len(sub):,}", fontsize=9)
            ax.set_aspect("equal")
            ax.grid(True, alpha=0.25)

    fig.suptitle("Muon Exit Position (x, y) at Each Depth", fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, outdir, "fig11_exit_position.png")


# ── Fig 12: Ratio plots (code / Geant4) ──────────────────────────────────────
def fig_comparison_ratios(results, outdir):
    """Fig 12: ratio of each code relative to Geant4 for transmission and exit KE."""
    g4 = _get(results, "Geant4")
    if g4 is None or len(results) < 2:
        return
    g4s = g4["stats"].set_index("MWE")
    others = [r for r in results if r["label"] != "Geant4"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.40)

    for ax, metric, title, ylabel in [
        (axes[0], "Transmission_%",  "Transmission ratio vs Geant4",    "ratio  (code / Geant4)"),
        (axes[1], "MeanExitKE_GeV", "Mean exit KE ratio vs Geant4",    "ratio  (code / Geant4)"),
    ]:
        ax.axhline(1.0, color="k", linestyle="--", linewidth=1.2, alpha=0.5, label="Geant4 (ref.)")
        ax.axhspan(0.9, 1.1, alpha=0.07, color="green")
        ax.text(0.01, 0.91, "±10 %", transform=ax.transAxes,
                fontsize=8, color="darkgreen", alpha=0.8)
        for j, res in enumerate(others):
            s = res["stats"].set_index("MWE")
            pts = []
            for mwe in sorted(set(s.index) & set(g4s.index)):
                g4v = g4s.loc[mwe, metric] if metric in g4s.columns else float("nan")
                v   = s.loc[mwe, metric]   if metric in s.columns  else float("nan")
                n_g4 = g4s.loc[mwe, "N_transmitted"] if "N_transmitted" in g4s.columns else 999
                if np.isnan(g4v) or g4v <= 0 or np.isnan(v) or n_g4 < 5:
                    continue
                pts.append({"MWE": mwe, "ratio": v / g4v})
            if pts:
                df = pd.DataFrame(pts)
                ax.plot(df["MWE"], df["ratio"],
                        label=res["label"], **line_kw(j + 1, res["color"]))
        ax.set_xlabel("Depth (m.w.e.)", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.set_ylim(0.6, 1.5)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    fig.suptitle("Code Comparison Relative to Geant4", fontsize=14, fontweight="bold")
    _save(fig, outdir, "fig12_comparison_ratios.png")


def _is_grid_scan(muons):
    """Detects discrete energy groups (benchmark grid) vs a continuous spectrum.

    Returns a list of cluster-centre energies (≤20 groups) if the dataset
    looks like a grid scan, otherwise returns None.

    Strategy: sort unique InitKE_GeV values and look for gaps where the ratio
    of consecutive values > 1.5 (factor-of-1.5 jump = new energy group).
    This correctly handles:
      • Exact grid codes (MUSIC/BB/UCMuon): energy values are identical within
        each group → ratio = 1.0 inside, > 1.5 between groups.
      • Geant4 benchmark: Ekin varies slightly per event (float storage) but
        all ~9.894 GeV, ~19.894 GeV, etc. → small ratios inside, ~2.0 between.
      • Real cosmic spectrum: energies span 2.89–365 GeV continuously → no
        large gaps, returns None.
    """
    if muons is None or "InitKE_GeV" not in muons.columns:
        return None
    sorted_e = np.sort(muons["InitKE_GeV"].unique())
    if len(sorted_e) == 0:
        return None
    if len(sorted_e) <= 20:
        return list(sorted_e)   # already a small discrete set

    # Ratio-based gap detection
    ratios      = sorted_e[1:] / np.maximum(sorted_e[:-1], 1e-9)
    big_gap_idx = np.where(ratios > 1.5)[0]
    if len(big_gap_idx) == 0 or len(big_gap_idx) >= 20:
        return None   # continuous spectrum or too fragmented

    # Build cluster centres
    starts  = [0]           + list(big_gap_idx + 1)
    ends    = list(big_gap_idx + 1) + [len(sorted_e)]
    centers = [float(np.mean(sorted_e[s:e])) for s, e in zip(starts, ends)]
    return centers if len(centers) <= 20 else None


# ── Fig 13: Code acceptability for subsurface muography ──────────────────────
def fig_acceptability(results, outdir):
    """Fig 13: From what depth can faster codes replace Geant4 for muography?

    Only compares codes that share the same source type (real spectrum).
    Grid-scan codes (fixed energies) are noted but excluded from the ratio,
    because their integrated transmission is not comparable to the real spectrum.
    See fig16 for the physics-valid comparison at fixed energies.
    """
    g4 = _get(results, "Geant4")
    if g4 is None or len(results) < 2:
        return
    g4s    = g4["stats"].set_index("MWE")

    # Separate real-spectrum vs grid-scan codes
    spectrum_others = [r for r in results
                       if r["label"] != "Geant4"
                       and (r["muons"] is None or _is_grid_scan(r["muons"]) is None)]
    grid_others     = [r for r in results
                       if r["label"] != "Geant4"
                       and r["muons"] is not None
                       and _is_grid_scan(r["muons"]) is not None]

    if not spectrum_others:
        # All non-G4 codes are grid-scan; just note it and skip
        print("\n  [fig13] All non-Geant4 codes used grid-scan (fixed-energy) input.")
        print("  Ratio-to-Geant4 is not meaningful; see fig16 for per-energy comparison.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.38)

    metric_defs = [
        ("Transmission_%",  "Transmission ratio (code / Geant4)"),
        ("MeanExitKE_GeV",  "Mean exit KE ratio (code / Geant4)"),
    ]

    accept_depths = {}

    for ax, (metric, ylabel) in zip(axes, metric_defs):
        ax.axhline(1.0, color="k", linestyle="--", linewidth=1.2,
                   alpha=0.6, label="Geant4 (ref.)")
        ax.axhspan(0.95, 1.05, alpha=0.15, color="green")
        ax.axhspan(0.90, 1.10, alpha=0.07, color="#FF9800")
        ax.text(0.02, 0.93, "±5 %",  transform=ax.transAxes, fontsize=8,
                color="green",  alpha=0.9)
        ax.text(0.02, 0.88, "±10 %", transform=ax.transAxes, fontsize=8,
                color="#E65100", alpha=0.9)

        for j, res in enumerate(spectrum_others):
            s   = res["stats"].set_index("MWE")
            pts = []
            for mwe in sorted(set(s.index) & set(g4s.index)):
                g4v  = g4s.loc[mwe, metric] if metric in g4s.columns else float("nan")
                v    = s.loc[mwe, metric]   if metric in s.columns  else float("nan")
                n_g4 = g4s.loc[mwe, "N_transmitted"] if "N_transmitted" in g4s.columns else 999
                if np.isnan(g4v) or g4v <= 0 or np.isnan(v) or n_g4 < 5:
                    continue
                pts.append({"MWE": mwe, "ratio": v / g4v})
            if not pts:
                continue
            df = pd.DataFrame(pts)
            cidx = results.index(res)
            ax.plot(df["MWE"], df["ratio"],
                    label=res["label"], **line_kw(cidx, res["color"]))
            ok = df[np.abs(df["ratio"] - 1.0) <= 0.10]
            if not ok.empty:
                accept_depths.setdefault(res["label"], {})[metric] = ok["MWE"].min()

        ax.set_xlabel("Depth (m.w.e.)", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_ylim(0.4, 2.2)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    axes[0].set_title("Transmission", fontsize=13)
    axes[1].set_title("Mean Exit KE", fontsize=13)

    fig.suptitle(
        "Code Acceptability for Subsurface Muography — real cosmic spectrum\n"
        "From what depth can faster codes replace Geant4?",
        fontsize=11, fontweight="bold")
    if grid_others:
        names = ", ".join(r["label"] for r in grid_others)
        fig.text(0.5, 0.01,
                 f"{names}: grid-scan source → see fig16 for per-energy comparison",
                 ha="center", va="bottom", fontsize=8, color="#555", style="italic")

    print("\n  ══ Muography Acceptability Summary (within 10% of Geant4) ══")
    print(f"  {'Code':<12}  {'Trans. OK from':>16}  {'KE OK from':>12}  {'Both OK from':>14}")
    print("  " + "─" * 62)
    for res in spectrum_others:
        code  = res["label"]
        d     = accept_depths.get(code, {})
        t_ok  = d.get("Transmission_%",  float("inf"))
        ke_ok = d.get("MeanExitKE_GeV", float("inf"))
        both  = max(t_ok, ke_ok)
        fmt   = lambda v: f"{v:.1f} mwe" if np.isfinite(v) else "never"
        print(f"  {code:<12}  {fmt(t_ok):>16}  {fmt(ke_ok):>12}  {fmt(both):>14}")
    if grid_others:
        print(f"\n  Grid-scan codes (see fig16): "
              f"{', '.join(r['label'] for r in grid_others)}")
    print()

    _save(fig, outdir, "fig13_acceptability.png")


# ── Fig 14: Input vs output energy distribution at each depth ────────────────
def fig_input_output_energy(results, outdir):
    """Fig 14: E_out/E_in distribution at each depth — how much energy is preserved?

    Each panel shows one depth; overlay histograms compare all codes.
    Muons that stopped before that depth are excluded (only survivors shown).
    """
    codes_data = [r for r in results
                  if r["muons"] is not None
                  and "InitKE_GeV" in r["muons"].columns
                  and "ExitKE_GeV" in r["muons"].columns]
    if not codes_data:
        return

    depth_cm_list = [d for d in [100, 1000, 2500, 5000, 10000, 20000]
                     if d not in EXCLUDE_DEPTH_CM]
    depth_labels  = [lbl for d, lbl in zip([100, 1000, 2500, 5000, 10000, 20000],
                                            ["1 m", "10 m", "25 m", "50 m", "100 m", "200 m"])
                     if d not in EXCLUDE_DEPTH_CM]

    ncols = min(3, len(depth_cm_list))
    nrows = (len(depth_cm_list) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).flatten()
    fig.subplots_adjust(hspace=0.42, wspace=0.35)

    for idx, (depth_cm, dlabel) in enumerate(zip(depth_cm_list, depth_labels)):
        ax = axes[idx]
        for j, res in enumerate(codes_data):
            m   = res["muons"]
            sub = m[(m["DepthCm"] == depth_cm) & (m["InitKE_GeV"] > 0.1)]
            if len(sub) < 20:
                continue
            frac = (sub["ExitKE_GeV"] / sub["InitKE_GeV"]).clip(0, 1.5)
            ax.hist(frac, bins=50, range=(0, 1.3), density=True, histtype="step",
                    label=f"{res['label']}  μ={frac.mean():.3f}",
                    color=res["color"], linewidth=1.6)
        ax.axvline(1.0, color="gray", linestyle=":", linewidth=1.2, alpha=0.7)
        ax.set_xlabel(r"$E_\mathrm{out}\,/\,E_\mathrm{in}$", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.set_title(f"Depth: {dlabel}", fontsize=11)
        ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.3)

    fig.suptitle(
        r"Energy Survival Fraction $E_\mathrm{out}/E_\mathrm{in}$ at Each Depth"
        "\n(surviving muons only; same source for all codes)",
        fontsize=13, fontweight="bold")
    _save(fig, outdir, "fig14_input_output_energy.png")


# ── Fig 15: Spatial blur for muography ───────────────────────────────────────
def fig_spatial_blur(results, outdir):
    """Fig 15: Position spread vs depth — key for muography image resolution.

    Left:  MCS transverse displacement (deviation from straight-line track).
    Right: Total horizontal displacement from entry position to exit position
           (combines slant drift + MCS; same source entry positions for all codes).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.subplots_adjust(wspace=0.40)

    # Left: MCS lateral displacement from straight-line projection
    plotted = 0
    for i, res in enumerate(results):
        s = res["stats"]
        if "MeanLatDisp_cm" not in s.columns:
            continue
        valid = s[s["MeanLatDisp_cm"].notna()]
        if valid.empty:
            continue
        axes[0].plot(valid["MWE"], valid["MeanLatDisp_cm"],
                     label=res["label"], **line_kw(i, res["color"]))
        plotted += 1
    axes[0].set_xlabel("Depth (m.w.e.)", fontsize=12)
    axes[0].set_ylabel("Mean MCS displacement (cm)", fontsize=12)
    axes[0].set_title("MCS Blur: Deviation from Straight-Line Track\n"
                      "(reconstruction accuracy for muography)", fontsize=11)
    if plotted:
        axes[0].legend(fontsize=9)
    else:
        axes[0].text(0.5, 0.5, "No MCS data available",
                     ha="center", va="center", transform=axes[0].transAxes,
                     fontsize=11, color="#777")
    axes[0].grid(True, alpha=0.3)

    # Right: total displacement from entry point (xs,ys) to exit (x,y)
    codes_entry = [r for r in results
                   if r["muons"] is not None
                   and "EntryX_cm" in r["muons"].columns
                   and "ExitX_cm"  in r["muons"].columns]
    if codes_entry:
        for i, res in enumerate(codes_entry):
            m    = res["muons"]
            rows = []
            for d in sorted(m["DepthCm"].unique()):
                if d in EXCLUDE_DEPTH_CM:
                    continue
                sub  = m[m["DepthCm"] == d]
                if len(sub) < 2:
                    continue
                disp = np.sqrt((sub["ExitX_cm"] - sub["EntryX_cm"])**2
                               + (sub["ExitY_cm"] - sub["EntryY_cm"])**2)
                rows.append({"MWE": d * ROCK_DENSITY / 100.0,
                             "Mean": disp.mean(), "RMS": disp.std()})
            if rows:
                df_r = pd.DataFrame(rows)
                cidx = results.index(res)
                axes[1].plot(df_r["MWE"], df_r["Mean"],
                             label=res["label"], **line_kw(cidx, res["color"]))
        axes[1].set_xlabel("Depth (m.w.e.)", fontsize=12)
        axes[1].set_ylabel("Mean total displacement from entry (cm)", fontsize=12)
        axes[1].set_title("Total Lateral Spread from Entry Position\n"
                          "(slant drift + MCS; same entry positions for all codes)", fontsize=11)
        axes[1].legend(fontsize=9)
        missing = [r["label"] for r in results
                   if r not in codes_entry
                   and r["muons"] is not None]
        if missing:
            axes[1].text(0.98, 0.02,
                         f"{', '.join(missing)}: entry position\nnot in output CSV",
                         transform=axes[1].transAxes, ha="right", va="bottom",
                         fontsize=8, color="#777", style="italic")
    else:
        axes[1].text(0.5, 0.5,
                     "Entry position not available\n(dat-format codes only)",
                     ha="center", va="center", transform=axes[1].transAxes,
                     fontsize=11, color="#777")
    axes[1].grid(True, alpha=0.3)

    fig.suptitle("Spatial Blur vs Depth — Impact on Muography Image Resolution",
                 fontsize=13, fontweight="bold")
    _save(fig, outdir, "fig15_spatial_blur.png")


# ── Fig 16: Transmission per energy (grid-scan comparison) ───────────────────
def _merge_energy_levels(grid_results):
    """Collect cluster-centre energies from all grid-scan results and merge
    energies that are within 30% of each other (same physics energy level)."""
    all_c = []
    for r in grid_results:
        cs = _is_grid_scan(r["muons"])
        if cs:
            all_c.extend(cs)
    if not all_c:
        return []
    all_c = sorted(all_c)
    merged, used = [], set()
    for i, c in enumerate(all_c):
        if i in used:
            continue
        group = [c]
        for j in range(i + 1, len(all_c)):
            if j not in used and all_c[j] / c < 1.30:
                group.append(all_c[j])
                used.add(j)
        used.add(i)
        merged.append(float(np.mean(group)))
    return sorted(merged)


def fig_transmission_per_energy(results, outdir):
    """Fig 16: Survival fraction vs depth at each fixed initial energy.

    All grid-scan codes (including Geant4 benchmark) are treated uniformly.
    Each code is matched to a canonical energy level within ±35%, which
    handles the small difference between e.g. Geant4's 9.894 GeV and
    MUSIC's 10.000 GeV at the '10 GeV' benchmark point.

    n_input per energy:
      • Grid codes (exact): n_input // n_energy_groups
      • Geant4 benchmark:   count at shallowest depth (proxy — all high-E
                            muons survive 1 m, so this ≈ actual input count)
    """
    grid_results = [r for r in results
                    if r["muons"] is not None
                    and _is_grid_scan(r["muons"]) is not None]
    if not grid_results:
        return

    energies = _merge_energy_levels(grid_results)
    if not energies:
        return
    n_e   = len(energies)
    ncols = 4
    nrows = (n_e + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows))
    axes = np.array(axes).flatten()
    fig.subplots_adjust(hspace=0.52, wspace=0.38)

    for ei, E_GeV in enumerate(energies):
        ax = axes[ei]
        ax.set_yscale("log")

        # Bounds: midpoints between adjacent canonical energies (±35% fallback)
        lo_e = energies[ei - 1] * 0.60 if ei > 0 else E_GeV * 0.65
        hi_e = energies[ei + 1] * 0.60 if ei < n_e - 1 else E_GeV * 1.35
        lo_e = max(lo_e, E_GeV * 0.65)
        hi_e = min(hi_e, E_GeV * 1.35)

        for j, res in enumerate(grid_results):
            m   = res["muons"]
            m_e = m[m["InitKE_GeV"].between(lo_e, hi_e)]
            if m_e.empty:
                continue

            n_clusters = len(_is_grid_scan(m))
            # For codes with exact energies, n_input is evenly distributed.
            # For Geant4 benchmark, use shallowest-depth count as proxy.
            if m["InitKE_GeV"].nunique() <= 20:
                n_input_e = res["n_input"] // max(n_clusters, 1)
            else:
                d_min = m_e["DepthCm"].min()
                n_input_e = len(m_e[m_e["DepthCm"] == d_min])

            rows = []
            for d in sorted(m_e["DepthCm"].unique()):
                if d in EXCLUDE_DEPTH_CM:
                    continue
                sub = m_e[m_e["DepthCm"] == d]
                rows.append({"MWE": d * ROCK_DENSITY / 100.0,
                             "Trans": 100.0 * len(sub) / max(n_input_e, 1)})
            if rows:
                cidx = results.index(res)
                df_e = pd.DataFrame(rows)
                ax.plot(df_e["MWE"], df_e["Trans"],
                        label=res["label"], **line_kw(cidx, res["color"]))

        ax.set_xlabel("Depth (m.w.e.)", fontsize=10)
        ax.set_ylabel("Survival fraction (%)", fontsize=10)
        ax.set_title(f"E₀ ≈ {E_GeV:.0f} GeV  (vertical)", fontsize=11)
        ax.set_ylim(5e-4, 200)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which="both")

    for ai in range(n_e, len(axes)):
        axes[ai].set_visible(False)

    fig.suptitle(
        "Transmission vs Depth at Fixed Initial Energy\n"
        "Vertical muons, same surface positions for all codes.\n"
        "All grid-scan datasets matched within ±35% energy window.",
        fontsize=12, fontweight="bold")
    _save(fig, outdir, "fig16_transmission_per_energy.png")


# ── Fig 17: Energy loss vs initial energy (Bethe-Bloch validation) ────────────
def fig_eloss_vs_energy(results, outdir):
    """Fig 17: Mean energy loss vs initial energy at each depth.

    For grid-scan codes each point is exact (100K muons per energy).
    For Geant4 the real spectrum is binned logarithmically.
    Black dashed line = Bethe-Bloch prediction (Groom 2001 parametrisation).
    This directly answers 'which code's dE/dx model agrees with Geant4?'
    """
    codes_ke = [r for r in results
                if r["muons"] is not None
                and "InitKE_GeV" in r["muons"].columns
                and "ExitKE_GeV" in r["muons"].columns]
    if not codes_ke:
        return

    depths_cm    = [d for d in [100, 1000, 2500, 5000, 10000, 20000]
                    if d not in EXCLUDE_DEPTH_CM]
    depth_labels = [lbl for d, lbl in zip([100, 1000, 2500, 5000, 10000, 20000],
                                           ["1 m", "10 m", "25 m", "50 m", "100 m", "200 m"])
                    if d not in EXCLUDE_DEPTH_CM]

    ncols = min(3, len(depths_cm))
    nrows = (len(depths_cm) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    axes = axes.flatten()
    fig.subplots_adjust(hspace=0.42, wspace=0.35)

    for idx, (depth_cm, dlabel) in enumerate(zip(depths_cm, depth_labels)):
        ax    = axes[idx]
        mwe   = depth_cm * ROCK_DENSITY / 100.0

        for j, res in enumerate(codes_ke):
            m   = res["muons"]
            m_d = m[m["DepthCm"] == depth_cm]
            if len(m_d) < 5:
                continue
            cidx = results.index(res)

            grid_e = _is_grid_scan(m)
            if grid_e:
                rows = []
                for E_GeV in grid_e:
                    # ±35% window; avoids fixed-0.5 GeV that fails at 10 TeV
                    sub = m_d[m_d["InitKE_GeV"].between(E_GeV * 0.65, E_GeV * 1.35)]
                    if len(sub) < 5:
                        continue
                    rows.append({"E": sub["InitKE_GeV"].mean(),
                                "dE": (sub["InitKE_GeV"] - sub["ExitKE_GeV"]).mean()})
                if not rows:
                    continue
                df_plot = pd.DataFrame(rows)
                ax.plot(df_plot["E"], df_plot["dE"], "o-",
                        label=res["label"], color=res["color"],
                        linewidth=1.5, markersize=6)
            else:
                # Real spectrum: bin by energy on log scale
                log_bins = np.linspace(np.log10(m_d["InitKE_GeV"].clip(lower=0.1).min()),
                                       np.log10(m_d["InitKE_GeV"].max()), 18)
                m_d = m_d.copy()
                m_d["Ebin"] = pd.cut(np.log10(m_d["InitKE_GeV"].clip(lower=0.1)),
                                     bins=log_bins)
                rows = []
                for _, grp in m_d.groupby("Ebin", observed=True):
                    if len(grp) < 5:
                        continue
                    rows.append({"E":  grp["InitKE_GeV"].mean(),
                                "dE": (grp["InitKE_GeV"] - grp["ExitKE_GeV"]).mean()})
                if not rows:
                    continue
                df_plot = pd.DataFrame(rows).sort_values("E")
                ax.plot(df_plot["E"], df_plot["dE"], ".",
                        label=res["label"], color=res["color"],
                        markersize=5, alpha=0.8)

        # Bethe-Bloch reference clipped to source energy range
        E_bb  = np.logspace(0.5, np.log10(400), 200)
        dE_bb = bethe_bloch_mwe(E_bb) * mwe
        ax.plot(E_bb, dE_bb, "k--", linewidth=1.5, alpha=0.7, label="Bethe-Bloch")

        ax.set_xscale("log")
        ax.set_xlim(3, 400)
        ax.set_xlabel("Initial KE (GeV)", fontsize=10)
        ax.set_ylabel("Mean ΔE = E_in − E_out (GeV)", fontsize=10)
        ax.set_title(f"Depth: {dlabel}", fontsize=11)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which="both")

    fig.suptitle(
        "Mean Energy Loss vs Initial Energy at Each Depth\n"
        "Grid codes (circles, exact) + Geant4 (dots, binned) vs Bethe-Bloch (dashed)\n"
        "Agreement with Bethe-Bloch validates the dE/dx model of each code",
        fontsize=12, fontweight="bold")
    _save(fig, outdir, "fig17_eloss_vs_energy.png")


# ── Quantitative comparison table ─────────────────────────────────────────────
def print_comparison_table(results, outdir):
    """Print and save a side-by-side comparison table: transmission, exit KE, MCS angle."""
    g4 = _get(results, "Geant4")
    if g4 is None or len(results) < 2:
        return
    g4s    = g4["stats"].set_index("MWE")
    codes  = [r["label"] for r in results]
    depths = sorted(g4s.index.tolist())
    all_stats = {r["label"]: r["stats"].set_index("MWE") for r in results}

    def _row(mwe, metric, fmt=".3f"):
        vals, ratios = [], []
        g4v = g4s.loc[mwe, metric] if metric in g4s.columns else float("nan")
        for code in codes:
            s   = all_stats[code]
            val = s.loc[mwe, metric] if (mwe in s.index and metric in s.columns) else float("nan")
            vals.append(val)
            if code != "Geant4" and not np.isnan(g4v) and g4v > 0 and not np.isnan(val):
                ratios.append(f"{val/g4v:+.2f}x")
            else:
                ratios.append("—")
        return vals, ratios

    sep = "─" * 88
    header = f"  {'MWE':>6}  " + "  ".join(f"{c:>12}" for c in codes)
    ratio_header = f"  {'MWE':>6}  " + "  ".join(f"{'Δ/G4 '+c:>12}" for c in codes[1:])

    rows_trans, rows_ke, rows_ang = [], [], []

    for metric_name, metric_key, metric_rows, fmt in [
        ("Transmission (%)",     "Transmission_%",   rows_trans, ".3f"),
        ("Mean Exit KE (GeV)",   "MeanExitKE_GeV",  rows_ke,    ".3f"),
        ("Mean MCS Angle (°)",   "MeanAngle_deg",    rows_ang,   ".3f"),
    ]:
        print(f"\n{'='*88}")
        print(f"  {metric_name}")
        print(sep)
        print(header)
        print(sep)
        for mwe in depths:
            g4n = g4s.loc[mwe, "N_transmitted"] if "N_transmitted" in g4s.columns else 999
            flag = " ⚠ low-N" if g4n < 5 else ""
            vals, ratios = _row(mwe, metric_key, fmt)
            line = f"  {mwe:>6.1f}  " + "  ".join(
                f"{v:>12{fmt}}" if not np.isnan(v) else f"  {'—':>10}" for v in vals)
            line += flag
            print(line)
            row = {"MWE": mwe}
            for code, v in zip(codes, vals):
                row[code] = v
            for code, r in zip(codes[1:], ratios[1:]):
                row[f"ratio_{code}"] = r
            note = ("dE/dx model diff: Geant4 5 GeV muon range ≈ 1000 cm plane"
                    if abs(mwe - 26.5) < 0.1 else "")
            row["Note"] = note
            metric_rows.append(row)
        print(sep)

    # Save CSV
    save_rows = []
    for metric_name, metric_rows in [
        ("Transmission_%",  rows_trans),
        ("MeanExitKE_GeV", rows_ke),
        ("MeanAngle_deg",  rows_ang),
    ]:
        for row in metric_rows:
            r2 = {"Quantity": metric_name}; r2.update(row)
            save_rows.append(r2)
    if save_rows:
        path = os.path.join(outdir, "comparison_metrics.csv")
        pd.DataFrame(save_rows).to_csv(path, index=False)
        print(f"\n  Quantitative comparison saved → {path}")


# ── Helpers ───────────────────────────────────────────────────────────────────
def _save(fig, outdir, fname):
    path = os.path.join(outdir, fname)
    fig.savefig(path)          # dpi=300, bbox="tight" set globally in rcParams
    plt.close(fig)
    print(f"  Saved {path}")

def _get(results, label):
    return next((r for r in results if r["label"] == label), None)

def _first_with(results, col):
    return next((r for r in results
                 if r["muons"] is not None and col in r["muons"].columns), None)

# ═════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="MuonRock benchmark — compare Geant4 / PHITS / MUSIC / PROPOSAL")
    parser.add_argument("run_dir",    nargs="?",
                        help="Geant4 run directory (positional, backward-compatible)")
    parser.add_argument("--geant4",   metavar="DIR")
    parser.add_argument("--phits",    metavar="DIR")
    parser.add_argument("--music",    metavar="DIR")
    parser.add_argument("--proposal", metavar="DIR")
    parser.add_argument("--bb",       metavar="DIR")
    parser.add_argument("--ucmuon",   metavar="DIR")
    parser.add_argument("--n",        type=int, metavar="N",
                        help="Override n_input for all codes")
    parser.add_argument("--outdir",   metavar="DIR", default=None,
                        help="Directory for output files  [default: <outputs>/<run_name>/]")
    args = parser.parse_args()

    g4_dir = args.geant4 or args.run_dir
    if not g4_dir:
        # Auto-detect: look for run_* in current dir or outputs/
        for candidate in [".", "outputs", "../outputs"]:
            runs = sorted(glob.glob(os.path.join(candidate, "run_*")))
            if runs:
                g4_dir = candidate
                break
        if not g4_dir:
            g4_dir = "."

    print("=" * 66)
    print("  MuonRock v6 — Benchmark Analysis  (UCLouvain / CP3)")
    print("=" * 66)

    results = []

    g4 = load_geant4(g4_dir, n_override=args.n)
    if g4 is None:
        sys.exit("[ERROR] Geant4 output required — aborting.")
    results.append(g4)

    # Default outdir: a subfolder named after the run, inside the outputs directory
    if args.outdir is None:
        base = g4.get("base", "")
        args.outdir = base + "_figures" if base else "."

    os.makedirs(args.outdir, exist_ok=True)

    for loader, flag in [(load_phits,    args.phits),
                          (load_music,    args.music),
                          (load_proposal, args.proposal),
                          (load_bb,       args.bb),
                          (load_ucmuon,   args.ucmuon)]:
        if flag:
            r = loader(flag, n_override=args.n)
            if r:
                results.append(r)

    # ── Print per-depth table ──────────────────────────────────────────────
    print(f"\n{'─'*66}")
    print(f"  Loaded {len(results)} code(s): "
          f"{', '.join(r['label'] for r in results)}")
    print(f"{'─'*66}")
    for res in results:
        s = res["stats"]
        print(f"\n  [{res['label']}]  n_input={res['n_input']:,}"
              f"  timing={res['timing_s']:.1f}s")
        cols = [c for c in ["DepthCm","MWE","N_transmitted","Transmission_%",
                             "MeanELoss_GeV","MeanAngle_deg","MeanLatDisp_cm"]
                if c in s.columns]
        print(s[cols].to_string(index=False))

    # ── Write combined CSV ─────────────────────────────────────────────────
    frames = []
    for res in results:
        s = res["stats"].copy()
        s.insert(0, "Code", res["label"])
        frames.append(s)
    csv_path = os.path.join(args.outdir, "benchmark_summary.csv")
    pd.concat(frames, ignore_index=True).to_csv(csv_path, index=False)
    print(f"\nSaved {csv_path}")

    # ── Quantitative comparison table ─────────────────────────────────────
    if len(results) > 1:
        print_comparison_table(results, args.outdir)

    # ── Generate figures ───────────────────────────────────────────────────
    print(f"\nGenerating figures → {args.outdir}/")
    fig_transmission(results, args.outdir)
    fig_energy_loss(results, args.outdir)
    fig_angular(results, args.outdir)
    fig_lateral(results, args.outdir)
    fig_exit_spectrum(results, args.outdir)
    fig_dedx(results, args.outdir)
    fig_charge_ratio(results, args.outdir)
    fig_secondaries(results, args.outdir)
    fig_stopped(results, args.outdir)
    fig_timing(results, args.outdir)
    fig_exit_position(results, args.outdir)
    if len(results) > 1:
        fig_summary(results, args.outdir)
        fig_comparison_ratios(results, args.outdir)
        fig_acceptability(results, args.outdir)
        fig_input_output_energy(results, args.outdir)
        fig_spatial_blur(results, args.outdir)
        fig_transmission_per_energy(results, args.outdir)
        fig_eloss_vs_energy(results, args.outdir)

    print("\n" + "=" * 66)
    print("  Done.  All outputs in: " + args.outdir)
    print("=" * 66)


if __name__ == "__main__":
    main()
