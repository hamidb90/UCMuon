#!/usr/bin/env python3
"""
make_fig06.py  —  Generate fig06_survival_curve.pdf for the UCMuon CPC paper.

Loads benchmark data from:
  - benchmark/geant4_muon_rock_v5/<ENGINE>/  (BB, MUSIC, PROPOSAL, UCMuon .dat files)
  - Geant4 CSV and PHITS summary (inside same scratch tree, override via CLI)

Output: manuscript/figs/fig06_survival_curve.pdf

Usage (from UCMuon repo root):
  python manuscript/make_fig06.py
  python manuscript/make_fig06.py \
      --geant4 /path/to/geant4/outputs/ \
      --phits  /path/to/phits/
"""
from __future__ import annotations
import argparse, glob, os, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Paths (relative to repo root) ─────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parent.parent.parent
BENCH_INDIR = REPO_ROOT / "benchmark" / "geant4_muon_rock_v5"
# G4/PHITS mono references come from the distilled v2 summary by default:
# the raw mono CSVs were pruned, and outputs/*_muons.csv + phits_summary.csv
# now hold the June REALISTIC-source runs — picking them by mtime silently
# plots realistic references against monoenergetic engines (wrong by design).
REF_SUMMARY = str(REPO_ROOT / "benchmark" / "geant4_muon_rock_v5" /
                  "figures_benchmark" / "benchmark_summary.csv")
G4_DEFAULT  = None
PHITS_DEFAULT = None
OUT_PDF     = Path(__file__).resolve().parent.parent / "figs/fig06_survival_curve.pdf"

# ── Physics constants ──────────────────────────────────────────────────────────
ROCK_DENSITY = 2.65   # g/cm³  Standard Rock
M_MU_GEV     = 0.10566
DEPTHS_M     = [1, 25, 50, 100, 200]
DEPTHS_MWE   = [d * ROCK_DENSITY for d in DEPTHS_M]
N_INPUT      = 600_000

# ── Style ──────────────────────────────────────────────────────────────────────
CODE_STYLE: dict[str, dict] = {
    "Geant4":   {"label": "Geant4 11.2 (FTFP_BERT)",
                 "color": "#1565C0", "marker": "o", "ls": "-",  "lw": 2.2, "ms": 6, "zorder": 5},
    "PHITS":    {"label": "PHITS 3.36",
                 "color": "#C62828", "marker": "s", "ls": "-.", "lw": 2.2, "ms": 6, "zorder": 5},
    "MUSIC":    {"label": "MUSIC (Engine 2)",
                 "color": "#FF8F00", "marker": "^", "ls": "-",  "lw": 1.6, "ms": 5, "zorder": 4},
    "PROPOSAL": {"label": "PROPOSAL (Engine 4)",
                 "color": "#2E7D32", "marker": "v", "ls": "--", "lw": 1.6, "ms": 5, "zorder": 4},
    "BB":       {"label": "Bethe–Bloch CSDA (Engine 3)",
                 "color": "#7B1FA2", "marker": "D", "ls": ":",  "lw": 1.6, "ms": 5, "zorder": 4},
    "UCMuon":   {"label": "UCMuon-MC (Engine 1)",
                 "color": "#00838F", "marker": "P", "ls": "--", "lw": 1.6, "ms": 5, "zorder": 4},
}

# ── Groom (2001) Standard Rock CSDA range table ────────────────────────────────
_GT = np.array([
    0.01,0.014,0.02,0.03,0.04,0.08,0.10,0.14,0.20,0.30,
    0.40,0.80,1.00,1.40,2.00,3.00,4.00,8.00,10.0,14.0,
    20.0,30.0,40.0,80.0,100.,140.,200.,300.,400.,800.,1000.,
])
_GR = np.array([
    0.8516,1.542,2.866,5.698,9.145,26.76,36.96,58.79,93.32,152.4,
    211.5,441.8,553.4,771.2,1088.,1599.,2095.,3998.,4920.,6724.,
    9360.,13620.,17760.,33430.,40840.,55460.,76650.,107900.,136100.,225300.,272200.,
])


def csda_survival(depths_m: np.ndarray, spectrum_gev: list[tuple[float, float]]) -> np.ndarray:
    """
    Semi-analytical CSDA survival fraction for a multi-energy source.
    spectrum_gev = [(E_GeV, fraction), ...], fractions must sum to 1.
    A muon of energy E survives depth d iff its CSDA range > rho*d.
    """
    survs = np.zeros(len(depths_m))
    for E0, frac in spectrum_gev:
        R0 = float(np.exp(np.interp(np.log(E0), np.log(_GT), np.log(_GR))))
        for i, dm in enumerate(depths_m):
            opacity = dm * ROCK_DENSITY * 100  # g/cm²  (cm → same units as _GR)
            survs[i] += frac * float(R0 > opacity)
    return survs


# Source beams: 5,10,20,50,100,300 GeV, 10^5 each → equal fractions
_SOURCE = [(e, 1/6) for e in [5, 10, 20, 50, 100, 300]]


# ── Loaders (subset of benchmark_analysis.py) ─────────────────────────────────
_DAT_COLS = ["EventID","xs","ys","zs","Es","ts","ps","charge","alive",
             "x","y","z","E","cx","cy","cz","theta_ug","phi_ug"]


def _find_dat(indir: Path, eng: str, dm: int) -> Path | None:
    for p in [indir / f"{eng}_bench_{dm}m.dat",
              indir / eng / f"{eng}_bench_{dm}m.dat"]:
        if p.exists():
            return p
    return None


def load_ucmuon_engines(indir: Path) -> dict[str, dict]:
    results: dict[str, dict] = {}
    for eng in ["BB", "MUSIC", "PROPOSAL", "UCMuon"]:
        depths, survs, n_tots = [], [], []
        for dm in DEPTHS_M:
            path = _find_dat(indir, eng, dm)
            if path is None:
                continue
            df = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                             names=_DAT_COLS, engine="python")
            n_total = len(df)
            n_alive = int((df["alive"] == 1).sum())
            depths.append(dm)
            survs.append(n_alive / n_total if n_total > 0 else 0.0)
            n_tots.append(n_total)
        if not depths:
            print(f"  [SKIP] {eng}: no files found in {indir}")
            continue
        results[eng] = {"depths": np.array(depths),
                        "survival": np.array(survs),
                        "n_total": np.array(n_tots)}
        print(f"  [OK]  {eng}: {depths} m")
    return results


def load_geant4(g4_dir: str) -> dict | None:
    if os.path.isdir(g4_dir):
        cands = sorted(glob.glob(os.path.join(g4_dir, "*_muons.csv")),
                       key=os.path.getmtime)
        if not cands:
            print(f"  [SKIP] No *_muons.csv in {g4_dir}")
            return None
        path = cands[-1]
    else:
        path = g4_dir
        if not os.path.exists(path):
            print(f"  [SKIP] {path} not found")
            return None
    print(f"  [G4]  {os.path.basename(path)}")
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    depths, survs, n_tots = [], [], []
    for dm in DEPTHS_M:
        sub = df[df["DepthCm"] == dm * 100]
        if sub.empty:
            continue
        depths.append(dm)
        survs.append(len(sub) / N_INPUT)
        n_tots.append(N_INPUT)
        print(f"  [G4]  {dm}m: {len(sub):,}/{N_INPUT:,} ({100*len(sub)/N_INPUT:.2f}%)")
    return {"depths": np.array(depths), "survival": np.array(survs),
            "n_total": np.array(n_tots)}


def load_phits(phits_dir: str) -> dict | None:
    summ = os.path.join(phits_dir, "phits_summary.csv")
    if not os.path.exists(summ):
        print(f"  [SKIP] {summ} not found")
        return None
    df = pd.read_csv(summ)
    df.columns = df.columns.str.strip()
    depths, survs, n_tots = [], [], []
    tr0 = df.iloc[0]["Transmission_%"]; n0 = df.iloc[0]["N_transmitted"]
    n_input = int(round(n0 / (tr0 / 100.0))) if tr0 > 0 else N_INPUT
    for _, row in df.iterrows():
        dm = int(row["DepthCm"]) // 100
        if dm not in DEPTHS_M:
            continue
        tr = row["Transmission_%"] / 100.0
        depths.append(dm)
        survs.append(tr)
        n_tots.append(n_input)
        print(f"  [PHITS] {dm}m: {100*tr:.2f}%")
    return {"depths": np.array(depths), "survival": np.array(survs),
            "n_total": np.array(n_tots)}


def load_from_summary(csv_path: str) -> dict[str, dict]:
    """Load per-depth survival for every code from a benchmark_summary CSV.

    Used when the raw per-event .dat files are not available (e.g. the 8 GB
    scratch tree was pruned). The distilled survival/N values are sufficient
    to reproduce the survival-vs-depth figure exactly.
    """
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()
    keep_cm = [d * 100 for d in DEPTHS_M]
    out: dict[str, dict] = {}
    for code in df["Code"].unique():
        sub = df[(df["Code"] == code) & (df["DepthCm"].isin(keep_cm))]
        sub = sub.sort_values("DepthCm")
        if sub.empty:
            continue
        surv = (sub["Transmission_%"].to_numpy() / 100.0)
        n    = sub["N_transmitted"].to_numpy()
        ntot = np.where(surv > 0, n / surv, N_INPUT)
        out[str(code)] = {"depths": (sub["DepthCm"].to_numpy() / 100.0),
                          "survival": surv, "n_total": ntot}
        print(f"  [SUMMARY] {code}: {list((sub['DepthCm']//100).astype(int))} m")
    return out


# ── Figure ─────────────────────────────────────────────────────────────────────

def make_fig06(all_codes: dict[str, dict], out_pdf: Path) -> None:
    plt.rcParams.update({
        "font.family":     "serif",
        "font.size":       11,
        "axes.labelsize":  12,
        "axes.titlesize":  12,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.grid":       True,
        "grid.alpha":      0.25,
        "grid.linestyle":  "--",
        "figure.dpi":      150,
        "savefig.dpi":     300,
        "savefig.bbox":    "tight",
    })

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7, 6.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        layout="constrained")

    mwe_arr = np.array(DEPTHS_MWE)
    # Fine curve for CSDA
    mwe_fine = np.linspace(mwe_arr[0] * 0.5, mwe_arr[-1] * 1.05, 200)
    dm_fine  = mwe_fine / ROCK_DENSITY
    csda_fine = csda_survival(dm_fine, _SOURCE)
    ax_top.plot(mwe_fine, csda_fine * 100, "k--", lw=1.4,
                label="CSDA (6-beam discrete spectrum)", zorder=2)

    geant4_ref = all_codes.get("Geant4")

    for name, style in CODE_STYLE.items():
        if name not in all_codes:
            continue
        code = all_codes[name]
        mwe  = code["depths"] * ROCK_DENSITY
        surv = code["survival"]
        n    = code["n_total"]
        err  = np.where(n > 0, np.sqrt(surv * (1 - surv) / np.maximum(n, 1)), 0.0)

        ax_top.errorbar(
            mwe, surv * 100, yerr=err * 100,
            fmt=f"{style['marker']}{style['ls']}",
            color=style["color"], lw=style["lw"], ms=style["ms"],
            capsize=3, label=style["label"], zorder=style["zorder"])

        # Bottom panel: % deviation from Geant4
        if geant4_ref is not None and name != "Geant4":
            ref_interp = np.interp(code["depths"],
                                   geant4_ref["depths"], geant4_ref["survival"],
                                   left=np.nan, right=np.nan)
            mask = (ref_interp > 1e-4) & ~np.isnan(ref_interp)
            if mask.any():
                dev = 100 * (surv[mask] - ref_interp[mask]) / ref_interp[mask]
                ax_bot.plot(mwe[mask], dev,
                            f"{style['marker']}{style['ls']}",
                            color=style["color"], lw=style["lw"], ms=style["ms"])

    ax_top.set_ylabel("Survival fraction [%]")
    ax_top.legend(loc="upper right", framealpha=0.9)
    ax_top.set_ylim(bottom=0)

    # Deviation panel decorations
    ax_bot.axhline(0,  color="#1565C0", lw=1.2, ls="-")
    ax_bot.axhspan(-1, 1, color="green",  alpha=0.10)
    ax_bot.axhspan(-3, 3, color="orange", alpha=0.07)
    ax_bot.axhline(+1, color="green",  lw=0.8, ls=":")
    ax_bot.axhline(-1, color="green",  lw=0.8, ls=":")
    ax_bot.set_ylabel("Δ vs Geant4 [%]")
    ax_bot.set_ylim(-5, 5)
    ax_bot.set_xlabel("Overburden [m.w.e.]")
    ax_bot.set_xticks(DEPTHS_MWE)
    ax_bot.set_xticklabels([f"{d:.0f}" for d in DEPTHS_MWE])

    # Secondary depth-in-metres axis on top panel
    secax = ax_top.secondary_xaxis("top", functions=(
        lambda x: x / ROCK_DENSITY, lambda x: x * ROCK_DENSITY))
    secax.set_xlabel("Vertical depth [m]")
    secax.set_xticks(DEPTHS_M)
    secax.set_xticklabels(DEPTHS_M)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"  → {out_pdf}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description="Generate fig06_survival_curve.pdf")
    ap.add_argument("--indir",  default=str(BENCH_INDIR))
    ap.add_argument("--geant4", default=G4_DEFAULT,
                    help="Raw Geant4 *_muons.csv (file or dir); default: "
                         "distilled reference summary")
    ap.add_argument("--phits",  default=PHITS_DEFAULT,
                    help="PHITS dir with phits_summary.csv; default: "
                         "distilled reference summary")
    ap.add_argument("--ref-summary", default=REF_SUMMARY,
                    help="benchmark_summary CSV supplying the Geant4/PHITS "
                         "reference rows (mono v2 run)")
    ap.add_argument("--out",    default=str(OUT_PDF))
    ap.add_argument("--summary", default=None,
                    help="Load all codes from a benchmark_summary CSV "
                         "(use when raw .dat are unavailable)")
    args = ap.parse_args()

    indir   = Path(args.indir)
    out_pdf = Path(args.out)

    print("Loading benchmark data …")
    all_codes: dict[str, dict] = {}

    if args.summary:
        all_codes = load_from_summary(args.summary)
        if not all_codes:
            print("No data loaded from summary — check the CSV.")
            sys.exit(1)
        print(f"\nGenerating {out_pdf} …")
        make_fig06(all_codes, out_pdf)
        print("Done.")
        return

    if args.geant4:
        g4 = load_geant4(args.geant4)
        if g4:
            all_codes["Geant4"] = g4

    if args.phits:
        ph = load_phits(args.phits)
        if ph:
            all_codes["PHITS"] = ph

    # Default: Geant4/PHITS reference rows from the distilled mono summary
    if args.ref_summary and ("Geant4" not in all_codes or "PHITS" not in all_codes):
        refs = load_from_summary(args.ref_summary)
        for code in ("Geant4", "PHITS"):
            if code not in all_codes and code in refs:
                all_codes[code] = refs[code]

    engines = load_ucmuon_engines(indir)
    all_codes.update(engines)

    if not all_codes:
        print("No data loaded — check paths.")
        sys.exit(1)

    print(f"\nGenerating {out_pdf} …")
    make_fig06(all_codes, out_pdf)
    print("Done.")


if __name__ == "__main__":
    main()
