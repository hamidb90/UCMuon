#!/usr/bin/env python3
"""
make_fig_realistic_survival.py — survival-vs-depth for the six benchmark codes
on the REALISTIC cosmic-ray source (not the monoenergetic 6-beam grid).

Reads the distilled six-code summary produced by
benchmark/analysis/analyze_sixcode.py:
    benchmark/run_realistic_20260614/benchmark_sixcode_summary.csv
(columns: code, depth_m, mwe, survival_pct, mean_exitKE_GeV)

Output: manuscript/figs/fig_survival_realistic.pdf

Companion of make_fig06.py (monoenergetic). Same styling so the two read as a
set. BB_fortran is excluded (duplicate of the CSDA Bethe-Bloch engine).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SUMMARY   = REPO_ROOT / "benchmark" / "run_realistic_20260614" / "benchmark_sixcode_summary.csv"
OUT_PDF   = Path(__file__).resolve().parent.parent / "figs" / "fig_survival_realistic.pdf"

ROCK_DENSITY = 2.65
N_SRC        = 100_000

# Same style/colours as make_fig06.py
CODE_STYLE = {
    "Geant4":   {"label": "Geant4 11.2 (FTFP_BERT)",      "color": "#1565C0", "marker": "o", "ls": "-",  "lw": 2.2, "ms": 6, "zorder": 5},
    "PHITS":    {"label": "PHITS 3.36",                   "color": "#C62828", "marker": "s", "ls": "-.", "lw": 2.2, "ms": 6, "zorder": 5},
    "MUSIC":    {"label": "MUSIC (Engine 2)",             "color": "#FF8F00", "marker": "^", "ls": "-",  "lw": 1.6, "ms": 5, "zorder": 4},
    "PROPOSAL": {"label": "PROPOSAL (Engine 4)",          "color": "#2E7D32", "marker": "v", "ls": "--", "lw": 1.6, "ms": 5, "zorder": 4},
    "BB":       {"label": "Bethe–Bloch CSDA (Engine 3)", "color": "#7B1FA2", "marker": "D", "ls": ":",  "lw": 1.6, "ms": 5, "zorder": 4},
    "UCMuon":   {"label": "UCMuon-MC (Engine 1)",         "color": "#00838F", "marker": "P", "ls": "--", "lw": 1.6, "ms": 5, "zorder": 4},
}


def main() -> None:
    df = pd.read_csv(SUMMARY)
    df.columns = df.columns.str.strip()

    codes = {}
    for name in CODE_STYLE:                      # ORDER + excludes BB_fortran
        sub = df[df["code"] == name].sort_values("depth_m")
        if sub.empty:
            print(f"  [skip] {name}: not in summary")
            continue
        codes[name] = {
            "mwe":  sub["mwe"].to_numpy(),
            "surv": sub["survival_pct"].to_numpy() / 100.0,
        }
        print(f"  [ok] {name}: {list(sub['depth_m'])} m")

    mwe_ticks = sorted(df["mwe"].unique())
    depths_m  = [round(m / ROCK_DENSITY) for m in mwe_ticks]

    plt.rcParams.update({
        "font.family": "serif", "font.size": 11, "axes.labelsize": 12,
        "legend.fontsize": 8.5, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7, 6.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        layout="constrained")

    g4 = codes.get("Geant4")
    for name, style in CODE_STYLE.items():
        if name not in codes:
            continue
        c = codes[name]
        surv = c["surv"]; mwe = c["mwe"]
        err = np.sqrt(surv * (1 - surv) / N_SRC)
        ax_top.errorbar(mwe, surv * 100, yerr=err * 100,
                        fmt=f"{style['marker']}{style['ls']}", color=style["color"],
                        lw=style["lw"], ms=style["ms"], capsize=3,
                        label=style["label"], zorder=style["zorder"])
        if g4 is not None and name != "Geant4":
            ref = np.interp(mwe, g4["mwe"], g4["surv"], left=np.nan, right=np.nan)
            m = (ref > 1e-4) & ~np.isnan(ref)
            if m.any():
                # Percentage points (matches the prose and the six-code
                # table); relative % blows past any sane axis at 530 mwe.
                ax_bot.plot(mwe[m], 100 * (surv[m] - ref[m]),
                            f"{style['marker']}{style['ls']}", color=style["color"],
                            lw=style["lw"], ms=style["ms"])

    ax_top.set_ylabel("Survival fraction [%]")
    ax_top.legend(loc="upper right", framealpha=0.9)
    ax_top.set_ylim(bottom=0)
    ax_top.set_title("Realistic cosmic-ray source ($10^5$ muons)", fontsize=11, loc="left")

    ax_bot.axhline(0, color="#1565C0", lw=1.2)
    ax_bot.axhspan(-1, 1, color="green", alpha=0.10)
    ax_bot.axhspan(-3, 3, color="orange", alpha=0.07)
    ax_bot.axhline(+1, color="green", lw=0.8, ls=":")
    ax_bot.axhline(-1, color="green", lw=0.8, ls=":")
    ax_bot.set_ylabel("Δ vs Geant4 [pp]")
    ax_bot.set_ylim(-9.5, 3)
    ax_bot.set_xlabel("Overburden [m.w.e.]")
    ax_bot.set_xticks(mwe_ticks)
    ax_bot.set_xticklabels([f"{m:.0f}" for m in mwe_ticks])

    secax = ax_top.secondary_xaxis("top", functions=(
        lambda x: x / ROCK_DENSITY, lambda x: x * ROCK_DENSITY))
    secax.set_xlabel("Vertical depth [m]")
    secax.set_xticks(depths_m)
    secax.set_xticklabels(depths_m)

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    plt.close(fig)
    print(f"  -> {OUT_PDF}")


if __name__ == "__main__":
    main()
