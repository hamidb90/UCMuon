#!/usr/bin/env python3
"""
make_fig_realistic_timing.py — CPU time (single-core-equivalent) of the six
benchmark codes on the REALISTIC cosmic-ray source (10^5 muons), as a
horizontal log-scale bar chart. Companion of make_fig10_timing.py
(monoenergetic 6x10^5 source); same styling and CPU-time convention.

Reads benchmark/run_realistic_20260614/realistic_timing.csv (code, seconds).
For the single-thread codes seconds = wall-clock; for PHITS it is the total CPU
summed over its 10 OpenMP threads (Case A, full physics).

Output: manuscript/figs/fig_timing_realistic.pdf
"""
from __future__ import annotations
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
TIMING    = REPO_ROOT / "benchmark" / "run_realistic_20260614" / "realistic_timing.csv"
OUT_PDF   = Path(__file__).resolve().parent.parent / "figs" / "fig_timing_realistic.pdf"

CODE_STYLE = {
    "Geant4":   {"label": "Geant4 11.2 (FTFP_BERT)",      "color": "#1565C0"},
    "PHITS":    {"label": "PHITS 3.36",                   "color": "#C62828"},
    "MUSIC":    {"label": "MUSIC (Engine 2)",             "color": "#FF8F00"},
    "PROPOSAL": {"label": "PROPOSAL (Engine 4)",          "color": "#2E7D32"},
    "BB":       {"label": "Bethe–Bloch CSDA (Engine 3)", "color": "#7B1FA2"},
    "UCMuon":   {"label": "UCMuon-MC (Engine 1)",         "color": "#00838F"},
}


def fmt_time(s: float) -> str:
    if s < 60:    return f"{s:.0f} s"
    if s < 3600:  return f"{s/60:.0f} min"
    return f"{s/3600:.1f} h"


def fmt_speedup(s: float, ref: float) -> str:
    if s < ref:  return f"{ref/s:.0f}× faster"
    if s > ref:  return f"{s/ref:.1f}× slower"
    return "—"


def main() -> None:
    timings = {}
    with open(TIMING) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            timings[row["code"]] = float(row["seconds"])

    geant4_t = timings["Geant4"]
    order  = sorted(timings, key=lambda c: timings[c], reverse=True)
    times  = [timings[c]             for c in order]
    colors = [CODE_STYLE[c]["color"] for c in order]
    ypos   = list(range(len(order)))

    plt.rcParams.update({
        "font.family": "serif", "font.size": 11, "axes.labelsize": 12,
        "axes.titlesize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    fig.subplots_adjust(left=0.20, right=0.66, top=0.90, bottom=0.16)

    ax.axhspan(-0.5, 1.5, alpha=0.07, color="#1565C0", zorder=0, lw=0)
    ax.axhspan( 1.5, 5.5, alpha=0.07, color="#2E7D32", zorder=0, lw=0)
    ax.barh(ypos, times, color=colors, height=0.62,
            edgecolor="white", linewidth=0.8, zorder=2)
    ax.axvline(geant4_t, color=CODE_STYLE["Geant4"]["color"],
               linestyle="--", linewidth=1.3, alpha=0.5, zorder=3)
    ax.text(geant4_t * 0.86, 2.6, "Geant4\nref.", ha="right", va="center",
            fontsize=8, color=CODE_STYLE["Geant4"]["color"], alpha=0.8, style="italic")

    ax.set_xscale("log")
    ax.set_xlim(5, 1_000_000)
    ax.set_ylim(-0.5, 5.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([CODE_STYLE[c]["label"] for c in order], fontsize=9.5)
    ax.set_xlabel("CPU time, single-core-equivalent  (log scale)", labelpad=6)

    ax.set_xticks([10, 60, 600, 3600, 36000, 360000])
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(
        ["10 s", "1 min", "10 min", "1 h", "10 h", "100 h"]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.tick_params(axis="x", labelsize=9, labelrotation=25, pad=3)
    ax.grid(True, axis="x", which="major", alpha=0.25, linestyle="--")
    ax.set_title("Simulation speed (realistic source, $10^{5}$ muons)",
                 fontsize=12, pad=8, loc="left")

    ax.text(1.04, 1.02, "Time", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, fontweight="bold", color="0.35", clip_on=False)
    ax.text(1.28, 1.02, "vs Geant4", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, fontweight="bold", color="0.35", clip_on=False)
    from matplotlib.transforms import blended_transform_factory
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for y, c in zip(ypos, order):
        ax.text(1.04, y, fmt_time(timings[c]), transform=trans, ha="left",
                va="center", fontsize=9.5, color=CODE_STYLE[c]["color"],
                fontweight="bold", clip_on=False)
        ax.text(1.28, y, fmt_speedup(timings[c], geant4_t), transform=trans,
                ha="left", va="center", fontsize=9, color="0.30", clip_on=False)

    ax.text(0.0, -0.22, "PHITS value is total CPU time over 10 OpenMP threads; "
            "all other codes single-thread (wall = CPU).",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="0.45", style="italic")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    plt.close(fig)
    print(f"[ok] wrote {OUT_PDF}")
    for c in order:
        print(f"     {c:9s} {timings[c]:>10.0f} s   {fmt_speedup(timings[c], geant4_t)}")


if __name__ == "__main__":
    main()
