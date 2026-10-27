#!/usr/bin/env python3
"""
make_fig10_timing.py — CPU time (single-core-equivalent) of the six benchmark
codes on the common 6x10^5-muon Standard-Rock source, rendered as a
horizontal log-scale bar chart (fig10_timing.pdf).

Every code runs on one core of the same machine, so wall time equals CPU
time.

All times are read from the shipped record benchmark/results/
engine_timing_20261027.csv: the TOTAL row of each UCMuon engine (sum over the
six depths, one process each) and the two reference rows, Geant4 and PHITS
(one core each, timed on 1/32 of the population and scaled; see the record's
header). Up to 1.4.1 the script read per-engine timing files from the
unshipped benchmark/geant4_muon_rock_v5/, and up to 2026-10-24 the PHITS
value was the CPU time summed over the 10 threads of an OpenMP run of
PHITS 3.36 (371 146 s).

Usage:  python3 manuscript/scripts/make_fig10_timing.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF    = SCRIPT_DIR.parent / "figs" / "fig10_timing.pdf"
TIMING_CSV = (SCRIPT_DIR.parent.parent / "benchmark" / "results"
              / "engine_timing_20261027.csv")
# CSV engine name -> key used below
CSV_KEYS = {"UCMuon-MC": "UCMuon", "MUSIC": "MUSIC", "BetheBloch": "BB",
            "PROPOSAL": "PROPOSAL", "Geant4": "Geant4", "PHITS": "PHITS"}

# Colours match CODE_STYLE in make_fig06.py so the two figures read as a set.
CODE_STYLE = {
    "Geant4":   {"label": "Geant4 11.4 (FTFP_BERT)",      "color": "#1565C0"},
    "PHITS":    {"label": "PHITS 3.37",                   "color": "#C62828"},
    "MUSIC":    {"label": "MUSIC (Engine 2)",             "color": "#FF8F00"},
    "PROPOSAL": {"label": "PROPOSAL (Engine 4)",          "color": "#2E7D32"},
    "BB":       {"label": "Bethe–Bloch CSDA (Engine 3)", "color": "#7B1FA2"},
    "UCMuon":   {"label": "UCMuon-MC (Engine 1)",         "color": "#00838F"},
}


def load_timings() -> dict[str, float]:
    """Seconds for the full 6x10^5-muon run of each code (TOTAL rows)."""
    t = {}
    with open(TIMING_CSV) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            key = CSV_KEYS.get(row["engine"])
            if key and row["depth_m"] == "TOTAL":     # UCMuon-MC-auto is not a key
                t[key] = float(row["seconds_6e5"])
    missing = set(CSV_KEYS.values()) - set(t)
    if missing:
        raise ValueError(f"{TIMING_CSV.name}: no TOTAL row for {sorted(missing)}")
    return t


def fmt_time(s: float) -> str:
    if s < 60:        return f"{s:.0f} s"
    if s < 600:       return f"{s / 60:.1f} min"
    if s < 3_600:     return f"{s / 60:.0f} min"
    return f"{s / 3_600:.1f} h"


def fmt_speedup(s: float, ref: float) -> str:
    if s < ref:   return f"{ref / s:.0f}× faster"
    if s > ref:   return f"{s / ref:.1f}× slower"
    return "—"


def make_fig10(timings: dict[str, float], out_pdf: Path) -> None:
    plt.rcParams.update({
        "font.family":     "serif",
        "font.size":       11,
        "axes.labelsize":  12,
        "axes.titlesize":  12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.dpi":      150,
        "savefig.dpi":     300,
        "savefig.bbox":    "tight",
    })

    geant4_t = timings["Geant4"]

    # Slowest -> fastest so the fastest engine sits at the top.
    order  = sorted(timings, key=lambda c: timings[c], reverse=True)
    times  = [timings[c]                  for c in order]
    colors = [CODE_STYLE[c]["color"]      for c in order]
    ypos   = list(range(len(order)))

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    fig.subplots_adjust(left=0.20, right=0.66, top=0.90, bottom=0.16)

    # Family bands: full-MC (Geant4+PHITS) at the bottom, fast engines above.
    ax.axhspan(-0.5, 1.5, alpha=0.07, color="#1565C0", zorder=0, lw=0)
    ax.axhspan( 1.5, 5.5, alpha=0.07, color="#2E7D32", zorder=0, lw=0)

    ax.barh(ypos, times, color=colors, height=0.62,
            edgecolor="white", linewidth=0.8, zorder=2)

    ax.axvline(geant4_t, color=CODE_STYLE["Geant4"]["color"],
               linestyle="--", linewidth=1.3, alpha=0.5, zorder=3)
    ax.text(geant4_t * 0.86, 2.6, "Geant4\nref.", ha="right", va="center",
            fontsize=8, color=CODE_STYLE["Geant4"]["color"],
            alpha=0.8, style="italic")

    ax.set_xscale("log")
    ax.set_xlim(20, 1_500_000)
    ax.set_ylim(-0.5, 5.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([CODE_STYLE[c]["label"] for c in order], fontsize=9.5)
    ax.set_xlabel("CPU time, single-core-equivalent  (log scale)", labelpad=6)

    tick_vals   = [60, 600, 3_600, 36_000, 360_000]   # Bethe-Bloch (48 s) sits left of 1 min
    tick_labels = ["1 min", "10 min", "1 h", "10 h", "100 h"]
    ax.set_xticks(tick_vals)
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(tick_labels))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.tick_params(axis="x", labelsize=9, labelrotation=25, pad=3)
    ax.grid(True, axis="x", which="major", alpha=0.25, linestyle="--")

    ax.set_title("Simulation speed ($6\\times10^{5}$ muons)",
                 fontsize=12, pad=8, loc="left")

    # Right-margin annotation columns: absolute time and speedup vs Geant4.
    ax.text(1.04, 1.02, "Time", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, fontweight="bold", color="0.35", clip_on=False)
    ax.text(1.28, 1.02, "vs Geant4", transform=ax.transAxes, ha="left",
            va="bottom", fontsize=9.5, fontweight="bold", color="0.35",
            clip_on=False)
    from matplotlib.transforms import blended_transform_factory
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for y, c in zip(ypos, order):
        ax.text(1.04, y, fmt_time(timings[c]), transform=trans, ha="left",
                va="center", fontsize=9.5, color=CODE_STYLE[c]["color"],
                fontweight="bold", clip_on=False)
        ax.text(1.28, y, fmt_speedup(timings[c], geant4_t), transform=trans,
                ha="left", va="center", fontsize=9, color="0.30", clip_on=False)

    ax.text(0.0, -0.22, "All codes on one core (wall = CPU); Geant4 and PHITS "
            "timed on 1/32 of the muons and scaled.",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="0.45", style="italic")

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"[ok] wrote {out_pdf}")
    for c in order:
        print(f"     {c:9s} {timings[c]:>9.0f} s   {fmt_speedup(timings[c], geant4_t)}")


def main() -> int:
    if not TIMING_CSV.is_file():
        print(f"[err] timing record not found: {TIMING_CSV}", file=sys.stderr)
        return 1
    make_fig10(load_timings(), OUT_PDF)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
