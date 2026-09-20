#!/usr/bin/env python3
"""
make_fig10_timing.py — CPU time (single-core-equivalent) of the six benchmark
codes on the common 6x10^5-muon Standard-Rock source, rendered as a
horizontal log-scale bar chart (fig10_timing.pdf).

For the single-thread codes this equals wall-clock time; the PHITS value is its
total CPU time summed over 10 OpenMP threads (the parser reads PHITS's
"total cpu time" into phits_timing.txt). The comparison is therefore a fair
single-core-equivalent CPU-time comparison, not a wall-clock one.

Timing data is read from the four-code benchmark folder so the figure stays
in sync with the survival table (tab:survival):

  benchmark/geant4_muon_rock_v5/
    {MUSIC,BB,PROPOSAL,UCMuon}/<code>_timing.txt   "Total : <s> s"   (all depths)
    phits/phits_timing.txt                         "Elapsed : <s> s"

Geant4 does not expose a machine-readable single-value timing file for the
full 6x10^5-muon run (the trimmed outputs/ keeps only a 10^5 sub-run), so its
12.5 h wall time is carried as a documented constant below. The MUSIC, BB,
PROPOSAL, UCMuon and PHITS values are parsed at runtime.

Both Geant4 (1 thread, SerialOnly) and the four fast engines run on a single
core (wall = CPU); PHITS ran with 10 OpenMP threads, and its total CPU time is
8.2x that of Geant4 -- i.e. it consumes 8.2x more CPU, which is the point of the
annotation.

Usage:  python3 manuscript/scripts/make_fig10_timing.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF    = SCRIPT_DIR.parent / "figs" / "fig10_timing.pdf"
BENCH_DIR  = (SCRIPT_DIR.parent.parent
              / "benchmark" / "geant4_muon_rock_v5")

# Geant4 full 6e5-muon run is not recoverable from the trimmed outputs/ folder;
# carried as a documented constant (12.5 h single thread, FTFP_BERT).
GEANT4_SECONDS = 45_134

# Colours match CODE_STYLE in make_fig06.py so the two figures read as a set.
CODE_STYLE = {
    "Geant4":   {"label": "Geant4 11.2 (FTFP_BERT)",      "color": "#1565C0"},
    "PHITS":    {"label": "PHITS 3.36",                   "color": "#C62828"},
    "MUSIC":    {"label": "MUSIC (Engine 2)",             "color": "#FF8F00"},
    "PROPOSAL": {"label": "PROPOSAL (Engine 4)",          "color": "#2E7D32"},
    "BB":       {"label": "Bethe–Bloch CSDA (Engine 3)", "color": "#7B1FA2"},
    "UCMuon":   {"label": "UCMuon-MC (Engine 1)",         "color": "#00838F"},
}


def _grep_seconds(path: Path, key: str) -> float:
    """Return the float following '<key> : <value>' in a timing file."""
    text = path.read_text()
    m = re.search(rf"{key}\s*:\s*([\d.]+)", text)
    if not m:
        raise ValueError(f"no '{key}' entry in {path}")
    return float(m.group(1))


def load_timings() -> dict[str, float]:
    t = {"Geant4": float(GEANT4_SECONDS)}
    for code in ("MUSIC", "BB", "PROPOSAL", "UCMuon"):
        t[code] = _grep_seconds(BENCH_DIR / code / f"{code}_timing.txt", "Total")
    t["PHITS"] = _grep_seconds(BENCH_DIR / "phits" / "phits_timing.txt", "Elapsed")
    return t


def fmt_time(s: float) -> str:
    if s < 60:        return f"{s:.0f} s"
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
    ax.set_xlim(60, 1_500_000)
    ax.set_ylim(-0.5, 5.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([CODE_STYLE[c]["label"] for c in order], fontsize=9.5)
    ax.set_xlabel("CPU time, single-core-equivalent  (log scale)", labelpad=6)

    tick_vals   = [60, 600, 3_600, 36_000, 360_000]
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

    ax.text(0.0, -0.22, "PHITS value is total CPU time over 10 OpenMP threads; "
            "all other codes single-thread (wall = CPU).",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="0.45", style="italic")

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"[ok] wrote {out_pdf}")
    for c in order:
        print(f"     {c:9s} {timings[c]:>9.0f} s   {fmt_speedup(timings[c], geant4_t)}")


def main() -> int:
    if not BENCH_DIR.is_dir():
        print(f"[err] benchmark folder not found: {BENCH_DIR}", file=sys.stderr)
        return 1
    make_fig10(load_timings(), OUT_PDF)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
