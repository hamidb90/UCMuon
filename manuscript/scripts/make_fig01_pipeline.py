#!/usr/bin/env python3
"""
make_fig01_pipeline.py — Generate fig01_pipeline.pdf: the UCMuon three-stage
simulation pipeline schematic, enriched with a real output thumbnail under
each stage (surface energy spectrum, survival-vs-depth curve, GUI screenshot).

Thumbnail data sources (all already in the repo):
  Stage 1  surface vertical intensity from ucmuon_flux_models.reyna_flux
  Stage 2  survival vs. depth from
           benchmark/geant4_muon_rock_v5/figures_benchmark/benchmark_summary.csv
  Stage 3  manuscript/figs/Screenshot_gui.png  (manual GUI capture)

Usage:  python3 manuscript/scripts/make_fig01_pipeline.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

SCRIPT_DIR = Path(__file__).resolve().parent
FIGS_DIR = SCRIPT_DIR.parent / "figs"
OUT_PDF = FIGS_DIR / "fig01_pipeline.pdf"
SURVIVAL_CSV = (SCRIPT_DIR.parent.parent / "benchmark" / "geant4_muon_rock_v5"
                / "figures_benchmark" / "benchmark_summary.csv")
GUI_PNG = FIGS_DIR / "Screenshot_gui.png"

sys.path.insert(0, str(SCRIPT_DIR))
import ucmuon_flux_models as fm  # noqa: E402

# Stage colours match make_fig06.py so the figure set reads as one family.
C1, C2, C3 = "#1565C0", "#C62828", "#2E7D32"

# Column centres (figure fraction) and box geometry.
CX = [0.16, 0.50, 0.84]
BOX_W, BOX_H = 0.225, 0.24
BOX_BOT, BOX_TOP = 0.60, 0.84

STAGES = [
    dict(color=C1, title="Stage 1 — Surface generator",
         body="ucmuon_gen_omp\nFortran 90 + OpenMP / MPI",
         items="8 spectra · 4 angular · 3 geometries",
         cap="Surface energy spectrum"),
    dict(color=C2, title="Stage 2 — Transport engine",
         body="UCMuon-MC | MUSIC | Bethe–Bloch\nPROPOSAL | BackMC | Terrain",
         items="slab / DEM overburden · energy loss + MCS",
         cap="Survival vs. depth"),
    dict(color=C3, title="Stage 3 — Analysis / GUI",
         body="Streamlit GUI or CLI scripts\nPython 3.9+",
         items="spectra · rates · transmission & density maps",
         cap="Streamlit GUI — Results tab"),
]
ARROWS = [("muons_surface.dat", "13/14-col ASCII"),
          ("muons_underground.dat", "18-col ASCII")]


def load_survival(preferred=("MUSIC", "UCMuon", "Geant4")):
    """Return (mwe, transmission%) for the first available engine, plus
    Geant4 as a faint reference if present."""
    if not SURVIVAL_CSV.exists():
        return None
    data: dict[str, list[tuple[float, float]]] = {}
    with open(SURVIVAL_CSV) as fh:
        for row in csv.DictReader(fh):
            try:
                mwe = float(row["MWE"]); tr = float(row["Transmission_%"])
            except (KeyError, ValueError):
                continue
            data.setdefault(row["Code"], []).append((mwe, tr))
    if not data:
        return None
    code = next((c for c in preferred if c in data), next(iter(data)))
    pts = sorted(data[code])
    out = {"code": code,
           "mwe": np.array([p[0] for p in pts]),
           "tr": np.array([p[1] for p in pts])}
    if "Geant4" in data and code != "Geant4":
        g = sorted(data["Geant4"])
        out["g4_mwe"] = np.array([p[0] for p in g])
        out["g4_tr"] = np.array([p[1] for p in g])
    return out


def thumb_spectrum(ax):
    p = np.geomspace(1.0, 2000.0, 200)
    iv = fm.reyna_flux(p, 1.0)
    ax.plot(p, iv, color=C1, lw=1.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$p_\mu$ [GeV/$c$]", fontsize=7, labelpad=1)
    ax.set_ylabel(r"$I_V$", fontsize=7, labelpad=1)
    ax.tick_params(labelsize=6, length=2, pad=1)
    ax.set_xticks([1, 100, 1000])
    ax.set_xticklabels(["1", "$10^2$", "$10^3$"])


def thumb_survival(ax):
    s = load_survival()
    if s is None:
        ax.text(0.5, 0.5, "survival data\nnot found", ha="center", va="center",
                fontsize=7, color="0.5", transform=ax.transAxes)
        return
    if "g4_mwe" in s:
        ax.plot(s["g4_mwe"], s["g4_tr"], color="0.55", lw=1.0, ls="--",
                marker="o", ms=2.5, label="Geant4", zorder=2)
    ax.plot(s["mwe"], s["tr"], color=C2, lw=1.8, marker="s", ms=3,
            label=s["code"], zorder=3)
    ax.set_yscale("log")
    ax.set_xlabel("depth [m.w.e.]", fontsize=7, labelpad=1)
    ax.set_ylabel("survival [%]", fontsize=7, labelpad=1)
    ax.tick_params(labelsize=6, length=2, pad=1)
    ax.legend(fontsize=5.5, loc="upper right", handlelength=1.3,
              borderpad=0.3, labelspacing=0.2, framealpha=0.85)


def thumb_screenshot(ax):
    if GUI_PNG.exists():
        ax.imshow(mpimg.imread(str(GUI_PNG)), aspect="equal")
    else:
        ax.text(0.5, 0.5, "Screenshot_gui.png\nnot found", ha="center",
                va="center", fontsize=7, color="0.5", transform=ax.transAxes)
    ax.axis("off")
    for spine in ax.spines.values():
        spine.set_visible(False)


def stage_box(ax, cx, color, title, body, items, **_):
    left = cx - BOX_W / 2
    for fc, ec, lw, alpha in [(color, "none", 0, 0.12), ("none", color, 1.6, 1)]:
        ax.add_patch(FancyBboxPatch((left, BOX_BOT), BOX_W, BOX_H,
                                    boxstyle="round,pad=0.012",
                                    fc=fc, ec=ec, lw=lw, alpha=alpha,
                                    transform=ax.transAxes, zorder=2))
    ax.text(cx, BOX_TOP - 0.045, title, transform=ax.transAxes,
            ha="center", va="center", fontsize=9, fontweight="bold", color=color)
    ax.text(cx, BOX_BOT + BOX_H / 2 - 0.005, body, transform=ax.transAxes,
            ha="center", va="center", fontsize=7.6, family="monospace")
    ax.text(cx, BOX_BOT + 0.028, items, transform=ax.transAxes,
            ha="center", va="center", fontsize=6.8, color="0.30")


def main():
    fig = plt.figure(figsize=(11.0, 6.2))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    for cx, s in zip(CX, STAGES):
        stage_box(ax, cx, **s)

    # Forward data-flow arrows with file-format labels.
    y_mid = (BOX_BOT + BOX_TOP) / 2
    for i, (fname, fmt) in enumerate(ARROWS):
        x1 = CX[i] + BOX_W / 2 + 0.005
        x2 = CX[i + 1] - BOX_W / 2 - 0.005
        ax.add_patch(FancyArrowPatch((x1, y_mid), (x2, y_mid),
                                     arrowstyle="-|>", mutation_scale=15,
                                     lw=1.5, color="0.15", transform=ax.transAxes))
        ax.text((x1 + x2) / 2, y_mid + 0.028, fname, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=6.2, family="monospace")
        ax.text((x1 + x2) / 2, y_mid - 0.028, fmt, transform=ax.transAxes,
                ha="center", va="top", fontsize=6.2, color="0.35")

    # Optional detector acceptance cut (dashed bar over the top).
    y_top, y_bar = BOX_TOP + 0.02, BOX_TOP + 0.085
    ax.plot([CX[2], CX[2], CX[0]], [y_top, y_bar, y_bar], ls="--", lw=1.2,
            color="0.4", transform=ax.transAxes)
    ax.add_patch(FancyArrowPatch((CX[0], y_bar), (CX[0], y_top),
                                 arrowstyle="-|>", mutation_scale=13, lw=1.2,
                                 ls="--", color="0.4", transform=ax.transAxes))
    ax.text((CX[0] + CX[2]) / 2, y_bar + 0.015,
            "Detector acceptance cut (optional): geometry module discards muons "
            "that cannot reach the detector", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8, color="0.4", style="italic")

    # Thumbnails + connectors + captions.
    th_bot, th_h = 0.13, 0.36
    th_top = th_bot + th_h
    builders = [thumb_spectrum, thumb_survival, thumb_screenshot]
    for cx, build, s in zip(CX, builders, STAGES):
        ax.add_patch(FancyArrowPatch((cx, BOX_BOT - 0.005), (cx, th_top + 0.01),
                                     arrowstyle="-|>", mutation_scale=12,
                                     lw=1.1, color="0.45", transform=ax.transAxes))
        tax = fig.add_axes([cx - BOX_W / 2, th_bot, BOX_W, th_h])
        build(tax)
        ax.text(cx, th_bot - 0.035, s["cap"], transform=ax.transAxes,
                ha="center", va="top", fontsize=8, style="italic",
                color=s["color"])

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    print(f"[OK] wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
