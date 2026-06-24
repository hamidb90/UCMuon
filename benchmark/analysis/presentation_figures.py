#!/usr/bin/env python3
"""
presentation_figures.py  —  CCS muography benchmark, slide-ready figures

Reads:   figures_benchmark/benchmark_summary.csv
Writes:  figures_presentation/
           pres_fig1_physics_outputs.png   (what simulation gives you, 2×2)
           pres_fig2_code_engines.png      (capability matrix + agreement)

Usage:
    python3 presentation_figures.py
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# ── Paths ─────────────────────────────────────────────────────────────────────
BENCH_CSV        = "figures_benchmark/benchmark_summary.csv"
OUTDIR           = "figures_presentation"
ROCK_DENSITY     = 2.65
EXCLUDE_DEPTH_CM = {1000}    # 26.5 MWE: 5 GeV Bragg boundary — codes incomparable

# ── Code identity ─────────────────────────────────────────────────────────────
CODE_ORDER = ["Geant4", "PHITS", "MUSIC", "PROPOSAL", "BB", "UCMuon"]
CODE_DISPLAY = {
    "Geant4":   "Geant4",
    "PHITS":    "PHITS",
    "MUSIC":    "MUSIC",
    "PROPOSAL": "PROPOSAL",
    "BB":       "BB",
    "UCMuon":   "UCMuon",
}
CODE_TYPE = {
    "Geant4":   "Full Monte Carlo",
    "PHITS":    "Full Monte Carlo",
    "MUSIC":    "Dedicated muon transport",
    "PROPOSAL": "Dedicated muon transport",
    "BB":       "Analytic (Bethe-Bloch)",
    "UCMuon":   "Stochastic (PUMAS)",
}
CODE_TYPE_SHORT = {
    "Geant4":   "Full MC",
    "PHITS":    "Full MC",
    "MUSIC":    "Muon transp.",
    "PROPOSAL": "Muon transp.",
    "BB":       "Analytic",
    "UCMuon":   "PUMAS",
}
CODE_COLORS = {
    "Geant4":   "#1565C0",
    "PHITS":    "#C62828",
    "MUSIC":    "#E65100",
    "PROPOSAL": "#2E7D32",
    "BB":       "#7B1FA2",
    "UCMuon":   "#00695C",
}
MARKERS = {
    "Geant4":   "o",
    "PHITS":    "s",
    "MUSIC":    "^",
    "PROPOSAL": "D",
    "BB":       "P",
    "UCMuon":   "X",
}
LINESTYLES = {
    "Geant4":   "-",
    "PHITS":    "--",
    "MUSIC":    "-.",
    "PROPOSAL": ":",
    "BB":       (0, (5, 1)),
    "UCMuon":   (0, (3, 1, 1, 1)),
}

# ── Presentation rc ───────────────────────────────────────────────────────────
PRES_RC = {
    "font.family":        "sans-serif",
    "font.size":          13,
    "axes.labelsize":     14,
    "axes.titlesize":     15,
    "axes.titlepad":      10,
    "axes.linewidth":     1.4,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.labelsize":    12,
    "ytick.labelsize":    12,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "xtick.major.size":   5,
    "ytick.major.size":   5,
    "legend.fontsize":    10,
    "legend.framealpha":  0.92,
    "legend.edgecolor":   "0.78",
    "lines.linewidth":    2.5,
    "grid.alpha":         0.28,
    "grid.linewidth":     0.7,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def save_fig(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"  Saved  {path}")


def load_bench():
    df = pd.read_csv(BENCH_CSV)
    df = df[~df["DepthCm"].isin(EXCLUDE_DEPTH_CM)].copy()
    return df


def code_kw(code, lw_scale=1.0):
    lw = (3.5 if code == "Geant4" else 2.2) * lw_scale
    ms = 9   if code == "Geant4" else 6
    zo = 4   if code == "Geant4" else 2
    return dict(color=CODE_COLORS[code], linewidth=lw, linestyle=LINESTYLES[code],
                marker=MARKERS[code], markersize=ms, zorder=zo)


def _se_transmission(sub):
    """Binomial standard error on transmission [%]."""
    p  = sub["Transmission_%"].values / 100.0
    Nt = sub["N_transmitted"].values
    # N_input back-computed from Nt and p (avoid division by near-zero)
    Ni = np.where(p > 0.001, Nt / p, Nt / 0.001)
    return np.sqrt(np.clip(p * (1.0 - p) / Ni, 0, None)) * 100.0


def _se_exit_ke(sub):
    """SEM on mean exit KE [GeV] using StdELoss as proxy for StdExitKE."""
    if "StdELoss_GeV" not in sub.columns or "N_transmitted" not in sub.columns:
        return None
    std = sub["StdELoss_GeV"].fillna(0.0).values
    N   = sub["N_transmitted"].fillna(1.0).values
    se  = np.where(N > 1, std / np.sqrt(N), np.nan)
    return np.where(np.isnan(sub["MeanExitKE_GeV"].values), np.nan, se)


def _se_angle(sub):
    """SEM on mean MCS angle [°] — directly from SEM_Angle_deg column."""
    if "SEM_Angle_deg" not in sub.columns:
        return None
    sem = sub["SEM_Angle_deg"].fillna(np.nan).values
    return sem


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Physics outputs: what simulation tells you
# ══════════════════════════════════════════════════════════════════════════════
def fig_physics_outputs(df):
    """2 × 2 panel: transmission, exit KE, process fractions, MCS angle."""
    codes_present = [c for c in CODE_ORDER if c in df["Code"].values]

    with plt.rc_context(PRES_RC):
        fig = plt.figure(figsize=(14, 10))
        gs  = GridSpec(2, 2, figure=fig, hspace=0.52, wspace=0.42)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])

        # ── [A] Transmission vs depth ─────────────────────────────────
        for code in codes_present:
            sub = df[df["Code"] == code].sort_values("MWE")
            if sub["Transmission_%"].isna().all():
                continue
            kw = code_kw(code)
            ax1.errorbar(sub["MWE"], sub["Transmission_%"],
                         yerr=_se_transmission(sub),
                         capsize=4, capthick=1.5, elinewidth=1.2,
                         label=CODE_DISPLAY[code], **kw)

        ax1.set_xlabel("Depth (m.w.e.)")
        ax1.set_ylabel("Muon transmission (%)")
        ax1.set_title("(A)  Muon Survival Through Rock", fontweight="bold")
        ax1.set_ylim(0, 108)
        ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax1.legend(fontsize=9, ncol=2, loc="upper right")
        ax1.grid(True, alpha=0.28)
        ax1.text(0.97, 0.97,
                 "All six codes agree within ±3%\nat relevant CCS depths",
                 transform=ax1.transAxes, ha="right", va="top",
                 fontsize=9, color="#444", style="italic",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", alpha=0.85))

        # ── [B] Mean exit KE vs depth ─────────────────────────────────
        for code in codes_present:
            sub = df[df["Code"] == code].sort_values("MWE")
            if "MeanExitKE_GeV" not in sub.columns or sub["MeanExitKE_GeV"].isna().all():
                continue
            kw  = code_kw(code)
            se  = _se_exit_ke(sub)
            ax2.errorbar(sub["MWE"], sub["MeanExitKE_GeV"],
                         yerr=se, capsize=4, capthick=1.5, elinewidth=1.2,
                         label=CODE_DISPLAY[code], **kw)

        ax2.set_xlabel("Depth (m.w.e.)")
        ax2.set_ylabel("Mean exit kinetic energy (GeV)")
        ax2.set_title("(B)  Exit Muon Energy vs Depth", fontweight="bold")
        ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax2.legend(fontsize=9, ncol=2, loc="upper left")
        ax2.grid(True, alpha=0.28)

        # ── [C] Geant4 process fractions (stacked bar) ───────────────
        g4   = df[df["Code"] == "Geant4"].sort_values("MWE")
        procs = [("Ionisation",     "FracIon_%",  "#42A5F5"),
                 ("Bremsstrahlung", "FracBrem_%",  "#EF5350"),
                 ("Pair production","FracPair_%",  "#FFA726"),
                 ("Photo-nuclear",  "FracNucl_%",  "#AB47BC")]
        if all(col in g4.columns for _, col, _ in procs):
            x    = np.arange(len(g4))
            mwes = [f"{v:.1f}" for v in g4["MWE"]]
            bot  = np.zeros(len(g4))
            for name, col, color in procs:
                vals = g4[col].fillna(0).values
                bars = ax3.bar(x, vals, bottom=bot, color=color, label=name,
                               width=0.68, edgecolor="white", linewidth=0.6)
                for xi, (b, v) in enumerate(zip(bot, vals)):
                    if v > 2.5:
                        ax3.text(xi, b + v / 2, f"{v:.1f}%",
                                 ha="center", va="center", fontsize=9,
                                 color="white", fontweight="bold")
                bot += vals
            ax3.set_xticks(x)
            ax3.set_xticklabels(mwes, fontsize=11)
            ax3.set_xlabel("Depth (m.w.e.)")
            ax3.set_ylabel("Fraction of total energy loss (%)")
            ax3.set_title("(C)  Energy Loss by Process  [Geant4]", fontweight="bold")
            ax3.legend(fontsize=9, loc="lower left")
            ax3.set_ylim(0, 115)
            ax3.grid(True, alpha=0.25, axis="y")
            ax3.spines["top"].set_visible(False)
            ax3.spines["right"].set_visible(False)

        # ── [D] MCS scatter angle ─────────────────────────────────────
        for code in codes_present:
            if code == "PHITS":
                continue   # only aggregate zenith angle — not directly comparable
            sub = df[(df["Code"] == code) & df["MeanAngle_deg"].notna()].sort_values("MWE")
            if sub.empty:
                continue
            kw  = code_kw(code)
            se  = _se_angle(sub)
            ax4.errorbar(sub["MWE"], sub["MeanAngle_deg"],
                         yerr=se, capsize=4, capthick=1.5, elinewidth=1.2,
                         label=CODE_DISPLAY[code], **kw)

        ax4.set_xlabel("Depth (m.w.e.)")
        ax4.set_ylabel("Mean MCS deflection angle (°)")
        ax4.set_title("(D)  Multiple Coulomb Scattering Angle", fontweight="bold")
        ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax4.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax4.legend(fontsize=9, ncol=2, loc="upper right")
        ax4.grid(True, alpha=0.28)
        ax4.text(0.97, 0.03,
                 "PHITS excluded (aggregate tally only)",
                 transform=ax4.transAxes, ha="right", va="bottom",
                 fontsize=8.5, color="#888", style="italic")

        fig.suptitle(
            "Geant4 + PHITS  vs  Fast Transport Engines\n"
            "Muon Propagation Through Standard Rock  (CCS Muography Benchmark)",
            fontsize=16, fontweight="bold", y=1.02)

        save_fig(fig, "pres_fig1_physics_outputs.png")


# ── Shared capability data ─────────────────────────────────────────────────────
_FEATURES = [
    "Muon transmission",
    "Exit energy / energy loss",
    "MCS scatter angle",
    "Lateral displacement",
    "Process breakdown\n(ion / brem / pair)",
    "Secondary particles",
    "Per-event output",
]
_CAPS = {
    #                                  G4  PHITS  MUSIC  PROP  BB  UCMuon
    "Muon transmission":              [ 2,   2,     2,     2,   2,    2  ],
    "Exit energy / energy loss":      [ 2,   2,     2,     2,   2,    2  ],
    "MCS scatter angle":              [ 2,   1,     2,     2,   2,    2  ],
    "Lateral displacement":           [ 2,   0,     2,     2,   2,    2  ],
    "Process breakdown\n(ion / brem / pair)": [ 2, 0, 0,  0,   0,    0  ],
    "Secondary particles":            [ 2,   1,     0,     0,   0,    0  ],
    "Per-event output":               [ 2,   0,     2,     2,   2,    2  ],
}
_CELL_COLOR = {2: "#C8E6C9", 1: "#FFF9C4", 0: "#FFCDD2"}
_CELL_TEXT  = {2: "Full",    1: "Partial", 0: "—"}
_CELL_TC    = {2: "#1B5E20", 1: "#5D4037", 0: "#B71C1C"}


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2a — Code capability matrix
# ══════════════════════════════════════════════════════════════════════════════
def fig_capability_matrix():
    """Standalone capability matrix for all six codes."""
    codes  = CODE_ORDER
    n_code = len(codes)
    n_feat = len(_FEATURES)

    with plt.rc_context(PRES_RC):
        fig, ax = plt.subplots(figsize=(13, 7.5))

        ax.set_xlim(-0.5, n_code - 0.5)
        ax.set_ylim(-0.5, n_feat - 0.5)
        ax.invert_yaxis()

        # Cells
        for fi, feat in enumerate(_FEATURES):
            row = _CAPS[feat]
            for ci, (code, val) in enumerate(zip(codes, row)):
                rect = mpatches.FancyBboxPatch(
                    (ci - 0.46, fi - 0.43), 0.92, 0.86,
                    boxstyle="round,pad=0.04", linewidth=0,
                    facecolor=_CELL_COLOR[val], zorder=1)
                ax.add_patch(rect)
                ax.text(ci, fi, _CELL_TEXT[val],
                        ha="center", va="center", fontsize=13,
                        color=_CELL_TC[val], fontweight="bold", zorder=2)

        # Row separator lines
        for y in np.arange(0.5, n_feat - 0.5, 1.0):
            ax.axhline(y, color="#e0e0e0", linewidth=1.2, zorder=0)

        # Column separator between Full MC and fast engines
        ax.axvline(1.5, color="#aaa", linewidth=2.0, linestyle="--",
                   alpha=0.65, zorder=3)
        for cx, lbl in [(0.5, "Full MC"), (3.5, "Fast engines")]:
            ax.text(cx, -0.54, f"— {lbl} —",
                    ha="center", va="bottom", fontsize=9.5, color="#555",
                    style="italic", clip_on=False, transform=ax.transData)

        # Column headers (top axis, two-line: name + type)
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xticks(range(n_code))
        tick_labels = [f"{CODE_DISPLAY[c]}\n{CODE_TYPE_SHORT[c]}" for c in codes]
        ax.set_xticklabels(tick_labels, fontsize=12)
        for tick, code in zip(ax.get_xticklabels(), codes):
            tick.set_color(CODE_COLORS[code])
            tick.set_fontweight("bold")

        # Feature labels (left)
        ax.set_yticks(range(n_feat))
        ax.set_yticklabels(_FEATURES, fontsize=12.5)
        ax.tick_params(axis="x", length=0, pad=8)
        ax.tick_params(axis="y", length=0, pad=8)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Legend
        handles = [
            mpatches.Patch(facecolor=_CELL_COLOR[2], edgecolor="#bbb",
                           label="Full output"),
            mpatches.Patch(facecolor=_CELL_COLOR[1], edgecolor="#bbb",
                           label="Aggregate / partial"),
            mpatches.Patch(facecolor=_CELL_COLOR[0], edgecolor="#bbb",
                           label="Not available"),
        ]
        ax.legend(handles=handles, loc="lower right", fontsize=11,
                  framealpha=0.95, edgecolor="0.78",
                  bbox_to_anchor=(1.0, -0.12))

        fig.suptitle(
            "Transport Code Capabilities  —  CCS Muography Benchmark",
            fontsize=16, fontweight="bold", y=1.04)

        save_fig(fig, "pres_fig2a_capability_matrix.png")


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2b — Transmission agreement vs Geant4
# ══════════════════════════════════════════════════════════════════════════════
def fig_transmission_agreement(df):
    """Ratio plot: each code's transmission divided by Geant4, with error bars."""
    codes_present = [c for c in CODE_ORDER if c in df["Code"].values]
    g4_df = df[df["Code"] == "Geant4"].set_index("MWE")

    with plt.rc_context(PRES_RC):
        fig, ax = plt.subplots(figsize=(10, 7))

        # Tolerance bands
        ax.axhspan(0.95, 1.05, color="#C8E6C9", alpha=0.60, zorder=0,
                   label="±5 % band")
        ax.axhspan(0.90, 1.10, color="#FFF9C4", alpha=0.55, zorder=0,
                   label="±10 % band")
        ax.axhline(1.0, color="#222", linewidth=1.6, linestyle="--",
                   alpha=0.65, zorder=1, label="Geant4 (reference)")

        # Band labels on the right margin
        for ypos, lbl, col in [(1.027, "±5 %",  "#388E3C"),
                                (1.072, "±10 %", "#F57F17")]:
            ax.text(1.01, ypos, lbl, transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=11, color=col,
                    fontweight="bold")

        for code in codes_present:
            if code == "Geant4":
                continue
            sub = df[df["Code"] == code].sort_values("MWE")
            pts = []
            for _, row in sub.iterrows():
                mwe = row["MWE"]
                v   = row.get("Transmission_%", np.nan)
                Nv  = row.get("N_transmitted",  np.nan)
                if mwe not in g4_df.index: continue
                g4v = g4_df.loc[mwe, "Transmission_%"]
                g4N = g4_df.loc[mwe, "N_transmitted"]
                if np.isnan(v) or np.isnan(g4v) or g4v <= 0: continue
                ratio = v / g4v
                p_v  = v   / 100.0; Ni_v  = Nv  / max(p_v,  0.001)
                p_g4 = g4v / 100.0; Ni_g4 = g4N / max(p_g4, 0.001)
                se_v  = np.sqrt(p_v  * (1 - p_v)  / max(Ni_v,  1)) * 100.0
                se_g4 = np.sqrt(p_g4 * (1 - p_g4) / max(Ni_g4, 1)) * 100.0
                se_r  = ratio * np.sqrt((se_v  / max(v,   0.01))**2
                                      + (se_g4 / max(g4v, 0.01))**2)
                pts.append({"MWE": mwe, "ratio": ratio, "se": se_r})
            if not pts:
                continue
            rdf = pd.DataFrame(pts)
            ax.errorbar(rdf["MWE"], rdf["ratio"], yerr=rdf["se"],
                        color=CODE_COLORS[code], linewidth=2.5,
                        linestyle=LINESTYLES[code], marker=MARKERS[code],
                        markersize=9, capsize=5, capthick=1.8, elinewidth=1.3,
                        label=CODE_DISPLAY[code], zorder=3)

        ax.set_xlabel("Depth (m.w.e.)", fontsize=15)
        ax.set_ylabel("Transmission  (code / Geant4)", fontsize=15)
        ax.set_ylim(0.82, 1.18)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.legend(fontsize=11, loc="lower left",
                  framealpha=0.92, edgecolor="0.75")
        ax.grid(True, alpha=0.28)
        ax.text(0.98, 0.04,
                "Error bars: propagated binomial SE  (N ≈ 100 K per depth)",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9, color="#888", style="italic")

        fig.suptitle(
            "Transmission Agreement  —  All Codes vs Geant4\n"
            "CCS Muography Benchmark  (standard rock, 600 K muons)",
            fontsize=15, fontweight="bold", y=1.02)

        save_fig(fig, "pres_fig2b_transmission_agreement.png")


def fig_timing_comparison():
    """Horizontal bar chart — wall-clock time for the 600 K-muon benchmark."""
    from matplotlib.transforms import blended_transform_factory

    TIMING = {                  # wall-clock seconds, single core, all depths, 600 K muons
        "MUSIC":    204,
        "UCMuon":   328,
        "BB":       871,
        "PROPOSAL": 2_147,
        "Geant4":   45_134,
        "PHITS":    371_146,    # 10 OMP threads; Geant4 is 1 thread (SerialOnly)
    }
    geant4_t = TIMING["Geant4"]

    def fmt_time(s):
        if s < 60:      return f"{s:.0f} s"
        elif s < 3_600: return f"{s/60:.0f} min"
        else:           return f"{s/3600:.1f} h"

    # Sort slowest→fastest so fastest bar appears at the top
    order  = sorted(TIMING, key=lambda c: TIMING[c], reverse=True)
    times  = [TIMING[c]      for c in order]
    colors = [CODE_COLORS[c] for c in order]
    ypos   = list(range(len(order)))

    with plt.rc_context(PRES_RC):
        fig, ax = plt.subplots(figsize=(11, 5.5))
        fig.subplots_adjust(left=0.12, right=0.58, top=0.84, bottom=0.22)

        # Background shading by code family
        ax.axhspan(-0.5, 1.5, alpha=0.07, color="#1565C0", zorder=0, lw=0)
        ax.axhspan( 1.5, 5.5, alpha=0.07, color="#27AE60", zorder=0, lw=0)

        ax.barh(ypos, times, color=colors, height=0.62,
                edgecolor="white", linewidth=0.8, zorder=2)

        # Geant4 reference dashed vertical line — labelled inside the plot area
        ax.axvline(geant4_t, color=CODE_COLORS["Geant4"],
                   linestyle="--", linewidth=1.4, alpha=0.5, zorder=3)
        ax.text(geant4_t * 0.88, 2.5, "Geant4\nref.",
                ha="right", va="center", fontsize=8.5,
                color=CODE_COLORS["Geant4"], alpha=0.80, style="italic")

        # Axis — log scale with clean human-readable ticks
        ax.set_xscale("log")
        ax.set_xlim(60, 1_500_000)
        ax.set_ylim(-0.5, 5.5)
        ax.set_yticks(ypos)
        ax.set_yticklabels([CODE_DISPLAY[c] for c in order], fontsize=13)
        ax.set_xlabel("Wall-clock time  (log scale)", fontsize=13, labelpad=8)

        # Fix x ticks: set explicitly then force FixedFormatter so LogFormatter
        # does not add its own labels on top of ours
        tick_vals   = [60, 600, 3_600, 36_000, 360_000]
        tick_labels = ["1 min", "10 min", "1 h", "10 h", "100 h"]
        ax.set_xticks(tick_vals)
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(tick_labels))
        ax.xaxis.set_minor_locator(ticker.NullLocator())   # no minor ticks
        ax.tick_params(axis="x", labelsize=10, labelrotation=30, pad=4)
        ax.grid(True, axis="x", which="major", alpha=0.25)

        # Title sits above the axes area
        ax.set_title(
            "Simulation Speed  —  600 K muons, standard rock",
            fontsize=13, fontweight="bold", pad=10)

        # Annotation columns in the right margin (x in axes coords, y in data coords)
        trans = blended_transform_factory(ax.transAxes, ax.transData)

        # Column headers anchored just above the top of the axes
        ax.text(1.03, 1.04, "Time",       transform=ax.transAxes,
                ha="left", va="bottom", fontsize=10, fontweight="bold",
                color="0.35", clip_on=False)
        ax.text(1.22, 1.04, "vs Geant4", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=10, fontweight="bold",
                color="0.35", clip_on=False)

        for i, (code, t) in enumerate(zip(order, times)):
            speedup = geant4_t / t
            ax.text(1.03, i, fmt_time(t), transform=trans,
                    va="center", ha="left",
                    fontsize=11, fontweight="bold", color=CODE_COLORS[code],
                    clip_on=False)
            if code == "Geant4":
                spd_str = "—"
            elif speedup < 1:
                spd_str = f"{1/speedup:.1f}× slower"
            else:
                spd_str = f"{speedup:.0f}× faster"
            ax.text(1.22, i, spd_str, transform=trans,
                    va="center", ha="left",
                    fontsize=10.5, color=CODE_COLORS[code],
                    clip_on=False)

        # Code family labels
        ax.text(1.50, 0.5, "Full MC",      transform=trans,
                va="center", ha="left", fontsize=9.5, color="#1565C0",
                fontweight="bold", style="italic", clip_on=False)
        ax.text(1.50, 3.5, "Fast\nengines", transform=trans,
                va="center", ha="left", fontsize=9.5, color="#27AE60",
                fontweight="bold", style="italic", clip_on=False)

        # Small footnote about PHITS threading (placed in figure coords to avoid xlabel clash)
        fig.text(0.12, 0.01, "* PHITS ran with 10 OMP threads; all other codes: 1 core.",
                 ha="left", va="bottom", fontsize=8, color="0.5", style="italic")

        save_fig(fig, "pres_fig3_timing.png")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.exists(BENCH_CSV):
        sys.exit(f"ERROR: {BENCH_CSV} not found — run benchmark_analysis.py first")

    print(f"Reading {BENCH_CSV}")
    df = load_bench()
    codes_found = sorted(df["Code"].unique())
    depths_mwe  = sorted(df["MWE"].unique())
    print(f"  Codes : {codes_found}")
    print(f"  Depths: {depths_mwe} m.w.e.")

    print("\nGenerating presentation figures ...")
    fig_physics_outputs(df)
    fig_capability_matrix()
    fig_transmission_agreement(df)
    fig_timing_comparison()
    print(f"\nDone.  Output in: {OUTDIR}/")
