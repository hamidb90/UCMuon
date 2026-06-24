#!/usr/bin/env python3
"""Compare the UCMuon source-spectrum parametrizations with original data.

Inputs (see README.md for provenance):
    data/cosmoaleph2013_table1.csv   CosmoALEPH measured spectrum + charge ratio
    data/reyna2006_fig3_data.csv     11 experiments digitized from Reyna Fig. 3
    data/reyna2006_fig3_bestfit.csv  Reyna's published best-fit curve, digitized

Outputs in figures/:
    fig_source_spectrum.{pdf,png}        vertical spectrum + ratio-to-Reyna panel
    fig_source_spectrum_p3.{pdf,png}     p^3-weighted spectrum (high-p detail)
    fig_charge_ratio.{pdf,png}           mu+/mu- charge ratio vs CosmoALEPH
    fig_zenith_dependence.{pdf,png}      zenith-angle dependence of the models
    fig_mode<N>_<name>.{pdf,png}         one figure per spectrum mode: that
                                         model alone vs the data, with a
                                         data/model ratio panel
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ucmuon_spectra as ucm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "legend.fontsize": 7.5,
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
})

EXPERIMENTS = {   # csv key -> (label, marker, colour)
    "nandi_sinha_0deg": ("Nandi & Sinha 0°", "*", "#777777"),
    "mars_0deg":        ("MARS 0°", "P", "#999999"),
    "okayama_0deg":     ("OKAYAMA 0°", "x", "#555555"),
    "okayama_30deg":    ("OKAYAMA 30°", "s", "#bb8888"),
    "kellogg_30deg":    ("Kellogg 30°", "d", "#cc6666"),
    "okayama_60deg":    ("OKAYAMA 60°", "^", "#88aa88"),
    "okayama_75deg":    ("OKAYAMA 75°", "v", "#8888bb"),
    "kellogg_75deg":    ("Kellogg 75°", "+", "#6666cc"),
    "kiel_desy_75deg":  ("Kiel–DESY 75°", "o", "#77aacc"),
    "okayama_80deg":    ("OKAYAMA 80°", "p", "#66bbbb"),
    "mutron_89deg":     ("MUTRON 89°", "D", "#bb88bb"),
}

# All seven source modes in numeric order, so the legend reads 1..7.
# Mode 2 is a sampling shape with no absolute normalization, so it is pinned
# to Mode 7 at 100 GeV/c (every other mode is absolutely normalized); it is
# kept in numeric position here rather than appended at the end.
MODES = [  # (function, label, colour, linestyle, linewidth)
    (ucm.mode1_cosmoaleph, "Mode 1: CosmoALEPH fit", "tab:purple", (0, (5, 2)), 1.4),
    (lambda p: ucm.mode2_powerlaw(p, norm_to=(100.0, ucm.mode7_reyna(100.0))),
     "Mode 2: $E^{-3.7}$ (shape, pinned 100 GeV)", "tab:brown",
     (0, (5, 1, 1, 1)), 1.6),
    (ucm.mode3_parma,      "Mode 3: PARMA/EXPACS (Sato 2015)", "tab:cyan",
     (0, (3, 1, 1, 1)), 1.4),
    (ucm.mode4_guan,       "Mode 4: Guan et al. 2015", "tab:blue", "-", 1.4),
    (ucm.mode5_frosin,     "Mode 5: Frosin et al. 2025", "tab:green", "-", 1.4),
    (ucm.mode6_gaisser,    "Mode 6: Gaisser 1990", "tab:orange", "-", 1.4),
    (ucm.mode7_reyna,      "Mode 7: Reyna 2006", "tab:red", "-", 1.4),
]


def load_compilation():
    data = {}
    with open(os.path.join(DATA, "reyna2006_fig3_data.csv")) as f:
        for row in csv.DictReader(f):
            data.setdefault(row["experiment"], []).append(
                (float(row["zeta_gev"]), float(row["intensity_cm2_s_sr_gev"]),
                 float(row["err_lo"] or 0), float(row["err_hi"] or 0)))
    return {k: np.array(v) for k, v in data.items()}


def load_cosmoaleph():
    rows = []
    with open(os.path.join(DATA, "cosmoaleph2013_table1.csv")) as f:
        for row in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
            rows.append([float(row[k]) for k in
                         ("p_gevc", "flux", "flux_err_stat", "flux_err_sys",
                          "charge_ratio", "cr_err_stat", "cr_err_sys")])
    return np.array(rows)


def load_bestfit():
    pts = []
    with open(os.path.join(DATA, "reyna2006_fig3_bestfit.csv")) as f:
        for row in csv.DictReader(f):
            pts.append((float(row["zeta_gev"]),
                        float(row["intensity_cm2_s_sr_gev"])))
    return np.array(pts)


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{stem}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  figures/{stem}.pdf  figures/{stem}.png")


def plot_compilation(ax, comp, weight=lambda p: 1.0, ms=3.5):
    for key, (label, marker, colour) in EXPERIMENTS.items():
        arr = comp[key]
        w = weight(arr[:, 0])
        ax.errorbar(arr[:, 0], arr[:, 1] * w,
                    yerr=[arr[:, 2] * w, arr[:, 3] * w],
                    fmt=marker, ms=ms, color=colour, mfc="none",
                    elinewidth=0.5, capsize=0, lw=0, label=label, zorder=2)


def plot_cosmoaleph(ax, ca, weight=lambda p: 1.0):
    p = ca[:, 0]
    err = np.hypot(ca[:, 2], ca[:, 3]) * weight(p)
    ax.errorbar(p, ca[:, 1] * weight(p), yerr=err, fmt="o", ms=4,
                color="black", elinewidth=0.8, capsize=1.5, lw=0,
                label="CosmoALEPH (2013)", zorder=3)


def fig_spectrum(comp, ca):
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.4, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [2.6, 1], "hspace": 0.06})

    p = np.geomspace(0.3, 3000, 400)
    plot_compilation(ax, comp)
    plot_cosmoaleph(ax, ca)
    for fn, label, colour, ls, lw in MODES:
        ax.plot(p, fn(p), ls=ls, color=colour, lw=lw, label=label, zorder=4)
    ax.plot(p, ucm.bugaev1998_vertical(p), ls=":", color="black", lw=1.2,
            label="Bugaev et al. 1998 (ref.)", zorder=4)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.3, 3000); ax.set_ylim(1e-11, 6e-3)
    ax.set_ylabel(r"$I_V$  [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV/$c$)$^{-1}$]")
    ax.legend(ncol=2, loc="lower left", framealpha=0.9)
    ax.set_title("UCMuon source spectra vs. sea-level vertical muon data")

    # ratio panel: everything / Reyna 2006
    ref = ucm.mode7_reyna
    for key, (label, marker, colour) in EXPERIMENTS.items():
        arr = comp[key]
        r = ref(arr[:, 0])
        axr.errorbar(arr[:, 0], arr[:, 1] / r,
                     yerr=[arr[:, 2] / r, arr[:, 3] / r],
                     fmt=marker, ms=3, color=colour, mfc="none",
                     elinewidth=0.5, lw=0, zorder=2)
    rca = ref(ca[:, 0])
    axr.errorbar(ca[:, 0], ca[:, 1] / rca,
                 yerr=np.hypot(ca[:, 2], ca[:, 3]) / rca,
                 fmt="o", ms=4, color="black", elinewidth=0.8, capsize=1.5,
                 lw=0, zorder=3)
    for fn, label, colour, ls, lw in MODES:
        axr.plot(p, fn(p) / ref(p), ls=ls, color=colour, lw=lw, zorder=4)
    axr.plot(p, ucm.bugaev1998_vertical(p) / ref(p), ls=":", color="black",
             lw=1.2, zorder=4)
    axr.axhline(1.0, color="tab:red", lw=0.8)
    axr.set_xscale("log")
    axr.set_xlim(0.3, 3000); axr.set_ylim(0.35, 1.85)
    axr.set_xlabel(r"$p_\mu$  or  $\zeta = p_\mu\cos\theta$  [GeV/$c$]")
    axr.set_ylabel("ratio to Mode 7\n(Reyna 2006)")
    save(fig, "fig_source_spectrum")


def fig_spectrum_p3(comp, ca, bestfit):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    w = lambda p: p**3
    plot_compilation(ax, comp, weight=w)
    plot_cosmoaleph(ax, ca, weight=w)
    p = np.geomspace(1, 3000, 400)
    for fn, label, colour, ls, lw in MODES:
        ax.plot(p, fn(p) * p**3, ls=ls, color=colour, lw=lw, label=label)
    ax.plot(p, ucm.bugaev1998_vertical(p) * p**3, ls=":", color="black",
            lw=1.2, label="Bugaev et al. 1998 (ref.)")
    ax.plot(bestfit[:, 0], bestfit[:, 1] * bestfit[:, 0]**3, "s", ms=2,
            mfc="none", color="tab:red", alpha=0.5,
            label="Reyna fit as drawn (digitized)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1, 3000); ax.set_ylim(2e-3, 0.5)
    ax.set_xlabel(r"$p_\mu$  or  $\zeta = p_\mu\cos\theta$  [GeV/$c$]")
    ax.set_ylabel(r"$p^3\, I_V$  [(GeV/$c$)$^2$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax.set_title(r"$p^3$-weighted vertical spectrum")
    ax.legend(ncol=2, loc="lower center", framealpha=0.9)
    save(fig, "fig_source_spectrum_p3")


def fig_single_mode(comp, ca, bestfit, fn, stem, title, note=None,
                    valid=None, show_bestfit=False):
    """One model against all data, with a data/model ratio panel.

    valid: (p_lo, p_hi) momentum range where the parametrization is meant
    to be used; shaded outside.
    """
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.4, 6.4), sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.06})

    p = np.geomspace(0.3, 3000, 400)
    plot_compilation(ax, comp)
    plot_cosmoaleph(ax, ca)
    ax.plot(p, fn(p), color="tab:red", lw=1.8, label=title, zorder=5)
    if show_bestfit:
        ax.plot(bestfit[:, 0], bestfit[:, 1], "s", ms=2.5, mfc="none",
                color="darkred", alpha=0.6, zorder=4,
                label="fit as drawn in the paper (digitized)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.3, 3000); ax.set_ylim(1e-11, 6e-3)
    ax.set_ylabel(r"$I_V$  [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV/$c$)$^{-1}$]")
    ax.set_title(title)
    ax.legend(ncol=2, loc="lower left", framealpha=0.9)
    if note:
        ax.text(0.975, 0.955, note, transform=ax.transAxes, ha="right",
                va="top", fontsize=7.5, style="italic",
                bbox=dict(fc="white", ec="0.7", alpha=0.85))

    # ratio panel: data / model
    for key, (label, marker, colour) in EXPERIMENTS.items():
        arr = comp[key]
        m = fn(arr[:, 0])
        axr.errorbar(arr[:, 0], arr[:, 1] / m,
                     yerr=[arr[:, 2] / m, arr[:, 3] / m],
                     fmt=marker, ms=3, color=colour, mfc="none",
                     elinewidth=0.5, lw=0, zorder=2)
    mca = fn(ca[:, 0])
    axr.errorbar(ca[:, 0], ca[:, 1] / mca,
                 yerr=np.hypot(ca[:, 2], ca[:, 3]) / mca,
                 fmt="o", ms=4, color="black", elinewidth=0.8, capsize=1.5,
                 lw=0, zorder=3)
    axr.axhline(1.0, color="tab:red", lw=1.0)
    axr.axhspan(0.9, 1.1, color="tab:red", alpha=0.10, lw=0)
    axr.set_xscale("log"); axr.set_yscale("log")
    axr.set_xlim(0.3, 3000); axr.set_ylim(0.2, 5)
    axr.set_yticks([0.25, 0.5, 1, 2, 4])
    axr.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
    axr.set_xlabel(r"$p_\mu$  or  $\zeta = p_\mu\cos\theta$  [GeV/$c$]")
    axr.set_ylabel("data / model")

    if valid:
        for a in (ax, axr):
            if valid[0] > 0.3:
                a.axvspan(0.3, valid[0], color="gray", alpha=0.12, lw=0)
            if valid[1] < 3000:
                a.axvspan(valid[1], 3000, color="gray", alpha=0.12, lw=0)
    save(fig, stem)


def fig_charge_ratio(ca):
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    p = np.geomspace(80, 3000, 800)
    ax.plot(p, ucm.charge_ratio(p), drawstyle="steps-post", color="tab:red",
            lw=1.4, label="UCMuon generator (step table)")
    err = np.hypot(ca[:, 5], ca[:, 6])
    ax.errorbar(ca[:, 0], ca[:, 4], yerr=err, fmt="o", ms=4, color="black",
                elinewidth=0.8, capsize=1.5, lw=0,
                label="CosmoALEPH (2013), stat ⊕ sys")
    ax.axhline(1.0, color="gray", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(80, 3000); ax.set_ylim(0.0, 2.4)
    ax.set_xlabel(r"$p_\mu$ [GeV/$c$]")
    ax.set_ylabel(r"$N_{\mu^+}/N_{\mu^-}$")
    ax.set_title("Muon charge ratio: generator table vs. CosmoALEPH")
    ax.legend(loc="lower left")
    save(fig, "fig_charge_ratio")


def fig_zenith():
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), sharey=True)
    theta = np.linspace(0, np.deg2rad(85), 300)
    deg = np.rad2deg(theta)
    for ax, p0 in zip(axes, (1.0, 10.0, 100.0)):
        models = [
            ("PARMA/EXPACS (Mode 3)", "tab:cyan",
             ucm.parma_angular(p0, deg) / ucm.parma_angular(p0, 0.0)),
            ("Guan et al. (Mode 4)", "tab:blue",
             ucm.mode4_guan(p0, np.cos(theta)) / ucm.mode4_guan(p0, 1.0)),
            ("Frosin et al. (Mode 5)", "tab:green",
             ucm.mode5_frosin(p0, np.cos(theta)) / ucm.mode5_frosin(p0, 1.0)),
            ("Gaisser 1990 (Mode 6)", "tab:orange",
             ucm.mode6_gaisser(p0, np.cos(theta)) / ucm.mode6_gaisser(p0, 1.0)),
            ("Reyna $\\cos^3\\!\\theta\\, I_V(p\\cos\\theta)$ (Mode 7)",
             "tab:red",
             ucm.reyna_angular(p0, theta) / ucm.reyna_angular(p0, 0.0)),
        ]
        for label, colour, y in models:
            ax.plot(deg, y, lw=1.4, color=colour, label=label)
        ax.plot(deg, np.cos(theta)**2, ls=":", color="gray", lw=1.0,
                label=r"$\cos^2\theta$ (ref.)")
        ax.set_title(f"$p_\\mu$ = {p0:g} GeV/$c$")
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_xlim(0, 85); ax.set_ylim(0, 1.6)
    axes[0].set_ylabel(r"$I(p,\theta)\;/\;I(p,0)$")
    axes[0].legend(loc="upper left")
    fig.suptitle("Zenith-angle dependence of the UCMuon source models", y=1.02)
    save(fig, "fig_zenith_dependence")


def main():
    os.makedirs(FIGS, exist_ok=True)
    comp = load_compilation()
    ca = load_cosmoaleph()
    bestfit = load_bestfit()
    print("writing figures:")
    fig_spectrum(comp, ca)
    fig_spectrum_p3(comp, ca, bestfit)
    fig_charge_ratio(ca)
    fig_zenith()

    # one figure per spectrum mode
    fig_single_mode(
        comp, ca, bestfit, ucm.mode1_cosmoaleph,
        "fig_mode1_cosmoaleph", "Mode 1: CosmoALEPH power-law fit",
        note="fit range 112–2239 GeV/c;\nextrapolation below is not physical",
        valid=(100, 2500))
    fig_single_mode(
        comp, ca, bestfit,
        lambda p: ucm.mode2_powerlaw(p, norm_to=(100.0, ucm.mode7_reyna(100.0))),
        "fig_mode2_powerlaw",
        "Mode 2: $E^{-3.7}$ power law (Kudryavtsev/MUSIC)",
        note="sampling shape only — no absolute normalization;\n"
             "pinned to Reyna 2006 at 100 GeV/c",
        valid=(50, 3000))
    fig_single_mode(
        comp, ca, bestfit, ucm.mode3_parma,
        "fig_mode3_parma", "Mode 3: PARMA/EXPACS (Sato 2015)",
        note="numerical atmospheric model (sea level, $W=0$);\n"
             "tracks the data to $\\sim$1 TeV/c, then\n"
             "underpredicts the high-energy tail",
        valid=(0.3, 1000))
    fig_single_mode(
        comp, ca, bestfit, ucm.mode4_guan,
        "fig_mode4_guan", "Mode 4: Guan et al. 2015")
    fig_single_mode(
        comp, ca, bestfit, ucm.mode5_frosin,
        "fig_mode5_frosin", "Mode 5: Frosin et al. 2025")
    fig_single_mode(
        comp, ca, bestfit, ucm.mode6_gaisser,
        "fig_mode6_gaisser", "Mode 6: Gaisser 1990 (pion+kaon)",
        note="no muon decay / energy-loss correction:\n"
             "valid only for $E \\gtrsim 100/\\cos\\theta$ GeV",
        valid=(30, 3000))
    fig_single_mode(
        comp, ca, bestfit, ucm.mode7_reyna,
        "fig_mode7_reyna", "Mode 7: Reyna 2006 log-polynomial",
        show_bestfit=True)


if __name__ == "__main__":
    main()
