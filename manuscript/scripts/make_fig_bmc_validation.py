#!/usr/bin/env python3
"""make_fig_bmc_validation.py: transmission of the backward CSDA flux
integrator (Engine 5) relative to the forward engines, against slant opacity.
(a) without the forward calibration (mode 1), at the calibration's nodes
(benchmark/results/backward_validation_20261019.csv); (b) with it (mode 2,
the default), at zenith angles and opacities that are not nodes
(benchmark/results/backward_validation_offnode_20261020.csv). Both written by
benchmark/analysis/validate_backward.py.

Usage:  python3 manuscript/scripts/make_fig_bmc_validation.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "benchmark" / "results"
OUT = REPO / "manuscript" / "figs" / "fig_bmc_validation.pdf"
MARK = {"UCMuon-MC": "o", "MUSIC": "s", "PROPOSAL": "^"}
CMAP = plt.get_cmap("viridis")


def load(path):
    rows = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            rows[float(r["theta_deg"])].append(r)
    return rows


def panel(ax, rows, column, csda, title):
    thetas = sorted(rows)
    for k, th in enumerate(thetas):
        rs = rows[th]
        col = CMAP(0.9 * k / max(len(thetas) - 1, 1))
        shift = 1.0 + 0.012 * (k - (len(thetas) - 1) / 2)    # spread the zeniths apart
        for eng, m in MARK.items():
            sel = [r for r in rs if r["engine"] == eng]
            mwe = np.array([float(r["mwe"]) for r in sel]) * shift
            tb = np.array([float(r[column]) for r in sel])
            tf = np.array([float(r["T_forward"]) for r in sel])
            err = np.array([float(r["T_forward_err"]) for r in sel]) / tf * tb / tf
            ax.errorbar(mwe, tb / tf, yerr=err, fmt=m, ms=4, color=col,
                        mfc=col if eng == "UCMuon-MC" else "none", lw=0.8, capsize=0,
                        label=f"{th:g}°" if eng == "UCMuon-MC" else None)
        if csda:
            ref = [r for r in rs if r["engine"] == "UCMuon-MC"]
            ax.plot(np.array([float(r["mwe"]) for r in ref]) * shift,
                    [float(r["T_backward_csda"]) / np.mean([float(q["T_forward"]) for q in rs
                                                            if q["X_gcm2"] == r["X_gcm2"]])
                     for r in ref], ls=":", color=col, lw=1.0)
    ax.axhline(1.0, color="gray", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Slant opacity [m w.e.]")
    ax.set_title(title, fontsize=9)
    ax.grid(alpha=0.3, which="both")
    leg = ax.legend(title="zenith", fontsize=7, title_fontsize=7, loc="lower left")
    ax.add_artist(leg)


def main() -> int:
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 4.0), sharey=True)
    panel(a, load(RES / "backward_validation_20261019.csv"), "T_backward", True,
          "(a) mode 1, at the calibration nodes")
    panel(b, load(RES / "backward_validation_offnode_20261020.csv"), "T_backward_calibrated",
          False, "(b) mode 2 (default), between the nodes")
    a.set_ylabel(r"$T_{\rm backward}\,/\,T_{\rm forward}$")
    a.set_ylim(0.80, 1.08)
    handles = [Line2D([], [], marker=m, ls="", color="k", mfc="k" if e == "UCMuon-MC" else "none")
               for e, m in MARK.items()]
    a.legend(handles + [Line2D([], [], ls=":", color="k")], list(MARK) + ["CSDA only (vs. mean)"],
             fontsize=7, loc="lower center")
    b.legend(handles, list(MARK), fontsize=7, loc="lower center")
    fig.tight_layout()
    fig.savefig(OUT)
    print(f"[ok] wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
