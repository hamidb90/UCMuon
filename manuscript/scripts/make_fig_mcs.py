#!/usr/bin/env python3
"""make_fig_mcs.py: multiple Coulomb scattering of every engine in the
monoenergetic benchmark (benchmark/results/mcs_widths_20261022.csv and
mcs_histograms_20261022.csv, written by benchmark/analysis/mcs_distributions.py).

(a) the effective constant E_s of theta0 = E_s/(beta c p) sqrt(x/X0) against
    depth, 50 and 300 GeV beams, with Highland capped at 1 X0 (up to 1.4.1),
    uncapped, and the thick-target limit (14.7 MeV);
(b, c) polar-angle distributions over the model's theta0, 300 GeV at 25 and
    200 m; (d) lateral displacement over the model's sigma_r, 300 GeV at 200 m.

Usage:  python3 manuscript/scripts/make_fig_mcs.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "benchmark" / "results"
OUT = REPO / "manuscript" / "figs" / "fig_mcs.pdf"
STYLE = {"Geant4": ("k", "o", "-"), "MUSIC": ("#1565C0", "s", "--"), "PROPOSAL": ("#2E7D32", "^", "-."),
         "UCMuon-MC": ("#C62828", "D", "-"), "Bethe-Bloch": ("#FF8F00", "v", ":")}
X0 = 26.54


def main() -> int:
    w = pd.read_csv(RES / "mcs_widths_20261022.csv", comment="#")
    h = pd.read_csv(RES / "mcs_histograms_20261022.csv", comment="#")
    fig, ax = plt.subplots(2, 2, figsize=(9.0, 7.0))
    a = ax[0, 0]
    d = np.logspace(0, np.log10(200), 100)
    a.plot(d, 13.6 * (1 + 0.038 * np.log(d * 265 / X0)), color="gray", lw=1, ls="--",
           label="Highland, log uncapped")
    a.axhline(13.6, color="gray", lw=1, ls=":", label="Highland, log capped at 1 $X_0$ (to 1.4.1)")
    a.axhline(14.7, color="gray", lw=2, alpha=0.5, label="thick-target limit, 14.7 MeV")
    for code, (c, m, _) in STYLE.items():
        for beam, fill in ((300, True), (50, False)):
            s = w[(w.code == code) & (w.beam_GeV == beam)].sort_values("depth_m")
            a.plot(s.depth_m, s.E_s_MeV, marker=m, color=c, mfc=c if fill else "none", ls="-", lw=0.8,
                   ms=5, label=code if beam == 300 else None)
    a.set_xscale("log"); a.set_xlabel("Depth [m of Standard Rock]"); a.set_ylabel(r"$E_s$ [MeV]")
    a.set_ylim(13.2, 17.8); a.set_title("(a) core width, 300 GeV (filled) and 50 GeV (open)", fontsize=9)
    a.legend(fontsize=6.5, ncol=2, loc="upper left"); a.grid(alpha=0.3, which="both")
    for axx, (beam, dep, q, title) in zip((ax[0, 1], ax[1, 0], ax[1, 1]),
                                          ((300, 25.0, "theta", "(b) polar angle, 300 GeV, 25 m"),
                                           (300, 200.0, "theta", "(c) polar angle, 300 GeV, 200 m"),
                                           (300, 200.0, "r", "(d) lateral displacement, 300 GeV, 200 m"))):
        for code, (c, _, ls) in STYLE.items():
            s = h[(h.code == code) & (h.beam_GeV == beam) & (h.depth_m == dep) & (h.quantity == q)]
            if s.empty:
                continue
            x = 0.5 * (s.lo + s.hi).to_numpy(); dens = s["count"].to_numpy() / s.n.iloc[0] / (s.hi - s.lo).to_numpy()
            axx.step(x, np.where(dens > 0, dens, np.nan), where="mid", color=c, ls=ls, lw=1.1, label=code)
        u = np.linspace(0, 8, 400)
        if q == "theta":
            axx.plot(u, u * np.exp(-u * u / 2), color="gray", lw=2, alpha=0.4, label="Gaussian core")
            axx.set_xlabel(r"$\theta\,/\,\theta_0$ (thick-target model)")
        else:
            axx.plot(u, 2 * u * np.exp(-u * u), color="gray", lw=2, alpha=0.4, label="Gaussian core")
            axx.axvline(2.5, color="gray", lw=0.8, ls=":")
            axx.text(2.55, 0.3, r"margin $2.5\,\sigma_r$", fontsize=7, color="gray")
            axx.set_xlabel(r"$r\,/\,\sigma_r$" + " (gui/mcs_margin.py)")
        axx.set_yscale("log"); axx.set_ylim(1e-5, 2); axx.set_xlim(0, 8)
        axx.set_ylabel("probability density"); axx.set_title(title, fontsize=9)
        axx.grid(alpha=0.3); axx.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT)
    print(f"[ok] wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
