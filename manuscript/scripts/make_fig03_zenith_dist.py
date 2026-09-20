#!/usr/bin/env python3
"""
make_fig03_zenith_dist.py — Generate fig03_zenith_dist.pdf:
normalised zenith-angle distributions sampled by the surface generator,
for the standard cos^2(theta) mode and the energy-dependent
self-consistent Guan (2015) CDF (Mode 4) at T = 10 and 100 GeV.

The curves are the exact sampling PDFs of sample_cos2 and
sample_guan_angle in src/generator/ucmuon_source_module.f90
(both include the sin(theta) phase-space factor; the Guan sampler
additionally weights by the projected-area factor cos(theta)).

Usage:  python3 manuscript/scripts/make_fig03_zenith_dist.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ucmuon_flux_models as fm

OUT_PDF = Path(__file__).resolve().parent.parent / "figs" / "fig03_zenith_dist.pdf"

THETA_MAX_DEG = 70.0


def main():
    th_max = np.radians(THETA_MAX_DEG)
    th = np.linspace(1e-4, th_max, 600)
    deg = np.degrees(th)

    N_COSN = 6   # representative exponent for the generalised cosⁿθ mode

    pdf_cos2 = fm.pdf_cos2_theta(th, th_max)
    pdf_cos3 = fm.pdf_cosn_theta(th, th_max, 3)
    pdf_cosn = fm.pdf_cosn_theta(th, th_max, N_COSN)
    # sample_guan_angle is called with the muon TOTAL energy
    pdf_g10 = fm.pdf_guan_theta(th, th_max, 10.0 + fm.MUON_MASS)
    pdf_g100 = fm.pdf_guan_theta(th, th_max, 100.0 + fm.MUON_MASS)

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(deg, pdf_cos2, color="#1565C0", lw=2.0,
            label=r"$\cos^2\theta$ (Modes 1–3, 6)")
    ax.plot(deg, pdf_cos3, color="#7B1FA2", lw=1.7, ls="-",
            label=r"$\cos^3\theta$ (Mode 7, Reyna)")
    ax.plot(deg, pdf_cosn, color="#FF8F00", lw=1.7, ls=(0, (5, 1)),
            label=rf"$\cos^n\theta$ generalised ($n={N_COSN}$)")
    ax.plot(deg, pdf_g10, color="#C62828", lw=1.7, ls="--",
            label=r"Guan (2015) CDF, $T = 10$ GeV (Modes 4–5)")
    ax.plot(deg, pdf_g100, color="#2E7D32", lw=1.7, ls="-.",
            label=r"Guan (2015) CDF, $T = 100$ GeV (Modes 4–5)")

    ax.axvline(THETA_MAX_DEG, color="0.5", ls=":", lw=1.0)
    ax.text(THETA_MAX_DEG - 1.0, ax.get_ylim()[1] * 0.55,
            r"$\theta_{\max} = 70^\circ$", rotation=90,
            ha="right", va="center", fontsize=8, color="0.4")

    ax.set_xlim(0, 75)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"Zenith angle $\theta$ [deg]")
    ax.set_ylabel(r"Normalised PDF $p(\theta)$ [rad$^{-1}$]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper left")

    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    print(f"[OK] wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
