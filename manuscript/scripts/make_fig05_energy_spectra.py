#!/usr/bin/env python3
"""
make_fig05_energy_spectra.py — Generate fig05_energy_spectra.pdf:
normalised muon energy spectra at the surface and after 100 m of
Standard Rock, for CosmoALEPH generation and MUSIC transport
(N = 1e6 muons, 10-2500 GeV, cos^2-theta angular mode,
theta_max = 85 deg).  The 10 GeV lower generation bound keeps a usable
number of survivors at 100 m (CSDA threshold 62 GeV): with the CosmoALEPH
p^-3.2 spectrum from 1 GeV only ~1e-4 of the muons would survive.

Runs bin/ucmuon_gen_omp and bin/ucmuon_transport_music_omp in a scratch
directory (manuscript/scripts/_fig05_work) and plots the resulting
13-column surface file and 18-column underground file (Appendix A
formats).  Re-running reuses existing output unless --rerun is given.

Usage:  python3 manuscript/scripts/make_fig05_energy_spectra.py [--rerun]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ucmuon_run_helpers as rh

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF = SCRIPT_DIR.parent / "figs" / "fig05_energy_spectra.pdf"
WORKDIR = SCRIPT_DIR / "_fig05_work"

N_MUONS = 1_000_000
E_MIN, E_MAX = 10.0, 2500.0
DEPTH_M = 100.0
T_MIN_CSDA = 62.0   # GeV — CSDA threshold for 100 m Standard Rock (Sec. 7.1)


def run_simulation() -> tuple[Path, Path]:
    surf = WORKDIR / "muons_surface.dat"
    underg = WORKDIR / "muons_underground.dat"
    if "--rerun" in sys.argv or not (surf.exists() and underg.exists()):
        rh.check_binaries()
        rh.prepare_workdir(WORKDIR)
        print(f"  [RUN] generator: CosmoALEPH, N={N_MUONS:,}")
        dt = rh.run_binary(rh.GEN_EXE,
                           rh.gen_stdin(N_MUONS, e_min=E_MIN, e_max=E_MAX,
                                        outfile=surf.name), WORKDIR)
        print(f"        done in {dt:.1f} s")
        print(f"  [RUN] MUSIC transport: {DEPTH_M:.0f} m Standard Rock")
        dt = rh.run_binary(rh.MUSIC_EXE,
                           rh.music_stdin(surf.name, underg.name,
                                          depth_m=DEPTH_M), WORKDIR)
        print(f"        done in {dt:.1f} s")
    else:
        print(f"  [SKIP] reusing {WORKDIR.name}/ output (--rerun to redo)")
    return surf, underg


def main():
    surf_file, underg_file = run_simulation()

    surf = np.loadtxt(surf_file)
    underg = np.loadtxt(underg_file)
    E_surf = surf[:, 10]                       # col 11: total energy [GeV]
    alive = underg[:, 8] == 1                  # col  9: alive flag
    E_und = underg[alive, 12]                  # col 13: underground energy
    surv = alive.mean()
    print(f"  surface muons: {len(E_surf):,}   "
          f"survivors at {DEPTH_M:.0f} m: {alive.sum():,} "
          f"({100 * surv:.2f}%)")

    bins = np.logspace(np.log10(E_MIN), np.log10(E_MAX), 50)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.hist(E_surf, bins=bins, histtype="step", density=True,
            color="#1565C0", lw=1.8, label="Surface (CosmoALEPH)")
    ax.hist(E_und, bins=bins, histtype="step", density=True,
            color="#FF8F00", lw=1.8,
            label=f"Underground ({DEPTH_M:.0f} m, MUSIC)")
    ax.axvline(T_MIN_CSDA, ls="--", color="gray", lw=1.2,
               label=rf"$T_{{\min}}$ CSDA = {T_MIN_CSDA:.0f} GeV")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(E_MIN, E_MAX)
    ax.set_xlabel("Total energy $E$ [GeV]")
    ax.set_ylabel("Normalised counts [GeV$^{-1}$]")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8.5, loc="upper right")

    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    print(f"[OK] wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
