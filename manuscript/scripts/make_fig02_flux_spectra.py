#!/usr/bin/env python3
"""
make_fig02_flux_spectra.py — Generate fig02_flux_spectra.pdf:
differential vertical muon intensity I(T, 0 deg) for the seven muon flux
parametrisations (Modes 1-7) implemented in UCMuon.

Modes 1, 2, 4-7 are evaluated analytically with the exact constants of
src/generator/ucmuon_source_module.f90 (see ucmuon_flux_models.py).
Mode 3 (PARMA/EXPACS) is evaluated by compiling and running the small
Fortran helper parma_vertical_spec.f90 against the PARMA subroutines in
src/parma/ and the data directory data/EXPACS/parma (sea level, zero
solar modulation, geomagnetic reference site).  If gfortran or the PARMA
data directory is unavailable the PARMA curve is skipped with a warning.

Usage:  python3 manuscript/scripts/make_fig02_flux_spectra.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ucmuon_flux_models as fm

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF = SCRIPT_DIR.parent / "figs" / "fig02_flux_spectra.pdf"
PARMA_DATA = REPO_ROOT / "data" / "EXPACS" / "parma"

# PARMA reference site: sea level, mid-latitude European site (UCLouvain),
# zero solar modulation (W = 0).
PARMA_SITE = dict(lat=50.67, lon=4.62, alt_km=0.0,
                  year=2026, month=1, day=1, W=0.0)

PDG_VERT_FLUX = 7.0e-3  # cm^-2 s^-1 sr^-1, integrated > 1 GeV (PDG 2022)


def parma_curve(t_min: float, t_max: float, npts: int = 120):
    """Compile (if needed) and run the PARMA helper; returns (T, I) or None."""
    exe = SCRIPT_DIR / "_parma_spec"
    if not PARMA_DATA.is_dir():
        print(f"  [SKIP] PARMA data not found at {PARMA_DATA}")
        return None
    if not exe.exists():
        gfortran = shutil.which("gfortran")
        if gfortran is None:
            print("  [SKIP] gfortran not found; PARMA curve omitted")
            return None
        cmd = [gfortran, "-O2", "-ffree-line-length-none", "-o", str(exe),
               str(REPO_ROOT / "src/parma/parma_path_module.f90"),
               str(REPO_ROOT / "src/parma/parma_subroutines.f90"),
               str(SCRIPT_DIR / "parma_vertical_spec.f90"),
               "-J", str(SCRIPT_DIR)]
        print("  [BUILD]", " ".join(cmd))
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  [SKIP] PARMA helper build failed:\n{r.stderr}")
            return None
    stdin = "\n".join([
        str(PARMA_DATA),
        str(PARMA_SITE["lat"]), str(PARMA_SITE["lon"]),
        str(PARMA_SITE["alt_km"]),
        str(PARMA_SITE["year"]), str(PARMA_SITE["month"]),
        str(PARMA_SITE["day"]), str(PARMA_SITE["W"]),
        str(t_min), str(t_max), str(npts),
    ]) + "\n"
    r = subprocess.run([str(exe)], input=stdin, capture_output=True,
                       text=True, cwd=REPO_ROOT)
    if r.returncode != 0:
        print(f"  [SKIP] PARMA helper run failed:\n{r.stderr}")
        return None
    data = np.array([[float(v) for v in line.split()]
                     for line in r.stdout.splitlines() if line.strip()])
    print(f"  [OK]  PARMA curve: {len(data)} points")
    return data[:, 0], data[:, 1]


def main():
    T = np.logspace(np.log10(0.2), np.log10(2.0e4), 400)  # kinetic energy GeV
    p = fm.p_of_T(T)
    E = T + fm.MUON_MASS
    jac = fm.jac_dp_dT(T)  # dN/dT = dN/dp * dp/dT

    # Modes 1 and 2 are sampling shapes with arbitrary absolute
    # normalisation in the generator; for display they are pinned to the
    # absolutely normalised Guan (2015) intensity at T = 100 GeV.
    T_ANCHOR = 100.0
    i_guan_anchor = float(fm.guan_flux(T_ANCHOR + fm.MUON_MASS, 1.0))
    cosmo = fm.cosmoaleph_dndp(p) * jac
    cosmo *= i_guan_anchor / float(
        fm.cosmoaleph_dndp(fm.p_of_T(T_ANCHOR)) * fm.jac_dp_dT(T_ANCHOR))
    plaw = fm.powerlaw_dndp(p) * jac
    plaw *= i_guan_anchor / float(
        fm.powerlaw_dndp(fm.p_of_T(T_ANCHOR)) * fm.jac_dp_dT(T_ANCHOR))

    curves = []  # (T, I, label, style)
    curves.append((T, cosmo,
                   "Mode 1 — CosmoALEPH (shape, norm. at 100 GeV)",
                   dict(color="#1565C0", ls="-", lw=2.0)))
    curves.append((T, plaw,
                   r"Mode 2 — Power-law $E^{-3.7}$ (shape, norm. at 100 GeV)",
                   dict(color="#C62828", ls="--", lw=1.6)))

    pc = parma_curve(0.2, 2.0e4)
    if pc is not None:
        curves.append((pc[0], pc[1],
                       "Mode 3 — PARMA/EXPACS (sea level, W=0)",
                       dict(color="#6A1B9A", ls="-", lw=1.6)))

    Tg = T[T >= 10.0]   # Guan/Frosin valid above ~10 GeV
    curves.append((Tg, fm.guan_flux(Tg + fm.MUON_MASS, 1.0,
                                    fm.GUAN_A, fm.GUAN_B),
                   "Mode 4 — Guan et al. (2015)",
                   dict(color="#2E7D32", ls="-", lw=1.6)))
    curves.append((Tg, fm.guan_flux(Tg + fm.MUON_MASS, 1.0,
                                    fm.FROSIN_A, fm.FROSIN_B),
                   "Mode 5 — Frosin et al. (2025)",
                   dict(color="#00838F", ls="--", lw=1.6)))
    Tb = T[T >= 1.0]    # Bugaev / Reyna-Bugaev shown from 1 GeV
    curves.append((Tb, fm.guan_flux(Tb + fm.MUON_MASS, 1.0,
                                    fm.BUGAEV_A, fm.BUGAEV_B),
                   "Mode 6 — Bugaev/Gaisser (1990)",
                   dict(color="#FF8F00", ls="-.", lw=1.6)))
    # Reyna-Bugaev log-polynomial fitted on 1-2000 GeV/c; it diverges
    # unphysically outside that range, so the curve is clipped there.
    Tr = T[(T >= 1.0) & (T <= 2000.0)]
    curves.append((Tr, fm.reyna_flux(fm.p_of_T(Tr), 1.0) * fm.jac_dp_dT(Tr),
                   "Mode 7 — Reyna–Bugaev (2006), 1–2000 GeV",
                   dict(color="#7B1FA2", ls=":", lw=1.9)))

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for x, y, label, style in curves:
        ok = y > 0
        ax.loglog(x[ok], y[ok], label=label, **style)

    ax.axhline(PDG_VERT_FLUX, color="0.4", ls="--", lw=1.0)
    ax.text(2.5e-1, PDG_VERT_FLUX * 1.35,
            r"PDG $\Phi_{\rm vert}(>1\,{\rm GeV}) = 7.0\times10^{-3}$"
            r"$\,{\rm cm^{-2}s^{-1}sr^{-1}}$",
            fontsize=7.5, color="0.35")

    ax.set_xlim(2e-1, 2e4)
    ax.set_ylim(1e-13, 3e-2)
    ax.set_xlabel(r"Kinetic energy $T$ [GeV]")
    ax.set_ylabel(r"$I(T,\,0^\circ)$ "
                  r"[cm$^{-2}$ s$^{-1}$ sr$^{-1}$ GeV$^{-1}$]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7.5, loc="lower left", framealpha=0.9)

    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    print(f"[OK] wrote {OUT_PDF}")

    # Sanity check: PARMA vertical integral above 1 GeV vs PDG
    if pc is not None:
        Tp, Ip = pc
        ok = Tp >= 1.0
        phi = np.trapz(Ip[ok], Tp[ok])
        print(f"  PARMA Phi_vert(>1 GeV) = {phi:.2e} cm^-2 s^-1 sr^-1 "
              f"(PDG: {PDG_VERT_FLUX:.1e})")


if __name__ == "__main__":
    main()
