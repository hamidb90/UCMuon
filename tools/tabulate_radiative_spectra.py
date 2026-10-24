#!/usr/bin/env python3
"""Tabulate the muon radiative-loss spectra UCMuon-MC samples its catastrophic
losses from: bremsstrahlung and e+e- pair production (Kelner, Kokoulin and
Petrukhin) and photonuclear interactions (Abramowicz et al. ALLM97), the
cross sections PROPOSAL uses by default for muons, evaluated with PROPOSAL.

Needs PROPOSAL (7.6.2 tested), so it is run once by the maintainers, not by
users:

    ~/venvs/ucmuon/bin/python3 tools/tabulate_radiative_spectra.py

and writes gui/ucmuon_radiative_spectra.npz, which UCMuon-MC loads with numpy
alone. For each material (Standard Rock, water, iron), process and total
energy E on a log grid, on a grid of ln v (v = energy transferred / E):

    ECUM[v] = int_0^v v' dsigma/dv' dv' / int_0^1 v' dsigma/dv' dv'
    NCUM[v] = int_{1e-3}^v dsigma/dv' dv' / int_0^1 v' dsigma/dv' dv'

(integrated on a ln v grid from 1e-7, stored from v = 1e-3, single precision)

so that, for a cut v_cut, the share of the process's mean loss carried by
transfers above v_cut is 1 - ECUM(v_cut) and the number of such transfers per
unit of (mean loss / E) is NCUM(1) - NCUM(v_cut). UCMuon-MC keeps the mean
losses of the PDG 2024 tables and takes only these shapes from here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "gui" / "ucmuon_radiative_spectra.npz"
E_GEV = np.logspace(0.0, 5.0, 41)                   # total energy, 1 GeV to 100 TeV
LNV = np.linspace(np.log(1e-7), 0.0, 1601)          # ln v grid of the integration
LNV_OUT = np.linspace(np.log(1e-3), 0.0, 301)       # ln v grid stored (v_cut >= 1e-3)


def main() -> int:
    import proposal as pp
    pdef = pp.particle.MuMinusDef()
    media = {"rock": pp.medium.StandardRock(), "water": pp.medium.Water(), "iron": pp.medium.Iron()}
    params = {
        "brems": pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(lpm=False),
        "pair": pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(lpm=False),
        "nucl": pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(
            pp.parametrization.photonuclear.ShadowButkevichMikheyev()),
    }
    v = np.exp(LNV)
    dlnv = np.diff(LNV)
    out = {"E_GeV": E_GEV, "lnv": LNV_OUT.astype(np.float32),
           "provenance": np.array(f"PROPOSAL {pp.__version__}: brems and epair KelnerKokoulinPetrukhin "
                                  f"(no LPM), photonuclear AbramowiczLevinLevyMaor97 with "
                                  f"ShadowButkevichMikheyev; MuMinus; media StandardRock, Water, Iron")}
    for mname, med in media.items():
        for proc, par in params.items():
            ecum = np.zeros((len(E_GEV), len(v)))
            ncum = np.zeros((len(E_GEV), len(v)))
            for i, e in enumerate(E_GEV * 1000.0):
                f = np.zeros(len(v))
                for comp, w in zip(med.components, _mass_weights(med)):
                    lim = par.kinematic_limits(pdef, comp, e)
                    ok = (v > lim.v_min) & (v < lim.v_max)
                    f[ok] += w * np.array([par.differential_crosssection(pdef, comp, e, float(x))
                                           for x in v[ok]])
                f = np.maximum(f, 0.0)
                g_e = f * v * v                                   # v dsigma/dv, per d ln v
                g_n = f * v                                       # dsigma/dv, per d ln v
                ce = np.concatenate(([0.0], np.cumsum(0.5 * (g_e[1:] + g_e[:-1]) * dlnv)))
                cn = np.concatenate(([0.0], np.cumsum(0.5 * (g_n[1:] + g_n[:-1]) * dlnv)))
                norm = ce[-1] if ce[-1] > 0 else 1.0
                ecum[i], ncum[i] = ce / norm, cn / norm
            e_out = np.array([np.interp(LNV_OUT, LNV, row) for row in ecum])
            n_out = np.array([np.interp(LNV_OUT, LNV, row) for row in ncum])
            out[f"{mname}_{proc}_ecum"] = e_out.astype(np.float32)
            out[f"{mname}_{proc}_ncum"] = (n_out - n_out[:, :1]).astype(np.float32)   # from v = 1e-3
            j = int(np.searchsorted(LNV_OUT, np.log(0.05)))
            print(f"{mname:5s} {proc:5s}: share of the mean loss above v = 0.05 at "
                  + ", ".join(f"{E_GEV[k]:g} GeV {100 * (1 - e_out[k, j]):.1f} %" for k in (8, 16, 24, 32)),
                  flush=True)
    np.savez_compressed(OUT, **out)
    print(f"[ok] {OUT} ({OUT.stat().st_size / 1024:.0f} kB)")
    return 0


def _mass_weights(med):
    """Weights per component making a per-gram cross section: atoms per molecule x
    1/(molecular mass), relative; PROPOSAL's dsigma/dv is per atom."""
    comps = med.components
    w = np.array([c.atoms_in_molecule for c in comps], float)
    return w / w.sum()


if __name__ == "__main__":
    sys.exit(main())
