#!/usr/bin/env python3
"""Multiple-scattering widths, tails and distributions of every engine in the
monoenergetic benchmark (six vertical beams, Standard Rock), per beam and depth.

For each code, beam and depth plane:
  theta0      core width of the projected angle, from the median polar angle
              (Rayleigh: median = theta0 sqrt(2 ln 2)), robust to the tails
  E_s         theta0 / sqrt(sum dx/X0/(beta c p)^2) along the CSDA momentum
              profile: the constant of theta0 = E_s/(beta c p) sqrt(x/X0)
              (13.6 MeV for Highland with its log capped at 1 X0, 14.7 MeV for
              the thick-target limit, growing with x for the uncapped log)
  p99, p999   99th and 99.9th percentiles of the polar angle over theta0
              (3.03 and 3.72 for a Gaussian core)
  sigma_r     RMS radial core of the lateral displacement, sqrt(2) times the
              per-projection width from its median (Rayleigh), and r99 / sigma_r
Writes benchmark/results/mcs_widths_<date>.csv and the binned distributions of
theta/theta0,model and r/sigma_r,model for the cases the paper shows to
benchmark/results/mcs_histograms_<date>.csv, theta0,model and sigma_r,model
being the thick-target model (gui/mcs_margin.py for sigma_r).

    python3 benchmark/analysis/mcs_distributions.py [GEANT4_RUN_PREFIX]
"""
from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "gui"))
sys.path.insert(0, str(REPO / "benchmark" / "analysis"))
import benchmark_analysis as ba                         # noqa: E402
import fast_flux_estimator as ff                        # noqa: E402
import mcs_margin                                       # noqa: E402
from ucmuon_stochastic_driver import ES_THICK_MEV, T_SAT_X0   # noqa: E402

V5 = REPO / "benchmark" / "geant4_muon_rock_v5"
G4_DEFAULT = Path.home() / "simulations/geant4/geant4_muon_rock_v5/outputs/run_stdrock_xy1km_20261022_mono"
SOURCE = V5 / "sources" / "benchmark_surface.dat"
OUT = REPO / "benchmark" / "results"
DATE = "20261022"
M_MU, X0, RHO = 0.1056584, 26.54, 2.65
CASES = [(300, 200.0), (300, 25.0), (50, 50.0), (20, 10.0)]          # (beam GeV, depth m)
CODES = ("Geant4", "MUSIC", "PROPOSAL", "UCMuon-MC", "Bethe-Bloch")


def load(g4_prefix):
    with contextlib.redirect_stdout(io.StringIO()):
        res = {"Geant4": ba.load_geant4(str(g4_prefix)), "MUSIC": ba.load_music(str(V5 / "MUSIC")),
               "PROPOSAL": ba.load_proposal(str(V5 / "PROPOSAL")),
               "Bethe-Bloch": ba.load_bb(str(V5 / "BB")), "UCMuon-MC": ba.load_ucmuon(str(V5 / "UCMuon"))}
    src = np.loadtxt(SOURCE, comments="#", usecols=(0, 1, 2))
    out = {}
    for code, r in res.items():
        m = (r["muons"] if isinstance(r, dict) else r[0]).copy()
        if code == "Geant4":           # its EventID k is source line k + 1; vertical beams
            x0, y0 = src[m["EventID"].to_numpy().astype(int), 1], src[m["EventID"].to_numpy().astype(int), 2]
            m["LatDisp_cm"] = np.hypot(m["ExitX_cm"] - x0, m["ExitY_cm"] - y0)
        m["Beam"] = np.round(m["InitKE_GeV"] + M_MU).astype(int)
        out[code] = m[["Beam", "DepthCm", "AngleScat_Deg", "LatDisp_cm"]]
    return out


def csda_profile(e_tot, depth_m, n=4000):
    """sum dx/X0/(beta p)^2 [GeV^-2] along a vertical path, CSDA momentum."""
    r0 = ff._R_of_T(e_tot - M_MU)
    x = np.linspace(0.0, depth_m * 100.0 * RHO, n + 1)
    xm = 0.5 * (x[1:] + x[:-1])
    e = ff._T_of_R(np.maximum(r0 - xm, 1e-3)) + M_MU
    p = np.sqrt(e * e - M_MU * M_MU)
    return float(np.sum(np.diff(x) / X0 / (p * p / e) ** 2))


def main() -> int:
    g4 = Path(sys.argv[1]) if len(sys.argv) > 1 else G4_DEFAULT
    data = load(g4)
    rows, hists = [], []
    for code in CODES:
        for (beam, dep), g in data[code].groupby(["Beam", "DepthCm"]):
            if len(g) < 2000:
                continue
            depth_m = dep / 100.0
            th = np.radians(g["AngleScat_Deg"].to_numpy())
            r = g["LatDisp_cm"].to_numpy()
            t0 = np.median(th) / np.sqrt(2 * np.log(2))
            s0 = np.median(r) / np.sqrt(2 * np.log(2))                     # per projection
            S = csda_profile(beam, depth_m)
            x_X0 = depth_m * 100 * RHO / X0
            sr_model = mcs_margin.sigma_r_cm(beam, depth_m * 100.0, rho=RHO, x0_gcm2=X0)
            rows.append(dict(code=code, beam_GeV=beam, depth_m=depth_m, n=len(g), x_X0=round(x_X0, 1),
                             theta0_mrad=round(1e3 * t0, 4), E_s_MeV=round(1e3 * t0 / np.sqrt(S), 3),
                             p99=round(np.percentile(th, 99) / t0, 3), p999=round(np.percentile(th, 99.9) / t0, 3),
                             sigma_r_cm=round(np.sqrt(2) * s0, 3),
                             sigma_r_model_cm=round(sr_model, 3) if sr_model else np.nan,
                             r99=round(np.percentile(r, 99) / (np.sqrt(2) * s0), 3)))
            if (beam, depth_m) in CASES:
                t_model = ES_THICK_MEV * 1e-3 * np.sqrt(S)
                for q, v, scale in (("theta", th, t_model), ("r", r, sr_model)):
                    edges = np.linspace(0.0, 8.0, 81)
                    cnt, _ = np.histogram(v / scale, bins=edges)
                    for lo, hi, c in zip(edges[:-1], edges[1:], cnt):
                        hists.append(dict(code=code, beam_GeV=beam, depth_m=depth_m, quantity=q,
                                          lo=lo, hi=hi, count=int(c), n=len(g), overflow=int(np.sum(v / scale >= 8.0))))
    w = pd.DataFrame(rows)
    hdr = (f"# Multiple-scattering widths of the monoenergetic benchmark ({DATE}): "
           f"mcs_distributions.py; Geant4 {g4.name}.\n"
           f"# E_s: 13.6 MeV = Highland with its log capped at 1 X0 (UCMuon-MC and Bethe-Bloch up to 1.4.1); "
           f"{ES_THICK_MEV} MeV = thick-target limit (log stopping at {T_SAT_X0:.1f} X0).\n")
    for name, df in ((f"mcs_widths_{DATE}.csv", w), (f"mcs_histograms_{DATE}.csv", pd.DataFrame(hists))):
        with open(OUT / name, "w") as fh:
            fh.write(hdr)
            df.to_csv(fh, index=False, lineterminator="\n")
        print(f"[ok] {OUT / name}")
    print(w.pivot_table(index=["beam_GeV", "depth_m"], columns="code", values="E_s_MeV").round(2).to_string())
    print(w[w.beam_GeV.isin([300, 50])].pivot_table(index=["beam_GeV", "depth_m"], columns="code",
                                                       values="sigma_r_cm").round(2).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
