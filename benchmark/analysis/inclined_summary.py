#!/usr/bin/env python3
"""Six-code summary of the inclined-beam benchmark (make_inclined_inputs.py):
survival, mean exit kinetic energy and mean deflection from the beam axis per
code, zenith angle, slant length L and beam, with the vertical benchmark as
theta = 0.

    python3 benchmark/analysis/inclined_summary.py [--theta 0,60,80,85]

Reads
  engines  benchmark/geant4_muon_rock_v5/{,inclined/th<theta>/}<ENG>/<ENG>_bench_<L>m.dat
  Geant4   <G4_OUT>/run_stdrock_xy1km_20261022_mono_muons.csv (theta = 0) and
           <G4_OUT>/run_stdrock_xy1km_20261027_incl<theta>_muons.csv
  PHITS    benchmark/codes/phits/phits_summary.csv (theta = 0) and
           benchmark/geant4_muon_rock_v5/inclined/th<theta>/phits/phits_summary.csv
           (beams summed only: PHITS has no per-muon output, and no deflection)
and writes benchmark/results/inclined_beams_<DATE>.csv. A missing input is
skipped with a note.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / "benchmark" / "geant4_muon_rock_v5"
G4_OUT = Path.home() / "simulations" / "geant4" / "geant4_muon_rock_v5" / "outputs"
DATE = "20261027"
M_MU = 0.1056584
L_M = (1, 10, 25, 50, 100, 200)
BEAMS = (5, 10, 20, 50, 100, 300)
N_BEAM = 100_000
ENGINES = {"UCMuon-MC": "UCMuon", "MUSIC": "MUSIC", "PROPOSAL": "PROPOSAL", "Bethe-Bloch": "BB"}


def axis(theta_deg: float) -> np.ndarray:
    t = math.radians(theta_deg)
    return np.array([math.sin(t), 0.0, -math.cos(t)])


def rows_for(code, theta, L, beam_of, alive, ke, defl):
    """One row per beam and one for all beams together (n = source muons)."""
    out = []
    for b in (*BEAMS, "all"):
        sel = np.ones_like(alive) if b == "all" else beam_of == b
        n = N_BEAM * (len(BEAMS) if b == "all" else 1)
        a = alive & sel
        k, p = int(a.sum()), a.sum() / n
        e = ke[a[alive]] if ke is not None else None
        d = defl[a[alive]] if defl is not None else None
        out.append(dict(code=code, theta_deg=theta, L_m=L, beam_GeV=b, n_source=n, n_pass=k,
                        survival_pct=round(100 * p, 4),
                        survival_err_pct=round(100 * math.sqrt(p * (1 - p) / n), 4),
                        exitKE_mean_GeV=round(float(e.mean()), 4) if e is not None and k else "",
                        exitKE_err_GeV=round(float(e.std() / math.sqrt(k)), 4) if e is not None and k > 1 else "",
                        defl_mean_deg=round(float(d.mean()), 5) if d is not None and k else ""))
    return out


def engine_rows(code: str, eng: str, theta: float) -> list[dict]:
    base = BENCH if theta == 0 else BENCH / "inclined" / f"th{theta:g}"
    rows = []
    for L in L_M:
        f = base / eng / f"{eng}_bench_{L}m.dat"
        if not f.exists():
            f = f.with_suffix(".dat.gz")
        if not f.exists():
            print(f"  [skip] {code} theta {theta:g} L {L}: no {f.name}")
            continue
        d = pd.read_csv(f, sep=r"\s+", comment="#", header=None, usecols=[4, 8, 12, 13, 14, 15]).to_numpy()
        beam_of = np.rint(d[:, 0]).astype(int)
        alive = d[:, 1] == 1
        ke = d[alive, 2] - M_MU
        u = d[alive, 3:6]
        cosd = np.clip((u @ axis(theta)) / np.linalg.norm(u, axis=1), -1, 1)
        rows += rows_for(code, theta, L, beam_of, alive, ke, np.degrees(np.arccos(cosd)))
    return rows


def geant4_rows(theta: float) -> list[dict]:
    name = ("run_stdrock_xy1km_20261022_mono" if theta == 0
            else f"run_stdrock_xy1km_{DATE}_incl{theta:g}")
    f = G4_OUT / f"{name}_muons.csv"
    if not f.exists():
        print(f"  [skip] Geant4 theta {theta:g}: no {f}")
        return []
    d = pd.read_csv(f, usecols=["EventID", "DepthCm", "InitKEGeV", "ExitKEGeV",
                                "ExitPxGeVc", "ExitPyGeVc", "ExitPzGeVc"])
    rows = []
    for L in L_M:
        s = d[np.isclose(d["DepthCm"], 100.0 * L)]
        # one entry per crossing muon (the first crossing of each plane)
        beam_x = np.rint(s["InitKEGeV"].to_numpy() + M_MU).astype(int)
        p = s[["ExitPxGeVc", "ExitPyGeVc", "ExitPzGeVc"]].to_numpy()
        cosd = np.clip((p @ axis(theta)) / np.linalg.norm(p, axis=1), -1, 1)
        # survivors only are recorded: build the per-muon arrays over survivors
        alive = np.ones(len(s), bool)
        r = rows_for("Geant4", theta, L, beam_x, alive, s["ExitKEGeV"].to_numpy(),
                     np.degrees(np.arccos(cosd)))
        rows += r
    return rows


def phits_rows(theta: float) -> list[dict]:
    f = (REPO / "benchmark" / "codes" / "phits" / "phits_summary.csv" if theta == 0
         else BENCH / "inclined" / f"th{theta:g}" / "phits" / "phits_summary.csv")
    if not f.exists():
        print(f"  [skip] PHITS theta {theta:g}: no {f}")
        return []
    rows = []
    for r in csv.DictReader(line for line in open(f) if not line.startswith("#")):
        L = round(float(r["DepthCm"]) / 100.0)
        n = N_BEAM * len(BEAMS)
        p = float(r["Transmission_%"]) / 100.0
        rows.append(dict(code="PHITS", theta_deg=theta, L_m=L, beam_GeV="all", n_source=n,
                         n_pass=round(p * n), survival_pct=round(100 * p, 4),
                         survival_err_pct=round(100 * math.sqrt(min(p, 1.0) * (1 - min(p, 1.0)) / n), 4),
                         exitKE_mean_GeV=round(float(r["MeanExitKE_GeV"]), 4), exitKE_err_GeV="",
                         defl_mean_deg=""))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--theta", default="0,60,80,85")
    ap.add_argument("--out", type=Path, default=REPO / "benchmark" / "results" / f"inclined_beams_{DATE}.csv")
    a = ap.parse_args()
    rows = []
    for th in (float(t) for t in a.theta.split(",")):
        rows += geant4_rows(th) + phits_rows(th)
        for code, eng in ENGINES.items():
            rows += engine_rows(code, eng, th)
        print(f"[ok] theta {th:g}")
    with open(a.out, "w") as fh:
        fh.write(f"# Inclined-beam benchmark ({DATE}): the six monoenergetic beams of the vertical\n"
                 "# benchmark (10^5 muons each, Standard Rock) tilted to theta (0 = vertical), scored\n"
                 "# at slant lengths L. survival = muons crossing the plane at L cos(theta) (Geant4,\n"
                 "# PHITS) or completing the slant path L (engines), per source muon; exitKE over the\n"
                 "# survivors; defl = mean angle between exit direction and beam axis. PHITS: beams\n"
                 "# summed only, mean energy from its 6 %-wide bins; it counts every downward\n"
                 "# crossing, so near the horizon a muon that crosses a plane twice counts twice\n"
                 "# (80 deg, 1 m: 100.001 %; dumping every crossing, at most 0.3 % of the muons\n"
                 "# that cross a plane cross it twice). Written by\n"
                 "# benchmark/analysis/inclined_summary.py.\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"[ok] {a.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
