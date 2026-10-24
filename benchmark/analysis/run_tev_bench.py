#!/usr/bin/env python3
"""TeV muons at km w.e.: UCMuon-MC against MUSIC and PROPOSAL beyond the
monoenergetic benchmark (which stops at 300 GeV and 530 m w.e.).

Vertical beams of 1, 3 and 10 TeV (total energy) through 1, 2 and 3 km w.e.
of Standard Rock (377, 755 and 1132 m at 2.65 g/cm3), 10^4 muons per cell
and engine, the engines' benchmark settings (validate_backward.run_engine).

    python3 benchmark/analysis/run_tev_bench.py [WORKDIR]

Writes benchmark/results/tev_beams_<date>.csv: survival with its binomial
error, mean and RMS exit kinetic energy, per engine, beam and depth.
"""
from __future__ import annotations

import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "benchmark" / "analysis"))
from validate_backward import run_engine, M_MU, RHO     # noqa: E402

DATE = "20261023"
BEAMS_GEV = (1000.0, 3000.0, 10000.0)
X_GCM2 = (1.0e5, 2.0e5, 3.0e5)
ENGINES = ("UCMuon-MC", "MUSIC", "PROPOSAL")
N = int(os.environ.get("NMU", 10_000))


def beam_file(path: Path, e_tot: float) -> None:
    p = np.sqrt(e_tot * e_tot - M_MU * M_MU)
    with open(path, "w") as fh:
        fh.write("# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV  theta_rad  "
                 "phi_rad  E_GeV  charge  hit_flag  det_mask\n")
        for i in range(N):       # spread in x and y: the engines take the source plane from it
            fh.write(f"{i + 1} {i % 100 * 10.0:.1f} {i // 100 % 100 * 10.0:.1f} 0.0 {p:.6f} 0.0 0.0 "
                     f"{-p:.6f} 0.0 0.0 {e_tot:.6f} 1 1 0\n")


def main() -> int:
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "benchmark" / "tev_bench_work"
    work.mkdir(parents=True, exist_ok=True)
    jobs = []
    for e in BEAMS_GEV:
        src = work / f"beam_{e:g}.dat"
        if not src.exists():
            beam_file(src, e)
        for x in X_GCM2:
            for eng in ENGINES:
                jobs.append((eng, src, work / f"{eng}_{e:g}_{x:g}.dat", x / (RHO * 100.0)))
    with ThreadPoolExecutor(max_workers=int(os.environ.get("WORKERS", "8"))) as pool:
        list(pool.map(lambda j: run_engine(*j), jobs))
    rows = []
    for eng, src, out, depth_m in jobs:
        a = np.loadtxt(out, comments="#", usecols=(8, 12), ndmin=2)
        alive = a[:, 0] == 1
        ek = a[alive, 1] - M_MU
        s = alive.mean()
        e = float(np.loadtxt(src, comments="#", usecols=(10,), max_rows=1))
        rows.append(dict(engine=eng, beam_GeV=e, X_gcm2=round(depth_m * RHO * 100.0), mwe=round(depth_m * RHO),
                         n=len(a), survival=round(s, 5), survival_err=round(np.sqrt(s * (1 - s) / len(a)), 5),
                         exitKE_mean_GeV=round(ek.mean(), 2), exitKE_rms_GeV=round(ek.std(), 2)))
        print(f"{eng:10s} {e:6.0f} GeV {depth_m * RHO:5.0f} m w.e.: survival {100 * s:6.2f} %, exit KE "
              f"{ek.mean():8.1f} +- {ek.std():7.1f} GeV", flush=True)
    out = REPO / "benchmark" / "results" / f"tev_beams_{DATE}.csv"
    with open(out, "w") as fh:
        fh.write(f"# TeV beams at km w.e. ({DATE}): benchmark/analysis/run_tev_bench.py, "
                 f"{N} muons per cell, vertical, Standard Rock\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"[ok] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
