#!/usr/bin/env python3
"""Does PHITS's larger energy loss account for its realistic-source deficit?

On the monoenergetic beams PHITS passes as many muons as Geant4 but loses more
energy: at 200 m, where only the 300 GeV beam survives, its mean exit energy
is below Geant4's. This script turns that difference into a fraction delta of
Geant4's mean energy loss, then makes every muon of Geant4's realistic-source
run lose (1 + delta) times its own loss, and reports how many survivors that
removes and how it shifts their mean energy, to compare with PHITS's own
realistic run.

    python3 benchmark/analysis/phits_energy_loss_check.py \\
        G4_MONO_PREFIX G4_REALISTIC_PREFIX

The prefixes are those of the Geant4 runs' _muons.csv, e.g.
.../outputs/run_stdrock_xy1km_20261022_mono and
.../outputs/run_stdrock_20261019_realistic. PHITS's numbers are read from
benchmark/codes/phits/phits_summary.csv (monoenergetic) and
benchmark/run_realistic_20260614/benchmark_sixcode_summary.csv (realistic).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    mono = pd.read_csv(sys.argv[1] + "_muons.csv", usecols=["DepthCm", "InitKEGeV", "ExitKEGeV"])
    real = pd.read_csv(sys.argv[2] + "_muons.csv", usecols=["DepthCm", "InitKEGeV", "ExitKEGeV"])
    ph_mono = {round(float(r["DepthCm"])): float(r["MeanExitKE_GeV"])
               for r in csv.DictReader(open(REPO / "benchmark" / "codes" / "phits" / "phits_summary.csv"))}
    six = list(csv.DictReader(open(REPO / "benchmark" / "run_realistic_20260614" / "benchmark_sixcode_summary.csv")))
    m = mono[np.isclose(mono.DepthCm, 20000)]
    loss = (m.InitKEGeV - m.ExitKEGeV).mean()
    delta = (m.ExitKEGeV.mean() - ph_mono[20000]) / loss
    print(f"monoenergetic, 200 m: Geant4 exit {m.ExitKEGeV.mean():.2f} GeV after losing {loss:.2f} GeV; "
          f"PHITS exit {ph_mono[20000]:.2f} GeV -> it loses {100 * delta:.1f} % more")
    for depth in (100, 200):
        s = real[np.isclose(real.DepthCm, 100 * depth)]
        lost = s.InitKEGeV - s.ExitKEGeV
        ok = s.ExitKEGeV > delta * lost
        g4 = next(r for r in six if r["code"] == "Geant4" and r["depth_m"] == str(depth))
        ph = next(r for r in six if r["code"] == "PHITS" and r["depth_m"] == str(depth))
        print(f"realistic, {depth} m: Geant4 with {100 * delta:.1f} % more loss: survivors "
              f"{100 * (ok.sum() / len(s) - 1):+.1f} %, their mean exit energy "
              f"{100 * ((s.ExitKEGeV - delta * lost)[ok].mean() / s.ExitKEGeV.mean() - 1):+.1f} %;  PHITS: "
              f"{100 * (float(ph['survival_pct']) / float(g4['survival_pct']) - 1):+.1f} % and "
              f"{100 * (float(ph['mean_exitKE_GeV']) / float(g4['mean_exitKE_GeV']) - 1):+.1f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main())
