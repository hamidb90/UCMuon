#!/usr/bin/env python3
"""Geant4 survival per beam and depth for the monoenergetic benchmark, the
small CSV manuscript/scripts/make_survival_table.py reads (the per-crossing
Geant4 output, 400 MB, is not shipped).

    python3 benchmark/analysis/geant4_survival_by_beam.py GEANT4_RUN_PREFIX

Writes benchmark/results/geant4_survival_by_beam.csv: muons crossing each
scoring plane / 10^5 per beam.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "benchmark" / "results" / "geant4_survival_by_beam.csv"
N_BEAM = 100_000
M_MU = 0.10566


def main() -> int:
    prefix = Path(sys.argv[1])
    d = np.loadtxt(f"{prefix}_muons.csv", delimiter=",", skiprows=1, usecols=(1, 4))
    beam = np.round(d[:, 1] + M_MU).astype(int)
    depth = np.round(d[:, 0] / 100.0, 3)
    rows = []
    for b in sorted(set(beam)):
        for dep in sorted(set(depth)):
            n = int(np.sum((beam == b) & (depth == dep)))
            rows.append(dict(beam_GeV=int(b), depth_m=f"{dep:g}", survival_pct=f"{100.0 * n / N_BEAM:.3f}"))
    with open(OUT, "w") as fh:
        fh.write("# Geant4 survival per beam and depth, monoenergetic benchmark (10^5 muons per\n"
                 f"# beam), Standard Rock: {prefix.name} (benchmark/codes/geant4).\n"
                 "# Muons crossing each scoring plane / 10^5. Read by\n"
                 "# manuscript/scripts/make_survival_table.py; written by\n"
                 "# benchmark/analysis/geant4_survival_by_beam.py.\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"[ok] {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
