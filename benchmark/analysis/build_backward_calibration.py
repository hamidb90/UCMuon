#!/usr/bin/env python3
"""Build the forward calibration of the backward CSDA flux integrator (mode 2)
from a calibration run of validate_backward.py.

    python3 benchmark/analysis/build_backward_calibration.py \
        benchmark/results/backward_calibration_run_20261020.csv

C(spectrum, theta, X) = mean over the forward engines of T_forward, divided by
the integrator's mode-1 transmission. Writes
benchmark/results/backward_calibration_20261020.csv (C with its statistical
uncertainty and the spread between engines) and replaces the _CAL table in
gui/ucmuon_backward_mc.py. The run (one row per cell and engine) is the
output of validate_backward.py with --spectra 1,2,3,4 --thetas 0,60,80,85
--x 1e4,2e4,5e4,1e5,2e5,3e5,4e5,5e5 --engines UCMuon-MC,MUSIC.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "benchmark" / "results" / "backward_calibration_20261020.csv"
MODULE = REPO / "gui" / "ucmuon_backward_mc.py"
THETAS = [0.0, 60.0, 80.0, 85.0]
XS = [1.0e4, 2.0e4, 5.0e4, 1.0e5, 2.0e5, 3.0e5, 4.0e5, 5.0e5]


def main() -> int:
    cells = defaultdict(list)
    with open(sys.argv[1]) as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            cells[(int(r["spectrum_mode"]), float(r["theta_deg"]), float(r["X_gcm2"]))].append(r)
    table, rows = {}, []
    for spec in sorted({k[0] for k in cells}):
        grid = np.full((len(THETAS), len(XS)), np.nan)
        for i, th in enumerate(THETAS):
            for j, X in enumerate(XS):
                rs = cells.get((spec, th, X))
                if not rs:
                    continue
                tf = np.array([float(r["T_forward"]) for r in rs])
                ef = np.array([float(r["T_forward_err"]) for r in rs])
                tb = float(rs[0]["T_backward"])
                c = tf.mean() / tb
                grid[i, j] = c
                rows.append(dict(spectrum_mode=spec, theta_deg=th, X_gcm2=X, mwe=X / 100,
                                 C=f"{c:.4f}",
                                 C_stat=f"{np.sqrt((ef ** 2).sum()) / len(tf) / tb:.4f}",
                                 engine_spread=f"{(tf.max() - tf.min()) / tf.mean():.4f}",
                                 engines="+".join(r["engine"] for r in rs)))
        if np.isnan(grid).any():
            sys.exit(f"spectrum {spec}: missing cells in the calibration grid")
        table[spec] = grid
    with open(OUT, "w", newline="") as fh:
        fh.write("# Forward calibration of the backward CSDA flux integrator (mode 2),\n"
                 "# 2026-10-20 (UCMuon-MC re-run 2026-10-23): C = mean T_forward / T_mode1 (benchmark/analysis/\n"
                 "# validate_backward.py with --engines UCMuon-MC,MUSIC, then\n"
                 "# build_backward_calibration.py). Standard Rock, flux above 1 GeV, 1-2500 GeV;\n"
                 "# 2x10^4 muons per cell and engine. C_stat: binomial; engine_spread:\n"
                 "# (max - min) / mean of the engines' T_forward.\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    lines = ["_CAL = {"]
    for spec, grid in table.items():
        lines.append(f"    {spec}: np.array([")
        for row in grid:
            lines.append("        [" + ", ".join(f"{v:.4f}" for v in row) + "],")
        lines.append("    ]),")
    lines.append("}")
    src = MODULE.read_text()
    new, n = re.subn(r"^_CAL = \{.*?^\}$|^_CAL = \{\}$", "\n".join(lines), src, count=1,
                     flags=re.M | re.S)
    if n != 1:
        sys.exit("could not find the _CAL table in the module")
    MODULE.write_text(new)
    print(f"[ok] {OUT}; _CAL updated for spectra {sorted(table)}")
    for spec, grid in table.items():
        print(f"spectrum {spec}: C from {grid.min():.3f} to {grid.max():.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
