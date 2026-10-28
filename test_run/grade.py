#!/usr/bin/env python3
"""Grade a test run against the reference statistics in expected/.

    python3 test_run/grade.py OUTPUT_DIR                 (used by run_test.sh)
    python3 test_run/grade.py OUTPUT_DIR --write-reference REF.csv

Six statistics are computed from the run's files: from Stage 1 (the muons
selected by the detector acceptance cut) the mean log10 of the energy, the
mu+/mu- ratio, the mean cos(theta) and the live time per selected muon; from Stage 2 the survival
fraction and the mean exit kinetic energy. Each comes with its statistical
error, estimated from the run itself.

A check passes when the run lies within TOL standard deviations of the
reference, the deviation combining the run's error with the reference's. The
reference (expected/reference_stats.csv) is a run 20 times larger with another
seed, made by make_reference.sh. With TOL = 4 a correct installation fails a
given check once in about 16 000 runs and one of the six once in about 2600,
while a change of the physics or of the sampling by more than about 4 sigma is
caught whatever the platform, compiler or random stream.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np

M_MU = 0.1056584
TOL = 4.0
NAMES = {
    "mean_log10_E": "mean log10(E/GeV), surface",
    "charge_ratio": "mu+/mu- ratio, surface",
    "mean_cos_theta": "mean cos(theta), surface",
    "live_time_ms": "live time per muon [ms]",
    "survival_pct": "survival fraction [%]",
    "mean_exit_KE_GeV": "mean exit kinetic energy [GeV]",
}


def stats(out: Path) -> dict[str, tuple[float, float]]:
    """{name: (value, statistical error)} from the files in OUT."""
    s = np.loadtxt(out / "muons_selected.dat", comments="#")
    n = len(s)
    loge, cost, q = np.log10(s[:, 10]), np.cos(s[:, 8]), s[:, 11]
    n_p, n_m = int(np.sum(q > 0)), int(np.sum(q < 0))
    r = n_p / n_m
    res = {
        "mean_log10_E": (loge.mean(), loge.std(ddof=1) / math.sqrt(n)),
        "charge_ratio": (r, r * math.sqrt(1 / n_p + 1 / n_m)),
        "mean_cos_theta": (cost.mean(), cost.std(ddof=1) / math.sqrt(n)),
    }
    log = (out / "generator.log").read_text() if (out / "generator.log").exists() else ""
    m = re.search(r"Tried:\s+(\d+)", log)
    t = re.search(r"Live time:\s+([0-9.Ee+-]+)", log)
    if m and t:
        tried, live = int(m.group(1)), 1000.0 * float(t.group(1)) / n
        # per selected muon, so that runs of any size compare; the number tried
        # until n were selected is negative-binomial
        res["live_time_ms"] = (live, live * math.sqrt((1 - n / tried) / n))
    u = np.loadtxt(out / "muons_underground.dat", comments="#")
    alive = u[:, 8] == 1
    p, k = alive.mean(), int(alive.sum())
    ke = u[alive, 12] - M_MU
    res["survival_pct"] = (100 * p, 100 * math.sqrt(p * (1 - p) / len(u)))
    res["mean_exit_KE_GeV"] = (ke.mean(), ke.std(ddof=1) / math.sqrt(k))
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--reference", type=Path, default=Path(__file__).with_name("expected") / "reference_stats.csv")
    ap.add_argument("--write-reference", type=Path, help="write this run's statistics as the reference")
    a = ap.parse_args()
    res = stats(a.out)

    if a.write_reference:
        with open(a.write_reference, "w") as fh:
            fh.write("# Reference statistics of the test run: written by test_run/grade.py "
                     "--write-reference from the run made by test_run/make_reference.sh.\n")
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["statistic", "value", "error"])
            for k, (v, e) in res.items():
                w.writerow([k, f"{v:.8g}", f"{e:.4g}"])
        print(f"[ok] {a.write_reference}")
        return 0

    ref = {r["statistic"]: (float(r["value"]), float(r["error"]))
           for r in csv.DictReader(line for line in open(a.reference) if not line.startswith("#"))}
    print(f"Checks against the reference ({a.reference.parent.name}/{a.reference.name}), "
          f"tolerance {TOL:g} sigma")
    ok = True
    for k, label in NAMES.items():
        if k not in res:
            print(f"  [ -- ] {label:32s} not checked (stage 1 did not run)")
            continue
        (v, e), (rv, re_) = res[k], ref[k]
        sig = math.hypot(e, re_)
        dev = (v - rv) / sig
        passed = abs(dev) <= TOL
        ok &= passed
        print(f"  [{'PASS' if passed else 'FAIL'}] {label:32s} {v:12.5g} +- {e:.2g}   "
              f"reference {rv:.5g} +- {re_:.2g}   ({dev:+.1f} sigma)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
