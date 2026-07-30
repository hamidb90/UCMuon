#!/usr/bin/env python3
"""Measure the validation suite's resolving power as a function of sample size.

This answers the only question that decides how big the reference samples need
to be: *how large an error would the suite have caught?*  Running it produces a
table of the smallest detectable distortion of each kind at each N, which is
both the basis for choosing N and a result worth publishing alongside the port.

Usage:
    python3 measure_power.py                 # default ladder
    python3 measure_power.py 50000 500000    # explicit sizes
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from ucmuref import GenConfig, derived, run
from ucmuref import stats

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "reference" / "resolving_power.json"

# Guan with a cos^2 angular model on a disk: the most-used configuration, and
# the one whose momentum/angle correlation is hardest to reproduce.
BASE = GenConfig(spectrum=4, angular_mode=2, source_mode=1)

DEFAULT_SIZES = (50_000, 200_000, 1_000_000)
CALIBRATION_RUNS = 12


def measure(n: int) -> dict:
    cfg = BASE.replace(nmuons=n)
    probes = stats.default_probes(cfg.angular_mode)

    t0 = time.time()
    thresholds, samples = stats.calibrate(
        lambda s: derived(run(cfg.replace(seed=s))),
        probes, n_runs=CALIBRATION_RUNS, return_samples=True,
    )
    t_cal = time.time() - t0

    # Reuse the calibration samples as disjoint (reference, other) pairs, so
    # measuring the resolving power costs no extra generator runs.
    pairs = [(samples[2 * i], samples[2 * i + 1])
             for i in range(len(samples) // 2)]

    null = stats.compare(pairs[0][0], pairs[0][1], f"null N={n}", thresholds,
                         probes, cfg.angular_mode)
    power, spread = stats.resolving_power(pairs, thresholds, cfg.angular_mode,
                                          return_spread=True)

    return {
        "n": n,
        "n_pairs": len(pairs),
        "null_passes": null.passed,
        "null_failures": [r.probe for r in null.failures],
        "thresholds": thresholds,
        "resolving_power": power,
        "resolving_power_spread": spread,
        "calibration_seconds": round(t_cal, 1),
    }


def main(argv: list[str]) -> int:
    sizes = [int(a) for a in argv[1:]] or list(DEFAULT_SIZES)
    results = []

    for n in sizes:
        print(f"\n=== N = {n:,} ===", flush=True)
        r = measure(n)
        results.append(r)
        status = "PASS" if r["null_passes"] else "FAIL " + ",".join(r["null_failures"])
        print(f"  null hypothesis: {status}   "
              f"(calibration {r['calibration_seconds']}s)")
        for name, mag in r["resolving_power"].items():
            print(f"  {name:<20s} detectable at >= {mag:.4%}"
                  if mag != float("inf") else
                  f"  {name:<20s} NOT detected even at 50%")

    print("\n" + "=" * 66)
    print(f"{'distortion':<20s}" + "".join(f"{n:>14,}" for n in sizes))
    print("-" * 66)
    names = sorted({k for r in results for k in r["resolving_power"]})
    for name in names:
        row = f"{name:<20s}"
        for r in results:
            mag = r["resolving_power"].get(name, float("inf"))
            row += f"{mag:>13.3%} " if mag != float("inf") else f"{'none':>14s}"
        print(row)
    print("=" * 66)
    print("Read as: the smallest error of this kind the suite would have caught.")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2))
    print(f"\nWritten: {OUT_JSON.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
