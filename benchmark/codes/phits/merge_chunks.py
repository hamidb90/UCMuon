#!/usr/bin/env python3
"""Merge the phits_summary.csv of independent PHITS runs on parts of one source
(for example the 32 chunks of run_phits_chunks.sh) into the summary of the
whole source.

    python3 merge_chunks.py OUT.csv CHUNK_DIR [CHUNK_DIR ...]

Each CHUNK_DIR holds a phits_summary.csv written by read_phits_output.py with
--n set to that chunk's number of source muons, which is read from its
n_source.txt. Transmitted muons add up; the mean and standard deviation of the
exit energy (and of the zenith angle) are combined with the transmitted
counts as weights, which is exact for the energy histogram's mean.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, dirs = Path(sys.argv[1]), [Path(d) for d in sys.argv[2:]]
    n_src = 0
    acc: dict[int, dict[str, float]] = {}
    for d in dirs:
        n_src += int((d / "n_source.txt").read_text().split()[0])
        for r in csv.DictReader(open(d / "phits_summary.csv")):
            depth = int(float(r["DepthCm"]))
            a = acc.setdefault(depth, dict(mwe=float(r["MWE"]), n=0.0, e1=0.0, e2=0.0, a1=0.0, a2=0.0))
            n = float(r["N_transmitted"])
            if n == 0:          # read_phits_output.py writes nan for the moments
                continue
            for key, mean, std in (("e", "MeanExitKE_GeV", "StdExitKE_GeV"),
                                   ("a", "MeanZenithAngle_deg", "StdZenithAngle_deg")):
                m, s = float(r[mean] or 0.0), float(r[std] or 0.0)
                a[key + "1"] += n * m
                a[key + "2"] += n * (s * s + m * m)
            a["n"] += n
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["DepthCm", "MWE", "N_transmitted", "Transmission_%", "MeanExitKE_GeV",
                    "StdExitKE_GeV", "MeanZenithAngle_deg", "StdZenithAngle_deg"])
        for depth in sorted(acc):
            a = acc[depth]
            n = a["n"]
            if n == 0:
                w.writerow([depth, f"{a['mwe']:.3f}", "0.0", "0.0000", "nan", "nan", "nan", "nan"])
                continue
            me, ma = a["e1"] / n, a["a1"] / n
            w.writerow([depth, f"{a['mwe']:.3f}", f"{n:.1f}", f"{100 * n / n_src:.4f}", f"{me:.4f}",
                        f"{math.sqrt(max(a['e2'] / n - me * me, 0)):.4f}", f"{ma:.4f}",
                        f"{math.sqrt(max(a['a2'] / n - ma * ma, 0)):.4f}"])
    print(f"[ok] {out}: {len(dirs)} chunks, {n_src} source muons")
    return 0


if __name__ == "__main__":
    sys.exit(main())
