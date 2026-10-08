#!/usr/bin/env python3
"""Rebuild benchmark/results/benchmark_summary_v2.csv from the v3 engine outputs.

The engine rows (MUSIC, PROPOSAL, BB, UCMuon) are recomputed from the per-event
files of the v3 benchmark run (benchmark/geant4_muon_rock_v5/<ENGINE>/, July
2026, run_mono_bench.sh) with the loaders and compute_stats of
benchmark_analysis.py. The Geant4 and PHITS rows are kept from the current file:
their raw monoenergetic outputs are no longer on disk, and they did not change.

Up to v1.3.2 the published file was a June 2026 copy with only the UCMuon rows
patched by hand: its BB rows predated the v3 fixes (mean exit energy 135.18
instead of 128.50 GeV at 200 m), and the UCMuon rows were not the v3 run the
paper's tables use.

Usage:  python3 benchmark/analysis/regen_summary_v2.py [output.csv]
        (default: overwrite benchmark/results/benchmark_summary_v2.csv)
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import benchmark_analysis as BA  # noqa: E402

V5 = os.path.join(BENCH, "geant4_muon_rock_v5")
PUBLISHED = os.path.join(BENCH, "results", "benchmark_summary_v2.csv")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else PUBLISHED
    pub = pd.read_csv(PUBLISHED)
    frames = [pub[pub["Code"].isin(["Geant4", "PHITS"])]]
    for loader, sub in [(BA.load_music, "MUSIC"), (BA.load_proposal, "PROPOSAL"),
                        (BA.load_bb, "BB"), (BA.load_ucmuon, "UCMuon")]:
        r = loader(os.path.join(V5, sub))
        s = r["stats"].copy()
        s.insert(0, "Code", r["label"])
        frames.append(s)
    new = pd.concat(frames, ignore_index=True).reindex(columns=pub.columns)
    new.to_csv(out, index=False)
    print(f"wrote {out} ({len(new)} rows)")


if __name__ == "__main__":
    main()
