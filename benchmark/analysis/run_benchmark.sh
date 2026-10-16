#!/bin/bash
# Six-code comparison figures and summary (benchmark/figures/v2_six_code/) from
# the raw per-event outputs, which are not in the repository (about 8 GB):
#
#   bash benchmark/analysis/run_benchmark.sh GEANT4_RUN_PREFIX [OUTDIR]
#
# GEANT4_RUN_PREFIX  the reference Geant4 run of the monoenergetic population
#                    (600 000 muons, run_file_20260522_210114; its _muons.csv
#                    is about 400 MB), e.g.
#                    ~/simulations/geant4/geant4_muon_rock_v5/outputs/run_file_20260522_210114
# OUTDIR             default: benchmark/figures/v2_six_code
#
# The UCMuon engines' outputs are read from benchmark/geant4_muon_rock_v5/<ENGINE>/
# (written by run_mono_bench.sh), PHITS from its summary in benchmark/codes/phits/
# (the per-event PHITS output was not kept). The summary it writes equals
# benchmark/results/benchmark_summary_v2.csv.
#
# Up to 1.4.0 this script took the most recently modified Geant4 CSV and the
# phits/ folder of the scratch tree, which by then held the realistic-source
# runs, and it worked only from inside that tree.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BENCH="$(dirname "$HERE")"
V5="$BENCH/geant4_muon_rock_v5"
G4="${1:?usage: run_benchmark.sh GEANT4_RUN_PREFIX [OUTDIR]}"
OUT="${2:-$BENCH/figures/v2_six_code}"

[[ -f "${G4}_muons.csv" ]] || { echo "ERROR: ${G4}_muons.csv not found" >&2; exit 1; }
echo "Using Geant4 run: $G4"

python3 "$HERE/benchmark_analysis.py" \
    --geant4   "$G4" \
    --phits    "$BENCH/codes/phits/" \
    --music    "$V5/MUSIC/" \
    --proposal "$V5/PROPOSAL/" \
    --bb       "$V5/BB/" \
    --ucmuon   "$V5/UCMuon/" \
    --outdir   "$OUT"
