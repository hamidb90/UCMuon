#!/bin/bash
# =============================================================================
#  Regenerate test_run/expected/ (a few minutes).
#
#  1. reference_stats.csv: the statistics run_test.sh grades against, from a
#     run 20 times larger than the test (2x10^5 selected muons) with other
#     seeds (generator 1, transport 7), so that the reference does not share
#     the test's random stream.
#  2. muons_selected.dat.gz and sha256.txt: the test-size run on this
#     platform, the first for machines without a Fortran compiler (Stage 2 then
#     runs on it), the second for the byte-identity report.
#
#  Usage:   bash test_run/make_reference.sh      (from the repository root)
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
REF=test_run/expected
BIG=$(mktemp -d)
trap 'rm -rf "$BIG"' EXIT
strip() { sed 's/#.*//' "$1" | sed 's/[[:space:]]*$//' | grep -v '^$'; }

echo "[1/2] reference run: 2x10^5 selected muons, seeds 1 and 7"
# value 15 of the generator input is nmuons, value 18 of the transport input the seed
{ strip test_run/input_gen.dat | awk 'NR == 15 { $0 = 200000 } { print }' \
      | sed "s|test_run/output/|${BIG}/|"; echo; } \
    | UCMUON_SEED=1 ./bin/ucmuon_gen_omp > "${BIG}/generator.log" 2>&1
strip test_run/input_transport.dat | awk 'NR == 18 { $0 = 7 } { print }' \
    | sed "s|test_run/output/|${BIG}/|" \
    | python3 gui/ucmuon_stochastic_driver.py > "${BIG}/transport.log" 2>&1
python3 test_run/grade.py "$BIG" --write-reference "${REF}/reference_stats.csv"

echo "[2/2] test-size run with the test's seeds"
bash test_run/run_test.sh > /dev/null || true
gzip -9 -n -c test_run/output/muons_selected.dat > "${REF}/muons_selected.dat.gz"
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }
{
    echo "# SHA-256 of the test run's output on the reference platform (make_reference.sh)"
    for f in muons_selected.dat muons_underground.dat; do echo "$(sha "test_run/output/$f") $f"; done
} > "${REF}/sha256.txt"
echo "[ok] ${REF}/: reference_stats.csv, muons_selected.dat.gz, sha256.txt"
