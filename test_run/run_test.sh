#!/bin/bash
# =============================================================================
#  UCMuon — comprehensive test run
#
#  Runs the full two-stage pipeline (Fortran surface generator with a detector
#  acceptance cut -> UCMuon-MC transport) on a small, seeded configuration and
#  grades the result statistically against expected/reference_stats.csv.
#
#  Usage:   bash test_run/run_test.sh        (from the repository root)
#
#  Exit status 0 means every statistic (four of the selected surface muons, two
#  of the transported ones) lies within 4 standard deviations of the
#  reference; see test_run/grade.py. Whether the output is also byte-identical
#  to the reference platform's is reported, but not required: a different
#  compiler, libm or random stream gives a statistically equivalent run.
#
#  Requirements: gfortran-built bin/ucmuon_gen_omp (run ./setup.sh) for Stage 1,
#  and Python 3.9+ with NumPy for Stage 2.  Stage 2 needs no compiler.
# =============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

TEST_DIR="test_run"
OUT_DIR="${TEST_DIR}/output"
REF_DIR="${TEST_DIR}/expected"

# Fixed seed and single-threaded execution: both stages are then deterministic
# and independent of how many cores the reviewer's machine has.
export UCMUON_SEED=20260807
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# 1 when stage 1 actually ran; 0 when its output was taken from expected/.
GEN_RAN=1

rm -rf "$OUT_DIR"; mkdir -p "$OUT_DIR"

strip() { sed 's/#.*//' "$1" | sed 's/[[:space:]]*$//' | grep -v '^$'; }

echo "============================================================"
echo "  UCMuon comprehensive test run"
echo "  seed = ${UCMUON_SEED}, single thread"
echo "============================================================"

# ---------------------------------------------------------------- Stage 1
echo
echo "[1/2] Surface generator (bin/ucmuon_gen_omp), detector acceptance cut"
if [[ ! -x bin/ucmuon_gen_omp ]]; then
    echo "      SKIP - bin/ucmuon_gen_omp not built. Run ./setup.sh first."
    echo "      (Stage 2 is still run, on the reference platform's selected muons.)"
    gunzip -c "${REF_DIR}/muons_selected.dat.gz" > "${OUT_DIR}/muons_selected.dat" || exit 1
    GEN_RAN=0
else
    # The generator ends with a "Press Enter to start" read, hence the blank line.
    { strip "${TEST_DIR}/input_gen.dat"; echo; } \
        | ./bin/ucmuon_gen_omp > "${OUT_DIR}/generator.log" 2>&1
    if [[ ! -s "${OUT_DIR}/muons_selected.dat" ]]; then
        echo "      FAIL - generator produced no output; see ${OUT_DIR}/generator.log"
        exit 1
    fi
    echo "      wrote ${OUT_DIR}/muons_selected.dat"
fi

# ---------------------------------------------------------------- Stage 2
echo
echo "[2/2] UCMuon-MC transport (gui/ucmuon_stochastic_driver.py), 25 m of rock"
strip "${TEST_DIR}/input_transport.dat" \
    | python3 gui/ucmuon_stochastic_driver.py > "${OUT_DIR}/transport.log" 2>&1
if [[ ! -s "${OUT_DIR}/muons_underground.dat" ]]; then
    echo "      FAIL - transport produced no output; see ${OUT_DIR}/transport.log"
    exit 1
fi
echo "      wrote ${OUT_DIR}/muons_underground.dat"

# ---------------------------------------------------------------- regression
#
# On the reference platform both stages reproduce the reference run's files
# exactly (expected/sha256.txt); elsewhere they need not, so this is reported
# but not graded.
echo
echo "Regression against the reference platform's output (not graded):"
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }
while read -r want name; do
    [[ -z "$want" || "$want" == \#* ]] && continue
    if [[ "$name" == muons_selected.dat && "$GEN_RAN" -eq 0 ]]; then
        echo "  not checked: ${name} (stage 1 skipped; this is the reference file)"
    elif [[ "$(sha "${OUT_DIR}/${name}")" == "$want" ]]; then
        echo "  byte-identical: ${name}"
    else
        echo "  not byte-identical: ${name} (expected off the reference platform)"
    fi
done < "${REF_DIR}/sha256.txt"

# ---------------------------------------------------------------- grade
echo
python3 "${TEST_DIR}/grade.py" "$OUT_DIR"
STATUS=$?

echo
if [[ $STATUS -eq 0 ]]; then
    echo "RESULT: PASS - every statistic agrees with the reference."
else
    echo "RESULT: FAIL - see the failing checks above."
fi
exit $STATUS
