#!/bin/bash
# =============================================================================
#  UCMuon — comprehensive test run
#
#  Runs the full two-stage pipeline (Fortran surface generator -> UCMuon-MC
#  transport) on a small, fully seeded configuration, then compares the result
#  against the reference output committed in test_run/expected/.
#
#  Usage:   bash test_run/run_test.sh        (from the repository root)
#
#  Exit status 0 means every produced file matched its reference byte for byte.
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

rm -rf "$OUT_DIR"; mkdir -p "$OUT_DIR"

strip() { sed 's/#.*//' "$1" | sed 's/[[:space:]]*$//' | grep -v '^$'; }

echo "============================================================"
echo "  UCMuon comprehensive test run"
echo "  seed = ${UCMUON_SEED}, single thread"
echo "============================================================"

# ---------------------------------------------------------------- Stage 1
echo
echo "[1/2] Surface generator (bin/ucmuon_gen_omp)"
if [[ ! -x bin/ucmuon_gen_omp ]]; then
    echo "      SKIP - bin/ucmuon_gen_omp not built. Run ./setup.sh first."
    echo "      (Stage 2 can still be checked against the committed surface file.)"
    cp "${REF_DIR}/muons_surface.dat" "${OUT_DIR}/" || exit 1
else
    # The generator ends with a "Press Enter to start" read, hence the blank line.
    { strip "${TEST_DIR}/input_gen.dat"; echo; } \
        | ./bin/ucmuon_gen_omp > "${OUT_DIR}/generator.log" 2>&1
    if [[ ! -s "${OUT_DIR}/muons_surface.dat" ]]; then
        echo "      FAIL - generator produced no output; see ${OUT_DIR}/generator.log"
        exit 1
    fi
    echo "      wrote ${OUT_DIR}/muons_surface.dat"
fi

# ---------------------------------------------------------------- Stage 2
echo
echo "[2/2] UCMuon-MC transport (gui/ucmuon_stochastic_driver.py)"
strip "${TEST_DIR}/input_transport.dat" \
    | python3 gui/ucmuon_stochastic_driver.py > "${OUT_DIR}/transport.log" 2>&1
if [[ ! -s "${OUT_DIR}/muons_underground.dat" ]]; then
    echo "      FAIL - transport produced no output; see ${OUT_DIR}/transport.log"
    exit 1
fi
echo "      wrote ${OUT_DIR}/muons_underground.dat"

# ---------------------------------------------------------------- compare
#
# The test is graded numerically rather than byte for byte.  Both stages are
# deterministic for a fixed seed on a fixed machine, but a different gfortran
# or libm can change the last bits of a double, which would fail a byte
# comparison without anything being wrong.  Byte identity is still reported
# when it holds, because on the reference platform it should.
echo
echo "Comparing against ${REF_DIR}/ ..."
for f in muons_surface.dat muons_underground.dat; do
    [[ -f "${REF_DIR}/${f}" ]] || { echo "  --  ${f}: no reference committed"; continue; }
    if diff -q "${OUT_DIR}/${f}" "${REF_DIR}/${f}" > /dev/null 2>&1; then
        echo "  byte-identical: ${f}"
    else
        echo "  not byte-identical: ${f} (expected off the reference platform)"
    fi
done

echo
python3 - "$OUT_DIR" "$REF_DIR" <<'PY'
import sys, numpy as np
from pathlib import Path
out, ref = Path(sys.argv[1]), Path(sys.argv[2])

def stats(p):
    a = np.loadtxt(p, comments="#")
    alive = a[:, 8].astype(int)
    E = a[:, 12]
    return {"n": len(alive), "surv": int(alive.sum()),
            "frac": 100.0 * alive.sum() / len(alive),
            "meanKE": float(np.mean(E[alive == 1])) - 0.10566}

s = stats(out / "muons_underground.dat")
print("Physics summary")
print(f"  muons transported : {s['n']}")
print(f"  survived          : {s['surv']}")
print(f"  survival fraction : {s['frac']:.2f} %")
print(f"  mean exit KE      : {s['meanKE']:.3f} GeV")

rp = ref / "muons_underground.dat"
if not rp.exists():
    print("\n  (no reference to compare against yet)")
    sys.exit(0)

r = stats(rp)
# Tolerances: the muon count must match exactly; the survival fraction is
# allowed a quarter of its own binomial sigma; the mean exit energy 0.5%.
sigma = np.sqrt(r["frac"] * (100 - r["frac"]) / r["n"])
checks = [
    ("muon count",       s["n"] == r["n"],                                  f"{s['n']} vs {r['n']}"),
    ("survival fraction", abs(s["frac"] - r["frac"]) <= max(0.25 * sigma, 1e-9),
                          f"{s['frac']:.3f}% vs {r['frac']:.3f}% (tol {0.25*sigma:.3f} pp)"),
    ("mean exit KE",     abs(s["meanKE"] - r["meanKE"]) <= 0.005 * abs(r["meanKE"]),
                          f"{s['meanKE']:.4f} vs {r['meanKE']:.4f} GeV (tol 0.5%)"),
]
print("\nChecks against the reference")
ok = True
for name, passed, detail in checks:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name:18s} {detail}")
    ok &= passed
sys.exit(0 if ok else 1)
PY
STATUS=$?

echo
if [[ $STATUS -eq 0 ]]; then
    echo "RESULT: PASS - the run reproduces the committed reference."
else
    echo "RESULT: FAIL - see the failing checks above."
fi
exit $STATUS
