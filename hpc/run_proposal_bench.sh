#!/bin/bash
# =============================================================================
#  UCMuon PROPOSAL benchmark — single-thread, 6 depths, SLURM
#  Lemaitre4 / CECI cluster
#
#  Usage (from UCMuon/ project root):
#    sbatch hpc/run_proposal_bench.sh
#
#  Outputs land in:
#    benchmark/geant4_muon_rock_v5/PROPOSAL/
#      PROPOSAL_bench_1m.dat   ...   PROPOSAL_bench_200m.dat
#      PROPOSAL_bench_*_timing.txt
#      PROPOSAL_bench_*_stopped.dat
#
#  Requirements on the cluster:
#    pip install proposal    (or module load python + pip install --user proposal)
#    source file at SOURCES path below (edit if needed)
# =============================================================================

#SBATCH --job-name=proposal_bench
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=logs/proposal_bench_%j.out
#SBATCH --error=logs/proposal_bench_%j.err

set -euo pipefail

# ─── paths ───────────────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GUI_DIR="${REPO_ROOT}/gui"
OUTDIR="${REPO_ROOT}/benchmark/geant4_muon_rock_v5/PROPOSAL"
SOURCES="${REPO_ROOT}/benchmark/geant4_muon_rock_v5/sources/benchmark_surface.dat"
PROP_TABLES="${HOME}/.proposal_tables"

# ─── single thread ───────────────────────────────────────────────────────────
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ─── setup ───────────────────────────────────────────────────────────────────
cd "${REPO_ROOT}"
mkdir -p "${OUTDIR}" logs

if [[ ! -f "${SOURCES}" ]]; then
    echo "ERROR: source file not found: ${SOURCES}"; exit 1
fi
if [[ ! -f "${GUI_DIR}/proposal_driver.py" ]]; then
    echo "ERROR: ${GUI_DIR}/proposal_driver.py not found"; exit 1
fi

# Use venv python if present (macOS ABI fix), else system python3
PROP_PY="${HOME}/venvs/ucmuon/bin/python3"
if [[ ! -x "${PROP_PY}" ]]; then PROP_PY="python3"; fi

# Check PROPOSAL imports on this node before wasting wall time
"${PROP_PY}" -c "import proposal as pp; print('PROPOSAL', pp.__version__, 'OK')" || {
    echo "ERROR: PROPOSAL not importable. Run: pip install proposal (or bash setup.sh)"; exit 1
}

echo "============================================================"
echo "  PROPOSAL benchmark (single thread, Moliere MCS)"
echo "============================================================"
echo "  Repo root  : ${REPO_ROOT}"
echo "  Source     : ${SOURCES}"
echo "  Output dir : ${OUTDIR}"
echo "  Tables     : ${PROP_TABLES}"
echo "  Started    : $(date)"
echo "============================================================"
echo ""

DEPTHS=(1 10 25 50 100 200)
t_total=0

for DEPTH in "${DEPTHS[@]}"; do
    OUTFILE="${OUTDIR}/PROPOSAL_bench_${DEPTH}m.dat"
    echo -n "  depth=${DEPTH}m ... "
    t0=$(date +%s)

    "${PROP_PY}" "${GUI_DIR}/proposal_driver.py" << EOF
${SOURCES}
${OUTFILE}
${DEPTH}
1
1
500
0.001
3
${PROP_TABLES}
EOF

    t1=$(date +%s)
    dt=$(( t1 - t0 ))
    t_total=$(( t_total + dt ))
    echo "${dt}s"
done

echo ""
echo "  PROPOSAL total wall time: ${t_total}s"
echo "  All outputs: ${OUTDIR}/"
echo "  Finished   : $(date)"
echo "============================================================"
