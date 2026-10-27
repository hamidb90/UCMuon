#!/bin/bash
# =============================================================================
# run_phits_mpi.sh: a PHITS benchmark deck under MPI on a Slurm cluster
# (written for CECI Lemaitre4 and PHITS 3.37).
#
#   sbatch --ntasks=61 --export=ALL,PHITSPATH=/path/to/phits \
#          run_phits_mpi.sh muon_rock.inp benchmark_phits.dat RUN_DIR
#
# PHITS under MPI has one master and N-1 workers; each worker reads its own
# source dump, SOURCE.001 .. SOURCE.<workers>, and maxcas * maxbch must equal
# the number of records exactly (otherwise "dump source file ended"). The
# script therefore splits SOURCE into equal per-worker files, sets
# maxbch = 10 x workers and maxcas = records / maxbch, and points file(1) at
# PHITSPATH, all in a copy of the deck inside RUN_DIR. The records are
# shuffled first (fixed seed): PHITS waits for every worker after each batch,
# and the six-beam source is sorted by energy, so unshuffled each worker got one
# beam and every batch took as long as the 300 GeV worker (the first run of
# 2026-10-22 timed out after 9 of 10 batches, its CPU time inflated by the
# waiting). Records are independent muons, so the shuffle changes nothing else.
#
# TZ is set because PHITS reads its CPU clock millions of times per run (for
# the time-by-process table in phits.out) and glibc's localtime, reached from
# that clock, re-reads /etc/localtime on every call when TZ is unset: on
# Lemaitre4 that made PHITS four times slower (211 s against 53 s for the same
# 50 muons). The remaining clock calls (getrusage, uname) still make the
# cluster build several times slower per muon than PHITS's own single-thread
# Mac binary, so CPU times from this script are not PHITS's speed; the paper
# times PHITS on one core of the M4.
#
# Sizes used for the benchmark (records must divide by 10 x workers):
#   six beams, 600 000 muons:  --ntasks=61  (60 workers, maxcas 1000)
#   realistic, 100 000 muons:  --ntasks=51  (50 workers, maxcas 200)
# =============================================================================
#SBATCH --job-name=phits_bench
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=08:00:00

set -euo pipefail
DECK=$(readlink -f "$1"); SOURCE=$(readlink -f "$2"); RUN_DIR=$3
: "${PHITSPATH:?set PHITSPATH to the PHITS installation (the directory holding bin/ and data/)}"
WORKERS=$(( ${SLURM_NTASKS:?run under sbatch} - 1 ))
PHITS_BIN="${PHITS_BIN:-${PHITSPATH}/bin/phits_LinIFX_OMP_MPI}"
[[ -x "${PHITS_BIN}" ]] || { echo "FATAL: ${PHITS_BIN} not found"; exit 1; }

NREC=$(grep -c . "${SOURCE}")
MAXBCH=$(( 10 * WORKERS ))
if (( NREC % MAXBCH != 0 )); then
    echo "FATAL: ${NREC} records do not divide into ${MAXBCH} batches (10 x ${WORKERS} workers)"; exit 1
fi
MAXCAS=$(( NREC / MAXBCH ))
SRCNAME=$(awk '/^[[:space:]]*\[ *Source *\]/{a=1} a && /^[[:space:]]*file *=/{print $3; exit}' "${DECK}")
[[ -n "${SRCNAME}" ]] || { echo "FATAL: no file = in the [Source] section of ${DECK}"; exit 1; }

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"
TAG=$(basename "${DECK}" .inp)
cp "${DECK}" "${TAG}.inp"
sed -i "s|^ *maxcas *=.*|    maxcas  = ${MAXCAS}|" "${TAG}.inp"
sed -i "s|^ *maxbch *=.*|    maxbch  = ${MAXBCH}|" "${TAG}.inp"
grep -q '^ *file(1) *=' "${TAG}.inp" \
    || sed -i "s|^ *maxbch *=.*|&\n    file(1) = ${PHITSPATH}|" "${TAG}.inp"

rm -f "${SRCNAME}".[0-9][0-9][0-9]
awk 'BEGIN { srand(20261022) } NF { printf "%.15f\t%s\n", rand(), $0 }' "${SOURCE}" \
    | sort -t "$(printf '\t')" -k1,1 | cut -f2- > shuffled_source.tmp
split -l "$(( NREC / WORKERS ))" -d -a 3 shuffled_source.tmp "${SRCNAME}."
rm -f shuffled_source.tmp
for ((i = WORKERS - 1; i >= 0; i--)); do
    mv "$(printf '%s.%03d' "${SRCNAME}" "$i")" "$(printf '%s.%03d' "${SRCNAME}" $((i + 1)))"
done
for f in "${SRCNAME}".[0-9][0-9][0-9]; do
    n=$(grep -c . "$f")
    (( n == NREC / WORKERS )) || { echo "FATAL: $f has $n records, expected $((NREC / WORKERS))"; exit 1; }
done

echo "file=${TAG}.inp" > phits.in
echo "PHITS ${PHITS_BIN}: ${TAG}, ${NREC} records, ${WORKERS} workers, maxcas ${MAXCAS}, maxbch ${MAXBCH}"
module load intel/2024a
export TZ="${TZ:-:/etc/localtime}"
START=$(date +%s)
mpiexec -n "${SLURM_NTASKS}" -genv OMP_NUM_THREADS 1 "${PHITS_BIN}"
echo "wall time: $(( $(date +%s) - START )) s on ${SLURM_NTASKS} tasks"
rm -f "${SRCNAME}".[0-9][0-9][0-9]
