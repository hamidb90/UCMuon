#!/bin/bash
# =============================================================================
#  UCMuon 4-engine benchmark on a REALISTIC generated source
#  Runs BB, MUSIC, PROPOSAL, UCMuon-MC through Standard Rock at several depths.
#
#  Differs from the canonical monoenergetic benchmark: the source here is a
#  realistic cosmic-ray distribution (output/muons_surface.dat), so Geant4/PHITS
#  references (made on the monoenergetic grid) are NOT comparable and excluded.
#
#  Usage:
#    bash run_realistic_bench.sh "1 10 25 50 100 200"     # depths (m), space-sep
#    bash run_realistic_bench.sh 1                        # smoke test, 1 m only
#    ENGINES="BB UCMuon" bash run_realistic_bench.sh "1 10 25 50 100 200"
#                                                         # only those engines
#    OUTDIR=/tmp/t UCMUON_WORKERS=1 bash run_realistic_bench.sh
#                       # timing run: outputs elsewhere, UCMuon-MC in one process
#
#  Every run's wall time goes to ${OUTDIR}/timing.csv (engine, depth_m,
#  seconds, processes). UCMUON_WORKERS is UCMuon-MC's process count (stdin
#  item 19): 0 (default) is automatic, one process per 2x10^4 muons up to the
#  core count, as the archived outputs were made; 1 is a single process, as for
#  the single-thread times in the paper. The other engines run single-threaded.
#
#  Driver stdin schemas verified against current gui/*.py (2026-06-14).
# =============================================================================
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GUI_DIR="${REPO_ROOT}/gui"
BIN_DIR="${REPO_ROOT}/bin"
# Dedicated benchmark source (CosmoALEPH 100-2500 GeV, disk r=200 m, cos^2,
# theta_max=85 deg, N=100k — same recipe as the 2026-06-14 run).  Kept inside
# benchmark/ so GUI runs that rewrite output/muons_surface.dat cannot silently
# change the benchmark input.
SOURCE_FILE="${REPO_ROOT}/benchmark/sources/muons_surface_real100.dat"
OUTDIR="${OUTDIR:-${REPO_ROOT}/benchmark/run_realistic_20260614}"
UCMUON_WORKERS="${UCMUON_WORKERS:-0}"

DEPTHS_STR="${1:-1 10 25 50 100 200}"
read -ra DEPTHS <<< "$DEPTHS_STR"

# single-thread for the Fortran/PROPOSAL engines (UCMuon-MC uses its own param)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PROP_PY="${HOME}/venvs/ucmuon/bin/python3"; [[ -x "$PROP_PY" ]] || PROP_PY="python3"
PROP_TABLES="${HOME}/.proposal_tables"

ENGINES="${ENGINES:-BB MUSIC PROPOSAL UCMuon}"
want() { [[ " ${ENGINES} " == *" $1 "* ]]; }

# Wall time of a run, to 0.01 s, printed and appended to timing.csv.
now() { perl -MTime::HiRes=time -e 'printf "%.3f", time'; }
lap() {   # lap ENGINE T0 PROCESSES
  local s; s=$(perl -e "printf '%.2f', $(now) - $2")
  echo "${s}s"
  echo "$1,${DEPTH},${s},$3" >> "${OUTDIR}/timing.csv"
}
mkdir -p "${OUTDIR}/BB" "${OUTDIR}/MUSIC" "${OUTDIR}/PROPOSAL" "${OUTDIR}/UCMuon"
[[ -f "${OUTDIR}/timing.csv" ]] || echo "engine,depth_m,seconds,processes" > "${OUTDIR}/timing.csv"

echo "============================================================"
echo "  Source : ${SOURCE_FILE}"
echo "  Outdir : ${OUTDIR}"
echo "  Depths : ${DEPTHS[*]} m"
echo "  Engines: ${ENGINES}"
echo "  Started: $(date)"
echo "============================================================"

for DEPTH in "${DEPTHS[@]}"; do
  echo ""
  echo "######################## depth = ${DEPTH} m ########################"

  # ── BB (Python CSDA driver) : infile,outfile,transport_all,ncols,depth,mat,ms
  if want BB; then
  echo -n "  [BB]       "; t0=$(now)
  python3 "${GUI_DIR}/ucmuon_bb_driver.py" > /dev/null 2>"${OUTDIR}/BB/BB_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/BB/BB_bench_${DEPTH}m.dat
1
14
${DEPTH}
1
1
EOF
  lap BB "$t0" 1
  fi

  # ── MUSIC (Fortran, run from bin/ for table files)
  if want MUSIC; then
  echo -n "  [MUSIC]    "; t0=$(now)
  ( cd "${BIN_DIR}" && "${BIN_DIR}/ucmuon_transport_music_omp" > /dev/null 2>"${OUTDIR}/MUSIC/MUSIC_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/MUSIC/MUSIC_bench_${DEPTH}m.dat
2.65
26.48
${DEPTH}
1
1
-30
1
1
1

EOF
  )
  lap MUSIC "$t0" 1
  fi

  # ── PROPOSAL (Moliere MCS) : infile,outfile,depth,mat,transport_all,e_cut,v_cut,scat,tables
  if want PROPOSAL; then
  echo -n "  [PROPOSAL] "; t0=$(now)
  "${PROP_PY}" "${GUI_DIR}/proposal_driver.py" > /dev/null 2>"${OUTDIR}/PROPOSAL/PROPOSAL_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/PROPOSAL/PROPOSAL_bench_${DEPTH}m.dat
${DEPTH}
1
1
500
0.001
3
${PROP_TABLES}
EOF
  lap PROPOSAL "$t0" 1
  fi

  # ── UCMuon-MC (current v2 per-process loss model; auto workers for speed)
  if want UCMuon; then
  echo -n "  [UCMuon]   "; t0=$(now)
  python3 "${GUI_DIR}/ucmuon_stochastic_driver.py" > "${OUTDIR}/UCMuon/UCMuon_${DEPTH}m.log" 2>&1 <<EOF
${SOURCE_FILE}
${OUTDIR}/UCMuon/UCMuon_bench_${DEPTH}m.dat
${DEPTH}
2.65
10.015
1
1
14
0
0.05
1
11
22
136.4
3.475e-6
1
2
42
${UCMUON_WORKERS}
1
EOF
  lap UCMuon "$t0" "$(sed -n 's/^ *Workers: *\([0-9]*\).*/\1/p' "${OUTDIR}/UCMuon/UCMuon_${DEPTH}m.log")"
  fi
done

echo ""
echo "============================================================"
echo "  Finished: $(date)"
echo "  Outputs : ${OUTDIR}/{BB,MUSIC,PROPOSAL,UCMuon}/"
echo "============================================================"
