#!/bin/bash
# =============================================================================
#  UCMuon 4-engine benchmark on the MONOENERGETIC v2 source
#  (6 vertical beams: 5/10/20/50/100/300 GeV x 100k muons each = 600k)
#
#  Regenerates benchmark/geant4_muon_rock_v5/{BB,MUSIC,PROPOSAL,UCMuon}/
#  engine outputs against the fixed Geant4/PHITS references in the same tree.
#  Feeds make_fig06.py and make_survival_table.py.
#
#  Usage:
#    bash run_mono_bench.sh "1 10 25 50 100 200"    # depths (m)
#    bash run_mono_bench.sh 1                       # smoke test
#    ENGINES="BB UCMuon" bash run_mono_bench.sh "1 10 25 50 100 200"
#                                                   # only those engines
#    OUTDIR=/tmp/t UCMUON_WORKERS=1 bash run_mono_bench.sh
#                       # timing run: outputs elsewhere, UCMuon-MC in one process
#    SOURCE_FILE=.../benchmark_surface_th85.dat INCL_DEG=85 OUTDIR=.../th85 \
#      bash run_mono_bench.sh   # inclined beams (make_inclined_inputs.py): the
#                       # depths are slant lengths L, the engines get L cos(85 deg)
#
#  Every run's wall time goes to ${OUTDIR}/timing.csv (engine, depth_m,
#  seconds, processes). UCMUON_WORKERS is UCMuon-MC's process count (stdin
#  item 19): 0 (default) is automatic, one process per 2x10^4 muons up to the
#  core count, as the archived outputs were made; 1 is a single process, as for
#  the single-thread times in the paper. The other engines run single-threaded.
# =============================================================================
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GUI_DIR="${REPO_ROOT}/gui"
BIN_DIR="${REPO_ROOT}/bin"
SOURCE_FILE="${SOURCE_FILE:-${REPO_ROOT}/benchmark/geant4_muon_rock_v5/sources/benchmark_surface.dat}"
INCL_DEG="${INCL_DEG:-0}"
OUTDIR="${OUTDIR:-${REPO_ROOT}/benchmark/geant4_muon_rock_v5}"
UCMUON_WORKERS="${UCMUON_WORKERS:-0}"

DEPTHS_STR="${1:-1 10 25 50 100 200}"
read -ra DEPTHS <<< "$DEPTHS_STR"

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
echo "  Depths : ${DEPTHS[*]} m$( [[ "$INCL_DEG" != 0 ]] && echo " (slant; beams inclined at ${INCL_DEG} deg)")"
echo "  Engines: ${ENGINES}"
echo "  Started: $(date)"
echo "============================================================"

for DEPTH in "${DEPTHS[@]}"; do
  echo ""
  echo "######################## depth = ${DEPTH} m ########################"
  # The engines take the vertical depth; DEPTH (file names, timing.csv) is the
  # length along the beam, the same for vertical beams.
  if [[ "$INCL_DEG" == 0 ]]; then VDEPTH="$DEPTH"
  else VDEPTH=$(perl -e "printf '%.6f', ${DEPTH} * cos(${INCL_DEG} * atan2(1, 1) / 45)"); fi

  if want BB; then
  echo -n "  [BB]       "; t0=$(now)
  python3 "${GUI_DIR}/ucmuon_bb_driver.py" > /dev/null 2>"${OUTDIR}/BB/BB_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/BB/BB_bench_${DEPTH}m.dat
1
14
${VDEPTH}
1
1
EOF
  lap BB "$t0" 1
  fi

  if want MUSIC; then
  echo -n "  [MUSIC]    "; t0=$(now)
  ( cd "${BIN_DIR}" && "${BIN_DIR}/ucmuon_transport_music_omp" > /dev/null 2>"${OUTDIR}/MUSIC/MUSIC_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/MUSIC/MUSIC_bench_${DEPTH}m.dat
2.65
26.48
${VDEPTH}
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

  if want PROPOSAL; then
  echo -n "  [PROPOSAL] "; t0=$(now)
  "${PROP_PY}" "${GUI_DIR}/proposal_driver.py" > /dev/null 2>"${OUTDIR}/PROPOSAL/PROPOSAL_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/PROPOSAL/PROPOSAL_bench_${DEPTH}m.dat
${VDEPTH}
1
1
500
0.001
3
${PROP_TABLES}
EOF
  lap PROPOSAL "$t0" 1
  fi

  if want UCMuon; then
  echo -n "  [UCMuon]   "; t0=$(now)
  python3 "${GUI_DIR}/ucmuon_stochastic_driver.py" > "${OUTDIR}/UCMuon/UCMuon_${DEPTH}m.log" 2>&1 <<EOF
${SOURCE_FILE}
${OUTDIR}/UCMuon/UCMuon_bench_${DEPTH}m.dat
${VDEPTH}
2.65
10.015
1
1
14
0
0.01
1
11
22
136.4
3.475e-6
1
3
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
echo "============================================================"
