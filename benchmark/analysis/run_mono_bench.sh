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
# =============================================================================
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GUI_DIR="${REPO_ROOT}/gui"
BIN_DIR="${REPO_ROOT}/bin"
SOURCE_FILE="${REPO_ROOT}/benchmark/geant4_muon_rock_v5/sources/benchmark_surface.dat"
OUTDIR="${REPO_ROOT}/benchmark/geant4_muon_rock_v5"

DEPTHS_STR="${1:-1 10 25 50 100 200}"
read -ra DEPTHS <<< "$DEPTHS_STR"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PROP_PY="${HOME}/venvs/ucmuon/bin/python3"; [[ -x "$PROP_PY" ]] || PROP_PY="python3"
PROP_TABLES="${HOME}/.proposal_tables"

echo "============================================================"
echo "  Source : ${SOURCE_FILE}"
echo "  Outdir : ${OUTDIR}"
echo "  Depths : ${DEPTHS[*]} m"
echo "  Started: $(date)"
echo "============================================================"

for DEPTH in "${DEPTHS[@]}"; do
  echo ""
  echo "######################## depth = ${DEPTH} m ########################"

  echo -n "  [BB]       "; t0=$(date +%s)
  python3 "${GUI_DIR}/ucmuon_bb_driver.py" > /dev/null 2>"${OUTDIR}/BB/BB_${DEPTH}m.log" <<EOF
${SOURCE_FILE}
${OUTDIR}/BB/BB_bench_${DEPTH}m.dat
1
14
${DEPTH}
1
1
EOF
  echo "$(( $(date +%s) - t0 ))s"

  echo -n "  [MUSIC]    "; t0=$(date +%s)
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
  echo "$(( $(date +%s) - t0 ))s"

  echo -n "  [PROPOSAL] "; t0=$(date +%s)
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
  echo "$(( $(date +%s) - t0 ))s"

  echo -n "  [UCMuon]   "; t0=$(date +%s)
  python3 "${GUI_DIR}/ucmuon_stochastic_driver.py" > /dev/null 2>"${OUTDIR}/UCMuon/UCMuon_${DEPTH}m.log" <<EOF
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
0
1
EOF
  echo "$(( $(date +%s) - t0 ))s"
done

echo ""
echo "============================================================"
echo "  Finished: $(date)"
echo "============================================================"
