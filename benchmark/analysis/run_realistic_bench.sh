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
#
#  Driver stdin schemas verified against current gui/*.py (2026-06-14).
# =============================================================================
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GUI_DIR="${REPO_ROOT}/gui"
BIN_DIR="${REPO_ROOT}/bin"
SOURCE_FILE="${REPO_ROOT}/output/muons_surface.dat"
OUTDIR="${REPO_ROOT}/benchmark/run_realistic_20260614"

DEPTHS_STR="${1:-1 10 25 50 100 200}"
read -ra DEPTHS <<< "$DEPTHS_STR"

# single-thread for the Fortran/PROPOSAL engines (UCMuon-MC uses its own param)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PROP_PY="${HOME}/venvs/ucmuon/bin/python3"; [[ -x "$PROP_PY" ]] || PROP_PY="python3"
PROP_TABLES="${HOME}/.proposal_tables"

mkdir -p "${OUTDIR}/BB" "${OUTDIR}/MUSIC" "${OUTDIR}/PROPOSAL" "${OUTDIR}/UCMuon"

echo "============================================================"
echo "  Source : ${SOURCE_FILE}"
echo "  Outdir : ${OUTDIR}"
echo "  Depths : ${DEPTHS[*]} m"
echo "  Started: $(date)"
echo "============================================================"

for DEPTH in "${DEPTHS[@]}"; do
  echo ""
  echo "######################## depth = ${DEPTH} m ########################"

  # ── BB (Python CSDA driver) : infile,outfile,transport_all,ncols,depth,mat,ms
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

  # ── MUSIC (Fortran, run from bin/ for table files)
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

  # ── PROPOSAL (Moliere MCS) : infile,outfile,depth,mat,transport_all,e_cut,v_cut,scat,tables
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

  # ── UCMuon-MC (current v2 per-process loss model; auto workers for speed)
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
echo "  Outputs : ${OUTDIR}/{BB,MUSIC,PROPOSAL,UCMuon}/"
echo "============================================================"
