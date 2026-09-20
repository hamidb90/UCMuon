#!/bin/bash
# =============================================================================
# run_scaling_hpc.sh — OpenMP scaling measurement for fig07 / tab:scaling
# on an HPC node (e.g. CECI Lemaitre4, 32-core AMD EPYC).
#
# Run from the UCMuon repo root on an exclusive compute node:
#   sbatch --exclusive --cpus-per-task=32 manuscript/scripts/run_scaling_hpc.sh
# or interactively:
#   bash manuscript/scripts/run_scaling_hpc.sh
#
# Produces manuscript/scripts/scaling_hpc.csv, which
# make_fig07_scaling.py picks up automatically.
# =============================================================================
#SBATCH --job-name=ucmuon_scaling
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00

set -euo pipefail

# Under `sbatch`, Slurm stages a copy of this script into its spool directory,
# so $0 no longer points into the repo; prefer SLURM_SUBMIT_DIR (the directory
# sbatch was invoked from, i.e. the repo root per this script's usage note).
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
SCRIPTS="$REPO_ROOT/manuscript/scripts"
WORK="$SCRIPTS/_hpc_scaling_work"
CSV="$SCRIPTS/scaling_hpc.csv"
THREADS="${THREADS:-1 2 4 8 16 32}"

mkdir -p "$WORK"
for t in music-eloss-rock.dat music-cross-sections-rock.dat \
         music-double-diff-rock.dat; do
  [ -e "$WORK/$t" ] || ln -s "$REPO_ROOT/bin/$t" "$WORK/$t"
done

# NOTE: each input ends with an explicit trailing newline (the `$'\n'` appended
# after the $() ) to feed the binaries' final "Press Enter to start..." read.
# Command substitution strips trailing newlines, so a trailing '' inside printf
# is NOT enough — the blank line must be re-added outside the substitution.
#
# Generator (timed): detector-filtered config of hpc/input_params.dat —
# 2000 SELECTED muons through two 4 cm cylinder detectors (compute-bound).
GEN_IN=$(printf '%s\n' 0 19.0 1500.0 1 1 1 0.0 0.0 500.0 0.0 0.0 0.0 \
                       2 85.0 2000 1 2 \
                       1 3.0 '990.0 0.0 -9000.0' '990.0 0.0 -3500.0' 4.0 \
                       1 3.0 '0.0 -1650.0 -9000.0' '0.0 -1650.0 -3500.0' 4.0 \
                       0 0 muons_surface_unused.dat muons_selected.dat)$'\n'
# MUSIC input prep (not timed): 1e6 muons, 10-2500 GeV, no detector
GEN_PREP=$(printf '%s\n' 0 10.0 2500.0 1 1 1 0.0 0.0 100.0 0.0 0.0 0.0 \
                         2 85.0 1000000 0 1 0 muons_surface.dat)$'\n'
# MUSIC stdin (timed): Standard Rock, d = 150 m, tables from disk
MUS_IN=$(printf '%s\n' muons_surface.dat muons_underground.dat \
                       2.65 26.48 150.0 1 1 -30 1 1 1)$'\n'

cd "$WORK"

# MUSIC input file + warm-up (cross-section table cache, page cache)
[ -f muons_surface.dat ] || \
  "$REPO_ROOT/bin/ucmuon_gen_omp" <<< "$GEN_PREP" > /dev/null
OMP_NUM_THREADS=2 "$REPO_ROOT/bin/ucmuon_gen_omp"            <<< "$GEN_IN" > /dev/null
OMP_NUM_THREADS=2 "$REPO_ROOT/bin/ucmuon_transport_music_omp" <<< "$MUS_IN" > /dev/null

echo "program,threads,seconds" > "$CSV"
for k in $THREADS; do
  for prog in generator music; do
    if [ "$prog" = generator ]; then
      exe="$REPO_ROOT/bin/ucmuon_gen_omp";            input="$GEN_IN"
    else
      exe="$REPO_ROOT/bin/ucmuon_transport_music_omp"; input="$MUS_IN"
    fi
    start=$(date +%s.%N)
    OMP_NUM_THREADS=$k "$exe" <<< "$input" > /dev/null
    end=$(date +%s.%N)
    secs=$(echo "$end - $start" | bc)
    echo "$prog,$k,$secs" >> "$CSV"
    echo "  $prog k=$k  ${secs}s"
  done
done

echo "[OK] wrote $CSV — re-run make_fig07_scaling.py to update the figure"
