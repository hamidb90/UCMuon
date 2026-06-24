#!/bin/bash
# =============================================================================
#  UCMuon Generator — MPI + OpenMP SLURM batch script
#  Lemaitre4 / CECI cluster
#
#  Usage (from UCMuon/ project root):
#    sbatch hpc/run_ucmuon_gen.sh hpc/input_params.dat
#
#  All per-rank files and the final merged file go into output_<JOBID>/.
#  Nothing lands in the project root during or after the run.
# =============================================================================

#SBATCH --job-name=ucmuon_gen
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=70             # set to free CPUs from `status`
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/ucmuon_gen_%j.out
#SBATCH --error=logs/ucmuon_gen_%j.err
##SBATCH --nodelist=lm4-w027

# -----------------------------------------------------------------------------
module load releases/2023b
module load foss/2023b

export GFORTRAN_UNBUFFERED_ALL=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMPI_MCA_btl_openib_allow_ib=1

cd $SLURM_SUBMIT_DIR
mkdir -p logs

# -----------------------------------------------------------------------------
INPUT=${1:-hpc/input_params.dat}
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file '$INPUT' not found."; exit 1
fi
if [[ ! -x "./bin/ucmuon_gen" ]]; then
    echo "ERROR: ./bin/ucmuon_gen not found. Run: make local"; exit 1
fi

# -----------------------------------------------------------------------------
#  Create output folder BEFORE srun
# -----------------------------------------------------------------------------
RUN_DIR="output_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR"

# Strip # comments from input
CLEAN=$(mktemp /tmp/ucmuon_gen_XXXXXX.dat)
sed 's/#.*//' "$INPUT" | grep -v '^\s*$' > "$CLEAN"

# The last two data lines of the input are the output file stems.
# Replace them unconditionally with paths inside $RUN_DIR so the Fortran
# writes all per-rank files directly into the output folder.
# (This works regardless of whether the stems have .dat or not, and
# regardless of what path the user wrote.)
NLINES=$(wc -l < "$CLEAN")
head -n $(( NLINES - 2 )) "$CLEAN" > "${CLEAN}.tmp"
echo "${RUN_DIR}/ucmuon_surface"  >> "${CLEAN}.tmp"
echo "${RUN_DIR}/ucmuon_selected" >> "${CLEAN}.tmp"
mv "${CLEAN}.tmp" "$CLEAN"

# -----------------------------------------------------------------------------
NRANKS=$SLURM_NTASKS
NTHREADS=$SLURM_CPUS_PER_TASK
echo "============================================================"
echo "  UCMuon Generator  (MPI + OpenMP)"
echo "============================================================"
echo "  Input      : $INPUT  ($(wc -l < "$CLEAN") data lines)"
echo "  Output dir : $RUN_DIR"
echo "  MPI ranks  : $NRANKS"
echo "  OMP/rank   : $NTHREADS"
echo "  Total CPUs : $(( NRANKS * NTHREADS ))"
echo "  Node(s)    : $SLURM_NODELIST"
echo "  Job ID     : $SLURM_JOB_ID"
echo "  Started    : $(date)"
echo "============================================================"
echo ""

# Verify last two lines look correct before launching
echo "  Output stems:"
tail -2 "$CLEAN" | while read line; do echo "    $line"; done
echo ""

# -----------------------------------------------------------------------------
time srun --mpi=pmix ./bin/ucmuon_gen < "$CLEAN"
EXIT_CODE=$?
rm -f "$CLEAN"

echo ""
echo "  Finished   : $(date)"
echo "  Exit code  : $EXIT_CODE"
echo ""

# -----------------------------------------------------------------------------
#  Merge per-rank files — they are in $RUN_DIR
# -----------------------------------------------------------------------------
echo "============================================================"
echo "  Merging output files in $RUN_DIR"
echo "============================================================"

merge_dir() {
    local DIR="$1" BASE="$2"
    local FILES=()
    for r in $(seq 0 $((NRANKS-1))); do
        local F; F=$(printf "${DIR}/${BASE}_%05d.dat" "$r")
        [[ -f "$F" ]] && FILES+=("$F")
    done
    if [[ ${#FILES[@]} -eq 0 ]]; then
        echo "  [skip] no files matching ${DIR}/${BASE}_*.dat"
        return
    fi
    cat "${FILES[@]}" > "${DIR}/${BASE}.dat"
    echo "  Merged ${#FILES[@]} files -> ${DIR}/${BASE}.dat  ($(wc -l < "${DIR}/${BASE}.dat") lines)"
    rm -f "${FILES[@]}"
}

merge_dir "$RUN_DIR" "ucmuon_surface"
merge_dir "$RUN_DIR" "ucmuon_selected"

# -----------------------------------------------------------------------------
#  PHITS conversion
# -----------------------------------------------------------------------------
echo ""
if [[ -x "./bin/ucmuon_to_phits" ]]; then
    for dat in ucmuon_selected ucmuon_surface; do
        if [[ -f "$RUN_DIR/${dat}.dat" ]]; then
            echo "  Converting ${dat}.dat -> PHITS format..."
            ./bin/ucmuon_to_phits gen < "$RUN_DIR/${dat}.dat" \
                              > "$RUN_DIR/${dat}_phits.dat"
            echo "  -> $RUN_DIR/${dat}_phits.dat"
        fi
    done
else
    echo "  NOTE: ucmuon_to_phits not found. Build: make ucmuon_to_phits"
fi

# -----------------------------------------------------------------------------
#  Geant4 conversion  (optional — set UCMUON_GEANT4=1 to enable)
# -----------------------------------------------------------------------------
if [[ "${UCMUON_GEANT4:-0}" == "1" ]]; then
    echo ""
    GEANT4_PY="src/converters/ucmuon_to_geant4.py"
    if [[ -f "$GEANT4_PY" ]]; then
        for dat in ucmuon_selected ucmuon_surface; do
            if [[ -f "$RUN_DIR/${dat}.dat" ]]; then
                echo "  Converting ${dat}.dat -> Geant4 ASCII format..."
                python3 "$GEANT4_PY" "$RUN_DIR/${dat}.dat" \
                                     "$RUN_DIR/${dat}_geant4.txt"
                echo "  -> $RUN_DIR/${dat}_geant4.txt"
            fi
        done
    else
        echo "  NOTE: $GEANT4_PY not found."
    fi
fi

cp "$INPUT" "$RUN_DIR/input_used.dat"

echo ""
echo "  All outputs in: $RUN_DIR/"
echo "  Next step — set line 1 of hpc/input_transport_*.dat to:"
echo "    $RUN_DIR/ucmuon_selected.dat"

# -----------------------------------------------------------------------------
#  Statistics
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Run statistics"
echo "============================================================"
OUTFILE="logs/ucmuon_gen_${SLURM_JOB_ID}.out"
if grep -q "UCMuon_gen  COMPLETE" "$OUTFILE" 2>/dev/null; then
    echo "  Status: COMPLETED NORMALLY"
    grep -E "Total saved|Total tried|Acceptance rate|Parallelisation" "$OUTFILE" | tail -6
else
    echo "  Status: DID NOT COMPLETE"
    grep "Saved .* / tried" "$OUTFILE" | tail -$((NRANKS * 2))
fi
echo "============================================================"
