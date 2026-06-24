#!/bin/bash
# =============================================================================
#  UCMuon Transport — MPI + OpenMP SLURM batch script
#  Works for MUSIC and Bethe-Bloch engines.
#  Lemaitre4 / CECI cluster
#
#  Usage (from UCMuon/ project root):
#    MUSIC:        sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_music.dat
#    Bethe-Bloch:  sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_bb.dat
#
#  All per-rank files and the final merged file go into output_<JOBID>/.
# =============================================================================

#SBATCH --job-name=ucmuon_transport
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=26             # set to free CPUs from `status`
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/ucmuon_transport_%j.out
#SBATCH --error=logs/ucmuon_transport_%j.err
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
INPUT=${1:-hpc/input_transport_music.dat}
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file '$INPUT' not found."; exit 1
fi

# Select binary
if [[ "$INPUT" == *"_bb"* ]]; then
    BINARY="./bin/ucmuon_transport_bb"; ENGINE="Bethe-Bloch"
else
    BINARY="./bin/ucmuon_transport_music"; ENGINE="MUSIC"
fi
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: $BINARY not found."; exit 1
fi

# -----------------------------------------------------------------------------
#  Create output folder BEFORE srun
# -----------------------------------------------------------------------------
RUN_DIR="output_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR"

# Strip # comments
CLEAN=$(mktemp /tmp/ucmuon_transport_XXXXXX.dat)
sed 's/#.*//' "$INPUT" | grep -v '^\s*$' > "$CLEAN"
NLINES=$(wc -l < "$CLEAN")

# Line 2 of the cleaned input is the output prefix.
# Replace it unconditionally with the full path inside $RUN_DIR.
OUTPUT_BASENAME=$(sed -n '2p' "$CLEAN" | tr -d '[:space:]')
OUTPUT_BASENAME=$(basename "${OUTPUT_BASENAME:-ucmuon_underground}")
OUTPUT_PREFIX="${RUN_DIR}/${OUTPUT_BASENAME}"
# Rewrite line 2
{
    sed -n '1p' "$CLEAN"
    echo "$OUTPUT_PREFIX"
    sed -n "3,${NLINES}p" "$CLEAN"
} > "${CLEAN}.tmp"
mv "${CLEAN}.tmp" "$CLEAN"

# MUSIC init_tables check (line 10 of cleaned input = init_tables)
if [[ "$ENGINE" == "MUSIC" ]]; then
    INIT_VAL=$(sed -n '10p' "$CLEAN" | tr -d '[:space:]')
    if [[ "$INIT_VAL" == "0" ]]; then
        echo "  init_tables=0: rank 0 will generate MUSIC tables (~1 min)."
        if [[ ! -f "data/music-double-diff-rock.dat" && ! -f "music-double-diff-rock.dat" ]]; then
            echo "ERROR: music-double-diff-rock.dat not found in ./ or ./data/"; exit 1
        fi
    else
        echo "  init_tables=1: reading pre-computed tables."
    fi
fi

# -----------------------------------------------------------------------------
NRANKS=$SLURM_NTASKS
NTHREADS=$SLURM_CPUS_PER_TASK
echo "============================================================"
echo "  UCMuon Transport  ($ENGINE — MPI + OpenMP)"
echo "============================================================"
echo "  Input      : $INPUT  ($(wc -l < "$CLEAN") data lines)"
echo "  Binary     : $BINARY"
echo "  Output dir : $RUN_DIR"
echo "  Output pfx : $OUTPUT_PREFIX"
echo "  MPI ranks  : $NRANKS"
echo "  OMP/rank   : $NTHREADS"
echo "  Total CPUs : $(( NRANKS * NTHREADS ))"
echo "  Node(s)    : $SLURM_NODELIST"
echo "  Job ID     : $SLURM_JOB_ID"
echo "  Started    : $(date)"
echo "============================================================"
echo ""

# -----------------------------------------------------------------------------
time srun --mpi=pmix "$BINARY" < "$CLEAN"
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
echo "  Merging: ${OUTPUT_PREFIX}_RRRRR.dat -> ${OUTPUT_PREFIX}.dat"
echo "============================================================"

FILES=()
for r in $(seq 0 $((NRANKS-1))); do
    F=$(printf "${OUTPUT_PREFIX}_%05d.dat" "$r")
    [[ -f "$F" ]] && FILES+=("$F")
done

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "  WARNING: no per-rank files found matching ${OUTPUT_PREFIX}_*.dat"
else
    cat "${FILES[@]}" > "${OUTPUT_PREFIX}.dat"
    echo "  Merged ${#FILES[@]} files -> ${OUTPUT_PREFIX}.dat  ($(wc -l < "${OUTPUT_PREFIX}.dat") lines)"
    rm -f "${FILES[@]}"
fi

# -----------------------------------------------------------------------------
#  PHITS conversion (survived muons only)
# -----------------------------------------------------------------------------
echo ""
if [[ -x "./bin/ucmuon_to_phits" ]]; then
    if [[ -f "${OUTPUT_PREFIX}.dat" ]]; then
        echo "  Converting underground muons (alive only) -> PHITS..."
        ./bin/ucmuon_to_phits transport < "${OUTPUT_PREFIX}.dat" \
                          > "${OUTPUT_PREFIX}_phits.dat"
        echo "  -> ${OUTPUT_PREFIX}_phits.dat"
    fi
else
    echo "  NOTE: ucmuon_to_phits not found. Build: make ucmuon_to_phits"
fi

# -----------------------------------------------------------------------------
#  Geant4 conversion  (optional — set UCMUON_GEANT4=1 to enable)
# -----------------------------------------------------------------------------
if [[ "${UCMUON_GEANT4:-0}" == "1" ]]; then
    echo ""
    GEANT4_PY="src/converters/ucmuon_to_geant4.py"
    if [[ -f "$GEANT4_PY" ]] && [[ -f "${OUTPUT_PREFIX}.dat" ]]; then
        echo "  Converting underground muons (alive only) -> Geant4 ASCII format..."
        python3 "$GEANT4_PY" "${OUTPUT_PREFIX}.dat" \
                             "${OUTPUT_PREFIX}_geant4.txt" --mode transport
        echo "  -> ${OUTPUT_PREFIX}_geant4.txt"
    elif [[ ! -f "$GEANT4_PY" ]]; then
        echo "  NOTE: $GEANT4_PY not found."
    fi
fi

cp "$INPUT" "$RUN_DIR/input_used.dat"
echo ""
echo "  All outputs in: $RUN_DIR/"

# -----------------------------------------------------------------------------
#  Statistics
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Run statistics"
echo "============================================================"
OUTFILE="logs/ucmuon_transport_${SLURM_JOB_ID}.out"
if grep -q "COMPLETE" "$OUTFILE" 2>/dev/null; then
    echo "  Status: COMPLETED NORMALLY"
    grep -E "Total transported|Survived|Survival rate|Wall time|Depth" "$OUTFILE" | tail -8
else
    echo "  Status: DID NOT COMPLETE"
    grep "Transported:" "$OUTFILE" | tail -$((NRANKS * 2))
fi
echo "============================================================"
