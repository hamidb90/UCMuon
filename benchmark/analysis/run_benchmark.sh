#!/bin/bash
# Pick the most recently modified Geant4 muons CSV and strip the _muons.csv suffix
LATEST=$(ls -t outputs/*_muons.csv 2>/dev/null | head -1)
if [ -z "$LATEST" ]; then
    echo "ERROR: no *_muons.csv found in outputs/" >&2
    exit 1
fi
GEANT4_PREFIX="${LATEST%_muons.csv}"
echo "Using Geant4 run: $GEANT4_PREFIX"

python3 benchmark_analysis.py \
    --geant4   "$GEANT4_PREFIX" \
    --phits    phits/ \
    --music    music/ \
    --proposal proposal/ \
    --bb       BB/ \
    --ucmuon   UCMuon/ \
    --outdir   figures_benchmark/
