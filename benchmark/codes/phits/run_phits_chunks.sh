#!/bin/bash
# =============================================================================
# run_phits_chunks.sh: a PHITS benchmark deck on one machine, as K independent
# single-thread runs of PHITS's own binary, P at a time, merged into one summary.
#
#   PHITSPATH=/path/to/phits bash run_phits_chunks.sh DECK SOURCE RUN_DIR [K] [P]
#
# Chunk c (0 ... K-1) gets the source records c+1, c+1+K, ..., so that every
# chunk sees the whole population, and the random seed rseed = 1000 + c: with
# one seed, chunks of a source sorted by beam would repeat the same histories.
# Each chunk is parsed by read_phits_output.py; merge_chunks.py combines them
# into RUN_DIR/phits_summary.csv. K defaults to 32 and P to 8. The record count
# per chunk must divide into 10 batches.
#
# On a cluster use run_phits_mpi.sh instead. The inclined-beam summaries of
# 2026-10-27 were made this way on an Apple M4 (K = 32, P = 8).
# =============================================================================
set -euo pipefail
DECK=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); SOURCE=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
RUN_DIR=$3; K=${4:-32}; P=${5:-8}
: "${PHITSPATH:?set PHITSPATH to the PHITS installation (the directory holding bin/ and data/)}"
HERE=$(cd "$(dirname "$0")" && pwd)
case "$(uname -s)" in Darwin) BIN=${PHITS_BIN:-$PHITSPATH/bin/phits_mac} ;; *) BIN=${PHITS_BIN:-$PHITSPATH/bin/phits_lin} ;; esac
[[ -x "$BIN" ]] || { echo "FATAL: $BIN not found (set PHITS_BIN)"; exit 1; }
SRCNAME=$(awk '/^[[:space:]]*\[ *Source *\]/{a=1} a && /^[[:space:]]*file *=/{print $3; exit}' "$DECK")
TAG=$(basename "$DECK")

mkdir -p "$RUN_DIR"
for c in $(seq 0 $((K - 1))); do
    d=$RUN_DIR/c$(printf %03d "$c"); mkdir -p "$d"
    awk -v c="$c" -v k="$K" 'NF && (NR - 1) % k == c' "$SOURCE" > "$d/$SRCNAME"
    n=$(grep -c . "$d/$SRCNAME"); echo "$n" > "$d/n_source.txt"
    (( n % 10 == 0 )) || { echo "FATAL: chunk $c has $n records, not a multiple of 10"; exit 1; }
    python3 - "$DECK" "$d/$TAG" "$((n / 10))" "$c" "$PHITSPATH" <<'EOF'
import re, sys
src, dst, cas, c, path = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
s = open(src).read()
s = re.sub(r"(?m)^ *maxcas *=.*$", f"    maxcas  = {cas}", s, count=1)
s = re.sub(r"(?m)^ *maxbch *=.*$", f"    maxbch  = 10\n    file(1) = {path}", s, count=1)
s, k = re.subn(r"(?m)^ *rseed *=.*$", f"    rseed   = {1000 + c}", s, count=1)
assert k == 1, "the deck has no rseed line"
s = re.sub(r"(?m)^\$OMP=.*$", "", s, count=1)
open(dst, "w").write(s)
EOF
    echo "file=$TAG" > "$d/phits.in"
done
echo "$TAG: $K chunks of $(cat "$RUN_DIR/c000/n_source.txt") records, $P at a time, $(date)"
export BIN
ls -d "$RUN_DIR"/c??? | xargs -P "$P" -n 1 bash -c 'cd "$1" && "$BIN" < phits.in > run.log 2>&1 || echo "  chunk $1 failed"' _
for d in "$RUN_DIR"/c???; do
    (cd "$d" && python3 "$HERE/read_phits_output.py" --n "$(cat n_source.txt)" > parse.log 2>&1)
done
python3 "$HERE/merge_chunks.py" "$RUN_DIR/phits_summary.csv" "$RUN_DIR"/c???
echo "done $(date)"
