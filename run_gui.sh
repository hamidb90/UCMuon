#!/usr/bin/env bash
# =============================================================================
#  run_gui.sh  —  UCMuon GUI launcher
#  UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
#
#  Always use this script instead of calling streamlit directly.
#  It automatically activates the correct Python environment so that
#  all six engines including PROPOSAL (Engine 3) are available.
#
#  Usage:
#    bash run_gui.sh
#    bash run_gui.sh --threads=4     override OMP_NUM_THREADS
# =============================================================================
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

VENV="$HOME/venvs/ucmuon"
THREADS=""
for a in "$@"; do case $a in --threads=*) THREADS="${a#*=}" ;; esac; done

# ── Activate system Python venv (required for PROPOSAL / Engine 3) ────────────
if [ -d "$VENV" ]; then
    source "$VENV/bin/activate"
    PROP_OK=0
    "$VENV/bin/python3" -c "import proposal" 2>/dev/null && PROP_OK=1 || true
    if [ "$PROP_OK" -eq 1 ]; then
        echo "  OK  Using system Python venv — all 6 engines available"
    else
        echo "  OK  Using system Python venv — Engine 3 (PROPOSAL) not installed"
        echo "      Run: bash setup.sh  to install it"
    fi
else
    echo "  NOTE: System Python venv not found at $VENV"
    echo "        Engine 3 (PROPOSAL) will be unavailable."
    echo "        Run: bash setup.sh  to create the venv automatically."
    echo ""
fi

# ── Set OMP thread count ──────────────────────────────────────────────────────
if [ -n "$THREADS" ]; then
    export OMP_NUM_THREADS="$THREADS"
elif [ -z "$OMP_NUM_THREADS" ]; then
    export OMP_NUM_THREADS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)
fi
echo "  OK  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo ""

# ── Launch ────────────────────────────────────────────────────────────────────
streamlit run gui/ucmuon_gui.py
