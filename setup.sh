#!/usr/bin/env bash
# =============================================================================
#  setup.sh  —  UCMuon build & environment check
#  UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
#  MIT License 2026
#
#  Builds the LOCAL (OMP-only) binaries used by the Streamlit GUI.
#  Also creates a system Python venv with PROPOSAL for Engine 3.
#  For HPC (MPI+OMP) builds see: hpc/README_HPC.md
#
#  Usage:
#    bash setup.sh                    full build + Python setup
#    bash setup.sh --no-python        skip all Python steps
#    bash setup.sh --no-proposal      skip PROPOSAL venv setup
#    bash setup.sh --python=python3   use a specific interpreter for core pkgs
# =============================================================================
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

NO_PY=0; NO_PROP=0; PYTHON="python3"
for a in "$@"; do case $a in
  --no-python)   NO_PY=1 ;;
  --no-proposal) NO_PROP=1 ;;
  --python=*)    PYTHON="${a#*=}" ;;
esac; done

VENV="$HOME/venvs/ucmuon"
SYS_PY="/usr/bin/python3"

echo ""
echo "======================================================="
echo "  UCMuon  —  setup & local build"
echo "  UCLouvain Muography Group"
echo "======================================================="
echo ""

# ── 1. gfortran ───────────────────────────────────────────────────────────────
echo "[1/6] Checking gfortran..."
if ! command -v gfortran &>/dev/null; then
    echo ""
    echo "  ERROR: gfortran not found. Install it first:"
    echo ""
    echo "    macOS (Homebrew) : brew install gcc"
    echo "    macOS (MacPorts) : sudo port install gcc14"
    echo "    Ubuntu / Debian  : sudo apt install gfortran"
    echo "    RHEL / Rocky     : sudo yum install gcc-gfortran"
    echo "    HPC (Lemaitre4)  : module load releases/2023b && module load foss/2023b"
    echo ""
    exit 1
fi
echo "  OK  $(gfortran --version | head -1)"

# ── 2. OpenMP ─────────────────────────────────────────────────────────────────
echo "[2/6] Checking OpenMP..."
printf 'program t\n  use omp_lib\n  implicit none\n  print*,omp_get_max_threads()\nend program\n' \
    > /tmp/_ucmuon_omp.f90
if ! gfortran -fopenmp /tmp/_ucmuon_omp.f90 -o /tmp/_ucmuon_omp 2>/dev/null; then
    echo ""
    echo "  ERROR: gfortran does not support -fopenmp."
    echo "  macOS: the Apple-bundled gcc is a clang alias with no OpenMP."
    echo "  Install Homebrew or MacPorts GCC:"
    echo "    brew install gcc   OR   sudo port install gcc14"
    echo ""
    rm -f /tmp/_ucmuon_omp.f90 /tmp/_ucmuon_omp; exit 1
fi
NCORES=$(/tmp/_ucmuon_omp 2>/dev/null | tr -d '[:space:]')
rm -f /tmp/_ucmuon_omp.f90 /tmp/_ucmuon_omp
echo "  OK  OpenMP — ${NCORES} logical core(s) available"

# ── 3. Core Python packages (conda/system — for Engines 1,2,4,5,6) ───────────
if [ "$NO_PY" -eq 0 ]; then
    echo "[3/6] Installing core Python packages..."
    if ! command -v "$PYTHON" &>/dev/null; then
        echo "  ERROR: '$PYTHON' not found. Try: bash setup.sh --python=python3.11"
        exit 1
    fi
    echo "  Using $($PYTHON --version 2>&1)"
    "$PYTHON" -m pip install -q -r requirements.txt \
        && echo "  OK  core requirements installed" \
        || echo "  WARN: pip install had warnings"

    # Optional: rasterio for Engine 6
    if "$PYTHON" -c "import rasterio" 2>/dev/null; then
        echo "  OK  rasterio  (Engine 6: UCMuon Terrain enabled)"
    else
        echo "  --  rasterio not installed  (Engine 6 disabled)"
        echo "      To enable: pip install rasterio"
    fi

    # Engine 4 self-test
    if [ -f "gui/ucmuon_stochastic_driver.py" ]; then
        "$PYTHON" -c "
import importlib.util
spec = importlib.util.spec_from_file_location('drv', 'gui/ucmuon_stochastic_driver.py')
drv  = importlib.util.module_from_spec(spec); spec.loader.exec_module(drv)
d = float(drv._dedx(1000.0))
assert 1.70 < d < 1.95, f'dE/dx out of range: {d:.4f}'
print(f'  OK  Engine 4 self-test  (dE/dx at 1 GeV = {d:.4f} MeV cm2/g)')
" 2>/dev/null || echo "  WARN: Engine 4 self-test failed"
    fi
else
    echo "[3/6] Core Python check skipped (--no-python)"
fi

# ── 4. PROPOSAL venv (Engine 3) ───────────────────────────────────────────────
# PROPOSAL segfaults under Anaconda/miniforge due to a pybind11 ABI mismatch.
# We create a separate venv based on the macOS system Python (/usr/bin/python3)
# which uses the OS C++ runtime that PROPOSAL was compiled against.
# The run_gui.sh launcher always activates this venv automatically.
if [ "$NO_PY" -eq 0 ] && [ "$NO_PROP" -eq 0 ]; then
    echo "[4/6] Setting up PROPOSAL venv (Engine 3)..."

    if [ ! -f "$SYS_PY" ]; then
        echo "  --  $SYS_PY not found — cannot create PROPOSAL venv"
        echo "      Install Xcode command line tools: xcode-select --install"
    else
        echo "  Using system Python: $($SYS_PY --version 2>&1)"

        # Create venv if it does not exist
        if [ ! -d "$VENV" ]; then
            "$SYS_PY" -m venv "$VENV"
            echo "  OK  Created venv: $VENV"
        else
            echo "  OK  Venv exists:  $VENV"
        fi

        # Install requirements into venv
        "$VENV/bin/pip" install -q --upgrade pip
        "$VENV/bin/pip" install -q -r requirements.txt \
            && echo "  OK  Core packages installed in venv"

        # Install PROPOSAL
        PROP_OK=0
        "$VENV/bin/pip" install -q proposal \
            && PROP_OK=1 \
            || echo "  WARN: PROPOSAL install failed — Engine 3 unavailable"

        if [ "$PROP_OK" -eq 1 ]; then
            # Test in subprocess to catch segfault safely
            "$VENV/bin/python3" -c \
                "import proposal; print(proposal.__version__)" \
                > /tmp/_prop_ver.txt 2>/dev/null && PROP_TEST=1 || PROP_TEST=0
            if [ "$PROP_TEST" -eq 1 ] && [ -s /tmp/_prop_ver.txt ]; then
                VER=$(cat /tmp/_prop_ver.txt | tr -d '[:space:]')
                echo "  OK  PROPOSAL ${VER} installed — Engine 3 enabled"
            else
                echo "  WARN: PROPOSAL installed but import test failed"
                echo "        Engine 3 may still be unavailable on this system"
            fi
            rm -f /tmp/_prop_ver.txt
        fi
    fi
else
    echo "[4/6] PROPOSAL venv skipped"
fi

# ── 5. Source file inventory ──────────────────────────────────────────────────
echo "[5/6] Checking source files..."

MISS_GEN=0; MISS_BB=0; MISS_MUSIC=0

for f in geom_module.f90 phits_module.f90 rng_parallel.f90 \
          ucmuon_source_module.f90 ucmuon_gen_omp.f90; do
    [ -f "src/generator/$f" ] \
        && echo "  OK  src/generator/$f" \
        || { echo "  --  src/generator/$f  (MISSING)"; MISS_GEN=1; }
done
[ -f "src/parma/cosmicray.f90" ] \
    && echo "  OK  src/parma/cosmicray.f90" \
    || { echo "  --  src/parma/cosmicray.f90  (MISSING)"; MISS_GEN=1; }
for f in ranlux.f ranlux_omp.f ranmar_omp.f rnorml.f corgen.f90 corset.f; do
    [ -f "src/common/$f" ] \
        && echo "  OK  src/common/$f" \
        || { echo "  --  src/common/$f  (MISSING)"; MISS_GEN=1; }
done
[ -f "src/transport/bethe_bloch/ucmuon_transport_bb_omp.f90" ] \
    && echo "  OK  src/transport/bethe_bloch/ucmuon_transport_bb_omp.f90" \
    || { echo "  --  ucmuon_transport_bb_omp.f90  (MISSING)"; MISS_BB=1; }
if [ -f "src/transport/music/music.f" ] && \
   [ -f "src/transport/music/music-crosssections.f" ] && \
   [ -f "src/transport/music/ucmuon_transport_music_omp.f90" ]; then
    echo "  OK  src/transport/music/  (Engine 1: MUSIC available)"
else
    echo "  --  src/transport/music/music.f  (Engine 1 disabled)"
    echo "      See docs/MUSIC_FILES.md"
    MISS_MUSIC=1
fi
for f in data/music-eloss-rock.dat data/music-double-diff-rock.dat; do
    [ -f "$f" ] && echo "  OK  $f" || echo "  --  $f  (MISSING)"
done

# ── 6. Build ──────────────────────────────────────────────────────────────────
echo "[6/6] Building local (OMP-only) binaries..."
echo ""
make data-links
[ "$MISS_GEN"   -eq 0 ] && make ucmuon_gen_omp              || echo "  SKIP  ucmuon_gen_omp"
[ "$MISS_MUSIC" -eq 0 ] && make ucmuon_transport_music_omp  || echo "  SKIP  ucmuon_transport_music_omp"
[ "$MISS_BB"    -eq 0 ] && make ucmuon_transport_bb_omp     || echo "  SKIP  ucmuon_transport_bb_omp"

# ── Summary ───────────────────────────────────────────────────────────────────
# Determine PROPOSAL status for summary (test in subprocess to catch segfault)
PROP_SUMMARY=0
if [ -d "$VENV" ]; then
    "$VENV/bin/python3" -c "import proposal; print('ok')" \
        2>/dev/null | grep -q ok && PROP_SUMMARY=1 || true
fi

echo ""
echo "======================================================="
echo "  ENGINE AVAILABILITY"
echo "======================================================="
echo ""
[ -f "gui/ucmuon_stochastic_driver.py" ] \
    && echo "  [x] Engine 1  UCMuon-MC (flagship)       (Python)" \
    || echo "  [ ] Engine 1  UCMuon-MC (flagship)       (gui/ucmuon_stochastic_driver.py missing)"
[ -f bin/ucmuon_transport_music_omp ] \
    && echo "  [x] Engine 2  MUSIC stochastic MC        (bin/ucmuon_transport_music_omp)" \
    || echo "  [ ] Engine 2  MUSIC                      (music.f not found — see docs/MUSIC_FILES.md)"
[ -f bin/ucmuon_transport_bb_omp ] \
    && echo "  [x] Engine 3  Bethe-Bloch + Highland MS  (bin/ucmuon_transport_bb_omp)" \
    || echo "  [ ] Engine 3  Bethe-Bloch                (build failed)"
[ "$PROP_SUMMARY" -eq 1 ] \
    && echo "  [x] Engine 4  PROPOSAL stochastic MC     (system Python venv)" \
    || echo "  [ ] Engine 4  PROPOSAL                   (venv at $VENV not ready)"
[ -f "gui/ucmuon_backward_mc.py" ] \
    && echo "  [x] Engine 5  Backward MC                (Python)" \
    || echo "  [ ] Engine 5  Backward MC                (gui/ucmuon_backward_mc.py missing)"
"${PYTHON}" -c "import rasterio" 2>/dev/null \
    && echo "  [x] Engine 6  UCMuon Terrain             (Python + rasterio)" \
    || echo "  [ ] Engine 6  UCMuon Terrain             (pip install rasterio)"
echo ""
echo "  Intermediate files (.o, .mod) are in build/"
echo ""
echo "======================================================="
echo "  LAUNCH THE GUI"
echo "======================================================="
echo ""
echo "    bash run_gui.sh"
echo ""
echo "  run_gui.sh activates the correct Python environment"
echo "  automatically — no manual venv activation needed."
echo ""
echo "  ── HPC (MPI+OMP) build ──────────────────────────────"
echo ""
echo "    module load releases/2023b && module load foss/2023b"
echo "    make hpc"
echo "    See hpc/README_HPC.md for full workflow."
echo ""
if [ "$MISS_MUSIC" -eq 1 ]; then
    echo "  NOTE: Engine 1 unavailable — music.f not found."
    echo "  Engines 2–6 are fully functional."
    echo "  See docs/MUSIC_FILES.md to enable Engine 1."
    echo ""
fi
