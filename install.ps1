# =============================================================================
#  install.ps1  --  UCMuon Windows installer
#  UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
#  MIT License 2026
#
#  Run from the project root in PowerShell:
#    powershell -ExecutionPolicy Bypass -File install.ps1
#
#  What this script does:
#    1. Checks Python >= 3.9
#    2. Installs all Python packages (streamlit, numpy, etc.)
#    3. Looks for gfortran (via MSYS2 or PATH)
#    4. If found: builds the Fortran OMP binaries with make
#    5. If not found: prints clear MSYS2 install instructions
#    6. Prints an engine availability summary
#
#  Engines 4 (UCMuon Stochastic) and 5 (Backward MC) are pure Python
#  and work on Windows without any Fortran compiler.
#  Engines 1 (MUSIC) and 2 (Bethe-Bloch) need gfortran via MSYS2.
#  Engine 3 (PROPOSAL) is not supported on Windows.
# =============================================================================

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host ""
Write-Host "======================================================="
Write-Host "  UCMuon  --  Windows installer"
Write-Host "  UCLouvain Muography Group"
Write-Host "======================================================="
Write-Host ""

# ── Helper ────────────────────────────────────────────────────────────────────
function Find-Exe($name, $extraPaths=@()) {
    $cmd = Get-Command $name -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    foreach ($p in $extraPaths) {
        if (Test-Path $p) { return $p }
    }
    return $null
}

# ── 1. Python ─────────────────────────────────────────────────────────────────
Write-Host "[1/4] Checking Python..."
$py = Find-Exe "python" @("C:\Python311\python.exe","C:\Python312\python.exe","C:\Python310\python.exe")
if (-not $py) {
    Write-Host ""
    Write-Host "  ERROR: Python not found."
    Write-Host ""
    Write-Host "  Install Python 3.11 or later from:"
    Write-Host "    https://www.python.org/downloads/"
    Write-Host ""
    Write-Host "  Or via winget (Windows 10/11):"
    Write-Host "    winget install Python.Python.3.11"
    Write-Host ""
    Write-Host "  Make sure to tick 'Add Python to PATH' during install."
    exit 1
}

# Check version
$ver_out = & $py --version 2>&1
if ($ver_out -match "(\d+)\.(\d+)") {
    $maj = [int]$Matches[1]; $min = [int]$Matches[2]
    if ($maj -lt 3 -or ($maj -eq 3 -and $min -lt 9)) {
        Write-Host "  ERROR: Python >= 3.9 required. Found: $ver_out"
        exit 1
    }
}
Write-Host "  OK  $ver_out  ($py)"

# ── 2. Python packages ────────────────────────────────────────────────────────
Write-Host "[2/4] Installing Python packages..."
& $py -m pip install -q --upgrade pip
& $py -m pip install -q -r requirements.txt
if ($LASTEXITCODE -ne 0) {
    Write-Host "  ERROR: pip install failed. Check your internet connection."
    exit 1
}
Write-Host "  OK  Core packages installed (streamlit, numpy, pandas, scipy, plotly, matplotlib)"

# Optional: rasterio for Engine 6
$rasterio_ok = & $py -c "import rasterio; print('ok')" 2>&1
if ($rasterio_ok -match "ok") {
    Write-Host "  OK  rasterio  (Engine 6: UCMuon Terrain enabled)"
} else {
    Write-Host "  --  rasterio not installed  (Engine 6 disabled)"
    Write-Host "      To enable: pip install rasterio"
}

# Engine 4 self-test
$e4_test = & $py -c @"
import importlib.util, sys, os
sys.path.insert(0, 'gui')
spec = importlib.util.spec_from_file_location('drv', 'gui/ucmuon_pumas_driver.py')
drv  = importlib.util.module_from_spec(spec); spec.loader.exec_module(drv)
d = float(drv._dedx(1000.0))
assert 1.70 < d < 1.95, f'bad dE/dx: {d}'
print(f'ok {d:.4f}')
"@ 2>&1
if ($e4_test -match "ok") {
    Write-Host "  OK  Engine 4 self-test passed"
} else {
    Write-Host "  WARN: Engine 4 self-test failed"
}

# ── 3. Fortran compiler (gfortran via MSYS2) ──────────────────────────────────
Write-Host "[3/4] Checking Fortran compiler (gfortran)..."

$gf_extra = @(
    "C:\msys64\ucrt64\bin\gfortran.exe",
    "C:\msys64\mingw64\bin\gfortran.exe",
    "C:\msys64\clang64\bin\gfortran.exe",
    "C:\msys64\usr\bin\gfortran.exe"
)
$gf = Find-Exe "gfortran" $gf_extra
$make_extra = @(
    "C:\msys64\ucrt64\bin\make.exe",
    "C:\msys64\mingw64\bin\make.exe",
    "C:\msys64\usr\bin\make.exe"
)
$mk = Find-Exe "make" $make_extra

$fortran_ok = $false

if ($gf -and $mk) {
    $gf_ver = & $gf --version 2>&1 | Select-Object -First 1
    Write-Host "  OK  $gf_ver"
    Write-Host "  OK  make: $mk"
    $fortran_ok = $true
} else {
    Write-Host ""
    Write-Host "  gfortran or make not found in PATH."
    Write-Host "  Engines 1 (MUSIC) and 2 (Bethe-Bloch) will be UNAVAILABLE."
    Write-Host "  Engines 4 and 5 (pure Python) will still work."
    Write-Host ""
    Write-Host "  To enable Fortran engines, install MSYS2:"
    Write-Host ""
    Write-Host "    Step 1 — Download MSYS2 from https://www.msys2.org/"
    Write-Host "             and run the installer (default path: C:\msys64)"
    Write-Host ""
    Write-Host "    Step 2 — Open 'MSYS2 UCRT64' from the Start menu and run:"
    Write-Host "             pacman -S mingw-w64-ucrt-x86_64-gcc-fortran make"
    Write-Host ""
    Write-Host "    Step 3 — Add C:\msys64\ucrt64\bin to your system PATH"
    Write-Host "             (Settings -> System -> Advanced -> Environment Variables)"
    Write-Host ""
    Write-Host "    Step 4 — Close and reopen PowerShell, then re-run install.ps1"
    Write-Host ""
}

# ── 4. Build Fortran binaries ─────────────────────────────────────────────────
Write-Host "[4/4] Building Fortran binaries..."

if ($fortran_ok) {
    # Pass the exact gfortran path to make in case it's not in PATH
    $env:FC = $gf
    & $mk local 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  OK  Fortran binaries built in bin\"
    } else {
        Write-Host "  WARN: make returned errors — some engines may be unavailable"
    }
} else {
    Write-Host "  SKIP (gfortran not found)"
}

# ── Summary ───────────────────────────────────────────────────────────────────
Write-Host ""
Write-Host "======================================================="
Write-Host "  ENGINE AVAILABILITY"
Write-Host "======================================================="
Write-Host ""

# Check for both .exe (MinGW) and no-extension (fallback)
function Test-Binary($stem) {
    return (Test-Path "bin\$stem.exe") -or (Test-Path "bin\$stem")
}

if (Test-Binary "ucmuon_transport_music_omp") {
    Write-Host "  [x] Engine 1  MUSIC stochastic MC        (bin\ucmuon_transport_music_omp)"
} else {
    Write-Host "  [ ] Engine 1  MUSIC                      (binary not built)"
    Write-Host "                                            music.f also required -- see docs\MUSIC_FILES.md"
}
if (Test-Binary "ucmuon_transport_bb_omp") {
    Write-Host "  [x] Engine 2  Bethe-Bloch + Highland MS  (bin\ucmuon_transport_bb_omp)"
} else {
    Write-Host "  [ ] Engine 2  Bethe-Bloch                (binary not built -- need gfortran)"
}
Write-Host "  [ ] Engine 3  PROPOSAL                   (not supported on Windows)"
Write-Host "  [x] Engine 4  UCMuon Stochastic          (Python -- always available)"
Write-Host "  [x] Engine 5  Backward MC                (Python -- always available)"
if ($rasterio_ok -match "ok") {
    Write-Host "  [x] Engine 6  UCMuon Terrain             (Python + rasterio)"
} else {
    Write-Host "  [ ] Engine 6  UCMuon Terrain             (pip install rasterio)"
}

Write-Host ""
Write-Host "======================================================="
Write-Host "  LAUNCH THE GUI"
Write-Host "======================================================="
Write-Host ""
Write-Host "    python -m streamlit run gui\ucmuon_gui.py"
Write-Host ""
Write-Host "  Or double-click run_gui.bat"
Write-Host ""
Write-Host "  Opens at http://localhost:8501 in your browser."
Write-Host ""
