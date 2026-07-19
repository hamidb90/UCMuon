# =============================================================================
#  install.ps1  --  UCMuon Windows installer
#  UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
#  MIT License 2026
#
#  Easiest: double-click install_windows.bat in the project root.
#  Or run from the project root in PowerShell:
#    powershell -ExecutionPolicy Bypass -File install.ps1
#  Scripted / CI use (no prompts, optional extras skipped):
#    powershell -ExecutionPolicy Bypass -File install.ps1 -NoPrompt
#
#  What this script does:
#    1. Checks Python >= 3.9
#    2. Installs all Python packages (streamlit, numpy, etc.)
#    3. Offers to install rasterio (Engine 6: UCMuon Terrain)
#    4. Looks for gfortran (via MSYS2 or PATH)
#    5. If found: builds the Fortran OMP binaries with make
#    6. If not found: offers to install MSYS2 + gfortran automatically
#       (winget), then builds; otherwise prints manual MSYS2 instructions
#    7. Prints an engine availability summary
#
#  Engines 1 (UCMuon-MC) and 5 (Backward MC) are pure Python
#  and work on Windows without any Fortran compiler.
#  Engines 2 (MUSIC) and 3 (Bethe-Bloch) need gfortran via MSYS2.
#  Engine 4 (PROPOSAL) is not supported on Windows.
# =============================================================================

param([switch]$NoPrompt)

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host ""
Write-Host "======================================================="
Write-Host "  UCMuon  --  Windows installer"
Write-Host "  UCLouvain Muography Group"
Write-Host "  Hamid Basiri  <hamid.basiri@uclouvain.be>"
Write-Host "  https://github.com/hamidb90/UCMuon"
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
# NOTE: with $ErrorActionPreference = "Stop", redirecting a native command's
# stderr (2>&1) turns any stderr line into a terminating NativeCommandError
# in Windows PowerShell 5.1. Every probe below that captures stderr must
# temporarily drop to "Continue".
$ErrorActionPreference = "Continue"
$ver_out = & $py --version 2>&1
$ErrorActionPreference = "Stop"
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
$ErrorActionPreference = "Continue"
$rasterio_ok = & $py -c "import rasterio; print('ok')" 2>&1
$ErrorActionPreference = "Stop"
if ($rasterio_ok -match "ok") {
    Write-Host "  OK  rasterio  (Engine 6: UCMuon Terrain enabled)"
} else {
    $do_install = $false
    if (-not $NoPrompt) {
        # Read-Host throws in a console-less host (e.g. CI); treat that as "skip"
        $ans = "n"
        try { $ans = Read-Host "  --  rasterio not found. Install it now to enable Engine 6 (UCMuon Terrain)? [Y/n]" } catch {}
        if ($ans -eq "" -or $ans -match "^[Yy]") { $do_install = $true }
    }
    if ($do_install) {
        Write-Host "      Installing rasterio..."
        # numpy<2.3 keeps pip from dragging numpy past what scipy pins allow
        & $py -m pip install -q rasterio "numpy<2.3"
        $ErrorActionPreference = "Continue"
        $rasterio_ok = & $py -c "import rasterio; print('ok')" 2>&1
        $ErrorActionPreference = "Stop"
        if ($rasterio_ok -match "ok") {
            Write-Host "  OK  rasterio installed  (Engine 6: UCMuon Terrain enabled)"
        } else {
            Write-Host "  WARN: rasterio install failed  (Engine 6 disabled)"
            Write-Host "      To enable manually: pip install rasterio"
        }
    } else {
        Write-Host "  --  rasterio not installed  (Engine 6 disabled)"
        Write-Host "      To enable: pip install rasterio"
    }
}

# Engine 1 (UCMuon-MC) self-test
$ErrorActionPreference = "Continue"
$e4_test = & $py -c @"
import importlib.util, sys, os
sys.path.insert(0, 'gui')
spec = importlib.util.spec_from_file_location('drv', 'gui/ucmuon_stochastic_driver.py')
drv  = importlib.util.module_from_spec(spec); spec.loader.exec_module(drv)
d = float(drv._dedx(1000.0))
assert 1.70 < d < 1.95, f'bad dE/dx: {d}'
print(f'ok {d:.4f}')
"@ 2>&1
$ErrorActionPreference = "Stop"
if ($e4_test -match "^ok") {
    Write-Host "  OK  Engine 1 (UCMuon-MC) self-test passed"
} else {
    Write-Host "  WARN: Engine 1 (UCMuon-MC) self-test failed:"
    $e4_test | ForEach-Object { Write-Host "      $_" }
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
    $ErrorActionPreference = "Continue"
    $gf_ver = & $gf --version 2>&1 | Select-Object -First 1
    $ErrorActionPreference = "Stop"
    Write-Host "  OK  $gf_ver"
    Write-Host "  OK  make: $mk"
    $fortran_ok = $true
} else {
    Write-Host ""
    Write-Host "  gfortran or make not found in PATH."
    Write-Host "  Without it, the muon GENERATOR (ucmuon_gen_omp) and Engines 2"
    Write-Host "  (MUSIC) and 3 (Bethe-Bloch) are unavailable."
    Write-Host ""

    # ── Offer automated MSYS2 + gfortran install via winget ──────────────────
    $msys_bash = "C:\msys64\usr\bin\bash.exe"
    $do_msys   = $false
    if (-not $NoPrompt) {
        if ((Test-Path $msys_bash) -or (Find-Exe "winget")) {
            $ans = "n"
            try { $ans = Read-Host "  Install MSYS2 + gfortran automatically now (~250 MB download, a few minutes)? [Y/n]" } catch {}
            if ($ans -eq "" -or $ans -match "^[Yy]") { $do_msys = $true }
        } else {
            Write-Host "  (winget not available -- cannot install MSYS2 automatically)"
        }
    }
    if ($do_msys) {
        if (-not (Test-Path $msys_bash)) {
            Write-Host "  Installing MSYS2 via winget..."
            $ErrorActionPreference = "Continue"
            & winget install --id MSYS2.MSYS2 -e --accept-source-agreements --accept-package-agreements
            $ErrorActionPreference = "Stop"
        } else {
            Write-Host "  OK  MSYS2 already present at C:\msys64"
        }
        if (Test-Path $msys_bash) {
            Write-Host "  Installing gcc-fortran and make (pacman)..."
            $ErrorActionPreference = "Continue"
            & $msys_bash -lc "pacman -Sy --noconfirm" 2>&1 | Out-Null
            & $msys_bash -lc "pacman -S --noconfirm --needed mingw-w64-ucrt-x86_64-gcc-fortran make" 2>&1 |
                ForEach-Object { Write-Host "      $_" }
            $ErrorActionPreference = "Stop"
            $gf = Find-Exe "gfortran" $gf_extra
            $mk = Find-Exe "make" $make_extra
            if ($gf -and $mk) {
                Write-Host "  OK  gfortran: $gf"
                Write-Host "  OK  make:     $mk"
                $fortran_ok = $true
            } else {
                Write-Host "  WARN: MSYS2 is installed but gfortran/make still not found."
            }
        } else {
            Write-Host "  WARN: MSYS2 installation did not complete."
        }
    }
}

if (-not $fortran_ok) {
    Write-Host ""
    Write-Host "  Pure-Python engines still work: Engine 5 (Backward MC) needs"
    Write-Host "  no generator; Engine 1 can transport existing muon files."
    Write-Host "  Muon generation also works without Fortran via the GUI's"
    Write-Host "  'Guaranteed-hit mode' workflow (pure Python)."
    Write-Host ""
    Write-Host "  To enable Fortran engines manually, install MSYS2:"
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

# ── 3b. Optional: PUMAS source for Engine 7 (downloaded, not bundled: LGPL-3.0) ──
$pumas_src = "external\pumas-master\src\pumas.c"
if ($fortran_ok -and (Test-Path $pumas_src)) {
    Write-Host "  OK  PUMAS source present  (Engine 7 will be built)"
} elseif ($fortran_ok) {
    $do_pumas = $false
    if (-not $NoPrompt) {
        $ans = "n"
        try { $ans = Read-Host "  --  Download PUMAS (~1 MB from github.com/niess/pumas) to enable Engine 7 (PUMAS backward MC)? [Y/n]" } catch {}
        if ($ans -eq "" -or $ans -match "^[Yy]") { $do_pumas = $true }
    }
    if ($do_pumas) {
        Write-Host "      Downloading PUMAS..."
        $ErrorActionPreference = "Continue"
        try {
            # PS 5.1 on older Windows may default to TLS 1.0; GitHub needs 1.2+
            [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
            $pumas_zip = Join-Path $env:TEMP "pumas-master.zip"
            Invoke-WebRequest -Uri "https://github.com/niess/pumas/archive/refs/heads/master.zip" `
                              -OutFile $pumas_zip -UseBasicParsing
            if (-not (Test-Path "external")) { New-Item -ItemType Directory -Path "external" | Out-Null }
            Expand-Archive -Path $pumas_zip -DestinationPath "external" -Force
            Remove-Item $pumas_zip -ErrorAction SilentlyContinue
        } catch {
            Write-Host "  WARN: PUMAS download failed: $($_.Exception.Message)"
        }
        $ErrorActionPreference = "Stop"
        if (Test-Path $pumas_src) {
            Write-Host "  OK  PUMAS source in external\pumas-master  (built in the next step)"
        } else {
            Write-Host "  WARN: PUMAS not available  (Engine 7 disabled)"
            Write-Host "      Manual: unzip the pumas repo to external\pumas-master (see docs\MUSIC_FILES.md)"
        }
    } else {
        Write-Host "  --  PUMAS skipped  (Engine 7 disabled; re-run install to enable later)"
    }
}

# ── 4. Build Fortran binaries ─────────────────────────────────────────────────
Write-Host "[4/4] Building Fortran binaries..."

if ($fortran_ok) {
    $ErrorActionPreference = "Continue"
    if (Test-Path "C:\msys64\usr\bin\bash.exe") {
        # Build inside an MSYS2 UCRT64 login shell: running make.exe directly
        # from PowerShell leaves its sub-shell without mkdir/cp/gfortran on PATH.
        $env:MSYSTEM        = "UCRT64"
        $env:CHERE_INVOKING = "1"    # keep the current directory in the login shell
        # -k: keep building other engines when one fails (e.g. MUSIC's
        #     ranmar_omp THREADPRIVATE common hits a gas .tls_common limit on PE)
        # -fno-common: works around that assembler limit with GCC on MinGW
        & "C:\msys64\usr\bin\bash.exe" -lc "make -k local FFLAGS_F77='-O2 -std=legacy -fno-common'" 2>&1 |
            ForEach-Object { Write-Host "  $_" }
    } else {
        # Non-MSYS2 toolchain already on PATH.
        # Forward slashes: backslashes would be eaten as escapes by make/sh.
        $env:FC = $gf -replace '\\', '/'
        & $mk local 2>&1
    }
    $ErrorActionPreference = "Stop"
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

if (Test-Binary "ucmuon_gen_omp") {
    Write-Host "  [x] Generator ucmuon_gen_omp            (bin\ucmuon_gen_omp)"
} else {
    Write-Host "  [ ] Generator ucmuon_gen_omp            (needs gfortran -- creates surface muons)"
}
Write-Host "  [x] Engine 1  UCMuon-MC (flagship)       (Python -- always available)"
if (Test-Binary "ucmuon_transport_music_omp") {
    Write-Host "  [x] Engine 2  MUSIC stochastic MC        (bin\ucmuon_transport_music_omp)"
} else {
    Write-Host "  [ ] Engine 2  MUSIC                      (binary not built)"
    Write-Host "                                            music.f also required -- see docs\MUSIC_FILES.md"
}
if (Test-Binary "ucmuon_transport_bb_omp") {
    Write-Host "  [x] Engine 3  Bethe-Bloch CSDA + MS      (bin\ucmuon_transport_bb_omp)"
} else {
    Write-Host "  [ ] Engine 3  Bethe-Bloch CSDA           (binary not built -- need gfortran)"
}
Write-Host "  [ ] Engine 4  PROPOSAL                   (not supported on Windows)"
Write-Host "  [x] Engine 5  Backward MC                (Python -- always available)"
if ($rasterio_ok -match "ok") {
    Write-Host "  [x] Engine 6  UCMuon Terrain             (Python + rasterio)"
} else {
    Write-Host "  [ ] Engine 6  UCMuon Terrain             (pip install rasterio)"
}
if (Test-Binary "ucmuon_transport_pumas") {
    Write-Host "  [x] Engine 7  PUMAS backward MC          (bin\ucmuon_transport_pumas)"
} else {
    Write-Host "  [ ] Engine 7  PUMAS backward MC          (re-run install and accept the PUMAS download)"
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
