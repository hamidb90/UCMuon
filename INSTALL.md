# UCMuon — Installation Guide

**UCLouvain Muography Group** · Hamid Basiri · [hamid.basiri@uclouvain.be](mailto:hamid.basiri@uclouvain.be) · MIT License

This guide covers installation on **Linux**, **macOS**, and **Windows**, plus the optional HPC (MPI) build. For a project overview see [README.md](README.md).

---

## Contents

1. [Requirements at a glance](#1-requirements-at-a-glance)
2. [Linux](#2-linux)
3. [macOS](#3-macos)
4. [Windows](#4-windows)
5. [HPC cluster (MPI+OMP)](#5-hpc-cluster-mpiomp)
6. [Optional components](#6-optional-components)
7. [Verifying the installation](#7-verifying-the-installation)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Requirements at a glance

| Component | Needed for | Linux | macOS | Windows |
|---|---|:---:|:---:|:---:|
| **Python ≥ 3.9** + pip | GUI and Engines ① ⑤ ⑥ | required | required | required |
| **gfortran with OpenMP** + make | Engines ② (MUSIC), ③ (Bethe-Bloch Fortran) | `apt install gfortran` | `brew install gcc` | MSYS2 |
| **gcc (C compiler)** | Engine ⑦ (PUMAS, optional) | ✓ | ✓ | MSYS2 |
| **rasterio** (pip, optional) | Engine ⑥ (UCMuon Terrain / DEM) | ✓ | ✓ | ✓ |
| **PROPOSAL** (pip, optional) | Engine ④ | ✓ | ✓ (system venv) | not supported |
| **mpif90 (OpenMPI/MPICH)** | HPC binaries only | ✓ | ✓ | — |

The flagship engine **① UCMuon-MC** and Engine **⑤ Backward MC** are pure Python — they work on every platform with no compiler at all. The installers (`setup.sh` / `install.ps1`) auto-detect what is available and only build what they can; missing optional components simply disable the corresponding engine.

Core Python packages (installed automatically): `streamlit`, `numpy`, `pandas`, `scipy`, `plotly`, `matplotlib`.

---

## 2. Linux

Tested on Ubuntu/Debian and RHEL/Rocky.

### 2.1 System dependencies

```bash
# Ubuntu / Debian
sudo apt update
sudo apt install gfortran make git python3 python3-pip python3-venv

# RHEL / Rocky / Fedora
sudo yum install gcc-gfortran make git python3 python3-pip
```

### 2.2 Clone and install

```bash
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon
bash setup.sh
```

`setup.sh` will:

1. Check `gfortran` and OpenMP support
2. Install the core Python packages (`pip install -r requirements.txt`)
3. Create a venv at `~/venvs/ucmuon` and install **PROPOSAL** (Engine ④) into it
4. Check the Fortran source inventory and build the local OMP binaries (`make local`)
5. Print an engine availability summary

Useful flags:

```bash
bash setup.sh --no-python              # skip all Python steps (build only)
bash setup.sh --no-proposal            # skip the PROPOSAL venv
bash setup.sh --python=python3.11      # use a specific interpreter
```

### 2.3 Launch

```bash
bash run_gui.sh        # opens http://localhost:8501 in your browser
bash run_gui.sh --threads=4    # optionally cap OMP threads
```

Always use `run_gui.sh` rather than calling `streamlit` directly — it activates the PROPOSAL venv and sets `OMP_NUM_THREADS` automatically.

---

## 3. macOS

Works on both Apple Silicon and Intel Macs.

### 3.1 System dependencies

The Apple-bundled `gcc` is a Clang alias **without OpenMP** — you need real GCC for the Fortran engines:

```bash
# Xcode command line tools (provides /usr/bin/python3, git, etc.)
xcode-select --install

# Real GCC with gfortran + OpenMP
brew install gcc                 # Homebrew (recommended)
# sudo port install gcc14        # MacPorts alternative
```

### 3.2 Clone and install

```bash
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon
bash setup.sh
```

Same steps and flags as on Linux (see §2.2).

> **PROPOSAL note (macOS):** PROPOSAL segfaults under Anaconda/miniforge due to a pybind11 ABI mismatch. `setup.sh` therefore creates a separate venv at `~/venvs/ucmuon` based on the **system Python** (`/usr/bin/python3`) and installs PROPOSAL there. `run_gui.sh` activates this venv automatically — no manual steps needed.

> **Compiler note:** If you have both Homebrew and MacPorts GCC installed, the Homebrew `mpif90` may wrap a different gfortran than the one in your PATH, causing a `.mod` file version mismatch when building the HPC binaries. Fix with:
> ```bash
> OMPI_FC=/opt/homebrew/bin/gfortran-15 make hpc   # adjust to your gfortran version
> ```

### 3.3 Launch

```bash
bash run_gui.sh        # opens http://localhost:8501
```

---

## 4. Windows

Engines ① (**UCMuon-MC**, flagship) and ⑤ (Backward MC) are pure Python and work immediately. Engines ② and ③ need gfortran via MSYS2. Engine ④ (PROPOSAL) is **not supported on Windows**.

### 4.1 Python

Install **Python 3.11+** from <https://www.python.org/downloads/> and **tick "Add Python to PATH"** during install. Or via winget:

```powershell
winget install Python.Python.3.11
```

### 4.2 Clone and install

Install [Git for Windows](https://git-scm.com/) if needed, then in **PowerShell**:

```powershell
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon
powershell -ExecutionPolicy Bypass -File install.ps1
```

Or, if you downloaded UCMuon as a zip (e.g. from Zenodo), simply **double-click
`install_windows.bat`** in the extracted folder — no terminal needed.

`install.ps1` will:

1. Check Python ≥ 3.9
2. Install all Python packages
3. Look for gfortran/make (in PATH or a standard MSYS2 location)
4. Build the Fortran OMP binaries if found — otherwise skip them and print instructions
5. Print an engine availability summary

### 4.3 Launch

Double-click **`run_gui.bat`**, or from a terminal:

```powershell
run_gui.bat
# equivalent to:  python -m streamlit run gui\ucmuon_gui.py
```

The GUI opens at <http://localhost:8501>.

### 4.4 Enabling the Fortran engines (② MUSIC, ③ Bethe-Bloch)

1. Install [MSYS2](https://www.msys2.org/) (accept the default path `C:\msys64`)
2. Open **MSYS2 UCRT64** from the Start menu and run:
   ```bash
   pacman -S mingw-w64-ucrt-x86_64-gcc-fortran make
   ```
3. Add `C:\msys64\ucrt64\bin` to your Windows PATH
   (Settings → System → About → Advanced system settings → Environment Variables)
4. Close and reopen PowerShell, then re-run `install.ps1` — it detects gfortran and builds the binaries automatically

---

## 5. HPC cluster (MPI+OMP)

The MPI binaries (`ucmuon_gen`, `ucmuon_transport_music`, `ucmuon_transport_bb`, `ucmuon_to_phits`) are for batch production runs on a cluster. On Lemaitre4 / CECI:

```bash
# One-time setup — add to ~/.bashrc:
echo "module load releases/2023b" >> ~/.bashrc
echo "module load foss/2023b"     >> ~/.bashrc
source ~/.bashrc

cd ~/UCMuon
make hpc

# Step 1 — generate surface muons:
sbatch hpc/run_ucmuon_gen.sh hpc/input_params.dat

# Step 2 — transport through rock:
sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_music.dat
```

On other clusters, any environment providing `gfortran` and `mpif90` (e.g. a `foss` toolchain) works. See [`hpc/README_HPC.md`](hpc/README_HPC.md) for the full workflow: MPI rank selection, MUSIC table initialisation, PHITS conversion, output formats, wall times, and troubleshooting.

---

## 6. Optional components

### Engine ② — MUSIC source files

`music.f` and `music-crosssections.f` (Kudryavtsev 2009) are **not redistributed** with UCMuon. If they are absent, the MUSIC targets are skipped automatically and all other engines remain fully functional. See [`docs/MUSIC_FILES.md`](docs/MUSIC_FILES.md) for how to obtain and place them, then re-run `bash setup.sh` (or `make local`).

### Engine ④ — PROPOSAL (Linux / macOS only)

Installed automatically by `setup.sh`. To do it manually:

```bash
/usr/bin/python3 -m venv ~/venvs/ucmuon      # system Python, NOT Anaconda
source ~/venvs/ucmuon/bin/activate
pip install -r requirements.txt
pip install proposal
```

### Engine ⑥ — UCMuon Terrain (DEM ray-tracing)

```bash
pip install rasterio
```

### Engine ⑦ — PUMAS (C binary)

Requires the PUMAS sources (LGPL, Niess et al. 2017) in `external/pumas-master/`:

```bash
git clone https://github.com/niess/pumas.git external/pumas-master
make pumas
```

If `external/pumas-master/` is absent the target is skipped with a note. The first backward-MC run builds a physics dump at `bin/pumas_StandardRock.pumas` (~10 s); subsequent runs reload it in under a second. PUMAS is local-only (not part of `make hpc`).

---

## 7. Verifying the installation

Both installers end with an **ENGINE AVAILABILITY** summary, e.g.:

```
  [x] Engine 1  UCMuon-MC (flagship)       (Python)
  [x] Engine 2  MUSIC stochastic MC        (bin/ucmuon_transport_music_omp)
  [x] Engine 3  Bethe-Bloch + Highland MS  (bin/ucmuon_transport_bb_omp)
  [x] Engine 4  PROPOSAL stochastic MC     (system Python venv)
  [x] Engine 5  Backward MC                (Python)
  [x] Engine 6  UCMuon Terrain             (Python + rasterio)
```

Re-run `bash setup.sh` (or `install.ps1`) at any time to re-check — both are idempotent.

Then launch the GUI and confirm it opens at <http://localhost:8501>:

```bash
bash run_gui.sh        # Linux / macOS
run_gui.bat            # Windows
```

A minimum working install only needs Engine ① checked — that is the flagship engine and the GUI default.

---

## 8. Troubleshooting

**`gfortran: command not found`**
Install it per the table in §1 and re-run the installer. On macOS make sure it is Homebrew/MacPorts GCC, not the Apple Clang alias.

**macOS: `gfortran does not support -fopenmp`**
You are picking up the Apple-bundled compiler. Run `brew install gcc` and make sure `/opt/homebrew/bin` (Apple Silicon) or `/usr/local/bin` (Intel) precedes `/usr/bin` in your PATH.

**PROPOSAL import segfaults / crashes**
You are likely under Anaconda/miniforge. PROPOSAL must live in a system-Python venv (see §6). `run_gui.sh` activates `~/venvs/ucmuon` automatically; if the venv is missing, re-run `bash setup.sh`.

**`make hpc` fails with "Fatal Error: Cannot read module file ..."**
`mpif90` wraps a different gfortran than the one that compiled the `.mod` files. Run `make clean`, then rebuild with `OMPI_FC=<path-to-matching-gfortran> make hpc` (see §3.2).

**Windows: `install.ps1` cannot be run ("running scripts is disabled")**
Run it exactly as shown — `powershell -ExecutionPolicy Bypass -File install.ps1` — which bypasses the policy for this one invocation only.

**Windows: GUI fails with "streamlit is not recognized"**
Python is not on PATH or packages were installed into a different interpreter. Re-install Python with "Add Python to PATH" ticked, reopen the terminal, and re-run `install.ps1`.

**Engine 2 (MUSIC) shows as unavailable**
The MUSIC source files are not bundled — see §6 and [`docs/MUSIC_FILES.md`](docs/MUSIC_FILES.md). All other engines work without it.

**GUI is slow / uses one core**
Set the thread count explicitly: `bash run_gui.sh --threads=8`, or set `OMP_NUM_THREADS` before launching.
