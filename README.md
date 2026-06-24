# UCMuon — You See Muon!

**A simulation suite for cosmic-ray muon generation, transport, and muography**

UCLouvain Muography Group · Hamid Basiri · [hamid.basiri@uclouvain.be](mailto:hamid.basiri@uclouvain.be)
MIT License · [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20826984.svg)](https://doi.org/10.5281/zenodo.20826984)

---

## What it does

UCMuon simulates cosmic muon flux from the surface through rock, water, or ice, with applications to muography of geological structures, CO₂ storage monitoring, and glacier/bedrock imaging. It can be run interactively through a browser-based GUI on any laptop, or in batch mode on an HPC cluster using MPI-parallelised Fortran executables.

- **Surface muon generator** — eight spectrum models (CosmoALEPH, power-law, PARMA/EXPACS, Guan 2015, Frosin 2025, Bugaev/Gaisser 1998, Reyna-Bugaev 2006, cosmic electrons), three source geometries, OpenMP and MPI parallelism
- **Six transport engines** — **UCMuon-MC** (the native stochastic MC introduced by UCMuon: PDG-table-anchored per-process sampling + δ-ray straggling), MUSIC, Bethe-Bloch+MS, PROPOSAL, Backward MC, UCMuon Terrain (DEM ray-tracing)
- **Streamlit GUI** — interactive simulation, live progress, 3D plots, PHITS/Geant4 export, density analysis tab, works on Windows / macOS / Linux
- **HPC workflow** — MPI+OMP Fortran executables, SLURM scripts, automatic PHITS conversion
- **Transmission maps** — terrain engine outputs `T_sim = Φ_rock / Φ_sky` per direction; run at multiple densities to build a T_sim library for inversion
- **Density analysis** — pixel-wise density inversion by three methods: forward-model T_sim library fitting, plus analytical opacity inversion (column + mean density) from a single open-sky and target measurement; double-ratio maps, chi-squared landscape, uncertainty propagation (GUI tab 🔬 Density)
- **Engine comparison tool** — `tools/compare_engines.py` runs all engines on the same input and produces a physics audit, survival table, and diagnostic plot

---

## Status & scope (v0.9.0)

This is the first public release. The **core simulation pipeline is validated**
against independent codes (Geant4, PHITS, MUSIC, PROPOSAL); other components are
included but still under validation and are marked accordingly.

**Validated core — use with confidence:**
- Surface muon generator (eight spectrum models; source spectra cross-checked
  against Reyna 2006 and CosmoALEPH 2013 data)
- **UCMuon-MC** stochastic transport engine (agrees with MUSIC/PROPOSAL to ≲1 %
  on survival fractions in standard rock)
- Bethe-Bloch (CSDA) transport engine
- PROPOSAL transport engine (wrapper)
- Streamlit GUI for the above

**Experimental — included but still under validation (use with caution):**
- **PUMAS** transport engine
- **UCMuon Terrain** (DEM ray-tracing) and the transmission-map workflow
- **Density-inversion** suite (GUI 🔬 Density tab)
- HPC MPI build targets

**Not redistributed:** the **MUSIC** engine (V. A. Kudryavtsev) is not included
in this repository for licensing reasons. UCMuon ships only the wrapper; to use
the MUSIC engine, obtain the source from the author — see
[docs/MUSIC_FILES.md](docs/MUSIC_FILES.md). All other engines run without it.

**Licensing:** UCMuon's own code is MIT (see [LICENSE](LICENSE)). It bundles
third-party components with their own terms — notably **PARMA/EXPACS**
(`src/parma/`, `data/EXPACS/`), which is **non-commercial use only** (JAEA;
cite Sato 2015/2016), and **PUMAS** (LGPL-3.0). See
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) before any commercial use.

---

## Installation

> Full step-by-step instructions for every platform, including optional engines and troubleshooting, are in [INSTALL.md](INSTALL.md).

### Linux (Ubuntu / Debian / RHEL)

```bash
# 1. System dependencies
sudo apt install gfortran make          # Ubuntu/Debian
# sudo yum install gcc-gfortran make    # RHEL/Rocky

# 2. Clone and install
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon
bash setup.sh          # installs Python packages, builds Fortran binaries, checks PROPOSAL

# 3. Launch
bash run_gui.sh        # opens http://localhost:8501
```

`setup.sh` auto-detects available source files and only builds what is present. Engines 1 (UCMuon-MC) and 5 (pure Python) work immediately without any Fortran compiler.

---

### macOS

```bash
# 1. Install gfortran with OpenMP support
#    The Apple-bundled gcc is a Clang alias with no OpenMP — you need the real GCC.
brew install gcc                        # Homebrew (recommended)
# sudo port install gcc14               # MacPorts alternative

# 2. Clone and install
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon
bash setup.sh          # installs Python packages, builds Fortran binaries, sets up PROPOSAL venv

# 3. Launch
bash run_gui.sh        # opens http://localhost:8501
```

> **PROPOSAL note (macOS):** `setup.sh` automatically creates a system Python venv at `~/venvs/ucmuon` and installs PROPOSAL there. This is required because PROPOSAL crashes under Anaconda/miniforge due to a pybind11 ABI mismatch. `run_gui.sh` activates this venv automatically — no manual steps needed.

> **Compiler note (macOS):** If you have both Homebrew and MacPorts GCC, the Homebrew `mpif90` may wrap a different gfortran version than the one in your PATH, causing a `.mod` file mismatch. Set `OMPI_FC=/opt/homebrew/bin/gfortran-15` (or whichever version Homebrew installed) when building the MPI binaries: `OMPI_FC=/opt/homebrew/bin/gfortran-15 make hpc`.

---

### Windows

Engines 1 (**UCMuon-MC**) and 5 (Backward MC) are pure Python and work immediately. Engines 2 and 3 require a Fortran compiler via MSYS2. Engine 4 (PROPOSAL) is not supported on Windows.

```powershell
# 1. Install Python 3.11+ from https://www.python.org/downloads/
#    Tick "Add Python to PATH" during install.

# 2. Clone the repository (Git for Windows: https://git-scm.com/)
git clone https://github.com/hamidb90/UCMuon.git
cd UCMuon

# 3. Run the installer (installs Python packages; builds Fortran binaries if gfortran found)
powershell -ExecutionPolicy Bypass -File install.ps1

# 4. Launch the GUI
run_gui.bat            # double-click, or run from a terminal
```

**To enable Engines 2 and 3 (Fortran) on Windows:**

1. Install [MSYS2](https://www.msys2.org/) (accept default path `C:\msys64`)
2. Open **MSYS2 UCRT64** from the Start menu and run:
   ```bash
   pacman -S mingw-w64-ucrt-x86_64-gcc-fortran make
   ```
3. Add `C:\msys64\ucrt64\bin` to your Windows PATH (Settings → System → Advanced system settings → Environment Variables)
4. Re-run `install.ps1` — it will detect gfortran and build the binaries automatically

---

### HPC cluster (Lemaitre4 / CECI)

```bash
# One-time setup — add to ~/.bashrc:
echo "module load releases/2023b" >> ~/.bashrc
echo "module load foss/2023b"     >> ~/.bashrc
source ~/.bashrc

cd ~/UCMuon
make hpc            # builds ucmuon_gen, ucmuon_transport_music,
                    # ucmuon_transport_bb, ucmuon_to_phits

# Step 1 — generate surface muons:
sbatch hpc/run_ucmuon_gen.sh hpc/input_params.dat

# Step 2 — transport through rock:
sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_music.dat
```

See [`hpc/README_HPC.md`](hpc/README_HPC.md) for the full cluster workflow: MPI rank selection, MUSIC table initialisation, PHITS conversion, output formats, wall times, and troubleshooting.

---

## Platform compatibility

| Feature | Linux | macOS | Windows |
|---|:---:|:---:|:---:|
| Engines 1 (**UCMuon-MC**), 5 (pure Python) | ✓ | ✓ | ✓ |
| Engine 6 (UCMuon Terrain) | ✓ | ✓ | ✓ |
| Engine 7 (PUMAS, C binary) | ✓ | ✓ | ✓ (via MSYS2 gcc) |
| Engines 2, 3 (Fortran OMP) | ✓ | ✓ | ✓ (via MSYS2) |
| Engine 4 (PROPOSAL) | ✓ | ✓ | — |
| MPI+OMP HPC binaries | ✓ | ✓ | — |
| Streamlit GUI | ✓ | ✓ | ✓ |
| One-command install | `bash setup.sh` | `bash setup.sh` | `install.ps1` |

---

## Requirements

### Python ≥ 3.9

All core packages are installed by `setup.sh` / `install.ps1`:

```bash
pip install -r requirements.txt
# installs: streamlit, pandas, numpy, plotly, scipy, matplotlib
```

Optional:
- `pip install rasterio` — Engine 6 (UCMuon Terrain / DEM ray-tracing)
- `pip install proposal` — Engine 4 (use a system Python venv, not Anaconda — see PROPOSAL section below)

### Fortran compiler with OpenMP

Required only for Engines 2 (MUSIC) and 3 (Bethe-Bloch).

| Platform | Command |
|---|---|
| Ubuntu / Debian | `sudo apt install gfortran` |
| RHEL / Rocky | `sudo yum install gcc-gfortran` |
| macOS (Homebrew) | `brew install gcc` |
| macOS (MacPorts) | `sudo port install gcc14` |
| Windows | MSYS2: `pacman -S mingw-w64-ucrt-x86_64-gcc-fortran make` |
| HPC (Lemaitre4) | `module load releases/2023b && module load foss/2023b` |

---

## Directory structure

```
UCMuon/
├── README.md / LICENSE / CITATION.cff / ROADMAP.md / .gitignore
├── requirements.txt              Python package list
├── Makefile                      build all local OMP and HPC MPI targets
├── setup.sh / install.ps1        installers (auto-detect what to build)
├── run_gui.sh / run_gui.bat      GUI launchers (Linux/macOS / Windows)
├── tools/
│   └── compare_engines.py        multi-engine diagnostic: survival, energy, plots
│
├── src/
│   ├── common/                   shared Fortran: RANLUX, RNORML, CORGEN
│   ├── generator/                surface muon generation (8 spectra)
│   │   ├── ucmuon_source_module.f90   spectrum sampling (CosmoALEPH → Frosin)
│   │   ├── ucmuon_gen_omp.f90         OMP-only binary — GUI
│   │   └── ucmuon_gen.f90             MPI+OMP binary — HPC
│   ├── transport/
│   │   ├── music/                     MUSIC stochastic MC (Engine 2)
│   │   │   ├── ucmuon_transport_music_omp.f90
│   │   │   ├── music.f                Kudryavtsev (not redistributed — see docs/)
│   │   │   └── music-crosssections.f  Kudryavtsev (not redistributed)
│   │   └── bethe_bloch/               Bethe-Bloch+MS engine (Engine 3)
│   │       └── ucmuon_transport_bb_omp.f90
│   ├── parma/                    PARMA/EXPACS model (Sato 2015)
│   └── converters/               PHITS/Geant4 source file converters
│
├── gui/
│   ├── ucmuon_gui.py             Streamlit entry point  ← streamlit run this
│   ├── fast_flux_estimator.py    analytical flux / rate estimator
│   ├── ucmuon_stochastic_driver.py   Engine 1: UCMuon-MC (flagship)
│   ├── proposal_driver.py        Engine 4: PROPOSAL driver
│   ├── ucmuon_backward_mc.py     Engine 5: Backward MC integrator
│   ├── ucmuon_terrain_driver.py  Engine 6: DEM overburden + transmission map
│   ├── ucmuon_density_analysis.py    density inversion (library fit + opacity)
│   ├── ucmuon_pumas_driver.py    Engine 7: PUMAS backward/forward MC
│   ├── ucmuon_bb_driver.py       Engine 3: Bethe-Bloch Python driver
│   └── gui_*.py                  GUI panel modules (terrain, density, stochastic, …)
│
├── bin/                          compiled binaries (built by make; git-ignored)
├── data/                         physics tables (PARMA, MUSIC tables, …)
├── docs/                         MUSIC_FILES.md, ENGINE6_USAGE_GUIDE.md
├── hpc/                          SLURM scripts + annotated input templates
├── output/                       simulation outputs (regenerable; git-ignored)
│
├── examples/
│   ├── vesuvius/                 MURAVES-style Vesuvius worked example
│   └── terrain/                  terrain engine example script
│
├── benchmark/                    6-code validation (Geant4, PHITS, 4 UCMuon engines)
│   ├── codes/                    per-code source + README (geant4, phits, music, …)
│   ├── analysis/                 cross-code scripts (benchmark_analysis.py, …)
│   ├── results/                  summary CSVs + timing files
│   ├── figures/                  canonical comparison plots (v2_six_code/)
│   ├── reports/                  BENCHMARK_FEEDBACK.md (v2 reference), write-ups
│   └── geant4_muon_rock_v5/      git-ignored scratch tree (~8 GB raw .dat)
│
├── external/                     vendored third-party libraries
│   └── pumas-master/             PUMAS C library (used by Engine 7 / Makefile)
│
├── manuscript/                   paper (ucmuon_cpc_paper.tex/.pdf, scripts/, figs/)
├── misc/                         loose regenerable scratch (validate_spectra.py, …)
└── references/                   literature PDFs (kept private; stripped at public release)
```

---

## Building

```bash
bash setup.sh                     # recommended: checks all deps + builds local + sets up PROPOSAL venv

make local                        # OMP-only binaries (GUI use)
make hpc                          # MPI+OMP binaries (HPC use, needs mpif90)
make pumas                        # PUMAS C binary — local only (requires pumas-master/ in project root)
make clean                        # remove build/
make veryclean                    # clean + remove all binaries + pumas physics dumps
```

Object files and `.mod` files go into `build/` — the project root stays clean. Compiled binaries go into `bin/`.

**PUMAS note:** `make pumas` is skipped automatically if `external/pumas-master/` is absent — it does not affect the rest of the build. PUMAS is **local-only** (single-threaded C binary); it is not included in `make hpc` / HPC cluster targets. The first backward-MC run builds a physics dump at `bin/pumas_StandardRock.pumas` (~10 s); subsequent runs reload it in under a second.

---

## Transport engines

Seven engines are available in Tab 2 (🪨 **Transport**):

| # | Name | Physics | Requirements | When to use |
|---|---|---|---|---|
| ① | **UCMuon-MC** ★ | PDG per-process radiative MC + δ-ray straggling + Highland MS + decay | No install needed (numpy only) | **Native flagship engine** — default choice; cross-platform, no Fortran |
| ② | **MUSIC** | Full stochastic MC (Kudryavtsev 2009) | Fortran binary + MUSIC source files | External reference; deep overburden |
| ③ | **Bethe-Bloch** | Deterministic CSDA + Highland MS | Fortran binary (OMP) or pure Python/NumPy (`ucmuon_bb_driver.py`) | Fast cross-check; well above threshold |
| ④ | **PROPOSAL** | Full stochastic MC (Koehne 2013, Alameddine 2024) | `pip install proposal` + system Python venv | Independent stochastic cross-check |
| ⑤ | **Backward MC** | CSDA backward inversion + stochastic P_surv | No install needed | Instant sensitivity estimate, no input file |
| ⑥ | **UCMuon Terrain** | DEM ray-tracing + Backward MC physics | `pip install rasterio` | Real field sites (volcano, CCS, alpine) |
| ⑦ | **PUMAS** | True backward MC (Niess 2017); also forward mode | `make pumas` (C binary, local only) | 100% detector-efficiency flux spectra; no muon input file |

**Key distinction of Engine ⑥:** All other engines assume a flat uniform slab. UCMuon Terrain computes the actual rock overburden per azimuth and zenith direction from a real DEM file — essential when terrain height varies with direction.

**Key distinction of Engine ⑦ (PUMAS):** Runs in reverse: particles start at the detector with a sampled final-state energy and are propagated backward through rock to the surface. Every event reaches the generation surface (100% efficiency), versus ~0.1% survival for forward MC at depth. This gives an unbiased underground flux spectrum without needing a surface muon input file.

**Survival rate ordering** (100 m Standard Rock, 6 × 10⁵ muons, 5–300 GeV):

```
UCMuon-MC (32.3%)  <  PROPOSAL (32.8%)  ≈  MUSIC (32.9%)  <  BB (33.3%)
```

BB overestimates survival because it uses deterministic mean energy loss — no stochastic catastrophic radiative events. MUSIC and PROPOSAL agree within 0.1 pp. UCMuon-MC is ~0.6 pp below MUSIC after switching to the Bethe-Heitler hard-event spectrum.

| Overburden | Recommended engine |
|---|---|
| < 200 m.w.e. | **UCMuon-MC** — all engines agree within ~5% |
| 200–800 m.w.e. | **UCMuon-MC**, cross-checked against MUSIC or PROPOSAL |
| > 800 m.w.e. | MUSIC or PROPOSAL (Landau fluctuations dominate); validate UCMuon-MC |
| Real field deployment | UCMuon Terrain (UCMuon-MC physics) + MUSIC/PROPOSAL |

---

## Engine comparison tool

`tools/compare_engines.py` runs all available engines on the same input file with the same depth and material, then reports survival rates, energy loss, lateral displacement, and a muon-by-muon breakdown:

```bash
python3 tools/compare_engines.py output/muons_surface.dat 90.0 2.65
#                          ^input file               ^depth(m) ^density(g/cm³)
```

Output includes:
- Source-plane auto-detection (XY / XZ / YZ) with a warning for non-XY planes
- Per-engine survival rate table with delta vs reference
- Muon-by-muon agreement: "both survive", "both stop", "only engine A", etc.
- Underground energy difference for co-surviving muons
- Physics parameter audit explaining systematic differences
- Comparison plot saved to `output/engine_comparison.png`

The script automatically finds a working Python+PROPOSAL installation — if the default Python crashes on `import proposal` (common under Anaconda/miniforge), it searches for an alternative Python in known venv and conda locations.

---

## Full benchmark suite

The curated 6-code benchmark lives in `benchmark/`. It cross-checks the four UCMuon engines against Geant4 11.2 and PHITS 3.36 using an identical 6 × 10⁵ muon source population (six monoenergetic beams, 5–300 GeV) in Standard Rock.

```bash
# Run all 4 UCMuon engines at 5 depths (1, 25, 50, 100, 200 m) — takes ~40 min
cd benchmark/analysis
bash run_benchmark.sh

# Then regenerate the paper's survival table:
python3 manuscript/scripts/make_survival_table.py > manuscript/tab_survival_matrix.tex
```

See [`benchmark/README.md`](benchmark/README.md) for the full reproduction roadmap (source generation, per-code instructions, analysis steps).

**Measured survival fractions at d = 100 m** (N = 6 × 10⁵, single thread, Standard Rock, Apple Silicon Mac, 2026):

| Engine | Survival fraction | Time per 10⁵ | Bias vs Geant4 |
|---|:---:|:---:|:---:|
| Geant4 (FTFP_BERT) | 32.7 ± 0.1% | — | Reference |
| PHITS 3.36 | 32.6 ± 0.1% | — | −0.4% |
| MUSIC | 32.9 ± 0.1% | ~9 s | +0.6% |
| PROPOSAL | 32.8 ± 0.1% | ~98 s | +0.3% |
| Bethe-Bloch | 33.3 ± 0.1% | ~28 s | +1.9% |
| **UCMuon-MC** | 32.3 ± 0.1% | ~12 s | −1.3% |

All UCMuon engines agree with Geant4 within ±2% at all depths ≥ 25 m. The 10m depth is excluded from Geant4 comparison due to a CSDA boundary effect for the 5 GeV group. UCMuon-MC bias vs MUSIC was reduced from −1.1% to −0.6% by switching to the Bethe-Heitler (1−v)/v hard-event spectrum (PDG 2024 range table). BB timing improved to ~28 s after adding a CSDA range pre-filter.

---

## Source-plane convention

The generator supports three source planes: **XY** (z = const, depth in Z), **XZ** (y = const, depth in Y), and **YZ** (x = const, depth in X). All transport engines (Fortran and Python) automatically detect the source plane from the input file by finding the coordinate with near-zero variance, then use the correct momentum component for slant path computation.

This matters for non-standard geometries: a tunnel detector (XZ plane) has depth in the Y direction. Using the Z-component instead would give wrong slant paths for ~50% of muons, causing large spurious discrepancies between engines.

---

## Engine 1: UCMuon-MC (flagship)

The native, self-developed pure-Python stochastic transport engine — the default and primary engine of UCMuon, introducing a **table-anchored stochastic decomposition**: every stochastic term subtracts its own mean from the continuous part, so the mean dE/dx equals the evaluated PDG table *exactly at all energies, by construction* — there is no cross-section integration that can drift from the evaluated data.

**Physics (v2):**

- Loss decomposition from the **PDG 2024 per-process tables** (Groom–Mokhov–Striganov, 2024 electronic revision; embedded): ionisation, bremsstrahlung, pair production, photonuclear `[MeV cm²/g]` vs energy
- **Native per-material tables** for Standard Rock, water, ice, and iron — no rescaling approximation (the old rock-rescale underestimated iron radiative losses by ~38% at 1 TeV); seawater uses the water table with sub-percent composition rescaling, custom materials rescale the rock table
- Hard radiative events Poisson-sampled **per process** above `v_cut`, each with its own energy-dependent rate and spectrum: brems `(1−v)/v`, pair `1/v³`, photonuclear `1/v`
- **δ-ray straggling**: knock-on electrons with `T > 10 MeV` sampled from the Rutherford `1/T²` spectrum; the continuous ionisation term is restricted accordingly (ionisation straggling — the dominant fluctuation below ~100 GeV)
- Multiple scattering: Highland (1979) per step; muon decay: Poisson per step (`p·cτ`, `p = √(E²−m²)`)
- **Per-muon adaptive stepping**: every muon gets `dx ≈ 5 g/cm²` regardless of zenith angle
- Deterministic-bound pre-filter: only muons that cannot survive even with zero stochastic losses are killed instantly (strict upper bound on penetration — does not bias survival)
- Multiprocess-parallel: muons split across worker processes with independent RNG streams (`0 = auto`, near-linear speedup; reproducible for a fixed seed + worker count)

**Validation:** mean energy loss matches the PDG table integral to < 0.3% (40k-muon monoenergetic beams at 25/50 m); legacy v1 spectra (`1/v`, Bethe–Heitler single-shape) remain selectable for reproducibility.

---

## Engine 5: Backward MC Flux Integrator

Computes expected muon flux at a detector **without any surface muon file** in seconds. Uses backward CSDA energy inversion: for each `(E_det, θ)` bin, finds the required surface energy `E_s` and weights by the surface spectrum and Jacobian `|dE_s/dE_det|`.

---

## Engine 6: UCMuon Terrain — DEM-aware flux integrator

### Install

```bash
pip install rasterio
```

### Get a DEM file

Option A — Auto-download in the GUI (no account needed):

> Tab 2 → UCMuon Terrain → Section 1 → "Auto-download" tab → set bounding box → Download

Option B — Command line with `eio`:

```bash
pip install elevation
eio clip -o site_dem.tif --bounds LON_W LAT_S LON_E LAT_N
# Example — Puy de Dôme:
eio clip -o puydedome.tif --bounds 2.8 45.5 3.2 46.0
```

Option C — [OpenTopography](https://portal.opentopography.org): select area → SRTM GL1 30 m → Export GeoTIFF.

### Configure and run

In the GUI:

1. Upload or auto-download a GeoTIFF DEM (Section 1)
2. Enter detector GPS coordinates — latitude, longitude, altitude a.s.l. (Section 2)
3. Set rock density and spectrum model (Section 3)
4. Click **Quick terrain preview** for the overburden shape (~30 s)
5. Click **▶ Run UCMuon Terrain** for the full overburden + flux maps

### Output files

| File | Columns | Description |
|---|---|---|
| `terrain_overburden.dat` | azimuth, zenith, overburden [g/cm²], open_sky | per-direction rock thickness |
| `terrain_flux.dat` | azimuth, zenith, flux [m⁻² s⁻¹ sr⁻¹] | muon flux map |
| `terrain_transmission.dat` | azimuth, elevation, T_sim | transmission ratio Φ_rock/Φ_sky ∈ [0,1] in elevation-angle convention |
| `terrain_summary.dat` | scalars | total rate, max/median overburden, peak flux direction |

`terrain_transmission.dat` is the key output for density analysis. For the library-fit method, run the terrain engine at N densities (e.g. 1.5, 2.0, 2.65, 3.0 g/cm³) to build a T_sim library, then load it in the 🔬 Density tab. For the direct opacity inversion, a single run suffices — one `T_sim` map (or an overburden map) provides the line-of-sight path length, and the measured open-sky/target transmission is inverted analytically.

### Recommended DEM sources

| Source | Resolution | Coverage | Account |
|---|---|---|---|
| [OpenTopography SRTM GL1](https://opentopography.org) | 30 m | Global 60°S–60°N | No |
| [Copernicus COP30](https://spacedata.copernicus.eu/) | 30 m | Global | No |
| [EU-DEM v1.1](https://www.eea.europa.eu/data-and-maps/) | 25 m | Europe | No |
| [USGS EarthExplorer SRTM](https://earthexplorer.usgs.gov/) | 30 m | Global | Yes (free) |

---

## Density analysis workflow

The 🔬 **Density** tab in the GUI offers **three selectable inversion methods**. Method 1 fits a forward-model library; Methods 2 and 3 invert a single open-sky + target measurement analytically — no density library required.

### Method 1 — library matching (forward model)

**Step 1 — Build a T_sim library.**
Run the terrain engine at several bulk densities (e.g. 1.5, 2.0, 2.65, 3.0 g/cm³). Each run produces a `terrain_transmission.dat` file with columns `azimuth[deg]  elevation[deg]  T_sim`, where `T_sim = Φ_rock / Φ_sky ∈ [0, 1]`.

**Step 2 — Load the library in the Density tab.**
Paste the file paths into the T_sim Library section and click Load. The tab validates that all files share the same angular grid and shows a preview heatmap.

**Step 3 — Provide measured transmission T_data.**
Upload a file in the same format, or use the synthetic generator (Poisson noise at a chosen true density and number of events) for testing and validation.

**Step 4 — Run inversion.**
For each pixel `(az, el)`, the code finds `ρ̂` such that `T_sim(ρ̂) = T_data` by monotone interpolation over the T_sim library. Uncertainty `σ_ρ = σ_T / |dT_sim/dρ|` is propagated analytically. Status codes flag pixels with no sensitivity (T_sim ≈ 0 at all densities, i.e. too-thick rock) or out-of-range results.

**Step 5 — Inspect results.**
Four interactive maps are shown: reconstructed density (g/cm³), statistical uncertainty, double ratio `D = T_data / T_sim(ρ_ref)` (D > 1 → lower density than reference), and status map. A chi-squared landscape at any selected pixel is also available. Results are downloadable as `density_map.dat`.

Most physically complete (the library carries the full transport response), but it needs 3–5 terrain runs to build the lookup table.

### Method 2 — direct opacity inversion (single open-sky + target)

You do **not** need a density library. Transmission depends on the rock only through the **opacity** `X = ∫ρ·dl` (column density, g/cm²) and zenith angle, via the minimum energy needed to punch through (Groom CSDA range). Since `T(X, θ) = I(>E_min(X), θ) / I(>0, θ)` is monotonic in `X`, one measured transmission map inverts analytically to a unique opacity:

1. **Load one measured transmission map** `T = Φ_target / Φ_open` (columns `az, angle, T`; no density header needed).
2. **Pick a flux model** (any of the 5 sea-level models) and detector altitude.
3. **Choose a path-length source** for `L(az, el)` — flexible, three options:
   - **From a T_sim map** at known `ρ_sim` → `L = X_sim / (100·ρ_sim)` (no extra files);
   - **From a terrain overburden map** → `L = X / (100·ρ)` (exact DEM geometry);
   - **Opacity-only** — skip `L` and report column density `X` [g/cm²] directly.
4. **Run** → mean density `ρ̄ = X̂ / (100·L)`, with `σ_ρ̄ = σ_X / (100·L)` when counting statistics are supplied.

This is the classical muography density estimate (Tanaka 2007; Lesparre 2012): one open-sky reference + one target, no repeated simulations. Results download as `opacity_map.dat`.

### Method 3 — two-flux-map ratio

Identical physics to Method 2, but you supply the **open-sky** and **target** flux maps separately and the tab forms `T = Φ_target / Φ_open` for you — matching the data products of a real campaign where the incident flux is measured directly.

> **Note.** A single detector resolves only the *mean* density along each line of sight. Recovering the density distribution *along* the ray (3-D tomography) requires multiple detector viewpoints.

**Command-line API** (no GUI needed):
```python
from ucmuon_density_analysis import (
    # Method 1 — library matching
    load_transmission_map, build_tsim_library,
    invert_density_map, compute_double_ratio, write_density_map,
    # Methods 2 & 3 — direct opacity inversion
    load_generic_map, invert_opacity_map, opacity_to_density,
    path_length_from_tsim, load_overburden_as_L,
    transmission_from_flux_maps, write_opacity_map,
)

# --- Method 1: library matching ---
tsim_lib = build_tsim_library(["t150.dat", "t200.dat", "t265.dat", "t300.dat"])
az_c, el_c, T_data, meta = load_transmission_map("measured.dat")
rho_map, sigma_rho, status = invert_density_map(T_data, tsim_lib)
D_map = compute_double_ratio(T_data, tsim_lib[2.65])
write_density_map(az_c, el_c, rho_map, sigma_rho, status, "density_map.dat", meta)

# --- Method 2: direct opacity inversion (single open-sky + target) ---
az_c, el_c, T_data, _ = load_generic_map("measured_transmission.dat")
opacity, sigma_X, status = invert_opacity_map(T_data, el_c, model="reyna_bugaev")
# path length L(az, el) from one T_sim map at known density:
L_map = path_length_from_tsim(load_generic_map("tsim_2.65.dat")[2], 2.65, el_c)
rho_map, sigma_rho = opacity_to_density(opacity, L_map, sigma_X)
write_opacity_map(az_c, el_c, opacity, sigma_X, status, "opacity_map.dat",
                  rho_map=rho_map, sigma_rho=sigma_rho)

# --- Method 3: two-flux-map ratio ---
_, _, F_target, _ = load_generic_map("flux_target.dat")
_, _, F_open,   _ = load_generic_map("flux_open.dat")
T_data = transmission_from_flux_maps(F_target, F_open)   # then invert as in Method 2
```

---

## Engine 4: PROPOSAL

PROPOSAL is incompatible with Anaconda/miniforge due to a pybind11 ABI mismatch. Use a **system Python venv** (created automatically by `setup.sh` on Linux/macOS):

```bash
# Manual setup if needed:
/usr/bin/python3 -m venv ~/venvs/ucmuon
source ~/venvs/ucmuon/bin/activate
pip install -r requirements.txt
pip install proposal
```

`run_gui.sh` activates this venv automatically. The `tools/compare_engines.py` script auto-discovers a working Python+PROPOSAL installation — it tests each candidate in a subprocess so a segfaulting install cannot crash the comparison.

PROPOSAL is **not supported on Windows** (no official Windows wheels).

---

## Getting the MUSIC source files (Engine 2)

`music.f` and `music-crosssections.f` are not redistributed. Contact **Prof. Vitaly Kudryavtsev** — v.kudryavtsev@sheffield.ac.uk.
Reference: Kudryavtsev (2009), CPC 180, 339.
Place both files in `src/transport/music/` and run `bash setup.sh`.
See [`docs/MUSIC_FILES.md`](docs/MUSIC_FILES.md) for full instructions.

---

## GUI tabs

| Tab | Content |
|---|---|
| 🌌 **Generator** | Surface muon generation — spectrum, geometry, detector filter, live OMP progress |
| 🪨 **Transport** | Seven engines with engine-specific controls and live progress |
| 🗺️ **Terrain** | UCMuon Terrain: DEM upload, overburden maps, flux polar plots, transmission map output |
| 🔬 **Density** | Density inversion: (1) T_sim library fit; (2) direct opacity inversion; (3) two-flux-map ratio |
| 📊 **Results** | Energy/angle distributions, survival rate, 3D trajectories, rate estimation |
| ⚙️ **Config** | Session autosave/restore, quick stats, CSV export |

---

## HPC cluster workflow

See [`hpc/README_HPC.md`](hpc/README_HPC.md) for:

- One-time module setup and binary build
- MPI rank selection strategy (`nranks ≤ N_muons / 1000`)
- Step-by-step: generate → transport → PHITS conversion
- Output column formats (14-column generator, 18-column transport)
- PHITS `s-type=17` source section
- Expected wall times (e.g. 1M muons transport: ~8 min on 26 ranks)
- Troubleshooting table

---

## References

| Reference | Used for |
|---|---|
| Kudryavtsev (2009), CPC 180, 339 | MUSIC engine |
| Koehne et al. (2013), CPC 184, 2070 | PROPOSAL cross-sections |
| Alameddine et al. (2024), CPC 305, 109243 | PROPOSAL v7 |
| Groom, Mokhov & Striganov (2001), ADNDT 78, 183 | dE/dx tables — Engines 1, 5, 6 |
| Sato (2015), PLOS ONE 10(12) | PARMA/EXPACS spectrum |
| Guan et al. (2015), arXiv:1509.06176 | Guan 2015 spectrum |
| Frosin et al. (2025), J. Phys. G 52, 035002 | Frosin 2025 spectrum |
| Highland (1979), NIM 129, 497 | Multiple Coulomb scattering |
| Lüscher (1994), CPC 79, 100 | RANLUX RNG |
| NASA/JPL SRTM (2000) | DEM source for Engine 6 |

---

## Citation

If you use UCMuon, please cite it using the metadata in [`CITATION.cff`](CITATION.cff):

```bibtex
@software{ucmuon2026,
  author    = {Basiri, Hamid},
  title     = {{UCMuon}: Open-source cosmic muon simulation suite for muography},
  year      = {2026},
  publisher = {UCLouvain Muography Group},
  url       = {https://github.com/hamidb90/UCMuon},
  doi       = {10.5281/zenodo.XXXXXXX}
}
```

---

## Contact

Hamid Basiri — [hamid.basiri@uclouvain.be](mailto:hamid.basiri@uclouvain.be)
UCLouvain Muography Group, Louvain-la-Neuve, Belgium
