# Geant4 — external reference code

Standalone Geant4 application (`MuonRock`) that propagates the benchmark muons
through a Standard-Rock slab and scores survival, energy, and scattering at each
depth plane. **Reference physics:** `FTFP_BERT` + `G4EmStandardPhysics_option4`,
Urban multiple-scattering model.

## Files here
| File / dir | Purpose |
|---|---|
| `MuonRock.cc` | main program |
| `src/`, `include/` | detector, primary generator, stepping/event/run actions |
| `include/EcoMug.h` | third-party EcoMug generator header (Pagano 2021, GPLv3) |
| `CMakeLists.txt` | build script |
| `run.mac`, `run_file.mac`, `run_powerlaw.mac`, `run_ecomug.mac`, `vis.mac` | macros |

## How to run the benchmark
```bash
# Build (needs a Geant4 11.x install with its environment sourced)
cmake -B build && cmake --build build -j

# File mode: read the shared benchmark source, 6e5 muons
./build/MuonRock -m file -f <benchmark_surface_as_geant4>.txt -q 600000
```
The Geant4 source file is produced from the common `benchmark_surface.dat` by
`../../codes/phits` conventions / the converter in the scratch tree
(`sources/convert_surface_to_geant4.py`). Output is a `*_muons.csv` per-event
file consumed by `../../analysis/benchmark_analysis.py` (`--geant4 <prefix>`).

> Per-event CSV and build artifacts are large and live only in the scratch tree
> `geant4_muon_rock_v5/` (git-ignored); only the code is kept here.
