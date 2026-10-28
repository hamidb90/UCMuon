# Geant4 — external reference code

Standalone Geant4 application (`MuonRock`) that propagates the benchmark muons
through a Standard-Rock slab and scores survival, energy, and scattering at each
depth plane. **Reference physics:** `FTFP_BERT` with its default electromagnetic
physics (`MuonRock.cc` registers `new FTFP_BERT()` and replaces nothing; up to
1.3.2 this file said `G4EmStandardPhysics_option4` with the Urban model).

**Rock.** Standard Rock as PDG defines it and as the UCMuon engines use it: one
nucleus with Z = 11, A = 22 (Z/A = 0.5), I = 136.4 eV, and the Sternheimer
density-effect parameters of Groom, Mokhov and Striganov (2001), Table IV-6
(`include/SimConfig.hh`). `DetectorConstruction.cc` builds it from a single
A = 22 isotope, because `G4Material(name, Z, A, density)` in Geant4 11.4 ignores
A and uses the natural element, here sodium: A = 22.99 (Z/A = 0.478),
I = 149 eV and generic density-effect parameters, 4 to 5 % less ionisation loss
than Standard Rock. That was the benchmark material up to UCMuon 1.4.1
(reference runs `run_file_20260522_210114`, six beams, and
`run_file_20260614_095725`, realistic source); the run log now prints
`Standard rock: Z/A = 0.5, I = 136.4 eV`.

## Files here
| File / dir | Purpose |
|---|---|
| `MuonRock.cc` | main program |
| `src/`, `include/` | detector, primary generator, stepping/event/run actions |
| `include/EcoMug_NOTICE.txt` | EcoMug (Pagano 2021, GPL-3.0) is not included; how to get it for `-m ecomug` |
| `run_parallel.py` | runs a FILE-mode source as independent chunks on several cores and merges them into one run |
| `CMakeLists.txt` | build script |
| `run.mac`, `run_file.mac`, `run_powerlaw.mac`, `run_ecomug.mac`, `vis.mac` | macros |

## How to run the benchmark
```bash
# Build (needs a Geant4 11.x install with its environment sourced)
cmake -B build && cmake --build build -j

# File mode: read the shared benchmark source, 6e5 muons, one core
./build/MuonRock -m file -f <benchmark_surface_as_geant4>.txt -s 20261019 -q 600000

# The same on 8 cores: 32 chunks, seeds 20261019 + chunk, merged with the event
# numbering and files of a single run (the reference of 2026-10-19)
python3 run_parallel.py --exe build/MuonRock --source <source>.txt \
    --out outputs --name run_stdrock_xy1km_20261022_mono --chunks 32 --workers 8 \
    --seed 20261019 --xy 100000

# The slab must be much wider than the source. The vertical beams start
# anywhere in +-25 m, the default -xy: in run_stdrock_20261019_mono, with that
# default, muons starting near the sides scattered out of the slab and were
# lost, up to 1.1 % of a beam at 100-200 m (100 GeV at 100 m: 96.52 % survival
# against 97.60 % with +-1 km). The realistic source starts up to 140 m from
# the axis.
python3 run_parallel.py --exe build/MuonRock --source source_realistic_g4.txt \
    --out outputs --name run_stdrock_20261019_realistic --chunks 32 --workers 8 \
    --seed 20261119 --xy 100000

# Inclined beams (the six beams tilted to 60, 80 and 85 deg by
# ../../analysis/make_inclined_inputs.py): the scoring plane of slant length L
# is placed at the depth L cos(theta) and keeps its name, so every depth in the
# output is the slant length L; the slab ends 5 cm below the deepest plane
for TH in 85 80 60; do
  python3 run_parallel.py --exe build/MuonRock --source benchmark_geant4_th${TH}.txt \
      --out outputs --name run_stdrock_xy1km_20261027_incl${TH} --chunks 32 --workers 8 \
      --seed 20261019 --xy 100000 --incl ${TH}
done
```
Options added after 1.4.1: `-s <seed>` (Geant4's default seed otherwise; recorded
in the run summary), `-xy <cm>` (slab half-width) and `-incl <deg>` (planes at
L cos(incl); without it the output is byte-identical to the build without the
option). The merged summary's
`Elapsed` is the CPU time summed over the chunks; its wall time is stated
separately.
The Geant4 source file is produced from the common `benchmark_surface.dat` by
`src/converters/ucmuon_to_geant4.py --mode gen`, which reproduces the reference
runs' source to the last printed digit of the kinetic energy (10^-9 relative;
the scratch tree's `sources/convert_surface_to_geant4.py` wrote it). Output is a `*_muons.csv` per-event
file consumed by `../../analysis/benchmark_analysis.py` (`--geant4 <prefix>`).

> Per-event CSV and build artifacts are large and live only in the scratch tree
> `geant4_muon_rock_v5/` (git-ignored); only the code is kept here.
