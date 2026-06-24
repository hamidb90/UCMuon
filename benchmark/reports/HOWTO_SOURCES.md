# MuonRock v6 — Switching Source Modes

Open `include/SimConfig.hh` and change **one line**:

```cpp
constexpr SourceMode SOURCE_MODE = SourceMode::FILE;      // Mode 1
constexpr SourceMode SOURCE_MODE = SourceMode::POWERLAW;  // Mode 2
constexpr SourceMode SOURCE_MODE = SourceMode::ECOMUG;    // Mode 3
```

Then recompile (`make`) and run the matching macro.

---

## Mode 1 — FILE

Reads every column from your text file. **z is used exactly as written**,
no override. All 8 columns must be present:

```
# PDG  x[mm]  y[mm]  z[mm]  px[MeV/c]  py[MeV/c]  pz[MeV/c]  Ekin[MeV]
 -13   -23295.8   11403.2   0.0   52.6   2949.2   -1083.8   3038.6
  13    17428.2  -24737.7   0.0  1062.5  -5036.7   -2094.4   5452.7
```

Configure column indices and units in `CSVColumnConfig` if your file
differs from the default layout.

Run: `./MuonRock run_file.mac`

---

## Mode 2 — POWERLAW

Built-in E^{-2.7} power-law spectrum between SPEC_EMIN_GEV and SPEC_EMAX_GEV.
Muons are shot **vertically downward** from a random XY position within
the slab. Good for quick testing and for comparing with MUSIC/PROPOSAL
under controlled conditions.

Key parameters in SimConfig.hh:
```cpp
constexpr double SPEC_EMIN_GEV = 100.0;   // 100 GeV minimum
constexpr double SPEC_EMAX_GEV = 10000.0; // 10 TeV maximum
constexpr double SPEC_GAMMA    = 2.7;      // spectral index
constexpr double POWERLAW_MU_MINUS_FRACTION = 0.43; // 43% mu-, 57% mu+
```

Run: `./MuonRock run_powerlaw.mac`

---

## Mode 3 — ECOMUG

Uses the EcoMug single-header library for a realistic cosmic-ray spectrum
including angular distribution. EcoMug generates from a hemispherical sky.

### Install EcoMug
```bash
git clone https://github.com/dr4kan/EcoMug.git
# No build needed — it is a single header: EcoMug/src/EcoMug.h
```

### Compile with EcoMug
```bash
mkdir build && cd build
cmake -DUSE_ECOMUG=ON -DECOMUG_INCLUDE=/path/to/EcoMug/src ..
make -j4
```

Key parameters in SimConfig.hh:
```cpp
constexpr double ECOMUG_SKY_RADIUS_MM = 25000.0; // 25 m sky dome radius
constexpr double ECOMUG_PMIN_MEV      = 100000.0; // 100 GeV/c
constexpr double ECOMUG_PMAX_MEV      = 10000000.0; // 10 TeV/c
constexpr double ECOMUG_MAX_THETA_RAD = 0.5236;  // 30 degrees
```

Run: `./MuonRock run_ecomug.mac`

---

## Summary

| Mode     | Macro              | Recompile? | Use case                        |
|----------|--------------------|------------|---------------------------------|
| FILE     | run_file.mac       | No         | Your own pre-generated source   |
| POWERLAW | run_powerlaw.mac   | No         | Quick tests, MUSIC benchmark    |
| ECOMUG   | run_ecomug.mac     | Only once  | Realistic angular+energy dist.  |
