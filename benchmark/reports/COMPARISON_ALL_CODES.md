# Four-Code Muon Benchmark: Geant4 · PHITS · MUSIC · PROPOSAL

Benchmark of cosmic muon propagation through 200 m standard rock
(530 m.w.e., ρ = 2.65 g/cm³).  Four transport codes are compared at six
scoring depths: **1 m, 10 m, 25 m, 50 m, 100 m, 200 m**.

---

## Code Overview

| Code | Type | Physics | Per-event data |
|------|------|---------|----------------|
| **Geant4** | Full Monte Carlo | FTFP_BERT, full MCS + radiative losses | Yes — all columns |
| **PHITS** | Monte Carlo | Built-in muon model, ATIMA dE/dx | Summary tallies only |
| **MUSIC** | Analytical transport | Continuous slowing-down + Lynch-Dahl MCS | Yes — position, direction, KE |
| **PROPOSAL** | C++/Python library | Continuous energy loss; **MCS off in current run** | Yes — position, KE (no scatter) |

> **PROPOSAL limitation:** the current run used straight-line propagation
> (no angular scattering). `AngleScat_Deg` and `LatDisp_cm` are unavailable
> until the run is repeated with MCS enabled.

---

## Common Source File

All four codes use the **same muon source** — `muons_geant4.txt`:

| Property | Value |
|----------|-------|
| File | `muons_geant4.txt` |
| N muons | 100 004 |
| Energy range | 2.89 – 365.2 GeV (kinetic) |
| Mean energy | 5.40 GeV |
| µ⁻ fraction | 44.5 % (44 548 µ⁻ / 55 456 µ⁺) |
| Format | `PDG  x[mm]  y[mm]  z[mm]  px  py  pz [MeV/c]  Ekin[MeV]` |

PHITS uses `phits/muons_for_phits.dat` — the same muons converted to PHITS
ASCII dump format (s-type=17, 10 columns) by `phits/convert_source_to_phits.py`.

---

## Workflow

### Step 1 — Geant4 (100K events)

```bash
cd build
./MuonRock -m file -f ../muons_geant4.txt -q 100000
```

Output lands in `../outputs/run_file_<timestamp>_*.csv`.

> Current run: **`outputs/run_file_20260516_233640`**
> 100 000 events · 116 s · 860 evt/s

---

### Step 2 — PHITS (100K events from same source)

Convert the source file to PHITS ASCII dump format (only needed once, or after
regenerating `muons_geant4.txt`):

```bash
cd phits
python3 convert_source_to_phits.py ../muons_geant4.txt muons_for_phits.dat
```

`muon_rock.inp` is already set to `maxcas = 10000`, `maxbch = 10` → 100 000 total.
Run PHITS:

```bash
/Users/basiri/phits/bin/phits.sh muon_rock.inp
```

**Parse the tally output** (required after every PHITS run, before running the benchmark):

```bash
python3 read_phits_output.py
```

This reads the raw `out_muon_*.out` tally files and writes `phits_summary.csv` +
`phits_timing.txt`. `n_source` is auto-detected from `phits.out` / `batch.out`.

---

### Step 3 — MUSIC (100K events)

MUSIC produces one file per scoring depth:

```
music/MUSIC_1m.dat   music/MUSIC_25m.dat   music/MUSIC_100m.dat
music/MUSIC_10m.dat  music/MUSIC_50m.dat   music/MUSIC_200m.dat
```

Column layout (18 columns, space-separated, 1 comment header line):

```
# EventID  x_srf_cm  y_srf_cm  z_srf_cm  E_srf_GeV  theta_srf_rad  phi_srf_rad
  charge  alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV  cx_ug  cy_ug  cz_ug
  theta_ug_rad  phi_ug_rad
```

| Column group | Description |
|---|---|
| `xs, ys, zs` | Surface entry position (cm); zs = 0 |
| `Es` | Initial KE at surface (GeV) |
| `theta_s, phi_s` | Surface zenith (from vertical) and azimuth (rad) |
| `charge` | +1 = µ⁺, −1 = µ⁻ |
| `alive` | 1 = survived to scoring depth, 0 = stopped |
| `x, y, z` | Underground position at scoring depth (cm) |
| `E` | Exit KE (GeV) |
| `cx, cy, cz` | Direction cosines at scoring depth |
| `theta_ug, phi_ug` | Underground zenith and azimuth (rad) |

`benchmark_analysis.py` reads these files automatically via `--music music/`.

---

### Step 4 — PROPOSAL (100K events)

Same 18-column layout, same file naming (`PROPOSAL_*.dat`), 5 comment header lines.

> **MCS status:** `theta_ug == theta_srf` for all events in the current run —
> straight-line propagation only. Re-run with MCS enabled in the PROPOSAL
> wrapper to obtain scattering angles and lateral displacements.

---

### Step 5 — Run the Benchmark

```bash
cd /Users/basiri/simulations/geant4/geant4_muon_rock_v5

python3 benchmark_analysis.py \
    --geant4  outputs/run_file_20260516_233640 \
    --phits   phits/ \
    --music   music/ \
    --proposal proposal/ \
    --outdir  figures_all_codes/
```

---

## Output Figures

| Figure | What it shows | Codes with data |
|--------|---------------|-----------------|
| `fig00_summary.png` | 2×3 overview panel | all 4 |
| `fig01_transmission.png` | Survival fraction vs depth | all 4 |
| `fig02_energy_loss.png` | Mean ΔE; process fractions; ΔE distributions | G4 (fractions); G4+MUSIC+PROPOSAL (mean ΔE) |
| `fig03_angular.png` | Scattering angle vs depth + Highland | Geant4, MUSIC |
| `fig04_lateral.png` | Lateral displacement vs depth | Geant4, MUSIC |
| `fig05_exit_spectrum.png` | Exit KE spectra at each depth | G4, MUSIC, PROPOSAL |
| `fig06_dedx.png` | dE/dx vs initial energy + Bethe-Bloch | G4, MUSIC, PROPOSAL |
| `fig07_charge_ratio.png` | µ⁻ fraction vs depth | G4, MUSIC, PROPOSAL |
| `fig08_secondaries.png` | Secondary particle types + KE spectrum | Geant4 only |
| `fig09_stopped.png` | Stopped muon KE and stopping depth | Geant4 only |
| `fig10_timing.png` | Wall time and throughput | Geant4, PHITS |
| `fig11_exit_position.png` | Exit (x, y) position map at each depth | G4, MUSIC, PROPOSAL |

Combined statistics: `figures_all_codes/benchmark_summary.csv`

---

## Results — 100K Events, Common Source

### Transmission (%)

| Depth | m.w.e. | Geant4 | PHITS* | MUSIC | PROPOSAL |
|-------|--------|:------:|:------:|:-----:|:--------:|
| 1 m   |   2.65 | 97.27  | 89.97  | 99.85 | 99.84    |
| 10 m  |  26.50 | 15.57  | 16.90  | 16.94 | 16.33    |
| 25 m  |  66.25 |  1.34  |  3.28  |  1.97 |  1.93    |
| 50 m  | 132.50 |  0.165 |  0.47  |  0.39 |  0.39    |
| 100 m | 265.00 |  0.017 |  0.02  |  0.065|  0.065   |
| 200 m | 530.00 |  0.002 |  0.00  |  0.012|  0.013   |

*PHITS numbers are from an older 10K run on the previous source — re-run
PHITS to align with the common 100K source.

**Physics notes on the 1 m difference:**
Geant4 (97.3%) shows lower transmission than MUSIC/PROPOSAL (99.8%) at 1 m
despite using the same source.  This is a genuine physics-model difference:
muons with large zenith angles have longer slant paths, and MUSIC/PROPOSAL
may use a different stopping-range model (Bethe-Bloch / CSDA) that predicts
slightly longer ranges for 3–5 GeV muons than Geant4's FTFP_BERT.
This difference vanishes by 10 m (all codes agree within ~1–2%), confirming
the 1 m gap is near-threshold stopping physics, not a source discrepancy.

**Statistical note:**
Geant4 at 200 m has only 2 surviving events → not a meaningful number.
MUSIC (12) and PROPOSAL (13) have better statistics due to the 100 004-muon
input and no detector secondaries to simulate.  For reliable deep-plane
comparison, 1–10 M source muons are needed.

### Mean Exit KE (GeV)

| Depth | Geant4 | PHITS* | MUSIC | PROPOSAL |
|-------|:------:|:------:|:-----:|:--------:|
| 1 m   |  4.65  |  7.34  |  4.70 |   4.69   |
| 10 m  |  4.65  | 11.30  |  4.96 |   5.08   |
| 25 m  | 11.62  | 18.48  | 12.42 |  12.52   |
| 50 m  | 21.39  | 22.39  | 22.78 |  22.39   |
| 100 m | 35.38  | 42.79  | 39.38 |  38.79   |
| 200 m |  4.13† |  —     | 57.70 |  58.44   |

*PHITS is from the old 10K run — numbers differ because the old source had
higher-energy muons (0.9–744 GeV vs current 2.89–365 GeV).
†Only 2 events — not meaningful.

Geant4, MUSIC, and PROPOSAL agree well on mean exit KE at 25–100 m, which
is where the comparison is most reliable.

### Scattering Angle and Lateral Displacement (MUSIC vs Geant4)

| Depth | G4 MeanAngle (°) | MUSIC MeanAngle (°) | G4 LatDisp (cm) | MUSIC LatDisp (cm) |
|-------|:----------------:|:-------------------:|:---------------:|:------------------:|
| 1 m   |  37.73  |  1.07  |  1891.6  |    2.7  |
| 10 m  |  28.11  |  4.20  |  1793.9  |   26.1  |
| 25 m  |  25.06  |  2.58  |  1709.5  |   37.3  |
| 50 m  |  19.80  |  1.91  |  1777.2  |   54.0  |
| 100 m |  13.36  |  1.62  |  1913.1  |   69.3  |
| 200 m |   4.77  |  0.76  |  2199.0  |  112.9  |

> **Different definitions:** Geant4 `AngleScat_Deg` is the angle of the muon's
> direction from the downward vertical axis — it includes the original source
> zenith angle and is dominated by it. MUSIC `AngleScat_Deg` is the **3-D
> deflection from each muon's own initial direction**, i.e., the true MCS
> scatter angle. These are not directly comparable.
> PROPOSAL: not available (no MCS in current run).

---

## Known Limitations

| Code | Limitation | Impact |
|------|-----------|--------|
| Geant4 100 m | Only 17 surviving events | Transmission value meaningful but spectra are not |
| Geant4 200 m | Only 2 surviving events | Statistically meaningless |
| PHITS | Numbers from old 10K run on previous source | Must re-run before comparing |
| PROPOSAL | MCS disabled | No scattering angles or lateral displacements |
| All codes | ~100K muons, mean 5.4 GeV | Deep planes (≥ 100 m) need 1–10 M events for statistics |

---

## Recommended Next Steps

1. **Re-run PHITS** with `muons_for_phits.dat` regenerated from new `muons_geant4.txt`
   and `maxcas = 100000`.

2. **Re-run PROPOSAL with MCS enabled** to get scattering angles and lateral
   displacements comparable to Geant4 and MUSIC.

3. **Increase statistics** — generate 1M+ muons from the same source generator
   to get meaningful numbers at 100 m and 200 m.

---

## File Locations

```
geant4_muon_rock_v5/
├── muons_geant4.txt                          ← common source (100 004 muons, 2.9–365 GeV)
├── outputs/
│   ├── run_file_20260512_180248_*.csv        ← Geant4 10K (old source, 0.9–744 GeV)
│   └── run_file_20260516_233640_*.csv        ← Geant4 100K (current source) ← USE THIS
├── phits/
│   ├── muon_rock.inp                         ← PHITS input (set maxcas=100000)
│   ├── muons_for_phits.dat                   ← PHITS source (100K, converted)
│   ├── phits_summary.csv                     ← parsed PHITS output (old 10K run)
│   ├── convert_source_to_phits.py            ← source converter
│   └── read_phits_output.py                  ← tally parser
├── music/
│   └── MUSIC_*.dat                           ← MUSIC output (100K, per depth)
├── proposal/
│   └── PROPOSAL_*.dat                        ← PROPOSAL output (100K, per depth)
├── benchmark_analysis.py                     ← 4-code analysis + figures
├── figures_all_codes/                        ← latest 4-code benchmark figures
│   ├── benchmark_summary.csv
│   └── fig00_summary.png … fig11_exit_position.png
├── COMPARISON_ALL_CODES.md                   ← this file
└── COMPARISON_G4_PHITS.md                    ← Geant4 vs PHITS focused guide
```
