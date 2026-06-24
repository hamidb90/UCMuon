# Geant4 vs PHITS Comparison

Focused guide for running and comparing **Geant4** and **PHITS** using the
common 100K muon source.

For the full four-code benchmark (+ MUSIC + PROPOSAL) see
**[COMPARISON_ALL_CODES.md](COMPARISON_ALL_CODES.md)**.

---

## Common Source File

`muons_geant4.txt` — 100 004 muons:

```
PDG   x[mm]   y[mm]   z[mm]   px[MeV/c]   py[MeV/c]   pz[MeV/c]   Ekin[MeV]
```

Energy range: 2.89 – 365.2 GeV · Mean: 5.40 GeV · µ⁻ fraction: 44.5 %

---

## Step 1 — Run Geant4

```bash
cd build
./MuonRock -m file -f ../muons_geant4.txt -q 100000
```

> Current run: **`outputs/run_file_20260516_233640`** — 100 000 events · 116 s · 860 evt/s

---

## Step 2 — Prepare PHITS Source

Convert `muons_geant4.txt` to PHITS ASCII dump format:

```bash
cd phits
python3 convert_source_to_phits.py ../muons_geant4.txt muons_for_phits.dat
```

`muon_rock.inp` is already set to `maxcas = 10000`, `maxbch = 10` → 100 000 total.
Run PHITS:

```bash
/Users/basiri/phits/bin/phits.sh muon_rock.inp
```

**Parse tally output** (required after every PHITS run, before running the benchmark):

```bash
python3 read_phits_output.py
```

This writes `phits_summary.csv` + `phits_timing.txt`. `n_source` is auto-detected.

---

## Step 3 — Run the Comparison

```bash
cd /Users/basiri/simulations/geant4/geant4_muon_rock_v5

python3 benchmark_analysis.py \
    --geant4 outputs/run_file_20260516_233640 \
    --phits  phits/ \
    --outdir figures_g4_vs_phits/
```

---

## Transmission Results

Numbers below use Geant4 100K (current source) and PHITS 10K (old source —
marked with *). Re-run PHITS for a proper apples-to-apples comparison.

| Depth | m.w.e. | Geant4 100K | PHITS 10K* |
|-------|--------|:-----------:|:----------:|
| 1 m   |   2.65 |  97.27 %    |  89.97 %   |
| 10 m  |  26.50 |  15.57 %    |  16.90 %   |
| 25 m  |  66.25 |   1.34 %    |   3.28 %   |
| 50 m  | 132.50 |   0.165 %   |   0.47 %   |
| 100 m | 265.00 |   0.017 %   |   0.02 %   |
| 200 m | 530.00 |   0.002 %   |   0.00 %   |

After PHITS is re-run with the same source, the agreement at 10–50 m depth
is expected to be within ~1–2 % (consistent with the previous 10K comparison).

---

## Output Figures (Geant4 + PHITS)

| Figure | What it shows | Both codes? |
|--------|---------------|-------------|
| `fig00_summary.png` | 2×3 overview panel | ✓ |
| `fig01_transmission.png` | Survival fraction vs depth | ✓ |
| `fig02_energy_loss.png` | Mean ΔE; process fractions; distributions | Geant4 only |
| `fig03_angular.png` | Scattering angle vs depth + Highland | Geant4 only |
| `fig04_lateral.png` | Lateral displacement vs depth | Geant4 only |
| `fig05_exit_spectrum.png` | Exit KE spectra at each depth | Geant4 only |
| `fig06_dedx.png` | dE/dx vs initial energy + Bethe-Bloch | Geant4 only |
| `fig07_charge_ratio.png` | µ⁻ fraction vs depth | Geant4 only |
| `fig08_secondaries.png` | Secondary particle types + KE spectrum | Geant4 only |
| `fig09_stopped.png` | Stopped muon KE and stopping depth | Geant4 only |
| `fig10_timing.png` | Wall time and throughput | ✓ |

PHITS provides only transmission and mean exit KE (from tallies, no per-event data).

---

## Known Limitations

| Issue | Impact |
|-------|--------|
| PHITS `phits_summary.csv` is from old 10K run / old source | Numbers not directly comparable until PHITS is re-run |
| Geant4 200 m: 2 surviving events | Statistically meaningless |
| PHITS tallies = integrated distributions | Only transmission and mean KE can be compared |
| PHITS `MeanZenithAngle_deg` ≠ Geant4 `AngleScat_Deg` | Different physical quantities — not comparable |

---

## File Locations

```
geant4_muon_rock_v5/
├── muons_geant4.txt                           ← common source (100 004 muons)
├── outputs/
│   └── run_file_20260516_233640_*.csv         ← Geant4 100K (use this)
├── phits/
│   ├── muon_rock.inp                          ← PHITS input (set maxcas=100000)
│   ├── muons_for_phits.dat                    ← PHITS source (re-generate from new muons_geant4.txt)
│   ├── phits_summary.csv                      ← PHITS tally output (old 10K — re-run)
│   ├── convert_source_to_phits.py             ← source converter
│   └── read_phits_output.py                   ← tally parser
├── benchmark_analysis.py                      ← comparison + figures
└── COMPARISON_ALL_CODES.md                    ← full 4-code guide
```
