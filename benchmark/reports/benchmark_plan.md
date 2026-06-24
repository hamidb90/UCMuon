# UCMuon Engine Benchmark Plan
**UCLouvain Muography Group | Hamid Basiri**  
Version 1.1 — Standard Rock Flat-Slab Benchmark

---

## Overview

Systematic, reproducible benchmark comparing the three primary UCMuon transport engines
against each other and against an analytical CSDA reference. Two phases:

- **Phase 1 (required):** Vertical muons (θ=0°), Standard Rock, five overburden depths.
- **Phase 2 (required):** Realistic cos²θ distribution, same material. Tests that the
  `depth / cos(θ)` path-length correction is applied consistently across all engines.

**Engines compared:**

| Label | Binary / Module | Physics |
|---|---|---|
| MUSIC | `bin/ucmuon_transport_music_omp` | Table-interpolated dE/dx + stochastic fluctuations (Kudryavtsev 2009) |
| BB | `bin/ucmuon_transport_bb_omp` or `gui/ucmuon_bb_driver.py` | CSDA + b·E radiative, Highland MS — deterministic, no fluctuations |
| UCMuon | `gui/ucmuon_stochastic_driver.py` | Groom 2001 range table, Poisson hard events, Highland MS, decay |
| CSDA | `benchmark_analysis.py` | Deterministic Groom 2001 range threshold — analytical upper bound |

---

## Benchmark Geometry Specification

| Parameter | Value | Notes |
|---|---|---|
| Material | Standard Rock | Z=11, A=22, ρ=2.65 g/cm³, I=136.4 eV, X₀=26.7 cm |
| Geometry | Infinite flat slab | One horizontal layer |
| Zenith angle (Phase 1) | θ = 0° | Path = depth exactly; cleanest comparison |
| Zenith angle (Phase 2) | cos²θ, θ_max=85° | Realistic sea-level distribution |
| E_min | **50 GeV** | See statistics note below |
| E_max | 2500 GeV | CosmoALEPH upper limit |
| Spectrum | CosmoALEPH | dN/dp ∝ p⁻³·¹⁹⁵² |
| N muons Phase 1 | **500,000** | Gives ≥ 650 survivors at D5 |
| N muons Phase 2 | **1,000,000** | Angular spread reduces per-bin statistics |
| Source shape | Circular disk, r=50 m | Avoids geometric edge effects |

> **Why E_min = 50 GeV, not 1 GeV?**
> The CSDA stopping threshold for 50 m of Standard Rock is ~29 GeV.
> Any muon below this energy is guaranteed to stop before reaching even the shallowest
> benchmark point. Generating them wastes CPU and dilutes statistics with zeroes.
> E_min = 50 GeV ensures every muon contributes meaningful physics information.

**Benchmark depth points — statistics with N_gen = 500,000:**

| Label | Depth | Opacity (g/cm²) | CSDA E_min | CSDA f_surv | Expected N_alive | σ_binomial |
|---|---|---|---|---|---|---|
| D1 | 50 m | 13,250 | 29 GeV | ~1.000 | ~500,000 | 0.1% |
| D2 | 100 m | 26,500 | 62 GeV | 0.623 | ~312,000 | 0.2% |
| D3 | 200 m | 53,000 | 133 GeV | 0.116 | ~58,000 | 0.4% |
| D4 | 500 m | 132,500 | 387 GeV | 0.011 | ~5,500 | 1.4% |
| D5 | 1000 m | 265,000 | 969 GeV | 0.0013 | ~650 | 3.9% |

D5 is at the edge of the Groom 2001 range table and has larger statistical uncertainty.
Include it in the paper with explicit error bars; note that the UCMuon Stochastic engine
is less reliable above ~800 m.w.e.

---

## Prerequisites

```bash
cd /Users/basiri/UCMoun_OMP
source ~/venvs/ucmuon/bin/activate

# Verify binaries
ls -lh bin/ucmuon_gen_omp bin/ucmuon_transport_music_omp bin/ucmuon_transport_bb_omp

# Verify MUSIC data files
ls -lh music-eloss-rock.dat music-double-diff-rock.dat

# Physics sanity check — UCMuon Stochastic dE/dx at 1 GeV
python -c "
import sys; sys.path.insert(0,'gui')
from ucmuon_stochastic_driver import _dedx
d = _dedx(1000.0)
print(f'dEdx at 1 GeV = {d:.4f} MeV cm2/g  (expected ~1.79–1.82)')
assert 1.70 < d < 1.95
print('OK')
"
```

If binaries are missing: `make omp` or `bash setup.sh`.

---

## Step 0 — Create Directory Structure

```bash
mkdir -p benchmark/phase1 benchmark/phase2 benchmark/results
```

---

## Step 1 — Generate Surface Muon File

One file, shared by all engines and all depths.

### Via GUI (Tab 1)

| Setting | Value |
|---|---|
| Spectrum model | ① CosmoALEPH |
| E_min | **50 GeV** |
| E_max | 2500 GeV |
| Source shape | 💿 Circular disk |
| Radius | 50 m |
| z plane | 0 m |
| N muons | **500,000** |
| Angular mode | **① Vertical only (θ=0°)** |
| Detector filter | **OFF** |
| Output file | `benchmark/phase1/muons_surface.dat` |

### Verify

```bash
wc -l benchmark/phase1/muons_surface.dat
# ~500,001 lines (1 header + 500,000 muons)

# Check all thetas are ~0 (vertical only)
awk 'NR>1 {print $9}' benchmark/phase1/muons_surface.dat | \
  python -c "import sys,math; t=[float(l) for l in sys.stdin]; \
  print(f'theta: min={min(t):.4f} max={max(t):.4f} rad  (all should be ~0)')"

# Check mean energy (should be ~100–200 GeV for steep spectrum from 50–2500 GeV)
awk 'NR>1 {sum+=$11; n++} END {printf "Mean E = %.1f GeV\n", sum/n}' \
  benchmark/phase1/muons_surface.dat
```

---

## Step 2 — Run Transport for All Engines and Depths

Save as `benchmark/run_phase1.sh` and run from the **project root**:

```bash
#!/usr/bin/env bash
#=============================================================
# UCMuon Phase 1 Benchmark  —  all engines x all depths
# Usage: bash benchmark/run_phase1.sh
# Run from: /Users/basiri/UCMoun_OMP
#=============================================================
set -euo pipefail

source ~/venvs/ucmuon/bin/activate
export OMP_NUM_THREADS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"

SURF="benchmark/phase1/muons_surface.dat"
RHO=2.65
X0=26.7       # Standard Rock radiation length [cm]

for DEPTH in 50 100 200 500 1000; do
  OPC=$(python -c "print(int($DEPTH * $RHO * 100))")
  echo ""
  echo "==================================================="
  echo "  DEPTH = ${DEPTH} m   (opacity = ${OPC} g/cm²)"
  echo "==================================================="

  # ── MUSIC ────────────────────────────────────────────────
  OUT="benchmark/phase1/ug_music_${DEPTH}m.dat"
  LOG="benchmark/phase1/log_music_${DEPTH}m.txt"
  echo "  → MUSIC ..."
  # stdin fields: infile / outfile / transport_all=1 / ncols_hint=13 /
  #               depth_m / mat_type=1(rock) / rho / rad / init=1 / [Enter]
  printf "%s\n%s\n1\n13\n%s\n1\n%s\n2.864e4\n1\n\n" \
      "$SURF" "$OUT" "$DEPTH" "$RHO" \
      | ./bin/ucmuon_transport_music_omp >"$LOG" 2>&1
  ALIVE=$(awk 'NR>1 && $9==1' "$OUT" | wc -l | tr -d ' ')
  TOTAL=$(awk 'NR>1' "$OUT" | wc -l | tr -d ' ')
  echo "     Survived $ALIVE / $TOTAL  ($(python -c "print(f'{$ALIVE/$TOTAL:.4f}')"))  → $OUT"

  # ── Bethe-Bloch ──────────────────────────────────────────
  OUT="benchmark/phase1/ug_bb_${DEPTH}m.dat"
  LOG="benchmark/phase1/log_bb_${DEPTH}m.txt"
  echo "  → Bethe-Bloch ..."
  # stdin: infile / outfile / transport_all=1 / ncols_hint=13 / depth_m / mat_type=1 / ms_enable=1
  printf "%s\n%s\n1\n13\n%s\n1\n1\n" \
      "$SURF" "$OUT" "$DEPTH" \
      | ./bin/ucmuon_transport_bb_omp >"$LOG" 2>&1
  ALIVE=$(awk 'NR>1 && $9==1' "$OUT" | wc -l | tr -d ' ')
  TOTAL=$(awk 'NR>1' "$OUT" | wc -l | tr -d ' ')
  echo "     Survived $ALIVE / $TOTAL  ($(python -c "print(f'{$ALIVE/$TOTAL:.4f}')"))  → $OUT"

  # ── UCMuon Stochastic ─────────────────────────────────────
  OUT="benchmark/phase1/ug_ucmuon_${DEPTH}m.dat"
  LOG="benchmark/phase1/log_ucmuon_${DEPTH}m.txt"
  echo "  → UCMuon Stochastic ..."
  # stdin: infile / outfile / depth_m / rho / X0_cm / mat_id=1 /
  #        transport_all=1 / ncols=13 / n_steps=0(auto) / v_cut=0.05 / ms_enable=1
  printf "%s\n%s\n%s\n%s\n%s\n1\n1\n13\n0\n0.05\n1\n" \
      "$SURF" "$OUT" "$DEPTH" "$RHO" "$X0" \
      | python gui/ucmuon_stochastic_driver.py >"$LOG" 2>&1
  ALIVE=$(awk 'NR>1 && $9==1' "$OUT" | wc -l | tr -d ' ')
  TOTAL=$(awk 'NR>1' "$OUT" | wc -l | tr -d ' ')
  echo "     Survived $ALIVE / $TOTAL  ($(python -c "print(f'{$ALIVE/$TOTAL:.4f}')"))  → $OUT"

done

echo ""
echo "==================================================="
echo "  ALL PHASE 1 RUNS COMPLETE"
echo "==================================================="
ls -lh benchmark/phase1/ug_*.dat
```

```bash
chmod +x benchmark/run_phase1.sh
bash benchmark/run_phase1.sh 2>&1 | tee benchmark/run_phase1.log
```

**If any run fails:** check `benchmark/phase1/log_{engine}_{depth}m.txt`.
The most common failure is a missing MUSIC data file or a wrong stdin field count.

> **Note on MUSIC `init` parameter:** Use `init=1` (use cached cross-section tables) after
> the very first run. If you are running MUSIC at a new depth for the first time and the
> cross-section cache does not yet exist, change `init=1` to `init=0` for that run, then
> revert to `init=1`. The cache file is `music-cross-sections-rock.dat`.

---

## Step 3 — Spot-Check Before Full Analysis

```bash
# Quick survival rate table
echo "Engine        50m    100m   200m   500m  1000m"
echo "----------------------------------------------------"
for ENGINE in music bb ucmuon; do
  ROW="$(printf '%-12s' "$ENGINE")"
  for DEPTH in 50 100 200 500 1000; do
    FILE="benchmark/phase1/ug_${ENGINE}_${DEPTH}m.dat"
    if [ -f "$FILE" ]; then
      N_A=$(awk 'NR>1 && $9==1' "$FILE" | wc -l | tr -d ' ')
      N_T=$(awk 'NR>1' "$FILE" | wc -l | tr -d ' ')
      F=$(python -c "print(f'{$N_A/$N_T:.3f}')")
      ROW="$ROW  $F"
    else
      ROW="$ROW   N/A"
    fi
  done
  echo "$ROW"
done
```

**Red flags that indicate a problem:**
- Any engine at 50 m gives f_surv < 0.95 → path-length or E_min error
- Any engine at 1000 m gives f_surv > 0.05 → depth unit error (metres vs cm?)
- MUSIC and BB differ by > 30% at 100 m → dEdx parameter mismatch
- All three engines give identical numbers → they may be reading the same output file

---

## Step 4 — Run the Analysis Script

```bash
cd /Users/basiri/UCMoun_OMP
python benchmark_analysis.py \
    --phase 1 \
    --indir benchmark/phase1 \
    --outdir benchmark/results
```

**Outputs:**

| File | Contents |
|---|---|
| `phase1_survival_fraction.pdf` | **Main figure** — f_surv vs opacity, all engines + CSDA ref |
| `phase1_relative_deviation.pdf` | (Engine − MUSIC) / MUSIC in % |
| `phase1_mean_energy.pdf` | Mean output energy vs depth |
| `phase1_spectrum_50m.pdf` ... | Output energy spectra at each depth |
| `phase1_benchmark_table.csv` | Numerical results |
| `phase1_benchmark_table.tex` | Ready-to-paste LaTeX table for the paper |

---

## Step 5 — Phase 2: cos²θ Distribution

### Generate Phase 2 surface file

Tab 1 changes from Phase 1:

| Setting | Phase 1 | Phase 2 |
|---|---|---|
| Angular mode | ① Vertical only | **② cos²θ (recommended)** |
| θ_max | — | **85°** |
| N muons | 500,000 | **1,000,000** |
| Output | `phase1/muons_surface.dat` | `phase2/muons_surface.dat` |

### Copy and modify the run script

```bash
sed 's|phase1|phase2|g' benchmark/run_phase1.sh > benchmark/run_phase2.sh
chmod +x benchmark/run_phase2.sh
bash benchmark/run_phase2.sh 2>&1 | tee benchmark/run_phase2.log

python benchmark_analysis.py \
    --phase 2 \
    --indir benchmark/phase2 \
    --outdir benchmark/results
```

### What Phase 2 tests

Each engine must compute the slant path `L = depth / cos(θ)` to apply the correct
overburden. A discrepancy not present in Phase 1 points to an angular correction bug in
one engine. To isolate it, extract near-vertical muons (θ < 3°) from the Phase 2 output
and compare them against Phase 1 results — they must match within statistical errors.

---

## Step 6 — Acceptance Criteria

| Comparison | Metric | Pass criterion | Physical reason |
|---|---|---|---|
| BB vs MUSIC | f_surv | ≤ 15% relative at D1–D4 | Known CSDA vs stochastic difference |
| UCMuon vs MUSIC | f_surv | ≤ 20% at D1–D3, ≤ 40% at D4–D5 | Simplified hard-event model |
| All engines | f_surv ≥ CSDA | Always true | Fluctuations allow survival past threshold |
| BB vs MUSIC | Mean E_out | ≤ 10% relative | Mean energy is robust to fluctuations |
| Phase 2 vs Phase 1 | Vertical sub-sample | Within Phase 1 tolerance | Angular correction consistency |

---

## Step 7 — Archive

```bash
# Copy outputs to docs/
cp benchmark/results/phase1_benchmark_table.csv docs/
cp benchmark/results/phase1_benchmark_table.tex docs/

# Tag the code at this benchmark version
git add benchmark/ docs/
git commit -m "Phase 1 benchmark: Standard Rock flat-slab (5 depths)"
git tag v1.0.0-benchmark
git push origin main --tags
```

**Paper methods template:**
> "The three transport engines were benchmarked against MUSIC (Kudryavtsev 2009,
> CPC 180, 339) using a Standard Rock flat-slab geometry (ρ = 2.65 g/cm³),
> 500,000 vertical muons from the CosmoALEPH spectrum (E ∈ [50, 2500] GeV),
> at five overburden depths (13,250–265,000 g/cm²). Survival fractions agreed
> to within X% for the Bethe-Bloch engine and Y% for UCMuon Stochastic for
> opacities up to Z g/cm²."

---

## Output File Format Reference (18 columns)

All three engines write an identical format:

```
Col  0: EventID
Col  1: x_srf_cm     Col  2: y_srf_cm     Col  3: z_srf_cm
Col  4: E_srf_GeV    ← surface (input) energy
Col  5: theta_srf    Col  6: phi_srf  [rad]
Col  7: charge       ← +1 / −1
Col  8: alive        ← 1 = survived,  0 = stopped   ← KEY COLUMN
Col  9: x_ug_cm      Col 10: y_ug_cm    Col 11: z_ug_cm
Col 12: E_ug_GeV     ← output energy (0.0 if stopped)
Col 13: cx_ug        Col 14: cy_ug      Col 15: cz_ug
Col 16: theta_ug     Col 17: phi_ug  [rad]
```

Survival fraction  = `sum(col[8] == 1) / N_total`  
Mean output energy = `mean(col[12][ col[8] == 1 ])`
