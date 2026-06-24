# MuonRock Benchmark — Analysis Feedback Report
**UCLouvain / CP3**
**Date:** 2026-05-23
**Geant4 run:** `outputs/run_file_20260522_210114`  (600,000 events, 18878 s)
**Reference code:** Geant4 11 + FTFP_BERT + G4EmStandardPhysics_option4

---

## 1. Source File Verification

All six codes were confirmed to use the correct, consistent input.

| Code | Source file | Format | Energy convention | Status |
|------|-------------|--------|-------------------|--------|
| MUSIC | `sources/benchmark_surface.dat` | 14-col ASCII | E_total (GeV) in col 11 | ✓ Correct |
| PROPOSAL | `sources/benchmark_surface.dat` | same | same | ✓ Correct |
| BB | `sources/benchmark_surface.dat` | same | same | ✓ Correct |
| UCMuon | `sources/benchmark_surface.dat` | same | same | ✓ Correct |
| PHITS | `sources/benchmark_phits.dat` | 10-col Fortran D-notation | KE (MeV) in col 8 | ✓ Correct |
| Geant4 | `sources/benchmark_geant4_new.txt` | Geant4 FILE-mode | KE (MeV) in col 8 | ✓ Correct |

**PHITS source detailed check:**
- 600,000 muons parsed; 6 bins × 100,000 each ✓
- KE bins (MeV): 4894, 9894, 19894, 49894, 99894, 299894 — each equals E_total − m_μ ✓
- Direction cosines: cx=cy=0, cz=−1 (all vertical downward) ✓
- Entry (x,y) positions: exact match to `benchmark_surface.dat` row-by-row ✓
- z=0 for all muons (surface injection) ✓
- PDG: −13 (mu+) and +13 (mu−) match surface.dat charge signs ✓

**`benchmark_surface.dat` is the common source** for all codes.
PHITS and Geant4 use kinematic energy (KE = E_total − 0.10566 GeV); MUSIC/PROPOSAL/BB/UCMuon
receive the total energy column, which their transport engines correctly interpret.

---

## 2. Bug Fixes Applied to benchmark_analysis.py

### 2.1 Energy convention in `.dat` loader (critical)
**Problem:** `_load_dat_files()` stored `Es` (total surface energy, GeV) and `E` (total underground energy, GeV) directly as `InitKE_GeV` and `ExitKE_GeV`.
These fields are total energies in the 18-column `.dat` format, not kinetic energies.
The error was +0.10566 GeV on every muon's stored KE — 2.2% at 5 GeV, negligible at 300 GeV.
**Fix:** `InitKE_GeV = Es − M_MU`, `ExitKE_GeV = max(E − M_MU, 0)`.
Energy *loss* was unaffected (mass cancels in the difference).

### 2.2 Geant4 ELoss accounting (critical for energy loss figures)
**Problem:** The Geant4 CSV column `ELossTotalGeV` = sum of `GetTotalEnergyDeposit()` over
muon steps. This is the **locally deposited** energy at the muon track — it does NOT include
energy carried away by bremsstrahlung photons and e⁺e⁻ pairs which are tracked as separate
particles. As a result:

| Depth | `ELossTotal` (local deposit) | `InitKE − ExitKE` (kinematic) | Ratio |
|-------|------------------------------|-------------------------------|-------|
| 2.65 MWE | 0.432 GeV | 0.649 GeV | 1.50× |
| 26.5 MWE | 4.313 GeV | 6.539 GeV | 1.52× |
| 66.25 MWE | 10.825 GeV | 17.389 GeV | 1.61× |
| 132.5 MWE | 21.681 GeV | 36.333 GeV | 1.68× |
| 265.0 MWE | 43.382 GeV | 75.871 GeV | 1.75× |
| 530.0 MWE | 86.797 GeV | 162.303 GeV | 1.87× |

The ratio grows with depth because radiative losses become more important at higher energies (the surviving muons at deep planes are increasingly energetic).
**Fix:** `ELoss_Total_GeV` replaced with `InitKE − ExitKE` for Geant4, making it directly comparable to all other codes. The per-process fractions (Ion/Brem/Pair/Nucl) retain their original local-deposit basis and are kept internally consistent by using the sum of process columns as their own denominator.

### 2.3 Error bars in MCS angle plot (cosmetic)
**Problem:** Error bars on mean MCS angle used population σ (spread of the distribution), not SEM (σ/√N). With ~50k–100k muons per point, σ ≈ 0.3–5° while SEM ≈ 0.001–0.007°. The large σ made the 26.5 MWE point look unreliable.
**Fix:** Error bars now show SEM.

### 2.4 fig04 right panel (misleading)
**Problem:** Plotted mean `sqrt(x²+y²)` from beam axis. This depends on the surface-sampling
geometry (±2500 cm slab for grid codes, same ±2500 cm for Geant4 but with finite-slab boundary
loss). Geant4's trend was opposite to UCMuon's because edge muons drift outside the 50m×50m slab
and are lost.
**Fix:** Changed to lateral displacement from entry position `sqrt(Δx²+Δy²)` for codes that record
entry positions (MUSIC, PROPOSAL, BB, UCMuon). Geant4 excluded from this panel (entry position
not in output CSV).

### 2.5 fig06 dE/dx left panel (functional)
**Problem:** Used deepest depth (530 MWE) where only one or two energy bins survive. Grid codes
silently produced zero points due to `pd.cut` failing on a degenerate bin (lo==hi).
**Fix:** Use shallowest depth (2.65 MWE, 1m) so all 6 energy bins contribute. Added guard for
lo≥hi edge case. All five codes now plot correctly.

### 2.6 fig13 caption overflow (cosmetic)
Grid-scan note moved from `suptitle` to `fig.text()` at figure bottom.

### 2.7 fig17 x-axis (cosmetic)
Bethe-Bloch reference clipped to source max energy (400 GeV). Consistent `xlim(3, 400)` across
all 6 subplots.

---

## 3. PHITS Input Review (`muon_rock.inp` — reviewed 2026-05-23)

### 3.1 Fixes applied to `muon_rock.inp`

| # | Parameter | Before | After | Impact |
|---|-----------|--------|-------|--------|
| 1 | `emin(5)` | `= 1.0` with comment "µ+ cutoff" | **removed** (emin(5) is pion−, defaults to 1.0 MeV) | None — pion− already at default |
| 2 | `emin(6)` | `= 1.0` with comment "µ− cutoff" | `= 1.0` with corrected comment "µ+ cutoff (itype=6)" | Comment fix only |
| 3 | `emin(7)` | missing → defaulted to **0.001 MeV** | `= 1.0` **µ− cutoff added** | Negligible (all muons > 5 GeV) |
| 4 | Tally 7 xy-map (200 m) | `output = current` | `output = f-curr` | Correct: forward-only flux |
| 5 | Tally 8 xy-map (1 m) | `output = current` | `output = f-curr` | Correct: forward-only flux |

### 3.2 Physics parameters confirmed correct

| Parameter | PHITS value | Geant4 equivalent | Status |
|-----------|------------|-------------------|--------|
| Muon bremsstrahlung | `imubrm = 1` | `G4MuBremsstrahlung` | ✓ |
| Muon pair production | `imuppd = 1` | `G4MuPairProduction` | ✓ |
| Muon nuclear reaction | `imuint = 1` | `G4MuonNuclearProcess` | ✓ |
| Photo-nuclear | `ipnint = 1` | `G4PhotoNuclearProcess` | ✓ |
| Stopping power | `ndedx = 3` (ATIMA) | NIST Bethe-Bloch | ✓ |
| Energy straggling | `nedisp = 1` (Landau-Vavilov) | `G4IonFluctuations` | ✓ |
| MCS | `nspred = -2` (Lynch-Molière) | `G4UrbanMscModel` | ⚠️ model differs |
| Muon capture | `imucap = 0` | not applicable (GeV muons) | ✓ |
| Geometry (XY, depths, density) | ±2500 cm, 6 planes, 2.65 g/cm³ | identical | ✓ |
| Total events | 60 000 × 10 batches = 600 000 | 600 000 | ✓ |

### 3.3 Known unfixable differences (PHITS architecture)

**a) No per-event output.**
All PHITS tallies are aggregate histograms (T-Cross, T-Deposit). There is no per-muon CSV equivalent to Geant4's `_muons.csv`. Per-muon InitKE, ExitKE, scatter angle, and position are unavailable. This is why MCS angles and lateral displacement show `—` for PHITS throughout the benchmark. A fundamentally different tally setup (PHITS dump output at each depth plane) would be required.

**b) MCS model: Lynch-Molière vs Urban MscModel.**
PHITS `nspred=-2` uses the Lynch-Molière parametrisation; Geant4 FTFP_BERT uses `G4UrbanMscModel` (Urban 2006). Both are derived from Molière multiple scattering theory but differ in step-limited accumulation and the handling of large-angle tails. This is a fundamental model difference, not a configuration error, and contributes to the ~8% MCS angle offset observed between Geant4 and MUSIC (which uses a similar formulation to Geant4).

**c) Material mean excitation energy I.**
Geant4 defines rock as a single-element pseudo-material (Z=11, A=22) → Geant4 assigns I from the NIST table for Z=11 (sodium): **I = 149 eV**. PHITS uses a Mg(75%)+O(25%) mixture → composite I ≈ **138 eV**. The 8% difference in I leads to < 0.3% change in ionisation dE/dx via the Bethe-Bloch logarithm at 300 GeV. This is **not** the cause of the −12% KE discrepancy.

**d) PHITS dE/dx discrepancy (−12.2% exit KE at 530 MWE).**
The gap grows monotonically with depth: −1.4% (66 MWE) → −3.1% (132 MWE) → −6.3% (265 MWE) → −12.2% (530 MWE). This signature — accelerating divergence at depth — points to higher radiative energy loss in PHITS for high-energy (> 100 GeV) muons. The most likely cause is different cross-section parametrisations for muon bremsstrahlung and pair production compared to Geant4's Kelner-Kokoulin-Petrukhin implementation. To confirm: run a single-layer PHITS test with monoenergetic 300 GeV muons at 200 m and compare the mean dE to Geant4.

---

## 4. Benchmark Results (run v2 — 2026-05-25)

> **Run v2 changes vs v1:** PROPOSAL MCS enabled (scattering was disabled in v1);
> UCMuon lateral displacement bug fixed (stochastic position kicks removed).
> All other codes unchanged. Numbers below are from the v2 run.

### 4.1 Transmission (%)

| MWE | Geant4 | PHITS | MUSIC | PROPOSAL | BB | UCMuon |
|----:|-------:|------:|------:|---------:|---:|-------:|
| 2.65 | **99.998** | 99.995 | 100.000 | 99.999 | 100.000 | 99.905 |
| **26.5** | **92.215** | 83.230 | 87.386 | 85.478 | 83.333 | 82.733 |
| 66.25 | **66.347** | 66.328 | 66.507 | 66.397 | 66.667 | 65.705 |
| 132.5 | **49.618** | 49.565 | 49.807 | 49.720 | 50.000 | 49.166 |
| 265.0 | **32.716** | 32.589 | 32.916 | 32.820 | 33.333 | 32.312 |
| 530.0 | **16.178** | 16.406 | 16.314 | 16.270 | 16.667 | 15.985 |

**Assessment:**
- ✅ All codes agree within ±2% at ≥66 MWE.
- ⚠️ **26.5 MWE outlier (Geant4 +7–9 pp above all other codes):** genuine dE/dx model difference.
  Geant4's 5 GeV muon range ≈ 996 cm; ~54% pass the 10m plane. All other codes stop all 5 GeV
  muons before 10m. Effect confined to this one depth; all codes agree at ≥66 MWE.
  This also affects MCS angles and lateral displacement at 26.5 MWE — that point should not be
  used for MCS comparisons.
- UCMuon −1.2% at 530 MWE (15.985% vs 16.178%) — improved from v1 (was −4.4%).

### 4.2 Mean Exit Kinetic Energy (GeV)

| MWE | Geant4 | PHITS | MUSIC | PROPOSAL | BB | UCMuon |
|----:|-------:|------:|------:|---------:|---:|-------:|
| 2.65 | **80.081** | 80.445 | 80.045 | 80.042 | 80.060 | 80.221 |
| 26.5 | **80.576** | 88.627 | 84.783 | 86.638 | 89.041 | 89.516 |
| 66.25 | **100.396** | 98.972 | 99.529 | 99.617 | 99.723 | 100.522 |
| 132.5 | **114.121** | 110.627 | 112.506 | 112.471 | 113.089 | 113.418 |
| 265.0 | **125.302** | 117.430 | 122.171 | 122.162 | 123.099 | 123.740 |
| 530.0 | **137.592** | 120.772 | 131.970 | 131.615 | 135.181 | 132.796 |

Deviation from Geant4 (code − G4) / G4:

| MWE | PHITS | MUSIC | PROPOSAL | BB | UCMuon |
|----:|------:|------:|---------:|---:|-------:|
| 2.65 | +0.5% | −0.0% | −0.0% | −0.0% | +0.2% |
| 66.25 | −1.4% | −0.9% | −0.8% | −0.7% | +0.1% |
| 132.5 | −3.1% | −1.4% | −1.4% | −0.9% | −0.6% |
| 265.0 | −6.3% | −2.5% | −2.5% | −1.8% | −1.2% |
| **530.0** | **−12.2%** | **−4.1%** | **−4.3%** | **−1.8%** | **−3.5%** |

**Assessment:**
- ✅ At 2.65 MWE all codes agree within 0.5%.
- ✅ **BB: ±2% at all depths** — best energy agreement.
- ✅ **UCMuon: ±3.5% at all depths** (v2 fixed; was +1.1% in v1 at 530 MWE, now −3.5%).
  UCMuon energy loss now consistent with MUSIC/PROPOSAL — the v1 inflated exit KE was tied
  to the same stochastic position bug that caused the lateral displacement excess.
- ⚠️ MUSIC and PROPOSAL both −4.1 to −4.3% at 530 MWE: consistent with each other,
  indicating similar radiative cross-section parametrisations (higher than Geant4).
- 🔴 **PHITS −12.2% at 530 MWE** — unchanged. Gap grows monotonically with depth.
  Most likely cause: different bremsstrahlung or pair-production cross sections at 100–300 GeV.
  Run a dedicated single-layer test (300 GeV muons, 200 m) to isolate the process.

### 4.3 Mean MCS Deflection Angle (°)

| MWE | Geant4 | MUSIC | PROPOSAL | BB | UCMuon | PHITS |
|----:|-------:|------:|---------:|---:|-------:|------:|
| 2.65 | **0.212** | 0.218 | 0.216 | 0.130 | 0.114 | — |
| 26.5 | **1.361** | 1.395 | 0.887* | 0.321 | 0.307 | — |
| 66.25 | **0.547** | 0.602 | 0.632 | 0.351 | 0.336 | — |
| 132.5 | **0.386** | 0.414 | 0.445 | 0.237 | 0.238 | — |
| 265.0 | **0.370** | 0.403 | 0.440 | 0.223 | 0.238 | — |
| 530.0 | **0.263** | 0.284 | 0.311 | 0.146 | 0.181 | — |

*26.5 MWE anomaly: PROPOSAL stops all 5 GeV muons before 10m; the surviving population
is higher-energy and less scattered than Geant4's surviving population, giving a spuriously
low mean angle. This point is excluded from physics comparisons (see §4.1).

MCS ratio (code / Geant4) at ≥66 MWE:

| MWE | MUSIC | PROPOSAL | BB | UCMuon |
|----:|------:|---------:|---:|-------:|
| 66.25 | +10% | **+16%** | −36% | −39% |
| 132.5 | +7% | **+15%** | −39% | −38% |
| 265.0 | +9% | **+19%** | −40% | −36% |
| 530.0 | +8% | **+18%** | −44% | −31% |

**Assessment:**
- ✅ **MUSIC: +3–10% of Geant4** — best match, within 10% at all depths.
- ⚠️ **PROPOSAL: +15–19% overestimate** at ≥66 MWE (v2, first run with MCS enabled).
  The overestimate is consistent across all deep depths, pointing to a systematic error in the
  radiation length X₀ used in PROPOSAL's scattering model. If X₀ is ~31% too small,
  θ ∝ 1/√X₀ would give a +18% overestimate.
  **Action: verify X₀ for standard rock in the PROPOSAL medium definition. Correct value:
  X₀ ≈ 26.7 g/cm² for Z=11, A=22, ρ=2.65 g/cm³.**
- ⚠️ **BB: 56–64% of Geant4** — consistent underestimate. Known limitation of single-plane
  Highland formula for thick absorbers (200 m). Do not tune to match.
- ⚠️ **UCMuon: 54–61% of Geant4** — similar to BB, consistent with PUMAS MCS implementation.
  MCS angles unchanged from v1 (only position tracking was fixed).
- ❌ **PHITS: no per-event angle data** (aggregate tallies only).
- **26.5 MWE peak (1.36° in Geant4):** expected physics — near-stopping 5 GeV muons
  scatter at large angles, inflating Geant4's mean. Not usable for code comparison.

### 4.4 Mean Lateral Displacement from Entry Position (cm)

*(Geant4/PHITS excluded — entry position not in their output)*

| MWE | MUSIC | PROPOSAL | BB | UCMuon |
|----:|------:|---------:|---:|-------:|
| 2.65 | 0.23 | 0.22 | 0.13 | 0.12 |
| 26.5 | 6.37 | 5.52 | 2.80 | 2.68 |
| 66.25 | 11.07 | 11.54 | 6.71 | 6.41 |
| 132.5 | 16.31 | 17.29 | 9.84 | 9.56 |
| 265.0 | 28.88 | 31.30 | 17.35 | 17.27 |
| 530.0 | 41.34 | 45.11 | 24.02 | 25.72 |

**Assessment:**
- Two groups emerge:
  - **MUSIC + PROPOSAL**: 41–45 cm at 530 MWE. Both have larger MCS angles (+8 to +18% vs Geant4),
    which drives more lateral spread. The two codes are consistent with each other.
  - **BB + UCMuon**: 24–26 cm at 530 MWE. Both have reduced MCS angles (−31 to −44% vs Geant4),
    giving proportionally less lateral drift. UCMuon is now internally consistent (v1 had 65 cm
    despite smaller MCS angles — the bug is fixed).
- None of the four codes can be compared directly to Geant4 for lateral displacement because
  Geant4 does not record entry positions in its CSV output.

### 4.5 Mean Energy Loss (GeV) — kinematic InitKE − ExitKE

| MWE | Geant4 | MUSIC | PROPOSAL | BB | UCMuon |
|----:|-------:|------:|---------:|---:|-------:|
| 2.65 | **0.649** | 0.683 | 0.686 | 0.668 | 0.572 |
| 26.5 | **6.539** | 6.889 | 6.969 | 6.853 | 6.902 |
| 66.25 | **17.389** | 18.076 | 18.161 | 17.671 | 18.019 |
| 132.5 | **36.333** | 37.693 | 37.892 | 36.805 | 37.748 |
| 265.0 | **75.871** | 78.670 | 78.983 | 76.796 | 78.519 |
| 530.0 | **162.303** | 167.924 | 168.279 | 164.714 | 167.098 |

Note: PHITS energy loss not available (summary-only, no per-muon InitKE).

**Assessment:**
- MUSIC, PROPOSAL, and UCMuon (v2) now all give nearly identical energy loss at deep depths
  (~167–168 GeV at 530 MWE), ~3–4% above Geant4.
- BB is closest to Geant4 at deep depths (+1.5% at 530 MWE).
- UCMuon at 2.65 MWE (0.572 GeV) is −12% below Geant4 (0.649 GeV). This is a UCMuon
  shallow-depth anomaly — possibly a threshold/cutoff effect. Warrants investigation but
  does not affect deep-depth results.

### 4.6 Geant4 Energy Loss Process Fractions (local-deposit basis)

| MWE | Ion (%) | Brem (%) | Pair (%) | Nucl (%) |
|----:|--------:|---------:|---------:|---------:|
| 2.65 | 97.70 | 0.12 | 2.14 | 0.03 |
| 26.5 | 97.11 | 0.15 | 2.70 | 0.04 |
| 66.25 | 96.33 | 0.16 | 3.45 | 0.05 |
| 132.5 | 95.71 | 0.17 | 4.07 | 0.05 |
| 265.0 | 95.07 | 0.17 | 4.70 | 0.06 |
| 530.0 | 94.09 | 0.18 | 5.66 | 0.07 |

**Caution:** Fractions from locally deposited energy at muon steps only. True radiative fraction
is substantially higher (~30% at 300 GeV) because brem photons and e⁺e⁻ pairs escape the step.
Use for relative depth comparisons only, not as absolute process fractions.

### 3.6 Geant4 Energy Loss Process Fractions (local-deposit basis)

| MWE | Ion (%) | Brem (%) | Pair (%) | Nucl (%) |
|----:|--------:|---------:|---------:|---------:|
| 2.65 | 97.70 | 0.12 | 2.14 | 0.03 |
| 26.5 | 97.11 | 0.15 | 2.70 | 0.04 |
| 66.25 | 96.33 | 0.16 | 3.45 | 0.05 |
| 132.5 | 95.71 | 0.17 | 4.07 | 0.05 |
| 265.0 | 95.07 | 0.17 | 4.70 | 0.06 |
| 530.0 | 94.09 | 0.18 | 5.66 | 0.07 |

**Caution:** These fractions are computed from locally deposited energy at muon track steps.
Bremsstrahlung photons and pair-produced e⁺e⁻ escape the current step and deposit elsewhere;
their energy is NOT counted here. The true radiative fraction is substantially higher (order 30%
for 300 GeV muons). These numbers are useful for relative comparisons across depths but should
not be used as absolute process fractions.

---

## 5. Code Acceptability for Subsurface Muography

### Transmission (within 5% of Geant4)
| Code | Acceptable from |
|------|----------------|
| PHITS | ≥66 MWE (26.5 MWE fails: 0.90×) |
| MUSIC | ≥66 MWE (26.5 MWE: 0.95×, model boundary) |
| PROPOSAL | ≥66 MWE (26.5 MWE: 0.93×) |
| BB | ≥66 MWE (26.5 MWE: 0.90×) |
| UCMuon | ≥66 MWE (26.5 MWE: 0.90×; −4.4% at 530 MWE) |

All codes agree at ≥66 MWE for relative muon counting.

### Mean Exit Energy (within 5% of Geant4)
| Code | Acceptable from | Notes |
|------|----------------|-------|
| UCMuon | All depths | ±1% |
| BB | All depths | ±2% |
| MUSIC | ≥66 MWE | −4.2% at 530 MWE |
| PROPOSAL | ≥66 MWE | −4.3% at 530 MWE |
| PHITS | ≤66 MWE only | −12.2% at 530 MWE — investigate dE/dx |

### MCS / Angle Reconstruction (v2)
| Code | Usable for angle reconstruction | Deviation vs Geant4 |
|------|--------------------------------|---------------------|
| Geant4 | ✅ Reference | — |
| MUSIC | ✅ Within 10% | +3 to +10% |
| PROPOSAL | ⚠️ Overestimates scatter | **+15 to +19%** — check X₀ in medium definition |
| BB | ⚠️ Underestimates scatter | −39 to −44% |
| UCMuon | ⚠️ Underestimates scatter | −31 to −46% |
| PHITS | ❌ No per-event angles | summary-only output |

---

## 6. Action Items (updated v2)

| Priority | Item | Code | Status |
|----------|------|------|--------|
| 🔴 High | PROPOSAL MCS +15–19% overestimate at ≥66 MWE. Verify radiation length X₀ in the PROPOSAL medium definition for standard rock. Correct value: X₀ ≈ 26.7 g/cm² (Z=11, A=22, ρ=2.65). Print `medium.radiation_length` in the wrapper. | PROPOSAL | Open |
| 🔴 High | PHITS −12.2% exit KE at 530 MWE. Run single-layer test: 300 GeV muons, 200 m depth. Compare mean dE/dx from T-Deposit vs Geant4. Likely cause: different brem/pair cross sections at 100–300 GeV. | PHITS | Open |
| 🟡 Medium | Add wall-clock timing files for all grid codes: `<Code>_timing.txt` with line `Elapsed : <seconds>`. Required for fig10 speed comparison. | MUSIC, PROPOSAL, BB, UCMuon | Open |
| 🟡 Medium | Add physics model headers to .dat files (brem model, pair model, MCS model, X₀). Needed for benchmark paper §3. | All grid codes | Open |
| 🟡 Medium | UCMuon energy loss at 2.65 MWE is −12% vs Geant4 (0.572 vs 0.649 GeV). Investigate threshold/cutoff effect at shallow depths in PUMAS. | UCMuon | Open |
| 🟡 Medium | Add Geant4 entry position to output CSV to enable lateral-displacement comparison in fig04 and fig15. | Geant4 | Open |
| 🟡 Medium | Harmonise rock I value: Geant4 Z=11 → I=149 eV (NIST Na); PHITS Mg/O → I≈138 eV. Consider ICRU-37 standard rock value (I=136 eV) in both. | Geant4/PHITS | Open |
| 🟢 Low | BB and UCMuon MCS 40–70% below Geant4. Document MCS model (Highland/Molière), X₀ value, and Z/A used for standard rock in file headers. Do not tune. | BB, UCMuon | Open |
| 🟢 Low | MUSIC/PROPOSAL/UCMuon energy loss ~4% above Geant4 at 530 MWE. Document which brem/pair cross-section parametrisation each uses. | MUSIC, PROPOSAL, UCMuon | Open |
| ✅ Done | PROPOSAL MCS enabled — scattering now recorded. cx/cy/cz reflect deflected exit direction. | PROPOSAL | v2 |
| ✅ Done | UCMuon lateral displacement fixed — stochastic position kicks removed. Now consistent with BB (25 vs 24 cm at 530 MWE). | UCMuon | v2 |
| ✅ Done | Fixed `muon_rock.inp`: emin(6)=µ+, emin(7)=µ− added, lateral map tallies changed to `f-curr`. | PHITS | v1 |

---

## 7. Summary Table (v2)

| Quantity | ≤26.5 MWE | 66–530 MWE | Notes |
|----------|-----------|------------|-------|
| Transmission | 0.90–0.95× all codes | ±2% all codes | 26.5 MWE: 5 GeV boundary — exclude from comparison |
| Mean exit KE (BB) | ≤0.5% | ≤2% | Best energy agreement |
| Mean exit KE (UCMuon) | +0.2% | ≤3.5% | Good; fixed from v1 (+1.1% to −3.5% at 530) |
| Mean exit KE (MUSIC, PROPOSAL) | ≤0.5% | ≤4.3% | Acceptable; shared dE/dx tables suspected |
| Mean exit KE (PHITS) | +0.5% | **−6 to −12%** | Requires investigation (radiative dE/dx) |
| MCS angle (MUSIC) | +2–3% | +7–10% | Best fast-code match |
| MCS angle (PROPOSAL) | +2% | **+15–19%** | Overestimates — check X₀ in medium definition |
| MCS angle (BB, UCMuon) | −39 to −77% | −31 to −44% | Consistent underestimate; model limitation |
| Lateral drift (MUSIC/PROPOSAL) | agree | 41–45 cm at 530 MWE | Consistent with their larger MCS angles |
| Lateral drift (BB/UCMuon) | agree | 24–26 cm at 530 MWE | Consistent with their smaller MCS angles |
| Speed (vs Geant4 18878 s) | — | — | PHITS ~371 000 s; grid codes <1 s (no timing files yet) |
