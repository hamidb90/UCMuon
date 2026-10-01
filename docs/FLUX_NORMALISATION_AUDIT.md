# Flux normalisation and angular dependence: audit of the surface muon spectra

2026-09-30, against v1.1.2 (private `main` = 439f18d). Requested after
absolute open-sky rates from UCMuGen were found to disagree by up to x100
between spectra.

**Changes results** (details in Section 5): Reyna in UCMuGen; the fast flux
estimator and everything built on it (Density tab, exposure panel); the
backward MC and the terrain engine (T_sim maps); the GUI rate and live time;
PARMA in multithreaded UCMuGen. Legacy generator events are bit-identical.

Part I (Sections 1-11) audits the spectra. Part II (Sections 12-17) covers
the second round: how the generator samples momentum and angle, the new joint
angular mode 6, one live-time formula across the code, and the detector
filter.

## 1. Scope and method

Every surface spectrum (CosmoALEPH, power law, PARMA, Guan, Frosin,
Gaisser/"Bugaev", Reyna, electrons) was checked in every code path that turns it
into a number:

| path | file | what it produces |
|---|---|---|
| UCMuGen | `ucmugen/include/UCMuGen.h` (`flux::`, `MomentumCdf`, `Generator`) | events, `rate()`, `liveTime()` |
| Fortran generator | `src/generator/ucmuon_source_module.f90` | events, the "Integrated flux" line the GUI reads |
| GUI rate | `gui/ucmuon_gui.py` (`_compute_flux`, exposure-time panel) | rates and live times shown in the GUI |
| fast estimator | `gui/fast_flux_estimator.py` | open-sky flux, transmission, exposure time, density inversion |
| backward MC | `gui/ucmuon_backward_mc.py` | flux and rate at depth |
| terrain engine | `gui/ucmuon_terrain_driver.py` (through the backward MC) | flux and T_sim maps |
| guaranteed-hit | `gui/ucmuon_dasrem_driver.py` | events (momentum shape only) |
| MPI generator | `src/generator/ucmuon_gen.f90` | events, "Total tried" |

Method: code comments and docstrings were not trusted. Each formula and
constant was read against the source paper (PDFs kept in the development
tree's `references/source/`, which is not published). Each
spectrum was then re-implemented independently from the paper, and compared
numerically with every code path and with PDG and published data. The
comparisons are now permanent tests (`tests/flux/`, Section 8).

Equation and table numbers used throughout:

| paper | file | used |
|---|---|---|
| Reyna 2006, arXiv:hep-ph/0604145 | `Reyna-B.pdf` | Eq. 1 (ζ = p cos θ), Eq. 2 (I = cos³θ I_V(ζ)), Eq. 3 (I_V), Sec. 4 (c1..c5, validity), Fig. 3 (data) |
| Guan et al. 2015, arXiv:1509.06176 | `Guan.pdf` | Eq. 1 (Gaisser), Eq. 2 (cos θ\*), Table 1 (P1..P5), Eq. 3, Sec. 1 (Gaisser validity conditions) |
| Frosin et al. 2025, J. Phys. G 52, 035002 | `Frosin.pdf` | Eqs. 5-6 (Guan form), Eqs. 7-8 (Bugaev, Reyna), Table 1 (7 datasets), Sec. 3.2 (fit range), Table 4 (a, b), Sec. 3.4.1 |
| Schmelling et al. 2013, Astropart. Phys. 49, 1 | `cosmoALEPH.pdf` | Table 1 (vertical spectrum, charge ratio, 11.2% normalisation) |
| Bugaev et al. 1998, Phys. Rev. D 58, 054001 | `Bugaev/Bugaev.pdf` | Eq. 3.4, Table II |
| Su et al. 2021, Front. Energy Res. 9, 750159 | `Bugaev/Bugaev-sui.pdf` | Eq. 5 (Reyna form, independent restatement) |
| Lesparre et al. 2010, GJI 183, 1348 | `Bugaev/Lesparre.pdf` | Eq. 11 (Reyna form, independent restatement) |
| PDG 2022, Review of Particle Physics, "Cosmic Rays" (Beatty, Matthews, Wakely) | `rpp2022-rev-cosmic-rays.pdf` | Sec. 30.3.1: I_v(>1 GeV/c) ≈ 70 m⁻²s⁻¹sr⁻¹ "with recent measurements favoring a lower normalization by 10-15%", ≈ 1 cm⁻²min⁻¹ horizontal, ∝ cos²θ at E ~ 3 GeV; Eq. 30.4 (Gaisser formula, valid E > 100/cos θ GeV, θ < 70°). The 2024 edition (`rpp2024-rev-cosmic-rays.pdf`) no longer has this section |
| Tang et al. 2006, Phys. Rev. D 74, 053007 | `PhysRevD.74.053007.pdf` | Sec. II.B, Eqs. 3-10 (modified Gaisser), Fig. 1 |

## 2. Findings at a glance

| # | finding | your hypothesis | status |
|---|---|---|---|
| 1 | UCMuGen `flux::reyna` has no cos³θ and uses cos θ\* instead of cos θ | **confirmed**, both parts | fixed |
| 2 | Fortran `reyna_flux` has the same form | **confirmed**, but it is only ever evaluated at cos θ = 1, so no Fortran output was affected. The docstring's 7.0e-3 claim is **true** | fixed for consistency; bit-identical |
| 2b | "the generated shape may be fine" | **partly wrong**: the Fortran samples p and θ independently for every spectrum, so the hardening with zenith angle is lost | new angular mode 6 (Section 13) |
| 3 | fast estimator Reyna formula is not Reyna's | **confirmed**. Its "~6.2e-3, 12% from PDG" claim is reproducible (6.18e-3), but it holds only through compensating errors | fixed |
| 3b | fast estimator Gaisser-formula models (`bugaev`, `gaisser_tang`) | new: prefactor 1.4e-2 instead of 0.14, **10x low** | fixed |
| 4 | CosmoALEPH extrapolated below its range; vertical-only | **confirmed**. Constants verified as a fit to Table 1. UCMuGen treats it as isotropic, which is defensible in its range; the Fortran/GUI path assumes cos²θ, which is wrong there | guardrails added; angular law unchanged |
| 5 | Guan constants, cos θ\*, Jacobian | **confirmed correct** in every path | no change |
| 6 | CITATION.cff titles | **confirmed**; also Frosin title and DOI, PARMA title, Highland year and title | fixed |
| 7 | backward MC spectra | new: none match their paper; x59 to x10⁸ off, and the terrain engine uses them | fixed |
| 8 | GUI exposure-time panel passes E where T is expected | new: 4% at a 1 GeV threshold, 1.5% at 10 GeV | fixed |
| 9 | GUI Results-tab rate formula | new: a fixed cos²θ law, no projection, and it grows with the source area in detector mode; up to x8 | replaced by T = N_tried / R (Section 14) |
| 10 | manuscript | new: Reyna equation, Frosin and Schmelling bib entries, stale statements | fixed; notes kept with the manuscript |
| 11 | angular mode 4 samples J(p,θ) jointly? | **no**: θ given E is right, E is weighted by the vertical spectrum; +7% to −30% | Section 12 |
| 12 | PARMA in the Fortran generator | new: total energy passed where PARMA expects kinetic (≤4% in shape) | fixed |
| 13 | PARMA in multithreaded UCMuGen / OpenMP | new: data race in PARMA's memoisation; 190 633 of 192 000 evaluations wrong with 8 threads | fixed (Section 16) |
| 14 | MPI generator angular modes | new: accepted 1-4 only, so mode 5 silently became cos²θ | fixed |
| 15 | detector filter intersection | correct for every source outside the inflated detector; caps misplaced for origins inside the margin zone | caps moved (Section 15) |
| 16 | "Total tried" under MPI + OpenMP | correct: no lost trials, at most N_threads−1 extra per rank | measured (Section 15) |
| 17 | GUI tried-count parser | new: took the first "tried N" (a progress line), not the total | fixed |

## 3. Spectrum by spectrum

Notation: I(p, θ) is the differential intensity per unit solid angle,
cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹; dN/dE and dN/dp are related by dN/dp = (p/E) dN/dE,
and dN/dT = dN/dE since T = E − m.

### 3.1 Reyna 2006 (spectrum 7)

**Paper.** Eq. 1: ζ = p cos θ. Eq. 2: I(p, θ) = cos³θ · I_V(ζ).
Eq. 3: I_V(p) = c1 p^−(c2 + c3 y + c4 y² + c5 y³), y = log10 p, in
cm⁻²s⁻¹sr⁻¹GeV⁻¹ (Sec. 1 and figure axes), with the Sec. 4 "Best Fit"
(c1..c5) = (0.00253, 0.2455, 1.288, −0.2555, 0.0209). Stated validity (Sec. 4):
p > 1 GeV/c and p < 2000 GeV/c / cos θ. The paper has three equations; the
"Eq. 6-7" cited in the code comments does not exist in it. Frosin Eq. 8, Su
Eq. 5 and Lesparre Eq. 11 restate the same form with plain cos θ.

| path | before | after |
|---|---|---|
| UCMuGen `flux::reyna` | I_V(p cos θ\*): no cos³θ, Guan cos θ\* | cos³θ · I_V(p cos θ) |
| Fortran `reyna_flux` | same as UCMuGen | same fix; only called at cos θ = 1, where cos θ\*(1) = 1.00000000033 is clamped to exactly 1, so the CDF is **bit-identical** |
| fast estimator `_reyna_bugaev` | log10(p³I) = 0.00253x³ − 0.2455x² + 1.288x − 4.25, ×100, × cos θ\*^1.85 with cos θ\* = √((c²+0.102²)/(1+0.102²)) | Eqs. 1-3, then × E/p for dN/dT |
| backward MC | no Reyna mode (its "Guan" used Reyna's constants, see 3.4) | n/a |

Constants c1..c5 are transcribed correctly in UCMuGen and the Fortran. The
fast estimator's polynomial is the same five numbers reshuffled, with sign
changes and an invented −4.25 and ×100. It is not algebraically equivalent to
Eq. 3.

Units and Jacobian: Eq. 3 is per GeV/c, so UCMuGen and the Fortran need no
Jacobian; the fast estimator returns dN/dT and multiplies by E/p. Correct in
all three after the fix.

Why the missing cos³θ mattered: I_V is steep, so I_V(p cos θ) > I_V(p) and the
intensity rose with zenith angle. The horizontal-surface rate was x4 (E > 1
GeV) to x12 (E > 100 GeV) too high.

The old docstring's validation claims:
- Fortran/UCMuGen "7.0e-3 cm⁻²s⁻¹sr⁻¹ above 1 GeV at cos θ = 1": **true**,
  7.02e-3 for p > 1 GeV/c (Section 4).
- fast estimator "~6.2e-3 above 1 GeV, 12% from PDG": reproducible (6.18e-3),
  but the formula was x0.48 of the data at 10 GeV/c, x1.71 at 1 TeV/c, and
  x0.24 at 10 GeV/c and 75°. The integral was right by accident.

### 3.2 Guan 2015 (spectrum 4)

**Paper.** Eq. 3:
dI/dE = 0.14 (E/GeV)^−2.7 [1 + 3.64 GeV/(E (cos θ\*)^1.29)]^−2.7
× [1/(1 + 1.1 E cos θ\*/115 GeV) + 0.054/(1 + 1.1 E cos θ\*/850 GeV)],
in cm⁻²s⁻¹sr⁻¹GeV⁻¹ (units under Eq. 1). Eq. 2:
cos θ\* = √[(cos²θ + P1² + P2 cos^P3 θ + P4 cos^P5 θ)/(1 + P1² + P2 + P4)],
with Table 1 (P1..P5) = (0.102573, −0.068287, 0.958633, 0.0407253, 0.817285).

| check | result |
|---|---|
| P1..P5 in all three implementations | match Table 1 |
| `GUAN_DENOM` = 0.99144315 | = √(1 + P1² + P2 + P4) = 0.991443 |
| E_eff = E(1 + a/(E cs^b)), then E_eff^−2.7 | equals E^−2.7 (1 + a/(E cs^b))^−2.7 |
| 0.14, 115, 850, 0.054, 1.1, −2.7, 3.64, 1.29 | match Eq. 3 |
| dN/dE → dN/dp | UCMuGen and Fortran multiply by p/E; the fast estimator returns dN/dT (no factor needed); the backward MC now uses the fast estimator |
| angular dependence | energy-dependent through cos θ\* in both the low-energy term and the meson terms |
| validity | Guan states no lower bound; Frosin Sec. 3.2 fits this form to data at 1 GeV-1 TeV and 0-90°. Warning below E = 1 GeV |

**Verdict: correct in UCMuGen, Fortran and fast estimator. No change.**
The GUI's claim that Guan "suppresses the surface flux 50-100x below 20 GeV"
was false: it compared Guan with the plain Gaisser formula as if that were the
truth. Guan is 9% below Reyna in the vertical at 1 GeV/c (Frosin 2%) and 11% below
the data at 10 GeV/c. The backward MC's "Guan" (mode 3) was an unrelated formula
(0.00253 E^−3.7 ..., Reyna's constants used as m⁻²); it now calls the real one.

### 3.3 Frosin 2025 (spectrum 5)

**Paper.** Eqs. 5-6 are the Guan form. Table 4 gives the refit
a = 3.512 ± 0.012, b = 1.388 ± 0.005. The fit used about 300 points from the 7
datasets of Table 1, at 1 GeV-1 TeV and 0-90° (Sec. 3.2). Frosin Table 6 also
refits Reyna's coefficients; those are not used anywhere in UCMuon.

**Verdict: constants correct in all paths; no change.** The GUI's "re-fitted on
304 datasets" is corrected to about 300 points from seven datasets.

### 3.4 Gaisser 1990 (spectrum 6, "Bugaev/Gaisser"; fast-estimator keys `bugaev`, `gaisser_tang`)

**Paper.** Guan Eq. 1 (= Frosin Eq. 2 = PDG):
dI/dE = 0.14 E^−2.7 [1/(1 + 1.1 E cos θ/115) + 0.054/(1 + 1.1 E cos θ/850)].
Guan Sec. 1 gives its validity: θ < 70° and E > 100/cos θ GeV (muon decay and
energy loss neglected).

| path | implementation | verdict |
|---|---|---|
| UCMuGen, Fortran (spectrum 6) | Guan Eq. 3 with a = 0, b = 1, i.e. Gaisser with cos θ\* for cos θ | a correct Gaisser formula with Guan's curvature fix. It contains **nothing from Bugaev 1998** despite its name |
| fast estimator `bugaev` | Gaisser with plain cos θ and **A = 1.4e-2** | 10x low. Now = spectrum 6 exactly |
| fast estimator `gaisser_tang` | Gaisser with plain cos θ, A = 1.4e-2, times (1 + 0.054E/800) | not Tang's formula: the factor is in no paper, and "Gaisser & Tang (1984)" does not exist. **Replaced** by Tang et al. 2006 Eqs. 3-10 (Section 3.8) |

Below 100 GeV this formula overestimates badly: 12x the PDG vertical integral
above 1 GeV, x2 at 10 GeV/c, x12 at 10 GeV/c and 75°. Inside its range it
agrees with CosmoALEPH Table 1 to +11% at 112 GeV/c and −28% at 1122 GeV/c
(the same as Guan). The old 1.4e-2 made its integral above 1 GeV look right
(84 vs 70) by being 10x low everywhere; the GUI's warning "10x too low below 10
GeV" described that bug as physics. Warnings now fire below E = 100/cos θ GeV.

### 3.5 CosmoALEPH (spectrum 1, default)

**Paper.** Schmelling et al. 2013 give Table 1: the vertical sea-level spectrum
in 14 bins centred at 112-2239 GeV/c, in (GeV/c s cm² sr)⁻¹, plus an 11.2%
global normalisation uncertainty. The paper gives **no fit formula**.

**Code.** dN/dp = 10^3.8467 p^−3.1952 m⁻²s⁻¹sr⁻¹(GeV/c)⁻¹, × 1e-4 for cm⁻².
Verified: an unweighted least-squares fit of log10 I against log10 p to all 14
bins of Table 1 (in m⁻²) returns exactly A = 3.8467, B = −3.1952. It agrees
with Table 1 to within 6% for the first 11 bins; the 1778 GeV/c bin is 31% off
(large errors there). The 1e-4 is correct: the fit gives 1.99e-7 cm⁻² at
112 GeV/c against 1.959e-7 in Table 1. The charge-ratio table matches Table 1.

**Validity.** Extrapolated below ~100 GeV/c the power law is far too high:
x3.4 against the vertical data at 10 GeV/c (not x2.5 as the GUI and paper
said), and x46 in the vertical integral above 1 GeV/c (3200 against 70
m⁻²s⁻¹sr⁻¹). That, and not a code bug, is why CosmoALEPH extrapolated to
E > 1 GeV gives 60 cm⁻²min⁻¹ through a horizontal surface (isotropic, as in
UCMuGen, over the whole upper hemisphere; 30 with cos²θ), against PDG's
I ≈ 1 cm⁻²min⁻¹ for horizontal detectors (PDG 2022 Sec. 30.3.1).

**Angular dependence.** A vertical measurement carries none, and the paths
disagree:

| path | assumption |
|---|---|
| UCMuGen | isotropic: I(p, θ) = I_V(p) |
| Fortran | whatever angular mode the user picks (cos²θ by default in the GUI) |
| GUI rate | cos²θ (Ω_cos²) |
| backward MC | cos²θ before; isotropic now, as UCMuGen |

At the momenta where the fit is valid, the angle-dependent models put
I(θ)/I(0) at: Guan 1.14 and Reyna 0.94 at 60° and 100 GeV/c; Guan 2.58 and
Reyna 1.42 at 75° and 1 TeV/c. Isotropic sits inside that spread. Two
alternatives were rejected. Reyna's scaling applied to the power law gives
cos^−0.195 θ, which diverges at the horizon and evaluates the fit below
100 GeV/c. cos²θ gives 0.07 at 75°, x12-38 below both models. **UCMuGen is
left isotropic.** The cos²θ that the Fortran/GUI path and the paper's Table 3
attach to Mode 1 is the wrong choice for thick targets, which is the use
Mode 1 is recommended for.

**Guardrail:** all three implementations warn when p_min < 99 GeV/c (or
p_max > 2.5 TeV/c). CosmoALEPH is not refused, because the GUI default range
(100-2500 GeV) is already valid and refusing would break existing configs.

### 3.6 PARMA/EXPACS (spectrum 3)

The PARMA functions are JAEA's code and were not re-derived; the C++ port
agrees with the Fortran to seven digits for muons (checked again here). What
was checked is the unit
handling in `UCMuGen_PARMA.h::intensity`: μ⁺ + μ⁻ spectrum per MeV, times the
per-sr angular factor, × 1e3 × p/E, which gives cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹.
Numerically (sea level, 3 GV cutoff, W = 0): vertical integral 63.0
m⁻²s⁻¹sr⁻¹ above 1 GeV/c; horizontal flux 0.61 cm⁻²min⁻¹ above 1 GeV and
0.79 cm⁻²min⁻¹ above 0.12 GeV. All checks in Section 4 pass. **Correct.** The
Fortran PARMA path passed total energy where PARMA takes kinetic energy, and
PARMA's routines are not thread-safe; both are fixed (Part II, Section 16).
The generator now prints a surface rate for PARMA too.

### 3.7 Power law E^−3.7 (spectrum 2) and cosmic electrons (spectrum 8)

Both are sampling shapes. UCMuGen returns p^−3.7 and p^−3; the Fortran samples
the same laws in p. The comments say dN/dE; for muons above 1 GeV and
electrons above 10 MeV, E/p differs from 1 by less than 0.6% and 0.2%.
Neither has an absolute normalisation, and neither path pretends otherwise:
UCMuGen `rate()` returns −1, the Fortran prints a "sampling-shape norm" rather
than "Integrated flux", and the GUI shows no rate. Messages saying so now
appear in all three. The backward MC's power law (mode 2) carries an arbitrary
95 m⁻²s⁻¹sr⁻¹GeV⁻¹ prefactor; its rates are now documented as meaningless.

### 3.8 Tang et al. 2006 (fast-estimator key `gaisser_tang`)

**Paper.** Sec. II.B, Eq. 3:
dN/dE dΩ = A · 0.14 E^−γ [1/(1 + 1.1 Ẽ cos θ*/115) + 0.054/(1 + 1.1 Ẽ cos θ*/850) + r_c]
(the printed "E^γ" is a sign slip; γ = 2.70), with cos θ* from Eq. 10 (Guan's
parametrisation, same P1..P5) and three segments:

| segment | A | Ẽ | r_c |
|---|---|---|---|
| E > 100/cos θ* | 1 | E | 0 |
| 1/cos θ* < E ≤ 100/cos θ* | 1.1 (90 √(cos θ + 0.001)/1030)^(4.5/(Ẽ cos θ*)) (Eq. 7) | E + Δ, Δ = 2.06·10⁻³ (950/cos θ* − 90) GeV (Eqs. 5-6) | 10⁻⁴ (Eq. 4) |
| E ≤ 1/cos θ* | as above, after E → (3E + 7 sec θ*)/10 (Eq. 9) | | |

**A misprint in Eq. 7.** The printed exponent is 4.5/(E cos θ*), with the
surface energy E. Evaluated that way the formula cannot be what Tang et al.
fitted: at θ = 0° and 1 GeV it gives E^2.7 dN/dE dΩ = 2.7·10⁻⁶ where their
Fig. 1 curve is ≈ 2·10⁻³, at 60° 1.6·10⁻⁸ against ≈ 4·10⁻⁴, at 87° and
10 GeV 3·10⁻⁸ against ≈ 3·10⁻⁴, and the vertical intensity above 1 GeV/c is
22 m⁻²s⁻¹sr⁻¹ (PDG ≈ 70). With Ẽ = E + Δ in the exponent (Δ is "the muon
energy loss in the atmosphere", Ẽ "the muon energy on top of the
atmosphere", so the decay factor naturally takes Ẽ) it reproduces Fig. 1 at
every angle shown: 3.0·10⁻³, 3.8·10⁻⁴, 5.3·10⁻⁴ at the same three points, and
0.2 GeV / 0° gives 5.7·10⁻⁵ against the data point at ≈ 5·10⁻⁵. The
implementation uses Ẽ. (Su et al. 2021, Eq. 3, restate the printed form, so
implementations that follow the printed equation are 10-10⁵× low below
~10 GeV.)

**Result.** Vertical intensity above 1 GeV/c 60.3 m⁻²s⁻¹sr⁻¹; all
Section 4 checks pass (0.86 of PDG, 0.88 of the 10 GeV/c data, 0.62-0.84 of
the OKAYAMA zenith data at 10 GeV/c); four points read off Fig. 1 agree
within a factor 2 (locked in `test_flux_reference.py`). Validity: the paper
fits 0.1 GeV-10 TeV at 0-87° and reports up to 40% disagreement at
E < 10 GeV and θ > 85°; the fast estimator warns there. Tang is not a
generator spectrum.

## 4. Numerical checks against independent references

Evaluated with the fixed code (`tests/flux/test_flux_reference.py`, which
prints this table). Each entry is model divided by reference. Reference
columns:
- PDG: I_v(p > 1 GeV/c) = 70 m⁻²s⁻¹sr⁻¹; horizontal flux (π/2)·70
  m⁻²s⁻¹ for E > 1 GeV, θ < 90°.
- R3: the θ = 0 experiments of Reyna Fig. 3 (Nandi & Sinha, MARS,
  OKAYAMA 0°), 7 points around 10 GeV/c.
- CA: CosmoALEPH Table 1.
- OKAYAMA: Reyna Fig. 3, unscaled by cos³θ, at p = 10 GeV/c.
- cos²θ: PDG's angular law, compared with ∫I(p > 1 GeV/c) dp at θ over the
  same integral at 0°.

Values in parentheses are outside the model's validity and are not tested.
Tang et al. 2006 (fast estimator only) is in `test_flux_reference.py` too: 0.862,
0.934, 0.877, 1.110, 0.723, 0.837, 0.780, 0.624, 1.058, 1.190, 1.274 in the
row order below, all within tolerance.

| check | tol. | Reyna (7) | Guan (4) | Frosin (5) | PARMA (3) | CosmoALEPH (1) | Gaisser (6) |
|---|---|---|---|---|---|---|---|
| I_v(>1 GeV/c) / PDG | ±20% | **1.002** | 0.862 | 0.905 | 0.900 | (45.7) | (12.0) |
| horizontal E > 1 GeV / (π/2)·PDG | ±20% | 1.076 | 0.909 | 0.922 | 0.929 | (92.6) | (24.5) |
| I_v(10 GeV/c) / R3 | ±15% | 0.967 | 0.888 | 0.911 | 0.932 | (3.41) | (2.05) |
| I_v(112 GeV/c) / CA | ±15% | 1.070 | 1.018 | 1.021 | 1.048 | 1.016 | 1.110 |
| I_v(1122 GeV/c) / CA | ±35% | 1.074 | 0.717 | 0.717 | 0.720 | 1.045 | 0.723 |
| I(10 GeV/c, 30°) / OKAYAMA | ±15% R, ±50% | 0.931 | 0.853 | 0.868 | 0.877 | (3.74) | (2.28) |
| I(10 GeV/c, 60°) / OKAYAMA | ±15% R, ±50% | 0.953 | 0.803 | 0.771 | 0.659 | (7.07) | (4.44) |
| I(10 GeV/c, 75°) / OKAYAMA | ±15% R, ±50% | 0.907 | 0.639 | 0.540 | 0.628 | (18.7) | (12.0) |
| ∫I(30°) / ∫I(0°) / cos²θ | ±30% | 1.049 | 1.037 | 1.019 | 1.031 | (1.33) | (1.34) |
| ∫I(60°) / ∫I(0°) / cos²θ | ±30% | 1.180 | 1.118 | 1.014 | 0.938 | (4.00) | (4.04) |
| ∫I(75°) / ∫I(0°) / cos²θ | ±30% | 1.241 | 1.185 | 0.969 | 1.113 | (14.9) | (15.2) |

Tolerance reasoning (also in the test header):
- PDG quotes its values as "≈" with no uncertainty.
- CosmoALEPH carries an 11.2% normalisation uncertainty and about 17% in
  total at 1122 GeV/c; the ±35% there is 2x that.
- Frosin Sec. 3.4.1 reports a factor-2 spread between models at 75°, hence
  ±50% for the non-Reyna models at large angle.
- Reyna was fitted to the OKAYAMA points, hence ±15% for it.

Absolute numbers behind the ratios:

| quantity | Reyna | Guan | Frosin | PARMA | PDG / data |
|---|---|---|---|---|---|
| I_v(p > 1 GeV/c) [m⁻²s⁻¹sr⁻¹] | 70.2 | 60.4 | 63.4 | 63.0 | ≈ 70 |
| horizontal, E > 1 GeV, θ < 90° [cm⁻²min⁻¹] | 0.71 | 0.60 | 0.61 | 0.61 | ≈ 1 for all energies (PARMA: 0.79 above 0.12 GeV) |
| p³I_V at 10 GeV/c [GeV² cm⁻²s⁻¹sr⁻¹] | 0.127 | 0.117 | 0.120 | 0.122 | 0.131 (R3) |
| p³I_V at 112 GeV/c | 0.295 | 0.280 | 0.281 | 0.288 | 0.275 (CA) |
| p³I_V at 1122 GeV/c | 0.183 | 0.122 | 0.122 | 0.123 | 0.171 (CA) |
| I_V(19-1500 GeV) [cm⁻²s⁻¹sr⁻¹] | 3.21e-4 | 3.11e-4 | | | |

The last row answers the fast-estimator question in the request. The old
fast-estimator Reyna gave 1.86e-4, 40% below Guan; it now gives 3.21e-4, 3%
above Guan. (The request compared Reyna's vertical integral with Guan's
*horizontal-surface* 8.45e-4 cm⁻²s⁻¹, which is a different quantity.)

Reading the table:
- Reyna reproduces PDG to 0.2% and every data point within 10%.
- Guan, Frosin, PARMA and Tang sit 9-14% below PDG's ≈ 70 in the vertical,
  which PDG itself qualifies: "recent measurements [favour] a lower
  normalization by 10-15%" (Sec. 30.3.1). They fall off too
  fast at large angle and 10 GeV/c, and all three are 28% below CosmoALEPH at
  1 TeV/c. That is model spread, not a code error.
- CosmoALEPH and plain Gaisser are fine at 112 GeV/c and useless at 1 GeV.

## 5. What changes for users

Rows marked **changes results** alter numbers users may have quoted.

| output | change | size |
|---|---|---|
| UCMuGen `Generator::rate()`, `liveTime()`, spectrum 7 (**changes results**) | Reyna fix | horizontal disk, θ ≤ 85°, E ≤ 1500 GeV: E > 1: 448.6 → 118.1 m⁻²s⁻¹ (×0.263); E > 10: 133.0 → 20.5 (×0.154); E > 19: 66.3 → 8.35 (×0.126); E > 100: 5.46 → 0.437 (×0.080). Live times ×3.8 to ×12.5 longer. `feature_tour` `spectrum_rate.reyna_bugaev` 657.1 → 380.4 (reference updated; no document quotes it) |
| UCMuGen events, spectrum 7 (**changes results**) | Reyna fix | the p-θ distribution now hardens with θ |
| UCMuGen, PARMA in multithreaded use (**changes results**) | thread-safety fix | single-threaded runs unchanged; with several threads, the old values were wrong (Section 16) |
| UCMuGen, all other spectra | none | rates and events bit-identical |
| Fortran generator, legacy angular modes 1-5, spectra ≠ 3 | none | events bit-identical (20/20 legacy configurations); new console lines: validity warnings, "Surface rate R", "Live time" |
| Fortran generator, PARMA (**changes results**, small) | kinetic energy | sampled energy density ×1.04 at 1 GeV, ×1.02 at 10 GeV, ×1.003 at 100 GeV |
| MPI generator, angular mode 5 (**changes results**) | mode accepted | runs asking for mode 5 got cos²θ before |
| GUI Results-tab rate and live time (**changes results**) | new formula (Section 14) | was Φ_v·A·Ω_cos²·N/N_gen: ×1.25-1.3 high for a plain horizontal source, up to ×8 with the detector filter |
| fast estimator, `reyna_bugaev` (GUI default, Density tab default) (**changes results**) | new formula | open-sky ∫ from T = 0.5 GeV: ×0.83 at 0°, ×0.79 at 60°, ×0.66 at 75°. Rock (ρ = 2.65), θ = 0: flux ×1.55 at X = 1000 g/cm², ×1.25 at 26 500 (100 m), ×0.80 at 100 000. θ = 60°: ×1.9, ×4.1, ×3.2. Transmission at 26 500 g/cm²: 0.30% → 0.45% vertical, 0.30% → 1.58% at 60°. **Density-tab opacity and density inversions shift accordingly** |
| fast estimator, `bugaev` (**changes results**) | prefactor | every flux ×10.000; exposure times ÷10; transmissions unchanged |
| fast estimator, `gaisser_tang` (**changes results**) | now Tang et al. 2006 | vertical intensity above 1 GeV/c 84.2 (v1.1.2) → 60.3 m⁻²s⁻¹sr⁻¹; shape changes at every energy below ~100/cos θ* GeV |
| GUI exposure-time panel (**changes results**) | E/T fix and R | now T = N / R with the same R as the generator; the E/T fix alone is +4.0% flux at E_min = 1 GeV, +1.5% at 10 GeV |
| backward MC (**changes results**) | new spectra | rate at depth, E 1-5000 GeV: CosmoALEPH ×59, Guan ×3.2e5-4.0e5, Frosin ×4.3e6-6.3e7 (10 m and 100 m rock). Old absolute values were meaningless |
| terrain engine flux and T_sim maps (**changes results**) | via the backward MC | backward-MC spectrum 3 (Guan, the MURAVES example's): T ×35-1000 (θ 60-85°, X = 1e4-3e5 g/cm²); spectrum 1 (CosmoALEPH, the committed Vesuvius library's): T ×1.48, nearly uniform |
| guaranteed-hit (DAS-REM) generator | shapes | spectrum 7 momentum shape at the batch-mean θ now hardens with θ |
| paper figures (`manuscript/scripts`) | Fig. 3 only | the non-existent cos⁶θ mode removed; Reyna calls are at cos θ = 1 |

## 6. The rate and live-time chain

What the code did before this audit: UCMuGen's `rateAndError` was correct
(full J(p,θ), projection). The Fortran printed only the vertical integral
Φ_v. The GUI Results tab computed R = N_rows · Φ_v · A_src · Ω / N_gen, with
Ω = 2π(1 − cos³θ_max)/3 without the filter and the cos²-weighted solid angle of
the detector seen from the source centre with it: a fixed cos²θ sky, no
projection, and (with the filter) growth with the source area. Errors:

| configuration | old GUI / correct |
|---|---|
| no filter, horizontal disk, θ ≤ 85°, E > 1 GeV: Reyna, Guan, Frosin | 1.25, 1.27, 1.32 |
| no filter, CosmoALEPH E > 100 GeV | 0.67 |
| filter, 1 m² detector 10 m below a disk of radius 5 / 10 / 20 m | 1.39 / 2.67 / 8.33 |
| Guan + mode 4, by the energy of the muons that matter | 1.61 at 1 GeV, 1.05 at 10, 0.63 at 100, 0.45 at 1000 GeV |

The exposure-time panel used a third formula (projection-correct for a cos²θ
sky, 4/3 off the Results tab). All three are replaced by one (Section 14).

## 7. Factorised momentum-angle sampling

Modes 1-3 and 5 draw p from the vertical spectrum and θ independently; mode 4
draws θ from P(θ|E) but E from the vertical spectrum. Real spectra harden
with zenith angle (Reyna at 100 GeV/c: I(p,θ)/I(p,0) = 1.00, 0.94, 0.78 at
30°, 60°, 75°, where cos²θ gives 0.75, 0.25, 0.07). Effective intensity of
spectrum 7 delivered to a detector over the true Reyna intensity:

| p [GeV/c] | 0° | 30° | 60° | 75° |
|---|---|---|---|---|
| cos³θ mode: 3 | 1.33 | 1.29 | 1.33 | 2.07 |
| cos³θ mode: 10 | 1.33 | 1.14 | 0.70 | 0.52 |
| cos³θ mode: 100 | 1.33 | 1.00 | 0.35 | 0.11 |
| cos³θ mode: 1000 | 1.33 | 0.95 | 0.27 | 0.06 |
| cos²θ mode: 100 | 1.00 | 0.87 | 0.53 | 0.33 |

The legacy modes are kept bit-exact; angular mode 6 (Section 13) samples the
joint density, and the GUI weights legacy events (Section 14).

## 8. Tests

| test | checks | result |
|---|---|---|
| `tests/flux/test_flux_consistency.py` | I(p,θ), vertical flux, horizontal rate identical in C++, Fortran, Python (fails above 1%); UCMuGen and Python warn alike | pointwise ≤ 1e-9; integrals < 0.1%; `rateAndError` within 0.1-0.23% of quadrature |
| `tests/flux/test_flux_reference.py` | each spectrum against PDG, CosmoALEPH Table 1, Reyna Fig. 3 data (Section 4) | all within the stated tolerances |
| `tests/flux/test_joint_sampling.py` | Fortran mode 6 against `ucmugen::Generator` in (E, θ) bins, 6 cases; legacy mode 4 as a control that must fail | all agree; mode 4 fails at 22σ (Section 13) |
| `tests/flux/test_live_time.py` | GUI live time against UCMuGen, plane/disk/hemisphere × filter on/off; weighted legacy modes | all within 1% (Section 14) |
| `tests/geometry/test_ray_intersections.f90` | detector-filter intersections against brute force, 399k rays | 0 mismatches (Section 15) |
| `ucmugen/validation/test_parma.cc` (new check) | 8 threads reproduce single-thread PARMA intensities exactly | 0 of 192 000 differ (old header: 190 633) |
| existing: `run_tests.sh`, `compare_legacy.py` | UCMuGen suite, bit-exact legacy streams | pass; 20/20 identical |

All but the MPI check run in CI (`.github/workflows/ucmugen.yml`, job `flux`).

## 9. Citations and the manuscript

**CITATION.cff**, corrected entries:

| entry | was | now |
|---|---|---|
| Kudryavtsev 2009 (MUSIC) | "Muon radiography of volcanoes: A review" | "Muon simulation codes MUSIC and MUSUN for underground physics" (from `MUSIC.pdf`) |
| Guan 2015 | Reyna's title, one author | "A parametrization of the cosmic-ray muon flux at sea-level", five authors |
| Frosin 2025 | "Muon flux parametrization for the muon depth-intensity relation", doi …/ada22d, start 35002 | "Optimizing sea level muon flux modelling: a 2D fit approach and integration into Geant4 generators", doi 10.1088/1361-6471/adb6c3, article 035002, seven authors (from `Frosin.pdf`) |
| Sato 2015 (PARMA) | "PARMA/EXPACS: A revision of PARMA" | "Analytical model for estimating terrestrial cosmic ray fluxes nearly anytime and anywhere in the world: Extension of PARMA/EXPACS", e0144679 |
| Highland | "Multiple Coulomb scattering", 1979 | "Some practical remarks on multiple scattering", 1975, pp. 497-499 (the DOI was already the 1975 paper's) |
| Lüscher 1994 | note "RANLUX RNG" | note says it is used by MUSIC and the serial transport; the generator's `par_ranlux` is an LCG |
| Groom 2001, Koehne 2013 | correct | end pages added |
| added | | Reyna 2006; Schmelling 2013 (doi 10.1016/j.astropartphys.2013.07.008, from the PDF footer); Bugaev 1998; Sato 2016 (PARMA zenith dependence) |

**.zenodo.json** carries no references; its title matches CITATION.cff.
Nothing to change.

**README**:
- the spectrum list said "Bugaev/Gaisser 1998" (it is Gaisser 1990);
- the references table lacked CosmoALEPH, Reyna, Bugaev and Gaisser;
- the BibTeX block had an old title and the placeholder DOI
  `10.5281/zenodo.XXXXXXX` (now the concept DOI 10.5281/zenodo.20826984).

All fixed, and a "Surface spectra: validity and absolute normalisation"
section was added.


Corrections to the CPC manuscript (not in this repository) were made
alongside: the Reyna equation, the Frosin and Schmelling bibliography entries,
and statements this audit showed to be wrong.

## 10. Open items and missing references

1. Done: PDG values now cited from PDG 2022 Sec. 30.3.1 / Eq. 30.4; the
   fast estimator's `gaisser_tang` replaced by Tang et al. 2006 (Section 3.8).
2. **Regenerate** what Section 17 lists (terrain T_sim libraries, MURAVES
   maps, Density-tab inversions made with the fast estimator).
3. The fast estimator's altitude factor exp(h/8500 m) has no reference.
4. The charge-ratio table applies each CosmoALEPH bin value up to the bin
   centre (a half-bin shift; negligible for rates).
5. The MPI generator writes surface rows for up to N_threads−1 attempts per
   rank after its target, when the filter is off (the OpenMP version gates
   them). Harmless for rates (the tried count is the reduced total), noted
   for completeness.
6. `tools/check_consistency.py` reports the public tree out of step until the
   mirror (prepared, not pushed) is committed.

## 11. Reproduce

```
python3 tests/flux/test_flux_consistency.py
python3 tests/flux/test_flux_reference.py
make ucmuon_gen_omp
python3 tests/flux/test_joint_sampling.py      # --quick for CI size
python3 tests/flux/test_live_time.py           # --quick for CI size
bash tests/geometry/run.sh
bash ucmugen/validation/run_tests.sh
(cd ucmugen/validation && python3 compare_legacy.py --n 20000)
```

---

# Part II: sampling, live time and the detector filter

## 12. Does angular mode 4 sample J(p, θ) jointly?

**No.** In both generators (`ucmuon_gen_omp`, `ucmuon_gen`), mode 4 draws p
from the vertical momentum CDF (`sample_momentum`, built from J(p, 0)) and
then θ from P(cos θ | E) ∝ J(E, θ) cos θ (`sample_guan_angle`). The sampled
density is therefore

    s(p, θ) = [J(p,0)/Φ_v] · J(p,θ) cos θ / N(p),   N(p) = ∫ J(p,θ) cos θ dΩ,

against the true through-surface density J(p,θ) cos θ / R. Their ratio,
J(p,0) R / (Φ_v N(p)), depends on p only: the zenith distribution at a given
energy is right, the energy weighting is not, because the angular integral
N(p) grows with energy (the sec θ enhancement) and the vertical spectrum
ignores that.

Measured for the production configuration (Guan, mode 4, 15-1500 GeV,
θ ≤ 85°, horizontal disk R = 500 m; 10⁶ events from `bin/ucmuon_gen_omp`,
against the exact quadrature of J(p,θ) cos θ). Sampled over true, with equal
totals:

| E \ θ | 0-20° | 20-40° | 40-60° | 60-85° |
|---|---|---|---|---|
| 15-30 GeV | 1.073 | 1.074 | 1.068 | 1.063 |
| 30-100 GeV | 0.897 | 0.892 | 0.894 | 0.891 |
| 100-1500 GeV | 0.685 | 0.721 | 0.697 | 0.702 |

Statistical errors 0.2-1.3%. With a cylinder filter (vertical cylinder,
r = 1.5 m, 20-24 m deep, margin 0.5 m, R = 40 m disk; against UCMuGen) the
pattern is the same within errors: 1.06-1.08, 0.89-0.91, 0.69-0.75. The
filter selects directions and positions, not energies, so it cannot change a
bias that depends on energy alone.

What it means for a run like the production source: muons above 100 GeV are
under-represented by 30%, 30-100 GeV by 11%, and 15-30 GeV over-represented by
7%, at every zenith angle. Deep detectors see mostly the high-energy part, so
their hit rates from such a sample are low by up to 30% before any rate
normalisation. Since v1.2.0 the GUI weights mode-4 events exactly (Section 14);
a new run with mode 6 needs no weights.

**§7 applies to all legacy modes**: 1-3 and 5 fully (θ independent of E,
wrong in both θ|E and E), 4 partly (E only), and the PARMA path in every
legacy mode (energy-averaged angular distribution).

## 13. Angular mode 6: joint J(p, θ) × projection

`generate_muon_joint` (Fortran module) and `parma_generate_joint` (PARMA, in
both main programs) draw (p, position, direction) from

    J(p, θ) · max(0, −n·d)

by the accept-reject of `ucmugen::Generator`: p from the vertical CDF,
position uniform on the surface, direction uniform in solid angle within
θ ≤ θ_max about the **world** zenith, accepted with probability
[J(p,θ)/J(p,0)] · projection / envelope. The envelope is scanned on a 128×128
(p, cos θ) grid with a 15% margin, as in UCMuGen. The sky no longer rotates
with a vertical source plane, and a tilted surface projects onto its own
normal. For PARMA, the angular factor is read from a table built serially
(Section 16). Legacy modes 1-5 are untouched and stay bit-exact. Modes 5 and
6 are now also accepted by the MPI generator.

Validation, Fortran mode 6 against UCMuGen (`test_joint_sampling.py`, full
size), fraction ratios in the requested bins:

| case | max \|pull\| | χ²/ndf |
|---|---|---|
| 1. Guan, disk R = 500 m, 15-1500 GeV, θ ≤ 85° (production) | 1.9 | 12.1/11 |
| 1b. same, legacy mode 4 (control, must fail) | 22.5 | 3649/11 |
| 2. Guan, R = 40 m, cylinder filter 20-24 m deep | 1.4 | 7.0/11 |
| 3. Reyna, hemisphere R = 20 m | 2.3 | 16.2/11 |
| 4. Frosin, vertical XZ rectangle (plus d_y, d_z components) | 1.6 | 9.3/11 |
| 5. CosmoALEPH, disk tilted 30° (plus d_x, d_y components) | 2.2 (3.0 in d_y) | 19.1/11 |
| 6. PARMA, Louvain sea level, disk | 2.7 | 17.2/11 |

Production case, Fortran mode 6 over the exact quadrature (10⁶ events):
0.989-1.004 in every bin, all within 1σ. Under MPI + OpenMP (3 ranks × 2
threads, cylinder filter), the distribution agrees (max pull 2.2).

Mode 6 is the default in the GUI and recommended in the docs.

## 14. One live-time formula

    T = N_tried / R,   R = ∫dp ∫dΩ J(p, θ) ∫_S max(0, −n·d) dA,

R being the rate of muons crossing the generation surface inside the energy
and zenith windows. It is computed identically by UCMuGen's `rate()`
(Monte Carlo), by the Fortran `surface_rate()` (quadrature; printed as
"Surface rate R", with "Live time" at the end, by both generators, PARMA
included) and by `gui/live_time.py` (the same quadrature, for planning before a
run). With mode 6 each tried muon is a draw from that flux, so n selected
rows correspond to a rate n / T. With legacy modes each row gets the exact
importance weight w = [J·proj / R] / [f_v(p) g(d|p) / A], reproducing the
Fortran samplers' densities (tabulated momentum CDF, 100-point Guan CDF), and
the rate is Σw / T. The GUI Results tab counts surviving rows only in
underground files and weights them with their surface kinematics. The
exposure-time panel uses the same R.

`test_live_time.py` (GUI path against UCMuGen `liveTime`, 120 000 hits per
filtered case):

| source | no filter | with filter |
|---|---|---|
| plane (60 m square, box detector) | +0.06% (±0.05%) | +0.50% (±0.29%) |
| disk (R = 30 m, cylinder) | +0.06% (±0.05%) | −0.24% (±0.29%) |
| hemisphere (R = 30 m, cylinder) | +0.03% (±0.07%) | −0.16% (±0.29%) |

R from Python and Fortran agree to 6 digits. Legacy modes with the filter,
weighted rate against UCMuGen: mode 2 +0.08%, mode 4 +0.30% (±0.29%);
unweighted they would be +3.6% and +0.8% in this geometry, and far more in inclined or
deep geometries (Section 7).

## 15. Detector filter and the tried count

**Intersections** (`src/generator/geom_module.f90`). Box: slab test on the box
inflated by the margin on every side; correct. Cylinder: the wall test uses
r + margin over the axial range [−margin, h + margin]; the caps sat at s = 0
and s = h (radius r + margin). For every ray from outside the inflated
cylinder that is equivalent to the full inflated volume: a ray entering
through the extended end face must leave through the wall inside the margin
zone or cross the original cap plane. `test_ray_intersections.f90` found 15
mismatches in 199 358 random rays, **all** with the origin inside the axial
margin zone, i.e. a source point inside the inflated detector. The caps now
sit at −margin and h + margin; 0 mismatches in 399k rays (tilted axes, rays
from inside, axial rays that only a cap catches). No run whose source lies
outside the inflated detector changes.

**Safety margin.** The margin is applied as stated: inflation of radius and
axial extent (cylinder) or of every face (box), per detector. It must exceed
the multiple-scattering displacement, because the filter only keeps muons
whose straight line crosses the inflated detector. `gui/mcs_margin.py`
integrates Highland (PDG Sec. 34.3) along the path with CSDA energy loss:

| 40 m of rock | σ_r (radial) | σ per projection |
|---|---|---|
| 50 GeV, ρ = 2.65 | 0.26 m | 0.18 m |
| 20 GeV, ρ = 2.65 | stops (range ≈ 35 m) | |
| 20 GeV, ρ = 2.0-2.2 | 0.68-0.79 m | 0.48-0.56 m |

PHITS gave 0.18 m at 50 GeV (agrees with the projected value) and 0.64 m at
20 GeV (between the radial and projected values for ρ ≈ 2.0-2.2; the PHITS
density and RMS definition would settle it). A point-like detector keeps a
fraction 1 − exp(−m²/σ_r²) of its hits: a 3 cm margin at 40 m keeps under 2%.
The GUI now warns below 2σ_r (98% kept) and the Helpers calculator suggests
2σ_r with energy loss included.

**Tried count.** Every attempt increments `ntry` once, by `!$OMP ATOMIC
CAPTURE`, before the only `cycle`, so no trial is lost. After the target is
reached each thread can add at most one more attempt, an overcount below
N_threads per rank (negligible against N_tried). Without the filter the OpenMP
generator resets `ntry = i`. The MPI generator sums the per-rank counts with
`MPI_Reduce`. Measured: 3 ranks × 2 threads, cylinder filter, 60 000 hits from
12.2 M tries: T = Tried / R = 205.8 s against UCMuGen's N_hits / rate_det =
206.8 s (−0.5%, stat ±0.4%).

## 16. Further bugs found and fixed

- **PARMA energy units (Fortran, both generators).** The energy CDF and the
  angular CDF evaluated PARMA at the total energy; PARMA takes kinetic energy
  (Sato 2008; `parma_subroutines.f90`). Effect on the sampled density ≤ 4%.
- **PARMA is not thread-safe.** Its angular routine memoises the last
  (energy, site) and derived coefficients in SAVE'd (Fortran) / static (C++)
  variables. Called from several threads it returns wrong values. In the
  Fortran this only matters for the new mode 6 (the legacy path calls PARMA in
  serial setup only), which reads a serially built table instead; measured
  +30% at 2-10 GeV, 60-80° with 8 threads before the table, correct with 1. In
  **UCMuGen** it affects every multithreaded Geant4 run with PARMA (one
  Generator per worker): 190 633 of 192 000 concurrent evaluations differed
  from single-thread values. The cache is now `thread_local`
  (`ucmugen/tools/make_parma_header.py`, which regenerates the header
  byte-for-byte from JAEA's source) and `parma::install()` loads the tables.
- **MPI generator angular modes.** Only 1-4 were accepted; mode 5 fell back to
  cos²θ silently.
- **GUI tried-count parser.** `re.search(r'Tried[:\s]+([0-9]+)')` matched the
  first progress line ("… / tried M"), not the final total.

Found during the v1.2.0 pre-release checks, when the examples were regenerated:

- **MURAVES example script** (`examples/vesuvius/ucmuon_vesuvius_muraves.py`).
  Its open-sky flux summed the spectrum over a log-energy grid without the
  factor E (∫f dE = ∫E f d ln E): 3-10× low, flat at ~10⁻³ m⁻²s⁻¹sr⁻¹ instead of
  1-58. Its through-rock flux called `backward_mc_flux` with `n_theta = 1`: a
  rate [m⁻²s⁻¹] over the cone 0..θ, evaluated at θ/2 and through
  X/(cosθ cos θ/2) instead of the slant opacity X. Both now call
  `directional_flux`, as `compute_flux_map` does; the open sky then agrees with
  the direct Guan integral to 4 % (the 50-point energy grid, which cancels in
  T_sim) and with PDG's cos²θ.
- **Terrain-engine summary rate** (`write_summary` in
  `gui/ucmuon_terrain_driver.py`, shown in the GUI's live panel). "Total
  expected rate" was `np.sum(flux_map)`, per-steradian fluxes summed over bins
  with no solid angle: ~100× too high on a 36×17 grid (2.3 × 10⁴ m⁻²s⁻¹) and
  dependent on the binning. It is now Σ Φ cosθ ΔΩ, the rate through a
  horizontal m²: 96.8 m⁻²s⁻¹ for an open sky with Guan above 1 GeV and
  θ < 85°, on both a 36×17 and a 360×85 grid, against 99.9 from the direct
  integral (the 3 % is again the energy grid).
- **MURAVES guide, Part 2.** The terrain-driver recipe had no spectrum line and
  four stray lines, so the driver read spectrum = 360, θ_max = 25°, a 1000 m
  ray step, and wrote to files named `1.0` and `2500.0`. Corrected, and checked
  by running it.

- **GUI Terrain tab, Section 4** (flux vs elevation) had the MURAVES script's
  bug: `backward_mc_flux` with `n_theta = 1`, divided by one pixel's solid
  angle. The open-sky curve was 190× to 1.4 × 10⁶× too high (more at larger θ
  and on finer grids). Now `directional_flux`.
- **Terrain-engine E_min** was 0.5 GeV, below Guan's and Frosin's fitted range;
  it is now 1 GeV, and the engine warns when a spectrum is used below its
  fitted range. CLI transmissions change ×1.00-1.17 (most near the vertical
  under thick rock), the open-sky flux ×0.82 at the vertical.
- **Range table.** The backward MC's range table ends at 2 TeV (3616 m w.e.).
  A detector energy whose surface energy is beyond it is left out, so under
  thick rock the flux is a lower bound: the cut removes about 0.6 % of the flux at 1000 m w.e., 6 % at 2000, 37 % at 3000 and 85 % at 3500 (θ = 0-60°, Guan; estimated by extending the table with dE/dX = a + bE fitted to its last points). It is exactly zero once every
  energy is beyond the table, where the true value is below 10⁻⁵ of the open
  sky. The table is not extended in this release; both integrators now raise
  a `RangeTableWarning` instead of cutting silently.
- `make_tsim_library.py` now reads the bundled `vesuvius_dem.tif` instead of
  `misc/dem_site.tif`, which is not in the tree: the two give a bit-identical
  overburden map.

## 17. What uses the changed values (to regenerate)

| where | uses | effect |
|---|---|---|
| `examples/vesuvius/tsim_library/*.dat`, `make_tsim_library.py` | terrain engine, backward-MC spectrum 1 | T ×1.48, nearly uniform. **Regenerated for v1.2.0** with Guan 1-2500 GeV (spectrum 1 extrapolated CosmoALEPH to 0.5 GeV); T in blocked directions rises ×300-19 000, 5th-95th percentile (×1.46 from the code, the rest from the spectrum), overburden identical |
| `examples/vesuvius/figs/fig_vesuvius_*.png`, `ucmuon_vesuvius_muraves.py` | backward MC spectrum 3 (Guan) | flux ×160-2600 at 60-85°, T ×35-1000; the script had its own flux bugs (Section 16). **Regenerated for v1.2.0** on the default cone model |
| `examples/vesuvius/MURAVES_GUIDE.md` | quotes "open sky flat ~10⁻³ m⁻²s⁻¹sr⁻¹" and "~10⁻⁸" floor, ρ̂ ≈ 1.98 | **Updated for v1.2.0** with the regenerated numbers (the old library gave no usable pixels in the guide's own inversion test) |
| T_sim libraries produced on HPC with the terrain engine (not in the tree) | backward-MC spectra | as the two rows above, by spectrum |
| Density tab: Methods 2-3 (`gui/ucmuon_density_analysis.py`, `gui/gui_density_analysis.py`) | fast estimator `reyna_bugaev` (default) | transmission curves change (Section 5): inferred opacities and densities shift |
| GUI exposure-time panel, Results-tab rates and live times | new R and weights | Section 14 |
| backward-MC tab outputs (rate, exposure) | backward MC | Section 5 |
| `benchmark/analysis/ucmuon_ffe_validate.py` | fast estimator | quick-reference numbers updated |
| UCMuGen PARMA Geant4 MT runs | PARMA | must be redone |
| UCMuGen runs with spectrum 7 | Reyna | rates ×0.08-0.26 |
| CPC manuscript | several | corrected in the development tree |

The transport benchmark (`benchmark/codes`, `benchmark/results`) and the
paper's survival tables use monoenergetic beams and are unaffected;
`benchmark/sources` evaluates Reyna only at cos θ = 1, which did not change.
