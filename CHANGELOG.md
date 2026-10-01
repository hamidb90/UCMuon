# Changelog

## [1.2.0] — 2026-10-01

Flux-normalisation audit: every surface spectrum in every code path checked
against its source paper, and against PDG and published data; then how the
generator samples momentum and angle, the rate and live-time chain, and the
detector filter. Full report with equation numbers, measurements and the
table of what changes: `docs/FLUX_NORMALISATION_AUDIT.md`.

### ⚠ Changes results
- **Reyna (spectrum 7) in UCMuGen**: horizontal-surface rates ×0.26 (E > 1 GeV)
  to ×0.08 (E > 100 GeV); the p-θ distribution now hardens with θ.
- **Fast flux estimator** (GUI default model, **Density tab**, exposure
  panel): Reyna replaced by Reyna's own formula (for example ×4.1 in flux at
  60° under 100 m of rock); plain Gaisser ×10; `gaisser_tang` is now the
  Tang et al. (2006) modified Gaisser formula (it was neither Tang's nor
  correctly normalised).
- **Backward MC and the terrain engine**: surface spectra replaced (backward-MC
  rates ×59 to ×10⁸; terrain transmissions ×1.48 with spectrum 1, ×35-1000
  with spectrum 3). The shipped Vesuvius T_sim library and MURAVES figures are
  regenerated; T_sim libraries made elsewhere must be.
- **Vesuvius T_sim library** (`examples/vesuvius/tsim_library/`): now Guan
  2015 over 1-2500 GeV; it was CosmoALEPH extrapolated down to 0.5 GeV, far below
  its ~100 GeV/c fitted range. Transmissions in blocked directions ×300-19 000; the
  MURAVES guide's own synthetic inversion goes from no usable pixels to
  ρ̂ = 2.007 for ρ_true = 2.0.
- **MURAVES example script**: open-sky flux 3-10× low (log-energy sum without
  the factor E), through-rock flux a cone-integrated rate at half the zenith
  angle; both now the backward MC's `directional_flux`. Figures and guide
  regenerated.
- **Terrain-engine summary** "Total expected rate" (summary file and the GUI
  live panel) summed per-steradian fluxes with no solid angle, ~100× high and
  grid-dependent; now the rate through a horizontal m² (Σ Φ cos θ ΔΩ).
- **GUI rates and live times**: one formula, T = N_tried / R (below). The old
  Results-tab rate was ×1.25-1.3 high for a plain horizontal source and up to
  ×8 with the detector filter.
- **PARMA in multithreaded UCMuGen**: concurrent evaluation was wrong (data
  race); multithreaded PARMA runs must be redone.
- **PARMA in the Fortran generator**: kinetic energy now passed to PARMA
  (≤ 4% change in the sampled spectrum).
- **MPI generator, angular mode 5**: previously replaced by cos²θ silently.

Legacy generator events (angular modes 1-5, all spectra but PARMA) are
bit-identical: 20/20 legacy configurations event-for-event equal.

### Added
- **UCMuGen 0.2.0** (`UCMUGEN_VERSION_*` in `UCMuGen.h`, was 0.1.0): a minor
  bump because its results change (Reyna, multithreaded PARMA) and its API
  grows (validity warnings); see `ucmugen/README.md`, "Versions".
- **Angular mode 6, joint sampling** (`generate_muon_joint`, both generators,
  PARMA included): draws (p, position, direction) from J(p, θ) × max(0, −n·d),
  the flux through the source surface, with the sky fixed in the world frame.
  Recommended, and the GUI default. Validated against `ucmugen::Generator` in
  six configurations (`tests/flux/test_joint_sampling.py`).
- **Live time.** The generators print "Surface rate R" (the flux through the
  source surface in the energy and zenith windows, the quantity UCMuGen's
  `rate()` computes) and "Live time = Tried / R". `gui/live_time.py` computes
  the same R and the exact importance weights of legacy-mode events; the GUI
  Results tab and exposure panel use it. Agrees with UCMuGen within 1% for
  plane, disk and hemisphere sources, with and without the detector filter
  (`tests/flux/test_live_time.py`).
- **Safety-margin guidance**: `gui/mcs_margin.py` (Highland with CSDA energy
  loss); the GUI warns when a detector margin is below 2σ_r of multiple
  scattering, and the Helpers calculator suggests one.
- **Validity guardrails** in all implementations, same limits (CosmoALEPH
  p ≳ 100 GeV/c; Gaisser E > 100/cos θ GeV; Guan, Frosin E > 1 GeV; Reyna
  1 < p < 2000/cos θ GeV/c), plus a notice for spectra with no absolute
  normalisation. UCMuGen: `flux::validity_warnings`,
  `flux::has_absolute_normalisation`, `Generator::warnings()`,
  `Generator::setPrintWarnings()`.
- Fortran `spectrum_intensity`, `surface_rate`, `geometry_h`,
  `source_normal_world`.
- Tests: `tests/flux/` (cross-implementation consistency, reference values,
  joint sampling, live time), `tests/geometry/` (detector intersections
  against brute force), a concurrency check in `ucmugen/validation/test_parma.cc`.
  All in CI.
- `docs/FLUX_NORMALISATION_AUDIT.md`.

### Fixed
- **GUI, a fresh install showed half the GUI empty.** With no output files yet,
  the Results tab called `st.stop()`, which ends the whole script: the Terrain,
  Config and Density tabs rendered empty and nothing was autosaved until the
  generator had run. The Results tab now only ends itself.
- **GUI without rasterio (optional, not installed by default).** The Terrain
  tab's DEM check called the terrain driver, which ran `sys.exit(1)` when
  rasterio was missing; a SystemExit is not caught by the GUI's error handling,
  so the script run ended there and every later tab, and the autosave, was
  lost. The driver now raises ImportError, which the GUI reports; its command
  line still exits with the same message.
- **GUI live time tied to its run.** The Results tab used the last generator
  run's tried count, surface rate and settings whatever file was loaded, and
  read the spectrum from the selector rather than the run (so a mono-energetic
  beam, run as spectrum 2 or 8, got a live time for the selected spectrum). A
  generator run now records its spectrum, start time and output files, a
  transport run the surface file each underground file came from; Results
  gives a live time only to a file from that run, and otherwise says why not.
  The run's live-time inputs are autosaved, so a restart keeps them
  (`gen_thetamax` was saved but never written; it is `gen_theta_max`).
- **GUI, Bethe-Bloch materials.** The menu labelled choices 2/3/4 Limestone /
  Water-Ice / Iron; both drivers run Ice / Water / Concrete. The menu now names
  what runs.
- **GUI, safety-margin warning** evaluated the depth at the detector's top
  face, so a detector reaching the surface (the default cylinder) never warned;
  it now uses the deepest face.
- **GUI, source too small.** New warning when a horizontal disk or rectangle
  source does not reach every straight path into the margin-inflated detector
  at θ ≤ θ_max. The default 200 m disk with the default 90 m-deep cylinder at
  85°: hit rate 11 % low (measured; 1 % at 600 m), and it needs R ≥ 1029 m.
- **GUI, Terrain altitude default** took the Generator's PARMA altitude in km
  as metres (three places).
- **GUI, Terrain cross-check legend** said "Guan" whatever spectrum was
  selected; the panel now also warns below a spectrum's fitted range.
- **GUI, Config "Reset autosave"** deleted the file, which the same run wrote
  straight back; it now also resets the session to the defaults.
- UCMuGen `flux::reyna` and Fortran `reyna_flux`: I(p, θ) = cos³θ · I_V(p cos θ)
  (Reyna 2006, Eqs. 1-3; was I_V(p cos θ*) with no prefactor). Bit-identical
  at cos θ = 1, so the momentum CDF and legacy streams are unchanged.
- `fast_flux_estimator`: Reyna formula; Gaisser prefactor 0.14 (was 1.4e-2);
  the `bugaev` key is now spectrum 6; `gaisser_tang` implements Tang et al.
  2006 (PRD 74, 053007) Eqs. 3-10, with Ẽ in the Eq. 7 exponent: the printed
  E cannot reproduce the paper's own Fig. 1 (10³-10⁵× low below a few GeV),
  Ẽ reproduces it at every angle shown.
- Backward MC surface spectra (none matched their paper).
- UCMuGen PARMA: memoisation cache and getPowCpp's scratch array `thread_local`
  (the array was caught by CI on Linux gcc; ThreadSanitizer is clean), tables loaded in
  `install()` (`tools/make_parma_header.py` regenerates the header).
- Fortran PARMA: kinetic energy; mode 6 reads PARMA's angular factor from a
  serially built table (PARMA's routines are not thread-safe).
- MPI generator accepts angular modes 5 and 6.
- Detector filter: cylinder caps at −margin and h + margin, consistent with
  the wall test. Identical for every source outside the inflated detector.
- GUI: exposure panel passed total energy where kinetic was expected; the
  tried count was parsed from the first progress line, not the total; spectrum
  and angular-mode help texts corrected (mode 6 was labelled "Bugaev 1998",
  Reyna "1-10000 GeV" and "~20%" PDG agreement, Guan "suppresses 50-100×",
  mode 4 "physically exact").
- MURAVES guide, Part 2: the terrain-driver recipe was missing the spectrum
  line and had four stray lines (the driver read spectrum = 360 and wrote
  files named `1.0` and `2500.0`); corrected and run. Part 2 now points to the
  shipped library.
- `examples/vesuvius/make_tsim_library.py` reads the bundled
  `vesuvius_dem.tif` (it needed `misc/dem_site.tif`, not in the tree); the
  overburden map is bit-identical, so the library now regenerates from the
  published tree alone.
- **GUI Terrain tab, Section 4** ("Integrated muon flux vs elevation") had the
  MURAVES script's bug: a cone-integrated rate at half the zenith angle,
  divided by one pixel's solid angle, 10²-10⁶× too high and grid-dependent.
  Now `directional_flux`.
- **Backward MC beyond the range table**: `directional_flux` and
  `backward_mc_flux` raise a `RangeTableWarning` when surface energies above
  2 TeV are left out, instead of returning a zero silently; the CLI summary,
  the results file, the backward-MC tab, the terrain map and the MURAVES
  script report it (`tests/flux/test_backward_mc_range.py`).
- **Terrain engine E_min**: `compute_flux_map` (and so the CLI) now integrates
  from 1 GeV, the lower end of Guan's and Frosin's fits, not 0.5 GeV: CLI
  transmissions rise by up to 17 % near the vertical (under 2 % at 80°). It
  warns when a spectrum is used below its fitted range, and for the power law,
  which has no absolute normalisation.
- `hpc/input_params.dat` is a neutral, valid demo (Guan, angular mode 6, one
  borehole-like detector, margin 2σ_r of multiple scattering, source radius
  from the detector's view cone), with each choice explained;
  `hpc/run_ucmuon_gen.sh` prints the surface rate and live time in its summary.
- `tests/flux`: pass thresholds set from the statistics (χ² p-value > 10⁻³,
  4σ): the first thresholds failed a correct build about one run in eight.
- THIRD_PARTY_LICENSES.md and `make_parma_header.py` list all four mechanical
  edits in `UCMuGen_PARMA.h`, including `thread_local`.
- CITATION.cff: titles of MUSIC, Guan, Frosin and PARMA entries, Frosin DOI,
  Highland year and title; added Reyna, Schmelling (CosmoALEPH), Bugaev,
  Sato 2016, Tang 2006 and the PDG 2022 cosmic-ray review. PDG values are
  cited from PDG 2022 Sec. 30.3.1 and Eq. 30.4 (the 2024 edition dropped the
  section). README references and BibTeX (placeholder DOI) corrected.

### Checked, unchanged
- Guan and Frosin constants, cos θ\*, the dN/dE → dN/dp Jacobian and the
  PARMA unit conversion in UCMuGen.
- CosmoALEPH's constants are exactly a log-log fit to Schmelling et al. 2013
  Table 1; its isotropic angular law in UCMuGen is kept.
- "Total tried" under MPI + OpenMP: no lost trials (measured against UCMuGen,
  −0.5% ± 0.4%).

### Known limitations
- The backward MC's range table ends at 2 TeV (3616 m w.e.), so under thick
  rock its flux is a lower bound, and zero once every energy is beyond it (the
  true value is then below 10⁻⁵ of the open sky). Measured: it removes about 0.6 % of the flux at 1000 m w.e., 6 % at 2000, 37 % at 3000 and 85 % at 3500 (θ = 0-60°, Guan; estimated by extending the table with dE/dX = a + bE fitted to its last points). This now raises a
  `RangeTableWarning`, and flux maps print one summary line; extending the
  table is left for a later release.
- The terrain-driver CLI has no energy-range input (1-5000 GeV).

## [1.1.2] — 2026-09-20

Publishes the benchmark suite and the paper's figure scripts. No code changes:
every engine, input file and GUI element is byte-for-byte what 1.1.1 shipped.

### Added
- **`benchmark/`**, the six-code validation campaign that the paper reports:
  the four UCMuon engines against Geant4 11.2 and PHITS 3.36 on one identical
  source population (six monoenergetic vertical beams, 5 to 300 GeV, 10^5 muons
  each, Standard Rock, scored at 1, 10, 25, 50, 100 and 200 m).
  - `codes/` one subdirectory per transport code, each with the setup and the
    instructions to run it; `analysis/` the cross-code scripts;
    `results/` the summary CSVs and the per-engine timing files;
    `figures/v2_six_code/` the canonical plots;
    `reports/BENCHMARK_FEEDBACK.md` the reference write-up of the v2 run.
  - `benchmark/README.md` is the roadmap for reproducing it end to end.
  - The roughly 8 GB of raw per-event output stays untracked and is regenerated
    by the per-code steps; what ships is the distilled product.
- **`manuscript/scripts/`**, one script per figure of the paper, so every figure
  can be regenerated from the code and the benchmark products. The paper itself
  stays out of the repository.

### Fixed
- `benchmark/README.md` numbered the engines wrongly, giving MUSIC and UCMuon-MC
  both as "Engine 1" and swapping Bethe-Bloch and PROPOSAL. The numbering now
  matches the suite: 1 UCMuon-MC, 2 MUSIC, 3 Bethe-Bloch, 4 PROPOSAL.
- `benchmark/README.md` listed the PHITS exit-energy disagreement as an open
  item to be fixed before submission. It is a reported result, not a defect,
  and it is in Section 8 of the paper and in `reports/BENCHMARK_FEEDBACK.md`.
- `reports/COMPARISON_G4_PHITS.md` linked to a superseded four-code write-up
  that is not published, and gave two commands as absolute paths on the author's
  machine.

### Note
Not everything in the development tree is published: the speaker notes and the
internal planning document are not, and neither is the superseded four-code
comparison and its figures, because publishing numbers we know to be superseded
next to the current ones invites exactly the wrong comparison.

## [1.1.1] — 2026-09-20

The version described in the Computer Physics Communications submission.

### Added
- **`test_run/`**, the comprehensive test run required by the CPC program
  submission. `bash test_run/run_test.sh` runs the full two-stage pipeline
  (Fortran surface generator to UCMuon-MC transport) on a seeded 2000-muon
  configuration in a few seconds and grades the result against committed
  reference output: 2000 CosmoALEPH muons, 10-2500 GeV, disk source R = 50 m,
  cos^2(theta) to a 70 degree zenith cut, through 25 m of Standard Rock, of
  which 574 survive (28.70%) with a mean exit kinetic energy of 13.562 GeV.
  - Graded numerically rather than byte for byte (exact muon count, survival
    fraction within a quarter of its binomial sigma, mean exit energy within
    0.5%), because a different `gfortran` or C library can change the last bits
    of a double without anything being wrong. Byte identity is still reported
    when it holds, as it does on the reference platform.
  - Runs without a Fortran compiler: stage 2 is pure Python, so when
    `bin/ucmuon_gen_omp` is absent the script falls back to the committed
    surface file and still grades the transport stage.
  - `test_run/expected/SUMMARY.txt` records the reference numbers, SHA-256
    sums, and the compiler, Python and NumPy versions they were produced with.

### Fixed
- `test_run/run_test.sh` reported `byte-identical: muons_surface.dat` even when
  stage 1 had been skipped for want of a Fortran compiler, in which case that
  file was the reference copied in and was only ever compared with itself. It
  now says the file was not checked, and grades stage 2 as before.
- The v1.1.0 changelog entry quoted the detector-directed sampling speed-up as
  5800x; the measured figure is 5700x, as reported in `ucmugen/README.md` and
  the paper.

### Changed
- README: the status heading still said v0.9.0 and "the first public release";
  the directory tree listed neither `ucmugen/` nor `test_run/`, and, in the
  public tree, listed a `references/` directory that is not published. A
  "Verify the installation" section now points at `test_run/run_test.sh`.

## [1.1.0] — 2026-08-05

### Added
- **UCMuGen** (`ucmugen/`), a single-file, dependency-free C++17 port of the
  surface generator for direct use inside Geant4. Integration is one `#include`
  and one call. It implements the seven muon flux parametrisations, the five
  angular models and the momentum-dependent charge ratio of the Fortran
  generator, and adds:
  - correct through-surface projection on flat, disk, hemispherical and
    cylindrical surfaces, so any spectrum works with any geometry;
  - an absolute rate, and hence a live time, for any spectrum and surface;
  - detector-directed sampling for a box, sphere or oriented cylinder,
    measured at 5700x faster for a 2 m detector under 200 m of rock;
  - optional location- and date-aware PARMA/EXPACS flux via a separate header
    (`UCMuGen_PARMA.h`, generated by `ucmugen/tools/make_parma_header.py`).
  Validated against the Fortran generator as *exact event equality* rather than
  a statistical tolerance: 20 of 23 configurations bit-identical at 20k events,
  3 differing by design. See `ucmugen/README.md`.
- Continuous integration for UCMuGen across four compiler and platform
  combinations, plus a bit-exactness job against the Fortran generator
  (`.github/workflows/ucmugen.yml`).
- `ucmugen/examples/geant4/`, a complete Geant4 application: geometry, physics
  list, macros and CMake, so the integration can be run and not only read. It
  runs multithreaded, and reports its predicted rate against the rate Geant4
  measures (325.517 Hz against 325.234 Hz, 0.09%), which checks the flux
  normalisation, the surface projection and the unit conversions end to end.
- `ucmugen/examples/features/`, a tour of the whole API in one file: all eight
  spectra, all four surfaces, all three detector shapes, the angular ranges,
  seeding, normalisation and the legacy driver. No Geant4 and no external data,
  about two seconds to run.
- `ucmugen/comparison/`, the two measured UCMuGen vs EcoMug programs, both
  driven by an identical differential flux so what is compared is the machinery
  rather than the default parametrisations.

### Fixed
- **MPI generator kept only rank 0's muons.** `output_all`, `output_sel` and
  `output_phits` are `character(512)` but were broadcast as 120 characters, so
  on every rank except 0 the tail of the string was uninitialised memory.
  `trim()` kept it, the result overflowed the per-rank filename buffer, and the
  `_RRRRR.dat` suffix was truncated away: ranks 1..n-1 all opened the same
  malformed file and overwrote each other. An 8-rank job wrote an eighth of the
  muons it reported, because the totals come from a reduction over counters
  rather than from the files. Affects every MPI run of `bin/ucmuon_gen`;
  single-process and OpenMP runs are unaffected.
- RNG seed collision: runs started within the same minute produced identical
  event streams.
- UCMuGen: `parma::install()` writes a process-wide provider, so calling it
  from a per-thread constructor was a data race. The Geant4 example installs it
  once before the workers start, and the constraint is now documented on the
  provider, on `install()` and in both example READMEs.
- UCMuGen: `UCMuGen.h` includes the Geant4 headers its `FireG4` block needs
  instead of assuming the including file arranged its own includes in a
  particular order. Defining `UCMUGEN_WITH_GEANT4` before the include used to
  fail a long way from the cause.
- UCMuGen: PARMA now supplies its own charge ratio. It models the two charges
  separately, so its ratio varies with atmospheric depth and cutoff, where the
  built-in table is a sea-level fit that is flat below 112 GeV/c. At 5 km the
  fit is high by about 12% at 1 GeV/c.
- The MUSIC data tables are now covered by `.gitignore` alongside the MUSIC
  sources. `docs/MUSIC_FILES.md` asks the user to paste them into `data/` to
  enable Engine 2, and the existing rule reached only the repository root.
- `hpc/run_ucmuon_gen.sh` reported the wrong build target when `bin/ucmuon_gen`
  was missing: that binary is MPI and comes from `make hpc`, not `make local`.
- The PUMAS "not built" note pointed at `docs/MUSIC_FILES.md` instead of
  `setup.sh`, which offers to download it.
- PARMA/EXPACS licence terms were recorded incorrectly. The published
  conditions **do** grant redistribution and modification for non-commercial
  use; the repository previously stated otherwise.
- Engine timings re-measured. The published figures predated the July engine
  work and overstated the two Python engines by large factors: UCMuon-MC is
  ~3.7 s per 1e5 muons at 100 m rather than ~14 s, Bethe-Bloch ~1.5 s rather
  than ~40 s. MUSIC, untouched since June, re-measures unchanged and serves as
  the control.

### Changed
- `data/EXPACS/` now ships only the PARMA data the engine reads
  (`parma/input/`), not the whole EXPACS distribution. Every file the engine
  opens lives under `input/`; the sample run output, the dose-calculator inputs
  and the EXPACS spreadsheet were unreachable. The model is unchanged, verified
  bit-identical against the full data. 21 MB to 3.4 MB.
- Terminology: "pre-filter" named two unrelated mechanisms and is replaced by
  **detector acceptance cut** (geometry) and **range cut** (CSDA). In
  muography "target" denotes the imaged body, so the detector volume is never
  called a target.

## [1.0.2] — 2026-07-21

### Added
- Bundled sample DEM `examples/vesuvius/vesuvius_dem.tif` (Mt. Vesuvius,
  SRTM GL1 30 m, public domain — provenance and citation in
  `examples/vesuvius/DEM_SOURCE.md`). The GUI Terrain tab now loads it by
  default when no DEM has been uploaded or downloaded, so the terrain
  workflow is runnable out of the box on a fresh install.

### Fixed (GUI)
- Terrain tab, Run & Results: switching the results view (e.g. selecting
  "3D Terrain") or touching any widget inside a view no longer snaps the
  view strip back to the first tab. The inner `st.tabs` (selection kept
  only client-side) was replaced with a keyed `st.segmented_control` whose
  selection persists across reruns.
- muRAvES / Vesuvius preset button now sets the documented MURAVES detector
  position (40.8271 °N, 14.4006 °E, 608 m) and loads the bundled DEM.

### Changed (GUI / examples)
- Terrain "MURAVES Comparison" result view renamed to "Literature cross-check"
  and redesigned: the 2D rock-thickness map now comes before the azimuth
  slice so it can be used to choose the target direction, the intro no longer
  advertises plots that do not exist, and the dead `ucmuon_mulder_crosscheck.py`
  reference was removed.
- Removed all references to an internal collaboration presentation (named
  colleague, meeting date, and slide numbers) and the reference data digitised
  from it, across the Terrain GUI, the Vesuvius example, and its guide.
- Corrected the Mt. Vesuvius reference: "Tioukov et al. 2019, Sci. Rep. 9, 6695"
  is the *Stromboli* muography paper, was not the source of the plotted data,
  and is replaced everywhere by the published MURAVES Vesuvius paper Hong et al.
  (2025), J. Appl. Phys. 138, doi:10.1063/5.0275078. The hardcoded thickness
  "reference" curve (actually unpublished preliminary simulation) was removed;
  the cross-check now shows only UCMuon's own curve and points to the paper.
- Citation audit: removed the "Lo Bue et al. 2023, JGR 128, e2022JB025446"
  reference (DOI does not resolve; the real R. Lo Bue paper is Etna seismic
  tomography, unrelated to Vesuvius muography). Corrected the Highland multiple-
  scattering reference year everywhere from 1979 to 1975 (NIM 129, 497 (1975)).
  The Frosin spectrum reference (J. Phys. G 52, 035002, 2025) and its fit
  parameters a=3.512, b=1.388 were verified against the paper (Table 4) and
  are correct. Corrected the PROPOSAL-update reference volume from CPC 305 to
  CPC 302 (Alameddine et al. 2024, CPC 302, 109243).
- All source-spectrum parameters were verified against the primary PDFs and
  match exactly: Guan P1-P5 + a,b (Guan Table 1), Gaisser constants, Frosin
  a,b (Table 4), Reyna c1-c5 (Reyna Eq. 3 best fit), CosmoALEPH charge-ratio
  table and power-law (Schmelling Table 1), and the Bugaev four-range Table II
  coefficients (the code correctly uses the 4.1625e5 breakpoint; the Reyna
  paper misprints it as 41625).

### Changed (installers)
- setup.sh: corrected the "Engines 2–6 are fully functional" note (stale
  range) to "All other engines (1, 3–7) are fully functional".

## [1.0.1] — 2026-07-19

Patch release: Windows support, PUMAS forward-mode fix, installer overhaul.
Verified end-to-end on Windows 11 (6 of 7 engines) and a fresh macOS
machine (all 7 engines).

### Fixed (critical)
- PUMAS engine, forward mode: the RNG was never seeded, so straggled runs
  hung above ~285 GeV (NaN energies) and mixed runs sampled hard losses
  with degenerate randomness. Forward transport now seeds the generator
  (new optional seed input; 0 = time-based). CSDA and backward-mode
  results were never affected.

### Windows support
- install.ps1 no longer crashes under PowerShell 5.1; it can now install
  MSYS2 + gfortran automatically (winget/pacman) and builds through the
  MSYS2 UCRT64 shell. New install_windows.bat double-click wrapper.
- RANMAR/RANLUX thread-private state moved from COMMON blocks to modules
  (bit-identical sequences), fixing the MinGW assembler failure that
  blocked the MUSIC and Bethe-Bloch builds on Windows.
- GUI resolves .exe binaries, adds the MSYS2 DLL directory to PATH, and
  disables Run buttons with a clear message when a binary is missing.
- run_gui.bat: thread count no longer depends on the removed wmic tool.

### Installers
- Both installers offer to download PUMAS (LGPL-3.0, github.com/niess/pumas)
  and build Engine 7 automatically instead of requiring manual steps.
- setup.sh: fixed a stale source-file check (cosmicray.f90) that caused
  the generator build to be skipped on every machine.
- Optional rasterio install is constrained to numpy<2.3 so pip cannot
  break pre-existing scipy installations.

### Changed (GUI)
- PUMAS defaults to forward transport; its underground detector filter is
  hidden in backward mode (flux output has nothing to filter).
- Bethe-Bloch Run button relabeled "Bethe-Bloch CSDA"; default minimum
  generator energy is 100 GeV to match the CosmoALEPH validity range.

## [1.0.0] — 2026-07-12

First stable release. Full pre-release verification passed (31-point check:
generator statistics, five-engine cross-validation against Geant4/PHITS,
GUI regression suite).

### Physics fixes
- Multiple scattering: polar deflection now drawn from a Rayleigh(θ₀)
  distribution (previous Gaussian draw gave √2-low RMS deflection) —
  Bethe–Bloch and UCMuon-MC engines; validated against the Highland
  expectation and Geant4.
- Bethe–Bloch engines: radiative-loss coefficient b(E) rebuilt on the
  PDG-2024 shape (`b_rad_shape`); Python and Fortran BB now agree to
  <0.4 GeV in mean exit energy at 200 m.
- Generator: removed the last-bin spike in the power-law spectrum (mode 2);
  uniform-cone angular mode now samples uniformly in solid angle;
  Reyna–Bugaev (mode 7) spectrum corrected.
- Energy conventions unified: 18-column output E is total energy for
  survivors and 0 for stopped muons in every engine; stopped muons are
  written to a companion file with initial kinetic energy and stopping depth.
- PROPOSAL driver: custom materials are now true custom media built from
  (Z, A, I, ρ) via Sternheimer density-effect parameters (previously
  transported as density-scaled Standard Rock).

### Added
- Vesuvius / MURAVES worked example with a shipped 5-density transmission
  (T_sim) library and DEM download recipe.
- MUSIC energy-loss table self-generation (driver and data included;
  the MUSIC source itself is obtained from its author — see
  docs/MUSIC_FILES.md).
- App version string in the GUI (v1.0.0, synced with CITATION.cff).

### Fixed (GUI)
- Data files containing only a header now produce a clear warning instead
  of a generic load error.
- Silenced spurious numpy warnings from masked divisions (axis-parallel
  cylinder intersection, terrain transmission map).

### Known limitations (see docs)
- MPI binaries require a cluster rebuild (all fixes are in the source).
- PROPOSAL writes ~0.05% of survivors with KE = 0 (boundary crawlers,
  cosmetic).
- The Terrain tab requires `pip install rasterio` (graceful error otherwise).

## [0.9.0] — 2026-06-24

Initial public release (concept DOI 10.5281/zenodo.20826984).
