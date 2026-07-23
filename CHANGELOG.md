# Changelog

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
