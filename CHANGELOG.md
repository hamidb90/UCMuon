# Changelog

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
