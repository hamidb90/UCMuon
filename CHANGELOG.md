# Changelog

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
