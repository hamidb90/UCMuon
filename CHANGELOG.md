# Changelog

## [1.4.0] — 2026-10-16

Fixes from a line-by-line review of the paper against the code (2026-10-08)
and from a six-part audit of the whole code (2026-10-14 to 16: generator and
flux, transport, terrain and density, GUI, UCMuGen, converters and tooling),
each finding reproduced by a run and checked by an independent second pass.
The data shipped with the code that the fixes change are regenerated.

### ⚠ Changes results
- **Multiple scattering** (UCMuon-MC, Bethe-Bloch Python and Fortran): the
  Highland variance accumulates along the path, so angles and lateral
  displacements no longer depend on the step size. UCMuon-MC's mean angle
  against Geant4 goes from −15.4 % to −2.8 % at 1 m of rock and is within 5 %
  to 200 m. Survival and energies are unchanged.
- **Random deviates** of the Fortran generators and UCMuGen carry 53 bits
  instead of 24: a fixed `UCMUON_SEED` gives different (statistically
  equivalent) files than 1.3.x, and the far momentum tail is no longer cut.
- **Terrain:** the forward transport put 1/cos²θ too much rock in front of
  every muon; the ray tracer interpolates the DEM bilinearly and bisects rock
  boundaries; projected DEMs are reprojected. The Vesuvius library gains 93
  rock directions.
- **Backward MC (Engine 5):** the range table goes to 10 TeV (was 2 TeV, which
  left out part of the flux under more than about 1000 m w.e.); water, ice,
  seawater and iron use their own PDG tables instead of a rescaled rock table.
- **Density inversion** interpolates ln T between library densities (bias
  between nodes +2.6 % → +0.5 %).
- **PARMA dates after the end of its table** got W = 0 (+6 % rate) in the
  generator and UCMuGen; they take the latest measured W.
- **UCMuGen 0.3.0:** the directed-sampling cone is exact (directed
  acceptances about 16 % higher), and the Plane's axes are the documented ones.
- Bethe-Bloch in water and ice, pencil beams, cosmic electrons through the
  command-line converters, MPI EventIDs: see Fixed.

### Fixed
- **Bethe-Bloch, Python driver, ice:** the radiation length was 33.10 g/cm²;
  it is 36.08 g/cm², as water, the Fortran engine and UCMuon-MC use (PDG).
  This changes mainly the multiple-scattering angles and lateral displacements
  of Bethe-Bloch runs in ice with multiple scattering on (θ₀ 4 to 5 % smaller).
- An autosave file holding valid JSON that is not an object (a list, say)
  crashed the GUI on start; it is now ignored.
- Config tab: the survival fraction is a caption under "Survived", no longer a
  green up arrow (the same fix as for the rate in 1.3.2).

Found by a code audit (2026-10-14):
- **PARMA, W index from the date:** any date after the last year of the
  neutron-monitor table (2026) got W = 0, solar minimum, without a warning:
  `getHP` took the 0.0 its table is initialised with as data. Such dates now
  use the latest day with data and say so, as missing days inside the table
  already did. The GUI's Advanced mode defaults to today's date, so from
  1 January 2027 every PARMA run would have been affected (+6 % surface rate
  at a site at 50.7° N, 0.1 km). The search back now reaches the table from
  any date the GUI accepts (up to 2050). Dates inside the table are unchanged.
- **MPI generator (`bin/ucmuon_gen`):** an E_min below 1 GeV was silently
  raised to 1 GeV, and spectrum 8 (cosmic electrons, listed in
  `hpc/input_params.dat`) silently became CosmoALEPH muons. It now keeps the
  energy window and generates electrons, as the OpenMP generator does.
  Without a detector, each rank's file now holds exactly its share of rows
  and "Total tried" equals "Total saved" (the extra muons the threads drew
  past the target were written before; the live time was consistent).
- **GUI, loading a PHITS dump:** a dump written with Fortran D exponents (the
  generator's own, the GUI's PHITS export, `ucmuon_to_phits.f90`) could not
  be loaded ("unsupported operand type(s)"); it now loads. Electrons (kf ±11)
  now get their own charge and mass, not the muon ones.
- **Bethe-Bloch in water and ice:** the Python driver rescaled the Standard
  Rock table by 1.046 (Z/A 0.523; water's is 0.555, and its radiative loss is
  25 % below rock's), and the Fortran engine (OpenMP and MPI) used Z/A 0.495
  and b = 3.40e-6 (PDG 2.28e-6). Muons lost too little energy: 100 GeV
  muons through 100 m of water, 5 % (Python) and 8 % (Fortran); 10 GeV
  through 10 m, 10 % (Fortran). The Python driver now uses the PDG 2024 ice
  and water tables (as UCMuon-MC), the Fortran engine H₂O's Z/A and the PDG
  b; both now agree with the PDG table (the Fortran one within its 2 %
  analytic-ionisation residual, as in rock). Standard Rock is unchanged.
- **Python Bethe-Bloch, MS off and range cut:** the analytical (MS off) path
  and the range cut multiplied the rock range by the material factor instead
  of dividing, so for a custom material the energy loss changed when MS was
  switched (10 GeV, Z/A 0.55: 17 %), and a material with Z/A below rock's
  (iron) had every muon cut. Concrete and custom materials now scale
  ionisation and radiative losses separately (as the Fortran engine and
  UCMuon-MC). Standard Rock is unchanged.
- **Exit point of survivors:** UCMuon-MC, the Python Bethe-Bloch and PROPOSAL
  wrote it at z = -depth, while they transport a slab of that depth below
  the source; for a source plane not at z = 0 (Advanced mode: plane
  coordinate W, tilted disks, hemisphere) the GUI's underground detector
  filter then shifted inclined muons by W tan θ. PROPOSAL also wrote
  z = -depth for vertical (XZ, YZ) source planes. All engines now write it
  depth below each muon's start, as the Fortran engines, MUSIC and PUMAS
  did. Sources at z = 0 are unchanged.
- **Pencil beams:** with every muon starting at one point, MUSIC and the
  Fortran Bethe-Bloch (OpenMP and MPI) took the depth along Y and PROPOSAL
  and PUMAS along X, so no muon survived. They now take it from the
  momentum, as UCMuon-MC already did. Generator files are not affected.
- **UCMuon-MC radiation length in the Terrain tab and the CSG detector step:**
  both passed the 2.65 g/cm³ value in cm with another density, so multiple
  scattering scaled as sqrt(2.65/ρ): +16 % in angle at ρ = 2.0, about 80x for
  the default air detector cell. Survival and transmission were not
  affected.
- **Multiple scattering (UCMuon-MC, Bethe-Bloch Python and Fortran):**
  Highland's θ₀ was computed for each step and the steps added in
  quadrature. Highland's formula describes a whole layer, and that sum depends
  on the step size (about 10 % in angle per tenfold change) and is too small:
  one radiation length of rock gave 22 to 30 % (UCMuon-MC) and 11 %
  (Fortran) less than Highland for that layer, and UCMuon-MC's mean
  deflection in the monoenergetic Geant4 benchmark was 15 % low at 1 m
  (its step there is 0.9 g/cm²). The variance is now accumulated along the
  path: Highland applied to the layer crossed so far up to one radiation
  length, growing linearly with the path beyond (the log term, a fit
  validated to 100 X₀, would put the angles 19 to 23 % above Geant4 in 25
  to 200 m of rock). The result no longer depends on the step (0.2 % over
  a hundredfold range), reproduces Highland for a single layer up to one X₀,
  and puts UCMuon-MC within 4 % of Geant4 from 1 to 200 m. UCMuon-MC also took βp from
  p = E instead of p = √(E² − m²) (angles 0.6 % low at a total energy of
  1 GeV, 7 % at 0.3 GeV). Survival and energies are unchanged; angles and
  lateral displacements change.
- **Random deviates (generators, UCMuGen.h):** only the top 24 bits of the
  64-bit generator were used, as a single-precision value, so every deviate
  was a multiple of 2^-24 and the far momentum tail was cut off and
  quantised: Guan 1-10000 GeV/c never exceeded 5923 GeV/c, and 36 muons
  above 1357 GeV/c in 4·10⁶ took only 28 distinct values. The deviates now
  carry 53 bits. The stream advances as before and each old deviate is the
  leading 24 bits of the new one, so seeded runs change only where a
  decision falls within 2^-24; UCMuGen.h still agrees with the Fortran bit
  for bit, and the test run still gives 574 of 2000 survivors (28.70 %,
  13.562 GeV), but its files are no longer byte-identical to the reference.
- **Terrain ray tracer (`gui/ucmuon_terrain_driver.py`):** it counted whole
  steps where the sample was below the DEM, with no boundary refinement
  (its docstring said it bisected). Each segment starting at the detector
  lost about half a step: a detector 100 m under a flat surface
  (underground deployment, default 50 m step) got 50 m of rock vertically
  and 100 m instead of 141 m at 45°, and short paths on a flank came out as
  0 or a multiple of the step. The ray is now sampled at a quarter of the
  DEM pixel (finer next to the detector) and every rock boundary is
  bisected; the flat slab is exact at any step. The DEM lookup was meant to
  be bilinear, but rasterio's `rowcol` returns integers, so every weight was
  zero and the surface was a staircase of pixels; with exact boundaries the
  steps put metres of spurious rock in front of a flank detector in every
  direction. It is now bilinear between pixel centres. On the shipped
  Vesuvius library this blocks 95 more bins (median 38 m of rock) and
  lengthens the paths by 3 % on average above
  2·10⁵ g/cm² (9 % at 0.5-2·10⁵); the library is to be regenerated.
- **DEM files in a projected CRS:** a GeoTIFF in UTM or Lambert had its
  metre transform read as degrees, and XYZ and ASC files were converted
  with a formula that ignored the UTM zone and false easting (Vesuvius in
  UTM 33N landed at 5.3° E): the detector fell outside the DEM and every
  direction was open sky. A projected GeoTIFF (or an .asc with its .prj) is
  now reprojected to longitude/latitude; XYZ and ASC files without a CRS
  must be in longitude/latitude and otherwise stop with an explanation.
  ASC files in `xllcenter` convention were shifted half a cell north-east,
  and one without the optional `NODATA_value` line lost its first row.
- **Backward MC, water/ice, seawater and iron:** the rock range was
  multiplied by the material factor a where it should be divided, so a
  material acted like X/a of rock instead of a·X, and one factor stood for
  all losses. Rates under water came out 11, 29 and 42 % too high at 10, 50
  and 200 m, under iron 13 to 27 % too low at 2 to 30 m. These materials now
  use their PDG 2024 tables (as UCMuon-MC), and the b of the survival
  factor follows them. Standard Rock, and so the terrain engine and its
  libraries, is unchanged.
- **Density analysis, loading a transmission map:** grid keys were the
  parsed angles but lookups rounded them to 8 decimals, so a map written at
  full precision (`np.savetxt`'s default, with the default 18 zenith bins of
  85/18°) loaded with 576 of 648 cells NaN. Through the GUI the library
  check usually rejected such a file; loaded as measured data first, or
  through the API, it went through silently.
- **STL geometries:** any error in the ray query came back as open sky; with
  trimesh installed as the loader asked (`pip install trimesh`) but without
  `rtree`, which its ray queries need, every direction was open sky. A
  detector inside the mesh (one exit per ray) also got open sky, and a ray
  through an edge shared by two triangles could count zero rock. Loading an
  STL now checks the ray engine and says to install `rtree`; inside
  detectors and shared-edge hits are handled. (An STL feeds the overburden
  preview and the single-ray check; it has no cells, so they trace from the
  geometry origin.)
- **Source-size panel, θ_max advice:** for the paper's underground example
  (r = 10 cm cylinder 90 m down, E_min 45 GeV) it said "Reduce to 0.7°", the
  angle beyond which scattering exceeds the detector size; cutting θ_max
  there would bias every rate low by orders of magnitude. It now only ever
  raises θ_max (for coverage) and explains the scattering limit instead. Its
  MCS margin was 1.5 σ_r at a constant surface momentum (102 cm) where the
  Generator's warning asks for 2 σ_r with energy loss (231 cm); it now uses
  the same rule. The CSG preview and single-ray check ignored the step and
  range set in the CSG panel (they read keys nothing writes); they use them.
- **UCMuGen (0.3.0), PARMA W index of a date:** `parma::w_index_from_date`
  returned W = 0 (solar minimum) for every date after the last measured day of
  the neutron-monitor table (2026-03-26, so also today) and for every later
  year, where the Fortran generator uses the latest measured W (107.4): the
  PARMA rate came out 5.9 % high (sea level, 50.7°N 4.4°E, 1-1000 GeV). It
  now falls back the same way; dates before 1951 use JAEA's yearly table
  (they returned 0), and W is clamped at -135.4 as in the Fortran (one
  ground-level-event day held -752, which made the spectrum unusable). A
  sixth mechanical edit in `tools/make_parma_header.py`; the header is
  regenerated from JAEA's source. The `Site` documentation called W = 0
  "mean activity"; it is a deep solar minimum (yearly means 6.6 in 2009,
  long-term mean 58.5).
- **UCMuGen, detector-directed sampling:** the largest cone the detector
  subtends over the generation surface, the bound of the acceptance test,
  was found by sampling 20 000 surface points. For a small detector close
  under a large surface it missed the true value (a 1 cm sphere 2 cm under a
  10 m disk: 0.05 sr instead of 0.84), so the events near the closest point
  were under-weighted: 6.5 % within 2 cm of the axis instead of 29 %. It is
  now exact for the built-in surfaces (`Surface::distanceTo`); the rate was
  not affected. `setSurface()` and `setSeed()` now rebuild the envelope
  (after a first run a new surface kept the old one: 14 000 envelope
  violations in 20 000 events).
- **UCMuGen, `Plane(half_x, half_y, ...)`:** half_x ran along y for the default
  normal, so `Plane(200, 20)` over a bar along x gave a quarter of the rate
  into it (8.4 instead of 31.0 /s). Every Plane in the repository is square.
- **UCMuGen, energy range:** E_min at or below the particle mass gave a NaN
  momentum table and a `generate()` that never returned; E_min > E_max became
  a beam at E_min. Both now throw. The acceptance envelope also scans
  log-spaced zenith cosines near the horizon, where Reyna's weight peaks above
  10 TeV (outside its validity), and `liveTime()` takes the number of
  integration points.
- **UCMuon-MC, where a range-cut muon stops (water, ice, seawater, iron):**
  a muon that cannot cross the slab is placed at its range, which was rock's
  range times the material's ionisation factor, the inverse correction: 10
  GeV muons stopped at 42 m in 45 m of water (stepped) but were placed at
  51 m in 200 m (cut); in iron at 7.0 and 5.8 m. It is now the material's
  own mean range; rock is byte-identical.
- **Density inversion between library densities:** ρ̂ interpolated T
  linearly between the library's densities, but T falls almost
  exponentially with ρ, so the chord lies above the curve and densities
  between library values came out high: +2.6, +1.5 and +0.4 % at 1.75,
  2.25 and 2.8 g/cm³ on the Vesuvius library. The inversion, the synthetic
  T_data and the χ² landscape interpolate ln T now (+0.5, +0.3, +0.04 %).
  At a library density the result is unchanged.
- **Terrain forward transport, rock in front of each muon:** the engine
  depth of a direction bin was its slant opacity divided by cos θ where it
  should be multiplied, so every engine (UCMuon-MC, MUSIC, Bethe-Bloch,
  PROPOSAL) put 1/cos²θ too much rock in front of the muons: through a flat
  20 m slab, 2582, 1180 and 0 of 3000 muons survived at 30°, 60° and 80°
  where 2820, 2391 and 1470 do (as the Terrain tab reported them; the
  transmission libraries come from the backward driver and are not
  affected). The terrain transport now agrees with the slab within its
  statistics; the 0.5 m depth floor, which added rock to thin
  near-horizontal paths, is 1 cm.
- **Terrain transport with MUSIC, Bethe-Bloch or PROPOSAL, bins of one
  muon:** each direction bin is written to a file for the engine, which
  takes the depth axis from the start positions' plane; a single muon has
  none, so the engines took it from the momentum, horizontal beyond 45°. A
  20 GeV muon at 81° then crossed 10 m of rock instead of 64 m and survived
  (Python Bethe-Bloch), or was sent away from the exit plane and lost (MUSIC,
  Fortran Bethe-Bloch). Such a muon is now written twice, 1 cm apart, and
  the copy dropped.
- **Geant4 and guaranteed-hit PHITS exports:** they read the run's file
  through the GUI's cached loader with the path as the only key, so a second
  run to the same file (the default name) exported the first run's muons
  (200 rows after a 300-muon run) and the success message gave the stale
  count. The cache belongs to the server, so a new browser session could
  export another session's run. The file's modification time is now part of
  the key.
- **MCS margin helper (Helpers, 🎯 MCS margin):** it took the slant-path
  threshold at 2.65 g/cm³ whatever the Transport density, so in water and ice
  it evaluated the scattering at too high an energy: through 100 m of water it
  suggested 96 cm where the Generator's warning asks 376 cm (a point-like
  detector keeps about a quarter of its slowest muons with 96 cm). It now
  uses the Transport density, as the warning and Basic do.
- **Safety margin, radiation length and range:** the margin rule (Basic's
  derived margin, the Generator's warning, the helper) used Standard Rock's
  X₀ and range table for every material. In iron (X₀ 13.84 g/cm², range per
  g/cm² 12 % longer) the margin came out 11-13 % below the scattering
  UCMuon-MC produces instead of above it; in water (range 13 % shorter) it
  was a third larger than needed (376 cm through 100 m, now 281 cm). It now
  takes the Transport material's X₀ and, for water, seawater, ice and iron,
  the PDG 2024 range table UCMuon-MC uses. Standard Rock margins are
  unchanged (the paper's 2.3 m).
- **MCS acceptance estimator (Generator, after a run with a detector
  filter):** it tested the scattered muons against the margin-inflated
  detector, evaluated Highland at the surface momentum and counted the
  muons that range out, so for the Basic default detector (r 50 cm, 30 m
  deep) it expected 5620 hits of 20 000 selected muons where UCMuon-MC
  transport gave 113, and its banner advised raising the margin. It now
  counts the muons that survive (CSDA) and land in the detector itself, with
  the margin rule's σ_r (energy loss included), and estimates the hits the
  margin loses from the muons that start inside the detector and scatter
  beyond the margin: 96 expected against 113 transported; at a 2σ_r margin
  98.9 % estimated retention against 99.4 % measured.
- **CSG (PHITS) cell expressions:** `#N` complements were dropped and
  parentheses flattened, so a cavity or a detector cut out of the rock with
  `#N` read as rock (39.5 m of rock above a point 40 m deep, through a 10 m
  void, instead of 30 m), `-3 (1 : 2)` meant `(-3 1) : 2` and `#(...)` an
  intersection. The CSG transport's overburden was affected, not only the
  preview. Cells are now parsed with the PHITS/MCNP precedence and
  complements. A body the parser does not support (BOX, REC, WED, ...) was
  left out of its cells, so a cell made of one filled every point no earlier
  cell claimed, the sky included; that part of the cell is now empty, and
  the loader lists such surfaces (and positive densities, atoms per barn-cm
  read as g/cm³) after loading.
- **CSG overburden preview, single-ray check and 3D marker:** they read a
  detector position that nothing set, so they traced from the geometry's
  origin (open sky above a detector 30 m deep). They trace from a point
  inside the first selected detector cell now. The step-march sampled each
  step at its end and lost half a step at boundaries on the step grid; it
  samples the midpoints.
- **CSG runs use the selected engine:** the run applied straight-line CSDA
  energy loss whatever engine was selected (the engine, stochastic-loss and
  multiple-scattering options gave identical outputs). Each muon whose
  straight line reaches a detector cell is now transported through the rock
  in front of it (Σρ·L of the cells crossed, as Standard Rock) by the
  selected engine, which decides survival and gives the energy, direction
  and lateral displacement at the cell (UCMuon-MC muon by muon; MUSIC,
  Bethe-Bloch and PROPOSAL in 2 % depth groups). Through a 30 m layer the
  CSG run matches each engine's slab transport (survivors 3676/3679
  UCMuon-MC, 3686/3681 Bethe-Bloch, 3677/3672 MUSIC; mean energies within
  0.2 GeV). The results views got a placeholder overburden map, zero
  everywhere and "open sky" wherever no muon survived, i.e. in the shielded
  directions; the map is now traced from the detector cell. The PHITS-style
  hit dumps were written to `gui/`; they go to `output/`.
- **Terrain transport, steep muons:** muons at or beyond the direction grid's
  maximum zenith (Angular grid settings) were neither transported nor
  written, without a word (45 of 1000 muons spread to 89° with the 85°
  default). The run now warns with their number.
- **Command-line converters, cosmic electrons:** `src/converters/ucmuon_to_geant4.py`,
  `ucmuon_to_phits.py` and `ucmuon_to_phits.f90` wrote every row as a muon
  (PDG ±13, E − m_μ). A generator file has no particle column, so a
  spectrum-8 file (cosmic e±) came out as muons, 99 % of them with zero
  kinetic energy (1982 of 2000 at 10 MeV to 1 GeV). This is the only PHITS
  writer of the HPC path. They now take `--electrons` (`electrons` for the
  Fortran one: PDG ±11, the electron mass) and stop with exit code 2 when a
  row's own E and p contradict the particle being written, either way;
  `hpc/run_ucmuon_gen.sh` passes it for spectrum 8. Muon files convert byte
  for byte as before. The GUI's own exports were right.
- **MPI generator, EventIDs:** every rank numbered its muons from 1, so the
  files `hpc/run_ucmuon_gen.sh` merges repeated EventIDs (1501 distinct in
  3000 rows on 2 ranks; 68 in 200 selected muons on 3). The selected file now
  runs 1..N over the ranks, as the OpenMP generator's does, and the surface
  file interleaves the ranks' try numbers. The GUI's CSG run with MUSIC,
  Bethe-Bloch or PROPOSAL matched the engine's output back by EventID and
  failed on such a file; it numbers each depth bin's muons itself now.
- **setup.sh, failed builds:** a compile error printed "SKIP" (the word for
  absent sources), the script exited 0 and the summary counted the binary
  left in `bin/` by an earlier build as built. A failed build is now reported
  as FAIL, marked `[!]` in the summary, which also lists the generator, and
  setup.sh exits 1.
- **setup.sh, optional venv:** under `set -e`, a venv that could not be made
  (Debian/Ubuntu without `python3-venv`) or had no pip stopped the whole
  install before anything was compiled (exit 127), on every re-run. The
  PROPOSAL venv is optional and its failure is now a warning; a half-made
  venv is removed.
- **setup.sh, Engine 6:** rasterio was installed only for `$PYTHON`, while
  the GUI runs in `~/venvs/ucmuon` (run_gui.sh), and the summary checked
  `$PYTHON`: with conda or miniforge as `python3` a fresh install reported
  Engine 6 while the GUI had no rasterio. It is installed in the venv too,
  and the summary checks the interpreter the GUI uses.
- **`tools/compare_engines.py`:** Bethe-Bloch and PROPOSAL were given the
  material without the density, so with a density other than 2.65 they ran
  at 2.65 while MUSIC and UCMuon-MC ran at the given one (ρ = 2.0, 10 m:
  40.4 and 41.0 % survival against 55.8 %; now 53.6 to 56.8 % for all four).
- `hpc/run_proposal_bench.sh` found the repository from `$0`, which under
  sbatch is the spool copy of the script; it uses the submit directory now.
- `src/converters/merge_sources.py --energies 1500` looked for
  `source_1TeV.dat` (1.5 TeV truncated); TeV is used for whole TeV only.
- The Makefile relinked the MUSIC binaries on every `make` (the phony
  `data-links` was a normal prerequisite; it is order-only now).
- `run_gui.bat` started the GUI with `python` from PATH, which is not the
  interpreter `install.ps1` used when Python was found only in `C:\Python31x`;
  it now takes the first of `python`, `C:\Python31x\python.exe` and `py -3`
  that has streamlit, and says so when none has. `install.ps1` passes `FC` on
  the make command line (the Makefile's `FC = gfortran` overrode the
  environment variable it set).
- `tools/check_consistency.py` compared and screened only the files tracked
  in the public tree; untracked files that `git add -A` would publish are now
  included. `tools/bump_version.py` rewrote every `vX.Y.Z` in the paper (a
  PROPOSAL or UCMuGen version too); it now moves only the version the code
  states and lists the others.
- Terrain transport (GUI) left two files per depth bin per run in the project
  root (the bin's muons and the engine's `_stopped.dat` and `_timing.txt`);
  they go to a temporary directory that is removed.

### Changed
- Both generators warn after 1e7, 1e8 and 1e9 tries if no muon has hit a
  detector yet (a detector that no straight ray from the source can reach
  made the run go on forever without a word).
- **Backward MC (Engine 5), range table to 10 TeV:** the standard rock CSDA
  range table ended at 2 TeV (3616 m w.e.), and flux needing a higher surface
  energy was left out: at Vesuvius 106 of the 462 density-sensitive
  directions of the shipped library lost more than 1 % (up to 33 %), which
  `MURAVES_GUIDE.md` said none did. Groom's table is now continued to 10 TeV
  (6889 m w.e.) with the PDG 2024 standard rock stopping power (UCMuon-MC's;
  its range agrees with Groom's to 0.03 % at 2 TeV), integrated on from
  Groom's last point. The cut is now below 1 % up to 3500 m w.e. (it was 37 %
  at 3000) and no density-sensitive Vesuvius direction is cut. Fluxes under
  less than 500 m w.e. rise by at most 0.25 % (detector energies of 2 to
  2.5 TeV are now inside the table), under 1000 to 2000 m w.e. by about 3 %.
- The Vesuvius T_sim library (`examples/vesuvius/tsim_library/`) is rebuilt
  with the corrected ray tracer and the 10 TeV table (1059 directions cross
  rock instead of 966); the guide's synthetic inversion gives 460
  density-sensitive pixels, ρ̂ = 2.002 g/cm³ (scatter 0.022), D = 1.68. The
  example's flux-vs-elevation and transmission figures are redrawn.
- `test_run/expected/` regenerated with the 53-bit deviates and the
  accumulated scattering variance: every number is unchanged (574 of 2000,
  28.70 %, 13.562 GeV), the files and their SHA-256 sums are new.
- `benchmark/results/benchmark_summary_v2.csv` rebuilt after re-running
  UCMuon-MC and Bethe-Bloch on the monoenergetic benchmark with the scattering
  variance accumulated along the path: survival and energies are unchanged;
  mean angles against Geant4 are now −2.8, −0.8, −3.4, −2.8 and −5.0 % for
  UCMuon-MC at 1, 25, 50, 100 and 200 m (were −15.4 to −1.9 %), and −2.8 to
  −14.1 % for Bethe-Bloch. `run_mono_bench.sh` takes `ENGINES="..."` to run a
  subset.
- setup.sh and install.ps1 download PUMAS at a fixed commit (529976ef, tag
  v1.2, the code UCMuon is tested with) instead of the master branch.

### Benchmark data
- `benchmark/results/benchmark_summary_v2.csv` is rebuilt from the v3 engine
  outputs that the paper's tables use, by the new
  `benchmark/analysis/regen_summary_v2.py`. It was a June 2026 copy with only
  the UCMuon rows patched by hand: the BB mean exit energy at 200 m was 135.18
  GeV (v3: 128.50, as in the paper) and the UCMuon rows were not the v3 run
  (survival at 200 m 16.369 %, v3 16.392 %). Geant4, PHITS, MUSIC and PROPOSAL
  are unchanged; nothing the README quotes changes.

### Documentation
- `docs/ENGINE6_USAGE_GUIDE.md`: the ray-trace step no longer changes the
  result; DEM interpolation and projected DEMs explained.
- UCMuGen README: the tilted-plane law (1 + cos α)/2 is for an isotropic sky
  over the upper hemisphere and cos α for a vertical beam only (the README
  also named a full-sphere field, where the rate does not depend on the
  tilt); the PARMA header lists six edits; the documented feature-tour
  numbers are regenerated (`validation/reference/tour_numbers.txt`: rates
  move within their Monte Carlo error with the Plane axes, directed
  acceptances rise about 16 % with the exact cone).
- `gui/mcs_margin.py`: σ_r is an upper envelope of the engines' scattering,
  measured with vertical beams through 40 and 90 m of rock: 8-10 % above
  MUSIC's displacement and 22 % above UCMuon-MC's, so 2σ_r keeps at least
  98 % of a point-like detector's hits with every engine.
- README and the UCMuon-MC driver docstrings: hard radiative events are
  sampled per process for bremsstrahlung and photonuclear interactions only;
  pair production stays in the continuous term (as the code has done since
  June 2026). They described pair production as sampled too.
- `src/generator/rng_parallel.f90`: a fixed `UCMUON_SEED` reproduces a run
  exactly with one OpenMP thread only (the comment said any run).
- `benchmark/codes/geant4/README.md` and `benchmark/reports/BENCHMARK_FEEDBACK.md`:
  the Geant4 reference uses `FTFP_BERT` with its default electromagnetic
  physics, not `G4EmStandardPhysics_option4` (`MuonRock.cc` replaces nothing).
- `benchmark/README.md` and `benchmark/codes/{music,bethe_bloch,proposal,ucmuon_stochastic}/README.md`
  gave a `compare_engines.py --engine ... --depths ...` command that no
  script accepts; the benchmark outputs come from
  `benchmark/analysis/run_mono_bench.sh`, which runs the Python Bethe-Bloch
  driver (as the paper says), not the Fortran binary.
- `hpc/README_HPC.md`: `logs/` has to exist before the first sbatch (SLURM
  opens the log files before the job script runs; it said "auto-created");
  the converter's `electrons` option.
- `run_gui.sh` called PROPOSAL Engine 3; it is Engine 4.
- `benchmark/codes/bethe_bloch/README.md` and `ucmuon_stochastic/README.md`:
  the notes on agreement with Geant4 were from older runs ("±2 % energy at all
  depths", "MCS angles −40 %", "Groom dE/dx", "lowest survivor"); they now
  give the measured values.

### Tests
- Two regression tests (`tests/gui/test_fixes_v133.py`), each failing on 1.3.2.
- `tests/flux/test_generator_fixes.py` (PARMA date, MPI parity; in CI, which
  now also runs on changes to `src/parma/`) and a PHITS-dump test in
  `tests/gui/test_fixes_v133.py`, each failing on 1.3.2.
- `tests/transport/test_transport_fixes.py` (Bethe-Bloch water, pencil beams;
  in the test-run CI job) and three transport tests in
  `tests/gui/test_fixes_v133.py`, each failing on 1.3.2.
- Multiple scattering against Highland for one radiation length, at two step
  sizes (UCMuon-MC and Python Bethe-Bloch in `tests/gui/test_fixes_v133.py`,
  Fortran OpenMP and MPI in `tests/transport/test_transport_fixes.py`), and
  the resolution of the random deviates of the Fortran generator and UCMuGen
  (`tests/flux/test_generator_fixes.py`), each failing on 1.3.2.
- Six terrain and density tests in `tests/gui/test_fixes_v133.py` (flat-slab
  ray tracing, projected and ASC DEMs, backward MC materials against the
  PDG tables, full-precision transmission maps, STL from inside the mesh,
  θ_max advice), each failing on 1.3.2.
- Ten GUI tests in `tests/gui/test_fixes_v133.py` (Geant4 export after a
  second run, the margin helper and Basic margin for water and iron, CSG
  complements, the preview position and the engine transport, the MCS
  estimator, the terrain zenith warning, the terrain depth against a flat
  slab, a terrain bin of one muon), one for the density inversion between
  library densities and one for range-cut stopping depths, each failing on
  1.3.2.
- UCMuGen: `validation/test_parma.cc` section 7 (W index of dates after the
  table), `validation/test_detector.cc` sections 5 and 6 (detector close under
  a large surface, setSurface, Plane axes, energy range), each failing on
  UCMuGen 0.2.
- `tests/converters/test_converters.py` (in the test-run CI job): cosmic
  electrons and muons through the three converters and both refusals; the
  electron checks fail on 1.3.2. MPI EventIDs in
  `tests/flux/test_generator_fixes.py` and a CSG run with repeated EventIDs in
  `tests/gui/test_fixes_v133.py`, each failing on 1.3.2.

### Benchmark analysis
- `benchmark/analysis/benchmark_analysis.py` took an engine's deflection
  angle as the arccos of a dot product of direction cosines written with 6
  decimals, so angles below about 1 mrad came out as 0 or 1.4 mrad: the
  engines' mean angles were 4 to 6 % low at 1 m and 0.4 to 1 % low at depth
  (Geant4's, from full-precision momenta, were not affected). It now uses
  atan2(|s × c|, s · c).

## [1.3.2] — 2026-10-08

A first run with the GUI's defaults now gives a meaningful result. No physics
or numerical result changes.

### Fixed
- **Default detector.** It was a 90 m tall column of 5 cm radius reaching up to
  the ground. In Basic mode, where the overburden is the depth of the
  detector's top face, a first run with the defaults transported through no
  rock (the Transport tab said so) after generating 100 000 hits for about
  13 minutes. The default is now a vertical cylinder of 0.5 m radius and 1 m
  height (box: 1 m x 1 m x 1 m) whose top face is 30 m deep, and the generator
  asks for 10 000 muons: a first run takes 3 to 4 minutes on 4 threads.
  Saved sessions keep their own detector.
- Advanced mode's guaranteed-hit workflow also defaults to 10 000 muons (it
  shares the field with the Standard workflow and had kept 100 000).
- A detector face at the ground printed as "-0 m".
- The Transport tab numbers PUMAS ⑦, as the README and the documentation do
  (⑥ is UCMuon Terrain, in its own tab).
- Results: a rate's statistical error is a caption under it ("± 1.24 /s
  (statistical, 1σ)"); it was the metric's delta, drawn as a green up arrow
  as if it were a change. Fig. 4 of the paper is retaken with it.
- The UCMuon-MC engine description says it runs one worker per ~20 000 muons,
  so a small run is serial (the Mode card said "Serial Python" next to
  "Multiprocess-parallel").
- `ucmugen/include/UCMuGen_PARMA.h` compiles under clang with `-Wpedantic
  -Werror`: its longest table string (139 000 characters) is above the length
  a compiler must accept, and the warning is now silenced around the tables
  (`ucmugen/tools/make_parma_header.py`; the regenerated header differs only
  by those lines). The UCMuGen CI now builds `test_parma.cc` with `-Werror`
  (clang warns only when the tables are compiled in, not on the include).

### Documentation
- README: a "Your first simulation (Basic mode)" walkthrough after the
  installation, with the times and numbers of measured runs and why only a
  few tens of the hits end in the detector; which tabs each mode shows; the
  engine count per tab and the engines in their ① to ⑦ order; which parts need
  gfortran (the surface generator, so Basic cannot generate without it) and
  which do not; `data/` does not ship the MUSIC tables; `tools/` and `docs/`
  list all their files; the Terrain tab's seven sections, section 4 (DEM
  check) included.
- Engine numbers corrected in `requirements.txt` (and `external/pumas-master/`
  for PUMAS), `INSTALL.md` (the installer's summary lists Engine 7) and the
  `benchmark/codes/` READMEs of MUSIC, Bethe-Bloch and PROPOSAL.
- `ROADMAP.md`, the working roadmap of June 2026 (it pointed to files that are
  not published and described the GUI as it was then), is no longer
  published; `tools/check_consistency.py` reports it if it comes back.

### Tests
- Four GUI regression tests, each failing on 1.3.1: fresh-install detector
  underground with 10 000 hits; no "-0"; guaranteed-hit default 10 000; the
  rate's error shown under it, not as a delta.

## [1.3.1] — 2026-10-08

Two GUI Results-tab bugs were found after 1.3.0, while re-measuring the paper's
examples and taking its GUI figure. A full audit followed (every number and
plot of the Results, Transport, Generator, Terrain and Density tabs recomputed
independently; every engine's input and output; the command-line, HPC and
documentation paths), and everything it found is fixed here.

### ⚠ Changes results
- **Rate at depth of the detector-hit file.** After transport, the
  underground detector filter tested each muon against the detector inflated
  by the generator's safety margin, so the detector-hit file, and the rate at
  depth that the Results tab gives it since 1.3.0, were those of the inflated
  volume. The margin belongs to the generator's straight-line cut (it keeps
  the muons that multiple scattering moves into the detector); underground,
  the real detector decides. Measured on a 10 cm × 50 cm cylinder with its top
  face 90 m deep and a 2.3 m margin (Guan, angular mode 6, MUSIC,
  4 × 10⁵ generated hits): 810 ± 3 muons/min before, 4.5 ± 0.2 muons/min now.
  Rerun the transport for detector-hit files made with 1.3.0 or earlier; the
  underground file itself is unchanged.
- **Fast flux estimator range table** (GUI fast-estimator panels, Density
  tab Methods 2 and 3, the margin helper): the CSDA ranges from 140 GeV on did
  not match Groom (2001) Table IV-6 (1 % long at 140 GeV rising to 11 % at
  1 TeV), and the table ended at 1 TeV, beyond which the flux was 0. Now the
  Groom values to 2 TeV (as the GUI's own table and the backward MC), then
  dE/dX = a + bE. The slab flux under 200 / 500 / 800 m of rock drops by
  about 4 / 12 / 20 %; Method 2 on the shipped library (ρ = 2.00) moved from
  ρ̂ = 2.115 to 2.049 (median; usable pixels 535 to 642). Below 100 GeV (≈ 150 m of rock) nothing changes.
- **Terrain tab, per-muon direction**: muons were binned by their direction of
  travel instead of their arrival direction, so each crossed the rock on the
  opposite side of the sky. Per-bin transmissions were right on average
  (every spectrum is uniform in azimuth); the per-muon alive flags in the
  terrain output, and anything using survivor directions, were mirrored.
- **Results, legacy angular modes (1-5)**: the mean energy, the survival
  fraction and the charge ratio were unweighted while the rate was weighted
  (Guan mode 2 at the surface: mean E 6.2 GeV shown, 8.3 weighted).
- **Results, survival vs depth**: the CSDA curve used the vertical opacity
  for every muon; it now uses each muon's slant path (all-sky, 20 m: 17.7 %
  before, 11.9 % now, equal to Bethe-Bloch's Monte Carlo).
- **Exposure panel**: R came from a separate flux-model selector (default
  Reyna) whatever the generator's spectrum (−15 % for a Guan run); it now
  uses the generator's spectrum (1, 4-7), the last run's printed R for PARMA,
  and the selector only for spectra 2 and 8.
- **PUMAS** (Engine 7): the backward flux spectrum dΦ/dE was Φ·E·ln(Emax/Emin)
  (each bin averaged its events' Φ/pdf instead of summing them over N·ΔE); the
  backward-event metrics (mean energies, μ⁺/μ⁻) were unweighted sampler
  statistics. Forward mode aborted on ordinary generator files ("bad norm for
  state direction": directions from 6-decimal momenta were not normalised),
  and read the detector mask of a 13-column hits file as its hit flag,
  dropping the hits of detectors 2, 3, ...
- **Bethe-Bloch, Python driver**: the Highland momentum treated the kinetic
  energy as total: NaN directions below 105.7 MeV (written alive), too much
  scattering above (about ×2 at 200 MeV). Survival and exit energies do not
  depend on it (unchanged, checked event by event).
- **Bethe-Bloch, MPI binary, custom material**: the radiation length was half
  Tsai's ((11.319 − ln Z) instead of ln(287/√Z)), so multiple scattering was
  √2 too large.
- **Geant4 and PHITS exports from the GUI**: cosmic e± (spectrum 8) were
  written as muons (PDG ±13, muon mass); and with the detector filter on and
  Save ALL off, the Geant4 export converted an all-muons file left by an
  earlier run.
- **Terrain tab**: Save T_sim, the 3D view and the cross-check used the
  density field's current value instead of the run's (a T_sim file could be
  saved under the wrong density); the flux-vs-elevation plot stayed on screen,
  relabelled, after the azimuth or the run changed.
- **Density tab**: library pixels that are NaN (bins of a Terrain-saved T_sim
  file with no muons) came out "OK" with ρ̂ = NaN, and the synthetic data
  generator crashed on them; downloads and the Methods 2-3 maps could be
  those of the previous inversion; a "synthetic validation" banner stayed
  after loading measured data.
- **Backward MC panel**: the mean surface energy and survival probability
  were averaged over the full solid angle, counting directions that do not
  contribute as 0 (at 1000 m: P_survival 0.39 shown, 1.00 over the
  contributing directions). The rate was right.
- **Guaranteed-hit mode, spectrum 2** sampled p^-2.7 instead of p^-3.7.
- **Source-size helper**: its "optimised radius" left out the depth × tan θ
  reach (127× too few hits in one case); it now gives the required source of
  Basic mode and the Generator's warning.

### Fixed
- Results, distributions of an underground file: the histogram of an
  at-depth quantity (θ, φ, x, y, z, E, direction cosines) and the θ-φ
  acceptance map included the stopped muons, which the engines write with
  zeros there, so for a deep run nearly every row sat in the first bin. They
  now use the survivors; the surface columns keep every row. The sea-level
  "cos²θ·sinθ peak" line, which held for none of the files, is gone.
- Results: the underground charge composition counts survivors; the detector
  solid angle is that of the detector, not of the margin-inflated volume; the
  rate of a generator detector-selection file says it is that of the inflated
  volume; a detector-hit file keeps the survivor count and the underground
  file it was cut from (a later transport or Terrain run overwrote them), and
  without them shows no survival curve; live times under a minute keep three
  digits (1.27 s, not "1 s"); the 3D title names the file kind.
- Transport: a failed MUSIC or Fortran Bethe-Bloch run returned exit status 0
  (now 1), and the underground filter then ran on the previous run's file (it
  now skips a file older than the run); the depth-mismatch warning only for
  the generator run's own output; the Iron preset says which engines use
  rock composition at iron density; the X₀ input says which engines use it;
  the filter caption no longer shows the margin.
- Generator: the directional-flux and MCS-margin helpers took the total E_min
  as kinetic (≤ 4 %); help texts (angular-mode bias, panel acceptance) and
  rounding.
- Terrain/Backward MC: "Underground detector", v_cut and MS are autosaved (a
  restart turned an underground detector into a surface one); outputs name
  the backward-MC spectrum (its numbering differs from the generator's); the
  backward-MC panel warns below a spectrum's fitted range; labels (slant, not
  vertical, path; nearest-pixel DEM lookup; the Density altitude factor
  cancels in T).
- Documentation: `hpc/` inputs (rectangle template 7 lines; X₀ in g/cm²;
  density note; these inputs are for the MPI binaries), `hpc/README_HPC.md`
  (paths, 13-column hits file), README and INSTALL (gfortran for the
  generator, engine ordering, Terrain walkthrough), `docs/ENGINE6_USAGE_GUIDE.md`
  (driver name, spectrum, typical fluxes, runtimes, GUI steps), the MURAVES
  guide, Makefile engine numbers, converter and driver docstrings.
- `hpc/input_params.dat` called the power-law spectrum E^-2.7; it is E^-3.7.
- Stopped muons in the 18-column output: every engine now writes where the
  muon stopped (x, y, z). Up to 1.3.0 MUSIC did, the Python engines
  (UCMuon-MC, Bethe-Bloch, PROPOSAL) wrote the surface x, y with the stopping
  depth, and PUMAS and the Fortran Bethe-Bloch range cut wrote the surface
  point. (Their direction columns are placeholders, as before.)
- PUMAS: live progress in the GUI (the driver held the binary's output until
  the end).
- Generator console: progress lines without the detector filter report the
  same tried count as the summary; "Integrated flux" gives its units
  (vertical, cm⁻² s⁻¹ sr⁻¹); the summary names the detector-hits file when one
  is written; "Plane at W" for the source plane position.
- GUI: the Terrain and Density inputs (method, file paths, download box,
  thresholds, detector layout, ...) are autosaved, so Basic-mode reruns no
  longer reset them; the survival-plot subtitle fits.

### Tests
- `tests/gui`: the underground filter keeps only muons that reach the real
  detector; the 3D title; the range table against Groom; the PUMAS spectrum
  integrates to its rate; Bethe-Bloch directions finite; Terrain arrival
  azimuth; NaN library pixels; stopped-muon positions (Bethe-Bloch,
  UCMuon-MC); Terrain and Density inputs surviving Basic reruns. Each fails
  on 1.3.0.

## [1.3.0] — 2026-10-06

The GUI gets a **Basic / Advanced** switch; before it, the GUI bugs found by an
inventory of every GUI control (B9-B15, and two found while checking them) are
fixed.

### Added
- **Basic / Advanced mode** (an "Advanced mode" toggle at the right of the
  page title). Basic has the Generator, Transport and Results tabs (Terrain,
  Density and Config are Advanced only) and the essential inputs: spectrum
  (all eight; PARMA with its site and date), energy range (choosing a
  spectrum fills in its recommended range, also θ_max 60° for Gaisser, and a
  button restores it), angular distribution (all six modes, ⑥ by default) and
  θ_max, source type, an optional detector (on by default), the output file
  names with the PHITS (s-type=17) and Geant4 (ASCII or HEPEvt) exports, the
  OpenMP threads, the transport input (the generator's output or the path of
  any muon file) and output, the transport engine (not the Backward MC flux
  integrator, which needs no muon file) and the overburden material; Results
  has every plot including the 3D trajectories (the overlay editor, which
  rewrites the run's detector record, stays Advanced). With the detector on, Basic derives and shows the
  safety margin (2σ_r of multiple scattering at the detector's deepest
  point), the source size (just large enough for every straight path at
  θ ≤ θ_max; measured: the hit rate does not change when it is doubled) and
  the overburden from the detector's top face (when the transport input is
  that run's output); with it off, or for any other input file, the user
  gives the source size and the overburden depth. Basic does not change the
  Advanced-only settings: the autosaved ones come back when switching (and
  after a restart), and after a Basic run Advanced shows the values it used,
  with a button to copy them. A fresh install opens in Basic, an existing
  autosave in Advanced.
- Terrain: "Bundled Mt. Vesuvius" DEM source (also sets the detector site),
  the default for a fresh session; "Altitude from the DEM" button. Density:
  the shipped Vesuvius T_sim library is the default when the project root has
  none of its own.
- `tests/gui/` (pytest + Streamlit AppTest) and the `GUI` CI workflow: Basic and
  Advanced mode, hidden settings surviving, the bug fixes below; latest
  Streamlit without rasterio, and Python 3.9 + Streamlit 1.49 (the minimum)
  with rasterio.

### Changed
- **Streamlit 1.49 or newer** (`requirements.txt` said 1.35): the GUI uses
  `width="stretch"` on buttons and tables: up to 1.47 `st.button` rejects
  it, 1.48 fails on it in `st.dataframe`. Measured by running `tests/gui`
  against 1.36, 1.40, 1.44, 1.46, 1.47, 1.48 (all fail) and 1.49 (passes).
- Terrain is one page: sections 1-7 from top to bottom (surface muon file,
  geometry, detector, DEM check, engine and material, overburden preview,
  run), then the results; up to v1.2.0 they were split over three sub-tabs
  (Setup, Overburden map, Run & Results) that referred to each other. The
  overburden preview and the DEM check's elevation profile and 3D view are
  collapsed; the DEM check's position and altitude verdict stays visible.
- Fresh sessions start on Guan 2015 over 1-2500 GeV (was CosmoALEPH from
  100 GeV, which drops every muon below 100 GeV).
- The Transport depth check (and the warning after the underground filter)
  compares the overburden with the detector's top face (the rock above it);
  it asked for the deepest face, which adds the detector's height as rock.
- Settings hidden by a condition (for example the PARMA site while another
  spectrum is chosen) keep their values; Streamlit used to drop them, and the
  autosave with them. Some settings that were never autosaved now are
  (transport output file, Bethe-Bloch menu, Geant4 options, custom
  compositions).

### ⚠ Changes results
- **PARMA (spectrum 3), solar modulation from the date**: the generators read
  the W index of the date (`getHP`) and discarded it, so every PARMA run used
  the typed W (GUI default 0, solar minimum) whatever its date. The GUI now
  sends "W from the date" by default (an explicit W behaves as before). After
  the end of PARMA's neutron-monitor table (2026-03-26) the latest earlier day
  is used, with a warning, instead of a silent W = 0. Measured at sea level
  (50.7°N, 4.4°E, 1-1000 GeV, θ ≤ 85°): R −5.0 % at 2024-10-01 (W 100.5),
  −5.6 % for dates after the table (W 107.4). UCMuGen's
  `parma::w_index_from_date` is the stock JAEA routine and returns 0 after the
  table end.
- **Bethe-Bloch and PROPOSAL density**: both ran at the table density of their
  own material menu (Bethe-Bloch 2.65 / 0.917 / 1.00 / 2.30, PROPOSAL its
  built-in medium) and ignored the Transport tab's ρ, and the Terrain tab's,
  which converts each bin's opacity to a depth with it. The material line of
  both stdin protocols takes an optional density (`1 2.71`); without it the
  output is unchanged (checked bit-identical).
- **Results plots for legacy angular modes 1-5**: energy, variable and θ-φ
  histograms were raw counts (the sampling distribution); they are now
  weighted to the true flux. PUMAS backward-event histograms are weighted by
  `flux_contribution`.

### Fixed
- Density σ (Gaussian prior, KDE) reaches no engine; the tab now says that
  only the mean ρ is transported. Switching to Gaussian no longer draws a
  zero-width pdf.
- Terrain: the "Surface muon file" choice was ignored (the run took the first
  candidate); a synthetic or CSG geometry kept overriding the DEM after
  switching mode or pressing Clear, and CSG Clear reset the DEM path; result
  captions named Puy de Dôme and "CosmoALEPH, 1293 m" for every site.
- Transport and Terrain preselected a remembered input file over the newest
  generator output.
- Results: underground and detector-hit files show rate and live time as
  headline metrics; the tab opens on the newest output of the session (the
  underground file after a transport run); the detector-hit file from the
  underground filter always said "No rate" (its source was not recorded);
  an unset session path could appear as an empty file choice.
- Results, survival vs depth: for the detector-hit underground file the
  CSDA curve was drawn from those hits only (all above threshold, so near
  100 %) and the plot said "MUSIC" for every engine; it now uses every
  transported muon and names the engine. The transport's own survival count
  is used only for that transport's files.
- Generator: the "Will write" line also names the PHITS and Geant4 files.
- Config: the downloaded config JSON restored nothing ("Restored 3
  settings"); it now carries the settings Restore reads, and an old
  summary-only file gets a clear message.
- Density: the plain Gaisser formula is no longer offered for the direct and
  two-flux inversions (it integrates from 0.5 GeV, ×12 high there); measured
  T_data without a `# Density:` header can be loaded; measured data can get
  σ_T from counting statistics, so the σ_ρ map and χ² work for them.
- Transport: the OpenMP slider is shown only for the engines it reaches
  (MUSIC, Fortran Bethe-Bloch), and is autosaved; UCMuon-MC's mode shows its
  worker processes, not "Single-thread".

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
