# UCMuon Terrain Engine — Complete Usage Guide
# Engine 6: DEM-aware muon flux for real field deployments
# UCLouvain Muography Group — Hamid Basiri <hamid.basiri@uclouvain.be>
# ─────────────────────────────────────────────────────────────────────────────

================================================================================
WHAT THIS ENGINE DOES AND WHY IT IS DIFFERENT
================================================================================

Every other engine in UCMuon (UCMuon-MC, MUSIC, Bethe-Bloch, PROPOSAL,
Backward MC) assumes a FLAT SLAB: one depth, one density, every direction from
the detector sees the same overburden.

This is a fine approximation when:
  - The detector is inside a mine or tunnel with known flat rock above
  - The site is flat and you are not close to any hills

It breaks down completely when:
  - The detector is on a volcano flank — one side sees 2000 m of rock,
    the opposite side sees open sky
  - The detector is near a glacier margin — varying ice thickness per direction
  - The detector is in a valley between two ridges — heavily asymmetric screening
  - You are planning a CCS monitoring array in hilly terrain
  - Any experiment where the background rate depends on azimuth

UCMuon Terrain solves this by:
  1. Reading a real terrain file (GeoTIFF DEM, e.g. SRTM 30m or 90m)
  2. For EVERY (azimuth, zenith) direction bin, tracing a ray backward
     from the detector through the terrain
  3. Computing the actual rock overburden along each ray
  4. Returning a full 2D directional flux map [m⁻² s⁻¹ sr⁻¹]

The output replaces the single "survival rate" number of other engines with
a complete picture of how terrain blocks muons from each direction.


================================================================================
STEP 1 — INSTALL RASTERIO
================================================================================

rasterio is the only new dependency. It is the standard Python library for
reading GeoTIFF files (the same format used by QGIS, GDAL, and all GIS tools).

  pip install rasterio

On Linux you may also need GDAL headers if pip cannot find a pre-built wheel:

  # Ubuntu / Debian:
  sudo apt install libgdal-dev python3-dev
  pip install rasterio

  # CentOS / Rocky:
  sudo yum install gdal-devel python3-devel
  pip install rasterio

  # macOS (Homebrew):
  brew install gdal
  pip install rasterio

  # macOS (MacPorts):
  sudo port install gdal
  pip install rasterio

Verify the installation:

  python -c "import rasterio; print(rasterio.__version__)"
  # Should print: 1.3.x or higher


================================================================================
STEP 2 — GET A DEM FILE (GeoTIFF)
================================================================================

You need a Digital Elevation Model file covering the area around your detector.
Three options, all free:

────────────────────────────────────────────────────────────────────────────────
OPTION A — Auto-download inside the GUI  (easiest)
────────────────────────────────────────────────────────────────────────────────

In the GUI under Tab 2 → UCMuon Terrain → Section 1 → "Auto-download" tab:

  1. Set the bounding box around your site (south/north/west/east lat/lon)
     — typically ±0.5° around the detector is enough for most sites
  2. Choose product: "SRTM GL1 (30m)" for best resolution
  3. Click "Download DEM"
  4. The file is saved locally and automatically loaded

This uses the OpenTopography public REST API (free, no account needed).
Rate limit: ~5 requests/day with the demo key.

────────────────────────────────────────────────────────────────────────────────
OPTION B — Download from OpenTopography website  (best control)
────────────────────────────────────────────────────────────────────────────────

  1. Go to https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1
  2. Draw a bounding box around your site
  3. Select "SRTM GL1 (30m resolution)" or "SRTM GL3 (90m, larger area)"
  4. Click "Export" → select "GeoTiff" format
  5. Download the .tif file

No account required for SRTM data. Registration is needed for some higher-res
products.

Recommended bounding box size:
  - Flat terrain:     ±0.3° around the detector (≈33 km × 33 km at 50°N)
  - Hilly terrain:    ±0.7° to capture all ridges that may block muons
  - Volcano flanks:   ±1.0° or larger to include the full edifice
  - Alpine glaciers:  ±0.5° typically sufficient

────────────────────────────────────────────────────────────────────────────────
OPTION C — Download from USGS EarthExplorer  (most complete archive)
────────────────────────────────────────────────────────────────────────────────

  1. Go to https://earthexplorer.usgs.gov/
  2. Create a free account (required for download)
  3. Enter your coordinates in the search
  4. Dataset tab → Digital Elevation → SRTM 1 Arc-Second Global
  5. Results → Download → GeoTIFF

This gives you the highest quality SRTM data directly from NASA/USGS.

────────────────────────────────────────────────────────────────────────────────
OPTION D — eio command-line tool  (scripted download)
────────────────────────────────────────────────────────────────────────────────

  pip install elevation gdal

  # Download SRTM 30m DEM for a 1°×1° box around Louvain-la-Neuve:
  eio clip -o louvain_dem.tif --bounds 4.0 50.2 5.2 51.0

  # Download SRTM 90m DEM for a larger area (e.g. Puy de Dôme, France):
  eio --product SRTM3 clip -o puydedome_dem.tif --bounds 2.5 45.5 3.5 46.5

  # For Etna volcano, Italy:
  eio clip -o etna_dem.tif --bounds 14.8 37.5 15.2 37.9

The resulting .tif file can then be uploaded or pointed to in the GUI.


================================================================================
STEP 3 — CONFIGURE THE DETECTOR POSITION
================================================================================

In Section 2 of the Terrain panel, enter:

  Latitude  [°N]     : decimal degrees, WGS84 (same system as GPS / Google Maps)
  Longitude [°E]     : decimal degrees, WGS84  (negative = West)
  Altitude  [m asl]  : metres above sea level (NOT above ground)

HOW TO GET COORDINATES:
  - Google Maps: right-click on the detector location → "What's here?"
                 The coordinates appear at the bottom (lat, lon).
  - GPS receiver: ensure WGS84 datum is selected
  - QGIS: open the DEM, hover the cursor over the detector site

IMPORTANT — ALTITUDE:
  The altitude must be the detector's actual elevation above sea level, NOT
  above the local ground. This is the same value shown by a GPS device.
  For a detector inside a tunnel or borehole, subtract the depth:
    alt_detector = alt_surface - depth_underground

EXAMPLE SITES:
  UCLouvain campus:          lat=50.6686, lon=4.6158,  alt=90 m
  Puy de Dôme summit:        lat=45.7716, lon=2.9645,  alt=1465 m
  Etna NE crater flank:      lat=37.7481, lon=15.0158, alt=2500 m
  Mont Blanc tunnel entry:   lat=45.8662, lon=6.8679,  alt=1275 m
  La Palma (CCS pilot site): lat=28.5740, lon=-17.8475, alt=400 m


================================================================================
STEP 4 — SET MATERIAL AND PHYSICS PARAMETERS
================================================================================

Rock density ρ [g/cm³]
  This density is applied uniformly along all rock paths. Use the mean
  density of the overburden material for your site.

  Standard Rock (benchmark):   2.65 g/cm³
  Granite:                      2.70 g/cm³
  Limestone / Chalk:            2.70 g/cm³
  Basalt (volcano):             2.85–3.00 g/cm³
  Ice (glacier):                0.917 g/cm³
  Saturated sandstone:          2.30 g/cm³
  Carbonates (CCS reservoir):   2.50–2.70 g/cm³

  For heterogeneous geology, use the density-weighted average:
    ρ_eff = Σ (ρᵢ × hᵢ) / Σ hᵢ

Surface spectrum model
  1 = CosmoALEPH (Schmelling 2013) — default; best for thick targets (E ≳ 50 GeV)
  3 = Guan et al. (2015)        — use for comparison / publication
  4 = Frosin et al. (2025)      — newest re-fitted parametrisation

Survival probability mode
  "CSDA + stochastic" (recommended) — includes the Poisson correction for
  catastrophic radiative losses. More accurate at overburden > 200 m.w.e.
  "CSDA only" — faster, slightly optimistic at large overburdens.

Grid settings
  Azimuth bins (n_az):  36 (10° steps) is standard. Use 72 for publication quality.
  Zenith bins  (n_ze):  18 ( 5° steps) is standard. Use 36 for publication quality.
  Max zenith (°):       75° is recommended. Beyond 80° the flat-Earth approximation
                        degrades and computation slows significantly.
  Ray-trace step (m):   50 m gives good accuracy. Use 100 m for fast preview,
                        20 m for high-precision near-horizontal directions.

Approximate runtimes (36 az × 18 ze = 648 directions):
  Preview  (n_az=36, n_ze=9,  step=100m): ~30 s
  Standard (n_az=36, n_ze=18, step=50m):  ~3 min
  Fine     (n_az=72, n_ze=36, step=20m):  ~25 min


================================================================================
STEP 5 — RUN THE ENGINE
================================================================================

IN THE GUI:
  1. Load DEM (Section 1)
  2. Set detector position (Section 2)
  3. Configure material & physics (Section 3)
  4. Optional: click "Quick terrain preview" to see the overburden shape
     immediately (coarse, ~30s) before committing to the full run
  5. Click "▶ Run UCMuon Terrain"
  6. Progress appears in the live console panel below the button
  7. Two polar heatmaps appear when the run completes

FROM THE COMMAND LINE (for HPC or scripting):

  echo "path/to/dem.tif
  50.6686
  4.6158
  90.0
  2.65
  1
  36
  18
  75.0
  50.0
  1
  terrain_overburden.dat
  terrain_flux.dat
  terrain_summary.dat" | python gui/cosmoaleph_terrain_driver.py


================================================================================
STEP 6 — INTERPRET THE RESULTS
================================================================================

The engine produces three output files:

────────────────────────────────────────────────────────────────────────────────
terrain_overburden.dat
────────────────────────────────────────────────────────────────────────────────
Columns: azimuth[deg]  zenith[deg]  overburden[g/cm²]  open_sky(0/1)

This is the core DEM-derived product. Each row gives the rock thickness
(in g/cm² = slant path × density) along one line of sight.

  overburden = 0          → open sky direction (no terrain in this direction)
  overburden = 26,500     → 100 m of standard rock (typical detector depth)
  overburden = 265,000    → 1000 m of rock (deep underground)

To convert to metres of rock:
  depth_m = overburden_gcm2 / (rho_gcm3 × 100)

To convert to metres water equivalent (m.w.e.):
  mwe = overburden_gcm2 / 100   [since water density = 1.0 g/cm³]

────────────────────────────────────────────────────────────────────────────────
terrain_flux.dat
────────────────────────────────────────────────────────────────────────────────
Columns: azimuth[deg]  zenith[deg]  flux[m⁻² s⁻¹ sr⁻¹]

Expected muon flux per solid angle at the detector for each direction.
The total rate in m⁻² s⁻¹ is the sum × solid angle element.

Typical values:
  Open sky, vertical:           ~170 m⁻² s⁻¹ sr⁻¹
  Open sky, ze=60°:             ~40 m⁻² s⁻¹ sr⁻¹
  100m rock overburden:         ~1×10⁻³ m⁻² s⁻¹ sr⁻¹
  1000m rock overburden:        ~1×10⁻⁸ m⁻² s⁻¹ sr⁻¹

────────────────────────────────────────────────────────────────────────────────
terrain_summary.dat
────────────────────────────────────────────────────────────────────────────────
Human-readable summary: number of rock/sky directions, median overburden,
max overburden direction, total expected rate, peak flux direction.

────────────────────────────────────────────────────────────────────────────────
GUI Polar Heatmaps
────────────────────────────────────────────────────────────────────────────────
Two interactive polar plots appear after the run:

  Overburden map  (azimuth × zenith, colour = g/cm²):
    Centre = vertical (ze=0°), edge = maximum zenith angle
    0° at top = North, 90° = East (geographic convention)
    Colours indicate rock thickness — dark = thin/no rock, bright = thick rock
    Masked areas = open sky

  Flux map  (azimuth × zenith, colour = log₁₀ flux):
    Same geometry. Bright = high flux (open sky or thin overburden)
    Dark = blocked directions (thick terrain)


================================================================================
WORKED EXAMPLE — Puy de Dôme volcano, France
================================================================================

Site: Detector on the south slope of Puy de Dôme, aimed at the summit.
Goal: Compute expected flux to plan a 2-week muography campaign.

Parameters:
  DEM:        SRTM GL1 30m, bounds: 2.8° W to 3.2° E, 45.5°N to 46.0°N
  Detector:   lat=45.760, lon=2.947, alt=1200 m  (south slope)
  ρ:          2.70 g/cm³ (volcanic basalt)
  Spectrum:   Guan et al. 2015 (mode 3)
  Grid:       72 az × 36 ze, ze_max=80°, step=30m  (fine, publication quality)

Expected result:
  - North sector (toward summit):   overburden 50,000–150,000 g/cm²
    → flux ~10⁻⁵–10⁻⁶ m⁻² s⁻¹ sr⁻¹
  - South sector (open valley):     overburden ~0 (open sky)
    → flux ~100–150 m⁻² s⁻¹ sr⁻¹  at ze<30°
  - Pronounced azimuthal asymmetry visible in the polar heatmap
  - Total rate in the "volcano direction" (±30° in azimuth, 40°–70° ze):
    → ~0.01–0.1 m⁻² s⁻¹  (highly dependent on exact geometry)

To estimate the number of muons detected in time T [s] with detector area A [m²]:
  N = Rate [m⁻² s⁻¹] × A × T × solid_angle_fraction


================================================================================
DIFFERENCES FROM A FLAT-SLAB CALCULATION
================================================================================

To understand the added value, compare:

Flat-slab (all other engines):
  Assumes all directions see the same overburden = depth_m × rho × 100 g/cm²
  → Single survival rate. Ignores asymmetry.
  → Overestimates flux from blocked directions.
  → Underestimates flux from open-sky directions on steep slopes.

UCMuon Terrain:
  Computes overburden per direction from actual DEM.
  → Correct azimuthal asymmetry.
  → Blocked directions (ridges, summits) get correct suppression.
  → Open directions get correct open-sky flux.
  → Total rate can differ from flat-slab by factors of 2–10× in hilly terrain.

Quantitative comparison (example, Puy de Dôme south slope):
  Flat-slab (assuming mean depth=200m):   total rate ≈ 2×10⁻³ m⁻² s⁻¹
  UCMuon Terrain (actual topography):     total rate ≈ 8×10⁻³ m⁻² s⁻¹
  → Flat-slab underestimates by 4× because it ignores the open south sector


================================================================================
KNOWN LIMITATIONS AND FUTURE IMPROVEMENTS
================================================================================

1. Flat-Earth approximation
   The ray tracing uses a simple flat-Earth geometry (ENU coordinates with
   constant lat/lon-per-metre scale factors). Valid for distances < ~50 km.
   For sites requiring rays > 50 km (very flat terrain, large zenith angles),
   accuracy degrades. Full spherical Earth ray tracing would improve this.

2. Single uniform density
   Rock density is applied uniformly along each ray. For sites with known
   density stratification (e.g. ice over rock, sediments over basement),
   a future multi-layer extension would improve accuracy.

3. DEM resolution vs. ray-trace step
   The ray-trace step should be ≥ 2× the DEM pixel size (typically 30m for
   SRTM GL1). Using step_m < 15m with a 30m DEM gains nothing and is slow.

4. No scattering through terrain
   The engine uses the backward CSDA flux formula which ignores Molière
   scattering through rock. At very large overburdens (>1000 m.w.e.), the
   spread of scattering angles causes muons from slightly different directions
   to mix. This effect is small for most applications.

5. Sea-level spectrum reference
   The surface spectrum (Guan, Frosin, etc.) is evaluated at sea level.
   For detectors at high altitude (>2000 m), the atmospheric depth is reduced
   and the actual muon flux is higher. An altitude correction to the surface
   spectrum would improve accuracy at alpine sites.

   Rough correction factor: multiply flux by exp((alt_m - 0) / 8500)
   At 2000m: factor ≈ 1.27   At 4000m: factor ≈ 1.60


================================================================================
FILES REFERENCE
================================================================================

gui/cosmoaleph_terrain_driver.py    Physics driver (subprocess target)
gui/gui_terrain_engine.py           Streamlit GUI panel
gui/cosmoaleph_backward_mc.py       Required: backward CSDA physics
PATCH_engine6_terrain.txt           4-edit patch for cosmoaleph_gui.py

DEM sources:
  https://portal.opentopography.org   (SRTM, free, no account for SRTM)
  https://earthexplorer.usgs.gov/     (SRTM+, free, account needed)
  https://spacedata.copernicus.eu/    (COP30, 30m, Europe, free)
