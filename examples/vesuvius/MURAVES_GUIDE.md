# Mt. Vesuvius — UCMuon/MURAVES Example

Standalone example reproducing the MURAVES muography simulation pipeline
(Rajan et al., Muographers 2026, Budapest) using UCMuon's backward-MC physics
and terrain engine. No GUI needed. Runs in **~3 seconds** at full resolution
and outputs four publication-style figures directly into `examples/vesuvius/figs/`.

**References**
- Hong et al. (2025) *J. Appl. Phys.* **138**, [doi:10.1063/5.0275078](https://doi.org/10.1063/5.0275078) — MURAVES project updates at Mt. Vesuvius
- Rajan et al. (2026) MURAVES/Muographers26 presentation

---

## Files in this folder

```
examples/vesuvius/
├── ucmuon_vesuvius_muraves.py   ← main simulation script
├── make_tsim_library.py         ← writes tsim_library/ (Part 2)
├── inversion_tests.py           ← the inversion beyond the closure test (Part 2)
├── parma_altitude.cc            ← PARMA at the site, for inversion_tests.py
├── MURAVES_GUIDE.md             ← this file
└── figs/                        ← output figures (auto-created on first run)
    ├── fig_vesuvius_thickness.png
    ├── fig_vesuvius_flux_elevation.png
    ├── fig_vesuvius_flux_maps.png
    └── fig_vesuvius_transmission.png
```

---

## What this example does

`ucmuon_vesuvius_muraves.py` runs the full UCMuon physics pipeline for
Mt. Vesuvius in a single script:

1. **Builds an overburden map** — ray-traces each (azimuth, elevation) direction
   from the MURAVES detector through a truncated-cone Vesuvius model (or a real
   DEM if supplied) to get rock thickness in g/cm².
2. **Computes the open-sky flux** — backward-MC integral over the Guan 2015
   spectrum, E = 1–2500 GeV, once per zenith angle.
3. **Computes the through-rock flux map** — backward-MC at each direction using
   the ray-traced path length.
4. **Computes the transmission map** — `T_sim = Φ_rock / Φ_sky` per direction.
5. **Generates four figures**: thickness map, flux vs elevation, flux maps, transmission ratio.

**Geometry** (MURAVES / Hong et al. 2025):

| | Value |
|---|---|
| Detector (SW flank, Osservatorio Vesuviano) | 40.8271°N, 14.4006°E, 608 m a.s.l. |
| Summit (Gran Cono) | 40.8218°N, 14.4265°E, 1281 m a.s.l. |
| Summit direction from detector | az = 105.1°, el = 16.6° above horizon |
| Cone model | truncated frustum, apex at summit, base radius 1900 m at 150 m a.s.l. |

**Spectrum**: Guan et al. 2015 (spectrum 3 in the backward-MC numbering), E_min = 1 GeV — required for the
low-elevation signal zone (el < 25°) where muons penetrate 100–1600 m of rock.

---

## Requirements

All standard packages — no extra install needed:

```bash
pip install numpy matplotlib scipy   # already in requirements.txt
```

Only needed for real DEM mode:
```bash
pip install rasterio
```

---

## Quick start

Run from the **UCMuon project root** (not from inside `examples/vesuvius/`):

```bash
cd /path/to/UCMuon
python examples/vesuvius/ucmuon_vesuvius_muraves.py
```

**Runtime: ~3 seconds** at full resolution (360 × 85 = 30 600 directions).  
Figures are saved to `examples/vesuvius/figs/` automatically.

---

## Command-line options

```
positional argument:
  [dem_file]        Path to DEM file (.tif / .xyz / .asc)
                    Omit to use the built-in synthetic cone model (default)

optional arguments:
  --rho   FLOAT     Bulk rock density [g/cm³]   default: 2.65
  --n-az  INT       Number of azimuth bins       default: 360  (1° per bin)
  --n-ze  INT       Number of zenith bins        default: 85   (1° per bin)
  --step  FLOAT     Ray-trace step size [m]      default: 25
  --out-dir PATH    Directory for output figures default: examples/vesuvius/figs
```

Examples:

```bash
# Full resolution with real DEM
python examples/vesuvius/ucmuon_vesuvius_muraves.py vesuvius_dem.tif

# Different density
python examples/vesuvius/ucmuon_vesuvius_muraves.py --rho 1.7

# Custom output directory
python examples/vesuvius/ucmuon_vesuvius_muraves.py --out-dir /tmp/myfigs

# Quick test at coarse resolution (~0.1 s)
python examples/vesuvius/ucmuon_vesuvius_muraves.py --n-az 72 --n-ze 30
```

---

## Output figures

### Fig 1 — Rock thickness map (`fig_vesuvius_thickness.png`)

![thickness](figs/fig_vesuvius_thickness.png)

**Left panel:** 2D map of rock thickness [m] as a function of azimuth and elevation.
The Vesuvius cone appears as a compact triangular structure centred at
az ≈ 105°, el ≈ 10–18° (white star = summit). Maximum thickness ≈ 1600 m at
the lowest elevations looking into the base of the cone. All other directions
see open sky (black = 0 m) — the Vesuvius cone subtends only ~40° × 12° of solid
angle from the detector.

**Right panel:** Slice through the summit azimuth (az = 105°). Rock thickness
drops sharply from ~1600 m at el = 5° to zero at el ≈ 16° — the upper edge of
the cone. Above that elevation the ray exits the cone entirely.

---

### Fig 2 — Flux vs elevation (`fig_vesuvius_flux_elevation.png`)

![flux_elevation](figs/fig_vesuvius_flux_elevation.png)

Integrated muon flux vs elevation angle at the summit azimuth (105°) for four
rock densities plus open sky.

- **Green dashed** (open sky): the Guan 2015 integral above 1 GeV, rising from
  1.1 m⁻² s⁻¹ sr⁻¹ at el = 5.5° to 5.6 at the summit elevation and 27 at
  el = 40°, roughly as cos²θ. At the zenith it is 58 (Guan's exact integral is
  60, the 50-point energy grid is 4 % low and cancels in T_sim); PDG quotes
  ≈ 70, with recent measurements 10–15 % lower.
- **Coloured curves** (ρ = 1.0, 2.0, 2.65, 3.0 g/cm³): one elevation bin below
  the summit (el = 15.5°) the flux is already 9 to 50 times below the open sky,
  at 0.55 (ρ = 1.0) and 0.10 (ρ = 2.65) m⁻² s⁻¹ sr⁻¹.
- **Density separation**: the curves stay apart all the way down the cone. At
  el = 10.5° they are 2.9 × 10⁻², 4.2 × 10⁻³, 1.6 × 10⁻³ and 9.7 × 10⁻⁴ for
  ρ = 1.0, 2.0, 2.65, 3.0: this is the density-sensitive zone used for
  inversion.
- **Near the cone base** the flux falls to 2.1 × 10⁻⁴ of the open sky at
  ρ = 2.65 and 1.1 × 10⁻⁴ at ρ = 3.0, under up to 4170 and 4730 m w.e. of rock.
  The backward MC's range table ends at 10 TeV (6889 m w.e.), and the flux
  needing more than that at the surface is left out: 2 % and 5 % here (the
  script prints a warning). Up to 1.3.2 the table ended at 2 TeV
  (3616 m w.e.): these curves dropped to zero below el ≈ 7° and 9°, and were
  6 % low at 2000 m w.e. and 37 % at 3000 m w.e.
- **Above el ≈ 17°**: all curves rejoin the open-sky reference — the ray has
  exited the cone completely.

The sharp upper edge is characteristic of the compact cone geometry. A real
DEM would show a broader, more gradual transition.

Up to v1.1.2 this figure showed the open sky flat at ~10⁻³ and the rock curves
dropping to ~10⁻⁸: the script summed the spectrum over a log-energy grid
without the factor E, and computed the rock flux as a rate over a whole cone of
directions at half the zenith angle. Both now use the backward MC's
`directional_flux`, as the terrain engine does.

---

### Fig 3 — Free-sky vs through-rock flux maps (`fig_vesuvius_flux_maps.png`)

![flux_maps](figs/fig_vesuvius_flux_maps.png)

**Left** (free-sky flux): smooth, azimuth-independent gradient — flux decreases
only with increasing zenith angle (lower elevation = higher zenith = lower flux),
from 58 m⁻² s⁻¹ sr⁻¹ at the zenith to about 1 near the horizon. On this linear
scale the whole band below el ≈ 20° is pale, so the cone shadow on the right is
faint; Fig 4 shows it clearly.

**Right** (through-rock flux, ρ = 2.65 g/cm³): identical to the left panel
everywhere except at the cone location (az ≈ 85–126°, el ≈ 5–16°) where the
flux collapses to near zero (white/yellow region = T_sim ≈ 0). The cone shadow
is the small bright triangular feature; bright here means near-zero flux because
the YlGn colormap runs white (low) → green (high).

> **Azimuth convention:** this script uses geographic convention (0 = North,
> 90 = East). MURAVES figures use a detector-centric convention where the
> summit is placed at az = 180°. Summit is at **az = 105°** here.
> To match MURAVES figures directly, rotate by 75° (180° − 105°).

---

### Fig 4 — Transmission ratio (`fig_vesuvius_transmission.png`)

![transmission](figs/fig_vesuvius_transmission.png)

`T_sim = Φ_rock / Φ_sky` per (az, el) pixel. Yellow = T_sim ≈ 1 (open sky,
~99% of the sky). Dark blue/purple = T_sim ≈ 0 (cone shadow). The cone
footprint is clearly triangular at az ≈ 85–126°, el ≈ 5–16°, with the apex
at the summit (white star).

Inside the cone at ρ = 2.65 g/cm³ the median transmission is 9 × 10⁻⁴: 51 % of
the blocked directions are below 0.001, and the rim reaches 0.16.
This map is the forward model input for density inversion: running at multiple
densities and comparing against measured T_data reveals the internal density
structure of the volcano.

---

## Part 2 — Density analysis (T_sim library + inversion)

To invert for density you need one `terrain_transmission.dat` file per bulk
density. A library for the MURAVES detector ships in `tsim_library/`
(ρ = 1.5, 2.0, 2.5, 2.65, 3.0 g/cm³; Guan 2015, 1–2500 GeV, real DEM), written
by `make_tsim_library.py` from the bundled `vesuvius_dem.tif` with the
integrator's forward calibration (mode 2, the default). Up to 1.4.1 the library
used mode 1, whose transmission in the pixels the inversion uses is 4.5 % below
the forward engines' (median; up to 13 %), which biased the synthetic test below
by -0.026 g/cm³; with mode 2 it is +0.010 g/cm³
(`benchmark/analysis/bmc_density_bias.py`). Directions beyond 5000 m w.e.,
the last opacity of the calibration, have no transmission (`nan`): 3 at
ρ = 1.5 and 350 at ρ = 3.0. Beyond it the correction is only held, and the
backward MC's range table (10 TeV) cuts more and more of the flux; up to 1.4.1
those directions were kept, as lower bounds, and 205 of them read exactly 0 at
ρ = 2.65. Transmissions are written with seven significant digits (six
decimals up to 1.4.1, which rounded T ≈ 10⁻⁶ by up to 40 %).

To write your own with the terrain driver CLI (stdin, one value per line: DEM,
latitude, longitude, altitude, density, spectrum, n_az, n_ze, ze_max, step,
integrator mode (2 = forward-calibrated), then the overburden, flux, summary
and transmission paths):

```bash
mkdir -p output/vesuvius_tsim

for rho in 1.5 2.0 2.65 3.0; do
    python gui/ucmuon_terrain_driver.py <<EOF
examples/vesuvius/vesuvius_dem.tif
40.8271
14.4006
608.0
${rho}
3
360
85
85.0
25.0
2
output/vesuvius_tsim/overburden_${rho}.dat
output/vesuvius_tsim/flux_${rho}.dat
output/vesuvius_tsim/summary_${rho}.dat
output/vesuvius_tsim/tsim_${rho}.dat
EOF
    echo "Done: rho=${rho}"
done
```

Spectrum 3 is Guan 2015. The CLI integrates the spectrum over 1–5000 GeV (1 GeV
is the lower end of Guan's fit); `make_tsim_library.py` uses 1–2500 GeV, which
gives the same T_sim: the median ratio is within 1 %, pixels with
T > 10⁻³ within 3 %. Each run ray-traces the DEM, which
takes about half a minute at this resolution.

### Inversion — GUI

1. `streamlit run gui/ucmuon_gui.py`
2. Open **🔬 Density** tab
3. Paste the four file paths into the T_sim Library box:
   ```
   examples/vesuvius/tsim_library/terrain_transmission_1.50.dat
   examples/vesuvius/tsim_library/terrain_transmission_2.00.dat
   examples/vesuvius/tsim_library/terrain_transmission_2.65.dat
   examples/vesuvius/tsim_library/terrain_transmission_3.00.dat
   ```
4. Click **Load Library**
5. Section 2 → **Generate synthetic**: True density = `2.0`, N events = `100000`
6. Section 3 → reference density = `2.65` → **▶ Run Inversion**

### Inversion — Python

```python
import sys; sys.path.insert(0, "gui")
from ucmuon_density_analysis import (
    load_transmission_map, build_tsim_library,
    invert_density_map, compute_double_ratio,
    generate_synthetic_tdata, write_density_map,
)

L = "examples/vesuvius/tsim_library/terrain_transmission"
tsim_lib = build_tsim_library([f"{L}_{r}.dat" for r in ("1.50", "2.00", "2.65", "3.00")])

# Synthetic test: true density = 2.0 g/cm³, 100k events
T_data, sigma_T = generate_synthetic_tdata(tsim_lib, true_rho=2.0, n_events=100_000)

rho_map, sigma_rho, status = invert_density_map(T_data, tsim_lib, sigma_T)
D_map = compute_double_ratio(T_data, tsim_lib[2.65])

az_c, el_c, _, meta = load_transmission_map(f"{L}_2.65.dat")
write_density_map(az_c, el_c, rho_map, sigma_rho, status,
                  "output/vesuvius_density.dat", meta)
```

### Expected inversion results

| Region | Status | Why |
|---|---|---|
| Most directions (el > 20°, no cone) | 1 (open sky) | T_sim ≈ 1 for all densities |
| Under the thickest rock (350 pixels, el 5.5–11.5°, median 7.5°) | 1 (no T_sim) | beyond 5000 m w.e. at ρ = 3.0, no library value |
| 12 pixels, el 8.5–12.5° | 4 (low sensitivity) | σ_ρ above 0.1 g/cm³ with 10⁵ events |
| The rest of the cone (697 pixels, el 5.5–16.5°, median 10.5°) | 0 (OK) | T_sim varies with ρ |

For the synthetic test with the shipped library (ρ_true = 2.0, N = 10⁵ open-sky
muons per pixel; the synthetic data are seeded, so these numbers reproduce
exactly):
- Recovered density in the 697 OK pixels: **ρ̂ = 2.004 g/cm³** (mean),
  pixel-to-pixel scatter 0.038 g/cm³, mean σ_ρ 0.033 g/cm³
- Double ratio D = 2.01 (median; D > 1 since ρ_true < ρ_ref = 2.65)

Status 4 comes from σ_ρ when the data have uncertainties (`max_sigma_rho`,
default 0.1 g/cm³; the GUI's "Max σ_ρ"), and from |dT/dρ| only when they do
not. Up to 1.4.1 it was always |dT/dρ| < 0.005, which ignores counting
statistics and flagged 591 pixels; with the binomial σ_T of the synthetic
data and a finite-difference slope, the pulls (ρ̂ − ρ_true)/σ_ρ had an RMS of
1.3, and 1.0 now.

### Beyond the closure test

The test above is a closure test: the data come from the library itself, at
one of its densities. `inversion_tests.py` (15 s; part 7 needs PARMA) measures
what it cannot show:

| Test | ρ_true [g/cm³] | Result |
|---|---|---|
| pull RMS, data from Engine 5, 20 seeds | 1.75, 2.0, 2.3, 2.8 | 0.99–1.09 (10⁴ events: 1.04–1.12) |
| interpolation, noise-free Engine 5 data | 1.75, 2.25, 2.8 | ρ̂ +0.5, +0.4, +0.1 % (interpolating T instead of ln T: +2.9, +3.0, +0.6 %) |
| data from the forward engines | 2.0, 2.3 | +0.009, +0.020 g/cm³ (Engine 5 data, same seeds: −0.001, +0.009) |
| data from Frosin 2025, Guan library | 2.0, 2.3 | −0.14, −0.16 g/cm³ |
| Method 2, flux models of `fast_flux_estimator.py` | 2.0 | ρ̂ 1.70 (Reyna), 2.00 (Guan, the default), 2.03 (Tang), 2.16 (Frosin), 0.15 (Gaisser, valid only above ~100 GeV: not offered in the GUI, and `invert_opacity_map` warns); the same at 0 and 608 m |
| sea-level library, data at 608 m (PARMA at the site) | 2.0 | read as 2.06–2.10 |
| exposure for 10⁵ open-sky muons per 1° × 1° pixel | | 1.9–9.4 m² yr of an ideal detector (median 4.1); with 10⁴ events 408 pixels invert, scatter 0.071 g/cm³ |

The flux model and the altitude dominate: the library uses a sea-level
spectrum (the terrain driver offers no PARMA), and the analytic models'
altitude factor cancels in T. Angular resolution, efficiency and the soft and
upward-going backgrounds, which matter most next to the skyline where these
pixels lie, are not modelled.

Up to v1.1.2 the shipped library gave no OK pixels in this test: its open-sky
reference integrated CosmoALEPH, a fit above ~100 GeV/c, from 0.5 GeV, which
put every T_sim 10²–10⁴ too low (`docs/FLUX_NORMALISATION_AUDIT.md`).

---

## Getting a real DEM

**A DEM is already bundled**: `vesuvius_dem.tif` in this directory (SRTM GL1
30 m, W=14.35, S=40.76, E=14.52, N=40.90; public domain, see
[`DEM_SOURCE.md`](DEM_SOURCE.md)). The GUI Terrain tab loads it by default,
and the script accepts it directly:

```bash
python examples/vesuvius/ucmuon_vesuvius_muraves.py examples/vesuvius/vesuvius_dem.tif
```

To re-download it or fetch a different area:

**Option A — eio (command line, no account):**
```bash
pip install elevation
eio clip -o vesuvius_dem.tif --bounds 14.35 40.76 14.52 40.90
```

**Option B — OpenTopography (browser, no account):**
1. Go to https://portal.opentopography.org
2. Select **SRTM GL1 (30 m)** or **Copernicus COP30 (30 m)**
3. Bounding box: W=14.35, S=40.76, E=14.52, N=40.90
4. Download GeoTIFF → save as `vesuvius_dem.tif`

**Option C — GUI Auto-download:**
Terrain tab (Advanced mode) → Section 2 "Geometry source" → auto-download → set the bounds above and an OpenTopography API key (or the demo key) → Download DEM.

Verify after downloading:
```bash
python - <<'EOF'
import rasterio, numpy as np
with rasterio.open("vesuvius_dem.tif") as src:
    d = src.read(1).astype(float)
    d[d == src.nodata] = np.nan
    print(f"Bounds: {src.bounds}")
    print(f"Elevation: {np.nanmin(d):.0f} – {np.nanmax(d):.0f} m")
    # Expected: ~50 – 1281 m
EOF
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: ucmuon_backward_mc` | Not running from project root | `cd /path/to/UCMuon` first |
| `ModuleNotFoundError: rasterio` | DEM mode without rasterio | `pip install rasterio` |
| Figures not found | Wrong working directory | Run from project root; figs go to `examples/vesuvius/figs/` |
| All T_sim ≈ 0 everywhere | Ray step too coarse for thin outer flanks | Use `--step 10` |
| Many pixels have status 4 (low sensitivity) | σ_ρ above "Max σ_ρ": too few events for the transmission there | Raise N events or Max σ_ρ; with 10⁵ events 697 pixels at el 5.5–16.5° invert (see Expected inversion results) |
| Density curves identical in Fig 2 | Cone edge narrower than bin width | Use default 360×85 grid (1° bins) |
| Summit at wrong azimuth vs MURAVES | Different convention | Add 75° to azimuth: MURAVES puts summit at 180°, here it is at 105° |
