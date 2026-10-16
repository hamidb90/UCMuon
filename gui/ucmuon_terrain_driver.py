#!/usr/bin/env python3
"""
cosmoaleph_terrain_driver.py  —  UCLouvain Muography Group
Engine 6: UCMuon Terrain  —  DEM-aware muon flux integrator.

What this does
──────────────
For a detector placed at a known GPS position, this driver:
  1. Reads a GeoTIFF Digital Elevation Model (SRTM 30m or 90m).
  2. For every (azimuth, zenith) direction bin, traces a ray backward
     from the detector through the terrain to find where it exits into
     open sky.  That intersection defines the muon entry point into rock.
  3. Computes the rock overburden (slant path × density) per direction.
  4. Applies the backward-CSDA flux formula (from cosmoaleph_backward_mc.py)
     to get the expected muon flux at the detector for each direction.
  5. Outputs:
       - A plain-text overburden map  (az × ze grid, g/cm²)
       - A plain-text flux map        (az × ze grid, m⁻² s⁻¹ sr⁻¹)
       - A summary statistics file

Why this is different from all other engines
─────────────────────────────────────────────
All other engines (UCMuon-MC, MUSIC, BB+MS, PROPOSAL, Backward MC)
assume a **flat uniform slab**: one depth, one density, all directions see
the same overburden.  This is a good approximation in the lab or for
borehole-mounted detectors, but fundamentally wrong for:
  - Alpine detectors (Puy de Dôme, Mont Blanc, alpine glaciers)
  - Volcano flank detectors (Etna, Merapi, Sakurajima)
  - CCS monitoring sites where surrounding hills contribute asymmetric screening
  - Any field deployment where the terrain height varies with azimuth

This engine replaces the fixed `depth_m` parameter with a full 2-D
overburden map derived from real topography.

Ray-tracing algorithm
──────────────────────
Uses a flat-Earth approximation valid for distances < ~100 km.  For each
(az, ze) direction:
  • Parametrise the ray in local East-North-Up (ENU) coordinates.
  • Sample the ray from the detector every min(step, DEM pixel / 4) metres
    (step configurable, default 50 m), converting the ENU offset to (lat, lon)
    and looking up the DEM elevation of the pixel.
  • Where the ray passes between rock (altitude below the DEM) and air from
    one sample to the next, locate the boundary by bisection.
  • Overburden = total slant path through rock × rho  [g/cm²] (all segments).
  • Directions that never hit terrain → open sky (overburden = 0).

Dependencies
─────────────
  rasterio ≥ 1.3     pip install rasterio
  numpy              (already required)
  scipy              (already required)
  requests           pip install requests  [optional, for auto-download only]

  cosmoaleph_backward_mc.py must be in the same directory (gui/).

Stdin protocol (one value per line):
   1  dem_file              path to GeoTIFF DEM  (.tif / .tiff)
   2  det_lat               detector latitude  [decimal degrees, WGS84]
   3  det_lon               detector longitude [decimal degrees, WGS84]
   4  det_alt_m             detector altitude above sea level [m]
   5  rho                   rock/material density [g/cm³]
   6  spectrum_mode         1=CosmoALEPH  2=Power-law  3=Guan  4=Frosin
   7  n_az                  number of azimuth bins  (default 36 = 10° steps)
   8  n_ze                  number of zenith bins   (default 18, i.e. ze_max/18 steps)
   9  ze_max_deg            maximum zenith angle    (default 85°)
  10  step_m                largest ray sampling step [m] (default 50; refined
                            to a quarter DEM pixel, boundaries bisected)
  11  mode                  0=CSDA only  1=+stochastic P_surv  (default 1)
  12  outfile_overburden    output overburden map        (default terrain_overburden.dat)
  13  outfile_flux          output flux map              (default terrain_flux.dat)
  14  outfile_summary       output summary stats         (default terrain_summary.dat)
  15  outfile_transmission  output transmission map      (default terrain_transmission.dat)
                            T_sim = Phi_rock / Phi_sky in elevation-angle convention

Author: Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026
"""

import sys
import importlib.util
import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Lazy-import rasterio (not in requirements until now; checked at runtime)
# ─────────────────────────────────────────────────────────────────────────────

def _import_rasterio():
    # Raises ImportError rather than sys.exit(): the GUI calls this module as a
    # library, and a SystemExit (not an Exception) ended the whole GUI script
    # run, blanking every tab after Terrain on an install without rasterio
    # (up to v1.2.0). The command-line entry point, main(), exits instead.
    try:
        import rasterio
        import rasterio.transform
        return rasterio
    except ImportError as exc:
        raise ImportError("rasterio is not installed; install it with: "
                          "pip install rasterio") from exc


# ─────────────────────────────────────────────────────────────────────────────
# DEM reader
# ─────────────────────────────────────────────────────────────────────────────

def _xyz_to_rasterio_transform(xs, ys, ncols, nrows):
    """Build a rasterio-compatible AffineTransform from XYZ grid extents."""
    import rasterio.transform as _rt
    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()
    res_x = (x_max - x_min) / (ncols - 1) if ncols > 1 else 1.0
    res_y = (y_max - y_min) / (nrows - 1) if nrows > 1 else 1.0
    # AffineTransform: (west, north) origin, pixel size, north-up (res_y negative)
    return _rt.from_bounds(
        x_min - res_x / 2, y_min - res_y / 2,
        x_max + res_x / 2, y_max + res_y / 2,
        ncols, nrows
    )


def load_dem_xyz(xyz_path):
    """
    Load an XYZ point-cloud DEM.

    Format: three whitespace-separated columns per line:
        X[m or °E]   Y[m or °N]   Z[m a.s.l.]
    Comment lines starting with # are skipped.
    No header is required.

    Coordinates must be geographic (lon/lat, WGS84): X in ±180, Y in ±90.
    An XYZ file carries no coordinate reference system, so projected
    coordinates (e.g. UTM metres) cannot be placed on the map and raise a
    ValueError; convert them to lon/lat, or use a GeoTIFF that records its
    CRS. (Up to 1.3.2 they were converted with a flat-earth formula that
    ignored the UTM zone and false easting: Vesuvius in UTM 33N landed at
    5.3° E instead of 14.4° E, every direction then saw open sky.)

    Returns:
        elev      : 2-D float32 array (nrows × ncols), rows sorted N→S
        transform : rasterio-compatible AffineTransform (lon/lat)
    """
    data = np.loadtxt(xyz_path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            f"XYZ file must have ≥ 3 columns (X Y Z); got shape {data.shape}. "
            "Check file format."
        )

    xs = data[:, 0].astype(float)
    ys = data[:, 1].astype(float)
    zs = data[:, 2].astype(float)

    # Detect geographic vs projected
    _geo = (xs.min() >= -180.01 and xs.max() <= 180.01 and
            ys.min() >=  -90.01 and ys.max() <=  90.01)

    if not _geo:
        raise ValueError(_projected_message(xyz_path, "an XYZ file"))
    lons = xs
    lats = ys

    # Build regular grid
    lon_unique = np.unique(np.round(lons, 8))
    lat_unique = np.unique(np.round(lats, 8))
    ncols = len(lon_unique)
    nrows = len(lat_unique)

    elev = np.full((nrows, ncols), np.nan, dtype=np.float32)
    lon_idx = {v: i for i, v in enumerate(lon_unique)}
    lat_idx = {v: i for i, v in enumerate(lat_unique)}

    for lon_v, lat_v, z_v in zip(
            np.round(lons, 8), np.round(lats, 8), zs):
        ci = lon_idx.get(lon_v)
        ri = lat_idx.get(lat_v)
        if ci is not None and ri is not None:
            elev[ri, ci] = z_v

    # rasterio expects rows ordered N→S (decreasing latitude)
    elev = elev[::-1, :]

    transform = _xyz_to_rasterio_transform(lon_unique, lat_unique, ncols, nrows)

    print(f"  XYZ DEM loaded: {Path(xyz_path).name}  "
          f"shape={elev.shape}  "
          f"min={np.nanmin(elev):.0f} m  max={np.nanmax(elev):.0f} m"
          f"  (geographic)", flush=True)
    return elev, transform


def _projected_message(path, what):
    return (f"{Path(path).name}: the coordinates look projected (metres, e.g. UTM), "
            f"and {what} has no coordinate reference system to place them on the "
            "map. Convert the DEM to longitude/latitude (WGS84), or to a GeoTIFF "
            "that records its CRS (for an .asc, a .prj file next to it is enough).")


def load_dem_asc(asc_path):
    """
    Load an Esri ASCII Grid (.asc) DEM.

    Header format (first 5–6 lines):
        ncols         <int>
        nrows         <int>
        xllcorner     <float>   or   xllcenter   <float>
        yllcorner     <float>   or   yllcenter   <float>
        cellsize      <float>
        NODATA_value  <float>   (optional)

    Values are elevation [m], stored row by row top-to-bottom (N→S).
    Geographic (lon/lat) grids are read directly. A projected grid needs its
    CRS from a .prj file next to it (it is then read through rasterio and
    reprojected, as a GeoTIFF); without one it raises a ValueError, as in
    load_dem_xyz.

    Returns:
        elev      : 2-D float32 array (nrows × ncols)
        transform : rasterio-compatible AffineTransform
    """
    import rasterio.transform as _rt

    with open(asc_path) as fh:
        header = {}
        first_data = []
        while len(header) < 6:
            ln = fh.readline()
            if not ln:
                break
            parts = ln.strip().split()
            if len(parts) == 2 and parts[0][:1].isalpha():
                header[parts[0].lower()] = parts[1]
            else:
                # hit the data rows: keep this one (up to 1.3.2 a header
                # without the optional NODATA_value line lost its first row)
                first_data = [ln]
                break
        data_lines = first_data + fh.readlines()

    ncols    = int(header["ncols"])
    nrows    = int(header["nrows"])
    cellsize = float(header["cellsize"])
    # xllcorner = lower-left corner; xllcenter = lower-left pixel centre, half
    # a cell further in (up to 1.3.2 the centre was used as the corner, which
    # shifted the grid half a cell north-east).
    xll      = float(header.get("xllcenter", header.get("xllcorner", 0)))
    yll      = float(header.get("yllcenter", header.get("yllcorner", 0)))
    if "xllcenter" in header:
        xll -= cellsize / 2.0
    if "yllcenter" in header:
        yll -= cellsize / 2.0
    nodata   = float(header.get("nodata_value", -9999))

    rows = []
    for ln in data_lines:
        ln = ln.strip()
        if ln:
            rows.append([float(v) for v in ln.split()])
    elev = np.array(rows, dtype=np.float32)
    elev[elev == nodata] = np.nan

    # Detect geographic vs projected
    x_max = xll + ncols * cellsize
    y_max = yll + nrows * cellsize
    _geo  = (xll >= -181 and x_max <= 181 and yll >= -91 and y_max <= 91)

    if not _geo:
        if Path(asc_path).with_suffix(".prj").exists():
            return _load_dem_geotiff(asc_path)      # GDAL reads the grid + .prj
        raise ValueError(_projected_message(asc_path, "an .asc file without a .prj"))
    cellsize_lon = cellsize_lat = cellsize

    # ASC stores top row first (N→S); transform origin is top-left
    transform = _rt.from_origin(
        xll, yll + nrows * cellsize_lat,   # west, north
        cellsize_lon, cellsize_lat          # pixel width, pixel height (both positive)
    )

    print(f"  ASC DEM loaded: {Path(asc_path).name}  "
          f"shape={elev.shape}  "
          f"min={np.nanmin(elev):.0f} m  max={np.nanmax(elev):.0f} m"
          f"  (geographic)", flush=True)
    return elev, transform


def load_dem(dem_path):
    """
    Load a DEM from any supported format and return (elevation_array, transform).

    Supported formats (auto-detected from file extension):
      .tif / .tiff  →  GeoTIFF via rasterio   (SRTM, Copernicus, etc.)
      .xyz          →  XYZ point cloud         (e.g. MURAVES 5-m Vesuvius LIDAR)
      .asc          →  Esri ASCII Grid         (common European survey format)

    Returns:
        elev      : 2-D float32 array [m a.s.l.], rows N→S
        transform : rasterio AffineTransform (col/row ↔ lon/lat)
    """
    ext = Path(dem_path).suffix.lower()
    if ext in (".xyz", ".txt"):   # .txt accepted too (users may keep original extension)
        return load_dem_xyz(dem_path)
    elif ext in (".asc",):
        return load_dem_asc(dem_path)
    elif ext in (".tif", ".tiff", ".geotiff"):
        return _load_dem_geotiff(dem_path)
    else:
        # Try GeoTIFF anyway (rasterio can handle many formats)
        try:
            return _load_dem_geotiff(dem_path)
        except Exception:
            raise ValueError(
                f"Unsupported DEM format: '{ext}'.  "
                "Supported: .tif/.tiff (GeoTIFF), .xyz/.txt (XYZ point cloud), .asc (Esri ASCII)."
            )


def _load_dem_geotiff(dem_path):
    """Internal: load a raster via rasterio, as lon/lat.

    The ray tracer works in longitude/latitude, so a raster in a projected
    CRS (UTM, Lambert, ...) is reprojected to WGS84 (bilinear). Up to 1.3.2 its
    metre transform was taken as degrees: the detector fell outside the DEM
    and every direction saw open sky. A raster with no CRS must be in lon/lat.
    """
    rio = _import_rasterio()
    note = ""
    with rio.open(dem_path) as ds:
        nodata = ds.nodata
        elev = ds.read(1).astype(np.float32)
        if nodata is not None:
            elev[elev == nodata] = np.nan
        transform = ds.transform
        if ds.crs is not None and not ds.crs.is_geographic:
            from rasterio.warp import calculate_default_transform, reproject, Resampling
            dst_t, w, h = calculate_default_transform(ds.crs, "EPSG:4326", ds.width,
                                                      ds.height, *ds.bounds)
            dst = np.full((h, w), np.nan, dtype=np.float32)
            reproject(elev, dst, src_transform=ds.transform, src_crs=ds.crs,
                      dst_transform=dst_t, dst_crs="EPSG:4326",
                      resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
            note = f"  (reprojected from {ds.crs.to_string()} to WGS84 lon/lat)"
            elev, transform = dst, dst_t
        elif ds.crs is None:
            b = ds.bounds
            if not (b.left >= -181 and b.right <= 181 and b.bottom >= -91 and b.top <= 91):
                raise ValueError(_projected_message(dem_path, "this raster"))
    print(f"  GeoTIFF DEM loaded: {Path(dem_path).name}  "
          f"shape={elev.shape}  "
          f"min={np.nanmin(elev):.0f} m  max={np.nanmax(elev):.0f} m{note}", flush=True)
    return elev, transform


def _dem_elevations(elev, transform, lat, lon):
    """DEM elevation [m] at arrays of points: bilinear interpolation between
    pixel centres, i.e. a continuous surface (NaN outside the grid; next to a
    NaN pixel, the weighted mean of the valid neighbours).

    Up to 1.3.2 the lookup was meant to be bilinear, but its fractional
    position came from rasterio's rowcol, which returns integers, so every
    weight was zero and the surface was a staircase of 30 m pixels. With
    boundaries located exactly (ray_overburden), the steps put metres of
    spurious rock in front of a detector on a slope in every direction.
    """
    colf, rowf = (~transform) * (np.asarray(lon, float), np.asarray(lat, float))
    u, v = colf - 0.5, rowf - 0.5                      # pixel-centre coordinates
    c0 = np.floor(u).astype(np.int64)
    r0 = np.floor(v).astype(np.int64)
    nr, nc = elev.shape
    ok = (r0 >= 0) & (r0 < nr - 1) & (c0 >= 0) & (c0 < nc - 1)
    out = np.full(np.shape(c0), np.nan)
    r, c = r0[ok], c0[ok]
    fu, fv = (u - c0)[ok], (v - r0)[ok]
    q = np.stack([elev[r, c], elev[r, c + 1], elev[r + 1, c], elev[r + 1, c + 1]]).astype(float)
    w = np.stack([(1 - fu) * (1 - fv), fu * (1 - fv), (1 - fu) * fv, fu * fv])
    good = ~np.isnan(q)
    num = np.sum(w * np.where(good, q, 0.0), axis=0)
    den = np.sum(w * good, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        val = np.where(good.all(axis=0), num, num / den)
    val[~good.any(axis=0) | (den <= 0)] = np.nan
    out[ok] = val
    return out


def dem_elevation_at(elev, transform, lat, lon):
    """
    DEM elevation [m] at (lat, lon): bilinear between pixel centres
    (_dem_elevations). Returns NaN if outside the DEM extent.
    """
    return float(_dem_elevations(elev, transform, np.array([lat], float),
                                 np.array([lon], float))[0])


# ─────────────────────────────────────────────────────────────────────────────
# Ray tracing
# ─────────────────────────────────────────────────────────────────────────────

# Flat-earth conversion constants
_M_PER_DEG_LAT = 111_320.0    # metres per degree of latitude (approximate)

def ray_overburden(elev, transform,
                   det_lat, det_lon, det_alt_m,
                   azimuth_deg, zenith_deg,
                   rho, step_m=50.0, max_dist_m=60_000.0,
                   underground=False):
    """
    Trace a ray from the detector outward and ACCUMULATE ALL underground path
    segments — i.e. integrate the full rock column along the slant path.

    Azimuth convention: 0° = North, 90° = East (geographic, clockwise).
    Zenith convention:  0° = straight up,  90° = horizontal.

    Returns (overburden_gcm2, total_slant_m, open_sky).
      overburden_gcm2  : total_slant_path_cm × rho  [g/cm²]
      total_slant_m    : total accumulated slant path through rock [m]
      open_sky         : True if zero rock encountered along the ray

    NOTE — why segment accumulation instead of first-intersection only:
    ──────────────────────────────────────────────────────────────────────
    For mountain muography the detector sits outside (or on the flank of) the
    topographic structure.  A ray traced backward from the detector can:
      1. Pass through a thin near-surface skin (GPS/DEM altitude mismatch)
      2. Travel through open air
      3. Enter the main mountain body from the far side
      4. Exit to open sky on the opposite flank
    The first-intersection approach captures only segment (1) and misses the
    dominant rock column (3).  It also fails catastrophically when the entered
    detector altitude is slightly below the DEM raster value (e.g. 12 m in the
    Vesuvius demo), producing entry_dist ≈ 0 → opacity ≈ 6 g/cm² for all
    blocked bins regardless of the true overburden.

    The accumulator approach sums every segment where the ray flies below
    the DEM surface, each bounded by bisection, giving the correct total
    column for any topology (flat slab, single peak, complex ridge, detector
    on flank, etc.).
    """
    az  = np.radians(azimuth_deg)
    ze  = np.radians(zenith_deg)

    # ENU direction unit vector (geographic convention: az from North, CW)
    sin_ze = np.sin(ze)
    dx_e   =  sin_ze * np.sin(az)   # East
    dx_n   =  sin_ze * np.cos(az)   # North
    dx_u   =  np.cos(ze)            # Up

    # Flat-earth lat/lon-per-metre factors
    lat_per_m = 1.0 / _M_PER_DEG_LAT
    lon_per_m = 1.0 / (_M_PER_DEG_LAT * max(np.cos(np.radians(det_lat)), 1e-6))

    # ------------------------------------------------------------------
    # Determine effective detector altitude.
    #
    # Surface deployment (underground=False, default):
    #   Clamp det_alt_eff upward to the DEM surface if the GPS altitude
    #   sits below the raster value. This absorbs GPS measurement error
    #   and 30-m pixel averaging artefacts without inflating near-field
    #   overburden for upward-going rays.
    #
    # Underground deployment (underground=True):
    #   Trust the supplied altitude exactly — DO NOT clamp upward.
    #   The detector is intentionally below the DEM surface (borehole,
    #   mine, cavern).  Clamping would teleport it to the surface and
    #   return zero overburden for every direction, which is wrong.
    #   The DEM still defines the rock boundary; the ray accumulates
    #   every step_m interval where ray_alt < DEM_alt, starting from
    #   the true underground position.
    # ------------------------------------------------------------------
    dem_at_det = dem_elevation_at(elev, transform, det_lat, det_lon)
    if np.isnan(dem_at_det):
        return 0.0, 0.0, True   # detector outside DEM — assume open sky
    if underground:
        det_alt_eff = float(det_alt_m)          # trust the given altitude
    else:
        det_alt_eff = max(float(det_alt_m),     # surface: absorb GPS/DEM mismatch
                          float(dem_at_det) + 0.5)

    # ------------------------------------------------------------------
    # Sample the ray, then locate every rock boundary by bisection.
    #
    # The ray is sampled every h = min(step_m, a quarter of the DEM pixel)
    # (64 times finer over the first four pixels) from the detector itself
    # (d = 0) to the edge of the DEM, and each change
    # between rock and air from one sample to the next is located by
    # bisection, so every segment is measured to ~1e-6 h whatever the step.
    # Up to 1.3.2 the walk counted whole steps of step_m where the sample was
    # below the DEM: a detector 100 m under a flat surface got 50 m of rock
    # vertically with the default 50 m step, and on a flank, rock within the
    # first step of the detector was missed (open sky).
    # ------------------------------------------------------------------
    pix_m = min(abs(transform.a) * _M_PER_DEG_LAT * max(np.cos(np.radians(det_lat)), 1e-6),
                abs(transform.e) * _M_PER_DEG_LAT)
    h = min(float(step_m), 0.25 * pix_m) if pix_m > 0 else float(step_m)
    h = max(h, 0.05)

    def _inside(dist):
        lat = det_lat + dx_n * dist * lat_per_m
        lon = det_lon + dx_e * dist * lon_per_m
        dem = _dem_elevations(elev, transform, lat, lon)
        with np.errstate(invalid="ignore"):
            return (det_alt_eff + dx_u * dist) < dem, np.isfinite(dem)

    # Finer sampling over the first four pixels, where every ray starts at
    # the ground and a pixel corner can put a few metres of rock in front
    # of the detector.
    near = min(4.0 * pix_m, max_dist_m) if pix_m > 0 else 0.0
    d = np.concatenate((np.arange(0.0, near, max(min(h, pix_m / 64.0), 0.05)),
                        np.arange(near, max_dist_m + h, h)))
    ins, valid = _inside(d)
    if not valid.all():                       # the ray leaves the DEM extent:
        n = int(np.argmin(valid))             # assume open sky beyond
        d, ins = d[:n], ins[:n]
    if len(d) == 0:
        return 0.0, 0.0, True

    # Boundaries between samples k-1 and k, refined by bisection
    k = np.nonzero(ins[1:] != ins[:-1])[0] + 1
    lo, hi = d[k - 1].copy(), d[k].copy()
    s_lo = ins[k - 1]
    for _ in range(25):
        mid = 0.5 * (lo + hi)
        m_in, _ = _inside(mid)
        same = m_in == s_lo
        lo = np.where(same, mid, lo)
        hi = np.where(same, hi, mid)
    cross = 0.5 * (lo + hi)

    # Rock segments: from an air->rock boundary (or the detector, if it is in
    # rock) to the next rock->air boundary (or the last sample in the DEM)
    starts = list(cross[~s_lo])
    ends = list(cross[s_lo])
    if ins[0]:
        starts.insert(0, 0.0)
    if ins[-1]:
        ends.append(float(d[-1]))
    total_underground_m = float(np.sum(np.asarray(ends) - np.asarray(starts))) if starts else 0.0

    if total_underground_m <= 0.0:
        return 0.0, 0.0, True

    # Overburden: accumulated slant path [cm] × density [g/cm³] → [g/cm²]
    opacity = total_underground_m * 100.0 * rho
    return opacity, total_underground_m, False


def compute_overburden_map(elev, transform,
                            det_lat, det_lon, det_alt_m,
                            rho, n_az=36, n_ze=18, ze_max_deg=85.0,
                            step_m=50.0, progress_cb=None,
                            underground=False):
    """
    Compute overburden [g/cm²] for every (azimuth, zenith) bin.

    underground : bool (default False)
        Set True when the detector is deployed below the terrain surface
        (borehole, mine drift, cavern).  Disables the GPS/DEM altitude
        clamping in ray_overburden() so the true underground position is
        used as the ray origin.  See ray_overburden() docstring for details.

    Returns:
        az_centres   : 1-D array of azimuth bin centres  [deg]
        ze_centres   : 1-D array of zenith bin centres   [deg]
        overburden   : 2-D array shape (n_az, n_ze)       [g/cm²]
        open_sky_map : 2-D bool array  (True = open sky / no rock above)
    """
    az_edges = np.linspace(0.0, 360.0, n_az + 1)
    ze_edges = np.linspace(0.0, ze_max_deg, n_ze + 1)
    az_c     = 0.5 * (az_edges[:-1] + az_edges[1:])
    ze_c     = 0.5 * (ze_edges[:-1] + ze_edges[1:])

    overburd  = np.zeros((n_az, n_ze), dtype=np.float64)
    open_sky  = np.zeros((n_az, n_ze), dtype=bool)

    total = n_az * n_ze
    done  = 0

    for ia, az in enumerate(az_c):
        for iz, ze in enumerate(ze_c):
            ob, _, sky = ray_overburden(
                elev, transform,
                det_lat, det_lon, det_alt_m,
                az, ze, rho, step_m,
                underground=underground,
            )
            overburd[ia, iz] = ob
            open_sky[ia, iz] = sky
            done += 1

        if progress_cb and (ia + 1) % max(1, n_az // 10) == 0:
            progress_cb(done, total)

    return az_c, ze_c, overburd, open_sky


# ─────────────────────────────────────────────────────────────────────────────
# Flux integration  (delegates to cosmoaleph_backward_mc.py physics)
# ─────────────────────────────────────────────────────────────────────────────

def _load_bmc(script_dir):
    """Import ucmuon_backward_mc from the gui/ directory."""
    bmc_path = Path(script_dir) / "ucmuon_backward_mc.py"
    if not bmc_path.exists():
        raise FileNotFoundError(f"ucmuon_backward_mc.py not found in {script_dir}")
    spec = importlib.util.spec_from_file_location("ucmuon_backward_mc",
                                                   str(bmc_path))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def compute_flux_map(az_c, ze_c, overburden, open_sky_map,
                     rho, spectrum_mode=1, mode=1,
                     E_min_GeV=1.0, E_max_GeV=5000.0,
                     n_E=40, script_dir=None, progress_cb=None):
    """
    For every (azimuth, zenith) bin compute the expected muon flux [m⁻² s⁻¹ sr⁻¹]
    using backward CSDA integration with the overburden from ray tracing.

    open_sky directions → assigned open-sky flux (max flux for that zenith angle).
    """
    bmc     = _load_bmc(script_dir or Path(__file__).parent)
    n_az, n_ze = overburden.shape
    # E_min defaults to 1 GeV, the lower end of Guan's and Frosin's fits (it was
    # 0.5 GeV up to v1.1.2); say so when a spectrum is used outside its range.
    note = bmc.spectrum_range_warning(spectrum_mode, E_min_GeV)
    if note:
        print(f"  WARNING: {note}", flush=True)
    flux_map = np.zeros((n_az, n_ze), dtype=np.float64)

    # Open-sky reference per zenith angle: same integrator with X = 0, so
    # transmission = flux_map / opensky_flux is consistent by construction.
    opensky_flux = np.zeros(n_ze)
    for iz, ze in enumerate(ze_c):
        opensky_flux[iz] = bmc.directional_flux(
            0.0, np.radians(ze), spectrum_mode,
            E_min_GeV=E_min_GeV, E_max_GeV=E_max_GeV, n_E=n_E, mode=mode,
        )

    total = n_az * n_ze
    done  = 0

    with bmc.collect_range_warnings() as beyond_table:
        for ia in range(n_az):
            for iz, ze in enumerate(ze_c):
                X = overburden[ia, iz]
                if open_sky_map[ia, iz] or X < 1.0:   # effectively open sky
                    flux_map[ia, iz] = opensky_flux[iz]
                else:
                    # Exact slant opacity at the exact zenith angle — no
                    # vertical-equivalent depth round trip, no cone average.
                    flux_map[ia, iz] = bmc.directional_flux(
                        X, np.radians(ze), spectrum_mode,
                        E_min_GeV=E_min_GeV, E_max_GeV=E_max_GeV,
                        n_E=n_E, mode=mode,
                    )
                done += 1

            if progress_cb and (ia + 1) % max(1, n_az // 10) == 0:
                progress_cb(done, total)

    note = bmc.range_summary(beyond_table, total)
    if note:
        print(f"  WARNING: {note}", flush=True)

    return flux_map, opensky_flux


def compute_transmission_map(flux_map, opensky_flux):
    """
    T_sim(az, ze) = Phi_rock / Phi_sky — transmission ratio in [0, 1].

    flux_map     : 2-D array (n_az, n_ze) from compute_flux_map — through-rock flux.
    opensky_flux : 1-D array (n_ze,)       from compute_flux_map — open-sky reference.

    Returns a 2-D array of shape (n_az, n_ze).  Open-sky bins where
    both numerator and denominator are equal give T ≈ 1.  Completely
    blocked directions give T ≈ 0.
    """
    opensky_2d = np.broadcast_to(opensky_flux[np.newaxis, :], flux_map.shape)
    with np.errstate(divide='ignore', invalid='ignore'):
        T = np.where(opensky_2d > 0.0, flux_map / opensky_2d, 0.0)
    return np.clip(T, 0.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Output writers
# ─────────────────────────────────────────────────────────────────────────────

def write_overburden_map(az_c, ze_c, overburden, open_sky_map, fpath,
                          det_lat, det_lon, det_alt_m, rho):
    with open(fpath, "w") as fh:
        fh.write("# UCMuon Terrain Engine — Overburden Map\n")
        fh.write(f"# Detector: lat={det_lat:.6f}  lon={det_lon:.6f}"
                 f"  alt={det_alt_m:.1f} m\n")
        fh.write(f"# Density: {rho:.3f} g/cm3\n")
        fh.write(f"# Grid: {len(az_c)} az bins  x  {len(ze_c)} ze bins\n")
        fh.write("# Cols: azimuth[deg]  zenith[deg]  overburden[g/cm2]  open_sky\n")
        for ia, az in enumerate(az_c):
            for iz, ze in enumerate(ze_c):
                fh.write(f"{az:8.2f}  {ze:7.2f}  {overburden[ia,iz]:14.2f}"
                         f"  {int(open_sky_map[ia,iz])}\n")


def write_flux_map(az_c, ze_c, flux_map, fpath, det_lat, det_lon, det_alt_m):
    with open(fpath, "w") as fh:
        fh.write("# UCMuon Terrain Engine — Directional Flux Map\n")
        fh.write(f"# Detector: lat={det_lat:.6f}  lon={det_lon:.6f}"
                 f"  alt={det_alt_m:.1f} m\n")
        fh.write("# Cols: azimuth[deg]  zenith[deg]  flux[m-2 s-1 sr-1]\n")
        for ia, az in enumerate(az_c):
            for iz, ze in enumerate(ze_c):
                fh.write(f"{az:8.2f}  {ze:7.2f}  {flux_map[ia,iz]:14.6e}\n")


def write_transmission_map(az_c, ze_c, T_sim, fpath,
                           det_lat, det_lon, det_alt_m, rho):
    """
    Write the transmission map T_sim in elevation-angle convention.

    Elevation = 90° − zenith so that 0° = horizon, 90° = vertical.
    This matches the MURAVES / muography community convention where
    the "interesting" rock directions sit at low elevations (5–30°),
    which would be large zenith angles (60–85°) in the zenith convention.

    Columns: azimuth[deg]  elevation[deg]  transmission[0-1]
    """
    el_c = 90.0 - ze_c
    with open(fpath, "w") as fh:
        fh.write("# UCMuon Terrain Engine — Transmission Map\n")
        fh.write("# T_sim(az, el) = Phi_rock / Phi_sky\n")
        fh.write(f"# Detector: lat={det_lat:.6f}  lon={det_lon:.6f}"
                 f"  alt={det_alt_m:.1f} m\n")
        fh.write(f"# Density: {rho:.3f} g/cm3\n")
        fh.write(f"# Grid: {len(az_c)} az bins  x  {len(ze_c)} el bins"
                 f"  (elevation = 90 - zenith)\n")
        fh.write("# Cols: azimuth[deg]  elevation[deg]  transmission\n")
        for ia, az in enumerate(az_c):
            for iz in range(len(ze_c)):
                fh.write(f"{az:8.2f}  {el_c[iz]:7.2f}  {T_sim[ia, iz]:.6f}\n")


def write_summary(az_c, ze_c, overburden, flux_map, open_sky_map,
                  fpath, det_lat, det_lon, det_alt_m, rho,
                  spectrum_mode, elapsed):
    # Rate through a horizontal 1 m² detector over the mapped sky:
    # Σ Φ(az, ze) cos(ze) ΔΩ. Up to v1.1.2 this was np.sum(flux_map), a sum of
    # m⁻² s⁻¹ sr⁻¹ values over bins with no solid angle, ~100x too high on the
    # default grid and dependent on the binning.
    d_az = 2.0 * np.pi / len(az_c)
    half = 0.5 * np.radians(ze_c[1] - ze_c[0]) if len(ze_c) > 1 else 0.0
    ze_r = np.radians(ze_c)
    proj_dom = d_az * 0.5 * (np.sin(ze_r + half) ** 2 - np.sin(ze_r - half) ** 2)
    total_rate = float(np.sum(flux_map * proj_dom[np.newaxis, :]))
    n_open = int(open_sky_map.sum())
    n_rock = open_sky_map.size - n_open
    ob_rock = overburden[~open_sky_map]
    ob_med  = float(np.median(ob_rock)) if ob_rock.size else 0.0
    ob_max  = float(np.max(ob_rock))    if ob_rock.size else 0.0

    # Direction of maximum overburden
    imax = np.unravel_index(np.argmax(overburden), overburden.shape)
    az_max_ob = az_c[imax[0]]
    ze_max_ob = ze_c[imax[1]]

    # Direction of maximum flux
    imax_f = np.unravel_index(np.argmax(flux_map), flux_map.shape)
    az_max_flux = az_c[imax_f[0]]
    ze_max_flux = ze_c[imax_f[1]]

    with open(fpath, "w") as fh:
        fh.write("# UCMuon Terrain Engine — Summary\n")
        fh.write(f"# Detector lat={det_lat:.6f} lon={det_lon:.6f}"
                 f" alt={det_alt_m:.1f} m  rho={rho:.3f} g/cm3\n")
        _names = {1: "CosmoALEPH", 2: "power law", 3: "Guan 2015", 4: "Frosin 2025"}
        fh.write(f"# Spectrum: {_names.get(int(spectrum_mode), spectrum_mode)} "
                 f"(backward-MC number {spectrum_mode})\n")
        fh.write(f"# Computed in {elapsed:.1f} s\n")
        fh.write(f"#\n")
        fh.write(f"# Open-sky directions : {n_open} / {open_sky_map.size}\n")
        fh.write(f"# Rock directions     : {n_rock} / {open_sky_map.size}\n")
        fh.write(f"# Median overburden   : {ob_med:.0f} g/cm2\n")
        fh.write(f"# Max overburden      : {ob_max:.0f} g/cm2"
                 f"  @ az={az_max_ob:.1f} deg  ze={ze_max_ob:.1f} deg\n")
        fh.write(f"# Rate, horizontal 1 m2, ze < {ze_c[-1] + np.degrees(half):.0f} deg :"
                 f" {total_rate:.4e} m-2 s-1\n")
        fh.write(f"# Peak flux direction : az={az_max_flux:.1f} deg"
                 f"  ze={ze_max_flux:.1f} deg\n")

    # Also print to stdout for GUI live panel
    print(f"  ═════════════════════════════════════", flush=True)
    print(f"  UCMuon Terrain — Summary", flush=True)
    print(f"  Open-sky directions : {n_open}/{open_sky_map.size}", flush=True)
    print(f"  Median overburden   : {ob_med:.0f} g/cm2", flush=True)
    print(f"  Max overburden      : {ob_max:.0f} g/cm2"
          f"  @ az={az_max_ob:.0f}° ze={ze_max_ob:.0f}°", flush=True)
    print(f"  Rate (horizontal m², ze < {ze_c[-1] + np.degrees(half):.0f}°) :"
          f" {total_rate:.4e} m-2 s-1", flush=True)
    print(f"  ═════════════════════════════════════", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    import time

    lines = [ln.strip() for ln in sys.stdin
             if ln.strip() and not ln.strip().startswith("#")]

    def _rd(i, default, typ=str):
        try:    return typ(lines[i]) if i < len(lines) else default
        except: return default

    dem_file     = _rd(0,  "")
    det_lat      = _rd(1,  50.668, float)   # default: UCLouvain
    det_lon      = _rd(2,   4.615, float)
    det_alt_m    = _rd(3,   90.0,  float)
    rho          = _rd(4,   2.65,  float)
    spec_mode    = _rd(5,   1,     int)
    n_az         = _rd(6,  36,     int)
    n_ze         = _rd(7,  18,     int)
    ze_max_deg   = _rd(8,  85.0,   float)
    step_m       = _rd(9,  50.0,   float)
    mode         = _rd(10,  1,     int)
    out_ob       = _rd(11, "terrain_overburden.dat")
    out_flux     = _rd(12, "terrain_flux.dat")
    out_summ     = _rd(13, "terrain_summary.dat")
    out_trans    = _rd(14, "terrain_transmission.dat")

    if not dem_file or not Path(dem_file).exists():
        print(f"  ERROR: DEM file not found: '{dem_file}'", flush=True)
        print(f"  Provide a GeoTIFF file (.tif) as the first stdin line.", flush=True)
        sys.exit(1)

    print(f"  UCMuon Terrain Engine v1.1", flush=True)
    print(f"  DEM:      {dem_file}", flush=True)
    print(f"  Detector: lat={det_lat:.6f}  lon={det_lon:.6f}  alt={det_alt_m:.1f} m", flush=True)
    print(f"  Grid:     {n_az} az × {n_ze} ze  (ze_max={ze_max_deg}°  step={step_m} m)", flush=True)
    print(f"  Material: rho={rho:.3f} g/cm3  spectrum={spec_mode}", flush=True)
    print(f"", flush=True)

    t0 = time.time()

    # 1. Load DEM
    try:
        elev, transform = load_dem(dem_file)
    except ImportError as exc:
        print(f"  ERROR: {exc}", flush=True)
        sys.exit(1)

    # 2. Overburden map
    print(f"  Ray tracing overburden map ({n_az} × {n_ze} = {n_az*n_ze} directions)...",
          flush=True)

    def _ob_progress(done, total):
        pct = 100 * done / total
        print(f"  Overburden: {done}/{total} ({pct:.0f}%)", flush=True)

    az_c, ze_c, overburden, open_sky = compute_overburden_map(
        elev, transform, det_lat, det_lon, det_alt_m,
        rho, n_az, n_ze, ze_max_deg, step_m,
        progress_cb=_ob_progress,
    )

    write_overburden_map(az_c, ze_c, overburden, open_sky, out_ob,
                          det_lat, det_lon, det_alt_m, rho)
    print(f"  Overburden map written: {out_ob}", flush=True)

    # 3. Flux map
    print(f"  Computing directional flux map...", flush=True)

    def _fl_progress(done, total):
        pct = 100 * done / total
        print(f"  Flux: {done}/{total} ({pct:.0f}%)", flush=True)

    script_dir = Path(__file__).resolve().parent
    flux_map, opensky_flux = compute_flux_map(
        az_c, ze_c, overburden, open_sky,
        rho, spec_mode, mode,
        script_dir=script_dir,
        progress_cb=_fl_progress,
    )

    write_flux_map(az_c, ze_c, flux_map, out_flux,
                   det_lat, det_lon, det_alt_m)
    print(f"  Flux map written: {out_flux}", flush=True)

    # 4. Transmission map
    T_sim = compute_transmission_map(flux_map, opensky_flux)
    write_transmission_map(az_c, ze_c, T_sim, out_trans,
                           det_lat, det_lon, det_alt_m, rho)
    print(f"  Transmission map written: {out_trans}", flush=True)

    # 5. Summary
    elapsed = time.time() - t0
    write_summary(az_c, ze_c, overburden, flux_map, open_sky,
                  out_summ, det_lat, det_lon, det_alt_m, rho,
                  spec_mode, elapsed)
    print(f"  Summary written: {out_summ}", flush=True)
    print(f"  Completed in {elapsed:.1f} s", flush=True)


if __name__ == "__main__":
    main()
