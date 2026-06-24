"""
gui_geometry_engine.py  —  UCLouvain Muography Group
Synthetic Terrain Geometry Builder

Generates analytical terrain geometries as synthetic DEMs compatible with
the existing cosmoaleph_terrain_driver ray-tracing infrastructure.

PUBLIC API
──────────
    SyntheticDEM              — dataclass holding (elev, transform, metadata)
    build_synthetic_dem(...)  — dispatcher: shape_id → SyntheticDEM
    render_geometry_builder() — Streamlit UI; returns SyntheticDEM | None

INTEGRATION
──────────────────────────────────────────────────────────────────────────────
In gui_terrain_engine.py, replace the GeoTIFF-load block with:

    from gui_geometry_engine import render_geometry_builder, SyntheticDEM

    dem_source = st.radio("DEM source", ["📂 File", "🔷 Synthetic Geometry"])
    if dem_source == "🔷 Synthetic Geometry":
        synth = render_geometry_builder(det_lat, det_lon, det_alt_m)
        if synth:
            elev      = synth.elev
            transform = synth.transform
    else:
        elev, transform = load_dem(dem_file_path)   # existing code

COORDINATE CONVENTION
──────────────────────
All shapes are centred on the detector GPS position (lat0, lon0).
The grid uses a flat-Earth approximation: 1° lat ≈ 111 320 m.
This is accurate to ≪1 m over the ~30 km extent used in practice.

Grid coordinate system:
    x  →  East  (metres from detector)
    y  →  North (metres from detector)
    z  →  elevation above sea level (metres)

Author : Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026
"""

from __future__ import annotations

import io
import tempfile
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Tuple

import numpy as np
import streamlit as st

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_M_PER_DEG_LAT = 111_320.0  # metres per degree of latitude (WGS84 mean)
_VERSION       = "1.0.0"

# Dark-theme Plotly layout (shared with gui_terrain_engine.py)
DARK = dict(
    paper_bgcolor="rgb(15,17,23)",
    plot_bgcolor="rgb(20,22,30)",
    font=dict(color="white", size=11),
)


# ─────────────────────────────────────────────────────────────────────────────
# Data container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SyntheticDEM:
    """
    Lightweight container for a synthetic DEM compatible with
    cosmoaleph_terrain_driver functions (dem_elevation_at, ray_overburden,
    compute_overburden_map).

    Fields
    ──────
    elev       : float32 array (nrows × ncols), elevation [m a.s.l.],
                 rows ordered North→South (rasterio convention).
    transform  : rasterio AffineTransform mapping (col, row) → (lon, lat).
    shape_id   : human-readable shape label.
    params     : dict of user-supplied parameters (for export / logging).
    det_lat    : detector latitude  used as grid centre [°].
    det_lon    : detector longitude used as grid centre [°].
    det_alt_m  : detector altitude [m] — used to clamp elev if needed.
    extent_km  : half-extent of the grid [km] (grid spans ±extent_km from det).
    res_m      : pixel resolution [m].
    """
    elev      : np.ndarray
    transform : object              # rasterio.transform.Affine
    shape_id  : str  = "unknown"
    params    : dict = field(default_factory=dict)
    det_lat   : float = 0.0
    det_lon   : float = 0.0
    det_alt_m : float = 0.0
    extent_km : float = 10.0
    res_m     : float = 50.0

    @property
    def nrows(self):  return self.elev.shape[0]
    @property
    def ncols(self):  return self.elev.shape[1]
    @property
    def z_min(self):  return float(np.nanmin(self.elev))
    @property
    def z_max(self):  return float(np.nanmax(self.elev))

    def summary(self) -> str:
        return (
            f"{self.shape_id}  {self.ncols}×{self.nrows} px  "
            f"res={self.res_m:.0f} m  extent=±{self.extent_km:.1f} km  "
            f"z=[{self.z_min:.0f}, {self.z_max:.0f}] m"
        )

    def export_xyz(self) -> str:
        """Return XYZ-format string (lon  lat  elev) for manual inspection."""
        import rasterio.transform as _rt
        lines = ["# lon[deg]  lat[deg]  elev[m]"]
        for row in range(self.nrows):
            for col in range(self.ncols):
                lon, lat = _rt.xy(self.transform, row, col, offset="center")
                lines.append(f"{lon:.6f}  {lat:.6f}  {self.elev[row, col]:.2f}")
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Grid factory
# ─────────────────────────────────────────────────────────────────────────────

def _make_grid(det_lat: float, det_lon: float,
               extent_km: float, res_m: float
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, object]:
    """
    Build a regular E-N grid centred on (det_lat, det_lon).

    Returns
    ───────
    x_e  : 1-D East offsets from detector [m], shape (ncols,)
    y_n  : 1-D North offsets from detector [m], shape (nrows,)
    XX   : 2-D East  grid [m], shape (nrows, ncols)
    YY   : 2-D North grid [m], shape (nrows, ncols)
    transform : rasterio AffineTransform (lon/lat ←→ col/row)
    """
    import rasterio.transform as _rt

    half = extent_km * 1_000.0          # metres
    x_e  = np.arange(-half, half + res_m, res_m)
    y_n  = np.arange(-half, half + res_m, res_m)
    XX, YY = np.meshgrid(x_e, y_n)     # (nrows, ncols)

    # Convert ENU offsets to lon/lat
    cos_lat  = max(np.cos(np.radians(det_lat)), 1e-6)
    lon_grid = det_lon + XX / (_M_PER_DEG_LAT * cos_lat)
    lat_grid = det_lat + YY / _M_PER_DEG_LAT

    # rasterio AffineTransform from the extent in lon/lat
    ncols = XX.shape[1]
    nrows = XX.shape[0]
    dx_lon = (lon_grid[0, -1] - lon_grid[0, 0]) / (ncols - 1)  if ncols > 1 else 1e-5
    dy_lat = (lat_grid[-1, 0] - lat_grid[0, 0]) / (nrows - 1)  if nrows > 1 else 1e-5

    # Origin: top-left corner (min lon, max lat) — rasterio N→S row ordering
    transform = _rt.from_origin(
        lon_grid[0, 0]  - dx_lon / 2,   # west edge
        lat_grid[-1, 0] + abs(dy_lat) / 2,  # north edge
        abs(dx_lon),                     # pixel width in lon
        abs(dy_lat),                     # pixel height in lat (positive)
    )

    # Flip Y so that row 0 = North (rasterio convention)
    return x_e, y_n, XX, YY[::-1], transform


def _wrap(elev_2d: np.ndarray, transform, shape_id: str, params: dict,
          det_lat, det_lon, det_alt_m, extent_km, res_m) -> SyntheticDEM:
    """Wrap a raw elevation array into a SyntheticDEM."""
    return SyntheticDEM(
        elev      = elev_2d.astype(np.float32),
        transform = transform,
        shape_id  = shape_id,
        params    = params,
        det_lat   = det_lat,
        det_lon   = det_lon,
        det_alt_m = det_alt_m,
        extent_km = extent_km,
        res_m     = res_m,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Shape generators
# Each function signature:
#   build_<shape>(det_lat, det_lon, det_alt_m, extent_km, res_m, **kwargs)
#                → SyntheticDEM
# ─────────────────────────────────────────────────────────────────────────────

def build_flat_slab(det_lat, det_lon, det_alt_m,
                    extent_km=15.0, res_m=50.0,
                    z_surface=0.0, **_) -> SyntheticDEM:
    """
    Infinite flat slab at constant elevation z_surface [m].
    Useful as a reference / sanity check against flat-earth engines.
    """
    _, _, XX, _, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    elev = np.full_like(XX, fill_value=z_surface, dtype=np.float64)
    return _wrap(elev, transform, "Flat Slab",
                 dict(z_surface=z_surface),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_cone(det_lat, det_lon, det_alt_m,
               extent_km=15.0, res_m=50.0,
               z_base=0.0, z_apex=1281.0, r_base_m=5_500.0,
               offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Right circular cone — classic idealised volcano model.

    Parameters
    ──────────
    z_base    : elevation of the cone base (flat surroundings) [m]
    z_apex    : elevation of the cone summit [m]
    r_base_m  : radius of the cone base [m]
    offset_e_m: East  offset of cone axis from detector [m]
    offset_n_m: North offset of cone axis from detector [m]
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    r = np.sqrt((XX - offset_e_m)**2 + (YY - offset_n_m)**2)
    h = (z_apex - z_base) * np.maximum(0.0, 1.0 - r / r_base_m)
    elev = z_base + h
    return _wrap(elev, transform, "Cone",
                 dict(z_base=z_base, z_apex=z_apex, r_base_m=r_base_m,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_truncated_cone(det_lat, det_lon, det_alt_m,
                         extent_km=15.0, res_m=50.0,
                         z_base=0.0, z_top=1281.0,
                         r_base_m=5_500.0, r_top_m=600.0,
                         offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Truncated cone (frustum) — caldera-forming volcano with flat top.

    Parameters
    ──────────
    z_base    : elevation of the surrounding plain [m]
    z_top     : elevation of the flat caldera floor / rim [m]
    r_base_m  : base radius [m]
    r_top_m   : caldera/summit radius [m]  (< r_base_m)
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    r = np.sqrt((XX - offset_e_m)**2 + (YY - offset_n_m)**2)
    # Three zones: flat top (r ≤ r_top), linear flank, flat base (r ≥ r_base)
    t = np.clip((r_base_m - r) / max(r_base_m - r_top_m, 1.0), 0.0, 1.0)
    elev = z_base + t * (z_top - z_base)
    return _wrap(elev, transform, "Truncated Cone",
                 dict(z_base=z_base, z_top=z_top,
                      r_base_m=r_base_m, r_top_m=r_top_m,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_gaussian_mound(det_lat, det_lon, det_alt_m,
                         extent_km=15.0, res_m=50.0,
                         z_base=0.0, amplitude_m=800.0, sigma_m=2_000.0,
                         offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Circularly-symmetric Gaussian hill: z(r) = z_base + A·exp(−r²/2σ²).

    This is the smoothest possible mound — useful for benchmarking because
    the cross-section at any zenith angle has an analytic closed form.

    Parameters
    ──────────
    amplitude_m : peak elevation above z_base [m]
    sigma_m     : Gaussian half-width (1σ) [m]
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    r2  = (XX - offset_e_m)**2 + (YY - offset_n_m)**2
    elev = z_base + amplitude_m * np.exp(-r2 / (2.0 * sigma_m**2))
    return _wrap(elev, transform, "Gaussian Mound",
                 dict(z_base=z_base, amplitude_m=amplitude_m, sigma_m=sigma_m,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_elliptical_mound(det_lat, det_lon, det_alt_m,
                           extent_km=15.0, res_m=50.0,
                           z_base=0.0, amplitude_m=800.0,
                           sigma_e_m=2_000.0, sigma_n_m=3_500.0,
                           rotation_deg=0.0,
                           offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Elliptical Gaussian mound — broken E/N symmetry.

    Useful for elongated geological structures (ridges, salt domes, anticlines).

    Parameters
    ──────────
    sigma_e_m    : East  Gaussian half-width [m]
    sigma_n_m    : North Gaussian half-width [m]
    rotation_deg : CCW rotation of the ellipse principal axes [°]
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    dx = XX - offset_e_m
    dy = YY - offset_n_m
    phi = np.radians(rotation_deg)
    # Rotate coordinates into ellipse frame
    xr = dx * np.cos(phi) + dy * np.sin(phi)
    yr = -dx * np.sin(phi) + dy * np.cos(phi)
    arg = (xr / sigma_e_m)**2 + (yr / sigma_n_m)**2
    elev = z_base + amplitude_m * np.exp(-0.5 * arg)
    return _wrap(elev, transform, "Elliptical Mound",
                 dict(z_base=z_base, amplitude_m=amplitude_m,
                      sigma_e_m=sigma_e_m, sigma_n_m=sigma_n_m,
                      rotation_deg=rotation_deg,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_hemisphere(det_lat, det_lon, det_alt_m,
                     extent_km=15.0, res_m=50.0,
                     z_base=0.0, radius_m=1_500.0,
                     offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Hemispherical dome: z(r) = z_base + sqrt(R² − r²)  for r ≤ R, else z_base.

    Useful for testing — the overburden along any chord through the sphere has
    an analytic solution: L = 2·sqrt(R² − d²) where d is the impact parameter.

    Parameters
    ──────────
    radius_m  : radius of the hemisphere [m]
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    r2  = (XX - offset_e_m)**2 + (YY - offset_n_m)**2
    R2  = radius_m**2
    inside = r2 <= R2
    elev   = np.where(inside, z_base + np.sqrt(np.maximum(R2 - r2, 0.0)), z_base)
    return _wrap(elev, transform, "Hemisphere",
                 dict(z_base=z_base, radius_m=radius_m,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_box(det_lat, det_lon, det_alt_m,
              extent_km=15.0, res_m=50.0,
              z_base=0.0, z_top=500.0,
              half_e_m=1_000.0, half_n_m=1_000.0,
              offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Rectangular box prism — flat-topped overburden block.

    Useful for industrial/lab scenarios (e.g. CCS storage reservoir overburden),
    underground cavern geometry, or as a worst-case flat-top test.

    Parameters
    ──────────
    z_base    : surrounding terrain elevation [m]
    z_top     : top surface of the box [m]
    half_e_m  : East  half-length of the box [m]
    half_n_m  : North half-length of the box [m]
    """
    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    in_box = (np.abs(XX - offset_e_m) <= half_e_m) & \
             (np.abs(YY - offset_n_m) <= half_n_m)
    elev   = np.where(in_box, z_top, z_base)
    return _wrap(elev, transform, "Box",
                 dict(z_base=z_base, z_top=z_top,
                      half_e_m=half_e_m, half_n_m=half_n_m,
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


def build_custom_radial(det_lat, det_lon, det_alt_m,
                        extent_km=15.0, res_m=50.0,
                        z_base=0.0,
                        r_pts: np.ndarray = None,
                        z_pts: np.ndarray = None,
                        offset_e_m=0.0, offset_n_m=0.0, **_) -> SyntheticDEM:
    """
    Custom rotationally symmetric profile defined by (r, z) tabulated points.

    The profile is interpolated with cubic splines and extrapolated flat
    (at z_base) beyond the last tabulated radius.

    Parameters
    ──────────
    r_pts  : 1-D array of radii [m], must start at r=0 and be strictly increasing
    z_pts  : 1-D array of elevations [m] at each radius
    z_base : elevation beyond the last tabulated radius [m]
    """
    from scipy.interpolate import interp1d

    if r_pts is None or z_pts is None or len(r_pts) < 2:
        raise ValueError("Custom profile requires at least 2 (r, z) points.")

    r_pts = np.asarray(r_pts, dtype=float)
    z_pts = np.asarray(z_pts, dtype=float)

    # Clamp: beyond max radius → z_base
    r_max = r_pts[-1]
    interp = interp1d(r_pts, z_pts, kind="cubic",
                      bounds_error=False, fill_value=z_base)

    _, _, XX, YY, transform = _make_grid(det_lat, det_lon, extent_km, res_m)
    r    = np.sqrt((XX - offset_e_m)**2 + (YY - offset_n_m)**2)
    elev = interp(r)
    elev[r > r_max] = z_base

    return _wrap(elev, transform, "Custom Radial Profile",
                 dict(z_base=z_base, r_pts=r_pts.tolist(), z_pts=z_pts.tolist(),
                      offset_e_m=offset_e_m, offset_n_m=offset_n_m),
                 det_lat, det_lon, det_alt_m, extent_km, res_m)


# ─────────────────────────────────────────────────────────────────────────────
# Dispatcher
# ─────────────────────────────────────────────────────────────────────────────

SHAPE_BUILDERS: dict[str, Callable] = {
    "Flat Slab"          : build_flat_slab,
    "Cone"               : build_cone,
    "Truncated Cone"     : build_truncated_cone,
    "Gaussian Mound"     : build_gaussian_mound,
    "Elliptical Mound"   : build_elliptical_mound,
    "Hemisphere"         : build_hemisphere,
    "Box"                : build_box,
    "Custom Radial Profile": build_custom_radial,
}

SHAPE_LABELS = list(SHAPE_BUILDERS.keys())


def build_synthetic_dem(shape_id: str, det_lat: float, det_lon: float,
                        det_alt_m: float, extent_km: float, res_m: float,
                        **kwargs) -> SyntheticDEM:
    """
    Build a SyntheticDEM by shape name.

    Parameters
    ──────────
    shape_id  : one of SHAPE_LABELS
    det_lat   : detector latitude  [°]
    det_lon   : detector longitude [°]
    det_alt_m : detector altitude  [m]
    extent_km : half-extent of grid [km] (grid spans ±extent_km)
    res_m     : pixel resolution [m]
    **kwargs  : shape-specific parameters (see individual build_* functions)
    """
    if shape_id not in SHAPE_BUILDERS:
        raise ValueError(f"Unknown shape '{shape_id}'. Choose from: {SHAPE_LABELS}")
    return SHAPE_BUILDERS[shape_id](
        det_lat, det_lon, det_alt_m, extent_km, res_m, **kwargs
    )


# ─────────────────────────────────────────────────────────────────────────────
# Plotly 3D preview
# ─────────────────────────────────────────────────────────────────────────────

def _plotly_3d_preview(synth: SyntheticDEM, det_alt_m: float,
                       subsample: int = 4):
    """
    Return a Plotly Figure with a 3-D surface preview of the synthetic DEM
    and a marker for the detector position.
    Uses the same dark theme as the rest of the GUI.

    subsample: stride for decimation (1 = full resolution, 4 = 1/4 pixels).
    """
    import plotly.graph_objects as go

    # Subsample for speed — full resolution is often overkill for preview
    s  = max(1, int(subsample))
    e  = synth.elev[::s, ::s]           # (nrows', ncols')

    # Build x, y in km from detector centre (for human-readable axis labels)
    half_px = synth.res_m * s / 2
    nx = e.shape[1]
    ny = e.shape[0]
    x_km = np.linspace(-synth.extent_km, synth.extent_km, nx)
    y_km = np.linspace(-synth.extent_km, synth.extent_km, ny)

    fig = go.Figure()

    fig.add_trace(go.Surface(
        z=e,
        x=x_km,
        y=y_km,
        colorscale="Viridis",
        colorbar=dict(title="Elev [m]", thickness=12),
        showscale=True,
        opacity=0.92,
        name=synth.shape_id,
    ))

    # Detector marker
    fig.add_trace(go.Scatter3d(
        x=[0.0], y=[0.0], z=[det_alt_m],
        mode="markers",
        marker=dict(size=8, color="red", symbol="diamond"),
        name="Detector",
    ))

    fig.update_layout(
        **DARK,
        height=450,
        margin=dict(l=0, r=0, t=30, b=0),
        scene=dict(
            xaxis_title="East [km]",
            yaxis_title="North [km]",
            zaxis_title="Elevation [m]",
            bgcolor="rgb(20,22,30)",
            xaxis=dict(gridcolor="rgba(255,255,255,0.1)"),
            yaxis=dict(gridcolor="rgba(255,255,255,0.1)"),
            zaxis=dict(gridcolor="rgba(255,255,255,0.1)"),
            camera=dict(eye=dict(x=1.4, y=-1.6, z=0.9)),
        ),
        title=dict(text=f"🔷 {synth.shape_id} — synthetic DEM preview", x=0.0,
                   font=dict(color="white", size=13)),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Streamlit UI
# ─────────────────────────────────────────────────────────────────────────────

def render_geometry_builder(det_lat: float, det_lon: float,
                             det_alt_m: float) -> Optional[SyntheticDEM]:
    """
    Render the geometry builder UI inside the current Streamlit container.

    Call this inside the terrain tab when the user selects "Synthetic Geometry"
    as the DEM source.  Returns a SyntheticDEM when the user has configured
    a valid shape and pressed "Build DEM", or None otherwise.

    Parameters
    ──────────
    det_lat   : detector latitude  (from Tab 5 GPS fields) [°]
    det_lon   : detector longitude (from Tab 5 GPS fields) [°]
    det_alt_m : detector altitude  (from Tab 5 GPS fields) [m]
    """

    st.markdown("#### 🔷 Synthetic Terrain Geometry")
    st.caption(
        f"Grid centred on detector at ({det_lat:.4f}°, {det_lon:.4f}°, "
        f"{det_alt_m:.0f} m a.s.l.) — flat-Earth approximation"
    )

    # ── Grid parameters ───────────────────────────────────────────────────────
    with st.expander("⚙️ Grid parameters", expanded=False):
        gc1, gc2 = st.columns(2)
        with gc1:
            extent_km = st.number_input(
                "Grid half-extent [km]", min_value=1.0, max_value=100.0,
                value=15.0, step=1.0, key="geom_extent_km",
                help="Grid spans ±extent_km around the detector. "
                     "Use ≥ 2× the largest feature radius."
            )
        with gc2:
            res_m = st.number_input(
                "Resolution [m/pixel]", min_value=10.0, max_value=500.0,
                value=50.0, step=10.0, key="geom_res_m",
                help="Ray-trace step in terrain_driver is typically 50 m; "
                     "DEM resolution should be ≤ step size."
            )
        n_px = int(2 * extent_km * 1000 / res_m) + 1
        st.info(f"Grid size: {n_px} × {n_px} = {n_px**2 / 1e6:.1f} M pixels  "
                f"({n_px**2 * 4 / 1e6:.0f} MB RAM)")

    # ── Shape selector ────────────────────────────────────────────────────────
    shape_id = st.selectbox(
        "Geometry shape", options=SHAPE_LABELS, key="geom_shape_id",
        help="Select the analytical terrain model."
    )

    # ── Per-shape parameter UI ─────────────────────────────────────────────────
    params: dict = {}

    # ── Cone ──────────────────────────────────────────────────────────────────
    if shape_id == "Cone":
        st.markdown("**Idealised volcano:** linear flanks, pointed summit.")
        ca, cb, cc = st.columns(3)
        with ca:
            params["z_base"]   = st.number_input("Base elevation [m]", value=0.0,  step=10.0, key="cone_zbase")
            params["z_apex"]   = st.number_input("Summit elevation [m]", value=1281.0, step=10.0, key="cone_zapex")
        with cb:
            params["r_base_m"] = st.number_input("Base radius [m]", value=5500.0, step=100.0, key="cone_rbase")
        with cc:
            params["offset_e_m"] = st.number_input("Axis offset East [m]",  value=0.0, step=100.0, key="cone_oe")
            params["offset_n_m"] = st.number_input("Axis offset North [m]", value=0.0, step=100.0, key="cone_on")
        # Sanity check
        if params["z_apex"] <= params["z_base"]:
            st.warning("Summit elevation must be above base elevation.")
        slope_deg = np.degrees(np.arctan((params["z_apex"] - params["z_base"]) /
                                          max(params["r_base_m"], 1.0)))
        st.caption(f"Flank slope ≈ {slope_deg:.1f}°  |  "
                   f"Volume ≈ {np.pi/3 * params['r_base_m']**2 * (params['z_apex']-params['z_base']) / 1e9:.2f} km³")

    # ── Truncated Cone ────────────────────────────────────────────────────────
    elif shape_id == "Truncated Cone":
        st.markdown("**Frustum / caldera volcano:** flat-topped, linear flanks.")
        ta, tb, tc = st.columns(3)
        with ta:
            params["z_base"]   = st.number_input("Base elevation [m]", value=0.0,  step=10.0, key="tc_zbase")
            params["z_top"]    = st.number_input("Top elevation [m]",  value=1281.0,step=10.0, key="tc_ztop")
        with tb:
            params["r_base_m"] = st.number_input("Base radius [m]", value=5500.0, step=100.0, key="tc_rbase")
            params["r_top_m"]  = st.number_input("Top (caldera) radius [m]", value=600.0, step=50.0, key="tc_rtop")
        with tc:
            params["offset_e_m"] = st.number_input("Axis East [m]",  value=0.0, step=100.0, key="tc_oe")
            params["offset_n_m"] = st.number_input("Axis North [m]", value=0.0, step=100.0, key="tc_on")
        if params["r_top_m"] >= params["r_base_m"]:
            st.warning("Top radius must be smaller than base radius.")

    # ── Gaussian Mound ────────────────────────────────────────────────────────
    elif shape_id == "Gaussian Mound":
        st.markdown("**Smooth hill:** z(r) = z₀ + A·exp(−r²/2σ²).  Good for analytic cross-checks.")
        ga, gb = st.columns(2)
        with ga:
            params["z_base"]       = st.number_input("Base elevation [m]",   value=0.0, step=10.0, key="gm_zbase")
            params["amplitude_m"]  = st.number_input("Amplitude A [m]",      value=800.0,step=10.0, key="gm_amp")
        with gb:
            params["sigma_m"]      = st.number_input("Width σ [m]",          value=2000.0, step=100.0, key="gm_sig")
            params["offset_e_m"]   = st.number_input("Centre East [m]",      value=0.0, step=100.0, key="gm_oe")
            params["offset_n_m"]   = st.number_input("Centre North [m]",     value=0.0, step=100.0, key="gm_on")
        fwhm = params["sigma_m"] * 2 * np.sqrt(2 * np.log(2))
        st.caption(f"FWHM ≈ {fwhm/1000:.2f} km  |  Peak at {params['z_base']+params['amplitude_m']:.0f} m a.s.l.")

    # ── Elliptical Mound ──────────────────────────────────────────────────────
    elif shape_id == "Elliptical Mound":
        st.markdown("**Asymmetric hill:** independent E/N widths + rotation.")
        ea, eb = st.columns(2)
        with ea:
            params["z_base"]      = st.number_input("Base elevation [m]",   value=0.0,   step=10.0,  key="em_zbase")
            params["amplitude_m"] = st.number_input("Amplitude A [m]",      value=800.0, step=10.0,  key="em_amp")
            params["sigma_e_m"]   = st.number_input("East width σ_E [m]",   value=2000.0,step=100.0, key="em_se")
        with eb:
            params["sigma_n_m"]   = st.number_input("North width σ_N [m]",  value=3500.0,step=100.0, key="em_sn")
            params["rotation_deg"]= st.slider("Rotation CCW [°]", 0.0, 180.0, 0.0, 5.0, key="em_rot")
            params["offset_e_m"]  = st.number_input("Centre East [m]",      value=0.0,   step=100.0, key="em_oe")
            params["offset_n_m"]  = st.number_input("Centre North [m]",     value=0.0,   step=100.0, key="em_on")

    # ── Hemisphere ────────────────────────────────────────────────────────────
    elif shape_id == "Hemisphere":
        st.markdown("**Hemispherical dome.**  Path length through a sphere: L = 2√(R²−d²)  "
                    "(analytic — useful for engine validation).")
        ha, hb = st.columns(2)
        with ha:
            params["z_base"]    = st.number_input("Base elevation [m]", value=0.0,    step=10.0,  key="hs_zbase")
            params["radius_m"]  = st.number_input("Radius R [m]",       value=1500.0, step=100.0, key="hs_rad")
        with hb:
            params["offset_e_m"] = st.number_input("Centre East [m]",   value=0.0, step=100.0, key="hs_oe")
            params["offset_n_m"] = st.number_input("Centre North [m]",  value=0.0, step=100.0, key="hs_on")
        vol = 2/3 * np.pi * params["radius_m"]**3 / 1e9
        st.caption(f"Height = {params['radius_m']:.0f} m  |  Volume = {vol:.2f} km³")

    # ── Box ───────────────────────────────────────────────────────────────────
    elif shape_id == "Box":
        st.markdown("**Rectangular prism** — flat-topped overburden block for CCS / lab scenarios.")
        ba, bb = st.columns(2)
        with ba:
            params["z_base"]    = st.number_input("Surrounding elevation [m]", value=0.0,   step=10.0,  key="bx_zbase")
            params["z_top"]     = st.number_input("Top elevation [m]",         value=500.0, step=10.0,  key="bx_ztop")
        with bb:
            params["half_e_m"]  = st.number_input("East  half-length [m]",  value=1000.0, step=100.0, key="bx_he")
            params["half_n_m"]  = st.number_input("North half-length [m]",  value=1000.0, step=100.0, key="bx_hn")
            params["offset_e_m"]= st.number_input("Centre East [m]",        value=0.0, step=100.0, key="bx_oe")
            params["offset_n_m"]= st.number_input("Centre North [m]",       value=0.0, step=100.0, key="bx_on")
        thickness = params["z_top"] - params["z_base"]
        vol = 4 * params["half_e_m"] * params["half_n_m"] * max(thickness, 0) / 1e9
        st.caption(f"Thickness = {thickness:.0f} m  |  Footprint = "
                   f"{2*params['half_e_m']/1000:.1f} × {2*params['half_n_m']/1000:.1f} km  |  "
                   f"Volume ≈ {vol:.3f} km³")

    # ── Flat Slab ─────────────────────────────────────────────────────────────
    elif shape_id == "Flat Slab":
        st.markdown("**Uniform flat terrain** — reference / sanity check.  "
                    "Results must match the Transport-tab flat-earth engines at the same depth.")
        params["z_surface"] = st.number_input("Surface elevation [m]", value=0.0, step=10.0, key="fs_z")

    # ── Custom Radial Profile ─────────────────────────────────────────────────
    elif shape_id == "Custom Radial Profile":
        st.markdown("**Tabulated radial profile** — define (r, z) pairs; "
                    "cubic spline interpolation; flat beyond max radius.")
        col_csv, col_hint = st.columns([2, 1])
        with col_hint:
            st.markdown(textwrap.dedent("""
                **CSV format:**
                ```
                # r[m]   z[m]
                0,       1281
                500,     1250
                2000,    900
                5000,    200
                8000,    0
                ```
            """))
        with col_csv:
            csv_text = st.text_area(
                "Profile table (r [m], z [m], one point per line)",
                value="0,1281\n500,1250\n2000,900\n5000,200\n8000,0",
                height=180, key="geom_csv_profile",
            )
        params["z_base"]    = st.number_input("Elevation beyond max radius [m]", value=0.0, step=10.0, key="cr_zbase")
        params["offset_e_m"]= st.number_input("Centre East [m]",  value=0.0, step=100.0, key="cr_oe")
        params["offset_n_m"]= st.number_input("Centre North [m]", value=0.0, step=100.0, key="cr_on")

        # Parse CSV
        r_pts, z_pts = [], []
        try:
            for line in csv_text.strip().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                vals = [float(v) for v in line.replace(";", ",").split(",")]
                if len(vals) >= 2:
                    r_pts.append(vals[0])
                    z_pts.append(vals[1])
            r_pts = np.asarray(r_pts)
            z_pts = np.asarray(z_pts)
            if len(r_pts) < 2:
                st.warning("Need at least 2 (r, z) points.")
            elif not np.all(np.diff(r_pts) > 0):
                st.warning("Radii must be strictly increasing.")
            else:
                params["r_pts"] = r_pts
                params["z_pts"] = z_pts
                st.success(f"Parsed {len(r_pts)} profile points: "
                           f"r ∈ [0, {r_pts[-1]:.0f}] m, "
                           f"z ∈ [{z_pts.min():.0f}, {z_pts.max():.0f}] m")
        except Exception as exc:
            st.error(f"CSV parse error: {exc}")

    # ── Build button ──────────────────────────────────────────────────────────
    st.divider()
    build_col, clear_col = st.columns([3, 1])
    with clear_col:
        if st.button("🗑 Clear", key="geom_clear", width='stretch'):
            st.session_state.pop("_synth_dem", None)
            st.rerun()
    with build_col:
        do_build = st.button("🔷 Build Synthetic DEM", type="primary",
                             key="geom_build", width='stretch')

    if do_build:
        # Validate minimum requirements
        if shape_id == "Custom Radial Profile" and (
                "r_pts" not in params or len(params.get("r_pts", [])) < 2):
            st.error("Fix the profile table before building.")
            return st.session_state.get("_synth_dem")

        with st.spinner("Generating synthetic DEM…"):
            try:
                synth = build_synthetic_dem(
                    shape_id, det_lat, det_lon, det_alt_m,
                    float(st.session_state.get("geom_extent_km", extent_km)),
                    float(st.session_state.get("geom_res_m",    res_m)),
                    **params,
                )
                st.session_state["_synth_dem"] = synth
            except Exception as exc:
                st.error(f"DEM generation failed: {exc}")
                return st.session_state.get("_synth_dem")

    # ── Preview ───────────────────────────────────────────────────────────────
    synth: Optional[SyntheticDEM] = st.session_state.get("_synth_dem")

    if synth is not None:
        st.success(f"✅  {synth.summary()}")

        # 3D preview
        try:
            sub = max(1, int(synth.nrows / 128))   # target ≤128² surface pts
            fig = _plotly_3d_preview(synth, det_alt_m, subsample=sub)
            st.plotly_chart(fig)
        except Exception as exc:
            st.warning(f"3D preview failed: {exc}")

        # Overburden profile along NS/EW transects
        with st.expander("📊 Elevation transects (QC)", expanded=False):
            _render_transects(synth)

        # Export
        with st.expander("💾 Export", expanded=False):
            _render_export(synth)

    return synth


# ─────────────────────────────────────────────────────────────────────────────
# QC helpers
# ─────────────────────────────────────────────────────────────────────────────

def _render_transects(synth: SyntheticDEM):
    """Plot N–S and E–W elevation cross-sections through the detector."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    mid_col = synth.ncols // 2
    mid_row = synth.nrows // 2

    # E–W transect (middle row in N→S ordering → centre)
    ew_row   = synth.nrows - 1 - mid_row   # flip because elev rows are N→S
    ew_elev  = synth.elev[ew_row, :]
    x_km     = np.linspace(-synth.extent_km, synth.extent_km, synth.ncols)

    # N–S transect (middle column)
    ns_elev  = synth.elev[::-1, mid_col]   # flip to get S→N
    y_km     = np.linspace(-synth.extent_km, synth.extent_km, synth.nrows)

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=["E–W transect (through detector)",
                                        "N–S transect (through detector)"])

    fig.add_trace(go.Scatter(x=x_km, y=ew_elev, mode="lines",
                             line=dict(color="#4fc3f7", width=2), name="E–W"),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=y_km, y=ns_elev, mode="lines",
                             line=dict(color="#81c784", width=2), name="N–S"),
                  row=1, col=2)

    # Detector altitude line
    for col_i in (1, 2):
        xr = x_km if col_i == 1 else y_km
        fig.add_hline(y=synth.det_alt_m, line_dash="dash",
                      line_color="red", row=1, col=col_i,
                      annotation_text="det", annotation_font_color="red")

    fig.update_layout(**DARK, height=300, showlegend=False,
                      margin=dict(l=0, r=0, t=30, b=0))
    fig.update_xaxes(title_text="Distance [km]", gridcolor="rgba(255,255,255,0.1)")
    fig.update_yaxes(title_text="Elevation [m]", gridcolor="rgba(255,255,255,0.1)")
    st.plotly_chart(fig)


def _render_export(synth: SyntheticDEM):
    """Download buttons: XYZ point cloud and parameter JSON."""
    import json as _json

    ec1, ec2 = st.columns(2)

    with ec1:
        xyz_str = synth.export_xyz()
        st.download_button(
            "📥 Download XYZ (lon lat elev)",
            data=xyz_str.encode(),
            file_name=f"synthetic_{synth.shape_id.replace(' ','_').lower()}.xyz",
            mime="text/plain",
            key="geom_dl_xyz",
            width='stretch',
        )

    with ec2:
        meta = {
            "shape_id" : synth.shape_id,
            "det_lat"  : synth.det_lat,
            "det_lon"  : synth.det_lon,
            "det_alt_m": synth.det_alt_m,
            "extent_km": synth.extent_km,
            "res_m"    : synth.res_m,
            "nrows"    : synth.nrows,
            "ncols"    : synth.ncols,
            "z_min_m"  : synth.z_min,
            "z_max_m"  : synth.z_max,
            "params"   : {k: (v.tolist() if hasattr(v, "tolist") else v)
                          for k, v in synth.params.items()},
        }
        st.download_button(
            "📥 Download parameters (JSON)",
            data=_json.dumps(meta, indent=2).encode(),
            file_name=f"synthetic_{synth.shape_id.replace(' ','_').lower()}.json",
            mime="application/json",
            key="geom_dl_json",
            width='stretch',
        )
