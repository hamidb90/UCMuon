"""
gui_csg_engine.py  —  UCLouvain Muography Group
CSG / Volumetric Geometry Engine for Underground Muography

Supports three geometry input formats:
  1. STL  (.stl)        — triangle mesh; loaded via trimesh
  2. PHITS input file   — [Surface] + [Cell] cards (CSG)
  3. MCNP input file    — surface + cell cards (CSG, similar syntax)

All formats produce a CSGGeometry object whose compute_overburden_map()
returns the same (az_c, ze_c, ob_map, sky_map) tuple consumed by
gui_terrain_engine.py — so all downstream transport engines (UCMuon-MC,
MUSIC, Bethe-Bloch) work without modification.

COORDINATE CONVENTION (internal)
──────────────────────────────────
  x = East  [m]    y = North  [m]    z = Up  [m]
  Azimuth: 0 = North, 90 = East (geographic, clockwise)
  Zenith:  0 = straight up, 90 = horizontal

  PHITS/MCNP files use cm.  All distances are converted to metres
  on load.  The user specifies the detector position in the geometry's
  own coordinate system (before conversion) using the GUI fields.

  STL files are assumed to be in metres unless the user ticks "cm".

DESIGN NOTES
─────────────
  • For STL: trimesh ray.intersects_location() gives exact mesh
    intersections analytically.  Path through the solid = sum of
    (exit−entry) distances for each intersection pair.

  • For PHITS/MCNP CSG: a step-march approach is used.  At each
    step_m interval along the ray we evaluate density_at(point) by
    testing point membership in each cell (sign of each surface
    function + boolean logic).  This is universal but slower than
    analytic ray-surface intersection.  Default step = 0.5 m inside
    small laboratory geometries; increase for large geological bodies.

  • The underground flag must be set in the GUI (Section 3 checkbox)
    so that cosmoaleph_terrain_driver.compute_overburden_map() does
    not clamp the detector altitude to the DEM surface.  For CSG
    geometries this is always True (no DEM is involved).

PUBLIC API
───────────
  CSGGeometry                      — main geometry container
  load_stl(path, density, scale)   → CSGGeometry
  load_phits(path)                 → CSGGeometry
  load_mcnp(path)                  → CSGGeometry
  render_csg_builder()             → CSGGeometry | None  (Streamlit UI)

Author : Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026
"""

from __future__ import annotations

import re
import io
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import streamlit as st

_VERSION = "1.0.0"

DARK = dict(
    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
    font=dict(color="white", size=11),
)

# ─────────────────────────────────────────────────────────────────────────────
# Surface algebra
# Every surface exposes  f(xyz) → float
# Convention (PHITS/MCNP): f < 0 → "inside" (negative region)
# ─────────────────────────────────────────────────────────────────────────────

class _SurfaceBase:
    """Abstract base for all analytical surfaces."""
    surf_id: int = 0
    def f(self, xyz: np.ndarray) -> float:
        raise NotImplementedError
    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.surf_id})"


class PlaneX(_SurfaceBase):
    def __init__(self, sid, x0): self.surf_id=sid; self.x0=x0
    def f(self, p): return p[0] - self.x0

class PlaneY(_SurfaceBase):
    def __init__(self, sid, y0): self.surf_id=sid; self.y0=y0
    def f(self, p): return p[1] - self.y0

class PlaneZ(_SurfaceBase):
    def __init__(self, sid, z0): self.surf_id=sid; self.z0=z0
    def f(self, p): return p[2] - self.z0

class Plane(_SurfaceBase):
    """General plane: ax + by + cz = d  (f = ax+by+cz-d)"""
    def __init__(self, sid, a, b, c, d):
        self.surf_id=sid; self.a=a; self.b=b; self.c=c; self.d=d
    def f(self, p): return self.a*p[0] + self.b*p[1] + self.c*p[2] - self.d

class Sphere(_SurfaceBase):
    def __init__(self, sid, x0, y0, z0, r):
        self.surf_id=sid; self.c=np.array([x0,y0,z0]); self.r2=r*r
    def f(self, p):
        d = p - self.c; return float(d@d) - self.r2

class CylinderZ(_SurfaceBase):
    """Infinite cylinder aligned with Z axis, radius r, axis at (x0,y0)."""
    def __init__(self, sid, x0, y0, r):
        self.surf_id=sid; self.x0=x0; self.y0=y0; self.r2=r*r
    def f(self, p): return (p[0]-self.x0)**2 + (p[1]-self.y0)**2 - self.r2

class CylinderX(_SurfaceBase):
    def __init__(self, sid, y0, z0, r):
        self.surf_id=sid; self.y0=y0; self.z0=z0; self.r2=r*r
    def f(self, p): return (p[1]-self.y0)**2 + (p[2]-self.z0)**2 - self.r2

class CylinderY(_SurfaceBase):
    def __init__(self, sid, x0, z0, r):
        self.surf_id=sid; self.x0=x0; self.z0=z0; self.r2=r*r
    def f(self, p): return (p[0]-self.x0)**2 + (p[2]-self.z0)**2 - self.r2

class ConeZ(_SurfaceBase):
    """Cone along Z: (x-x0)^2+(y-y0)^2 - t^2*(z-z0)^2 = 0"""
    def __init__(self, sid, x0, y0, z0, t2, sheet=0):
        self.surf_id=sid; self.x0=x0; self.y0=y0; self.z0=z0
        self.t2=t2; self.sheet=sheet
    def f(self, p):
        return (p[0]-self.x0)**2+(p[1]-self.y0)**2 - self.t2*(p[2]-self.z0)**2

class Quadric(_SurfaceBase):
    """General quadric: Ax²+By²+Cz²+Dxy+Eyz+Fxz+Gx+Hy+Iz+J=0"""
    def __init__(self, sid, coeffs):
        self.surf_id=sid
        self.A,self.B,self.C,self.D,self.E,self.F,self.G,self.H,self.I_,self.J = coeffs
    def f(self, p):
        x,y,z = p
        return (self.A*x*x + self.B*y*y + self.C*z*z +
                self.D*x*y + self.E*y*z + self.F*x*z +
                self.G*x   + self.H*y   + self.I_*z  + self.J)


# ─────────────────────────────────────────────────────────────────────────────
# Cell (region) — boolean combination of half-spaces
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CSGCell:
    """
    A material region defined by boolean combinations of surfaces.

    tokens : list of (sign: int, surf_id: int) OR the string 'OR'
             sign = -1 → inside (f < 0), +1 → outside (f > 0)
             'OR' is the union operator between groups
    density : g/cm³ (0 = void / outside world)
    mat_id  : material identifier (for multi-material geometries)
    cell_id : integer cell number from the input file
    """
    cell_id : int
    mat_id  : int
    density : float                      # g/cm³; 0 = void
    tokens  : list                       # [(sign, surf_id), ...] with 'OR' separators
    label   : str = ""

    def contains(self, point: np.ndarray,
                 surf_map: Dict[int, _SurfaceBase]) -> bool:
        """
        Return True if point is inside this cell.

        Each token group (separated by 'OR') is evaluated as an AND of
        half-space tests. The cell contains the point if ANY group is True.
        """
        group_result = True
        any_group_true = False

        for tok in self.tokens:
            if tok == 'OR':
                if group_result:
                    return True
                any_group_true |= group_result
                group_result = True
                continue
            sign, sid = tok
            surf = surf_map.get(abs(sid))
            if surf is None:
                continue
            fval = surf.f(point)
            # sign=-1 → inside (f<0); sign=+1 → outside (f>0)
            half_ok = (fval < 0) if sign < 0 else (fval > 0)
            group_result = group_result and half_ok

        return group_result


# ─────────────────────────────────────────────────────────────────────────────
# Main geometry container
# ─────────────────────────────────────────────────────────────────────────────

class CSGGeometry:
    """
    Volumetric geometry for underground muon transport.

    Contains either:
      A) A triangle mesh (from STL)  → uses trimesh ray intersection
      B) A set of CSG cells          → uses step-march density evaluation

    Public methods mirror cosmoaleph_terrain_driver so gui_terrain_engine
    can call them transparently.
    """

    def __init__(self):
        self.source_format: str = "unknown"   # "STL", "PHITS", "MCNP"
        self.source_path  : str = ""
        self.label        : str = "CSG Geometry"

        # STL path
        self._mesh = None            # trimesh.Trimesh or None
        self._stl_density: float = 2.65    # g/cm³ — uniform for STL

        # CSG path
        self._surf_map: Dict[int, _SurfaceBase] = {}
        self._cells   : List[CSGCell]           = []

        # Bounding box in metres (for step-march termination)
        self.bbox_min = np.array([-1e4, -1e4, -1e4])
        self.bbox_max = np.array([ 1e4,  1e4,  1e4])

        # Coordinate scale applied on load
        self._cm_to_m: bool = True   # True if source file uses cm

    # ── Info ─────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        if self._mesh is not None:
            return (f"STL  {self.source_path}  "
                    f"{len(self._mesh.faces)} triangles  "
                    f"ρ={self._stl_density:.2f} g/cm³")
        return (f"{self.source_format}  {self.source_path}  "
                f"{len(self._surf_map)} surfaces  "
                f"{len(self._cells)} cells")

    def materials(self) -> List[Tuple[int, float]]:
        """Return list of (mat_id, density) for non-void cells."""
        seen = {}
        for c in self._cells:
            if c.density > 0 and c.mat_id not in seen:
                seen[c.mat_id] = c.density
        return sorted(seen.items())

    # ── Point-density query ───────────────────────────────────────────────────

    def density_at(self, point_m: np.ndarray) -> float:
        """
        Return the material density [g/cm³] at point_m [m, ENU].
        Returns 0 for void / outside all cells.
        Cells are tested in order; first match wins.
        """
        if self._mesh is not None:
            # STL: uniform density inside the closed mesh
            try:
                if self._mesh.contains([point_m])[0]:
                    return self._stl_density
            except Exception:
                pass
            return 0.0

        for cell in self._cells:
            if cell.density <= 0:
                continue
            if cell.contains(point_m, self._surf_map):
                return cell.density
        return 0.0

    # ── Ray overburden ────────────────────────────────────────────────────────

    def ray_overburden(self,
                       az_deg: float, ze_deg: float,
                       det_pos_m: np.ndarray,
                       step_m: float = 0.5,
                       max_dist_m: float = 10_000.0) -> Tuple[float, float, bool]:
        """
        Compute rock overburden along one ray from the detector.

        Parameters
        ──────────
        az_deg    : azimuth [°], geographic convention (0=N, 90=E, CW)
        ze_deg    : zenith  [°] from vertical (0=up, 90=horiz)
        det_pos_m : detector ENU position [m] in the geometry frame
        step_m    : march step [m] (default 0.5 m — suitable for lab/mine)
        max_dist_m: maximum ray length [m]

        Returns (overburden_gcm2, slant_m, open_geometry)
        """
        az = np.radians(az_deg)
        ze = np.radians(ze_deg)
        sin_ze = np.sin(ze)
        # ENU unit direction
        direction = np.array([
            sin_ze * np.sin(az),   # East
            sin_ze * np.cos(az),   # North
            np.cos(ze),            # Up
        ], dtype=float)

        # ── STL: exact ray-mesh intersection ─────────────────────────────
        if self._mesh is not None:
            try:
                locs, ray_idx, _ = self._mesh.ray.intersects_location(
                    ray_origins    = [det_pos_m],
                    ray_directions = [direction],
                )
                if len(locs) < 2:
                    return 0.0, 0.0, True
                # Sort intersections by distance from detector
                dists = np.linalg.norm(locs - det_pos_m, axis=1)
                order = np.argsort(dists)
                dists = dists[order]
                # If detector is already inside the mesh, first intersection
                # is an exit; handle by prepending dist=0 as entry
                inside_start = self._mesh.contains([det_pos_m])[0]
                if inside_start:
                    dists = np.concatenate([[0.0], dists])
                # Pair (entry, exit): 0-1, 2-3, ...
                total_m = 0.0
                for i in range(0, len(dists) - 1, 2):
                    seg = dists[i+1] - dists[i]
                    if seg > 0:
                        total_m += seg
                if total_m <= 0:
                    return 0.0, 0.0, True
                ob = total_m * 100.0 * self._stl_density
                return ob, total_m, False
            except Exception:
                return 0.0, 0.0, True

        # ── CSG: step-march density integration ──────────────────────────
        total_m  = 0.0
        total_ob = 0.0
        pos = det_pos_m.copy().astype(float)

        for _ in range(int(max_dist_m / step_m) + 1):
            pos += direction * step_m
            # Terminate if outside bounding box
            if np.any(pos < self.bbox_min) or np.any(pos > self.bbox_max):
                break
            rho = self.density_at(pos)
            if rho > 0:
                total_m  += step_m
                total_ob += step_m * 100.0 * rho   # cm × g/cm³ = g/cm²

        if total_m <= 0:
            return 0.0, 0.0, True
        return total_ob, total_m, False

    # ── Overburden map (same API as cosmoaleph_terrain_driver) ───────────────

    def compute_overburden_map(self,
                                det_pos_m: np.ndarray,
                                n_az: int = 36, n_ze: int = 18,
                                ze_max_deg: float = 85.0,
                                step_m: float = 0.5,
                                max_dist_m: float = 10_000.0,
                                progress_cb=None
                                ) -> Tuple[np.ndarray, np.ndarray,
                                           np.ndarray, np.ndarray]:
        """
        Compute overburden [g/cm²] for every (az, ze) bin.

        Parameters
        ──────────
        det_pos_m : detector position [m] in the geometry ENU frame
        n_az      : number of azimuth bins (default 36 = 10° steps)
        n_ze      : number of zenith  bins (default 18 = 5° steps)
        ze_max_deg: maximum zenith angle [°] (default 85°)
        step_m    : march step [m]
        max_dist_m: ray max length [m]
        progress_cb: optional callable(done, total) for progress reporting

        Returns (az_centres, ze_centres, overburden[n_az,n_ze], open_sky[n_az,n_ze])
        """
        az_c  = np.linspace(0, 360, n_az, endpoint=False) + 180/n_az
        ze_c  = np.linspace(0, ze_max_deg, n_ze, endpoint=False) + ze_max_deg/(2*n_ze)
        ob    = np.zeros((n_az, n_ze), dtype=np.float64)
        sky   = np.zeros((n_az, n_ze), dtype=bool)
        total = n_az * n_ze
        done  = 0
        for ia, az in enumerate(az_c):
            for iz, ze in enumerate(ze_c):
                ob[ia, iz], _, sky[ia, iz] = self.ray_overburden(
                    az, ze, det_pos_m, step_m, max_dist_m
                )
                done += 1
            if progress_cb and (ia+1) % max(1, n_az//10) == 0:
                progress_cb(done, total)
        return az_c, ze_c, ob, sky

    # ── Analytical per-cell bounding box ─────────────────────────────────────

    def cell_bbox_analytical(self, cell: "CSGCell") -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the exact bounding box [m] for a CSGCell from its surface tokens,
        without any grid sampling.

        Handles PlaneZ, CylinderZ, and Sphere surface types.
        Exclusion tokens (#cell) are ignored for bbox purposes (they shrink the
        region but never expand it).

        Returns (bbox_min, bbox_max) as np.ndarray([x, y, z]) pairs [m].
        """
        lo = self.bbox_min.copy()   # start from global bounds
        hi = self.bbox_max.copy()

        for sign, sid in cell.tokens:
            surf = self._surf_map.get(sid)
            if surf is None:
                continue

            stype = type(surf).__name__

            if stype == "PlaneZ":
                z0 = surf.z0
                if sign == -1:          # z < z0
                    hi[2] = min(hi[2], z0)
                else:                   # z > z0
                    lo[2] = max(lo[2], z0)

            elif stype == "CylinderZ":
                import math
                r = math.sqrt(surf.r2)
                xc, yc = surf.x0, surf.y0
                if sign == -1:          # inside cylinder → tight xy bbox
                    lo[0] = max(lo[0], xc - r)
                    hi[0] = min(hi[0], xc + r)
                    lo[1] = max(lo[1], yc - r)
                    hi[1] = min(hi[1], yc + r)
                # sign = +1 (outside cylinder) → no constraint on bbox

            elif stype == "Sphere":
                import math
                r = math.sqrt(surf.r2)
                cx, cy, cz = surf.c[0], surf.c[1], surf.c[2]
                if sign == -1:          # inside sphere
                    lo[0] = max(lo[0], cx - r)
                    hi[0] = min(hi[0], cx + r)
                    lo[1] = max(lo[1], cy - r)
                    hi[1] = min(hi[1], cy + r)
                    lo[2] = max(lo[2], cz - r)
                    hi[2] = min(hi[2], cz + r)

        # Clamp to global bbox
        lo = np.maximum(lo, self.bbox_min)
        hi = np.minimum(hi, self.bbox_max)
        return lo, hi

    # ── Plotly 3D preview ─────────────────────────────────────────────────────

    def plotly_preview(self, det_pos_m: np.ndarray = None,
                       highlighted_cells: List[int] = None):
        """
        Return a Plotly 3D figure showing the PHITS/STL geometry.

        STL: renders the actual mesh surface.
        PHITS/CSG: uses `cell_bbox_analytical()` to obtain exact per-cell
        bounding boxes, then draws filled Mesh3d boxes.  This works correctly
        for all cell sizes — from 3 cm detector cylinders to 265 m geology
        layers — without grid sampling.

        Parameters
        ──────────
        det_pos_m          : optional detector position marker [m]
        highlighted_cells  : list of cell_ids to draw with a thick red outline
        """
        import plotly.graph_objects as go
        import math

        if highlighted_cells is None:
            highlighted_cells = []

        fig = go.Figure()

        # ── STL: render the actual mesh ───────────────────────────────────────
        if self._mesh is not None:
            v = self._mesh.vertices
            f = self._mesh.faces
            fig.add_trace(go.Mesh3d(
                x=v[:,0], y=v[:,1], z=v[:,2],
                i=f[:,0], j=f[:,1], k=f[:,2],
                opacity=0.35, color="#4fc3f7", name="STL mesh",
                flatshading=True, lighting=dict(ambient=0.6, diffuse=0.8),
            ))
            if det_pos_m is not None:
                fig.add_trace(go.Scatter3d(
                    x=[det_pos_m[0]], y=[det_pos_m[1]], z=[det_pos_m[2]],
                    mode="markers", name="Detector",
                    marker=dict(size=10, color="red", symbol="diamond"),
                ))
            self._apply_3d_layout(fig)
            return fig

        # ── PHITS/CSG: analytical bboxes + filled Mesh3d boxes ───────────────
        # Colour palette indexed by density (bright = high density, cool = low/void)
        DENSITY_COLOURS = {
            "detector": "#ff1744",      # red — detector (highlighted)
            "void":     "rgba(60,80,120,0.15)",   # near-transparent blue
            "low":      "#29b6f6",      # light blue
            "mid_low":  "#66bb6a",      # green
            "mid":      "#ffca28",      # amber
            "mid_high": "#ffa726",      # orange
            "high":     "#ef5350",      # red-orange (dense rock)
        }
        def _rho_colour(rho, is_det=False, is_void=False):
            if is_det: return "#ff1744"
            if is_void or rho <= 0: return DENSITY_COLOURS["void"]
            if rho < 1.2:  return DENSITY_COLOURS["low"]
            if rho < 1.6:  return DENSITY_COLOURS["mid_low"]
            if rho < 1.8:  return DENSITY_COLOURS["mid"]
            if rho < 2.0:  return DENSITY_COLOURS["mid_high"]
            return DENSITY_COLOURS["high"]

        def _opacity(rho, is_det=False, is_void=False):
            if is_det:  return 0.85
            if is_void: return 0.05
            # Scale opacity 0.25–0.65 by density (denser = more opaque)
            return 0.20 + min(rho / 3.0, 1.0) * 0.45

        def _mesh3d_box(lo, hi, color, opacity, name):
            """Return a filled Mesh3d box from lo to hi."""
            x0, y0, z0 = lo; x1, y1, z1 = hi
            # 8 corners
            xs = [x0,x1,x1,x0, x0,x1,x1,x0]
            ys = [y0,y0,y1,y1, y0,y0,y1,y1]
            zs = [z0,z0,z0,z0, z1,z1,z1,z1]
            # 12 triangles (2 per face × 6 faces)
            ii = [0,0,1,1,2,2,4,4,0,0,1,1]
            jj = [1,3,2,0,3,1,5,7,4,5,5,6]
            kk = [2,2,3,3,0,0,6,6,5,1,6,2]
            return go.Mesh3d(
                x=xs, y=ys, z=zs,
                i=ii, j=jj, k=kk,
                color=color, opacity=opacity,
                name=name, showlegend=True,
                flatshading=True,
                hovertemplate=(
                    f"<b>{name}</b><br>"
                    f"x: {x0:.2f}–{x1:.2f} m<br>"
                    f"y: {y0:.2f}–{y1:.2f} m<br>"
                    f"z: {z0:.2f}–{z1:.2f} m<extra></extra>"
                ),
            )

        def _wireframe_box(lo, hi, color, name, width=3):
            """Wireframe overlay for highlighted cells."""
            x0,y0,z0 = lo; x1,y1,z1 = hi
            e = lambda a,b: [a, b, None]
            xs = e(x0,x1)+e(x0,x1)+e(x0,x1)+e(x0,x1)+e(x0,x0)+e(x1,x1)+e(x0,x0)+e(x1,x1)
            ys = e(y0,y0)+e(y1,y1)+e(y0,y0)+e(y1,y1)+e(y0,y1)+e(y0,y1)+e(y0,y1)+e(y0,y1)
            zs = e(z0,z0)+e(z0,z0)+e(z1,z1)+e(z1,z1)+e(z0,z0)+e(z0,z0)+e(z1,z1)+e(z1,z1)
            return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines",
                                line=dict(color=color, width=width),
                                name=f"{name} [selected]", showlegend=True)

        cells = [c for c in self._cells if c.mat_id not in (-1,)]

        # --- Find the "display range" = the actual subsurface region ---
        # We hide the large outer sphere/cylinder void cells (901, 902, 999)
        # They are air/void and would dominate the view
        SKIP_VOID_GLOBAL = {c.cell_id for c in cells
                            if c.mat_id in (0,) and c.density == 0}

        for cell in cells:
            cmin, cmax = self.cell_bbox_analytical(cell)
            is_det  = cell.cell_id in highlighted_cells
            is_void = (cell.mat_id == 0 or cell.density == 0)

            # Skip cells that are just the outer void (global extent) — too large
            if cell.cell_id in SKIP_VOID_GLOBAL:
                # Still show them, but as faint outline only
                color = "rgba(80,100,150,0.08)"
                label = f"Cell {cell.cell_id} [air/void]"
                # Only show if explicitly highlighted
                if not is_det:
                    continue

            # Expand tiny cells for visibility (min 0.5m in any dimension)
            MIN_VIS = 0.5  # m
            for ax in range(3):
                span = cmax[ax] - cmin[ax]
                if span < MIN_VIS:
                    mid = (cmin[ax] + cmax[ax]) / 2
                    cmin[ax] = mid - MIN_VIS / 2
                    cmax[ax] = mid + MIN_VIS / 2

            color   = _rho_colour(cell.density, is_det, is_void)
            opacity = _opacity(cell.density, is_det, is_void)
            rho_str = f"ρ={cell.density:.2f} g/cm³" if cell.density > 0 else "void/air"
            mat_str = f"mat {cell.mat_id}" if cell.mat_id > 0 else "void"
            label   = f"Cell {cell.cell_id}  [{mat_str}]  {rho_str}"
            if is_det:
                label = f"📡 {label}  ← DETECTOR"

            fig.add_trace(_mesh3d_box(cmin, cmax, color, opacity, label))

            # Thick wireframe outline for highlighted detector cells
            if is_det:
                fig.add_trace(_wireframe_box(cmin, cmax, "#ff1744", label, width=4))

        # ── Detector position marker (if a single-point position is given) ────
        if det_pos_m is not None:
            fig.add_trace(go.Scatter3d(
                x=[det_pos_m[0]], y=[det_pos_m[1]], z=[det_pos_m[2]],
                mode="markers+text",
                marker=dict(size=10, color="red", symbol="diamond",
                            line=dict(color="white", width=1)),
                text=["📡 Detector"], textfont=dict(color="white", size=10),
                textposition="top center",
                name="Detector position", showlegend=True,
            ))

        self._apply_3d_layout(fig)
        return fig

    def _apply_3d_layout(self, fig):
        """Apply standard dark-theme 3D layout."""
        fig.update_layout(
            paper_bgcolor="rgb(15,17,23)",
            font=dict(color="white", size=11),
            height=500,
            margin=dict(l=0, r=0, t=40, b=0),
            scene=dict(
                xaxis_title="East [m]", yaxis_title="North [m]", zaxis_title="Up [m]",
                bgcolor="rgb(20,22,30)",
                xaxis=dict(gridcolor="rgba(255,255,255,0.1)", color="white"),
                yaxis=dict(gridcolor="rgba(255,255,255,0.1)", color="white"),
                zaxis=dict(gridcolor="rgba(255,255,255,0.1)", color="white"),
                camera=dict(eye=dict(x=1.5, y=-1.5, z=0.9)),
            ),
            legend=dict(
                bgcolor="rgba(0,0,0,0.5)", font=dict(color="white", size=10),
                x=1.0, y=0.98, xanchor="right",
            ),
            title=dict(text=f"<b>{self.summary()}</b>",
                       x=0.02, font=dict(color="white", size=12)),
        )


# ─────────────────────────────────────────────────────────────────────────────
# STL loader
# ─────────────────────────────────────────────────────────────────────────────

def load_stl(path: str, density: float = 2.65, units_cm: bool = False) -> CSGGeometry:
    """
    Load an STL file as a watertight triangle mesh.

    Parameters
    ──────────
    path      : path to .stl file (ASCII or binary)
    density   : uniform material density [g/cm³]
    units_cm  : True if the STL file uses cm; False (default) = metres

    Returns CSGGeometry with mesh stored for ray intersection.
    """
    try:
        import trimesh
    except ImportError:
        raise ImportError(
            "trimesh is required to load STL files.\n"
            "Install with:  pip install trimesh"
        )

    mesh = trimesh.load(path, force="mesh")

    if not isinstance(mesh, trimesh.Trimesh):
        # Scene or multiple meshes → merge
        mesh = trimesh.util.concatenate(list(mesh.geometry.values()))

    if units_cm:
        mesh.apply_scale(0.01)   # cm → m

    if not mesh.is_watertight:
        # Attempt to heal
        trimesh.repair.fix_normals(mesh)
        trimesh.repair.fill_holes(mesh)

    geom = CSGGeometry()
    geom.source_format = "STL"
    geom.source_path   = str(path)
    geom.label         = Path(path).stem
    geom._mesh         = mesh
    geom._stl_density  = float(density)

    bb = mesh.bounds           # [[xmin,ymin,zmin],[xmax,ymax,zmax]] in metres
    geom.bbox_min = bb[0] - 1.0
    geom.bbox_max = bb[1] + 1.0

    print(f"  STL loaded: {Path(path).name}  "
          f"{len(mesh.faces)} faces  "
          f"bounds=[{bb[0]}]–[{bb[1]}]  "
          f"watertight={mesh.is_watertight}  "
          f"ρ={density} g/cm³", flush=True)
    return geom


# ─────────────────────────────────────────────────────────────────────────────
# PHITS / MCNP surface card parser
# ─────────────────────────────────────────────────────────────────────────────

_FLOAT = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

def _tok_floats(s: str) -> List[float]:
    """Extract all floats from a string."""
    return [float(x) for x in re.findall(_FLOAT, s)]


def _eval_phits_param_expr(s: str, c_vars: Dict[int, float]) -> str:
    """
    Evaluate PHITS set:c[N] variable references in a surface parameter string.

    Handles all patterns found in real PHITS input files:
      NUM+cN   e.g. -3900.0+c5  → -3900.0+12000.0 → 8100.0
      NUM-cN   e.g. -3900.0-c5  → -3900.0-12000.0 → -15900.0
      cN+NUM   e.g. c5-100.0    → 12000.0-100.0   → 11900.0
      bare cN  e.g. so c1       → so 16000.0

    Variables are substituted in descending order of N so that
    c10 is replaced before c1, avoiding partial matches.

    Safe to call on parameter-only strings (after the surface type token
    has been stripped), so c/z surface types are never affected.
    """
    if not c_vars:
        return s

    for n in sorted(c_vars.keys(), reverse=True):
        val  = c_vars[n]
        cn   = f'c{n}'
        esc  = re.escape(cn)
        guard = r'(?!\d)'   # not followed by digit (avoids c1 matching c12)

        # Pattern 1: float_expr  +/-  cN
        def _repl1(m: re.Match, v: float = val) -> str:
            base = float(m.group(1))
            return f'{base + v:.10g}' if m.group(2) == '+' else f'{base - v:.10g}'
        s = re.sub(
            r'([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*([+-])\s*' + esc + guard,
            _repl1, s)

        # Pattern 2: cN  +/-  float_expr
        def _repl2(m: re.Match, v: float = val) -> str:
            num = float(m.group(2))
            return f'{v + num:.10g}' if m.group(1) == '+' else f'{v - num:.10g}'
        s = re.sub(
            esc + guard + r'\s*([+-])\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)',
            _repl2, s)

        # Pattern 3: bare cN (not preceded by alphanumeric / slash)
        s = re.sub(r'(?<![/\w])' + esc + guard, f'{val:.10g}', s)

    return s


def _parse_phits_set_vars(text: str) -> Dict[int, float]:
    """
    Scan a PHITS input file text and collect all set:c[N][value] variable
    definitions.  Returns {N: value} dict.

    PHITS syntax: set:c1[16000]   or   set:c5[12000.0]
    These may appear in [Parameters], [Surface], or other sections.
    """
    c_vars: Dict[int, float] = {}
    for line in text.splitlines():
        clean = line.strip()
        if clean.startswith('$') or clean.startswith('!'):
            continue
        # Strip inline comments
        clean = re.split(r'\$', clean)[0].strip()
        m = re.match(r'set\s*:\s*c(\d+)\s*\[\s*(' + _FLOAT + r')\s*\]',
                     clean, re.IGNORECASE)
        if m:
            n   = int(m.group(1))
            val = float(m.group(2))
            c_vars[n] = val
    return c_vars


def _parse_surface_card(sid: int, stype: str, vals: List[float]) -> Optional[_SurfaceBase]:
    """
    Parse one surface card and return a _SurfaceBase subclass.

    Supports: PX PY PZ P SO S SX SY SZ
              CX CY CZ  C/X C/Y C/Z
              KX KY KZ
              RPP RCC SPH TRC GQ
    All input values in metres (caller converts from cm if needed).
    """
    t = stype.upper()

    # ── Planes ──────────────────────────────────────────────────────────────
    if t == "PX":
        return PlaneX(sid, vals[0])
    if t == "PY":
        return PlaneY(sid, vals[0])
    if t == "PZ":
        return PlaneZ(sid, vals[0])
    if t == "P":
        a,b,c,d = vals[0],vals[1],vals[2],vals[3]
        return Plane(sid, a,b,c,d)

    # ── Spheres ─────────────────────────────────────────────────────────────
    if t in ("SO",):
        return Sphere(sid, 0,0,0, vals[0])
    if t == "S":
        return Sphere(sid, vals[0],vals[1],vals[2], vals[3])
    if t == "SX":
        return Sphere(sid, vals[0],0,0, vals[1])   # PHITS: SX x r
    if t == "SY":
        return Sphere(sid, 0,vals[0],0, vals[1])
    if t == "SZ":
        return Sphere(sid, 0,0,vals[0], vals[1])

    # ── Cylinders ────────────────────────────────────────────────────────────
    if t == "CX":
        return CylinderX(sid, 0,0, vals[0])
    if t == "CY":
        return CylinderY(sid, 0,0, vals[0])
    if t == "CZ":
        return CylinderZ(sid, 0,0, vals[0])
    if t in ("C/X",):
        return CylinderX(sid, vals[0],vals[1], vals[2])
    if t in ("C/Y",):
        return CylinderY(sid, vals[0],vals[1], vals[2])
    if t in ("C/Z",):
        return CylinderZ(sid, vals[0],vals[1], vals[2])

    # ── Cones ────────────────────────────────────────────────────────────────
    if t == "KX":
        # KX x t^2 [sheet]  — (y^2+z^2) - t^2*(x-x0)^2 = 0
        x0=vals[0]; t2=vals[1]
        class _CX(_SurfaceBase):
            def __init__(self,sid,x0,t2): self.surf_id=sid; self.x0=x0; self.t2=t2
            def f(self,p): return (p[1]**2+p[2]**2) - self.t2*(p[0]-self.x0)**2
        return _CX(sid, x0, t2)
    if t == "KY":
        y0=vals[0]; t2=vals[1]
        class _CY(_SurfaceBase):
            def __init__(self,sid,y0,t2): self.surf_id=sid; self.y0=y0; self.t2=t2
            def f(self,p): return (p[0]**2+p[2]**2) - self.t2*(p[1]-self.y0)**2
        return _CY(sid, y0, t2)
    if t == "KZ":
        z0=vals[0]; t2=vals[1]
        return ConeZ(sid, 0,0, z0, t2)

    # ── Macrobodies ──────────────────────────────────────────────────────────
    if t == "RPP":
        # RPP xmin xmax ymin ymax zmin zmax → 6 bounding planes as a compound
        # We expand into 6 planes and store as a pseudo-surface using Quadric
        xmin,xmax,ymin,ymax,zmin,zmax = vals[:6]
        # Use a helper compound class
        class _RPP(_SurfaceBase):
            def __init__(self, sid, vals):
                self.surf_id=sid
                self.xmin,self.xmax=vals[0],vals[1]
                self.ymin,self.ymax=vals[2],vals[3]
                self.zmin,self.zmax=vals[4],vals[5]
            def f(self, p):
                # Inside RPP ↔ all six inequalities hold → return negative value inside
                dx = max(self.xmin-p[0], p[0]-self.xmax)
                dy = max(self.ymin-p[1], p[1]-self.ymax)
                dz = max(self.zmin-p[2], p[2]-self.zmax)
                # Convention: f<0 inside, f>0 outside
                # Inside box: all d<0 → f = max(dx,dy,dz) < 0
                return max(dx, dy, dz)
        return _RPP(sid, vals[:6])

    if t == "SPH":
        # SPH x y z r  (macrobody sphere)
        return Sphere(sid, vals[0],vals[1],vals[2], vals[3])

    if t == "RCC":
        # RCC vx vy vz  hx hy hz  r
        # Right circular cylinder: vxyz = base centre, hxyz = axis vector, r = radius
        vx,vy,vz,hx,hy,hz,r = vals[:7]
        class _RCC(_SurfaceBase):
            def __init__(self,sid,v,h,r):
                self.surf_id=sid
                self.v=np.array(v,dtype=float)
                self.h=np.array(h,dtype=float)
                self.h_hat=self.h/np.linalg.norm(self.h)
                self.h_len=np.linalg.norm(self.h)
                self.r2=r*r
            def f(self, p):
                d  = np.array(p,dtype=float) - self.v
                s  = float(d @ self.h_hat)
                if s < 0 or s > self.h_len: return 1.0   # outside caps
                perp2 = float(d@d) - s*s
                return perp2 - self.r2
        return _RCC(sid, [vx,vy,vz],[hx,hy,hz], r)

    if t == "TRC":
        # TRC vx vy vz  hx hy hz  r_base r_top
        vx,vy,vz,hx,hy,hz,rb,rt = vals[:8]
        class _TRC(_SurfaceBase):
            def __init__(self,sid,v,h,rb,rt):
                self.surf_id=sid
                self.v=np.array(v,dtype=float)
                self.h=np.array(h,dtype=float)
                self.h_hat=self.h/np.linalg.norm(self.h)
                self.h_len=np.linalg.norm(self.h)
                self.rb=rb; self.rt=rt
            def f(self, p):
                d = np.array(p,dtype=float) - self.v
                s = float(d @ self.h_hat)
                if s < 0 or s > self.h_len: return 1.0
                r_at_s = self.rb + (self.rt - self.rb) * s / self.h_len
                perp2  = float(d@d) - s*s
                return perp2 - r_at_s**2
        return _TRC(sid, [vx,vy,vz],[hx,hy,hz], rb, rt)

    if t == "GQ":
        # GQ A B C D E F G H I J  (10 coefficients)
        return Quadric(sid, vals[:10])

    return None   # unsupported surface type


def _parse_cell_tokens(expr: str) -> List:
    """
    Parse a PHITS/MCNP geometry expression string into tokens.

    Input:  "-1 2 -3 : -4 5"  (space=AND, :=OR, #=complement not yet supported)
            "(-1 2) : (-3 4)"  (parentheses are stripped)
    Output: [(-1,1),(+1,2),(-1,3), 'OR', (-1,4),(+1,5)]
    Each non-OR token is (sign, surf_id_abs).

    Cell ordering note: in PHITS/MCNP, cells are evaluated in the order they
    appear in the input deck.  More specific regions (e.g. an anomaly sphere
    inside the bulk rock) must be listed BEFORE the enclosing region, or the
    enclosing cell matches first and the anomaly is never reached.
    """
    tokens = []
    # Remove parentheses (used for visual grouping only in PHITS/MCNP)
    expr = expr.replace('(', ' ').replace(')', ' ')
    # Replace : and | with a canonical OR separator
    expr = re.sub(r'[:|]', ' OR ', expr)
    for tok in expr.split():
        if tok.upper() == 'OR':
            if tokens and tokens[-1] != 'OR':
                tokens.append('OR')
            continue
        try:
            n = int(tok)
            sign = -1 if n < 0 else +1
            tokens.append((sign, abs(n)))
        except ValueError:
            pass   # skip unknown tokens (like region keywords)
    return tokens


# ─────────────────────────────────────────────────────────────────────────────
# PHITS input file parser
# ─────────────────────────────────────────────────────────────────────────────

def _strip_comments(line: str) -> str:
    """Remove PHITS inline comments ($ ...) and return stripped line."""
    idx = line.find('$')
    if idx >= 0:
        line = line[:idx]
    return line.strip()


def load_phits(path: str) -> CSGGeometry:
    """
    Parse a PHITS input file and extract surface + cell definitions.

    Supported sections: [Surface], [Cell]
    All dimensions converted from cm (PHITS default) to metres.

    Variable substitution: PHITS set:c[N][value] variables (e.g. set:c5[12000])
    are collected in a first pass and substituted into every surface parameter
    expression (e.g. ``pz -3900.0+c5`` → ``pz 8100.0``) before float parsing.
    Without this step the geometry is misplaced by the variable's value, causing
    zero hits even when muons are correctly targeted.
    """
    text = Path(path).read_text(errors="replace")
    geom = CSGGeometry()
    geom.source_format = "PHITS"
    geom.source_path   = str(path)
    geom.label         = Path(path).stem
    geom._cm_to_m      = True

    # ── Pass 0: collect set:c[N][value] variables ─────────────────────────────
    c_vars = _parse_phits_set_vars(text)
    if c_vars:
        _cv_str = ', '.join(f'c{k}={v:g}' for k, v in sorted(c_vars.items()))
        print(f"  PHITS variables: {_cv_str}", flush=True)

    surf_raw: Dict[int, Tuple[str, List[float]]] = {}   # sid → (type, vals_m)
    cell_raw: List[Dict] = []

    # ── Pass 1: section-by-section parse ─────────────────────────────────────
    lines = text.splitlines()
    section = None

    # continuation: PHITS uses & at line end
    merged: List[str] = []
    for line in lines:
        clean = _strip_comments(line)
        if clean.endswith('&'):
            merged.append(clean[:-1].rstrip())
            continue
        if merged:
            clean = ' '.join(merged) + ' ' + clean
            merged.clear()
        merged_line = clean

        sec_match = re.match(r'^\s*\[(.+?)\]', merged_line, re.IGNORECASE)
        if sec_match:
            section = sec_match.group(1).strip().lower()
            continue

        if not merged_line or merged_line.startswith('c ') or merged_line.lower().startswith('c\t'):
            continue

        # ── [Surface] section ─────────────────────────────────────────────
        if section in ('surface', 'surfaces'):
            # Format: surf_id  [transform]  surf_type  params...
            parts = merged_line.split()
            if len(parts) < 2: continue
            try:
                sid = int(parts[0])
            except ValueError:
                continue
            # Skip optional transform ID (integer after sid before alpha token)
            idx = 1
            if idx < len(parts):
                try:
                    int(parts[idx])
                    idx += 1   # skip transform number
                except ValueError:
                    pass
            if idx >= len(parts): continue
            stype = parts[idx]
            # ── Apply variable substitution on parameter string only ───────
            # stype is always alphabetic (e.g. "pz", "c/z") — never matches cN
            param_raw = ' '.join(parts[idx+1:])
            param_sub = _eval_phits_param_expr(param_raw, c_vars)
            vals_cm   = _tok_floats(param_sub)
            vals_m    = [v * 0.01 for v in vals_cm]   # cm → m
            surf_raw[sid] = (stype, vals_m)

        # ── [Cell] section ────────────────────────────────────────────────
        elif section in ('cell', 'cells'):
            parts = merged_line.split()
            if len(parts) < 3: continue
            try:
                cid = int(parts[0])
            except ValueError:
                continue
            try:
                mat_id = int(parts[1])
            except ValueError:
                continue
            # Density (negative = g/cm³ for mass, positive = at/cm³ for number)
            # Void cell: mat_id = -1 or 0 with no density
            rho = 0.0
            geom_start = 2
            if mat_id not in (-1, 0):
                try:
                    rho = abs(float(parts[2]))
                    geom_start = 3
                except ValueError:
                    pass
            geom_expr = ' '.join(parts[geom_start:])
            # Strip any trailing keywords (imp:n=1 etc.)
            geom_expr = re.sub(r'\s+imp:\w+=\S+', '', geom_expr)
            geom_expr = re.sub(r'\s+imp=\S+', '', geom_expr)
            tokens = _parse_cell_tokens(geom_expr)
            cell_raw.append(dict(cid=cid, mat=mat_id, rho=rho, tokens=tokens))

    # ── Build surface map ─────────────────────────────────────────────────────
    for sid, (stype, vals_m) in surf_raw.items():
        s = _parse_surface_card(sid, stype, vals_m)
        if s is not None:
            geom._surf_map[sid] = s

    # ── Build cell list ───────────────────────────────────────────────────────
    for cd in cell_raw:
        cell = CSGCell(
            cell_id = cd['cid'],
            mat_id  = cd['mat'],
            density = cd['rho'],
            tokens  = cd['tokens'],
            label   = f"Cell {cd['cid']} mat={cd['mat']} ρ={cd['rho']:.2f}",
        )
        geom._cells.append(cell)

    # ── Estimate bounding box from surface extents ────────────────────────────
    _estimate_bbox(geom)

    mats = geom.materials()
    print(f"  PHITS loaded: {Path(path).name}  "
          f"{len(geom._surf_map)} surfaces  "
          f"{len(geom._cells)} cells  "
          f"materials={mats}", flush=True)
    return geom


# ─────────────────────────────────────────────────────────────────────────────
# MCNP input file parser
# ─────────────────────────────────────────────────────────────────────────────

def load_mcnp(path: str) -> CSGGeometry:
    """
    Parse an MCNP5/6 input file and extract surface + cell cards.

    MCNP and PHITS surface card syntax is nearly identical.
    Key differences handled here:
      • No [Section] headers — blocks separated by blank lines
        (cell block first, then surface block, then data block)
      • Comments: lines starting with 'c ' or 'C '
      • Continuation lines start with 5+ spaces
    All dimensions converted from cm to metres.
    """
    text = Path(path).read_text(errors="replace")
    geom = CSGGeometry()
    geom.source_format = "MCNP"
    geom.source_path   = str(path)
    geom.label         = Path(path).stem
    geom._cm_to_m      = True

    # MCNP input: title card (line 1), cell block, blank, surface block, blank, data block
    lines = text.splitlines()

    # Merge continuation lines (5+ leading spaces)
    merged_lines: List[str] = []
    for i, line in enumerate(lines):
        if i == 0:
            merged_lines.append('')   # skip title card
            continue
        # Comment line — bare 'c' acts as blank-line separator; inline keeps block context
        stripped_lo = line.strip().lower()
        if stripped_lo == 'c' or stripped_lo == '':
            merged_lines.append('')   # treat as blank separator
            continue
        if re.match(r'^[cC][\s$]', line):
            continue   # inline comment — skip but don't add blank
        # Inline comment
        line = re.sub(r'\s*\$.*', '', line)
        # Continuation
        if merged_lines and re.match(r'^ {5,}', line):
            merged_lines[-1] = merged_lines[-1].rstrip() + ' ' + line.strip()
        else:
            merged_lines.append(line.rstrip())

    # Split into blocks by blank lines
    blocks: List[List[str]] = []
    current: List[str] = []
    for line in merged_lines:
        if line.strip() == '':
            if current:
                blocks.append(current)
                current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)

    # block[0] = cells, block[1] = surfaces, block[2+] = data
    # (title card becomes a blank line → no block created, so index starts at 0)
    cell_lines  = blocks[0] if len(blocks) > 0 else []
    surf_lines  = blocks[1] if len(blocks) > 1 else []

    # ── Parse surface block ───────────────────────────────────────────────────
    surf_raw: Dict[int, Tuple[str, List[float]]] = {}
    for line in surf_lines:
        parts = line.split()
        if not parts: continue
        try:
            sid = int(parts[0])
        except ValueError:
            continue
        idx = 1
        # Optional transform integer
        if idx < len(parts):
            try:
                int(parts[idx]); idx += 1
            except ValueError:
                pass
        if idx >= len(parts): continue
        stype = parts[idx]
        vals_cm = _tok_floats(' '.join(parts[idx+1:]))
        vals_m  = [v * 0.01 for v in vals_cm]
        surf_raw[sid] = (stype, vals_m)

    # ── Parse cell block ──────────────────────────────────────────────────────
    cell_raw: List[Dict] = []
    for line in cell_lines:
        parts = line.split()
        if len(parts) < 3: continue
        try:
            cid = int(parts[0])
        except ValueError:
            continue
        try:
            mat_id = int(parts[1])
        except ValueError:
            continue
        rho = 0.0
        geom_start = 2
        if mat_id not in (0, -1):
            try:
                rho = abs(float(parts[2])); geom_start = 3
            except ValueError:
                pass
        geom_expr = ' '.join(parts[geom_start:])
        # Strip MCNP keywords: imp:n=, imp:p=, vol=, etc.
        geom_expr = re.sub(r'\s+\w+[:=]\S+', '', geom_expr)
        tokens = _parse_cell_tokens(geom_expr)
        cell_raw.append(dict(cid=cid, mat=mat_id, rho=rho, tokens=tokens))

    for sid, (stype, vals_m) in surf_raw.items():
        s = _parse_surface_card(sid, stype, vals_m)
        if s is not None:
            geom._surf_map[sid] = s

    for cd in cell_raw:
        geom._cells.append(CSGCell(
            cell_id=cd['cid'], mat_id=cd['mat'], density=cd['rho'],
            tokens=cd['tokens'],
            label=f"Cell {cd['cid']} mat={cd['mat']} ρ={cd['rho']:.2f}",
        ))

    _estimate_bbox(geom)

    mats = geom.materials()
    print(f"  MCNP loaded: {Path(path).name}  "
          f"{len(geom._surf_map)} surfaces  "
          f"{len(geom._cells)} cells  "
          f"materials={mats}", flush=True)
    return geom


# ─────────────────────────────────────────────────────────────────────────────
# Bounding box estimator
# ─────────────────────────────────────────────────────────────────────────────

def _estimate_bbox(geom: CSGGeometry, margin_m: float = 50.0):
    """
    Estimate geometry bounding box from PlaneX/Y/Z surfaces and Sphere radii.
    Falls back to ±1 km if no planes found.
    """
    extremes = []
    for s in geom._surf_map.values():
        if isinstance(s, PlaneX): extremes += [s.x0]
        elif isinstance(s, PlaneY): extremes += [s.y0]
        elif isinstance(s, PlaneZ): extremes += [s.z0]
        elif isinstance(s, Sphere):
            r = np.sqrt(s.r2)
            for dim in range(3):
                extremes += [s.c[dim]-r, s.c[dim]+r]
        elif isinstance(s, CylinderZ):
            r = np.sqrt(s.r2)
            extremes += [s.x0-r, s.x0+r, s.y0-r, s.y0+r]

    if extremes:
        lo = min(extremes) - margin_m
        hi = max(extremes) + margin_m
        geom.bbox_min = np.array([lo, lo, lo])
        geom.bbox_max = np.array([hi, hi, hi])
    else:
        geom.bbox_min = np.array([-1000.0, -1000.0, -1000.0])
        geom.bbox_max = np.array([ 1000.0,  1000.0,  1000.0])


# ─────────────────────────────────────────────────────────────────────────────
# Streamlit UI
# ─────────────────────────────────────────────────────────────────────────────

def render_csg_builder() -> Optional[CSGGeometry]:
    """
    Render the PHITS / STL geometry loader inside the current Streamlit container.

    Returns a CSGGeometry when the user has successfully loaded a file, or None.
    MCNP support has been removed — use PHITS input for CSG geometries.
    """
    st.caption(
        "Upload your PHITS input file (or an STL mesh).  "
        "The geometry is parsed to extract cells, surfaces, and material densities.  "
        "After loading, select the detector cell in **Section 3b** below."
    )

    # ── Format selector (PHITS or STL only) ──────────────────────────────────
    fmt = st.radio(
        "Format",
        ["🔩 PHITS input (.inp)", "🧊 STL mesh (.stl)"],
        horizontal=True, key="csg_format",
    )
    is_phits = fmt.startswith("🔩")

    # ── File upload ───────────────────────────────────────────────────────────
    accept_label = "PHITS input file (.inp, .i, .phits, .txt)" if is_phits else "STL mesh file (.stl)"
    uploaded = st.file_uploader(accept_label, type=None, key="csg_file_upload",
                                help="PHITS [Surface] / [Cell] deck or binary/ASCII STL.")

    # ── STL-specific density setting ──────────────────────────────────────────
    if not is_phits:
        _sc1, _sc2 = st.columns(2)
        stl_rho = _sc1.number_input("Material density ρ [g/cm³]", 0.1, 20.0, 2.65, 0.05,
                                     key="csg_stl_rho",
                                     help="Uniform density applied to the STL solid.")
        _sc2.checkbox("File units in cm (not m)", value=False, key="csg_stl_cm",
                      help="Tick if the STL was exported in centimetres.")

    # ── Load / Clear buttons ──────────────────────────────────────────────────
    st.divider()
    load_col, clear_col = st.columns([3, 1])
    with clear_col:
        if st.button("🗑 Clear", key="csg_clear", width='stretch'):
            for k in ("_csg_geom", "terrain_dem_path"):
                st.session_state.pop(k, None)
            st.rerun()
    with load_col:
        do_load = st.button("🔩 Load Geometry", type="primary",
                            key="csg_load", width='stretch')

    if do_load:
        if uploaded is None:
            st.error("Please upload a geometry file first.")
        else:
            import tempfile, os
            suffix = Path(uploaded.name).suffix or ".inp"
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(uploaded.read())
                tmp_path = tmp.name
            try:
                with st.spinner("Parsing geometry…"):
                    if is_phits:
                        geom = load_phits(tmp_path)
                    else:
                        geom = load_stl(
                            tmp_path,
                            density  = float(st.session_state.get("csg_stl_rho", 2.65)),
                            units_cm = bool(st.session_state.get("csg_stl_cm", False)),
                        )
                    geom.source_path = uploaded.name
                st.session_state["_csg_geom"] = geom
            except Exception as exc:
                st.error(f"❌ Load failed: {exc}")
                import traceback; st.code(traceback.format_exc())
            finally:
                try: os.unlink(tmp_path)
                except Exception: pass

    # ── Preview when geometry is loaded ──────────────────────────────────────
    geom: Optional[CSGGeometry] = st.session_state.get("_csg_geom")

    if geom is None:
        return None

    # ── Success banner ────────────────────────────────────────────────────────
    st.success(f"✅  {geom.summary()}", icon="🔩")

    # ── 3D geometry viewer ────────────────────────────────────────────────────
    st.markdown("#### 🗺️ Geometry Preview")
    st.caption(
        "Each coloured wireframe box shows the approximate bounding region of one cell.  "
        "The detector cell you select in Section 3b will be highlighted there."
    )
    with st.spinner("Rendering geometry…"):
        try:
            # Get any detector cells already selected (from Section 3b)
            _hl = st.session_state.get("_csg_selected_det_cells", [])
            fig = geom.plotly_preview(highlighted_cells=_hl)
            st.plotly_chart(fig)
        except Exception as exc:
            st.warning(f"3D preview unavailable: {exc}")

    # ── Materials table ───────────────────────────────────────────────────────
    mats = geom.materials()
    cells_nv = [c for c in geom._cells if c.density > 0]
    cells_v  = [c for c in geom._cells if c.density == 0 and c.mat_id != -1]

    if mats or cells_v:
        st.markdown("#### 🧱 Geometry Cells & Materials")

        # Build rows for ALL non-outside cells
        all_cells = [c for c in geom._cells if c.mat_id != -1]
        if all_cells:
            # Colour scale: map density → a colour band
            max_rho = max((c.density for c in all_cells if c.density > 0), default=1.0)

            def _density_bar(rho: float, max_rho: float) -> str:
                """Return an inline CSS colour based on density value."""
                if rho <= 0:
                    return "#2e3a4e"          # void / air — dark blue-grey
                frac = min(rho / max_rho, 1.0)
                # Interpolate steel-blue → amber → red
                if frac < 0.5:
                    r = int(79  + frac * 2 * (255 - 79))
                    g = int(195 - frac * 2 * (195 - 160))
                    b = int(247 - frac * 2 * 247)
                else:
                    r = 255
                    g = int(160 - (frac - 0.5) * 2 * 160)
                    b = 0
                return f"rgb({r},{g},{b})"

            rows_html = ""
            for c in sorted(all_cells, key=lambda x: x.cell_id):
                bg   = _density_bar(c.density, max_rho)
                void = c.density == 0
                rho_str = "— void / air" if void else f"{c.density:.3f} g/cm³"
                mat_str = f"mat {c.mat_id}" if c.mat_id > 0 else "void"
                label   = c.label if c.label else ""
                text_c  = "#dde" if void else "#fff"
                rows_html += (
                    f"<tr>"
                    f"<td style='padding:5px 10px;font-weight:600;color:#fff;'>{c.cell_id}</td>"
                    f"<td style='padding:5px 10px;color:#aac;'>{mat_str}</td>"
                    f"<td style='padding:5px 10px;'>"
                    f"  <span style='display:inline-block;width:14px;height:14px;"
                    f"  border-radius:3px;background:{bg};vertical-align:middle;"
                    f"  margin-right:6px;border:1px solid rgba(255,255,255,0.2);'></span>"
                    f"  <span style='color:{text_c};font-family:monospace;'>{rho_str}</span>"
                    f"</td>"
                    f"<td style='padding:5px 10px;color:#8a9;font-size:0.85em;'>{label}</td>"
                    f"</tr>"
                )

            table_html = f"""
<div style="background:rgb(22,26,36);border-radius:8px;padding:10px;
            border:1px solid rgba(255,255,255,0.08);overflow-x:auto;">
  <table style="width:100%;border-collapse:collapse;font-size:0.9em;">
    <thead>
      <tr style="border-bottom:1px solid rgba(255,255,255,0.15);">
        <th style="padding:6px 10px;text-align:left;color:#7a8aaa;font-weight:500;">Cell ID</th>
        <th style="padding:6px 10px;text-align:left;color:#7a8aaa;font-weight:500;">Material</th>
        <th style="padding:6px 10px;text-align:left;color:#7a8aaa;font-weight:500;">Density</th>
        <th style="padding:6px 10px;text-align:left;color:#7a8aaa;font-weight:500;">Label</th>
      </tr>
    </thead>
    <tbody>{rows_html}</tbody>
  </table>
</div>"""
            st.markdown(table_html, unsafe_allow_html=True)

    # ── Quick single-ray overburden test ──────────────────────────────────────
    with st.expander("🔬  Single-ray overburden check", expanded=False):
        _t1, _t2 = st.columns(2)
        test_az = _t1.number_input("Azimuth [°]", 0.0, 360.0, 0.0, 5.0, key="csg_test_az")
        test_ze = _t2.number_input("Zenith [°]",  0.0,  90.0, 0.0, 5.0, key="csg_test_ze")
        det_pos = st.session_state.get("_csg_det_pos", np.zeros(3))
        if st.button("▶ Trace ray", key="csg_test_ray"):
            ob, slant, sky = geom.ray_overburden(
                test_az, test_ze, det_pos,
                step_m    = float(st.session_state.get("csg_step", 0.5)),
                max_dist_m= float(st.session_state.get("csg_maxdist", 5000.0)),
            )
            if sky:
                st.info(f"Az={test_az:.0f}°  Ze={test_ze:.0f}°  →  open geometry (no rock)")
            else:
                st.success(
                    f"Az={test_az:.0f}°  Ze={test_ze:.0f}°  →  "
                    f"slant={slant:.1f} m  |  OB={ob:.0f} g/cm²  |  ≈{ob/263:.0f} m w.e."
                )

    return geom


def _render_ob_heatmap(az_c, ze_c, ob_map, sky_map):
    """Render a simple heatmap of the overburden map for quick QC."""
    import plotly.graph_objects as go
    disp = np.where(sky_map, 0.0, ob_map)
    fig  = go.Figure(go.Heatmap(
        z=disp, x=ze_c, y=az_c,
        colorscale="Viridis",
        colorbar=dict(title="OB [g/cm²]", thickness=12),
        hovertemplate="Az=%{y:.0f}°  Ze=%{x:.0f}°  OB=%{z:.0f} g/cm²<extra></extra>",
    ))
    fig.update_layout(
        **DARK, height=320,
        xaxis_title="Zenith [°]", yaxis_title="Azimuth [°]",
        title=dict(text="Overburden map preview [g/cm²]", font=dict(color="white",size=12)),
        margin=dict(l=60, r=20, t=40, b=50),
    )
    st.plotly_chart(fig)
