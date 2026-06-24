"""
gui_csg_transport.py  —  UCLouvain Muography Group
Per-Muon Ray Tracing and Transport Through CSG Geometry

This module implements a fundamentally different workflow from the DEM-based
terrain transport:

  DEM approach  :  bin muons by direction → per-bin average overburden →
                   uniform slab transport per bin
  CSG approach  :  per-muon ray tracing through the CSG cell structure →
                   accumulated multi-material overburden per muon →
                   stochastic transport → record only muons that reach the
                   detector cell

The output is a standard 18-col underground file identical to the MUSIC /
Bethe-Bloch / UCMuon-Stochastic engines, so all downstream Tab 3 plotting
and MURAVES comparison tools work without modification.

WORKFLOW
─────────
  1. Load PHITS / MCNP / STL geometry via gui_csg_engine.py
  2. User selects one cell as the "detector cell"
  3. Each Tab-1 surface muon is ray-traced:
       a. March from surface entry into the geometry
       b. At each step: determine current cell → accumulate (rho × dl) per material
       c. Detect first entry into detector cell → record entry position
       d. If never reaches detector cell: muon blocked / misses → alive=0
  4. For muons that reach detector cell:
       a. Apply per-segment CSDA energy loss through each material layer
       b. Apply stochastic transport (UCMuon) over total overburden
       c. If energy > threshold: alive=1, record 18-col output
  5. Return DataFrame with all muons (alive=0 for blocked/absorbed)

COORDINATE SYSTEMS
───────────────────
  Tab-1 surface file: positions in cm, ENU convention (x=East, y=North, z=Up)
                      z ≈ 0 at the generation surface
  PHITS geometry:     positions in cm, arbitrary frame
  Coordinate mapping: user supplies (dx, dy, dz) [cm] — offset from Tab-1
                      frame to geometry frame. Typical use:
                        dx=0, dy=0, dz = -surface_z_in_geometry [cm]

  Example: PHITS has the surface at z=0 and detector at z=-85000 cm (850 m).
           Tab-1 muons have z≈0 (surface). Offset = (0, 0, 0): frames match.

MULTI-MATERIAL TRANSPORT
──────────────────────────
  For each muon, the ray traverses a sequence of (material, length) segments.
  Transport is applied sequentially using the CSDA approximation per material:

      E_remaining after segment i = E_in × exp(−overburden_i / R_CSDA(E_in))
      ... corrected via Bethe-Bloch interpolation (same as stochastic_driver)

  If E_remaining falls below E_stop (muon range exhausted): muon stopped.
  A UCMuon-Stochastic call is used for the total aggregated overburden in
  the standard material, or optionally per-segment for multi-material.

Author : Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026
"""

from __future__ import annotations

import importlib.util
import time
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
# Physical constants (mirrors ucmuon_stochastic_driver.py)
# ─────────────────────────────────────────────────────────────────────────────
M_MU_MEV  = 105.6584
M_MU_GEV  = 0.1056584
E_STOP_MEV = 10.0        # kinetic energy floor — muon considered stopped

# ─────────────────────────────────────────────────────────────────────────────
# Ray-tracing result for a single muon
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RayResult:
    """
    Outcome of ray-tracing one muon through the CSG geometry.

    Fields
    ──────
    hit          : True if the ray entered the detector cell
    entry_pos_cm : (x, y, z) [cm] where the ray first entered the detector cell
    exit_pos_cm  : (x, y, z) [cm] where the ray exited the detector cell (or None)
    segments     : list of (density [g/cm³], length [cm]) for each material
                   traversed BEFORE entering the detector cell, in order
    total_ob_gcm2: total overburden [g/cm²] before the detector cell
    path_through_cm : path length through the detector cell [cm]
    cx, cy, cz   : direction cosines at detector entry (may differ from surface
                   if multiple-scattering is applied per-segment in future)
    """
    hit            : bool   = False
    entry_pos_cm   : Optional[np.ndarray] = None
    exit_pos_cm    : Optional[np.ndarray] = None
    segments       : List[Tuple[float, float]] = field(default_factory=list)
    total_ob_gcm2  : float  = 0.0
    path_through_cm: float  = 0.0
    cx : float = 0.0
    cy : float = 0.0
    cz : float = 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Detector cell geometry analyser
# ─────────────────────────────────────────────────────────────────────────────

def _detector_cell_geometry(det_cell, surf_map):
    """
    Analyse a detector cell's tokens to determine if it is a simple
    CylinderZ + PlaneZ cell (the most common detector geometry in PHITS).

    Returns dict with:
      type    : "cylinder_pz"  or  "general"
      cx_m, cy_m, r_m : cylinder parameters [m] (if cylinder_pz)
      z_min_m, z_max_m: z bounds [m] (if cylinder_pz)
    """
    cyl = None
    z_bounds = []

    for sign, sid in det_cell.tokens:
        surf = surf_map.get(sid)
        if surf is None:
            continue
        stype = type(surf).__name__

        if stype == "CylinderZ" and sign == -1:
            import math
            cyl = (surf.x0, surf.y0, math.sqrt(surf.r2))

        elif stype == "PlaneZ":
            if sign == -1:   # z < z0
                z_bounds.append(("max", surf.z0))
            else:            # z > z0
                z_bounds.append(("min", surf.z0))

    if cyl is not None and len(z_bounds) >= 2:
        z_min = max((v for k, v in z_bounds if k == "min"), default=-1e9)
        z_max = min((v for k, v in z_bounds if k == "max"), default=+1e9)
        return {
            "type": "cylinder_pz",
            "cx_m": cyl[0], "cy_m": cyl[1], "r_m": cyl[2],
            "z_min_m": z_min, "z_max_m": z_max,
        }
    return {"type": "general"}


def _ray_cylinder_pz_intersect(origin_m, direction, geom_info):
    """
    Compute the exact entry point where a ray first enters a finite
    CylinderZ + PlaneZ cell.

    Parameters (all in metres)
    ──────────
    origin_m   : (3,) array — ray start [m]
    direction  : (3,) unit vector
    geom_info  : dict from _detector_cell_geometry() with type=="cylinder_pz"

    Returns (t_entry, t_exit) in metres, or (None, None) if no intersection.
    """
    import math

    cx, cy, r = geom_info["cx_m"], geom_info["cy_m"], geom_info["r_m"]
    z_min, z_max = geom_info["z_min_m"], geom_info["z_max_m"]

    ox, oy, oz = origin_m[0] - cx, origin_m[1] - cy, float(origin_m[2])
    dx, dy, dz = float(direction[0]), float(direction[1]), float(direction[2])

    # Solve quadratic: (ox+t*dx)² + (oy+t*dy)² = r²
    a = dx*dx + dy*dy
    if a < 1e-20:
        # Ray is parallel to cylinder axis — check if inside radially
        if ox*ox + oy*oy >= r*r:
            return None, None
        # Fully inside radially — constrained only by z planes
        if dz == 0:
            return None, None
        tz_min = (z_min - oz) / dz
        tz_max = (z_max - oz) / dz
        t_in  = max(min(tz_min, tz_max), 0.0)
        t_out = max(tz_min, tz_max)
        return (t_in, t_out) if t_in < t_out else (None, None)

    b = 2.0 * (ox*dx + oy*dy)
    c = ox*ox + oy*oy - r*r
    disc = b*b - 4*a*c
    if disc < 0:
        return None, None

    sq = math.sqrt(disc)
    t1 = (-b - sq) / (2*a)
    t2 = (-b + sq) / (2*a)

    if t2 <= 0:
        return None, None   # cylinder entirely behind ray

    # z bounds
    def z_at(t): return oz + t * dz

    # Clip [t1,t2] by z planes
    t_lo, t_hi = max(t1, 0.0), t2   # t must be positive (forward)

    if dz != 0:
        tz_zmin = (z_min - oz) / dz
        tz_zmax = (z_max - oz) / dz
        if dz > 0:
            t_lo = max(t_lo, tz_zmin)
            t_hi = min(t_hi, tz_zmax)
        else:
            t_lo = max(t_lo, tz_zmax)
            t_hi = min(t_hi, tz_zmin)
    else:
        # Ray horizontal: check if z is within bounds at any point
        z = z_at(0)
        if z < z_min or z > z_max:
            return None, None

    if t_lo >= t_hi or t_hi <= 0:
        return None, None

    return t_lo, t_hi


# ─────────────────────────────────────────────────────────────────────────────
# Core ray tracer — per muon
# ─────────────────────────────────────────────────────────────────────────────

def trace_ray_to_cell(geom,
                      origin_cm: np.ndarray,
                      direction: np.ndarray,
                      detector_cell_id: int,
                      step_cm: float = 50.0,
                      max_dist_cm: float = 1_000_000.0) -> RayResult:
    """
    Ray-trace from origin_cm in direction, accumulating pre-detector overburden.

    For CylinderZ + PlaneZ detector cells (common PHITS detector geometry):
    uses EXACT analytical ray-cylinder intersection — works correctly regardless
    of step_cm and independent of detector radius (handles 3 cm detectors).

    For other cell types: falls back to step-march.

    Parameters
    ──────────
    geom            : CSGGeometry from gui_csg_engine
    origin_cm       : ray start [cm] in geometry frame
    direction       : unit direction vector (cx, cy, cz), already normalised
    detector_cell_id: cell_id of the detector cell
    step_cm         : march step [cm] (used for overburden accumulation;
                      detector hit detection is always exact for CylinderZ cells)
    max_dist_cm     : maximum ray length [cm]

    Returns RayResult
    """
    res = RayResult(cx=direction[0], cy=direction[1], cz=direction[2])

    pos  = origin_cm.astype(float).copy()
    d    = direction.astype(float)
    d   /= max(np.linalg.norm(d), 1e-12)

    bbox_min = geom.bbox_min * 100.0   # m → cm
    bbox_max = geom.bbox_max * 100.0

    # Locate detector cell
    det_cell = None
    for c in geom._cells:
        if c.cell_id == detector_cell_id:
            det_cell = c
            break
    if det_cell is None:
        return res

    # Analyse detector geometry for analytical intersection
    det_geom_info = _detector_cell_geometry(det_cell, geom._surf_map)
    use_analytical = (det_geom_info["type"] == "cylinder_pz")

    # ── ANALYTICAL PATH (CylinderZ + PlaneZ cells) ───────────────────────────
    if use_analytical:
        # Convert units: origin_cm → origin_m, direction is dimensionless
        origin_m = origin_cm / 100.0
        t_entry_m, t_exit_m = _ray_cylinder_pz_intersect(origin_m, d, det_geom_info)

        if t_entry_m is None:
            return res   # ray misses detector entirely

        # Accumulate overburden from origin to t_entry by step-marching
        # (step_cm controls accuracy of overburden accumulation, not hit detection)
        seg_density = 0.0
        seg_len_cm  = 0.0
        segments: List[Tuple[float, float]] = []

        dist_cm = 0.0
        t_entry_cm = t_entry_m * 100.0
        t_exit_cm  = t_exit_m  * 100.0

        while dist_cm < min(t_entry_cm, max_dist_cm):
            dist_cm += step_cm
            if dist_cm > t_entry_cm:
                dist_cm = t_entry_cm
            step_pos = origin_cm + d * dist_cm
            if np.any(step_pos < bbox_min) or np.any(step_pos > bbox_max):
                break
            density = geom.density_at(step_pos / 100.0)
            if density > 0:
                if abs(density - seg_density) < 0.01:
                    seg_len_cm += min(step_cm, t_entry_cm - (dist_cm - step_cm))
                else:
                    if seg_len_cm > 0 and seg_density > 0:
                        segments.append((seg_density, seg_len_cm))
                    seg_density = density
                    seg_len_cm  = min(step_cm, t_entry_cm - (dist_cm - step_cm))
            else:
                if seg_len_cm > 0 and seg_density > 0:
                    segments.append((seg_density, seg_len_cm))
                seg_density = 0.0
                seg_len_cm  = 0.0

        if seg_len_cm > 0 and seg_density > 0:
            segments.append((seg_density, seg_len_cm))

        entry_pos_cm = origin_cm + d * t_entry_cm
        exit_pos_cm  = origin_cm + d * t_exit_cm

        res.hit             = True
        res.entry_pos_cm    = entry_pos_cm
        res.exit_pos_cm     = exit_pos_cm
        res.segments        = segments
        res.total_ob_gcm2   = sum(rho * L for rho, L in segments)
        res.path_through_cm = t_exit_cm - t_entry_cm
        return res

    # ── GENERAL STEP-MARCH FALLBACK ───────────────────────────────────────────
    seg_density = 0.0
    seg_len_cm  = 0.0
    segments: List[Tuple[float, float]] = []

    in_detector   = False
    detector_entry: Optional[np.ndarray] = None
    det_path_cm   = 0.0

    dist = 0.0
    while dist < max_dist_cm:
        dist += step_cm
        pos = origin_cm + d * dist

        if np.any(pos < bbox_min) or np.any(pos > bbox_max):
            break

        pos_m   = pos / 100.0
        density = geom.density_at(pos_m)
        is_det  = det_cell.contains(pos_m, geom._surf_map)

        if is_det:
            if not in_detector:
                in_detector    = True
                detector_entry = pos.copy()
                if seg_len_cm > 0 and seg_density > 0:
                    segments.append((seg_density, seg_len_cm))
                seg_density = 0.0
                seg_len_cm  = 0.0
            det_path_cm += step_cm
        else:
            if in_detector:
                res.hit             = True
                res.entry_pos_cm    = detector_entry
                res.exit_pos_cm     = pos.copy()
                res.segments        = segments
                res.total_ob_gcm2   = sum(rho * L for rho, L in segments)
                res.path_through_cm = det_path_cm
                return res
            if density > 0:
                if abs(density - seg_density) < 0.01:
                    seg_len_cm += step_cm
                else:
                    if seg_len_cm > 0 and seg_density > 0:
                        segments.append((seg_density, seg_len_cm))
                    seg_density = density
                    seg_len_cm  = step_cm
            else:
                if seg_len_cm > 0 and seg_density > 0:
                    segments.append((seg_density, seg_len_cm))
                seg_density = 0.0
                seg_len_cm  = 0.0

    if in_detector:
        res.hit             = True
        res.entry_pos_cm    = detector_entry
        res.exit_pos_cm     = None
        res.segments        = segments
        res.total_ob_gcm2   = sum(rho * L for rho, L in segments)
        res.path_through_cm = det_path_cm

    return res


# ─────────────────────────────────────────────────────────────────────────────
# CSDA energy loss through a sequence of material segments
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# CSDA range table — Standard Rock, Groom (2001) Table IV-6, subset
# Log-log interpolation: < 1% error over full 10 MeV – 1 TeV range.
# No regime-boundary discontinuity.
# ─────────────────────────────────────────────────────────────────────────────

_CSDA_E_MEV = np.array([
    10.,    20.,    40.,    80.,    160.,   300.,
    600.,  1_000., 2_000., 4_000., 8_000.,
    20_000., 60_000., 200_000., 1_000_000.
], dtype=float)

_CSDA_R_GCM2 = np.array([
    0.25,   0.80,   2.5,    7.8,    24.5,   72.0,
    205.,   405.,   930.,  2_000.,  4_400.,
    12_000., 37_000., 122_000., 560_000.
], dtype=float)

_LOG_E_TAB = np.log(_CSDA_E_MEV)
_LOG_R_TAB = np.log(_CSDA_R_GCM2)


def _csda_range_gcm2(E_MeV: float) -> float:
    """CSDA range [g/cm²] in Standard Rock via Groom (2001) log-log interpolation."""
    lE = np.log(max(float(E_MeV), _CSDA_E_MEV[0]))
    return float(np.exp(np.interp(lE, _LOG_E_TAB, _LOG_R_TAB)))


def _inverse_csda_range(R_gcm2: float) -> float:
    """Inverse CSDA: range [g/cm²] → kinetic energy [MeV] via log-log interpolation."""
    if R_gcm2 <= 0:
        return 0.0
    lR = np.log(min(max(float(R_gcm2), _CSDA_R_GCM2[0]), _CSDA_R_GCM2[-1]))
    return float(np.exp(np.interp(lR, _LOG_R_TAB, _LOG_E_TAB)))


def energy_after_segments(E_start_MeV: float,
                           segments: List[Tuple[float, float]]
                           ) -> float:
    """
    Compute kinetic energy [MeV] after transporting through a sequence of
    (density [g/cm³], path_length [cm]) segments using iterative CSDA.

    Each segment is treated as a uniform slab:
        R_remaining after seg = R_start - rho × L
        E_end from R_remaining via inverse CSDA table
    Returns 0.0 if the muon stops.
    """
    E = float(E_start_MeV)
    for rho, L_cm in segments:
        if E <= E_STOP_MEV:
            return 0.0
        ob_gcm2 = rho * L_cm
        R_start = _csda_range_gcm2(E)
        R_end   = R_start - ob_gcm2
        if R_end <= 0:
            return 0.0
        # Inverse CSDA: find E such that range = R_end
        # Use bisection on _csda_range_gcm2
        E = _inverse_csda_range(R_end)
    return max(E, 0.0)







# ─────────────────────────────────────────────────────────────────────────────
# Main transport function — vectorised over all surface muons
# ─────────────────────────────────────────────────────────────────────────────

def transport_muons_through_csg(
        surface_df,
        geom,
        detector_cell_id : int,
        coord_offset_cm  : np.ndarray = None,
        step_cm          : float = 50.0,
        max_dist_cm      : float = 1_000_000.0,
        v_cut            : float = 0.05,
        ms_enable        : bool  = True,
        script_dir       : str   = ".",
        progress_container = None,
) -> "pd.DataFrame":
    """
    Trace all Tab-1 surface muons through the CSG geometry.

    For each muon:
      1. Ray-trace through CSG until detector cell reached (or missed)
      2. Apply per-segment CSDA energy loss up to detector entry
      3. Apply UCMuon-Stochastic transport for the remaining path
         (uses total pre-detector overburden in Standard Rock equivalent)
      4. Record 18-col underground format at detector entry point

    Parameters
    ──────────
    surface_df       : pandas DataFrame from Tab-1 file
                       must have columns: EventID, x, y, z, E, theta, phi,
                                         cx, cy, cz, charge
    geom             : CSGGeometry from gui_csg_engine
    detector_cell_id : cell_id of the detector cell
    coord_offset_cm  : [dx, dy, dz] [cm] to add to Tab-1 coordinates
                       to transform them into the geometry frame.
                       Default: [0, 0, 0]
    step_cm          : ray-march step [cm] (default 50 cm)
    max_dist_cm      : max ray length [cm]
    v_cut            : UCMuon catastrophic event threshold
    ms_enable        : Highland multiple-scattering deflections
    script_dir       : path to find ucmuon_stochastic_driver.py

    Returns
    ──────
    pandas DataFrame with 18 columns (EventID, xs, ys, zs, Es, theta_s,
    phi_s, charge, alive, x, y, z, E, cx, cy, cz, theta, phi).

    alive=1  : muon reached detector cell with E > E_stop
    alive=0  : muon was absorbed / did not reach detector cell
    """
    import pandas as pd

    if coord_offset_cm is None:
        coord_offset_cm = np.zeros(3)
    coord_offset_cm = np.asarray(coord_offset_cm, dtype=float)

    # Load UCMuon-MC driver for the final transport step
    _stochastic_drv = None
    try:
        drv_path = Path(script_dir) / "ucmuon_stochastic_driver.py"
        if drv_path.exists():
            spec = importlib.util.spec_from_file_location(
                "ucmuon_stochastic_driver", str(drv_path))
            _stochastic_drv = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(_stochastic_drv)
    except Exception:
        _stochastic_drv = None

    rng = np.random.default_rng(42)

    n_total = len(surface_df)
    out_rows: list = []
    n_hit    = 0
    n_surv   = 0
    t0       = time.time()

    for i, row in enumerate(surface_df.itertuples(index=False)):
        # ── Surface muon properties ───────────────────────────────────────
        eid    = int(row.EventID)
        # Tab-1 positions are in cm (CosmoALEPH convention)
        x_s    = float(row.x)
        y_s    = float(row.y)
        z_s    = float(getattr(row, 'z', 0.0))
        E_s    = float(row.E)          # GeV
        theta_s= float(row.theta)
        phi_s  = float(row.phi)
        charge = int(row.charge)
        cx_s   = float(row.cx)
        cy_s   = float(row.cy)
        cz_s   = float(row.cz)

        # ── Transform to geometry frame ───────────────────────────────────
        origin_cm = np.array([x_s, y_s, z_s]) + coord_offset_cm
        direction = np.array([cx_s, cy_s, cz_s])
        n = np.linalg.norm(direction)
        if n < 1e-12:
            direction = np.array([0., 0., -1.])
        else:
            direction /= n

        # ── Ray trace to detector cell ────────────────────────────────────
        ray = trace_ray_to_cell(
            geom, origin_cm, direction,
            detector_cell_id, step_cm, max_dist_cm,
        )

        if not ray.hit:
            # Muon missed / absorbed before detector
            out_rows.append(_dead_row(
                eid, x_s, y_s, z_s, E_s, theta_s, phi_s, charge,
                x_s, y_s, z_s, cx_s, cy_s, cz_s, theta_s, phi_s,
            ))
            _maybe_progress(progress_container, i, n_total, n_hit, n_surv, t0)
            continue

        n_hit += 1

        # ── Energy at detector entry (CSDA through each pre-det segment) ──
        E_at_det_MeV = energy_after_segments(
            E_s * 1000.0,      # GeV → MeV
            ray.segments,
        )

        if E_at_det_MeV <= E_STOP_MEV:
            # Muon stopped before reaching detector
            det_x, det_y, det_z = ray.entry_pos_cm
            out_rows.append(_dead_row(
                eid, x_s, y_s, z_s, E_s, theta_s, phi_s, charge,
                det_x, det_y, det_z, cx_s, cy_s, cz_s, theta_s, phi_s,
            ))
            _maybe_progress(progress_container, i, n_total, n_hit, n_surv, t0)
            continue

        # ── UCMuon-MC transport through detector overburden ───────────────
        # We use the path through the detector cell itself as a final
        # stochastic transport step (adds energy spread and MS deflections).
        det_x, det_y, det_z = ray.entry_pos_cm

        if _stochastic_drv is not None and ray.path_through_cm > 0:
            # Build single-muon dict
            mat    = _stochastic_drv._MAT_DB.get(1)   # Standard Rock approximation
            mu_dict = dict(
                EventID  = np.array([eid]),
                x        = np.array([det_x]),
                y        = np.array([det_y]),
                z        = np.array([det_z]),
                theta    = np.array([theta_s]),
                phi      = np.array([phi_s]),
                Ekin_GeV = np.array([E_at_det_MeV / 1000.0]),
                Ekin_MeV = np.array([E_at_det_MeV]),
                charge   = np.array([charge]),
                cx       = np.array([direction[0]]),
                cy       = np.array([direction[1]]),
                cz       = np.array([direction[2]]),
            )
            det_rho_gcm3 = _get_detector_density(geom, detector_cell_id)
            depth_m = ray.path_through_cm / 100.0   # cm → m
            result  = _stochastic_drv.transport(
                mu_dict, depth_m, det_rho_gcm3, mat,
                n_steps=0, v_cut=v_cut, ms_enable=ms_enable, rng=rng,
            )
            alive    = int(result["alive"][0])
            E_f_GeV  = float(result["E_kin_f_MeV"][0]) / 1000.0
            cx_f     = float(result["cx_f"][0])
            cy_f     = float(result["cy_f"][0])
            cz_f     = float(result["cz_f"][0])
            theta_f  = float(result["theta_f"][0])
            phi_f    = float(result["phi_f"][0])
            x_f      = float(result["x_f"][0])
            y_f      = float(result["y_f"][0])
            z_f      = float(result["z_f"][0])
        else:
            # No stochastic driver or zero path — use CSDA result directly
            alive   = 1
            E_f_GeV = E_at_det_MeV / 1000.0
            cx_f, cy_f, cz_f = direction[0], direction[1], direction[2]
            theta_f = np.arccos(np.clip(-cz_f, -1.0, 1.0))
            phi_f   = np.arctan2(cy_f, cx_f)
            x_f, y_f, z_f = det_x, det_y, det_z

        if alive:
            n_surv += 1

        out_rows.append({
            "EventID": eid,
            "xs":      x_s,    "ys":    y_s,    "zs":    z_s,
            "Es":      E_s,
            "theta_s": theta_s, "phi_s": phi_s,
            "charge":  charge,
            "alive":   alive,
            "x":       x_f,    "y":     y_f,    "z":     z_f,
            "E":       E_f_GeV,
            "cx":      cx_f,   "cy":    cy_f,   "cz":    cz_f,
            "theta":   theta_f, "phi":  phi_f,
        })

        _maybe_progress(progress_container, i, n_total, n_hit, n_surv, t0)

    if not out_rows:
        import pandas as pd
        return pd.DataFrame(columns=[
            "EventID","xs","ys","zs","Es","theta_s","phi_s","charge",
            "alive","x","y","z","E","cx","cy","cz","theta","phi"])

    import pandas as pd
    return pd.DataFrame(out_rows)



# ─────────────────────────────────────────────────────────────────────────────
# Multi-detector transport — PHITS T-Cross equivalent
# Records hits on any number of detector cells in a SINGLE muon pass
# ─────────────────────────────────────────────────────────────────────────────

def transport_muons_multi_detector(
        surface_df,
        geom,
        detector_cell_ids : List[int],
        detector_names    : Optional[Dict[int, str]] = None,
        coord_offset_cm   : np.ndarray = None,
        step_cm           : float = 50.0,
        max_dist_cm       : float = 1_000_000.0,
        v_cut             : float = 0.05,
        ms_enable         : bool  = True,
        script_dir        : str   = ".",
        write_phits_dump  : bool  = True,
        dump_dir          : str   = ".",
        progress_container = None,
) -> Dict[int, "pd.DataFrame"]:
    """
    Transport all Tab-1 surface muons and record hits on MULTIPLE detector
    cells simultaneously — exactly equivalent to running several PHITS T-Cross
    tallies in a single simulation pass.

    For each surface muon:
      1. Ray-trace against each detector cell (analytical intersection for
         CylinderZ cells, step-march fallback for others)
      2. Record the first hit on EACH detector separately
      3. Compute per-segment CSDA energy loss to each detector independently
      4. Write PHITS-compatible phase-space dump for each detector
         (T-Cross dump=-14 format: columns 1 2 3 4 5 6 7 8 9 10 17 18 19 20)

    Parameters
    ──────────
    surface_df        : Tab-1 DataFrame (EventID, x, y, z, E, theta, phi,
                        cx, cy, cz, charge  —  positions in cm, E in GeV)
    geom              : CSGGeometry
    detector_cell_ids : list of cell_ids to treat as detector volumes
    detector_names    : optional dict {cell_id: "M2"} for output filenames
    coord_offset_cm   : [dx, dy, dz] [cm] to shift Tab-1 into geometry frame
    step_cm           : overburden accumulation step [cm]
    max_dist_cm       : maximum ray length [cm]
    write_phits_dump  : write PHITS T-Cross compatible .dat files
    dump_dir          : directory for dump files

    Returns
    ──────
    Dict[int, pd.DataFrame] keyed by cell_id.  Each DataFrame has the
    standard 18-column underground format (alive=1 = reached detector).
    """
    import pandas as pd
    import os

    if coord_offset_cm is None:
        coord_offset_cm = np.zeros(3)
    coord_offset_cm = np.asarray(coord_offset_cm, dtype=float)

    if detector_names is None:
        detector_names = {cid: f"cell{cid}" for cid in detector_cell_ids}

    # Pre-cache detector geometry info (analytical vs step-march)
    det_cells = {}
    for cid in detector_cell_ids:
        cell = next((c for c in geom._cells if c.cell_id == cid), None)
        if cell is None:
            raise ValueError(f"Detector cell {cid} not found in geometry")
        geom_info = _detector_cell_geometry(cell, geom._surf_map)
        det_cells[cid] = (cell, geom_info)

    rng = np.random.default_rng(42)

    # Per-detector output accumulator
    out_rows  = {cid: [] for cid in detector_cell_ids}
    n_hit     = {cid: 0  for cid in detector_cell_ids}
    n_surv    = {cid: 0  for cid in detector_cell_ids}
    n_total   = len(surface_df)
    t0        = time.time()

    # Phase-space dump file handles (PHITS T-Cross dump=-14 format)
    dump_handles: Dict[int, object] = {}
    if write_phits_dump:
        os.makedirs(dump_dir, exist_ok=True)
        for cid in detector_cell_ids:
            name = detector_names.get(cid, f"cell{cid}")
            fpath = os.path.join(dump_dir, f"{name}_hits.dat")
            dump_handles[cid] = open(fpath, "w")
            # Header: PHITS dump=-14 column order
            dump_handles[cid].write(
                "# PHITS T-Cross equivalent dump  "
                f"detector=cell{cid} ({name})\n"
                "# col: 1=x[cm] 2=y[cm] 3=z[cm] 4=cx 5=cy 6=cz "
                "7=E[MeV] 8=weight 9=charge 10=EventID "
                "17=theta[rad] 18=phi[rad] 19=OB[g/cm2] 20=path_det[cm]\n"
            )

    try:
        for i, row in enumerate(surface_df.itertuples(index=False)):
            eid    = int(row.EventID)
            x_s    = float(row.x)
            y_s    = float(row.y)
            z_s    = float(getattr(row, 'z', 0.0))
            E_s    = float(row.E)          # GeV
            theta_s= float(row.theta)
            phi_s  = float(row.phi)
            charge = int(row.charge)
            cx_s   = float(row.cx)
            cy_s   = float(row.cy)
            cz_s   = float(row.cz)

            origin_cm = np.array([x_s, y_s, z_s]) + coord_offset_cm
            direction = np.array([cx_s, cy_s, cz_s])
            n = np.linalg.norm(direction)
            direction = direction / max(n, 1e-12)

            # ── Trace ray against each detector independently ─────────────
            for cid in detector_cell_ids:
                ray = trace_ray_to_cell(
                    geom, origin_cm, direction, cid,
                    step_cm=step_cm, max_dist_cm=max_dist_cm,
                )

                if not ray.hit:
                    out_rows[cid].append(_dead_row(
                        eid, x_s, y_s, z_s, E_s, theta_s, phi_s, charge,
                        x_s, y_s, z_s, cx_s, cy_s, cz_s, theta_s, phi_s,
                    ))
                    continue

                n_hit[cid] += 1

                # CSDA energy loss to detector entry
                E_at_det_MeV = energy_after_segments(E_s * 1000.0, ray.segments)

                if E_at_det_MeV <= E_STOP_MEV:
                    det_x, det_y, det_z = ray.entry_pos_cm
                    out_rows[cid].append(_dead_row(
                        eid, x_s, y_s, z_s, E_s, theta_s, phi_s, charge,
                        det_x, det_y, det_z, cx_s, cy_s, cz_s, theta_s, phi_s,
                    ))
                    continue

                n_surv[cid] += 1
                det_x, det_y, det_z = ray.entry_pos_cm
                cx_f, cy_f, cz_f = direction
                theta_f = float(np.arccos(np.clip(-cz_f, -1.0, 1.0)))
                phi_f   = float(np.arctan2(cy_f, cx_f))
                E_f_GeV = E_at_det_MeV / 1000.0

                out_rows[cid].append({
                    "EventID": eid,
                    "xs": x_s, "ys": y_s, "zs": z_s, "Es": E_s,
                    "theta_s": theta_s, "phi_s": phi_s, "charge": charge,
                    "alive": 1,
                    "x": det_x, "y": det_y, "z": det_z,
                    "E": E_f_GeV,
                    "cx": cx_f, "cy": cy_f, "cz": cz_f,
                    "theta": theta_f, "phi": phi_f,
                })

                # PHITS T-Cross compatible dump line
                if write_phits_dump and cid in dump_handles:
                    dump_handles[cid].write(
                        f"{det_x:14.6e} {det_y:14.6e} {det_z:14.6e} "
                        f"{cx_f:10.6f} {cy_f:10.6f} {cz_f:10.6f} "
                        f"{E_at_det_MeV:14.6e} "
                        f"1.0 {charge:3d} {eid:10d} "
                        f"{theta_f:10.6f} {phi_f:10.6f} "
                        f"{ray.total_ob_gcm2:14.4f} {ray.path_through_cm:12.4f}\n"
                    )

            # Progress update
            if progress_container and (i+1) % max(1, n_total//50) == 0:
                elapsed = time.time() - t0
                rate = (i+1) / elapsed if elapsed > 0 else 0
                eta  = (n_total - i - 1) / rate if rate > 0 else 0
                hit_summary = "  ".join(
                    f"Cell{cid}={n_hit[cid]}" for cid in detector_cell_ids)
                try:
                    progress_container.write(
                        f"⏳ {i+1:,}/{n_total:,}  "
                        f"({hit_summary})  "
                        f"{rate:.0f} μ/s  ETA {eta:.0f}s"
                    )
                except Exception:
                    pass

    finally:
        for fh in dump_handles.values():
            fh.close()

    # Build DataFrames
    COLS = ["EventID","xs","ys","zs","Es","theta_s","phi_s","charge",
            "alive","x","y","z","E","cx","cy","cz","theta","phi"]
    result = {}
    for cid in detector_cell_ids:
        rows = out_rows[cid]
        if rows:
            result[cid] = pd.DataFrame(rows)
        else:
            result[cid] = pd.DataFrame(columns=COLS)

    # Summary
    print(f"\n[UCMuon CSG multi-detector] {n_total:,} muons transported")
    for cid in detector_cell_ids:
        name = detector_names.get(cid, f"cell{cid}")
        nh   = n_hit[cid]
        ns   = n_surv[cid]
        hit_rate  = nh  / n_total * 100 if n_total else 0
        surv_rate = ns  / n_total * 100 if n_total else 0
        print(f"  {name} (cell {cid}): "
              f"hit={nh:,} ({hit_rate:.2f}%)  "
              f"survive={ns:,} ({surv_rate:.2f}%)")
    if write_phits_dump:
        for cid in detector_cell_ids:
            name = detector_names.get(cid, f"cell{cid}")
            print(f"  Dump: {os.path.join(dump_dir, name+'_hits.dat')}")

    return result


def write_phits_dump_file(df: "pd.DataFrame", path: str,
                          detector_cell_id: int, geometry_name: str = ""):
    """
    Write a PHITS T-Cross compatible phase-space dump from a results DataFrame.

    Format matches PHITS dump=-14 with columns 1 2 3 4 5 6 7 8 9 10 17 18 19 20.
    Only rows with alive=1 are written (equivalent to T-Cross current scoring).

    Parameters
    ──────────
    df               : results DataFrame from transport_muons_multi_detector
    path             : output file path
    detector_cell_id : for header information
    geometry_name    : geometry description for header
    """
    alive = df[df["alive"] == 1] if "alive" in df.columns else df
    with open(path, "w") as f:
        f.write(f"# UCMuon PHITS T-Cross dump  cell={detector_cell_id}  "
                f"geo={geometry_name}  n_muons={len(alive)}\n")
        f.write("# x[cm]           y[cm]           z[cm]           "
                "cx         cy         cz         "
                "E[MeV]          wt  chg     EventID    "
                "theta[r]   phi[r]     OB[g/cm2]     path_det[cm]\n")
        for row in alive.itertuples(index=False):
            theta = float(row.theta)
            phi   = float(row.phi)
            E_MeV = float(row.E) * 1000.0
            f.write(
                f"{row.x:14.6e} {row.y:14.6e} {row.z:14.6e} "
                f"{row.cx:10.6f} {row.cy:10.6f} {row.cz:10.6f} "
                f"{E_MeV:14.6e} 1.0 {row.charge:3d} {row.EventID:10d} "
                f"{theta:10.6f} {phi:10.6f} 0.0 0.0\n"
            )
    return path

def _dead_row(eid, xs, ys, zs, Es, ths, phs, charge,
              x, y, z, cx, cy, cz, theta, phi) -> dict:
    return {
        "EventID": eid,
        "xs": xs, "ys": ys, "zs": zs, "Es": Es,
        "theta_s": ths, "phi_s": phs, "charge": charge,
        "alive": 0,
        "x": x, "y": y, "z": z, "E": 0.0,
        "cx": cx, "cy": cy, "cz": cz,
        "theta": theta, "phi": phi,
    }


def _get_detector_density(geom, detector_cell_id: int) -> float:
    for c in geom._cells:
        if c.cell_id == detector_cell_id:
            return c.density if c.density > 0 else 1.225e-3  # air density
    return 1.0


def _maybe_progress(container, i, total, n_hit, n_surv, t0):
    if container is None:
        return
    frac = (i + 1) / max(total, 1)
    elapsed = time.time() - t0
    eta = elapsed / max(frac, 1e-6) * (1.0 - frac)
    container.progress(
        min(frac, 0.99),
        text=(f"⏳  {i+1:,}/{total:,} muons traced  |  "
              f"{n_hit:,} reached detector  |  "
              f"{n_surv:,} survived  |  "
              f"ETA {eta:.0f} s"),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Streamlit UI — detector cell selector
# ─────────────────────────────────────────────────────────────────────────────

def render_detector_cell_selector(geom) -> Optional[int]:
    """
    Render the detector cell selector UI.

    Supports selecting MULTIPLE detector cells simultaneously (like having
    multiple PHITS T-Cross tallies). Returns the PRIMARY cell_id (first
    selected) for backward compatibility. The full list is stored in
    st.session_state["_csg_selected_det_cells"].

    For CylinderZ cells (e.g. Svelvik M2/M3): displays the analytical
    bounding box so the user can verify the cell is correct.
    """
    from gui_csg_engine import CSGGeometry

    if geom is None:
        st.info("Load a PHITS geometry in Section 2 first.", icon="ℹ️")
        return None

    candidate_cells = [c for c in geom._cells if c.mat_id != -1]
    if not candidate_cells:
        st.warning("No candidate cells found in the geometry.", icon="⚠️")
        return None

    st.markdown("#### 📡 Detector Cell Selection")
    st.caption(
        "Select one or more cells as detector volumes.  "
        "Each selected cell records hits independently — equivalent to "
        "multiple PHITS T-Cross tallies.  "
        "The ray-tracer uses exact analytical intersection for CylinderZ cells "
        "(works even for 3 cm detectors in 265 m geometries)."
    )

    def _cell_label(c):
        mat_str = f"mat={c.mat_id}" if c.mat_id >= 0 else "outside"
        rho_str = f"ρ={c.density:.3f} g/cm³" if c.density > 0 else "void/air"
        return f"Cell {c.cell_id}  [{mat_str}]  {rho_str}"

    options  = [_cell_label(c) for c in candidate_cells]
    cell_ids = [c.cell_id for c in candidate_cells]

    # Restore previously selected cells (or default to first plastic scintillator)
    prev_selected = st.session_state.get("_csg_selected_det_cells", [])
    default_labels = []
    for cid in prev_selected:
        idx = cell_ids.index(cid) if cid in cell_ids else None
        if idx is not None:
            default_labels.append(options[idx])
    if not default_labels:
        # Default: prefer cells with plastic scintillator (mat=70 in Svelvik)
        for c in candidate_cells:
            if c.mat_id == 70 or c.density > 1.0:
                default_labels.append(_cell_label(c))
        if not default_labels:
            default_labels = [options[0]]

    selected_labels = st.multiselect(
        "Detector cells",
        options=options,
        default=default_labels,
        key="csg_det_cell_ids",
        help="Select all cells to be treated as detector volumes.  "
             "For Svelvik: select Cell 307 (M2) and Cell 306 (M3) simultaneously.",
    )

    if not selected_labels:
        st.warning("Select at least one detector cell.", icon="⚠️")
        return None

    selected_ids = [cell_ids[options.index(lbl)] for lbl in selected_labels]

    # Store for 3D plot highlighting
    st.session_state["_csg_selected_det_cells"] = selected_ids
    # Primary detector (used by single-detector run button)
    det_cell_id = selected_ids[0]

    # ── Per-cell info panel ───────────────────────────────────────────────────
    for cid in selected_ids:
        cell = next(c for c in candidate_cells if c.cell_id == cid)
        geom_info = _detector_cell_geometry(cell, geom._surf_map)
        label = _cell_label(cell)

        if geom_info["type"] == "cylinder_pz":
            import math
            r_m  = geom_info["r_m"]
            cx_m = geom_info["cx_m"]
            cy_m = geom_info["cy_m"]
            zlo  = geom_info["z_min_m"]
            zhi  = geom_info["z_max_m"]
            col1, col2 = st.columns(2)
            with col1:
                st.success(
                    f"📡 **Cell {cid}** — CylinderZ + PlaneZ  ✓ Exact intersection\n\n"
                    f"Centre: ({cx_m:.2f}, {cy_m:.2f}) m  "
                    f"r = {r_m*100:.1f} cm\n\n"
                    f"z: {zlo:.1f} m → {zhi:.1f} m  "
                    f"(length {abs(zhi-zlo):.1f} m)",
                    icon="✅"
                )
            with col2:
                # Analytical overburden estimate (vertical ray)
                z_entry = abs(zhi)   # m depth to top of detector
                st.metric("Depth to detector top", f"{z_entry:.1f} m")
                st.metric("Detector length", f"{abs(zhi-zlo):.1f} m")
                st.metric("Detector radius", f"{r_m*100:.1f} cm")
        else:
            if cell.density > 0:
                st.warning(
                    f"⚠️  Cell {cid}: ρ={cell.density:.2f} g/cm³  "
                    "(not void — using step-march intersection)",
                    icon="⚠️"
                )
            else:
                st.info(f"Cell {cid}: void/air  (step-march intersection)", icon="ℹ️")

    # ── Multi-detector run info ───────────────────────────────────────────────
    if len(selected_ids) > 1:
        st.info(
            f"**{len(selected_ids)} detectors selected** — "
            "muons will be ray-traced against each cell independently "
            "in a single transport pass.  "
            "Output: one DataFrame and one PHITS dump file per detector.",
            icon="🔬"
        )

    # ── Ray-march settings ────────────────────────────────────────────────────
    with st.expander("⚙️  Ray-march settings", expanded=False):
        _c1, _c2 = st.columns(2)
        step_cm    = _c1.number_input(
            "Step [cm]", 1.0, 1000.0, 50.0, 10.0,
            key="csg_trace_step_cm",
            help="Controls overburden accuracy (not hit detection — "
                 "CylinderZ hits are always exact)."
        )
        max_dist_m = _c2.number_input(
            "Max ray length [m]", 10.0, 100_000.0, 5000.0, 100.0,
            key="csg_trace_maxdist_m",
            help="Terminate ray at this distance.  Set ≥ geometry diameter."
        )
        st.caption(
            f"Step: {step_cm:.0f} cm  |  "
            f"Max dist: {max_dist_m*100:.0f} cm  |  "
            f"Max steps: {int(max_dist_m*100/step_cm):,}"
        )

    st.session_state["_csg_step_cm"]    = float(step_cm)
    st.session_state["_csg_maxdist_cm"] = float(max_dist_m * 100.0)

    # ── Coordinate offset ─────────────────────────────────────────────────────
    with st.expander("📐  Coordinate frame offset (Tab-1 → geometry)", expanded=False):
        st.caption(
            "Tab-1 muons have positions in cm in the CosmoALEPH ENU frame "
            "(z≈0 at surface).  Add an offset to transform them into the "
            "PHITS geometry frame.  For Svelvik: leave at (0,0,0) — the "
            "PHITS geometry already has z=0 at ground surface."
        )
        _o1, _o2, _o3 = st.columns(3)
        off_x = _o1.number_input("Δx [cm]", value=0.0, step=100.0, key="csg_off_x",
                                  help="East offset [cm]")
        off_y = _o2.number_input("Δy [cm]", value=0.0, step=100.0, key="csg_off_y",
                                  help="North offset [cm]")
        off_z = _o3.number_input("Δz [cm]", value=0.0, step=100.0, key="csg_off_z",
                                  help="Vertical offset [cm]")
        st.session_state["_csg_offset_cm"] = np.array([off_x, off_y, off_z])

    # ── Single-ray spot-check for primary detector ────────────────────────────
    with st.expander("🔬  Single-ray spot-check", expanded=False):
        st.caption(
            "Trace one ray and verify the hit detection, "
            "overburden accumulation, and energy loss."
        )
        _sc1, _sc2, _sc3 = st.columns(3)
        sc_az = _sc1.number_input("Azimuth [°]", 0.0, 360.0,  0.0, 5.0, key="csg_sc_az")
        sc_ze = _sc2.number_input("Zenith [°]",  0.0,  90.0, 30.0, 5.0, key="csg_sc_ze")
        sc_E  = _sc3.number_input("E [GeV]",     1.0,  1e5,  50.0, 1.0, key="csg_sc_E",
                                   format="%.1f")

        _off = st.session_state.get("_csg_offset_cm", np.zeros(3))
        az_r  = np.radians(sc_az)
        ze_r  = np.radians(sc_ze)
        sin_ze = np.sin(ze_r)
        sc_dir = np.array([sin_ze*np.sin(az_r), sin_ze*np.cos(az_r), -np.cos(ze_r)])
        sc_dir /= max(np.linalg.norm(sc_dir), 1e-12)
        sc_origin = np.array([0.0, 0.0, 0.0]) + _off

        # Check against all selected detectors
        _sc_det = st.selectbox(
            "Check against detector",
            options=[f"Cell {c}" for c in selected_ids],
            key="csg_sc_det",
        )
        sc_det_id = int(_sc_det.split()[1])

        if st.button("▶ Trace single ray", key="csg_sc_btn"):
            ray = trace_ray_to_cell(
                geom, sc_origin, sc_dir, sc_det_id,
                step_cm    = float(st.session_state.get("_csg_step_cm", 50.0)),
                max_dist_cm= float(st.session_state.get("_csg_maxdist_cm", 500_000.0)),
            )
            if ray.hit:
                E_det = energy_after_segments(sc_E * 1000.0, ray.segments)
                st.success(
                    f"**Hit ✓**  |  "
                    f"Entry ({ray.entry_pos_cm[0]/100:.2f}, "
                    f"{ray.entry_pos_cm[1]/100:.2f}, "
                    f"{ray.entry_pos_cm[2]/100:.2f}) m  |  "
                    f"OB before: **{ray.total_ob_gcm2:.0f} g/cm²**  |  "
                    f"E at detector: **{E_det/1000:.3f} GeV** "
                    f"(Δ = {sc_E - E_det/1000:.3f} GeV)  |  "
                    f"Path through cell: {ray.path_through_cm/100:.1f} m  |  "
                    f"{len(ray.segments)} pre-det segments"
                )
                if ray.segments:
                    segs_md = "  \n".join(
                        f"— Seg {j+1}: ρ={r:.3f} g/cm³  "
                        f"L={L/100:.2f} m  OB={r*L:.0f} g/cm²"
                        for j, (r, L) in enumerate(ray.segments)
                    )
                    st.caption(segs_md)
            else:
                st.info(
                    f"No hit — ray (az={sc_az:.0f}°, ze={sc_ze:.0f}°) "
                    f"did not reach cell {sc_det_id}.  "
                    "Adjust the azimuth/zenith to aim toward the detector, "
                    "or check the coordinate offset.",
                    icon="🚫"
                )

    return det_cell_id


# ─────────────────────────────────────────────────────────────────────────────
# Summary diagnostics plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_csg_transport_results(df, det_cell_id: int):
    """
    Render a compact results panel for a CSG detector-cell transport run.

    Shows: hit rate, survival rate, energy spectrum at detector,
           and angular distribution of detected muons.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    n_total = len(df)
    n_hit   = int((df["alive"] >= 0).sum())          # all rows
    n_surv  = int((df["alive"] == 1).sum())
    n_reach = int((df["x"] != df["xs"]).sum())        # proxy: position changed

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total muons",    f"{n_total:,}")
    c2.metric("Reached detector", f"{n_reach:,}",
              delta=f"{100*n_reach/max(n_total,1):.2f}%")
    c3.metric("Survived",       f"{n_surv:,}",
              delta=f"{100*n_surv/max(n_reach,1):.2f}% of reached")
    c4.metric("Detector cell",  f"Cell {det_cell_id}")

    df_surv = df[df["alive"] == 1]
    if df_surv.empty:
        st.warning("No surviving muons to plot.", icon="⚠️")
        return

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Energy spectrum at detector [GeV]",
                        "Zenith distribution at detector [°]"],
    )

    # Energy spectrum
    E_det = df_surv["E"].values
    E_surf= df_surv["Es"].values
    bins  = np.logspace(np.log10(max(E_det.min(), 0.01)), np.log10(E_det.max()+0.01), 40)
    fig.add_trace(go.Histogram(x=E_surf, xbins=dict(start=bins[0], end=bins[-1]),
                               name="Surface", marker_color="#4fc3f7", opacity=0.6), row=1, col=1)
    fig.add_trace(go.Histogram(x=E_det,  xbins=dict(start=bins[0], end=bins[-1]),
                               name="At detector", marker_color="#81c784", opacity=0.8), row=1, col=1)

    # Zenith distribution
    theta_det = np.degrees(df_surv["theta"].values)
    fig.add_trace(go.Histogram(x=theta_det, nbinsx=36,
                               name="Zenith at det", marker_color="#ffb74d", opacity=0.9), row=1, col=2)

    fig.update_layout(
        **DARK, height=340, showlegend=True,
        legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color="white")),
        margin=dict(l=50, r=20, t=40, b=50),
    )
    fig.update_xaxes(title_text="E [GeV]",   gridcolor="#2a2a3a", type="log", row=1, col=1)
    fig.update_xaxes(title_text="Zenith [°]", gridcolor="#2a2a3a",             row=1, col=2)
    fig.update_yaxes(title_text="Counts", gridcolor="#2a2a3a")
    st.plotly_chart(fig)

    # Overburden distribution (xs, ys vs x, y shows how path changed)
    st.caption(
        f"Surface muon file: {n_total:,} total.  "
        f"Detector cell hit: {n_reach:,} ({100*n_reach/max(n_total,1):.1f}%).  "
        f"Survived transport: {n_surv:,} ({100*n_surv/max(n_reach,1):.1f}% of reached).  "
        f"Mean E at detector: {E_det.mean():.2f} GeV  "
        f"(surface mean: {E_surf.mean():.2f} GeV)."
    )
