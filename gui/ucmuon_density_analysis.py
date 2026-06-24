#!/usr/bin/env python3
"""
ucmuon_density_analysis.py  —  UCLouvain Muography Group
Density inversion module for muon tomography.

Three inversion methods are provided:

  1. Library matching (forward model)
     Recover ρ̂(az, el) by interpolating a library of simulated T_sim(ρ) maps
     produced at several discrete densities.  Robust but needs 3–5 full
     terrain forward simulations.
        → invert_density_map()

  2. Direct opacity inversion (analytical)
     Recover the column density (opacity) ϱ = ∫ρ·dl directly from a single
     measured transmission map T_data = Φ_target / Φ_open, using an analytic
     sea-level flux model + CSDA range (fast_flux_estimator).  Needs only ONE
     open-sky + ONE target measurement.  Mean density ρ̄ = ϱ / L follows once a
     line-of-sight path length L(az, el) is supplied.
        → invert_opacity_map() , opacity_to_density()

  3. Two-flux-map ratio
     Same physics as (2) but the transmission map is formed here from two raw
     flux maps (open-sky + target) instead of a precomputed T file.
        → transmission_from_flux_maps()

Path length L(az, el) for methods (2)/(3) may be obtained from:
  • a T_sim map at known ρ_sim      → path_length_from_tsim()
  • a terrain overburden map        → load_overburden_as_L()
  • or skipped entirely (report opacity only)

Author: Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026
"""

import numpy as np
from pathlib import Path

# Optional scipy — graceful fallback to numpy-only interpolation
try:
    from scipy.interpolate import interp1d as _scipy_interp1d
    _SCIPY_OK = True
except ImportError:
    _SCIPY_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_OPENSKY_THRESHOLD = 0.999   # T_sim values above this → open sky for all ρ
_GRADIENT_THRESHOLD = 1e-6   # |dT/dρ| below this → low-sensitivity flag


# ─────────────────────────────────────────────────────────────────────────────
# File I/O
# ─────────────────────────────────────────────────────────────────────────────

def load_transmission_map(filepath):
    """
    Parse a terrain_transmission.dat file into arrays and metadata.

    Returns (az_c, el_c, T_2d, metadata_dict).
      az_c         : 1-D float64 array (n_az,) — azimuth bin centres [deg]
      el_c         : 1-D float64 array (n_el,) — elevation bin centres [deg]
      T_2d         : 2-D float64 array (n_az, n_el) — transmission in [0, 1]
      metadata_dict: dict with keys 'density', 'det_lat', 'det_lon', 'det_alt',
                     'n_az', 'n_el', 'filepath'
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise ValueError(f"Transmission file not found: {filepath}")

    metadata = {
        "density":  None,
        "det_lat":  None,
        "det_lon":  None,
        "det_alt":  None,
        "n_az":     None,
        "n_el":     None,
        "filepath": str(filepath),
    }

    rows = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                _parse_header_line(line, metadata)
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                az, el, T = float(parts[0]), float(parts[1]), float(parts[2])
                rows.append((az, el, T))
            except ValueError:
                continue

    if not rows:
        raise ValueError(f"No data rows found in {filepath}")

    data = np.array(rows, dtype=np.float64)
    az_vals = data[:, 0]
    el_vals = data[:, 1]
    T_vals  = data[:, 2]

    az_c = np.unique(az_vals)
    el_c = np.unique(el_vals)
    n_az = len(az_c)
    n_el = len(el_c)

    # Validate grid completeness
    expected = n_az * n_el
    if len(rows) != expected:
        raise ValueError(
            f"Expected {expected} data rows ({n_az} az × {n_el} el) "
            f"but found {len(rows)} in {filepath}. File may be truncated."
        )

    # Build az/el index maps
    az_idx = {v: i for i, v in enumerate(az_c)}
    el_idx = {v: i for i, v in enumerate(el_c)}

    T_2d = np.full((n_az, n_el), np.nan, dtype=np.float64)
    for az, el, T in rows:
        ia = az_idx.get(np.round(az, 8))
        ie = el_idx.get(np.round(el, 8))
        if ia is not None and ie is not None:
            T_2d[ia, ie] = T

    # Update metadata
    if metadata["n_az"] is None:
        metadata["n_az"] = n_az
    if metadata["n_el"] is None:
        metadata["n_el"] = n_el

    if metadata["density"] is None:
        raise ValueError(
            f"Could not parse density from header in {filepath}. "
            "Expected a line like:  # Density: 2.65 g/cm3"
        )

    return az_c, el_c, T_2d, metadata


def _parse_header_line(line, metadata):
    """Extract metadata from a comment line in-place."""
    line_low = line.lower()

    # Density: 2.65 g/cm3
    if "density:" in line_low:
        parts = line.split(":")
        if len(parts) >= 2:
            try:
                metadata["density"] = float(parts[1].split()[0])
            except (ValueError, IndexError):
                pass

    # Detector: lat=40.827100  lon=14.400600  alt=608.0 m
    if "detector:" in line_low:
        import re
        m_lat = re.search(r"lat\s*=\s*([-\d.]+)", line, re.IGNORECASE)
        m_lon = re.search(r"lon\s*=\s*([-\d.]+)", line, re.IGNORECASE)
        m_alt = re.search(r"alt\s*=\s*([\d.]+)", line, re.IGNORECASE)
        if m_lat:
            try:
                metadata["det_lat"] = float(m_lat.group(1))
            except ValueError:
                pass
        if m_lon:
            try:
                metadata["det_lon"] = float(m_lon.group(1))
            except ValueError:
                pass
        if m_alt:
            try:
                metadata["det_alt"] = float(m_alt.group(1))
            except ValueError:
                pass

    # Grid: 72 az bins  x  30 el bins
    if "grid:" in line_low and "az bins" in line_low:
        import re
        m = re.search(r"(\d+)\s+az\s+bins\s+x\s+(\d+)\s+el\s+bins", line, re.IGNORECASE)
        if m:
            metadata["n_az"] = int(m.group(1))
            metadata["n_el"] = int(m.group(2))


def build_tsim_library(file_list):
    """
    Load multiple transmission maps and return a density-keyed dict.

    Returns dict {rho_float: T_2d_array}, sorted ascending by density.
    All maps must share the same az/el grid.
    """
    if not file_list:
        raise ValueError("file_list is empty — supply at least one T_sim file.")

    lib = {}
    ref_az = None
    ref_el = None

    for fpath in file_list:
        fpath = str(fpath).strip()
        if not fpath:
            continue
        az_c, el_c, T_2d, meta = load_transmission_map(fpath)
        rho = float(meta["density"])

        if ref_az is None:
            ref_az = az_c
            ref_el = el_c
        else:
            if not (np.allclose(ref_az, az_c, atol=1e-4) and
                    np.allclose(ref_el, el_c, atol=1e-4)):
                raise ValueError(
                    f"Grid mismatch: file '{fpath}' has az grid "
                    f"{az_c[:3]}... el grid {el_c[:3]}... "
                    f"but reference grid is az {ref_az[:3]}... el {ref_el[:3]}..."
                )

        if rho in lib:
            import warnings
            warnings.warn(
                f"Duplicate density {rho} g/cm³ in library — "
                f"overwriting with {fpath}", stacklevel=2
            )
        lib[rho] = T_2d

    if not lib:
        raise ValueError("No valid files were loaded into the library.")

    # Return sorted by ascending density
    return dict(sorted(lib.items()))


# ─────────────────────────────────────────────────────────────────────────────
# Core inversion
# ─────────────────────────────────────────────────────────────────────────────

def invert_density_map(T_data, tsim_lib, sigma_T=None, min_sensitivity=0.005):
    """
    Recover density ρ̂(az, el) by interpolating T_sim(ρ) = T_data pixel-by-pixel.

    Returns (rho_map, sigma_rho, status_map), all shape (n_az, n_el).

    Status codes:
      0 = OK
      1 = open sky (T_sim ≈ 1 for all densities)
      2 = T_data above library range (ρ < ρ_min, under-dense)
      3 = T_data below library range (ρ > ρ_max, over-dense)
      4 = low sensitivity (|dT/dρ| < min_sensitivity)
    """
    if not tsim_lib:
        raise ValueError("tsim_lib is empty.")

    rho_vec = np.array(sorted(tsim_lib.keys()), dtype=np.float64)
    n_rho   = len(rho_vec)

    # Stack T_sim arrays: shape (n_rho, n_az, n_el)
    T_stack = np.stack([tsim_lib[r] for r in rho_vec], axis=0)

    n_az, n_el = T_data.shape
    rho_map    = np.full((n_az, n_el), np.nan, dtype=np.float64)
    sigma_rho  = np.full((n_az, n_el), np.nan, dtype=np.float64)
    status_map = np.zeros((n_az, n_el), dtype=np.int8)

    for ia in range(n_az):
        for ie in range(n_el):
            td = float(T_data[ia, ie])

            if np.isnan(td):
                status_map[ia, ie] = 1
                continue

            T_vec = T_stack[:, ia, ie]   # T_sim at each ρ, shape (n_rho,)

            # Check for open sky: T_sim ≈ 1 for all densities
            if np.all(T_vec >= _OPENSKY_THRESHOLD):
                rho_map[ia, ie]    = np.nan
                status_map[ia, ie] = 1
                continue

            T_min = float(np.min(T_vec))
            T_max = float(np.max(T_vec))

            # T_sim decreases with ρ, so T_max is at ρ_min and T_min is at ρ_max.
            # Check range
            if td > T_max:
                # T_data higher than the highest T_sim → density below ρ_min
                rho_map[ia, ie]    = rho_vec[0]   # clamp to minimum
                status_map[ia, ie] = 2
            elif td < T_min:
                # T_data lower than the lowest T_sim → density above ρ_max
                rho_map[ia, ie]    = rho_vec[-1]  # clamp to maximum
                status_map[ia, ie] = 3
            else:
                # Interpolate: T_sim is decreasing in ρ, so reverse for np.interp
                # np.interp requires xp to be increasing
                rho_hat = np.interp(td, T_vec[::-1], rho_vec[::-1])
                rho_map[ia, ie] = rho_hat

                # Estimate gradient dT/dρ at ρ̂ using finite differences on T_vec
                grad_T = np.gradient(T_vec, rho_vec)   # (n_rho,)
                dT_drho = float(np.interp(rho_hat, rho_vec, grad_T))

                if abs(dT_drho) < min_sensitivity:
                    status_map[ia, ie] = 4
                else:
                    status_map[ia, ie] = 0

                # Propagate uncertainty
                if sigma_T is not None:
                    sT = float(sigma_T[ia, ie])
                    if abs(dT_drho) > _GRADIENT_THRESHOLD and sT > 0.0:
                        sigma_rho[ia, ie] = sT / abs(dT_drho)

    return rho_map, sigma_rho, status_map


# ─────────────────────────────────────────────────────────────────────────────
# Double ratio
# ─────────────────────────────────────────────────────────────────────────────

def compute_double_ratio(T_data, T_sim_ref):
    """
    Compute D = T_data / T_sim_ref (element-wise, clipped to [0, 5]).

    D = 1 → density matches reference; D > 1 → lower density; D < 1 → higher density.
    """
    T_data   = np.asarray(T_data,   dtype=np.float64)
    T_sim_ref = np.asarray(T_sim_ref, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        D = np.where(T_sim_ref > 0.0, T_data / T_sim_ref, np.nan)

    return np.clip(D, 0.0, 5.0)


# ─────────────────────────────────────────────────────────────────────────────
# Chi-squared landscape (single pixel diagnostic)
# ─────────────────────────────────────────────────────────────────────────────

def compute_chi2_landscape(T_data, tsim_lib, sigma_T, az_idx, el_idx,
                            n_rho_fine=200):
    """
    Compute χ²(ρ) curve for a single pixel — useful for diagnostic visualisation.

    Returns (rho_fine, chi2_curve).
    """
    if not tsim_lib:
        raise ValueError("tsim_lib is empty.")

    rho_vec  = np.array(sorted(tsim_lib.keys()), dtype=np.float64)
    T_stack  = np.stack([tsim_lib[r] for r in rho_vec], axis=0)

    td  = float(T_data[az_idx, el_idx])
    sT  = float(sigma_T[az_idx, el_idx])

    if sT <= 0.0 or np.isnan(td) or np.isnan(sT):
        raise ValueError(
            f"Pixel ({az_idx}, {el_idx}): T_data={td:.4f}, sigma_T={sT:.4g} — "
            "cannot compute chi2 (zero or NaN uncertainty)."
        )

    T_pixel = T_stack[:, az_idx, el_idx]   # shape (n_rho,)

    rho_fine = np.linspace(rho_vec[0], rho_vec[-1], n_rho_fine)

    if _SCIPY_OK and len(rho_vec) >= 4:
        _fn      = _scipy_interp1d(rho_vec, T_pixel, kind="cubic",
                                   bounds_error=False, fill_value="extrapolate")
        T_interp = np.clip(_fn(rho_fine), 0.0, 1.0)
    else:
        T_interp = np.interp(rho_fine, rho_vec, T_pixel)

    chi2_curve = (td - T_interp) ** 2 / sT ** 2
    return rho_fine, chi2_curve


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic T_data generator
# ─────────────────────────────────────────────────────────────────────────────

def generate_synthetic_tdata(tsim_lib, true_rho, n_events=10_000, seed=42):
    """
    Generate a synthetic T_data map with Poisson noise for testing.

    true_rho: scalar (uniform density) or (n_az, n_el) array.
    Returns (T_data, sigma_T) both shape (n_az, n_el).
    """
    if not tsim_lib:
        raise ValueError("tsim_lib is empty.")

    rho_vec = np.array(sorted(tsim_lib.keys()), dtype=np.float64)
    T_stack = np.stack([tsim_lib[r] for r in rho_vec], axis=0)

    # Determine shape from the library
    _, n_az, n_el = T_stack.shape

    # Broadcast true_rho to 2D
    true_rho = np.asarray(true_rho, dtype=np.float64)
    if true_rho.ndim == 0:
        rho_map = np.full((n_az, n_el), float(true_rho))
    elif true_rho.shape == (n_az, n_el):
        rho_map = true_rho
    else:
        raise ValueError(
            f"true_rho shape {true_rho.shape} is incompatible with "
            f"library grid ({n_az}, {n_el})."
        )

    rng     = np.random.default_rng(seed)
    T_data  = np.zeros((n_az, n_el), dtype=np.float64)
    sigma_T = np.zeros((n_az, n_el), dtype=np.float64)

    for ia in range(n_az):
        for ie in range(n_el):
            rho_true_pix = float(rho_map[ia, ie])
            T_vec        = T_stack[:, ia, ie]

            # Interpolate true transmission
            T_true = float(np.interp(rho_true_pix, rho_vec, T_vec))
            T_true = float(np.clip(T_true, 0.0, 1.0))

            if T_true >= _OPENSKY_THRESHOLD:
                # Open sky
                T_data[ia, ie]  = 1.0
                sigma_T[ia, ie] = 1.0 / np.sqrt(max(n_events, 1))
            else:
                # Poisson sampling
                N_rock = int(rng.poisson(n_events * T_true))
                T_data[ia, ie]  = N_rock / n_events
                sigma_T[ia, ie] = float(np.sqrt(
                    T_true * max(1.0 - T_true, 0.0) / max(n_events, 1)
                ))

    return T_data, sigma_T


# ─────────────────────────────────────────────────────────────────────────────
# Output writer
# ─────────────────────────────────────────────────────────────────────────────

def write_density_map(az_c, el_c, rho_map, sigma_rho, status_map, fpath,
                       metadata=None):
    """
    Write the inverted density map to a text file.

    Columns: azimuth[deg]  elevation[deg]  density[g/cm3]  sigma_rho[g/cm3]  status
    """
    fpath = Path(fpath)
    fpath.parent.mkdir(parents=True, exist_ok=True)

    metadata = metadata or {}
    with open(fpath, "w") as fh:
        fh.write("# UCMuon Density Analysis — Inverted Density Map\n")

        if metadata.get("det_lat") is not None:
            fh.write(
                f"# Detector: lat={metadata['det_lat']:.6f}"
                f"  lon={metadata.get('det_lon', 0):.6f}"
                f"  alt={metadata.get('det_alt', 0):.1f} m\n"
            )
        if metadata.get("rho_ref") is not None:
            fh.write(f"# Reference density: {metadata['rho_ref']:.3f} g/cm3\n")
        if metadata.get("n_events") is not None:
            fh.write(f"# N_events (synthetic): {metadata['n_events']}\n")
        if metadata.get("extra_comment"):
            fh.write(f"# {metadata['extra_comment']}\n")

        fh.write(f"# Grid: {len(az_c)} az bins  x  {len(el_c)} el bins\n")
        fh.write(
            "# Status codes: 0=OK 1=open_sky 2=below_rho_min "
            "3=above_rho_max 4=low_sensitivity\n"
        )
        fh.write(
            "# Cols: azimuth[deg]  elevation[deg]  density[g/cm3]"
            "  sigma_rho[g/cm3]  status\n"
        )

        for ia, az in enumerate(az_c):
            for ie, el in enumerate(el_c):
                rho_val = rho_map[ia, ie]
                sig_val = sigma_rho[ia, ie]
                sta_val = int(status_map[ia, ie])
                rho_str = f"{rho_val:.4f}" if not np.isnan(rho_val) else "nan"
                sig_str = f"{sig_val:.4f}" if not np.isnan(sig_val) else "nan"
                fh.write(
                    f"{az:8.2f}  {el:7.2f}  {rho_str:>10}  {sig_str:>10}  {sta_val}\n"
                )


# ═════════════════════════════════════════════════════════════════════════════
#  DIRECT OPACITY INVERSION  (analytical, single open-sky + target)
#
#  Physics
#  ───────
#  The transmission T(az,el) = Φ_target / Φ_open depends on the rock only through
#  the opacity  ϱ = ∫ρ·dl  [g/cm²] and the zenith angle θ:
#
#       T(ϱ, θ) = I(>E_min(ϱ), θ) / I(>0, θ)
#
#  where I(>E, θ) is the integrated muon flux above kinetic energy E and
#  E_min(ϱ) is the minimum energy needed to traverse opacity ϱ (CSDA range).
#  T is strictly monotonically *decreasing* in ϱ, so for each angle the curve
#  T(ϱ) can be inverted to give a unique ϱ̂ for every measured T.  No density
#  library is required.  Mean density follows from the geometry:  ρ̄ = ϱ̂ / L.
# ═════════════════════════════════════════════════════════════════════════════

# Opacity sampling grid [g/cm²] — 0 (open sky) plus log-spaced up to ~CSDA max.
_OPACITY_GRID = np.concatenate((
    [0.0],
    np.logspace(np.log10(1.0), np.log10(2.7e5), 512),
))


def _load_flux_engine():
    """Lazy-load the sibling fast_flux_estimator.py module."""
    import importlib.util
    here = Path(__file__).resolve().parent
    fpath = here / "fast_flux_estimator.py"
    if not fpath.exists():
        raise FileNotFoundError(
            f"fast_flux_estimator.py not found next to ucmuon_density_analysis.py "
            f"({here}) — required for the direct opacity inversion."
        )
    spec = importlib.util.spec_from_file_location("fast_flux_estimator", str(fpath))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def flux_model_labels():
    """Return {model_key: human label} for the GUI dropdown."""
    try:
        return dict(_load_flux_engine().MODEL_LABELS)
    except Exception:
        return {"reyna_bugaev": "Reyna–Bugaev (2006)"}


def transmission_vs_opacity(theta_deg, model="reyna_bugaev", altitude_m=0.0,
                            opacity_grid=None):
    """
    Tabulate T(ϱ) = I(>E_min(ϱ), θ) / I(>0, θ) for a single zenith angle.

    Returns (opacity_grid, T_curve), with T_curve monotonically decreasing
    from ~1 (ϱ=0) toward 0 (deep rock).

    The differential flux dΦ/dT is evaluated once for this angle and integrated
    cumulatively, so the whole curve costs a single flux evaluation.
    """
    ff = _load_flux_engine()
    if opacity_grid is None:
        opacity_grid = _OPACITY_GRID

    T_energy = ff._T_GRID                                  # energy grid [GeV]
    phi      = ff.differential_flux(T_energy, theta_deg, model, altitude_m)

    # Cumulative integral I(>E): area from each node to the top of the grid.
    dE  = np.diff(T_energy)
    seg = 0.5 * (phi[1:] + phi[:-1]) * dE                 # trapezoid segments
    rev = np.cumsum(seg[::-1])[::-1]                       # integral node→top
    I_above = np.concatenate((rev, [0.0]))                # len == len(T_energy)
    I_open  = float(I_above[0])

    if I_open <= 0.0:
        return opacity_grid, np.full(len(opacity_grid), np.nan)

    T_curve = np.empty(len(opacity_grid))
    for i, x in enumerate(opacity_grid):
        emin = ff.emin_from_opacity(float(x))
        if emin is None:                                  # beyond CSDA table
            T_curve[i] = 0.0
        else:
            emin_c     = max(float(emin), float(T_energy[0]))
            I_x        = float(np.interp(emin_c, T_energy, I_above))
            T_curve[i] = I_x / I_open
    return opacity_grid, np.clip(T_curve, 0.0, 1.0)


def invert_opacity_map(T_data, el_c, model="reyna_bugaev", altitude_m=0.0,
                       sigma_T=None, opacity_grid=None, opensky_threshold=0.999):
    """
    Invert a measured transmission map into an opacity map ϱ̂(az, el) [g/cm²].

    Parameters
    ----------
    T_data : (n_az, n_el) array     measured transmission T = Φ_target / Φ_open
    el_c   : (n_el,) array          elevation bin centres [deg]; zenith = 90 − el
    model  : str                    sea-level flux model key (fast_flux_estimator)
    altitude_m : float              detector altitude [m a.s.l.]
    sigma_T : (n_az, n_el) | None   1-σ uncertainty on T (for error propagation)

    Returns
    -------
    opacity_map  : (n_az, n_el)  ϱ̂ [g/cm²]   (NaN where open sky)
    sigma_opacity: (n_az, n_el)  σ_ϱ [g/cm²]  (NaN if sigma_T not given)
    status_map   : (n_az, n_el) int8
        0 = OK
        1 = open sky (T ≳ 1, no rock)
        3 = beyond range (T below deepest tabulated → ϱ clamped to max)

    Notes
    -----
    The flux models here are azimuth-symmetric, so the T(ϱ) curve depends only
    on elevation; it is tabulated once per elevation bin and reused across all
    azimuths for speed.
    """
    T_data = np.asarray(T_data, dtype=np.float64)
    n_az, n_el = T_data.shape
    if opacity_grid is None:
        opacity_grid = _OPACITY_GRID

    opacity_map   = np.full((n_az, n_el), np.nan, dtype=np.float64)
    sigma_opacity = np.full((n_az, n_el), np.nan, dtype=np.float64)
    status_map    = np.zeros((n_az, n_el), dtype=np.int8)

    for ie in range(n_el):
        theta = 90.0 - float(el_c[ie])
        _, T_curve = transmission_vs_opacity(theta, model, altitude_m, opacity_grid)

        if np.all(np.isnan(T_curve)):
            status_map[:, ie] = 1
            continue

        # Build increasing-T arrays for np.interp (T decreases with ϱ).
        T_inc   = T_curve[::-1]
        x_inc   = opacity_grid[::-1]
        # dT/dϱ for uncertainty propagation.
        with np.errstate(divide="ignore", invalid="ignore"):
            dTdx = np.gradient(T_curve, opacity_grid)

        T_top = float(T_curve[0])      # transmission at ϱ=0  (≈1)
        T_bot = float(T_curve[-1])     # transmission at deepest ϱ

        for ia in range(n_az):
            td = float(T_data[ia, ie])
            if np.isnan(td):
                status_map[ia, ie] = 1
                continue

            if td >= min(T_top, opensky_threshold):
                # No measurable attenuation → open sky.
                status_map[ia, ie] = 1
                continue

            if td <= T_bot:
                opacity_map[ia, ie] = float(opacity_grid[-1])
                status_map[ia, ie]  = 3
                continue

            x_hat = float(np.interp(td, T_inc, x_inc))
            opacity_map[ia, ie] = x_hat
            status_map[ia, ie]  = 0

            if sigma_T is not None:
                sT = float(sigma_T[ia, ie])
                g  = float(np.interp(x_hat, opacity_grid, dTdx))
                if sT > 0.0 and abs(g) > 0.0:
                    sigma_opacity[ia, ie] = sT / abs(g)

    return opacity_map, sigma_opacity, status_map


def opacity_to_density(opacity_map, L_map, sigma_opacity=None):
    """
    Convert opacity ϱ [g/cm²] and line-of-sight path length L [m] into mean
    density ρ̄ = ϱ / (100·L)  [g/cm³].   (100 converts L from m to cm.)

    Pixels with L ≤ 0 (open sky / no rock) become NaN.
    Returns (rho_map, sigma_rho).
    """
    opacity_map = np.asarray(opacity_map, dtype=np.float64)
    L_map       = np.asarray(L_map,       dtype=np.float64)
    if L_map.shape != opacity_map.shape:
        raise ValueError(
            f"Path-length grid {L_map.shape} does not match opacity grid "
            f"{opacity_map.shape}."
        )

    rho_map   = np.full_like(opacity_map, np.nan)
    sigma_rho = np.full_like(opacity_map, np.nan)

    valid = np.isfinite(opacity_map) & np.isfinite(L_map) & (L_map > 0.0)
    rho_map[valid] = opacity_map[valid] / (100.0 * L_map[valid])

    if sigma_opacity is not None:
        sigma_opacity = np.asarray(sigma_opacity, dtype=np.float64)
        vsig = valid & np.isfinite(sigma_opacity)
        sigma_rho[vsig] = sigma_opacity[vsig] / (100.0 * L_map[vsig])

    return rho_map, sigma_rho


def path_length_from_tsim(T_sim, rho_sim, el_c, model="reyna_bugaev",
                          altitude_m=0.0, opacity_grid=None):
    """
    Recover a line-of-sight path length map L(az, el) [m] from a single
    simulated transmission map at known density ρ_sim.

    The forward sim encodes opacity ϱ_sim = ρ_sim·L; inverting T_sim with the
    analytic flux model recovers ϱ_sim, hence L = ϱ_sim / (100·ρ_sim).

    Open-sky pixels are returned as L = 0.
    """
    if rho_sim <= 0.0:
        raise ValueError(f"rho_sim must be positive, got {rho_sim}.")
    opac_sim, _, status = invert_opacity_map(
        T_sim, el_c, model=model, altitude_m=altitude_m,
        opacity_grid=opacity_grid,
    )
    L_map = np.where(np.isfinite(opac_sim), opac_sim / (100.0 * rho_sim), 0.0)
    L_map[status == 1] = 0.0
    return L_map


# ─────────────────────────────────────────────────────────────────────────────
#  Generic map I/O for measured data (transmission / flux / overburden)
# ─────────────────────────────────────────────────────────────────────────────

def load_generic_map(filepath, angle_is_zenith=False, value_col=2):
    """
    Load a 3+ column az / angle / value map without requiring a density header.

    Parameters
    ----------
    filepath : str
    angle_is_zenith : bool
        If True the second column is zenith [deg] and is converted to
        elevation (el = 90 − zenith); rows are reordered so el ascends.
        If False (default) the second column is already elevation.
    value_col : int
        Zero-based index of the data column (default 2 = third column).

    Returns
    -------
    az_c  : (n_az,) elevation-sorted azimuth centres [deg]
    el_c  : (n_el,) ascending elevation centres [deg]
    V_2d  : (n_az, n_el) values (NaN where missing)
    meta  : dict (may contain 'density', 'det_lat', 'det_lon', 'det_alt')
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise ValueError(f"Map file not found: {filepath}")

    meta = {"density": None, "det_lat": None, "det_lon": None,
            "det_alt": None, "filepath": str(filepath)}

    rows = []
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                _parse_header_line(s, meta)
                continue
            parts = s.split()
            if len(parts) <= value_col:
                continue
            try:
                a   = float(parts[0])
                ang = float(parts[1])
                v   = float(parts[value_col])
            except ValueError:
                continue
            rows.append((a, ang, v))

    if not rows:
        raise ValueError(f"No data rows found in {filepath}")

    data = np.array(rows, dtype=np.float64)
    az_vals  = data[:, 0]
    ang_vals = data[:, 1]
    v_vals   = data[:, 2]

    el_vals = (90.0 - ang_vals) if angle_is_zenith else ang_vals

    az_c = np.unique(az_vals)
    el_c = np.unique(el_vals)
    az_idx = {round(float(v), 6): i for i, v in enumerate(az_c)}
    el_idx = {round(float(v), 6): i for i, v in enumerate(el_c)}

    V_2d = np.full((len(az_c), len(el_c)), np.nan, dtype=np.float64)
    for a, el, v in zip(az_vals, el_vals, v_vals):
        ia = az_idx.get(round(float(a), 6))
        ie = el_idx.get(round(float(el), 6))
        if ia is not None and ie is not None:
            V_2d[ia, ie] = v

    return az_c, el_c, V_2d, meta


def transmission_from_flux_maps(flux_target, flux_open, clip=True):
    """
    Form a transmission map T = Φ_target / Φ_open from two raw flux maps.

    Open-sky pixels with Φ_open ≤ 0 become NaN.  When clip is True the result
    is clipped to [0, 1] (transmission is physically bounded).
    """
    flux_target = np.asarray(flux_target, dtype=np.float64)
    flux_open   = np.asarray(flux_open,   dtype=np.float64)
    if flux_target.shape != flux_open.shape:
        raise ValueError(
            f"Flux map shapes differ: target {flux_target.shape} vs "
            f"open-sky {flux_open.shape}."
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(flux_open > 0.0, flux_target / flux_open, np.nan)
    return np.clip(T, 0.0, 1.0) if clip else T


def load_overburden_as_L(filepath):
    """
    Load a terrain-engine overburden map and convert it to a path-length map.

    The overburden file stores columns  az / zenith / overburden[g/cm²] / open_sky
    and a  '# Density: <ρ> g/cm3'  header.  Path length L = ϱ / (100·ρ)  [m].

    Returns (az_c, el_c, L_map, meta).  Open-sky pixels → L = 0.
    """
    az_c, el_c, ob_2d, meta = load_generic_map(
        filepath, angle_is_zenith=True, value_col=2
    )
    rho = meta.get("density")
    if rho is None or rho <= 0.0:
        raise ValueError(
            f"Overburden file {filepath} has no usable '# Density:' header — "
            "cannot convert opacity to path length."
        )
    L_map = np.where(np.isfinite(ob_2d), ob_2d / (100.0 * float(rho)), 0.0)
    return az_c, el_c, L_map, meta


def write_opacity_map(az_c, el_c, opacity_map, sigma_opacity, status_map,
                      fpath, rho_map=None, sigma_rho=None, metadata=None):
    """
    Write the direct-inversion result (opacity, optional density) to a .dat file.

    Columns: azimuth[deg] elevation[deg] opacity[g/cm2] sigma_opacity[g/cm2]
             status [density[g/cm3] sigma_rho[g/cm3]]
    """
    fpath = Path(fpath)
    fpath.parent.mkdir(parents=True, exist_ok=True)
    metadata = metadata or {}
    has_rho  = rho_map is not None

    with open(fpath, "w") as fh:
        fh.write("# UCMuon Density Analysis — Direct Opacity Inversion\n")
        if metadata.get("model"):
            fh.write(f"# Flux model: {metadata['model']}\n")
        if metadata.get("altitude_m") is not None:
            fh.write(f"# Altitude: {metadata['altitude_m']:.1f} m\n")
        if metadata.get("L_source"):
            fh.write(f"# Path length source: {metadata['L_source']}\n")
        fh.write(f"# Grid: {len(az_c)} az bins  x  {len(el_c)} el bins\n")
        fh.write("# Status codes: 0=OK 1=open_sky 3=beyond_range\n")
        cols = "# Cols: azimuth[deg]  elevation[deg]  opacity[g/cm2]  sigma_opacity[g/cm2]  status"
        if has_rho:
            cols += "  density[g/cm3]  sigma_rho[g/cm3]"
        fh.write(cols + "\n")

        for ia, az in enumerate(az_c):
            for ie, el in enumerate(el_c):
                ov = opacity_map[ia, ie]
                sv = sigma_opacity[ia, ie]
                ov_s = f"{ov:.3f}" if np.isfinite(ov) else "nan"
                sv_s = f"{sv:.3f}" if np.isfinite(sv) else "nan"
                line = (f"{az:8.2f}  {el:7.2f}  {ov_s:>12}  {sv_s:>12}  "
                        f"{int(status_map[ia, ie])}")
                if has_rho:
                    rv = rho_map[ia, ie]
                    rs = sigma_rho[ia, ie] if sigma_rho is not None else np.nan
                    rv_s = f"{rv:.4f}" if np.isfinite(rv) else "nan"
                    rs_s = f"{rs:.4f}" if np.isfinite(rs) else "nan"
                    line += f"  {rv_s:>10}  {rs_s:>10}"
                fh.write(line + "\n")
