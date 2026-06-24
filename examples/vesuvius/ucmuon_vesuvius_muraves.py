#!/usr/bin/env python3
"""
ucmuon_vesuvius_muraves.py  —  UCLouvain Muography Group
─────────────────────────────────────────────────────────────────────────────
Reproduces the four key MURAVES simulation outputs from the Muographers 2026
presentation (Rajan et al., Budapest, 4 Jun 2026):

  Fig 1  — Rock-thickness map   (az × el heatmap)         ↔ slide 9
  Fig 2  — Flux vs elevation at az=summit                 ↔ slides 11, 20
  Fig 3  — Free-sky vs through-rock flux maps             ↔ slide 12
  Fig 4  — Transmission ratio T_sim = Φ_rock / Φ_sky     ↔ slide 14

DEM handling
────────────
If a DEM file is passed on the command line the full UCMuon terrain engine is
used (requires rasterio).  Without a DEM, a smooth truncated-cone model of
Vesuvius is used for the geometry — good enough to reproduce the slide shapes
and all physics comparisons.

Usage:
    python ucmuon_vesuvius_muraves.py                      # cone model
    python ucmuon_vesuvius_muraves.py vesuvius_dem.tif     # real DEM
    python ucmuon_vesuvius_muraves.py vesuvius_ingv.xyz    # INGV 5 m DEM

Output:  four PNG figures saved to the current directory.

References:
  Tioukov et al. (2019) Sci. Rep. 9, 6695
  Lo Bue   et al. (2023) J. Geophys. Res. 128, e2022JB025446
  Rajan    et al. (2026) MURAVES/Muographers26 presentation
"""

from __future__ import annotations
import sys
import os
import importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Paths: gui/ must be on sys.path to import terrain driver and backward MC
# ─────────────────────────────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent
_GUI  = _ROOT / "gui"
if str(_GUI) not in sys.path:
    sys.path.insert(0, str(_GUI))


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_bmc = _load_module("ucmuon_backward_mc",   _GUI / "ucmuon_backward_mc.py")
_drv = _load_module("ucmuon_terrain_driver", _GUI / "ucmuon_terrain_driver.py")

# ─────────────────────────────────────────────────────────────────────────────
# MURAVES geometry (Rajan 2026 / Tioukov 2019)
# ─────────────────────────────────────────────────────────────────────────────
DET_LAT  = 40.8271      # MURAVES detector, SW flank
DET_LON  = 14.4006
DET_ALT  = 608.0        # m a.s.l.

SUM_LAT  = 40.8218      # Gran Cono summit
SUM_LON  = 14.4265
SUM_ALT  = 1281.0       # m a.s.l.

# Summit in detector-centric East-North-Up [m]
_M_PER_DEG_LAT = 111_320.0
_SUM_E = (SUM_LON - DET_LON) * _M_PER_DEG_LAT * np.cos(np.radians(DET_LAT))
_SUM_N = (SUM_LAT - DET_LAT) * _M_PER_DEG_LAT
_SUM_U = SUM_ALT - DET_ALT

# Geographic azimuth to summit (0=N, 90=E) and elevation
AZ_TO_SUMMIT = float(np.degrees(np.arctan2(_SUM_E, _SUM_N))) % 360.0
EL_TO_SUMMIT = float(np.degrees(np.arctan2(_SUM_U, np.sqrt(_SUM_E**2 + _SUM_N**2))))

# Cone model parameters (used when no DEM is supplied)
_CONE_BASE_ALT  = 150.0     # m a.s.l. — cone base (lava apron boundary)
_CONE_BASE_RAD  = 1900.0    # m — cone base radius at _CONE_BASE_ALT
_CONE_APEX      = (_SUM_E, _SUM_N, _SUM_U)   # apex in ENU relative to detector

# Angular grid (1° bins over the MURAVES-relevant range)
N_AZ     = 360    # 1° resolution
N_ZE     = 85     # 1° steps, ze = 0.5° … 84.5°
ZE_MAX   = 85.0

# Simulation parameters
RHO_REF  = 2.65   # g/cm³ — reference density used in MURAVES T_sim slide
SPEC_MODE = 3     # Guan et al. 2015 — used by MURAVES
E_MIN    = 1.0    # GeV
E_MAX    = 2500.0 # GeV
N_E      = 50     # energy integration points

# ─────────────────────────────────────────────────────────────────────────────
# Synthetic cone overburden (no DEM)
# ─────────────────────────────────────────────────────────────────────────────

def cone_overburden(az_c: np.ndarray, ze_c: np.ndarray,
                    rho: float, step_m: float = 25.0) -> np.ndarray:
    """
    Compute overburden [g/cm²] through a truncated-cone Vesuvius model.

    The cone:
      apex  = summit (ENU relative to detector)
      base  = circle of radius _CONE_BASE_RAD at altitude _CONE_BASE_ALT
      slope = linear taper from apex down to base

    A ray from the detector in direction (az, ze) accumulates path length
    wherever it flies below the cone surface.
    """
    E_s, N_s, U_s = _CONE_APEX
    h_sum   = DET_ALT + U_s        # absolute summit altitude [m]
    h_base  = _CONE_BASE_ALT       # absolute base altitude [m]
    R_base  = _CONE_BASE_RAD
    max_d   = 10_000.0             # maximum ray length [m]

    n_az, n_ze = len(az_c), len(ze_c)
    ob = np.zeros((n_az, n_ze), dtype=np.float64)

    for ia, az in enumerate(az_c):
        az_r = np.radians(az)
        sin_az, cos_az = np.sin(az_r), np.cos(az_r)
        for iz, ze in enumerate(ze_c):
            ze_r   = np.radians(ze)
            sin_ze = np.sin(ze_r)
            dE     =  sin_ze * sin_az   # East unit vector
            dN     =  sin_ze * cos_az   # North unit vector
            dU     =  np.cos(ze_r)      # Up unit vector

            underground_m = 0.0
            dist = step_m / 2.0         # start at half-step to avoid origin issues
            while dist < max_d:
                E_pos = dE * dist       # position in ENU
                N_pos = dN * dist
                U_pos = dU * dist
                alt   = DET_ALT + U_pos

                if alt < h_base or alt > h_sum:
                    dist += step_m
                    continue

                # Cone radius at this altitude (linear taper)
                frac      = (h_sum - alt) / (h_sum - h_base)
                r_cone    = R_base * frac
                r_horiz   = np.sqrt((E_pos - E_s)**2 + (N_pos - N_s)**2)

                if r_horiz < r_cone:
                    underground_m += step_m

                dist += step_m

            ob[ia, iz] = underground_m * 100.0 * rho   # [g/cm²]

    return ob


# ─────────────────────────────────────────────────────────────────────────────
# Flux integration (delegates to backward MC physics module)
# ─────────────────────────────────────────────────────────────────────────────

def _opensky_flux(ze_c: np.ndarray, spec_mode: int,
                  E_min: float, E_max: float, n_E: int) -> np.ndarray:
    """Integrated open-sky flux [m⁻² s⁻¹ sr⁻¹] for each zenith angle."""
    E_grid = np.logspace(np.log10(E_min), np.log10(E_max), n_E)
    dln_E  = (np.log(E_max) - np.log(E_min)) / n_E
    result = np.zeros(len(ze_c))
    for iz, ze in enumerate(ze_c):
        result[iz] = float(np.sum(
            _bmc._flux_surface(E_grid, np.radians(ze), spec_mode)
        )) * dln_E
    return result


def _rock_flux(depth_m: float, rho: float, ze: float,
               spec_mode: int, E_min: float, E_max: float, n_E: int) -> float:
    """Flux after traversing depth_m of rock at zenith ze."""
    res = _bmc.backward_mc_flux(
        depth_m       = depth_m,
        rho           = rho,
        mat_id        = 1,
        spectrum_mode = spec_mode,
        E_min_GeV     = E_min,
        E_max_GeV     = E_max,
        theta_max_deg = ze,
        n_E           = n_E,
        n_theta       = 1,
        mode          = 1,
    )
    return float(res["rate_m2_s"])


def flux_maps(az_c: np.ndarray, ze_c: np.ndarray, overburden: np.ndarray,
              rho: float, sky: np.ndarray,
              spec_mode: int = SPEC_MODE,
              E_min: float = E_MIN, E_max: float = E_MAX,
              n_E: int = N_E) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute (through-rock flux map, transmission map) given a pre-computed
    open-sky flux array `sky` (shape n_ze).
    Returns: (flux_map (n_az, n_ze), T_sim (n_az, n_ze)).
    """
    n_az, n_ze = overburden.shape
    flux = np.zeros((n_az, n_ze))

    total = n_az * n_ze
    done  = 0
    report_every = max(1, n_az // 10)

    for ia in range(n_az):
        for iz, ze in enumerate(ze_c):
            X = overburden[ia, iz]
            if X < 1.0:
                flux[ia, iz] = sky[iz]
            else:
                cos_ze  = max(np.cos(np.radians(ze)), 0.02)
                depth_m = max(X / (rho * 100.0 * cos_ze), 1.0)
                flux[ia, iz] = _rock_flux(depth_m, rho, ze,
                                          spec_mode, E_min, E_max, n_E)
            done += 1

        if (ia + 1) % report_every == 0:
            print(f"  Flux: {done}/{total} ({100*done/total:.0f}%)", flush=True)

    T_sim = _drv.compute_transmission_map(flux, sky)
    return flux, T_sim


def flux_slice_at_azimuth(ob_slice: np.ndarray, ze_c: np.ndarray,
                           rho: float, sky: np.ndarray,
                           spec_mode: int = SPEC_MODE,
                           E_min: float = E_MIN, E_max: float = E_MAX,
                           n_E: int = N_E) -> np.ndarray:
    """
    Compute through-rock flux for a single azimuth slice (1D, n_ze values).
    ob_slice : overburden [g/cm²] at fixed azimuth, shape (n_ze,).
    Returns flux array (n_ze,).
    """
    result = np.zeros(len(ze_c))
    for iz, ze in enumerate(ze_c):
        X = ob_slice[iz]
        if X < 1.0:
            result[iz] = sky[iz]
        else:
            cos_ze  = max(np.cos(np.radians(ze)), 0.02)
            depth_m = max(X / (rho * 100.0 * cos_ze), 1.0)
            result[iz] = _rock_flux(depth_m, rho, ze, spec_mode, E_min, E_max, n_E)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

_CMAP_THICK  = "viridis"
_CMAP_FLUX   = "YlGn"
_CMAP_TRANS  = "plasma"

COLORS_DENS  = {
    1.0:  ("#00bcd4", "ρ = 1.0 g/cm³"),
    2.0:  ("#9c27b0", "ρ = 2.0 g/cm³"),
    2.65: ("#e91e63", "ρ = 2.65 g/cm³"),
    3.0:  ("#3f51b5", "ρ = 3.0 g/cm³"),
}


def _az_el_extent(az_c, ze_c):
    """imshow extent in (azimuth, elevation) space."""
    el_c   = 90.0 - ze_c
    daz    = az_c[1] - az_c[0] if len(az_c) > 1 else 1.0
    del_   = el_c[1] - el_c[0] if len(el_c) > 1 else 1.0
    return [az_c[0]  - daz/2, az_c[-1] + daz/2,
            el_c[-1] - del_/2, el_c[0]  + del_/2]


def _summit_indicator(ax, az_to_sum: float, el_to_sum: float, label=True):
    ax.axvline(az_to_sum, color="white", lw=0.8, ls="--", alpha=0.7)
    ax.axhline(el_to_sum, color="white", lw=0.8, ls=":",  alpha=0.7)
    if label:
        ax.scatter([az_to_sum], [el_to_sum], s=50, c="white",
                   marker="*", zorder=5)
        ax.text(az_to_sum + 0.8, el_to_sum + 0.5, "summit",
                color="white", fontsize=7, va="bottom")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Rock-thickness map  (slide 9)
# ─────────────────────────────────────────────────────────────────────────────

def plot_thickness_map(az_c, ze_c, overburden, rho, outfile="fig_vesuvius_thickness.png"):
    thickness_m = overburden / (rho * 100.0)   # [m]
    el_c = 90.0 - ze_c

    # Slice at az ≈ az_to_summit
    ia_sum = int(np.argmin(np.abs(az_c - AZ_TO_SUMMIT)))
    slice_el = el_c
    slice_th = thickness_m[ia_sum, :]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                              facecolor="#0f1117")
    for ax in axes:
        ax.set_facecolor("#141620")

    # Left: 2-D thickness map
    ax = axes[0]
    extent = _az_el_extent(az_c, ze_c)
    im = ax.imshow(thickness_m.T, origin="upper", aspect="auto",
                   extent=extent, cmap=_CMAP_THICK,
                   vmin=0, vmax=np.nanmax(thickness_m))
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("Rock thickness [m]", color="white", fontsize=9)
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
    _summit_indicator(ax, AZ_TO_SUMMIT, EL_TO_SUMMIT)
    ax.set_xlabel("Azimuth [deg]", color="white")
    ax.set_ylabel("Elevation [deg]", color="white")
    ax.set_title("Rock thickness map — view from MURAVES detector", color="white")
    ax.tick_params(colors="white")
    for sp in ax.spines.values():
        sp.set_edgecolor("white")

    # Right: slice at summit azimuth
    ax2 = axes[1]
    ax2.plot(slice_el, slice_th, color="#4fc3f7", lw=2)
    ax2.axvline(EL_TO_SUMMIT, color="white", lw=0.8, ls="--", alpha=0.7,
                label=f"Summit el = {EL_TO_SUMMIT:.1f}°")
    ax2.set_xlabel("Elevation [deg]", color="white")
    ax2.set_ylabel("Rock thickness [m]", color="white")
    az_lab = f"{AZ_TO_SUMMIT:.0f}°"
    ax2.set_title(f"Rock thickness vs elevation  (az = {az_lab})", color="white")
    ax2.tick_params(colors="white")
    ax2.legend(framealpha=0.3, labelcolor="white", fontsize=8)
    ax2.grid(alpha=0.2)
    for sp in ax2.spines.values():
        sp.set_edgecolor("white")

    fig.suptitle("UCMuon / MURAVES — Mt. Vesuvius  (slide 9 equivalent)",
                 color="white", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {outfile}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Flux vs elevation for multiple densities  (slides 11, 20)
# ─────────────────────────────────────────────────────────────────────────────

def plot_flux_vs_elevation(el_c, sky_flux, rock_fluxes_by_rho,
                           outfile="fig_vesuvius_flux_elevation.png"):
    """
    rock_fluxes_by_rho : dict {rho_value: flux_1d_array (n_ze,)}
                         pre-computed 1D rock flux at summit azimuth per density.
    """
    fig, ax = plt.subplots(figsize=(9, 6), facecolor="#0f1117")
    ax.set_facecolor("#141620")

    # Open-sky reference
    ax.plot(el_c, sky_flux, color="#4caf50", lw=2, ls="--",
            label=f"Open sky ({DET_ALT:.0f} m a.s.l.)")

    # Through-rock for each density
    for rho, rock_flux in rock_fluxes_by_rho.items():
        col, lab = COLORS_DENS.get(rho, ("#aaa", f"ρ = {rho} g/cm³"))
        ax.plot(el_c, rock_flux, color=col, lw=2, label=lab)

    ax.set_xlabel("Elevation [deg]", color="white")
    ax.set_ylabel("Integrated muon flux  [m⁻² s⁻¹ sr⁻¹]", color="white")
    ax.set_title(f"Azimuth = {AZ_TO_SUMMIT:.0f}°  —  transport modes & densities\n"
                 f"(Guan spectrum, E_min = {E_MIN} GeV)", color="white")
    ax.set_yscale("log")
    ax.axvline(EL_TO_SUMMIT, color="white", lw=0.8, ls=":", alpha=0.6,
               label=f"Summit el = {EL_TO_SUMMIT:.1f}°")
    ax.set_xlim(el_c.min() - 0.5, min(el_c.max() + 0.5, 40.0))
    ax.tick_params(colors="white")
    ax.legend(framealpha=0.25, labelcolor="white", fontsize=9)
    ax.grid(alpha=0.2)
    for sp in ax.spines.values():
        sp.set_edgecolor("white")

    fig.suptitle("UCMuon / MURAVES — Flux vs elevation (slides 11, 20)",
                 color="white", fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {outfile}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Free-sky and through-rock flux maps side by side  (slide 12)
# ─────────────────────────────────────────────────────────────────────────────

def plot_flux_maps(az_c, ze_c, flux_map, sky_flux,
                   outfile="fig_vesuvius_flux_maps.png"):
    el_c       = 90.0 - ze_c
    sky_map_2d = np.tile(sky_flux, (len(az_c), 1))
    extent     = _az_el_extent(az_c, ze_c)
    v_max      = np.nanmax(sky_map_2d)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5),
                              facecolor="#0f1117", sharey=True)
    for ax in axes:
        ax.set_facecolor("#141620")

    titles = ["Free-sky muon flux", "Muon flux through Vesuvius"]
    data   = [sky_map_2d, flux_map]

    for ax, d, title in zip(axes, data, titles):
        im = ax.imshow(d.T, origin="upper", aspect="auto",
                       extent=extent, cmap=_CMAP_FLUX,
                       vmin=0, vmax=v_max)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("m⁻² sr⁻¹ s⁻¹", color="white", fontsize=8)
        cb.ax.yaxis.set_tick_params(color="white")
        plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
        _summit_indicator(ax, AZ_TO_SUMMIT, EL_TO_SUMMIT)
        ax.set_xlabel("Azimuth [deg]", color="white")
        ax.set_title(title, color="white")
        ax.tick_params(colors="white")
        for sp in ax.spines.values():
            sp.set_edgecolor("white")

    axes[0].set_ylabel("Elevation [deg]", color="white")

    note = ("Note: MURAVES azimuth convention centres the summit at 180°.  "
            f"In geographic convention (0=N) the summit is at az={AZ_TO_SUMMIT:.0f}°.")
    fig.text(0.5, -0.04, note, ha="center", color="grey", fontsize=7)
    fig.suptitle(f"UCMuon / MURAVES — Simulated muon flux  ρ = {RHO_REF} g/cm³  (slide 12)",
                 color="white", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {outfile}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Transmission ratio map  (slide 14)
# ─────────────────────────────────────────────────────────────────────────────

def plot_transmission_map(az_c, ze_c, T_sim,
                          outfile="fig_vesuvius_transmission.png"):
    el_c   = 90.0 - ze_c
    extent = _az_el_extent(az_c, ze_c)

    fig, ax = plt.subplots(figsize=(9, 5), facecolor="#0f1117")
    ax.set_facecolor("#141620")

    im = ax.imshow(T_sim.T, origin="upper", aspect="auto",
                   extent=extent, cmap=_CMAP_TRANS, vmin=0, vmax=1)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("Simulated Transmission  T_sim", color="white", fontsize=9)
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")

    _summit_indicator(ax, AZ_TO_SUMMIT, EL_TO_SUMMIT)
    ax.set_xlabel("Azimuth [deg]", color="white")
    ax.set_ylabel("Elevation [deg]", color="white")
    ax.tick_params(colors="white")
    for sp in ax.spines.values():
        sp.set_edgecolor("white")

    fig.suptitle(f"UCMuon / MURAVES — Transmission ratio  ρ_sim = {RHO_REF} g/cm³  (slide 14)",
                 color="white", fontsize=11)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {outfile}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# DEM-based overburden (delegates to terrain driver)
# ─────────────────────────────────────────────────────────────────────────────

def dem_overburden(dem_file: str, rho: float, n_az: int, n_ze: int,
                   ze_max: float, step_m: float = 25.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load DEM, ray-trace overburden.  Returns (az_c, ze_c, overburden)."""
    elev, transform = _drv.load_dem(dem_file)
    az_c, ze_c, ob, _ = _drv.compute_overburden_map(
        elev, transform,
        DET_LAT, DET_LON, DET_ALT,
        rho, n_az, n_ze, ze_max, step_m,
    )
    return az_c, ze_c, ob


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    import argparse
    import time

    p = argparse.ArgumentParser(
        description="UCMuon/MURAVES Vesuvius simulation — reproduces slides 9, 11, 12, 14")
    p.add_argument("dem", nargs="?", default=None,
                   help="DEM file (.tif / .xyz / .asc).  Omit to use synthetic cone model.")
    p.add_argument("--rho",    type=float, default=RHO_REF,
                   help=f"Rock density g/cm³ (default {RHO_REF})")
    p.add_argument("--n-az",   type=int,   default=N_AZ,
                   help=f"Azimuth bins (default {N_AZ})")
    p.add_argument("--n-ze",   type=int,   default=N_ZE,
                   help=f"Zenith bins (default {N_ZE})")
    p.add_argument("--step",   type=float, default=25.0,
                   help="Ray-trace step [m] (default 25)")
    _default_out = str(_HERE / "figs")
    p.add_argument("--out-dir", default=_default_out, help="Output directory for PNG files")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("═" * 65, flush=True)
    print("  UCMuon / MURAVES — Mt. Vesuvius simulation", flush=True)
    print("═" * 65, flush=True)
    print(f"  Detector:  lat={DET_LAT}°N  lon={DET_LON}°E  alt={DET_ALT} m", flush=True)
    print(f"  Summit:    lat={SUM_LAT}°N  lon={SUM_LON}°E  alt={SUM_ALT} m", flush=True)
    print(f"  Summit direction (geographic): az={AZ_TO_SUMMIT:.1f}°  el={EL_TO_SUMMIT:.1f}°", flush=True)
    print(f"  Density: {args.rho} g/cm³   Spectrum: Guan 2015 (mode 3)", flush=True)
    print(f"  Grid:  {args.n_az} az × {args.n_ze} ze  (ze_max={ZE_MAX}°)", flush=True)

    az_edges = np.linspace(0.0, 360.0, args.n_az + 1)
    ze_edges = np.linspace(0.0, ZE_MAX, args.n_ze + 1)
    az_c     = 0.5 * (az_edges[:-1] + az_edges[1:])
    ze_c     = 0.5 * (ze_edges[:-1] + ze_edges[1:])

    # ── 1. Overburden map ────────────────────────────────────────────────────
    t0 = time.time()
    if args.dem:
        print(f"\n  Loading DEM: {args.dem}", flush=True)
        print(f"  Ray-tracing overburden map …", flush=True)
        az_c, ze_c, ob_ref = dem_overburden(args.dem, args.rho,
                                             args.n_az, args.n_ze,
                                             ZE_MAX, args.step)
    else:
        print("\n  No DEM supplied — using synthetic truncated-cone model.", flush=True)
        print(f"  Cone apex: E={_SUM_E:.0f} m  N={_SUM_N:.0f} m  U={_SUM_U:.0f} m  (ENU)", flush=True)
        print(f"  Base radius: {_CONE_BASE_RAD} m at alt={_CONE_BASE_ALT} m a.s.l.", flush=True)
        print(f"  Computing overburden … ({args.n_az * args.n_ze} directions)", flush=True)
        ob_ref = cone_overburden(az_c, ze_c, args.rho, step_m=args.step)

    print(f"  Overburden done in {time.time()-t0:.1f} s", flush=True)

    # ── 2. Open-sky flux (computed once, reused for all densities) ───────────
    print("\n  Computing open-sky flux …", flush=True)
    sky_flux = _opensky_flux(ze_c, SPEC_MODE, E_MIN, E_MAX, N_E)

    # ── 3. Through-rock flux + transmission (reference density, full 2D) ─────
    print(f"\n  Computing flux maps (ρ = {args.rho} g/cm³) …", flush=True)
    flux_ref, T_sim_ref = flux_maps(az_c, ze_c, ob_ref, args.rho, sky_flux)

    # ── 4. 1D flux slices at summit azimuth for multi-density Fig 2 ──────────
    ia_sum = int(np.argmin(np.abs(az_c - AZ_TO_SUMMIT)))
    densities_all = sorted({1.0, 2.0, 2.65, 3.0, args.rho})
    rock_slices = {}
    print(f"\n  Computing 1D summit-azimuth flux slices for density comparison …",
          flush=True)
    for rho_i in densities_all:
        # Scale path length by density ratio (same geometry, different density)
        ob_slice = ob_ref[ia_sum, :] * (rho_i / args.rho)
        rock_slices[rho_i] = flux_slice_at_azimuth(ob_slice, ze_c, rho_i, sky_flux)
        print(f"    ρ = {rho_i} g/cm³  done", flush=True)

    el_c = 90.0 - ze_c

    # ── 5. Plots ─────────────────────────────────────────────────────────────
    print("\n  Generating figures …", flush=True)

    plot_thickness_map(az_c, ze_c, ob_ref, args.rho,
                       str(out_dir / "fig_vesuvius_thickness.png"))

    plot_flux_vs_elevation(el_c, sky_flux, rock_slices,
                           str(out_dir / "fig_vesuvius_flux_elevation.png"))

    plot_flux_maps(az_c, ze_c, flux_ref, sky_flux,
                  str(out_dir / "fig_vesuvius_flux_maps.png"))

    plot_transmission_map(az_c, ze_c, T_sim_ref,
                          str(out_dir / "fig_vesuvius_transmission.png"))

    elapsed = time.time() - t0
    print(f"\n  ═════════════════════════════════════════", flush=True)
    print(f"  Done in {elapsed:.0f} s", flush=True)
    print(f"  Figures written to: {out_dir.resolve()}", flush=True)
    print(f"  ─────────────────────────────────────────", flush=True)
    print(f"  Azimuth note: geographic convention (0=North, 90=East).", flush=True)
    print(f"  MURAVES slides centre the summit at az=180°.", flush=True)
    print(f"  Summit is at az={AZ_TO_SUMMIT:.0f}° in this output.", flush=True)
    print(f"  Rotate/remap azimuth by ({180 - AZ_TO_SUMMIT:.0f}°) for direct comparison.", flush=True)


if __name__ == "__main__":
    main()
