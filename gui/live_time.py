"""
live_time.py: one live-time formula for the GUI, the Fortran and UCMuGen
UCLouvain Muography Group | MIT License 2026

The live time of a generator run is

    T = N_tried / R,

with R the absolute rate of muons crossing the generation surface inside the
energy and zenith windows,

    R = ∫dp ∫dΩ J(p, θ) ∫_S max(0, -n·d) dA,

the quantity ``ucmugen::Generator::rate()`` estimates and
``surface_rate()`` in src/generator/ucmuon_source_module.f90 computes (the
Fortran prints it as "Surface rate R").  The same quadrature is repeated here
so that the GUI can estimate R before a run and for runs that did not print
it; tests/flux/test_live_time.py holds the three together.

With angular mode 6 the generator draws every tried muon from that flux, so a
subset of n rows (detector hits, survivors) corresponds to a rate n / T.  The
legacy modes 1-5 draw p from the vertical spectrum and θ from a fixed law
with no projection, so each row carries an importance weight

    w = [J(p,θ) max(0,-n·d) / R] / [f_v(p) g(d|p) / A],

the ratio of the true flux density to the density the generator sampled.  It
averages to 1 over the tried muons, and a subset corresponds to Σw / T.  The
densities reproduce the Fortran samplers exactly (the tabulated momentum CDF,
the 100-point Guan angular CDF), see docs/FLUX_NORMALISATION_AUDIT.md.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np

import fast_flux_estimator as _ffe

M_MU = 0.10566            # GeV, as in the Fortran and UCMuGen
PI = 3.141592654          # the Fortran's PI

# Spectrum -> dN/dT function of fast_flux_estimator (identical to the Fortran
# spectrum_intensity and ucmugen::flux::intensity, see tests/flux).
_DNDT = {1: _ffe._cosmoaleph, 4: _ffe._guan_2015, 5: _ffe._frosin_2025,
         6: _ffe._bugaev, 7: _ffe._reyna_bugaev}


@dataclass
class Source:
    """The generation surface, in the Fortran generator's conventions."""
    mode: int = 1               # 1 disk, 2 rectangle, 3 hemisphere
    plane: int = 1              # 1 XY, 2 XZ, 3 YZ (disk and rectangle)
    radius_cm: float = 100.0    # disk or hemisphere radius
    half_lx_cm: float = 100.0   # rectangle half-extents
    half_ly_cm: float = 100.0
    tilt_deg: float = 0.0
    tilt_az_deg: float = 0.0
    centre_z_cm: float = 0.0    # hemisphere centre (world z)

    def area(self) -> float:
        if self.mode == 3:
            return 2.0 * PI * self.radius_cm ** 2
        if self.mode == 2:
            return 4.0 * self.half_lx_cm * self.half_ly_cm
        return PI * self.radius_cm ** 2

    def flat_normal(self) -> np.ndarray:
        """World normal of a disk/rectangle (source_normal_world in the Fortran)."""
        a, f = math.radians(self.tilt_deg), math.radians(self.tilt_az_deg)
        if a > 1e-9:
            n = [math.sin(a) * math.cos(f), math.sin(a) * math.sin(f), math.cos(a)]
        else:
            n = [0.0, 0.0, 1.0]
        if self.plane == 2:
            n = [n[0], n[2], n[1]]
        elif self.plane == 3:
            n = [n[2], n[0], n[1]]
        return np.array(n)

    def normals(self, x, y, z) -> np.ndarray:
        """Outward normals at world positions, shape (N, 3)."""
        if self.mode == 3:
            v = np.stack([x, y, z - self.centre_z_cm], axis=1)
            return v / np.linalg.norm(v, axis=1, keepdims=True)
        return np.broadcast_to(self.flat_normal(), (len(np.atleast_1d(x)), 3))


def intensity_p(spectrum, p, cos_theta: float) -> np.ndarray:
    """J(p, θ) = dN/dp [cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹]; zero for θ ≥ 90°.

    ``spectrum`` is a generator spectrum number, or a fast_flux_estimator
    model key (e.g. "reyna_bugaev") for the GUI's planning panels."""
    p = np.atleast_1d(np.asarray(p, float))
    if cos_theta <= 0.0:
        return np.zeros_like(p)
    E = np.sqrt(p * p + M_MU * M_MU)
    theta = math.degrees(math.acos(min(1.0, cos_theta)))
    fn = _ffe._MODELS[spectrum] if isinstance(spectrum, str) else _DNDT[spectrum]
    return fn(E - M_MU, theta) * p / E


def geometry_h(src: Source, c: float) -> float:
    """∫dφ ∫_S max(0, -n·d) dA for zenith cosine c (geometry_h in the Fortran)."""
    if src.mode == 3:
        return 2.0 * PI * PI * src.radius_cm ** 2 * (1.0 + c) / 2.0
    n = src.flat_normal()
    a = n[2] * c
    b = math.hypot(n[0], n[1]) * math.sqrt(max(0.0, 1.0 - c * c))
    if b <= abs(a):
        return src.area() * 2.0 * PI * max(a, 0.0)
    psi0 = math.acos(a / b)
    return src.area() * (2.0 * a * (PI - psi0) + 2.0 * b * math.sin(psi0))


def surface_rate(spectrum, e_min: float, e_max: float,
                 theta_max_deg: float, src: Source) -> float | None:
    """R [s⁻¹], by the Fortran surface_rate() quadrature.

    ``spectrum``: generator spectrum number, or a fast_flux_estimator model
    key. None for spectra this module cannot evaluate: 2 and 8 (no absolute
    normalisation) and 3 (PARMA, which only the Fortran and UCMuGen carry;
    use the rate the generator prints).
    """
    known = spectrum in _ffe._MODELS if isinstance(spectrum, str) else spectrum in _DNDT
    if not known or e_max <= e_min or theta_max_deg <= 0.0:
        return None
    p_min = math.sqrt(max(e_min ** 2 - M_MU ** 2, 0.0))
    p_max = math.sqrt(max(e_max ** 2 - M_MU ** 2, 0.0))
    npts, nc = 400, 400
    p = p_min * np.exp(math.log(p_max / p_min) * np.arange(npts) / (npts - 1))
    c_lo = math.cos(math.radians(theta_max_deg))
    total = 0.0
    for k in range(nc):
        c = c_lo + (k + 0.5) * (1.0 - c_lo) / nc
        hc = geometry_h(src, c)
        if hc <= 0.0:
            continue
        f = intensity_p(spectrum, p, c)
        total += float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(p))) * hc * (1.0 - c_lo) / nc
    return total


def live_time(n_tried: int, rate: float) -> float | None:
    """T = N_tried / R [s]."""
    if not rate or rate <= 0.0 or not n_tried:
        return None
    return float(n_tried) / rate


# ---------------------------------------------------------------------------
#  The densities the legacy Fortran samplers draw from
# ---------------------------------------------------------------------------
def _vertical_momentum_density(spectrum: int, p_min: float, p_max: float, p):
    """f_v(p) of sample_momentum: analytic power law, or the 300-point
    tabulated CDF with linear interpolation (piecewise-constant density)."""
    p = np.asarray(p, float)
    alpha = None
    if spectrum == 1 and p_min < 10.0:
        alpha = -3.1952 + 1.0
    elif spectrum == 2:
        alpha = -2.7
    elif spectrum == 8:
        alpha = -2.0
    if alpha is not None:
        norm = (p_max ** alpha - p_min ** alpha) / alpha
        return p ** (alpha - 1.0) / norm
    grid = p_min * np.exp(np.arange(300) / 299.0 * math.log(p_max / p_min))
    f = intensity_p(spectrum, grid, 1.0)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(grid))])
    dens = np.diff(cdf) / np.diff(grid) / cdf[-1]
    j = np.clip(np.searchsorted(grid, p, side="right") - 1, 0, 298)
    return dens[j]


def _cos_star_vec(c):
    """Guan 2015 Eq. 2, vectorised (fast_flux_estimator._guan_cos_star)."""
    c = np.clip(np.asarray(c, float), 0.0, 1.0)
    num = (c ** 2 + _ffe._GUAN_P1 ** 2 + _ffe._GUAN_P2 * c ** _ffe._GUAN_P3
           + _ffe._GUAN_P4 * c ** _ffe._GUAN_P5)
    return np.sqrt(np.maximum(num, 0.0)) / _ffe._GUAN_DENOM


def _guan_vec(E, c, a, b):
    """Guan 2015 Eq. 3 dN/dE, vectorised over E and c (a=0, b=1: Gaisser)."""
    cs = _cos_star_vec(c)
    with np.errstate(divide="ignore", invalid="ignore"):
        e_eff = np.where(cs > 0.0, E * (1.0 + a / (E * cs ** b)), E)
        phi = (_ffe._GUAN_PRE * e_eff ** _ffe._GUAN_IDX
               * (1.0 / (1.0 + 1.1 * E * cs / _ffe._GUAN_EPI)
                  + _ffe._GUAN_KF / (1.0 + 1.1 * E * cs / _ffe._GUAN_EK)))
    return np.where(np.isfinite(phi) & (phi > 0.0), phi, 0.0)


_GUAN_AB = {4: (3.64, 1.29), 5: (3.512, 1.388), 6: (0.0, 1.0)}


def intensity_vec(spectrum: int, p, c) -> np.ndarray:
    """J(p, θ) for arrays of p and cos θ; same function as intensity_p."""
    p = np.asarray(p, float)
    c = np.asarray(c, float)
    E = np.sqrt(p * p + M_MU * M_MU)
    if spectrum in _GUAN_AB:
        j = _guan_vec(E, c, *_GUAN_AB[spectrum]) * p / E
    elif spectrum == 7:
        zeta = np.maximum(p * np.clip(c, 0.0, None), 1e-300)
        y = np.log10(zeta)
        c1, c2, c3, c4, c5 = _ffe._REYNA_C
        j = np.clip(c, 0.0, None) ** 3 * c1 * zeta ** -(c2 + c3 * y + c4 * y ** 2 + c5 * y ** 3)
    elif spectrum == 1:
        j = 10.0 ** 3.8467 * 1.0e-4 * p ** -3.1952 * np.ones_like(c)
    else:
        raise ValueError(f"no intensity for spectrum {spectrum}")
    return np.where(c > 0.0, j, 0.0)


def _guan_angular_density(spectrum: int, E, c, theta_max_deg: float):
    """Per-steradian density of sample_guan_angle: the pdf guan(E, c)·c on a
    100-point grid, trapezoid CDF, linear inversion (so piecewise-constant
    density in cos θ). Guan's (a, b) unless the spectrum is Frosin, as in the
    Fortran."""
    a, b = (3.512, 1.388) if spectrum == 5 else (3.64, 1.29)
    cmin = math.cos(math.radians(theta_max_deg))
    grid = cmin + np.arange(100) / 99.0 * (1.0 - cmin)
    E = np.asarray(E, float)[:, None]
    pdf = _guan_vec(E, grid[None, :], a, b) * grid[None, :]
    cdf = np.concatenate([np.zeros((len(E), 1)),
                          np.cumsum(0.5 * (pdf[:, 1:] + pdf[:, :-1]) * np.diff(grid), axis=1)],
                         axis=1)
    k = np.clip(np.searchsorted(grid, np.asarray(c, float), side="right") - 1, 0, 98)
    rows = np.arange(len(E))
    return ((cdf[rows, k + 1] - cdf[rows, k]) / (grid[k + 1] - grid[k])
            / cdf[:, -1] / (2.0 * PI))


def _legacy_direction_density(angular_mode, spectrum, E, c_canon, theta_max_deg):
    cm = math.cos(math.radians(theta_max_deg))
    c = np.asarray(c_canon, float)
    inside = (c >= cm - 1e-12) & (c <= 1.0 + 1e-12)
    if angular_mode == 2:
        g = c ** 2 / (2.0 * PI * (1.0 - cm ** 3) / 3.0)
    elif angular_mode == 3:
        g = np.full_like(c, 1.0 / (2.0 * PI * (1.0 - cm)))
    elif angular_mode == 5:
        g = c ** 3 / (2.0 * PI * (1.0 - cm ** 4) / 4.0)
    elif angular_mode == 4:
        g = _guan_angular_density(spectrum, E, c, theta_max_deg)
    else:
        raise ValueError("vertical beam (mode 1) has no density to weight against")
    return np.where(inside, g, 0.0)


def event_weights(spectrum: int, angular_mode: int, e_min: float, e_max: float,
                  theta_max_deg: float, src: Source,
                  p, direction, position, rate: float) -> np.ndarray | None:
    """Importance weights of generated events, averaging to 1 over N_tried.

    ``direction`` and ``position`` are world-frame (N, 3) arrays as the
    generator writes them. Mode 6 returns ones. None when a weight cannot be
    formed: vertical beam, PARMA with a legacy mode, or a spectrum without an
    absolute normalisation.
    """
    p = np.atleast_1d(np.asarray(p, float))
    if angular_mode == 6:
        return np.ones_like(p)
    if angular_mode == 1 or spectrum not in _DNDT or not rate:
        return None
    d = np.asarray(direction, float).reshape(-1, 3)
    x = np.asarray(position, float).reshape(-1, 3)
    E = np.sqrt(p * p + M_MU * M_MU)
    c_world = -d[:, 2]
    # Legacy modes rotate the sky with a vertical source plane; undo that to
    # find the zenith cosine the sampler actually drew.
    if src.mode != 3 and src.plane == 2:
        c_canon = -d[:, 1]
    elif src.mode != 3 and src.plane == 3:
        c_canon = -d[:, 0]
    else:
        c_canon = c_world
    n = src.normals(x[:, 0], x[:, 1], x[:, 2])
    proj = np.clip(-(n * d).sum(axis=1), 0.0, None)
    p_min = math.sqrt(max(e_min ** 2 - M_MU ** 2, 0.0))
    p_max = math.sqrt(max(e_max ** 2 - M_MU ** 2, 0.0))
    J = intensity_vec(spectrum, p, c_world)
    q = (_vertical_momentum_density(spectrum, p_min, p_max, p)
         * _legacy_direction_density(angular_mode, spectrum, E, c_canon, theta_max_deg)
         / src.area())
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(q > 0.0, J * proj / rate / q, 0.0)
    return w


# ---------------------------------------------------------------------------
#  Reading what the generator printed
# ---------------------------------------------------------------------------
_RE_TRIED = re.compile(r"^\s*(?:Total\s+tried\s*\(all ranks\)|Tried)\s*:\s*(\d+)", re.I)
_RE_RATE = re.compile(r"Surface rate R\s*:\s*([0-9.Ee+\-]+)")


def parse_generator_output(lines) -> dict:
    """Tried count and surface rate from ucmuon_gen / ucmuon_gen_omp stdout."""
    out = {"tried": None, "surface_rate": None}
    for line in lines:
        m = _RE_TRIED.search(line)
        if m:
            out["tried"] = int(m.group(1))
        m = _RE_RATE.search(line)
        if m:
            try:
                out["surface_rate"] = float(m.group(1))
            except ValueError:
                pass
    return out
