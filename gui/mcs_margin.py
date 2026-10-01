"""
mcs_margin.py: how large the detector-filter safety margin must be
UCLouvain Muography Group | MIT License 2026

The generator's detector filter keeps a muon when its STRAIGHT line from the
source crosses the detector inflated by the safety margin. Multiple Coulomb
scattering in the overburden then moves it sideways by σ_r (RMS, radial). A
muon whose straight line misses the inflated detector is never generated,
even if it would have scattered into the detector, so a margin much smaller
than σ_r silently loses hits and biases every rate and live time low.

σ_r is the Highland–Lynch–Dahl lateral displacement (PDG, Passage of
particles through matter, Sec. 34.3), integrated along the path with the
momentum falling by CSDA energy loss (Groom et al. 2001 range table,
Standard Rock):

    σ_y² = ∫₀ᴸ (L − s)² [13.6 MeV / (β c p(s))]² [1 + 0.038 ln(ρL/X₀)]² ds / (X₀/ρ)
    σ_r  = √2 σ_y

At constant momentum this is L·θ₀/√3 per projection, the formula the GUI's
older helper used; with energy loss it grows sharply for muons that arrive
with little energy left. For a point-like detector, a Gaussian displacement
keeps a fraction 1 − exp(−m²/σ_r²) of the hits with margin m: 1.5σ_r keeps
89%, 2σ_r 98%, 3σ_r 99.99%.
"""
from __future__ import annotations

import math

import numpy as np

import fast_flux_estimator as _ffe

M_MU = 0.10566
X0_STANDARD_ROCK = 26.54      # g/cm², PDG Standard Rock


def sigma_r_cm(e_surface_total_gev: float, slant_path_cm: float,
               rho: float = 2.65, x0_gcm2: float = X0_STANDARD_ROCK,
               n: int = 2000) -> float | None:
    """RMS radial displacement [cm] after slant_path_cm of rock, with energy
    loss. None when the muon stops before the end of the path."""
    t0 = e_surface_total_gev - M_MU
    if t0 <= 0.0 or slant_path_cm <= 0.0:
        return None
    r0 = float(_ffe._R_of_T(t0))
    s = np.linspace(0.0, slant_path_cm, n)
    rem = r0 - rho * s
    if rem[-1] <= _ffe._GROOM_R_GCM2[0]:
        return None
    t = np.asarray(_ffe._T_of_R(rem), float)
    e = t + M_MU
    p = np.sqrt(e * e - M_MU * M_MU)
    beta = p / e
    corr = 1.0 + 0.038 * math.log(max(rho * slant_path_cm / x0_gcm2, 1e-12))
    theta2_per_cm = (0.0136 / (beta * p)) ** 2 * corr ** 2 / (x0_gcm2 / rho)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    sy2 = float(trapz((slant_path_cm - s) ** 2 * theta2_per_cm, s))
    return math.sqrt(2.0 * sy2)


def arriving_energy_floor(slant_path_cm: float, rho: float = 2.65,
                          t_arrive_gev: float = 1.0) -> float:
    """Total surface energy [GeV] of a muon reaching the end of the path with
    kinetic energy t_arrive_gev."""
    r = float(_ffe._R_of_T(t_arrive_gev)) + rho * slant_path_cm
    return float(_ffe._T_of_R(r)) + M_MU


def suggested_margin(depth_m: float, rho: float, e_min_total_gev: float,
                     theta_deg: float = 0.0, x0_gcm2: float = X0_STANDARD_ROCK,
                     k: float = 2.0) -> dict:
    """Margin recommendation for a detector at vertical depth depth_m.

    The displacement is evaluated for the slowest muon that matters: the
    generated E_min, or, if that cannot cross the slant path, the surface
    energy that arrives with 1 GeV (muons arriving with less scatter even
    more). k = 2 keeps ~98% of the hits of a point-like detector.
    """
    c = max(math.cos(math.radians(theta_deg)), 0.02)
    path_cm = depth_m / c * 100.0
    e_floor = arriving_energy_floor(path_cm, rho)
    e_eff = max(e_min_total_gev, e_floor)
    sig = sigma_r_cm(e_eff, path_cm, rho, x0_gcm2)
    return {"sigma_r_cm": sig, "margin_cm": None if sig is None else k * sig,
            "e_eval_gev": e_eff, "e_floor_gev": e_floor, "slant_m": path_cm / 100.0,
            "k": k}


def retained_fraction(margin_cm: float, sigma_cm: float | None) -> float | None:
    """Fraction of hits a point-like detector keeps with this margin."""
    if not sigma_cm:
        return None
    return 1.0 - math.exp(-(margin_cm / sigma_cm) ** 2)
