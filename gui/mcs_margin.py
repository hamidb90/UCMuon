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

Taking the log term over the whole path makes σ_r an upper envelope of the
engines' scattering, on purpose: measured with vertical beams (50 GeV through
40 m, 86 and 171 GeV through 90.5 m of Standard Rock), σ_r is 8-10 % above
MUSIC's RMS radial displacement and about 22 % above UCMuon-MC's (whose
angles are within 4 % of Geant4's; MUSIC's are about 9 % above). So 2σ_r
keeps at least 98 % of a point-like detector's hits whichever engine
transports the muons (UCMuon-MC 99.7 %, MUSIC 99 %). A σ_r matched to
UCMuon-MC would keep only about 95 % with MUSIC.
"""
from __future__ import annotations

import math

import numpy as np

import fast_flux_estimator as _ffe

M_MU = 0.10566
X0_STANDARD_ROCK = 26.54      # g/cm², PDG Standard Rock

# UCMuon-MC material numbers of the materials with their own PDG 2024 table
_MAT_KEYS = {"water": 2, "seawater": 3, "iron": 4, "ice": 6}
_TABLES: dict = {}


def _range_fns(material=None):
    """(R_of_T, T_of_R, R_min) for the CSDA range [g/cm²] against kinetic
    energy [GeV]. Standard Rock (None or "rock") uses the Groom 2001 table of
    fast_flux_estimator; water, seawater, ice and iron use the PDG 2024 tables
    UCMuon-MC transports with. Up to 1.3.2 every material used rock's table:
    water's range per g/cm² is 13 % shorter (its margins came out 12 %
    larger than needed), iron's 12 % longer (11 % smaller)."""
    if material not in _MAT_KEYS:
        return _ffe._R_of_T, _ffe._T_of_R, float(_ffe._GROOM_R_GCM2[0])
    if material not in _TABLES:
        import ucmuon_stochastic_driver as _u
        mat = _u._MAT_DB[_MAT_KEYS[material]]
        T = _u._PROC_T                                   # MeV
        a, Lb, Lp, Ln = _u._loss_components(T, mat)
        tab, s_a, _ = _u._select_proc(mat)
        inv = 1.0 / (a + Lb + Lp + Ln)
        R = tab["R10"] / s_a + np.concatenate(
            ([0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(T))))
        lT, lR = np.log(T / 1000.0), np.log(R)
        _TABLES[material] = (
            lambda t: np.exp(np.interp(np.log(np.maximum(np.asarray(t, float), 1e-12)), lT, lR)),
            lambda r: np.exp(np.interp(np.log(np.maximum(np.asarray(r, float), 1e-12)), lR, lT)),
            float(R[0]))
    return _TABLES[material]


def threshold_kinetic_gev(opacity_gcm2: float, material=None) -> float:
    """Kinetic energy [GeV] whose CSDA range is opacity_gcm2."""
    return float(_range_fns(material)[1](opacity_gcm2))


def sigma_r_cm(e_surface_total_gev: float, slant_path_cm: float,
               rho: float = 2.65, x0_gcm2: float = X0_STANDARD_ROCK,
               n: int = 2000, material=None) -> float | None:
    """RMS radial displacement [cm] after slant_path_cm of rock (or of the
    material, see _range_fns), with energy loss. None when the muon stops
    before the end of the path."""
    R_of_T, T_of_R, R_min = _range_fns(material)
    t0 = e_surface_total_gev - M_MU
    if t0 <= 0.0 or slant_path_cm <= 0.0:
        return None
    r0 = float(R_of_T(t0))
    s = np.linspace(0.0, slant_path_cm, n)
    rem = r0 - rho * s
    if rem[-1] <= R_min:
        return None
    t = np.asarray(T_of_R(rem), float)
    e = t + M_MU
    p = np.sqrt(e * e - M_MU * M_MU)
    beta = p / e
    corr = 1.0 + 0.038 * math.log(max(rho * slant_path_cm / x0_gcm2, 1e-12))
    theta2_per_cm = (0.0136 / (beta * p)) ** 2 * corr ** 2 / (x0_gcm2 / rho)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    sy2 = float(trapz((slant_path_cm - s) ** 2 * theta2_per_cm, s))
    return math.sqrt(2.0 * sy2)


def sigma_r_cm_many(e_surface_total_gev, slant_path_cm, rho: float = 2.65,
                    x0_gcm2: float = X0_STANDARD_ROCK, n: int = 400,
                    material=None) -> np.ndarray:
    """sigma_r_cm for arrays of muons (the same integral, vectorised); NaN
    where the muon stops before the end of its path."""
    R_of_T, T_of_R, R_min = _range_fns(material)
    e0 = np.atleast_1d(np.asarray(e_surface_total_gev, float))
    L = np.broadcast_to(np.asarray(slant_path_cm, float), e0.shape)
    out = np.full(e0.shape, np.nan)
    ok = (e0 > M_MU) & (L > 0.0)
    if not ok.any():
        return out
    r0 = np.asarray(R_of_T(e0[ok] - M_MU), float)
    u = np.linspace(0.0, 1.0, n)
    s = L[ok][:, None] * u[None, :]
    rem = r0[:, None] - rho * s
    live = rem[:, -1] > R_min
    t = np.asarray(T_of_R(np.maximum(rem, R_min)), float)
    e = t + M_MU
    p = np.sqrt(e * e - M_MU * M_MU)
    beta = p / e
    corr = 1.0 + 0.038 * np.log(np.maximum(rho * L[ok] / x0_gcm2, 1e-12))
    theta2_per_cm = (0.0136 / (beta * p)) ** 2 * (corr ** 2)[:, None] / (x0_gcm2 / rho)
    f = (L[ok][:, None] - s) ** 2 * theta2_per_cm
    sy2 = np.sum(0.5 * (f[:, 1:] + f[:, :-1]) * np.diff(s, axis=1), axis=1)
    res = np.sqrt(2.0 * sy2)
    res[~live] = np.nan
    out[ok] = res
    return out


def arriving_energy_floor(slant_path_cm: float, rho: float = 2.65,
                          t_arrive_gev: float = 1.0, material=None) -> float:
    """Total surface energy [GeV] of a muon reaching the end of the path with
    kinetic energy t_arrive_gev."""
    R_of_T, T_of_R, _ = _range_fns(material)
    r = float(R_of_T(t_arrive_gev)) + rho * slant_path_cm
    return float(T_of_R(r)) + M_MU


def suggested_margin(depth_m: float, rho: float, e_min_total_gev: float,
                     theta_deg: float = 0.0, x0_gcm2: float = X0_STANDARD_ROCK,
                     k: float = 2.0, material=None) -> dict:
    """Margin recommendation for a detector at vertical depth depth_m.

    The displacement is evaluated for the slowest muon that matters: the
    generated E_min, or, if that cannot cross the slant path, the surface
    energy that arrives with 1 GeV (muons arriving with less scatter even
    more). k = 2 keeps ~98% of the hits of a point-like detector.
    """
    c = max(math.cos(math.radians(theta_deg)), 0.02)
    path_cm = depth_m / c * 100.0
    e_floor = arriving_energy_floor(path_cm, rho, material=material)
    e_eff = max(e_min_total_gev, e_floor)
    sig = sigma_r_cm(e_eff, path_cm, rho, x0_gcm2, material=material)
    return {"sigma_r_cm": sig, "margin_cm": None if sig is None else k * sig,
            "e_eval_gev": e_eff, "e_floor_gev": e_floor, "slant_m": path_cm / 100.0,
            "k": k}


def retained_fraction(margin_cm: float, sigma_cm: float | None) -> float | None:
    """Fraction of hits a point-like detector keeps with this margin."""
    if not sigma_cm:
        return None
    return 1.0 - math.exp(-(margin_cm / sigma_cm) ** 2)
