"""
fast_flux_estimator.py  —  UCMuon analytical muon flux utilities
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026

Provides fast analytical estimates for muon flux, transmission, exposure time,
and minimum detectable energy through rock overburdens.

Five sea-level flux models are implemented:
  "reyna_bugaev"  — Reyna (2006) / Bugaev et al. (1998) [recommended]
  "bugaev"        — Bugaev et al. (1998) / Gaisser (1990)
  "gaisser_tang"  — Gaisser & Tang (1984) / PDG 2022 §30
  "guan_2015"     — Guan et al. (2015) arXiv:1509.06176 — correct cosθ* (Earth curvature)
  "frosin_2025"   — Frosin et al. (2025) J. Phys. G 52, 035002 — re-fitted on 304 datasets

All five models are azimuth-symmetric at sea level.  Azimuth dependence
requires the PARMA interface (spectrum mode ③ in the generator tab).

CSDA range table:
  Groom, Mokhov & Striganov (2001) ADNDT 78, 183  (Standard Rock)

Public API (matches cosmoaleph_gui.py imports):
  integrated_flux(opacity_gcm2, theta_deg, model, altitude_m)
      → (I_flux [cm⁻²sr⁻¹s⁻¹], E_min [GeV] | None)

  flux_vs_depth(depths_m, rho, theta_deg, model, altitude_m)
      → (I_arr, T_arr, Emin_arr)  all numpy arrays

  differential_flux(T_GeV, theta_deg, model, altitude_m)
      → ndarray  dΦ/dT  [cm⁻²s⁻¹sr⁻¹GeV⁻¹]

  angular_profile(theta_arr_deg, E_min_GeV, model, altitude_m)
      → (I_arr [cm⁻²sr⁻¹s⁻¹], T_arr normalised)

  exposure_time(n_muons, flux_cm2_sr_s, acceptance_cm2_sr)
      → t [s]

  emin_from_opacity(opacity_gcm2)
      → E_min [GeV]

  MODEL_LABELS : dict[str, str]
      Human-readable model names for use in GUI dropdowns.

  RHO_STANDARD_ROCK = 2.65  (g/cm³)
"""
from __future__ import annotations
import math
import numpy as np

# ---------------------------------------------------------------------------
RHO_STANDARD_ROCK: float = 2.65   # g/cm³  (PDG Standard Rock)
M_MU_GEV: float = 0.10566         # muon rest mass [GeV]

# ---------------------------------------------------------------------------
#  Groom (2001) CSDA range table for Standard Rock
#  T [GeV] → R [g/cm²]
# ---------------------------------------------------------------------------
_GROOM_T_GEV = np.array([
    0.01, 0.014, 0.02, 0.03, 0.04, 0.08, 0.10, 0.14, 0.20, 0.30,
    0.40, 0.80, 1.00, 1.40, 2.00, 3.00, 4.00, 8.00,
    10.0, 14.0, 20.0, 30.0, 40.0, 80.0, 100.0,
    140.0, 200.0, 300.0, 400.0, 800.0, 1000.0,
])
# Values from Groom (2001) Table 26.1 Standard Rock, divided by 100 for correct g/cm² units.
_GROOM_R_GCM2 = np.array([
    0.8516, 1.542, 2.866, 5.698, 9.145, 26.76, 36.96, 58.79, 93.32, 152.4,
    211.5, 441.8, 553.4, 771.2, 1088., 1599., 2095., 3998.,
    4920., 6724., 9360., 13620., 17760., 33430., 40840.,
    55460., 76650., 107900., 136100., 225300., 272200.,
])  # g/cm²  — Groom (2001) Table 26.1 Standard Rock (already in g/cm²)

# Log-log interpolator: T → R
_log_T_tab = np.log(_GROOM_T_GEV)
_log_R_tab = np.log(_GROOM_R_GCM2)


def _R_of_T(T_GeV: float | np.ndarray) -> float | np.ndarray:
    """CSDA range [g/cm²] for kinetic energy T [GeV] — log-log interpolation."""
    scalar = np.ndim(T_GeV) == 0
    T = np.atleast_1d(np.asarray(T_GeV, dtype=float))
    logT = np.log(np.clip(T, _GROOM_T_GEV[0], _GROOM_T_GEV[-1]))
    logR = np.interp(logT, _log_T_tab, _log_R_tab)
    R = np.exp(logR)
    return float(R[0]) if scalar else R


def _T_of_R(R_gcm2: float | np.ndarray) -> float | np.ndarray:
    """Inverse CSDA: opacity [g/cm²] → minimum kinetic energy [GeV]."""
    scalar = np.ndim(R_gcm2) == 0
    R = np.atleast_1d(np.asarray(R_gcm2, dtype=float))
    logR = np.log(np.clip(R, _GROOM_R_GCM2[0], _GROOM_R_GCM2[-1]))
    logT = np.interp(logR, _log_R_tab, _log_T_tab)
    T = np.exp(logT)
    return float(T[0]) if scalar else T


# ---------------------------------------------------------------------------
#  Sea-level flux models
#  All return differential flux dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹] at (T [GeV], θ [deg])
# ---------------------------------------------------------------------------

def _reyna_bugaev(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Reyna (2006) / Bugaev (1998) vertical muon spectrum.
    Reyna Eq. 3: log₁₀(p³ × I_vert) = c₀x³+c₁x²+c₂x+c₃  (x = log₁₀(p))
    I_vert = 10^polynomial / p³ × 100   [cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹]

    Derivation of the ×100 factor:
      Reyna fitted coefficients with I_vert in 10⁻² m⁻²sr⁻¹s⁻¹(GeV/c)⁻¹.
      ×10⁻² (implicit) × ×10⁴ (m²→cm²) = ×100.
    Validation: integrated above 1 GeV → ~6.2×10⁻³ cm⁻²sr⁻¹s⁻¹
                PDG reference: ~7×10⁻³ cm⁻²sr⁻¹s⁻¹ (12% agreement) ✓
    """
    cos_th = math.cos(math.radians(theta_deg))
    p = np.maximum(np.sqrt((T_GeV + M_MU_GEV)**2 - M_MU_GEV**2), 0.5)  # [GeV/c]
    log10_p = np.log10(p)
    c0, c1, c2, c3 = 0.00253, -0.2455, 1.288, -4.25
    log10_p3I = c0*log10_p**3 + c1*log10_p**2 + c2*log10_p + c3
    # I_vert [cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹]
    phi_vert = 10.0**log10_p3I / (p**3) * 100.0
    # Angular dependence: cos²θ* approximation (Reyna Eq. 4)
    cos_th_star = np.sqrt((cos_th**2 + 0.102**2) / (1.0 + 0.102**2))
    phi_p = phi_vert * cos_th_star**1.85
    # Convert dΦ/dp → dΦ/dT via Jacobian dp/dT = E/p (T=E−m → dT=dE, dp/dE=E/p)
    E = T_GeV + M_MU_GEV
    return phi_p * (E / p)   # [cm⁻²s⁻¹sr⁻¹GeV⁻¹]


def _bugaev(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Bugaev et al. (1998) / Gaisser parametrisation.
    dΦ/dE [cm⁻²s⁻¹sr⁻¹GeV⁻¹]
    """
    cos_th = math.cos(math.radians(theta_deg))
    E = T_GeV + M_MU_GEV
    # Gaisser approximation with pion and kaon terms
    # Gaisser (1990) Eq. 15.3; A in cm⁻²s⁻¹sr⁻¹GeV⁻¹
    A = 1.4e-2
    phi = (A * E**(-2.7) *
           (1.0 / (1.0 + 1.1*E*cos_th/115.0) +
            0.054 / (1.0 + 1.1*E*cos_th/850.0)))
    return phi


def _gaisser_tang(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Gaisser & Tang (1984) / PDG 2022 §30.3 Eq. 30.6
    dΦ/dE [cm⁻²s⁻¹sr⁻¹GeV⁻¹]
    """
    cos_th = math.cos(math.radians(theta_deg))
    E = T_GeV + M_MU_GEV
    A = 1.4e-2
    phi = (A * E**(-2.7) *
           (1.0 / (1.0 + 1.1*E*cos_th/115.0) +
            0.054 / (1.0 + 1.1*E*cos_th/850.0)) *
           (1.0 + 0.054*E/800.0))
    return phi


# ---------------------------------------------------------------------------
#  Guan / Frosin models — exact Python transcription of cosmoaleph_module_omp.f90
# ---------------------------------------------------------------------------
# cosθ* parameters (Guan 2015, Table 1 / arXiv:1509.06176)
_GUAN_P1    =  0.102573
_GUAN_P2    = -0.068287
_GUAN_P3    =  0.958633
_GUAN_P4    =  0.0407253
_GUAN_P5    =  0.817285
_GUAN_DENOM =  0.99144315   # = sqrt(1 + P1² + P2 + P4) normalisation factor
_GUAN_EPI   =  115.0        # [GeV]  effective pion critical energy
_GUAN_EK    =  850.0        # [GeV]  effective kaon critical energy
_GUAN_KF    =    0.054      # kaon / (pion + kaon) fraction
_GUAN_PRE   =    0.14       # normalisation [cm⁻²s⁻¹sr⁻¹GeV^{1.7}] — matches Fortran GUAN_PRE
_GUAN_IDX   =   -2.7        # Gaisser spectral index


def _guan_cos_star(cos_th: float) -> float:
    """
    Effective cosθ* that accounts for Earth's curvature at large zenith angles.
    Exact match to the Fortran `guan_cos_star` in cosmoaleph_module_omp.f90.

    Parameters
    ----------
    cos_th : float
        cos(zenith angle), clamped to [0, 1].

    Returns
    -------
    float
        cosθ* ∈ [0, 1].
    """
    cos_th = max(0.0, min(1.0, cos_th))
    numer = (cos_th**2 + _GUAN_P1**2
             + _GUAN_P2 * cos_th**_GUAN_P3
             + _GUAN_P4 * cos_th**_GUAN_P5)
    return math.sqrt(max(0.0, numer)) / _GUAN_DENOM


def _guan_frosin(T_GeV: np.ndarray, theta_deg: float, a: float, b: float) -> np.ndarray:
    """
    General Guan / Frosin formula.

    dΦ/dT = PRE × E_eff^{-2.7} × (pion_term + kaon_term)

    where
        E_eff = E × (1 + a / (E × cosθ*^b))   [energy-shift from low-E correction]
        pion_term = 1 / (1 + 1.1 E cosθ* / 115)
        kaon_term = 0.054 / (1 + 1.1 E cosθ* / 850)

    Parameters a, b:
        Guan 2015 (arXiv:1509.06176) : a = 3.64,  b = 1.29
        Frosin 2025 (JPG 52, 035002) : a = 3.512, b = 1.388

    Returns
    -------
    ndarray  dΦ/dT  [cm⁻²s⁻¹sr⁻¹GeV⁻¹]
    """
    cos_th = math.cos(math.radians(theta_deg))
    cs = _guan_cos_star(cos_th)                 # scalar effective cosθ*
    E  = T_GeV + M_MU_GEV                       # total energy [GeV]
    # Low-energy correction term from Guan (2015): shifts E slightly upward at low E
    if cs > 0.0:
        E_eff = E * (1.0 + a / (E * cs**b))
    else:
        E_eff = E
    pion_t = 1.0        / (1.0 + 1.1 * E * cs / _GUAN_EPI)
    kaon_t = _GUAN_KF   / (1.0 + 1.1 * E * cs / _GUAN_EK)
    phi = _GUAN_PRE * np.power(E_eff, _GUAN_IDX) * (pion_t + kaon_t)
    return np.where(np.isfinite(phi) & (phi >= 0), phi, 0.0)


def _guan_2015(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """Guan et al. (2015) arXiv:1509.06176.  a=3.64, b=1.29."""
    return _guan_frosin(T_GeV, theta_deg, a=3.64, b=1.29)


def _frosin_2025(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """Frosin et al. (2025) J. Phys. G 52, 035002.  a=3.512, b=1.388."""
    return _guan_frosin(T_GeV, theta_deg, a=3.512, b=1.388)


_MODELS = {
    "reyna_bugaev": _reyna_bugaev,
    "bugaev":       _bugaev,
    "gaisser_tang": _gaisser_tang,
    "guan_2015":    _guan_2015,
    "frosin_2025":  _frosin_2025,
}

# Human-readable labels (used in the GUI dropdown)
MODEL_LABELS: dict[str, str] = {
    "reyna_bugaev": "Reyna–Bugaev (2006) ← recommended",
    "bugaev":       "Bugaev (1998) / Gaisser",
    "gaisser_tang": "Gaisser–Tang (1984)",
    "guan_2015":    "Guan et al. (2015)  [a=3.64, b=1.29]",
    "frosin_2025":  "Frosin et al. (2025) [a=3.512, b=1.388]",
}

_ALTITUDE_FACTOR = {
    # Approximate flux scaling vs altitude [m a.s.l.] — from PDG §30
    # sea level = 1.0; multiplicative correction
    # Uses exp(h/h_scale) with h_scale ≈ 8500 m for muons
}

def _altitude_correction(altitude_m: float) -> float:
    """
    Approximate muon flux correction for altitude above sea level.
    Based on exponential atmosphere model: φ(h) ≈ φ₀ · exp(h / 8500).
    """
    return math.exp(altitude_m / 8500.0)


# ---------------------------------------------------------------------------
#  Integration grid
# ---------------------------------------------------------------------------
_T_GRID = np.logspace(np.log10(0.5), np.log10(1.5e4), 600)  # 0.5 GeV → 15 TeV


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------

def emin_from_opacity(opacity_gcm2: float) -> float | None:
    """
    Minimum muon kinetic energy [GeV] to traverse opacity_gcm2 [g/cm²].
    Returns None if opacity exceeds the maximum tabulated CSDA range.
    """
    if opacity_gcm2 <= 0.0:
        return 0.0
    if opacity_gcm2 >= _GROOM_R_GCM2[-1]:
        return None   # above table maximum (~27,000 km.w.e.)
    return float(_T_of_R(opacity_gcm2))


def integrated_flux(
    opacity_gcm2: float,
    theta_deg: float = 0.0,
    model: str = "reyna_bugaev",
    altitude_m: float = 0.0,
) -> tuple[float, float | None]:
    """
    Integrated muon flux after traversing opacity_gcm2 [g/cm²] of rock.

    Parameters
    ----------
    opacity_gcm2 : float
        Rock opacity X = ρ·L [g/cm²]. Pass 0 for open-sky (no rock).
    theta_deg : float
        Zenith angle [°].
    model : str
        Flux model: "reyna_bugaev" (recommended), "bugaev", "gaisser_tang".
    altitude_m : float
        Altitude above sea level [m] for flux correction.

    Returns
    -------
    I_flux : float
        Integrated flux [cm⁻²sr⁻¹s⁻¹] above E_min.
    E_min : float | None
        Minimum kinetic energy [GeV] to traverse opacity. None if too deep.
    """
    E_min = emin_from_opacity(opacity_gcm2)
    if E_min is None:
        return 0.0, None

    flux_fn = _MODELS.get(model, _reyna_bugaev)
    alt_corr = _altitude_correction(altitude_m)

    # Integrate dΦ/dT from E_min to T_max
    mask = _T_GRID >= E_min
    T = _T_GRID[mask]
    if len(T) < 2:
        return 0.0, E_min

    phi = flux_fn(T, theta_deg) * alt_corr
    phi = np.where(np.isfinite(phi) & (phi > 0), phi, 0.0)
    _integrate = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
    I = float(_integrate(phi, T))
    return I, float(E_min)


def flux_vs_depth(
    depths_m: np.ndarray,
    rho: float = RHO_STANDARD_ROCK,
    theta_deg: float = 0.0,
    model: str = "reyna_bugaev",
    altitude_m: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute muon flux, transmission and E_min vs depth.

    Parameters
    ----------
    depths_m : array-like [m]
        Rock thickness values.
    rho : float [g/cm³]
        Rock density.
    theta_deg : float [°]
        Zenith angle.
    model : str
        Flux model.
    altitude_m : float [m]
        Altitude above sea level.

    Returns
    -------
    I_arr   : ndarray  integrated flux [cm⁻²sr⁻¹s⁻¹]
    T_arr   : ndarray  transmission = I(L) / I(0)
    Emin_arr: ndarray  E_min [GeV]  (0.0 where below table)
    """
    depths_m = np.asarray(depths_m, dtype=float)
    opacities = rho * depths_m * 100.0   # [g/cm²]

    I_open, _ = integrated_flux(0.0, theta_deg, model, altitude_m)

    I_arr    = np.zeros(len(opacities))
    T_arr    = np.zeros(len(opacities))
    Emin_arr = np.zeros(len(opacities))

    for i, X in enumerate(opacities):
        I, Emin = integrated_flux(float(X), theta_deg, model, altitude_m)
        I_arr[i]    = I
        T_arr[i]    = (I / I_open) if I_open > 0 else 0.0
        Emin_arr[i] = Emin if Emin is not None else 0.0

    return I_arr, T_arr, Emin_arr


def differential_flux(
    T_GeV: np.ndarray,
    theta_deg: float = 0.0,
    model: str = "reyna_bugaev",
    altitude_m: float = 0.0,
) -> np.ndarray:
    """
    Differential muon flux dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹] at sea level (or altitude).

    Parameters
    ----------
    T_GeV : array-like
        Muon kinetic energy [GeV].
    theta_deg : float
        Zenith angle [°]. All five models are azimuth-symmetric at sea level;
        azimuth dependence requires the PARMA interface (spectrum mode ③).
    model : str
        One of: "reyna_bugaev", "bugaev", "gaisser_tang", "guan_2015", "frosin_2025".
    altitude_m : float
        Altitude above sea level [m].

    Returns
    -------
    ndarray  dΦ/dT  [cm⁻²s⁻¹sr⁻¹GeV⁻¹]
    """
    T = np.atleast_1d(np.asarray(T_GeV, dtype=float))
    flux_fn  = _MODELS.get(model, _reyna_bugaev)
    alt_corr = _altitude_correction(altitude_m)
    phi = flux_fn(T, theta_deg) * alt_corr
    return np.where(np.isfinite(phi) & (phi >= 0), phi, 0.0)


def angular_profile(
    theta_arr_deg: np.ndarray,
    E_min_GeV: float = 1.0,
    model: str = "reyna_bugaev",
    altitude_m: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Integrated muon flux I(θ) [cm⁻²sr⁻¹s⁻¹] and normalised ratio I(θ)/I(0°)
    for an array of zenith angles, integrated above E_min.

    Parameters
    ----------
    theta_arr_deg : array-like
        Zenith angles [°].
    E_min_GeV : float
        Lower energy cut-off [GeV].
    model : str
        Flux model key.
    altitude_m : float
        Altitude [m].

    Returns
    -------
    I_arr : ndarray  integrated flux [cm⁻²sr⁻¹s⁻¹]
    T_arr : ndarray  normalised ratio I(θ)/I(0°)  (transmission analogue)
    """
    theta_arr = np.asarray(theta_arr_deg, dtype=float)
    flux_fn   = _MODELS.get(model, _reyna_bugaev)
    alt_corr  = _altitude_correction(altitude_m)
    _integrate = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)

    # Energy grid clipped to user E_min
    mask = _T_GRID >= max(E_min_GeV, _T_GRID[0])
    T    = _T_GRID[mask]
    if len(T) < 2:
        z = np.zeros(len(theta_arr))
        return z, z

    I_arr = np.zeros(len(theta_arr))
    for i, th in enumerate(theta_arr):
        phi = flux_fn(T, float(th)) * alt_corr
        phi = np.where(np.isfinite(phi) & (phi >= 0), phi, 0.0)
        I_arr[i] = float(_integrate(phi, T))

    I_0 = I_arr[0] if I_arr[0] > 0 else 1.0
    T_arr = I_arr / I_0
    return I_arr, T_arr


def exposure_time(
    n_muons: int,
    flux_cm2_sr_s: float,
    acceptance_cm2_sr: float,
) -> float:
    """
    Equivalent measurement time [s] for n_muons simulated muons.

    Parameters
    ----------
    n_muons : int
        Number of simulated muons above E_min.
    flux_cm2_sr_s : float
        Integrated muon flux [cm⁻²sr⁻¹s⁻¹].
    acceptance_cm2_sr : float
        Detector acceptance [cm²·sr].

    Returns
    -------
    t : float  [s]  (inf if rate = 0)
    """
    rate = flux_cm2_sr_s * acceptance_cm2_sr
    if rate <= 0.0:
        return float('inf')
    return float(n_muons) / rate
