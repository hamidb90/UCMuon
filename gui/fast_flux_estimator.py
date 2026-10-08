"""
fast_flux_estimator.py  —  UCMuon analytical muon flux utilities
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>
MIT License 2026

Provides fast analytical estimates for muon flux, transmission, exposure time,
and minimum detectable energy through rock overburdens.

Five sea-level flux models are implemented (the keys are historical; see
docs/FLUX_NORMALISATION_AUDIT.md for the audit against the source papers):
  "reyna_bugaev"  Reyna (2006) Eqs. 1-3, = generator spectrum 7 [recommended]
  "bugaev"        Gaisser (1990) formula (PDG 2022 Eq. 30.4) with Guan's cosθ*,
                    = generator spectrum 6.  Not Bugaev (1998).  Valid E > 100/cosθ
                    GeV only
  "gaisser_tang"  Tang et al. (2006) modified Gaisser formula, PRD 74, 053007
  "guan_2015"     Guan et al. (2015) arXiv:1509.06176 Eq. 3, = spectrum 4
  "frosin_2025"   Frosin et al. (2025) J. Phys. G 52, 035002, = spectrum 5

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
    140.0, 200.0, 300.0, 400.0, 800.0, 1000.0, 1400.0, 2000.0,
])
# Groom, Mokhov & Striganov, ADNDT 78 (2001), Table IV-6, Standard Rock, in
# g/cm²; the same values as _GROOM_R_GCM2 in ucmuon_gui.py and the backward
# MC. Up to 1.3.0 the entries from 140 GeV on were wrong (1 % long at 140 GeV
# rising to 11 % at 1 TeV) and the table stopped at 1 TeV.
_GROOM_R_GCM2 = np.array([
    0.8516, 1.542, 2.866, 5.698, 9.145, 26.76, 36.96, 58.79, 93.32, 152.4,
    211.5, 441.8, 553.4, 771.2, 1088., 1599., 2095., 3998.,
    4920., 6724., 9360., 13620., 17760., 33430., 40840.,
    54950., 74590., 104000., 130200., 212900., 245300., 299000., 361600.,
])

# Beyond the table (2 TeV, 3616 m w.e.) the range is continued with
# dE/dX = a + b E, fitted to the slopes of its last three points:
# a = 2.31e-3 GeV cm²/g, b = 4.28e-6 cm²/g (standard-rock values).
_A_LOSS, _B_LOSS = 2.31e-3, 4.28e-6

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
    hi = T > _GROOM_T_GEV[-1]       # continue with dE/dX = a + bE
    if hi.any():
        R[hi] = _GROOM_R_GCM2[-1] + np.log((_A_LOSS + _B_LOSS * T[hi]) /
                                           (_A_LOSS + _B_LOSS * _GROOM_T_GEV[-1])) / _B_LOSS
    return float(R[0]) if scalar else R


def _T_of_R(R_gcm2: float | np.ndarray) -> float | np.ndarray:
    """Inverse CSDA: opacity [g/cm²] → minimum kinetic energy [GeV]."""
    scalar = np.ndim(R_gcm2) == 0
    R = np.atleast_1d(np.asarray(R_gcm2, dtype=float))
    logR = np.log(np.clip(R, _GROOM_R_GCM2[0], _GROOM_R_GCM2[-1]))
    logT = np.interp(logR, _log_R_tab, _log_T_tab)
    T = np.exp(logT)
    hi = R > _GROOM_R_GCM2[-1]      # continue with dE/dX = a + bE
    if hi.any():
        e0 = _GROOM_T_GEV[-1] + _A_LOSS / _B_LOSS
        T[hi] = e0 * np.exp(_B_LOSS * (R[hi] - _GROOM_R_GCM2[-1])) - _A_LOSS / _B_LOSS
    return float(T[0]) if scalar else T


# ---------------------------------------------------------------------------
#  Sea-level flux models
#  All return differential flux dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹] at (T [GeV], θ [deg])
# ---------------------------------------------------------------------------

# Reyna (2006), arXiv:hep-ph/0604145: Eq. 3 with the Sec. 4 "Best Fit"
# coefficients, I_V in cm^-2 s^-1 sr^-1 (GeV/c)^-1.  Same constants as
# REYNA_C0..C4 in src/generator/ucmuon_source_module.f90 and UCMuGen.h.
_REYNA_C = (0.00253, 0.2455, 1.288, -0.2555, 0.0209)


def _reyna_bugaev(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Reyna (2006) differential intensity, arXiv:hep-ph/0604145:

        I(p, θ) = cos³θ · I_V(ζ),   ζ = p cosθ                    (Eqs. 1-2)
        I_V(x)  = c1 · x^-(c2 + c3 y + c4 y² + c5 y³),  y = log10 x   (Eq. 3)

    with (c1..c5) = (0.00253, 0.2455, 1.288, -0.2555, 0.0209) in
    cm⁻²s⁻¹sr⁻¹(GeV/c)⁻¹, the same function as generator spectrum 7.
    Vertical integral above 1 GeV/c: 7.02e-3 cm⁻²s⁻¹sr⁻¹ (PDG: ≈7e-3).
    Validity (Sec. 4): 1 GeV/c < p < 2000 GeV/c / cosθ.

    Returns dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹], i.e. dΦ/dp · dp/dT with dp/dT = E/p.
    """
    cos_th = math.cos(math.radians(theta_deg))
    T = np.asarray(T_GeV, dtype=float)
    if cos_th <= 0.0:
        return np.zeros_like(T)
    E = T + M_MU_GEV
    p = np.sqrt(np.maximum(E**2 - M_MU_GEV**2, 0.0))          # [GeV/c]
    zeta = np.maximum(p * cos_th, 1e-300)
    y = np.log10(zeta)
    c1, c2, c3, c4, c5 = _REYNA_C
    phi_p = cos_th**3 * c1 * zeta**-(c2 + c3*y + c4*y**2 + c5*y**3)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p > 0.0, phi_p * E / p, 0.0)


def _bugaev(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Gaisser (1990) pion+kaon formula (PDG 2022, Cosmic Rays, Eq. 30.4; Guan
    2015 Eq. 1):

        dΦ/dE = 0.14 E^-2.7 [1/(1 + 1.1 E cosθ*/115) + 0.054/(1 + 1.1 E cosθ*/850)]

    in cm⁻²s⁻¹sr⁻¹GeV⁻¹, with Guan's cosθ* (Eq. 2) for θ, i.e. Guan Eq. 3 with
    a = 0, b = 1 and the same function as generator spectrum 6.  The key is
    historical: this is not the Bugaev et al. (1998) fit.
    Valid only for E > 100/cosθ GeV (muon decay and energy loss neglected);
    below that it overestimates the flux (x12 the PDG vertical integral above
    1 GeV).  The prefactor was 1.4e-2 before the 2026-09 audit, 10x too low.
    """
    return _guan_frosin(T_GeV, theta_deg, a=0.0, b=1.0)


def _gaisser_tang(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """
    Tang et al. (2006), Phys. Rev. D 74, 053007, Sec. II.B: the modified
    Gaisser parametrisation, dΦ/dE [cm⁻²s⁻¹sr⁻¹GeV⁻¹],

        dN/dE dΩ = A · 0.14 E^-2.70 [1/(1 + 1.1 Ẽ cosθ*/115)
                                     + 0.054/(1 + 1.1 Ẽ cosθ*/850) + r_c]   (Eq. 3)

    with cosθ* from Eq. 10 (Guan's parametrisation, same P1..P5), in three
    segments of the surface energy E:
      E > 100/cosθ*          A = 1, Ẽ = E, r_c = 0 (plain Gaisser)
      1/cosθ* < E ≤ 100/cosθ*
                             r_c = 1e-4                                      (Eq. 4)
                             Ẽ = E + Δ, Δ = 2.06e-3 (950/cosθ* − 90) GeV    (Eqs. 5-6)
                             A = 1.1 (90 √(cosθ + 0.001) / 1030)^(4.5/(Ẽ cosθ*))
                             (plain cosθ inside the root)                   (Eq. 7)
      E ≤ 1/cosθ*            E → (3E + 7 secθ*)/10, then as above           (Eq. 9)
    Fitted to the world data set; the paper reports agreement within 40% in
    the worst case (E < 10 GeV, θ > 85°).

    Eq. 7 as printed has E, not Ẽ, in the exponent. Evaluated that way the
    formula falls orders of magnitude below the paper's own fitted curves in
    Fig. 1 (θ = 0° at 1 GeV: 2.7e-6 against ≈2e-3 for E^2.7 dN/dE dΩ;
    θ = 60° at 1 GeV: 1.6e-8 against ≈4e-4) and gives a vertical intensity of
    22 m⁻²s⁻¹sr⁻¹ above 1 GeV/c (PDG: ≈70). With Ẽ in the exponent it
    reproduces Fig. 1 at every angle shown (3.0e-3 and 3.8e-4 at those two
    points) and gives 60 m⁻²s⁻¹sr⁻¹. The printed E is taken to be a typo for
    Ẽ; see docs/FLUX_NORMALISATION_AUDIT.md. Before the 2026-09 audit this key
    held the Gaisser formula times an unexplained (1 + 0.054E/800) with a
    prefactor 10x low.
    """
    cos_th = math.cos(math.radians(theta_deg))
    cs = _guan_cos_star(cos_th)
    E0 = np.asarray(T_GeV, dtype=float) + M_MU_GEV
    E = np.where(E0 <= 1.0 / cs, (3.0 * E0 + 7.0 / cs) / 10.0, E0)       # Eq. 9
    high = E0 > 100.0 / cs
    delta = 2.06e-3 * (950.0 / cs - 90.0)                                  # Eq. 5
    E_t = np.where(high, E, E + delta)                                     # Eq. 6
    r_c = np.where(high, 0.0, 1.0e-4)                                      # Eq. 4
    A = np.where(high, 1.0,
                 1.1 * (90.0 * math.sqrt(max(cos_th, 0.0) + 0.001) / 1030.0)
                 ** (4.5 / (E_t * cs)))                 # Eq. 7, Ẽ: see above
    phi = A * 0.14 * E ** (-2.70) * (1.0 / (1.0 + 1.1 * E_t * cs / 115.0)
                                     + 0.054 / (1.0 + 1.1 * E_t * cs / 850.0)
                                     + r_c)                                # Eq. 3
    return np.where(np.isfinite(phi) & (phi > 0.0), phi, 0.0)


# CosmoALEPH: power-law fit 10^3.8467 p^-3.1952 [m^-2 s^-1 sr^-1 (GeV/c)^-1] to
# the vertical spectrum of Schmelling et al. (2013), Table 1, 112-2239 GeV/c.
# Vertical-only, returned isotropic as in the generator and UCMuGen.  Not
# valid below ~100 GeV/c.  Not offered in MODEL_LABELS (the GUI dropdowns);
# used by the backward MC for its spectrum 1.
def _cosmoaleph(T_GeV: np.ndarray, theta_deg: float) -> np.ndarray:
    """CosmoALEPH dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹]; isotropic, valid p ≳ 100 GeV/c."""
    E = np.asarray(T_GeV, dtype=float) + M_MU_GEV
    p = np.sqrt(np.maximum(E**2 - M_MU_GEV**2, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p > 0.0,
                        10.0**3.8467 * 1.0e-4 * p**-3.1952 * E / p, 0.0)


# ---------------------------------------------------------------------------
#  Guan / Frosin models — exact Python transcription of cosmoaleph_module_omp.f90
# ---------------------------------------------------------------------------
# cosθ* parameters (Guan 2015, Table 1 / arXiv:1509.06176)
_GUAN_P1    =  0.102573
_GUAN_P2    = -0.068287
_GUAN_P3    =  0.958633
_GUAN_P4    =  0.0407253
_GUAN_P5    =  0.817285
_GUAN_DENOM =  0.99144315   # = sqrt(1 + P1² + P2 + P4), Guan Eq. 2
_GUAN_EPI   =  115.0        # [GeV]  effective pion critical energy
_GUAN_EK    =  850.0        # [GeV]  effective kaon critical energy
_GUAN_KF    =    0.054      # kaon / (pion + kaon) fraction
_GUAN_PRE   =    0.14       # [cm⁻²s⁻¹sr⁻¹GeV⁻¹ at 1 GeV], Guan Eq. 1/3, = Fortran GUAN_PRE
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
    "reyna_bugaev": "Reyna (2006) ← recommended",
    "bugaev":       "Gaisser (1990), cosθ*  [E > 100/cosθ GeV only]",
    "gaisser_tang": "Tang et al. (2006) modified Gaisser",
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


def validity_warning(model: str, T_min_GeV: float, theta_deg: float = 0.0) -> str | None:
    """
    A one-line warning when `model` is evaluated below the energy its source
    fitted, or None.  Same limits as ucmugen::flux::validity_warnings:
    Reyna p > 1 GeV/c (Reyna 2006 Sec. 4); Guan/Frosin E > 1 GeV (Frosin 2025
    Sec. 3.2); plain Gaisser E > 100/cosθ GeV (PDG 2022 Eq. 30.4, Guan 2015
    Sec. 1); Tang et al. (2006) fitted 0.1 GeV-10 TeV, 0-87° (their Fig. 1),
    up to 40% off at E < 10 GeV and θ > 85°.  1% tolerance, so E_min = 1 GeV
    does not trip a 1 GeV/c limit.
    """
    E = T_min_GeV + M_MU_GEV
    p = math.sqrt(max(E * E - M_MU_GEV**2, 0.0))
    c = max(math.cos(math.radians(theta_deg)), 1e-3)
    if model == "reyna_bugaev" and p < 0.99:
        return f"Reyna (2006) is fitted for p > 1 GeV/c; p_min = {p:.3g} GeV/c."
    if model in ("guan_2015", "frosin_2025") and E < 0.99:
        return f"Guan/Frosin are fitted to data above 1 GeV; E_min = {E:.3g} GeV."
    if model == "gaisser_tang" and theta_deg > 85.0 and E < 10.0:
        return ("Tang et al. (2006) report up to 40% disagreement with data at "
                "E < 10 GeV and θ > 85° (Sec. II.B).")
    if model == "bugaev" and E < 0.99 * 100.0 / c:
        return (f"The Gaisser (1990) formula is valid only for E > 100/cosθ = "
                f"{100.0 / c:.0f} GeV; at E_min = {E:.3g} GeV it overestimates "
                f"the flux (x12 the PDG vertical integral above 1 GeV).")
    return None


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
    Beyond the table (3616 m w.e.) the range is continued with a + bE.
    """
    if opacity_gcm2 <= 0.0:
        return 0.0
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
        Flux model: "reyna_bugaev" (recommended), "guan_2015", "frosin_2025",
        "gaisser_tang" (Tang et al. 2006), "bugaev" (plain Gaisser, E > 100/cosθ).
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
