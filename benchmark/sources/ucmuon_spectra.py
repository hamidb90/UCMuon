"""Python ports of the UCMuon source-spectrum parametrizations.

Each function mirrors, line for line, the corresponding code in
src/generator/ucmuon_source_module.f90 so the figures produced from this
module show exactly what the generator samples.  All differential
intensities are returned in cm^-2 s^-1 sr^-1 (GeV/c)^-1 as a function of
muon momentum p [GeV/c] (the unit used by every reference data set).

Spectrum modes (build_cosmoaleph_cdf):
    1  CosmoALEPH power-law fit          (Schmelling et al. 2013)
    2  Power-law E^-3.7, shape only      (Kudryavtsev/MUSIC convention)
    3  PARMA/EXPACS                      (Sato 2015, 2016) -- numerical model;
                                         loaded from a cached CSV produced by
                                         make_parma_spectrum.py, not a closed form
    4  Guan et al. 2015                  (arXiv:1509.06176)
    5  Frosin et al. 2025 refit of Guan  (J. Phys. G 52, 035002)
    6  Gaisser 1990 pion+kaon formula    ("Bugaev/Gaisser" in the GUI)
    7  Reyna 2006 log-polynomial         (arXiv:hep-ph/0604145, Eq. 6-7)

Reference curves (not UCMuon modes):
    bugaev1998_vertical  Bugaev et al., PRD 58 (1998) 054001, Table II
"""

import os

import numpy as np

MUON_MASS = 0.10566          # GeV/c^2

_DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_PARMA_CSV = os.path.join(_DATADIR, "parma_vertical_spectrum.csv")
_PARMA_ANG_CSV = os.path.join(_DATADIR, "parma_angular.csv")
_parma_table = None          # lazily loaded (log p, log I_V) for interpolation
_parma_ang = None            # lazily loaded {p0: (theta_deg, F_ang)}

# --- Guan/Gaisser constants (ucmuon_source_module.f90 lines 34-48) ---------
GUAN_P = (0.102573, -0.068287, 0.958633, 0.0407253, 0.817285)
GUAN_DENOM = 0.99144315
GUAN_EPI, GUAN_EK, GUAN_KF = 115.0, 850.0, 0.054
GUAN_PRE, GUAN_IDX = 0.14, -2.7
GUAN_A, GUAN_B = 3.64, 1.29
FROSIN_A, FROSIN_B = 3.512, 1.388

# --- CosmoALEPH fit (lines 31-32); flux in m^-2 s^-1 sr^-1 (GeV/c)^-1 ------
A_COSMO, B_COSMO = 3.8467, -3.1952

# --- Reyna 2006 best fit (lines 54-58) --------------------------------------
REYNA_C = (0.00253, 0.2455, 1.288, -0.2555, 0.0209)


def guan_cos_star(cos_th):
    """Atmospheric curvature correction cos(theta*), Guan et al. Eq. 4."""
    cos_th = np.asarray(cos_th, dtype=float)
    p1, p2, p3, p4, p5 = GUAN_P
    numer = cos_th**2 + p1**2 + p2 * cos_th**p3 + p4 * cos_th**p5
    return np.clip(np.sqrt(np.maximum(0.0, numer)) / GUAN_DENOM, 0.0, 1.0)


def guan_flux(E, cos_th, a_par, b_par):
    """Modified Gaisser formula, dN/dE [cm^-2 s^-1 sr^-1 GeV^-1].

    a_par/b_par: (3.64, 1.29) Guan, (3.512, 1.388) Frosin, (0, 1) plain
    Gaisser 1990.
    """
    E = np.asarray(E, dtype=float)
    cs = guan_cos_star(cos_th)
    E_eff = E * (1.0 + a_par / (E * cs**b_par))
    pion = 1.0 / (1.0 + 1.1 * E * cs / GUAN_EPI)
    kaon = GUAN_KF / (1.0 + 1.1 * E * cs / GUAN_EK)
    return GUAN_PRE * E_eff**GUAN_IDX * (pion + kaon)


def _per_momentum(p, flux_per_energy):
    """dN/dp = dN/dE * dE/dp with dE/dp = p/E."""
    E = np.sqrt(np.asarray(p, dtype=float) ** 2 + MUON_MASS**2)
    return flux_per_energy * p / E


def mode1_cosmoaleph(p):
    """Mode 1: dN/dp = 10^3.8467 * p^-3.1952 m^-2 ... -> converted to cm^-2.

    Power-law fit to the CosmoALEPH measurement (valid 100 GeV/c - 2.5 TeV/c;
    the generator extrapolates it over the requested momentum window).
    """
    p = np.asarray(p, dtype=float)
    return 1.0e-4 * 10.0**A_COSMO * p**B_COSMO


def mode2_powerlaw(p, norm_to=None):
    """Mode 2: sampling density ~ E^-3.7 (Kudryavtsev/MUSIC convention).

    The generator never normalizes this shape to an absolute flux; for
    plotting, pass norm_to=(p0, I0) to pin the curve to a reference point.
    """
    p = np.asarray(p, dtype=float)
    E = np.sqrt(p**2 + MUON_MASS**2)
    shape = _per_momentum(p, E**-3.7)
    if norm_to is not None:
        p0, I0 = norm_to
        E0 = np.sqrt(p0**2 + MUON_MASS**2)
        shape = shape * I0 / (_per_momentum(np.array([p0]), E0**-3.7)[0])
    return shape


def _load_parma():
    """Load the cached PARMA vertical spectrum as (log10 p, log10 I_V)."""
    global _parma_table
    if _parma_table is None:
        if not os.path.exists(_PARMA_CSV):
            raise FileNotFoundError(
                f"{_PARMA_CSV} not found; run make_parma_spectrum.py first")
        arr = np.loadtxt(_PARMA_CSV, delimiter=",", comments="#")
        _parma_table = (np.log10(arr[:, 0]), np.log10(arr[:, 1]))
    return _parma_table


def mode3_parma(p):
    """Mode 3: PARMA/EXPACS vertical muon spectrum (Sato 2015, 2016).

    The only UCMuon source mode without a closed form: it is the numerical
    PARMA model (src/parma/parma_subroutines.f90).  This loads the cached
    vertical intensity I_V(p) = (mu+ + mu-) tabulated by
    make_parma_spectrum.py and log-log interpolates onto p.  Returns 0
    outside the tabulated momentum range (no extrapolation).
    Units: cm^-2 s^-1 sr^-1 (GeV/c)^-1, as a function of p [GeV/c].
    """
    logp, logiv = _load_parma()
    p = np.asarray(p, dtype=float)
    out = 10.0 ** np.interp(np.log10(np.maximum(p, 1e-300)), logp, logiv,
                            left=np.nan, right=np.nan)
    return np.where(np.isnan(out), 0.0, out)


def _load_parma_angular():
    """Load PARMA angular factors as {p0: (theta_deg, F_ang)} per momentum."""
    global _parma_ang
    if _parma_ang is None:
        if not os.path.exists(_PARMA_ANG_CSV):
            raise FileNotFoundError(
                f"{_PARMA_ANG_CSV} not found; run make_parma_spectrum.py first")
        arr = np.loadtxt(_PARMA_ANG_CSV, delimiter=",", comments="#")
        _parma_ang = {}
        for p0 in np.unique(arr[:, 0]):
            m = arr[:, 0] == p0
            _parma_ang[round(float(p0), 6)] = (arr[m, 1], arr[m, 2])
    return _parma_ang


def parma_angular(p0, theta_deg):
    """Mode 3 (PARMA) zenith dependence: angular factor F_ang(E(p0), theta).

    p0 must be one of the tabulated momenta in data/parma_angular.csv
    (1, 10, 100 GeV/c).  theta_deg may be a scalar or array.  The zenith
    figure plots parma_angular(p0, theta) / parma_angular(p0, 0).
    """
    table = _load_parma_angular()
    key = round(float(p0), 6)
    if key not in table:
        raise KeyError(f"PARMA angular data has no p0={p0} GeV/c; "
                       f"available: {sorted(table)}")
    th, ang = table[key]
    return np.interp(np.asarray(theta_deg, dtype=float), th, ang)


def mode4_guan(p, cos_th=1.0):
    """Mode 4: Guan et al. 2015, dN/dp [cm^-2 s^-1 sr^-1 (GeV/c)^-1]."""
    E = np.sqrt(np.asarray(p, dtype=float) ** 2 + MUON_MASS**2)
    return _per_momentum(p, guan_flux(E, cos_th, GUAN_A, GUAN_B))


def mode5_frosin(p, cos_th=1.0):
    """Mode 5: Frosin et al. 2025 refit (a=3.512, b=1.388) of the Guan model."""
    E = np.sqrt(np.asarray(p, dtype=float) ** 2 + MUON_MASS**2)
    return _per_momentum(p, guan_flux(E, cos_th, FROSIN_A, FROSIN_B))


def mode6_gaisser(p, cos_th=1.0):
    """Mode 6: plain Gaisser 1990 pion+kaon formula (a=0, b=1)."""
    E = np.sqrt(np.asarray(p, dtype=float) ** 2 + MUON_MASS**2)
    return _per_momentum(p, guan_flux(E, cos_th, 0.0, 1.0))


def mode7_reyna(p, cos_th=1.0):
    """Mode 7: Reyna 2006 vertical fit, with the Guan cos(theta*) substitution
    used by the generator (reyna_flux in the Fortran).  At cos_th=1 this is
    exactly Reyna Eq. 6-7: I_V(p) = c1 * p^-(c2 + c3*z + c4*z^2 + c5*z^3),
    z = log10(p).  Already per (GeV/c)."""
    p = np.asarray(p, dtype=float)
    c1, c2, c3, c4, c5 = REYNA_C
    p_eff = p * guan_cos_star(cos_th)
    z = np.log10(np.maximum(p_eff, 1e-300))
    n = c2 + c3 * z + c4 * z**2 + c5 * z**3
    return np.where(p_eff > 0, c1 * p_eff**-n, 0.0)


def reyna_angular(p, theta_rad):
    """Reyna's full angular prescription, Eq. 2 of the paper:
    I(p, theta) = cos^3(theta) * I_V(p * cos(theta)).
    This is what UCMuon reproduces by sampling cos^3 zenith angles
    (angular mode 5) together with spectrum mode 7."""
    c = np.cos(theta_rad)
    c1, c2, c3, c4, c5 = REYNA_C
    zeta = np.maximum(np.asarray(p, dtype=float) * c, 1e-300)
    z = np.log10(zeta)
    n = c2 + c3 * z + c4 * z**2 + c5 * z**3
    return c**3 * c1 * zeta**-n


def charge_ratio(p):
    """CosmoALEPH charge ratio step function (charge_ratio_from_p)."""
    edges = [112, 141, 178, 224, 282, 355, 447, 562, 708, 891, 1122, 1413, 1778]
    ratios = [1.252, 1.293, 1.259, 1.271, 1.239, 1.348, 1.541, 1.373, 1.243,
              1.547, 1.785, 1.361, 0.648, 1.495]
    p = np.asarray(p, dtype=float)
    return np.array(ratios)[np.searchsorted(edges, p, side="left")]


def bugaev1998_vertical(p):
    """Bugaev et al., PRD 58 (1998) 054001, Eq. 3.4 + Table II.

    Vertical sea-level momentum spectrum of conventional muons,
    cm^-2 s^-1 sr^-1 (GeV/c)^-1, fitted to 2% accuracy for p > 1 GeV/c.
    This is the parametrization Reyna 2006 reshaped; shown as an
    independent reference curve.
    """
    p = np.asarray(p, dtype=float)
    z = np.log10(p)
    out = np.empty_like(p)
    pieces = [
        (p < 927.65, 2.950e-3, (0.3061, 1.2743, -0.2630, 0.0252)),
        ((p >= 927.65) & (p < 1587.8), 1.781e-2, (1.7910, 0.3040, 0.0, 0.0)),
        ((p >= 1587.8) & (p < 4.1625e5), 14.35, (3.6720, 0.0, 0.0, 0.0)),
        (p >= 4.1625e5, 1.0e3, (4.0, 0.0, 0.0, 0.0)),
    ]
    for mask, C, (g0, g1, g2, g3) in pieces:
        n = g0 + g1 * z[mask] + g2 * z[mask] ** 2 + g3 * z[mask] ** 3
        out[mask] = C * p[mask] ** -n
    return out
