#!/usr/bin/env python3
"""
cosmoaleph_backward_mc.py  —  UCLouvain Muography Group
Self-developed backward Monte Carlo muon flux integrator.

Purpose
───────
Compute the expected muon flux and count rate at a detector sited at
depth d [m] without requiring a surface-generated muon file.
Instead of tracking forward from the atmosphere, we trace backward from
the detector using the CSDA inversion relation and the known surface spectrum.

Algorithm
─────────
Expected differential flux at the detector:

    dΦ_det/dE_det dΩ = Φ_surf(E_s(E_det,θ), θ) × |dE_s/dE_det| × P_surv(E_s, X_slant)

where
    E_s(E_det, θ) = CSDA backward inversion: surface energy that reaches E_det
                    after traversing X_slant = depth × rho × 100 / cosθ [g/cm²]
    |dE_s/dE_det| = dE/dx(E_s) / dE/dx(E_det)  (Jacobian of the energy mapping)
    P_surv        = survival probability including stochastic correction

Integrated expected rate [muons m⁻² s⁻¹]:
    Rate = ∫dE_det ∫dΩ dΦ_det/dE_det cosθ

Stochastic correction for P_surv
─────────────────────────────────
CSDA gives P_surv = 1 for all muons above threshold. The stochastic correction
accounts for muons that are stopped by a single catastrophic radiative event
even though their mean energy loss (CSDA) would let them survive.

Poisson approximation:
    P(no fatal catastrophic event) = exp(−λ_fatal × X_slant)
where λ_fatal = b_total × ln(1/v_stop) / ln(1/v_cut)  [events/g/cm²]
and v_stop is the minimum energy fraction that would stop the muon.

Spectrum models (backward-MC numbering, NOT the generator's 1-8; e.g. 3 is
Guan here, PARMA in the generator):
    1 = CosmoALEPH (Schmelling 2013)  dN/dp ∝ p^{-3.195}
    2 = Power-law                   dN/dE ∝ E^{-3.7}
    3 = Guan et al. (2015)          arXiv:1509.06176
    4 = Frosin et al. (2025)        J. Phys. G 52, 035002

Stdin input:
     1  depth_m        overburden depth [m]
     2  rho            density [g/cm³]
     3  X0_cm          radiation length [cm]
     4  mat_id         1=rock 2=water 3=seawater 4=iron
     5  spectrum_mode  1–4 (see above)
     6  E_min_GeV      minimum detector energy
     7  E_max_GeV      maximum detector energy
     8  theta_max_deg  maximum zenith angle
     9  n_E            energy grid bins
    10  n_theta        zenith angle bins
    11  mode           0=CSDA only  1=CSDA + stochastic P_surv correction
    12  v_cut          catastrophic threshold (default 0.05)
    13  b_rad          radiative coefficient override [cm2/g] (0=use database)
    14  outfile        output ASCII file

Author: Hamid Basiri <hamid.basiri@uclouvain.be>
"""

import contextlib
import sys
import warnings
import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants
# ─────────────────────────────────────────────────────────────────────────────
M_MU     = 105.6584         # muon rest mass [MeV]
M_MU_GEV = 0.1056584

# ─────────────────────────────────────────────────────────────────────────────
# Groom (2001) Standard Rock CSDA range table (33 entries, same as GUI)
# ─────────────────────────────────────────────────────────────────────────────
_GROOM_T = np.array([
    1.0e+01, 1.4e+01, 2.0e+01, 3.0e+01, 4.0e+01,
    8.0e+01, 1.0e+02, 1.4e+02, 2.0e+02, 3.0e+02,
    4.0e+02, 8.0e+02, 1.0e+03, 1.4e+03, 2.0e+03,
    3.0e+03, 4.0e+03, 8.0e+03, 1.0e+04, 1.4e+04,
    2.0e+04, 3.0e+04, 4.0e+04, 8.0e+04, 1.0e+05,
    1.4e+05, 2.0e+05, 3.0e+05, 4.0e+05, 8.0e+05,
    1.0e+06, 1.4e+06, 2.0e+06,
])
_GROOM_R = np.array([
    8.516e-01, 1.542e+00, 2.866e+00, 5.698e+00, 9.145e+00,
    2.676e+01, 3.696e+01, 5.879e+01, 9.332e+01, 1.524e+02,
    2.115e+02, 4.418e+02, 5.534e+02, 7.712e+02, 1.088e+03,
    1.599e+03, 2.095e+03, 3.998e+03, 4.920e+03, 6.724e+03,
    9.360e+03, 1.362e+04, 1.776e+04, 3.343e+04, 4.084e+04,
    5.495e+04, 7.459e+04, 1.040e+05, 1.302e+05, 2.129e+05,
    2.453e+05, 2.990e+05, 3.616e+05,
])

# Above 2 TeV the table is continued to 10 TeV with the PDG 2024 standard rock
# stopping power (gui/ucmuon_stochastic_driver.py, the table UCMuon-MC uses; its
# range agrees with Groom's to 0.03 % up to 2 TeV), integrated from Groom's
# last point so the range is continuous. Up to 1.3.2 the table ended at 2 TeV
# (3616 m w.e.): under thick rock that cut the flux (at Vesuvius, 106 of the
# 462 density-sensitive directions lost more than 1 %, up to 33 %).
def _pdg_rock_continuation(T_max=1.0e7):
    gui_dir = str(Path(__file__).resolve().parent)
    if gui_dir not in sys.path:
        sys.path.insert(0, gui_dir)
    import ucmuon_stochastic_driver as _usd
    t, T = _usd._PROC_TABLES["rock"], _usd._PROC_T
    dedx = t["ION"] + t["BREMS"] + t["PAIR"] + t["NUCL"]          # MeV cm2/g
    sel = (T >= _GROOM_T[-1]) & (T <= T_max * 1.000001)
    Ts, inv = T[sel], 1.0 / dedx[sel]
    R = _GROOM_R[-1] + np.concatenate(([0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(Ts))))
    return Ts[1:], R[1:]


_EXT_T, _EXT_R = _pdg_rock_continuation()
_TAB_T = np.concatenate([_GROOM_T, _EXT_T])       # 10 MeV - 10 TeV
_TAB_R = np.concatenate([_GROOM_R, _EXT_R])

# Pre-compute dE/dx fine grid once at import: the 1000 points up to 2 TeV as
# before, continued above it with the same logarithmic step.
_N_FINE   = 1000
_T_FINE   = np.logspace(np.log10(_GROOM_T[0]), np.log10(_GROOM_T[-1]), _N_FINE)
_LSTEP    = np.log10(_T_FINE[1] / _T_FINE[0])
_T_FINE   = np.concatenate([_T_FINE, 10.0 ** np.linspace(
    np.log10(_GROOM_T[-1]), np.log10(_TAB_T[-1]),
    int(np.ceil(np.log10(_TAB_T[-1] / _GROOM_T[-1]) / _LSTEP)) + 1)[1:]])
_R_FINE   = np.exp(np.interp(np.log(_T_FINE), np.log(_TAB_T), np.log(_TAB_R)))
_DRDT     = np.gradient(_R_FINE, _T_FINE)
_DEDX_STD = 1.0 / np.maximum(_DRDT, 1e-12)        # [MeV cm2/g] Standard Rock


def _dedx(E_MeV, tab=None):
    """Stopping power [MeV cm2/g] of material table `tab` (None = Standard
    Rock, Groom 2001)."""
    T, dedx = (_T_FINE, _DEDX_STD) if tab is None else (tab["T"], tab["dedx"])
    E = np.clip(np.asarray(E_MeV, float), T[0], T[-1])
    return np.interp(E, T, dedx)


def _range(E_MeV, tab=None):
    """CSDA range [g/cm2] in material table `tab` (None = Standard Rock)."""
    T, R = (_T_FINE, _R_FINE) if tab is None else (tab["T"], tab["R"])
    E = np.clip(np.asarray(E_MeV, float), T[0], T[-1])
    return np.interp(E, T, R)


# ─────────────────────────────────────────────────────────────────────────────
# End of the range table
#
# The table stops at 10 TeV, a range of 6889 m w.e. in standard rock (2 TeV and
# 3616 m w.e. up to 1.3.2). A detector energy whose surface energy would lie
# beyond it cannot be mapped back and is left out of the integral, so under
# thick rock the flux is a lower bound, and exactly zero once every detector
# energy is beyond the table. Up to v1.1.2 this was silent. Both integrators now
# estimate what is left out, by continuing the table with dE/dX = a + bE fitted
# to its last four points, and raise a RangeTableWarning when that is more than
# 1 % of the flux or the flux is zero. The estimate only sizes the warning; the
# flux returned is unchanged. (Guan, 1-2500 GeV at the detector, θ = 0-60°: the
# cut is below 1 % up to 3500 m w.e., 1.3 % at 4000, 7 % at 5000, 30 % at 6000
# and 59 % at 6500; with the 2 TeV table it was 0.6 % at 1000 m w.e., 6 % at
# 2000 and 37 % at 3000.) Detector energies above 10 TeV are outside the table
# altogether and are left out under any rock.
# ─────────────────────────────────────────────────────────────────────────────
RANGE_TABLE_E_MAX_GEV  = float(_TAB_T[-1]) / 1000.0       # 10 000 GeV
RANGE_TABLE_R_MAX_GCM2 = float(_TAB_R[-1])                # 688 900 g/cm2
RANGE_WARN_FRACTION    = 0.01

_EXT_B, _EXT_A = np.polyfit(_TAB_T[-4:], 1.0 / np.gradient(_TAB_R, _TAB_T)[-4:], 1)


class RangeTableWarning(UserWarning):
    """Part of a flux was left out because it lies beyond the range table."""

    def __init__(self, note, X_gcm2=0.0, fraction=0.0, zero=False):
        super().__init__(note)
        self.X_gcm2, self.fraction, self.zero = X_gcm2, fraction, zero


def _beyond_table_contribution(E_det_MeV, X_gcm2, theta_rad, spectrum_mode,
                               tab, b_total, v_cut, mode):
    """φ_det [m⁻² s⁻¹ sr⁻¹ GeV⁻¹] of a detector energy whose surface energy is
    beyond the table, with the table continued as dE/dX = a + bE. For the
    warning only."""
    R_det = float(_range(E_det_MeV, tab))
    R_surf = R_det + X_gcm2
    if tab is None:
        a, b, T0, R0 = _EXT_A, _EXT_B, float(_T_FINE[-1]), float(_R_FINE[-1])
    else:
        a, b, T0, R0 = tab["ext_a"], tab["ext_b"], float(tab["T"][-1]), float(tab["R"][-1])
    E_s = ((a + b * T0) * np.exp(b * (R_surf - R0)) - a) / b
    J = (a + b * E_s) / max(float(_dedx(E_det_MeV, tab)), 1e-12)
    Ps = 1.0
    if mode == 1:
        v_stop = float(np.clip(1.0 - R_det / R_surf, v_cut, 1.0 - 1e-6))
        Ps = float(np.exp(-b_total * np.log(1.0 / v_stop) / np.log(1.0 / v_cut) * X_gcm2))
    return float(_flux_surface(E_s / 1000.0, theta_rad, spectrum_mode)) * J * Ps


def range_table_note(X_gcm2, flux, missing):
    """The warning text, or None when what is left out is below 1 % of the flux."""
    if missing <= 0.0:
        return None
    total = flux + missing
    frac = missing / total
    if flux > 0.0 and frac <= RANGE_WARN_FRACTION:
        return None
    head = (f"Backward MC: at a slant opacity of {X_gcm2 / 100.0:.0f} m w.e., muons "
            f"that need a surface energy above {RANGE_TABLE_E_MAX_GEV / 1000.0:g} TeV "
            f"(the end of the range table) are left out")
    if flux <= 0.0:
        return (head + ", so the flux returned is 0, which is not a result: the "
                "true flux is small but not zero.")
    return (head + f", about {100.0 * frac:.0f} % of the flux (estimated by "
            f"continuing the table); the flux returned is a lower bound.")


def _warn_off_table(X_gcm2, flux, missing):
    note = range_table_note(X_gcm2, flux, missing)
    if note:
        warnings.warn(RangeTableWarning(note, X_gcm2, missing / (flux + missing),
                                        flux <= 0.0), stacklevel=3)
    return note


@contextlib.contextmanager
def collect_range_warnings():
    """Gather the RangeTableWarnings raised inside the block into a list instead
    of printing one per call (a flux map makes thousands); summarise them with
    range_summary(). Other warnings are passed on."""
    got = []
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        yield got
    for w in rec:
        # by name: callers load this module with importlib, one copy each
        if w.category.__name__ == "RangeTableWarning":
            got.append(w.message)
        else:
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)


def range_summary(got, n_directions):
    """One line for a map or a scan, from collect_range_warnings(); None if empty."""
    if not got:
        return None
    n_zero = sum(m.zero for m in got)
    X_min = min(m.X_gcm2 for m in got)
    worst = max((m.fraction for m in got if not m.zero), default=0.0)
    return (f"Backward MC: in {len(got)} of {n_directions} directions (slant opacity "
            f"from {X_min / 100.0:.0f} m w.e.) more than 1 % of the flux needs surface "
            f"energies above {RANGE_TABLE_E_MAX_GEV / 1000.0:g} TeV, the end of the "
            f"range table, and is left out: those fluxes are lower bounds"
            + (f" (up to {100.0 * worst:.0f} % low)" if worst > 0 else "")
            + (f", and {n_zero} are 0, where the true flux is small but not zero."
               if n_zero else "."))


SPECTRUM_NAMES = {1: "CosmoALEPH", 2: "power law", 3: "Guan 2015", 4: "Frosin 2025"}


def spectrum_name(spectrum_mode):
    """Name of a backward-MC spectrum number (its numbering differs from the
    generator's, so outputs print the name)."""
    return SPECTRUM_NAMES.get(int(spectrum_mode), f"spectrum {spectrum_mode}")


def spectrum_range_warning(spectrum_mode, E_min_GeV):
    """A one-line warning when E_min [GeV, kinetic] is below the fitted range of
    the surface spectrum, or when the spectrum has no absolute normalisation;
    None otherwise. Same limits as ucmugen::flux::validity_warnings."""
    if spectrum_mode == 1:
        p = float(np.sqrt((E_min_GeV + M_MU_GEV) ** 2 - M_MU_GEV ** 2))
        if p < 99.0:
            return (f"CosmoALEPH is a fit to vertical data at 112-2239 GeV/c; from "
                    f"p_min = {p:.3g} GeV/c it extrapolates (x46 in the vertical "
                    f"integral above 1 GeV/c), and so does the open-sky flux that "
                    f"every T = Phi_rock / Phi_sky is divided by.")
    elif spectrum_mode == 2:
        return ("The E^-3.7 power law has no absolute normalisation: transmissions "
                "are meaningful, absolute fluxes and rates are not.")
    elif spectrum_mode in (3, 4):
        return _ffe.validity_warning("guan_2015" if spectrum_mode == 3 else "frosin_2025",
                                     E_min_GeV)
    return None


def _backward_energy(E_det_MeV, X_gcm2, tab=None):
    """
    Backward CSDA mapping: given E_det_MeV at detector, find surface energy E_s
    and Jacobian |dE_s/dE_det| = dE/dx(E_s) / dE/dx(E_det), in material table
    `tab` (None = Standard Rock).
    Returns (E_s_MeV, jacobian).  E_s = inf if unreachable (off-table).
    """
    T, R = (_T_FINE, _R_FINE) if tab is None else (tab["T"], tab["R"])
    R_det  = _range(E_det_MeV, tab)
    R_surf = R_det + X_gcm2
    if np.ndim(R_surf) == 0:
        if R_surf > R[-1]:
            return np.inf, 1.0
        E_s = float(np.interp(R_surf, R, T))
        J   = float(_dedx(E_s, tab) / max(_dedx(E_det_MeV, tab), 1e-12))
        return E_s, J
    # Vectorised
    E_s  = np.full_like(np.asarray(E_det_MeV, float), np.inf)
    J    = np.ones_like(E_s)
    ok   = R_surf <= R[-1]
    E_s[ok]  = np.interp(R_surf[ok], R, T)
    J[ok]    = _dedx(E_s[ok], tab) / np.maximum(_dedx(E_det_MeV[ok], tab), 1e-12)
    return E_s, J


# ─────────────────────────────────────────────────────────────────────────────
# Material database
# ─────────────────────────────────────────────────────────────────────────────
# Standard Rock uses the Groom (2001) table above. The other materials use
# their PDG 2024 tables (gui/ucmuon_stochastic_driver.py, as UCMuon-MC and the
# Python Bethe-Bloch do); seawater is water with the same composition
# rescaling as UCMuon-MC. Up to 1.3.2 they rescaled the rock table by one
# factor a (water 1.046, iron 0.930) and multiplied the RANGE by it where it
# should divide, so a material acted like X/a of rock instead of a*X: water
# gave rates 11, 29 and 42 % too high at 10, 50 and 200 m, iron 13 to 27 %
# too low. Their b for P_surv (water 3.20e-6, iron 4.06e-6) is now the rock
# value scaled by the tables' radiative ratio at 1 TeV (0.75 and 1.87).
_MAT_DB = {
    1: {"name": "Standard Rock", "b_rad": 3.475e-6, "X0_cm": 26.48},
    2: {"name": "Water/Ice",     "proc": "water", "X0_cm": 36.08},
    3: {"name": "Seawater",      "proc": "water", "X0_cm": 35.75,
        "s_ion": 0.55250 / 0.55509, "s_rad": 3.22e-6 / 3.20e-6},
    4: {"name": "Iron",          "proc": "iron",  "X0_cm": 1.757},
}
_TABLES = {}


def _material(mat_id):
    """(material dict, range/stopping-power table or None for Standard Rock,
    b for P_surv) of a material id."""
    mat = _MAT_DB.get(mat_id, _MAT_DB[1])
    key = mat.get("proc")
    if key is None:
        return mat, None, mat["b_rad"]
    if mat_id not in _TABLES:
        import ucmuon_stochastic_driver as _usd
        t, T = _usd._PROC_TABLES, _usd._PROC_T
        s_ion, s_rad = mat.get("s_ion", 1.0), mat.get("s_rad", 1.0)
        rad = lambda k: t[k]["BREMS"] + t[k]["PAIR"] + t[k]["NUCL"]
        dedx = t[key]["ION"] * s_ion + rad(key) * s_rad
        d = np.exp(np.interp(np.log(_T_FINE), np.log(T), np.log(dedx)))
        inv = 1.0 / d
        R = t[key]["R10"] / s_ion + np.concatenate(
            ([0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(_T_FINE))))
        top = (T >= 2.0e6) & (T <= _T_FINE[-1] * 1.000001)   # 2 - 10 TeV
        ext_b, ext_a = np.polyfit(T[top], dedx[top], 1)
        lt = np.log(1.0e6)                                    # 1 TeV
        b_ratio = (np.interp(lt, np.log(T), rad(key)) * s_rad
                   / np.interp(lt, np.log(T), rad("rock")))
        _TABLES[mat_id] = (dict(T=_T_FINE, R=R, dedx=d, ext_a=ext_a, ext_b=ext_b),
                           _MAT_DB[1]["b_rad"] * b_ratio)
    tab, b = _TABLES[mat_id]
    return mat, tab, b


# ─────────────────────────────────────────────────────────────────────────────
# Surface flux models
#
# Taken from fast_flux_estimator so that the backward MC, the fast estimator
# and the generator share one audited implementation of each spectrum (see
# docs/FLUX_NORMALISATION_AUDIT.md).  Before the 2026-09 audit this module had
# its own formulas, none of which matched its paper: mode 3 "Guan" was built
# from Reyna's constants, and the CosmoALEPH normalisation was 27x low at
# 100 GeV.
# ─────────────────────────────────────────────────────────────────────────────
_gui_dir = str(Path(__file__).resolve().parent)
if _gui_dir not in sys.path:
    sys.path.insert(0, _gui_dir)
import fast_flux_estimator as _ffe                      # noqa: E402


def _flux_surface(E_GeV, theta_rad, mode):
    """
    Differential surface intensity Φ(T, θ) [muons m⁻² s⁻¹ GeV⁻¹ sr⁻¹] at
    kinetic energy T = E_GeV (the backward CSDA map works in kinetic energy).

    mode 1  CosmoALEPH power-law fit (vertical-only, isotropic as in the
            generator); valid p ≳ 100 GeV/c
    mode 2  power law ∝ E^-3.7: a shape with NO absolute normalisation; the
            95 m⁻²s⁻¹sr⁻¹GeV⁻¹ prefactor is arbitrary and so is any rate or
            exposure time computed from it
    mode 3  Guan et al. (2015) Eq. 3          (generator spectrum 4)
    mode 4  Frosin et al. (2025) Table 4 fit  (generator spectrum 5)
    """
    E = np.asarray(E_GeV, float)
    theta_deg = float(np.degrees(theta_rad))
    if mode == 1:
        return 1.0e4 * _ffe._cosmoaleph(E, theta_deg)
    if mode == 2:
        return 95.0 * E**(-3.7) * np.cos(theta_rad)**2
    if mode == 3:
        return 1.0e4 * _ffe._guan_2015(E, theta_deg)
    if mode == 4:
        return 1.0e4 * _ffe._frosin_2025(E, theta_deg)
    return np.zeros_like(E)


# ─────────────────────────────────────────────────────────────────────────────
# Stochastic survival probability (Poisson analytical approximation)
# ─────────────────────────────────────────────────────────────────────────────

def _P_surv_stochastic(E_s_MeV, E_det_MeV, X_gcm2, b_total, v_cut, tab=None):
    """
    Survival probability correction for catastrophic radiative events.

    P_surv = exp(−λ_fatal × X_slant)
    where λ_fatal [events/g/cm²] = b_total × ln(1/v_stop) / ln(1/v_cut)
    and v_stop = min fraction that would stop the muon:
        v_stop = 1 − R_det / R_surf  (fraction of E_s that can be lost)

    Reduces to 1 in CSDA limit (v_stop → 0) and to 0 when every catastrophic
    event above v_cut would stop the muon.
    """
    R_det  = float(_range(E_det_MeV, tab))
    R_surf = float(_range(E_s_MeV,   tab))
    if R_surf <= 0:
        return 0.0
    v_stop = 1.0 - R_det / R_surf
    v_stop = float(np.clip(v_stop, v_cut, 1.0 - 1e-6))
    ln_ivc    = np.log(1.0 / v_cut)
    ln_ivstop = np.log(1.0 / v_stop)
    lam_fatal = b_total * ln_ivstop / ln_ivc
    return float(np.exp(-lam_fatal * X_gcm2))


# ─────────────────────────────────────────────────────────────────────────────
# Backward MC flux integration
# ─────────────────────────────────────────────────────────────────────────────

def backward_mc_flux(depth_m, rho, mat_id, spectrum_mode,
                     E_min_GeV=0.1, E_max_GeV=5000.0,
                     theta_max_deg=70.0, n_E=80, n_theta=30,
                     mode=1, v_cut=0.05, b_rad_override=0.0,
                     progress_cb=None):
    """
    Compute expected muon flux at depth by backward CSDA integration.

    Returns dict with:
        E_det_GeV    : detector energy grid [GeV]
        flux_det     : dΦ/dE integrated over θ [m⁻² s⁻¹ GeV⁻¹]
        flux_surf    : surface reference flux [m⁻² s⁻¹ GeV⁻¹]
        E_surf_GeV   : mean required surface energy per E_det bin [GeV]
        P_survival   : mean stochastic survival probability
        rate_m2_s    : total expected rate [muons m⁻² s⁻¹]
        X_vert_gcm2  : vertical opacity [g/cm²]
        info         : metadata dict
    """
    mat, tab, b_mat = _material(mat_id)
    b_total  = b_rad_override if b_rad_override > 0 else b_mat

    X_vert   = depth_m * 100.0 * rho           # vertical opacity [g/cm²]

    # Zenith angle grid
    th_edges = np.linspace(0.0, np.deg2rad(theta_max_deg), n_theta + 1)
    th_mid   = 0.5 * (th_edges[:-1] + th_edges[1:])
    dth      = th_edges[1:] - th_edges[:-1]
    dOmega   = 2.0 * np.pi * np.sin(th_mid) * dth     # solid angle [sr]
    cos_mid  = np.cos(th_mid)

    # Slant opacity per zenith bin [g/cm²]
    X_slant  = X_vert / np.maximum(cos_mid, 0.02)

    # Detector energy grid [MeV]
    E_det_MeV = np.logspace(np.log10(E_min_GeV * 1000.0),
                             np.log10(E_max_GeV * 1000.0), n_E)
    # Actual node spacing of the log grid: dE_i = E_i * dlnE (trapezoid
    # end-weights), so the integrated rate is independent of n_E.
    # (A hard-coded 5% bin width made the rate scale linearly with n_E.)
    dlnE = (np.log(E_max_GeV) - np.log(E_min_GeV)) / max(n_E - 1, 1)

    # Output arrays
    flux_det   = np.zeros(n_E)     # dΦ/dE_det integrated over Ω [m⁻² s⁻¹ GeV⁻¹]
    flux_surf  = np.zeros(n_E)     # surface reference dΦ/dE [m⁻² s⁻¹ GeV⁻¹]
    E_surf_out = np.zeros(n_E)     # mean surface energy [GeV]
    Ps_out     = np.zeros(n_E)     # mean P_survival
    rate_total = 0.0               # [m⁻² s⁻¹]

    dOmega_total = float(dOmega.sum())
    missing_rate = 0.0             # left out beyond the range table [m⁻² s⁻¹]

    for iE, E_det in enumerate(E_det_MeV):
        w_trap = 0.5 if iE in (0, n_E - 1) else 1.0
        dE_GeV = (E_det / 1000.0) * dlnE * w_trap   # energy bin width [GeV]

        sum_flux = 0.0; sum_surf = 0.0
        sum_Es   = 0.0; sum_Ps  = 0.0
        sum_dO   = 0.0      # solid angle of the zenith bins that contribute

        for ith, (th, dO, cos_t, X_sl) in enumerate(zip(th_mid, dOmega, cos_mid, X_slant)):

            # Backward CSDA: what surface energy E_s is required?
            E_s, Jacob = _backward_energy(E_det, X_sl, tab)
            if E_s == np.inf or np.isnan(E_s):
                if E_det <= _T_FINE[-1]:       # surface energy beyond the table
                    missing_rate += _beyond_table_contribution(
                        E_det, X_sl, th, spectrum_mode, tab, b_total, v_cut,
                        mode) * dO * cos_t * dE_GeV
                continue
            E_s_GeV = E_s / 1000.0
            if E_s_GeV < E_min_GeV or E_s_GeV > E_max_GeV:
                continue

            # Surface flux
            phi_s = float(_flux_surface(E_s_GeV, th, spectrum_mode))
            if phi_s <= 0:
                continue

            # Survival probability
            if mode == 1:
                Ps = _P_surv_stochastic(E_s, E_det, X_sl, b_total, v_cut, tab)
            else:
                Ps = 1.0

            # Detector flux contribution  [m⁻² s⁻¹ GeV⁻¹ sr⁻¹ → m⁻² s⁻¹ GeV⁻¹]
            phi_det  = phi_s * Jacob * Ps       # [m⁻² s⁻¹ GeV⁻¹ sr⁻¹]
            contrib  = phi_det * dO * cos_t     # [m⁻² s⁻¹ GeV⁻¹]

            sum_flux += contrib
            sum_surf += phi_s * dO * cos_t
            sum_Es   += E_s_GeV * dO
            sum_Ps   += Ps * dO
            sum_dO   += dO

            rate_total += contrib * dE_GeV

        flux_det[iE]   = sum_flux
        flux_surf[iE]  = sum_surf
        # Means over the zenith bins that contribute: directions beyond the
        # range table or outside the energy window carry no E_s or P_surv
        # (up to 1.3.0 they counted as 0 through the full solid angle).
        E_surf_out[iE] = sum_Es  / sum_dO if sum_dO > 0 else np.nan
        Ps_out[iE]     = sum_Ps  / sum_dO if sum_dO > 0 else np.nan

        if progress_cb and iE % max(1, n_E // 10) == 0:
            progress_cb(f"E bin {iE+1}/{n_E}  E_det={E_det/1000:.2f} GeV"
                        f"  rate so far {rate_total:.3e} m-2 s-1")

    range_note = _warn_off_table(float(X_slant.max()), rate_total, missing_rate)

    return dict(
        E_det_GeV   = E_det_MeV / 1000.0,
        flux_det    = flux_det,
        flux_surf   = flux_surf,
        E_surf_GeV  = E_surf_out,
        P_survival  = Ps_out,
        rate_m2_s   = rate_total,
        X_vert_gcm2 = X_vert,
        info        = dict(depth_m=depth_m, rho=rho, mat=mat["name"],
                           spectrum_mode=spectrum_mode,
                           n_E=n_E, n_theta=n_theta,
                           theta_max_deg=theta_max_deg,
                           mode="CSDA+stochastic" if mode == 1 else "CSDA only",
                           v_cut=v_cut, range_note=range_note),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Single-direction integrated flux (used by the terrain flux map)
# ─────────────────────────────────────────────────────────────────────────────

def directional_flux(X_slant_gcm2, theta_rad, spectrum_mode,
                     E_min_GeV=0.5, E_max_GeV=5000.0, n_E=40,
                     mode=1, v_cut=0.05, mat_id=1, b_rad_override=0.0):
    """
    Through-rock integrated directional flux [m⁻² s⁻¹ sr⁻¹] at ONE zenith
    angle and ONE slant opacity:

        Φ(θ, X) = ∫_{E_min}^{E_max} dE_det
                  Φ_surf(E_s(E_det, X), θ) · |dE_s/dE_det| · P_surv

    The slant opacity X is used as-is (no vertical-equivalent conversion),
    and the surface flux is evaluated at the exact zenith angle θ.  The
    required surface energy E_s may exceed E_max — the spectrum
    parametrisations remain valid there, and truncating would bias the
    through-rock integral low.

    With X = 0 this reduces to ∫ Φ_surf dE over the same grid, so the
    numerator and denominator of a transmission map T = Φ(X)/Φ(0) use
    identical integration by construction.
    """
    mat, tab, b_mat = _material(mat_id)
    b_total = b_rad_override if b_rad_override > 0 else b_mat

    E_det_MeV = np.logspace(np.log10(E_min_GeV * 1000.0),
                             np.log10(E_max_GeV * 1000.0), n_E)
    dlnE = (np.log(E_max_GeV) - np.log(E_min_GeV)) / max(n_E - 1, 1)

    flux = 0.0
    missing = 0.0                      # left out beyond the range table
    for iE, E_det in enumerate(E_det_MeV):
        w_trap = 0.5 if iE in (0, n_E - 1) else 1.0
        dE_GeV = (E_det / 1000.0) * dlnE * w_trap
        if X_slant_gcm2 > 0.0:
            E_s, Jacob = _backward_energy(E_det, X_slant_gcm2, tab)
            if not np.isfinite(E_s):
                if E_det <= _T_FINE[-1]:   # surface energy beyond the table
                    missing += _beyond_table_contribution(
                        E_det, X_slant_gcm2, theta_rad, spectrum_mode, tab,
                        b_total, v_cut, mode) * dE_GeV
                continue
            Ps = (_P_surv_stochastic(E_s, E_det, X_slant_gcm2,
                                     b_total, v_cut, tab)
                  if mode == 1 else 1.0)
        else:
            E_s, Jacob, Ps = E_det, 1.0, 1.0
        phi_s = float(_flux_surface(E_s / 1000.0, theta_rad, spectrum_mode))
        if phi_s > 0.0:
            flux += phi_s * Jacob * Ps * dE_GeV
    _warn_off_table(X_slant_gcm2, flux, missing)
    return flux


# ─────────────────────────────────────────────────────────────────────────────
# Output writer
# ─────────────────────────────────────────────────────────────────────────────

def _write_results(res, fpath):
    info = res["info"]
    with open(fpath, "w") as fh:
        fh.write("# UCMuon Backward MC — UCLouvain Muography Group\n")
        fh.write(f"# Depth: {info['depth_m']:.1f} m  "
                 f"rho: {info['rho']:.3f} g/cm3  mat: {info['mat']}\n")
        fh.write(f"# Vertical opacity: {res['X_vert_gcm2']:.1f} g/cm2\n")
        fh.write(f"# Mode: {info['mode']}  spectrum: {spectrum_name(info['spectrum_mode'])}"
                 f" (backward-MC number {info['spectrum_mode']})\n")
        fh.write(f"# Total expected rate: {res['rate_m2_s']:.5e} m-2 s-1\n")
        if info.get("range_note"):
            fh.write(f"# WARNING: {info['range_note']}\n")
        fh.write("# Cols: E_det[GeV]  dPhi/dE[m-2 s-1 GeV-1]  E_surf[GeV]"
                 "  P_survival  dPhi_surf/dE[m-2 s-1 GeV-1]\n")
        for i in range(len(res["E_det_GeV"])):
            fh.write(f"{res['E_det_GeV'][i]:12.4f}"
                     f"  {res['flux_det'][i]:14.6e}"
                     f"  {res['E_surf_GeV'][i]:12.4f}"
                     f"  {res['P_survival'][i]:10.6f}"
                     f"  {res['flux_surf'][i]:14.6e}\n")


def _print_summary(res):
    info = res["info"]
    print("", flush=True)
    print("  ══════════════════════════════════════════════", flush=True)
    print("   UCMuon Backward MC — Flux at depth", flush=True)
    print("  ══════════════════════════════════════════════", flush=True)
    print(f"   Depth    : {info['depth_m']:.1f} m  ({info['mat']})", flush=True)
    print(f"   Opacity  : {res['X_vert_gcm2']:.1f} g/cm2 (vertical)", flush=True)
    print(f"   Mode     : {info['mode']}", flush=True)
    print(f"   Rate     : {res['rate_m2_s']:.4e} m-2 s-1", flush=True)
    imax = int(np.argmax(res["flux_det"]))
    print(f"   Peak flux: E_det = {res['E_det_GeV'][imax]:.2f} GeV", flush=True)
    if res["E_surf_GeV"][imax] > 0:
        print(f"   E_surf at peak: {res['E_surf_GeV'][imax]:.2f} GeV", flush=True)
    alive_ps = res["P_survival"][res["P_survival"] > 1e-4]
    if len(alive_ps):
        print(f"   Mean P_surv: {float(np.mean(alive_ps)):.4f}", flush=True)
    if info.get("range_note"):
        print(f"   WARNING: {info['range_note']}", flush=True)
    print("  ══════════════════════════════════════════════", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    lines = [ln.strip() for ln in sys.stdin if ln.strip() and not ln.strip().startswith("#")]

    def _rd(i, default, typ=str):
        try:    return typ(lines[i]) if i < len(lines) else default
        except: return default

    depth_m       = _rd(0,  500.0,   float)
    rho           = _rd(1,  2.65,    float)
    X0_cm         = _rd(2,  26.48,   float)
    mat_id        = _rd(3,  1,       int)
    spec_mode     = _rd(4,  1,       int)
    E_min_GeV     = _rd(5,  0.5,     float)
    E_max_GeV     = _rd(6,  5000.0,  float)
    theta_max_deg = _rd(7,  70.0,    float)
    n_E           = _rd(8,  60,      int)
    n_theta       = _rd(9,  25,      int)
    mode          = _rd(10, 1,       int)
    v_cut         = _rd(11, 0.05,    float)
    b_rad_ov      = _rd(12, 0.0,     float)
    outfile       = _rd(13, "backward_mc_results.dat")

    mat = _MAT_DB.get(mat_id, _MAT_DB[1])
    print(f"  UCMuon Backward MC v1.0 — self-developed CSDA+stochastic integrator",
          flush=True)
    print(f"  {mat['name']}  rho={rho:.3f} g/cm3  depth={depth_m:.1f} m"
          f"  opacity={depth_m*100*rho:.1f} g/cm2", flush=True)
    print(f"  Spectrum {spectrum_name(spec_mode)}  theta_max={theta_max_deg:.0f} deg"
          f"  E={E_min_GeV}–{E_max_GeV} GeV  bins={n_E}x{n_theta}", flush=True)

    res = backward_mc_flux(
        depth_m=depth_m, rho=rho, mat_id=mat_id, spectrum_mode=spec_mode,
        E_min_GeV=E_min_GeV, E_max_GeV=E_max_GeV,
        theta_max_deg=theta_max_deg, n_E=n_E, n_theta=n_theta,
        mode=mode, v_cut=v_cut, b_rad_override=b_rad_ov,
        progress_cb=lambda msg: print(f"  {msg}", flush=True),
    )

    _print_summary(res)
    _write_results(res, outfile)


if __name__ == "__main__":
    main()
