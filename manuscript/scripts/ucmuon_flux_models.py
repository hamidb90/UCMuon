#!/usr/bin/env python3
"""
ucmuon_flux_models.py — Python transcription of the surface flux
parametrisations in src/generator/ucmuon_source_module.f90, used by the
manuscript figure scripts (fig02, fig03).

All constants are copied verbatim from the Fortran module so that the
curves shown in the paper are exactly the spectra the generator samples.

Differential intensities are returned in cm^-2 s^-1 sr^-1 GeV^-1:
  * cosmoaleph(p)      — dN/dp  (per GeV/c)
  * powerlaw(p)        — dN/dp, arbitrary normalisation (see docstring)
  * guan_flux(E, cos)  — dN/dE  (per GeV), Gaisser-type formulas
  * reyna_flux(p, cos) — dN/dp  (per GeV/c)
"""
from __future__ import annotations

import numpy as np

# ── Constants (ucmuon_source_module.f90) ──────────────────────────────────────
MUON_MASS = 0.10566          # GeV/c^2

A_COSMO = 3.8467
B_COSMO = -3.1952

GUAN_P1 = 0.102573
GUAN_P2 = -0.068287
GUAN_P3 = 0.958633
GUAN_P4 = 0.0407253
GUAN_P5 = 0.817285
GUAN_DENOM = 0.99144315
GUAN_EPI = 115.0
GUAN_EK = 850.0
GUAN_KF = 0.054
GUAN_PRE = 0.14
GUAN_IDX = -2.7
GUAN_A, GUAN_B = 3.64, 1.29
FROSIN_A, FROSIN_B = 3.512, 1.388
BUGAEV_A, BUGAEV_B = 0.0, 1.0

REYNA_C0 = 0.00253
REYNA_C1 = 0.2455
REYNA_C2 = 1.288
REYNA_C3 = -0.2555
REYNA_C4 = 0.0209


# ── Kinematics helpers ─────────────────────────────────────────────────────────
def p_of_T(T):
    """Momentum [GeV/c] from kinetic energy T [GeV]."""
    E = np.asarray(T) + MUON_MASS
    return np.sqrt(E**2 - MUON_MASS**2)


def jac_dp_dT(T):
    """dp/dT = E/p — converts dN/dp to dN/dT (= dN/dE)."""
    E = np.asarray(T) + MUON_MASS
    return E / p_of_T(T)


# ── Models ─────────────────────────────────────────────────────────────────────
def cosmoaleph_dndp(p):
    """Mode 1 — CosmoALEPH: dN/dp = 10^3.8467 * p^-3.1952 (absolute units
    follow the ALEPH fit used by the generator)."""
    return 10.0**A_COSMO * np.asarray(p, float) ** B_COSMO


def powerlaw_dndp(p, norm_at=50.0):
    """Mode 2 — Kudryavtsev/MUSIC power-law dN/dE ∝ E^-3.7.
    The generator treats this spectrum as shape-only (arbitrary absolute
    normalisation); here it is pinned to the CosmoALEPH intensity at
    p = `norm_at` GeV/c so both curves can share one absolute axis."""
    p = np.asarray(p, float)
    shape = p**-3.7
    return shape * cosmoaleph_dndp(norm_at) / (norm_at**-3.7)


def guan_cos_star(cos_th):
    """Atmospheric curvature correction cos θ* (Chirkin / Guan Eq. 4)."""
    c = np.asarray(cos_th, float)
    numer = c**2 + GUAN_P1**2 + GUAN_P2 * c**GUAN_P3 + GUAN_P4 * c**GUAN_P5
    cs = np.sqrt(np.maximum(0.0, numer)) / GUAN_DENOM
    return np.clip(cs, 0.0, 1.0)


def guan_flux(E_GeV, cos_th, a_par=GUAN_A, b_par=GUAN_B):
    """Modified Gaisser formula (Guan 2015 / Frosin 2025 / Bugaev-Gaisser).
    Returns dN/dE [cm^-2 s^-1 sr^-1 GeV^-1] at total energy E [GeV].
    Mode 4: (a,b)=(3.64,1.29); Mode 5: (3.512,1.388); Mode 6: (0,1)."""
    E = np.asarray(E_GeV, float)
    cs = guan_cos_star(cos_th)
    E_eff = E * (1.0 + a_par / (E * cs**b_par))
    pion_t = 1.0 / (1.0 + 1.1 * E * cs / GUAN_EPI)
    kaon_t = GUAN_KF / (1.0 + 1.1 * E * cs / GUAN_EK)
    return GUAN_PRE * E_eff**GUAN_IDX * (pion_t + kaon_t)


def reyna_flux(p_GeV, cos_th):
    """Mode 7 — Reyna–Bugaev (2006) dΦ/dp [cm^-2 s^-1 sr^-1 (GeV/c)^-1].
    I_V(p) = C0 * p^-(C1 + C2 z + C3 z^2 + C4 z^3), z = log10(p)
    (hep-ph/0604145 Eq. 6-7); integrates to 7.0e-3 cm^-2 s^-1 sr^-1
    above 1 GeV at cos(theta) = 1, matching the PDG reference value."""
    p = np.asarray(p_GeV, float)
    cs = guan_cos_star(cos_th)
    p_eff = p * cs
    out = np.zeros_like(p)
    ok = p_eff > 0
    lp = np.log10(p_eff[ok])
    n_exp = REYNA_C1 + REYNA_C2 * lp + REYNA_C3 * lp**2 + REYNA_C4 * lp**3
    out[ok] = REYNA_C0 * p_eff[ok] ** (-n_exp)
    return np.maximum(out, 0.0)


# ── Generator zenith-angle PDFs (sample_cos2 / sample_guan_angle) ─────────────
def pdf_cos2_theta(theta, theta_max):
    """PDF in θ sampled by the generator's cos²θ mode: ∝ cos²θ sinθ."""
    th = np.asarray(theta, float)
    pdf = np.where(th <= theta_max, np.cos(th) ** 2 * np.sin(th), 0.0)
    return pdf / np.trapz(pdf, th)


def pdf_cosn_theta(theta, theta_max, n):
    """PDF in θ sampled by the generator's cosⁿθ mode: ∝ cosⁿθ sinθ.
    n=2 reproduces pdf_cos2_theta; n=3 is the Reyna (Mode 7) recommendation."""
    th = np.asarray(theta, float)
    pdf = np.where(th <= theta_max, np.cos(th) ** n * np.sin(th), 0.0)
    return pdf / np.trapz(pdf, th)


def pdf_guan_theta(theta, theta_max, E_GeV, a_par=GUAN_A, b_par=GUAN_B):
    """PDF in θ sampled by sample_guan_angle: ∝ Φ(E, cosθ)·cosθ in cosθ,
    i.e. ∝ Φ(E, cosθ)·cosθ·sinθ in θ (projected-area weighting)."""
    th = np.asarray(theta, float)
    pdf = np.where(
        th <= theta_max,
        guan_flux(E_GeV, np.cos(th), a_par, b_par) * np.cos(th) * np.sin(th),
        0.0,
    )
    return pdf / np.trapz(pdf, th)
