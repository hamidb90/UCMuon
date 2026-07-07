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

Spectrum models (same as Tab 1 GUI):
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

import sys
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

# Pre-compute dE/dx fine grid once at import
_N_FINE   = 1000
_T_FINE   = np.logspace(np.log10(_GROOM_T[0]), np.log10(_GROOM_T[-1]), _N_FINE)
_R_FINE   = np.exp(np.interp(np.log(_T_FINE), np.log(_GROOM_T), np.log(_GROOM_R)))
_DRDT     = np.gradient(_R_FINE, _T_FINE)
_DEDX_STD = 1.0 / np.maximum(_DRDT, 1e-12)        # [MeV cm2/g] Standard Rock


def _dedx(E_MeV, a_scale=1.0):
    E = np.clip(np.asarray(E_MeV, float), _T_FINE[0], _T_FINE[-1])
    return np.interp(E, _T_FINE, _DEDX_STD) * a_scale


def _range(E_MeV, a_scale=1.0):
    E = np.clip(np.asarray(E_MeV, float), _T_FINE[0], _T_FINE[-1])
    return np.interp(E, _T_FINE, _R_FINE) * a_scale


def _backward_energy(E_det_MeV, X_gcm2, a_scale=1.0):
    """
    Backward CSDA mapping: given E_det_MeV at detector, find surface energy E_s
    and Jacobian |dE_s/dE_det| = dE/dx(E_s) / dE/dx(E_det).
    Returns (E_s_MeV, jacobian).  E_s = inf if unreachable (off-table).
    """
    R_det  = _range(E_det_MeV, a_scale)
    R_surf = R_det + X_gcm2
    if np.ndim(R_surf) == 0:
        if R_surf > _R_FINE[-1] * a_scale:
            return np.inf, 1.0
        E_s = float(np.interp(R_surf / a_scale, _R_FINE, _T_FINE))
        J   = float(_dedx(E_s, a_scale) / max(_dedx(E_det_MeV, a_scale), 1e-12))
        return E_s, J
    # Vectorised
    E_s  = np.full_like(np.asarray(E_det_MeV, float), np.inf)
    J    = np.ones_like(E_s)
    ok   = R_surf <= _R_FINE[-1] * a_scale
    E_s[ok]  = np.interp(R_surf[ok] / a_scale, _R_FINE, _T_FINE)
    J[ok]    = _dedx(E_s[ok], a_scale) / np.maximum(_dedx(E_det_MeV[ok], a_scale), 1e-12)
    return E_s, J


# ─────────────────────────────────────────────────────────────────────────────
# Material database
# ─────────────────────────────────────────────────────────────────────────────
_MAT_DB = {
    1: {"name": "Standard Rock", "b_rad": 3.475e-6, "a_scale": 1.000, "X0_cm": 26.48},
    2: {"name": "Water/Ice",     "b_rad": 3.20e-6,  "a_scale": 1.046, "X0_cm": 36.08},
    3: {"name": "Seawater",      "b_rad": 3.22e-6,  "a_scale": 1.028, "X0_cm": 35.75},
    4: {"name": "Iron",          "b_rad": 4.06e-6,  "a_scale": 0.930, "X0_cm": 1.757},
}


# ─────────────────────────────────────────────────────────────────────────────
# Surface flux models
# Approximate absolute normalisation from Gaisser & Stanev (PDG 2022)
# ─────────────────────────────────────────────────────────────────────────────

def _flux_surface(E_GeV, theta_rad, mode):
    """
    Differential surface flux Φ(E, θ) [muons m⁻² s⁻¹ GeV⁻¹ sr⁻¹].
    Normalisation is self-consistent across modes for relative comparisons.
    """
    cos_t = np.cos(theta_rad)
    E     = np.asarray(E_GeV, float)

    if mode == 1:
        # CosmoALEPH (Schmelling 2013):  dN/dp ∝ p^{-3.195}
        # Rough absolute: ~170 m⁻² s⁻¹ sr⁻¹ GeV⁻¹ at 1 GeV vertical
        return 170.0 * E**(-3.195) * cos_t**2

    elif mode == 2:
        # Power-law:  dN/dE ∝ E^{-3.7}
        return 95.0 * E**(-3.7) * cos_t**2

    elif mode == 3:
        # Guan et al. (2015) arXiv:1509.06176
        # Full parametrisation with energy-zenith coupling
        Ec   = 854.0    # critical energy [GeV]
        A    = 0.00253
        a    = 0.2455;  b = 1.288;  c = -0.2455;  d = 0.2949
        flux = (A * E**(-3.7) / (1.0 + E/Ec)**1.55
                * (1.0 + a * E**b * cos_t**(c * E**d + 1.0)))
        return np.maximum(flux, 0.0)

    elif mode == 4:
        # Frosin et al. (2025) J. Phys. G 52, 035002
        Ec   = 902.0
        return np.maximum(0.00245 * E**(-3.72) / (1.0 + E/Ec)**1.52 * cos_t**2, 0.0)

    return np.zeros_like(E)


# ─────────────────────────────────────────────────────────────────────────────
# Stochastic survival probability (Poisson analytical approximation)
# ─────────────────────────────────────────────────────────────────────────────

def _P_surv_stochastic(E_s_MeV, E_det_MeV, X_gcm2, b_total, v_cut, a_scale=1.0):
    """
    Survival probability correction for catastrophic radiative events.

    P_surv = exp(−λ_fatal × X_slant)
    where λ_fatal [events/g/cm²] = b_total × ln(1/v_stop) / ln(1/v_cut)
    and v_stop = min fraction that would stop the muon:
        v_stop = 1 − R_det / R_surf  (fraction of E_s that can be lost)

    Reduces to 1 in CSDA limit (v_stop → 0) and to 0 when every catastrophic
    event above v_cut would stop the muon.
    """
    R_det  = float(_range(E_det_MeV, a_scale))
    R_surf = float(_range(E_s_MeV,   a_scale))
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
    mat      = _MAT_DB.get(mat_id, _MAT_DB[1])
    b_total  = b_rad_override if b_rad_override > 0 else mat["b_rad"]
    a_scale  = mat["a_scale"]

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

    for iE, E_det in enumerate(E_det_MeV):
        w_trap = 0.5 if iE in (0, n_E - 1) else 1.0
        dE_GeV = (E_det / 1000.0) * dlnE * w_trap   # energy bin width [GeV]

        sum_flux = 0.0; sum_surf = 0.0
        sum_Es   = 0.0; sum_Ps  = 0.0

        for ith, (th, dO, cos_t, X_sl) in enumerate(zip(th_mid, dOmega, cos_mid, X_slant)):

            # Backward CSDA: what surface energy E_s is required?
            E_s, Jacob = _backward_energy(E_det, X_sl, a_scale)
            if E_s == np.inf or np.isnan(E_s):
                continue                       # direction unreachable
            E_s_GeV = E_s / 1000.0
            if E_s_GeV < E_min_GeV or E_s_GeV > E_max_GeV:
                continue

            # Surface flux
            phi_s = float(_flux_surface(E_s_GeV, th, spectrum_mode))
            if phi_s <= 0:
                continue

            # Survival probability
            if mode == 1:
                Ps = _P_surv_stochastic(E_s, E_det, X_sl, b_total, v_cut, a_scale)
            else:
                Ps = 1.0

            # Detector flux contribution  [m⁻² s⁻¹ GeV⁻¹ sr⁻¹ → m⁻² s⁻¹ GeV⁻¹]
            phi_det  = phi_s * Jacob * Ps       # [m⁻² s⁻¹ GeV⁻¹ sr⁻¹]
            contrib  = phi_det * dO * cos_t     # [m⁻² s⁻¹ GeV⁻¹]

            sum_flux += contrib
            sum_surf += phi_s * dO * cos_t
            sum_Es   += E_s_GeV * dO
            sum_Ps   += Ps * dO

            rate_total += contrib * dE_GeV

        flux_det[iE]   = sum_flux
        flux_surf[iE]  = sum_surf
        E_surf_out[iE] = sum_Es  / dOmega_total if dOmega_total > 0 else 0.0
        Ps_out[iE]     = sum_Ps  / dOmega_total if dOmega_total > 0 else 0.0

        if progress_cb and iE % max(1, n_E // 10) == 0:
            progress_cb(f"E bin {iE+1}/{n_E}  E_det={E_det/1000:.2f} GeV"
                        f"  rate so far {rate_total:.3e} m-2 s-1")

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
                           v_cut=v_cut),
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
    mat     = _MAT_DB.get(mat_id, _MAT_DB[1])
    b_total = b_rad_override if b_rad_override > 0 else mat["b_rad"]
    a_scale = mat["a_scale"]

    E_det_MeV = np.logspace(np.log10(E_min_GeV * 1000.0),
                             np.log10(E_max_GeV * 1000.0), n_E)
    dlnE = (np.log(E_max_GeV) - np.log(E_min_GeV)) / max(n_E - 1, 1)

    flux = 0.0
    for iE, E_det in enumerate(E_det_MeV):
        w_trap = 0.5 if iE in (0, n_E - 1) else 1.0
        dE_GeV = (E_det / 1000.0) * dlnE * w_trap
        if X_slant_gcm2 > 0.0:
            E_s, Jacob = _backward_energy(E_det, X_slant_gcm2, a_scale)
            if not np.isfinite(E_s):
                continue
            Ps = (_P_surv_stochastic(E_s, E_det, X_slant_gcm2,
                                     b_total, v_cut, a_scale)
                  if mode == 1 else 1.0)
        else:
            E_s, Jacob, Ps = E_det, 1.0, 1.0
        phi_s = float(_flux_surface(E_s / 1000.0, theta_rad, spectrum_mode))
        if phi_s > 0.0:
            flux += phi_s * Jacob * Ps * dE_GeV
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
        fh.write(f"# Mode: {info['mode']}  spectrum: {info['spectrum_mode']}\n")
        fh.write(f"# Total expected rate: {res['rate_m2_s']:.5e} m-2 s-1\n")
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
    print(f"  Spectrum mode {spec_mode}  theta_max={theta_max_deg:.0f} deg"
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
