#!/usr/bin/env python3
"""
ucmuon_stochastic_driver.py  —  UCLouvain Muography Group
UCMuon-MC: the native forward stochastic muon transport engine (pure Python/NumPy).

Physics model (v2 — table-anchored stochastic decomposition)
─────────────────────────────────────────────────────────────
Total energy loss per step dx [g/cm²]:

    dE = a_res(E) × dx  +  ΔE_δ  +  v_cut × L_rad(E) × dx  +  ΔE_hard

  a_res    : restricted ionisation = a_ion(E) − k_δ·ln(T_max/T_cut)
  ΔE_δ     : δ-rays (knock-on electrons), Poisson-sampled above T_cut=10 MeV
             from the Rutherford 1/T² close-collision spectrum
  L_rad(E) : L_brems(E) + L_pair(E) + L_photonuc(E)  — per-process mean
             radiative losses from the PDG 2024 evaluated table for
             Standard Rock (embedded; pdg.lbl.gov muE table)
  ΔE_hard  : catastrophic radiative events, Poisson-sampled above v_cut
             per process: brems φ(v)∝(1−v)/v, pair φ(v)∝1/v³,
             photonuclear φ(v)∝1/v, each with rate
             λ_i = (1−v_cut)·L_i(E)/(E_tot·<v>_i)

  Conservation: every stochastic term subtracts its own mean from the
  continuous part, so <dE/dx> equals the evaluated PDG table EXACTLY at
  all energies, by construction.  Single-process legacy spectra
  (1/v, Bethe-Heitler) remain selectable for reproducibility.

Additional physics:
  · Multiple scattering  : Highland (1979) projected angles per step
  · Muon decay           : Poisson probability per step (p·cτ, p = √(E²−m²))
  · Pre-filter           : muons whose *deterministic-loss* range (a_ion +
                           v_cut·L_rad only — a strict upper bound on
                           penetration) is below the slant path are marked
                           dead instantly.  Muons between the mean-loss CSDA
                           range and this bound are transported stochastically
                           (they survive if they avoid hard radiative events).
  · Adaptive stepping    : per-muon step count targets dx ≈ 5 g/cm²
                           (near-horizontal muons no longer get oversized steps)

Position update (lateral displacement bug fix)
──────────────────────────────────────────────
Previously, lateral exit position was computed as xs + cx_final × slant_cm,
i.e. the FINAL (post-all-scattering) direction times the FULL path length.
This over-estimated lateral displacement by ~2× because the muon only
occupies the final direction for the last step, not the entire path.

The correct formula integrates position step-by-step:
  x(t+dt) = x(t) + cx(t) × (dx_gcm2 / rho)
where cx(t) is the current direction before the MCS update for that step.
This gives the physical random-walk lateral displacement.

Progress output (parsed by GUI parse_progress):
    'Transported: N  Survived: M  Total: T'

Output: 18-column underground file (same format as MUSIC)
    EventID xs ys zs Es[GeV] θs φs charge alive x y z E[GeV] cx cy cz θ φ

Stdin input:
    1  infile        e.g. muons_selected.dat
    2  outfile       e.g. muons_stochastic.dat
    3  depth_m       overburden depth [m]
    4  rho           density [g/cm³]
    5  X0_cm         radiation length [cm]
    6  mat_id        1=rock  2=water  3=seawater  4=iron  5=custom  6=ice
                     (1, 2, 4, 6 use native embedded PDG per-process tables;
                      3 rescales the water table; 5 rescales the rock table)
    7  transport_all 0=hit_flag only  1=all
    8  ncols         informational (reader auto-detects)
    9  n_steps       transport steps (0 = auto)
   10  v_cut         catastrophic threshold (default 0.05)
   11  ms_enable     1=Highland MS on  0=off
   12  Z_eff         [custom only]
   13  A_eff         [custom only]
   14  I_eV          [custom only]
   15  b_rad         [custom only, cm²/g]
   16  range_table   0=groom2001  1=pdg2024 (legacy; v2 loss model is PDG-anchored)
   17  hard_spectrum 0=groom (1/v)  1=bh ((1-v)/v)  2=per-process (default)
   18  seed          RNG seed (int, default 42; reproducible runs by default)
   19  n_workers     parallel worker processes (default 1 = serial;
                     0 = auto: one worker per ~20k muons, up to all cores).
                     Results are reproducible for a given (seed, n_workers).
   20  delta_rays    1=explicit δ-ray straggling (default)  0=off

Author: Hamid Basiri <hamid.basiri@uclouvain.be>
"""

import sys
import time
import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Physical constants
# ─────────────────────────────────────────────────────────────────────────────
M_MU      = 105.6584        # muon rest mass  [MeV]
M_MU_GEV  = 0.1056584       # muon rest mass  [GeV]
C_CM_S    = 2.99792458e10   # speed of light  [cm/s]
TAU_MU    = 2.1969811e-6    # muon lifetime   [s]
PCTAU_CM  = C_CM_S * TAU_MU / M_MU_GEV  # cτ/m → dec_len=(p/m)cτ=βγcτ [cm/(GeV/c)]

REPORT_EVERY = 500          # progress print every N transported muons

# ─────────────────────────────────────────────────────────────────────────────
# Groom (2001) Standard Rock CSDA range table  (33 entries, Table IV-6)
# Standard Rock: Z=11, A=22, rho=2.65 g/cm3, I=136.4 eV
# Same values as _GROOM_T_MEV / _GROOM_R_GCM2 in cosmoaleph_gui_v17_omp.py
# ─────────────────────────────────────────────────────────────────────────────
_GROOM_T = np.array([
    1.0e+01, 1.4e+01, 2.0e+01, 3.0e+01, 4.0e+01,
    8.0e+01, 1.0e+02, 1.4e+02, 2.0e+02, 3.0e+02,
    4.0e+02, 8.0e+02, 1.0e+03, 1.4e+03, 2.0e+03,
    3.0e+03, 4.0e+03, 8.0e+03, 1.0e+04, 1.4e+04,
    2.0e+04, 3.0e+04, 4.0e+04, 8.0e+04, 1.0e+05,
    1.4e+05, 2.0e+05, 3.0e+05, 4.0e+05, 8.0e+05,
    1.0e+06, 1.4e+06, 2.0e+06,
])                          # kinetic energy T [MeV]

_GROOM_R = np.array([       # CSDA range R [g/cm²]
    8.516e-01, 1.542e+00, 2.866e+00, 5.698e+00, 9.145e+00,
    2.676e+01, 3.696e+01, 5.879e+01, 9.332e+01, 1.524e+02,
    2.115e+02, 4.418e+02, 5.534e+02, 7.712e+02, 1.088e+03,
    1.599e+03, 2.095e+03, 3.998e+03, 4.920e+03, 6.724e+03,
    9.360e+03, 1.362e+04, 1.776e+04, 3.343e+04, 4.084e+04,
    5.495e+04, 7.459e+04, 1.040e+05, 1.302e+05, 2.129e+05,
    2.453e+05, 2.990e+05, 3.616e+05,
])

# ─────────────────────────────────────────────────────────────────────────────
# Pre-compute dE/dx on a 1000-point log-grid at import time (one-time cost ~1 ms)
# dE/dx(T) [MeV cm2/g] = 1 / (dR/dT) by numerical differentiation
# ─────────────────────────────────────────────────────────────────────────────
_N_FINE   = 1000
_T_FINE   = np.logspace(np.log10(_GROOM_T[0]), np.log10(_GROOM_T[-1]), _N_FINE)
_R_FINE   = np.exp(np.interp(np.log(_T_FINE), np.log(_GROOM_T), np.log(_GROOM_R)))
_DRDT     = np.gradient(_R_FINE, _T_FINE)          # [g/cm² / MeV]
_DEDX_STD = 1.0 / np.maximum(_DRDT, 1e-12)         # [MeV cm2/g]  Standard Rock

# ─────────────────────────────────────────────────────────────────────────────
# PDG 2024 Standard Rock table  (56 entries, T in MeV)
# Source: PDG 2024 stopping power table, Standard Rock Z=11 A=22 I=136.4 eV
# Direct dE/dx column avoids range-differentiation noise.
# Values are essentially identical to Groom 2001 at overlapping energies.
# ─────────────────────────────────────────────────────────────────────────────
_PDG24_T = np.array([
    1.0e+01, 1.4e+01, 2.0e+01, 3.0e+01, 4.0e+01,
    8.0e+01, 1.0e+02, 1.4e+02, 2.0e+02, 3.0e+02,
    4.0e+02, 8.0e+02, 1.0e+03, 1.4e+03, 2.0e+03,
    3.0e+03, 4.0e+03, 8.0e+03, 9.0e+03, 1.0e+04,
    1.2e+04, 1.4e+04, 1.7e+04, 2.0e+04, 2.5e+04,
    3.0e+04, 3.5e+04, 4.0e+04, 4.5e+04, 5.0e+04,
    5.5e+04, 6.0e+04, 7.0e+04, 8.0e+04, 9.0e+04,
    1.0e+05, 1.2e+05, 1.4e+05, 1.7e+05, 2.0e+05,
    2.5e+05, 3.0e+05, 3.5e+05, 4.0e+05, 4.5e+05,
    5.0e+05, 5.5e+05, 6.0e+05, 7.0e+05, 8.0e+05,
    9.0e+05, 1.0e+06, 1.2e+06, 1.4e+06, 1.7e+06,
    2.0e+06,
])  # kinetic energy T [MeV]

_PDG24_DEDX = np.array([
    6.619, 5.180, 4.057, 3.157, 2.702,
    2.029, 1.904, 1.779, 1.710, 1.688,
    1.698, 1.775, 1.808, 1.862, 1.922,
    1.990, 2.038, 2.152, 2.171, 2.188,
    2.218, 2.244, 2.277, 2.305, 2.347,
    2.383, 2.416, 2.447, 2.476, 2.503,
    2.530, 2.556, 2.605, 2.653, 2.700,
    2.747, 2.836, 2.925, 3.056, 3.187,
    3.399, 3.610, 3.823, 4.035, 4.249,
    4.463, 4.675, 4.888, 5.315, 5.745,
    6.177, 6.612, 7.471, 8.336, 9.642,
    1.096e+01,
])  # total dE/dx [MeV cm²/g]

_PDG24_R = np.array([
    8.400e-01, 1.530e+00, 2.854e+00, 5.687e+00, 9.133e+00,
    2.675e+01, 3.695e+01, 5.878e+01, 9.331e+01, 1.523e+02,
    2.114e+02, 4.418e+02, 5.534e+02, 7.712e+02, 1.088e+03,
    1.599e+03, 2.095e+03, 3.998e+03, 4.461e+03, 4.920e+03,
    5.828e+03, 6.724e+03, 8.051e+03, 9.360e+03, 1.151e+04,
    1.362e+04, 1.571e+04, 1.776e+04, 1.979e+04, 2.180e+04,
    2.379e+04, 2.576e+04, 2.963e+04, 3.343e+04, 3.717e+04,
    4.084e+04, 4.801e+04, 5.495e+04, 6.498e+04, 7.459e+04,
    8.978e+04, 1.041e+05, 1.175e+05, 1.302e+05, 1.423e+05,
    1.538e+05, 1.647e+05, 1.752e+05, 1.948e+05, 2.129e+05,
    2.297e+05, 2.453e+05, 2.738e+05, 2.991e+05, 3.325e+05,
    3.617e+05,
])  # CSDA range R [g/cm²]

_PDG24_T_FINE    = np.logspace(np.log10(_PDG24_T[0]), np.log10(_PDG24_T[-1]), _N_FINE)
_PDG24_DEDX_FINE = np.exp(np.interp(
    np.log(_PDG24_T_FINE), np.log(_PDG24_T), np.log(_PDG24_DEDX)))
_PDG24_R_FINE    = np.exp(np.interp(
    np.log(_PDG24_T_FINE), np.log(_PDG24_T), np.log(_PDG24_R)))

# ─────────────────────────────────────────────────────────────────────────────
# Per-process mean muon energy losses [MeV cm²/g] — embedded evaluated tables
# Source: Groom, Mokhov & Striganov, At. Data Nucl. Data Tables 78 (2001) 183,
#         in the PDG 2024 electronic revision (post-Born pair production):
#         pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_<material>.txt
# Materials embedded natively (97 entries each, common T grid 10 MeV–10 TeV):
#   rock  = standard rock          water = water (liquid)
#   ice   = water (ice)            iron  = iron (Fe)
# Seawater uses the water table with sub-percent composition rescaling;
# custom materials rescale the rock table by (Z/A)/0.5 and b_rad/b_rock.
# Verified: ION + BREMS + PAIR + NUCL = dE/dx column to < 0.07% in all files.
# ─────────────────────────────────────────────────────────────────────────────
_PROC_T = np.array([
    1.000e+01, 1.200e+01, 1.400e+01, 1.700e+01, 2.000e+01,
    2.500e+01, 3.000e+01, 3.500e+01, 4.000e+01, 4.500e+01,
    5.000e+01, 5.500e+01, 6.000e+01, 7.000e+01, 8.000e+01,
    9.000e+01, 1.000e+02, 1.200e+02, 1.400e+02, 1.700e+02,
    2.000e+02, 2.500e+02, 3.000e+02, 3.500e+02, 4.000e+02,
    4.500e+02, 5.000e+02, 5.500e+02, 6.000e+02, 7.000e+02,
    8.000e+02, 9.000e+02, 1.000e+03, 1.200e+03, 1.400e+03,
    1.700e+03, 2.000e+03, 2.500e+03, 3.000e+03, 3.500e+03,
    4.000e+03, 4.500e+03, 5.000e+03, 5.500e+03, 6.000e+03,
    7.000e+03, 8.000e+03, 9.000e+03, 1.000e+04, 1.200e+04,
    1.400e+04, 1.700e+04, 2.000e+04, 2.500e+04, 3.000e+04,
    3.500e+04, 4.000e+04, 4.500e+04, 5.000e+04, 5.500e+04,
    6.000e+04, 7.000e+04, 8.000e+04, 9.000e+04, 1.000e+05,
    1.200e+05, 1.400e+05, 1.700e+05, 2.000e+05, 2.500e+05,
    3.000e+05, 3.500e+05, 4.000e+05, 4.500e+05, 5.000e+05,
    5.500e+05, 6.000e+05, 7.000e+05, 8.000e+05, 9.000e+05,
    1.000e+06, 1.200e+06, 1.400e+06, 1.700e+06, 2.000e+06,
    2.500e+06, 3.000e+06, 3.500e+06, 4.000e+06, 4.500e+06,
    5.000e+06, 5.500e+06, 6.000e+06, 7.000e+06, 8.000e+06,
    9.000e+06, 1.000e+07,
])
_PROC_TABLES = {
    "rock": {
        "ION": np.array([
    6.619e+00, 5.787e+00, 5.180e+00, 4.524e+00, 4.057e+00,
    3.520e+00, 3.157e+00, 2.897e+00, 2.701e+00, 2.550e+00,
    2.430e+00, 2.331e+00, 2.249e+00, 2.121e+00, 2.028e+00,
    1.958e+00, 1.904e+00, 1.828e+00, 1.779e+00, 1.734e+00,
    1.710e+00, 1.691e+00, 1.688e+00, 1.691e+00, 1.698e+00,
    1.707e+00, 1.716e+00, 1.726e+00, 1.736e+00, 1.756e+00,
    1.774e+00, 1.792e+00, 1.808e+00, 1.836e+00, 1.861e+00,
    1.893e+00, 1.920e+00, 1.957e+00, 1.986e+00, 2.011e+00,
    2.033e+00, 2.051e+00, 2.067e+00, 2.082e+00, 2.095e+00,
    2.118e+00, 2.138e+00, 2.155e+00, 2.170e+00, 2.195e+00,
    2.215e+00, 2.241e+00, 2.262e+00, 2.289e+00, 2.311e+00,
    2.329e+00, 2.344e+00, 2.358e+00, 2.369e+00, 2.380e+00,
    2.389e+00, 2.406e+00, 2.420e+00, 2.433e+00, 2.444e+00,
    2.463e+00, 2.479e+00, 2.499e+00, 2.515e+00, 2.538e+00,
    2.557e+00, 2.573e+00, 2.586e+00, 2.598e+00, 2.609e+00,
    2.619e+00, 2.628e+00, 2.644e+00, 2.658e+00, 2.670e+00,
    2.681e+00, 2.700e+00, 2.717e+00, 2.737e+00, 2.755e+00,
    2.779e+00, 2.798e+00, 2.815e+00, 2.830e+00, 2.843e+00,
    2.855e+00, 2.865e+00, 2.875e+00, 2.892e+00, 2.908e+00,
    2.921e+00, 2.933e+00,
        ]),
        "BREMS": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    5.446e-07, 1.383e-06, 2.246e-06, 4.044e-06, 5.935e-06,
    7.912e-06, 9.973e-06, 1.433e-05, 1.896e-05, 2.641e-05,
    3.438e-05, 4.871e-05, 6.418e-05, 8.064e-05, 9.800e-05,
    1.162e-04, 1.350e-04, 1.546e-04, 1.747e-04, 2.168e-04,
    2.608e-04, 3.066e-04, 3.540e-04, 4.531e-04, 5.572e-04,
    7.212e-04, 8.933e-04, 1.195e-03, 1.512e-03, 1.843e-03,
    2.184e-03, 2.536e-03, 2.898e-03, 3.273e-03, 3.656e-03,
    4.442e-03, 5.253e-03, 6.086e-03, 6.939e-03, 8.708e-03,
    1.054e-02, 1.338e-02, 1.631e-02, 2.139e-02, 2.666e-02,
    3.208e-02, 3.763e-02, 4.330e-02, 4.906e-02, 5.489e-02,
    6.080e-02, 7.283e-02, 8.512e-02, 9.764e-02, 1.103e-01,
    1.361e-01, 1.624e-01, 2.027e-01, 2.440e-01, 3.132e-01,
    3.839e-01, 4.559e-01, 5.290e-01, 6.030e-01, 6.778e-01,
    7.520e-01, 8.268e-01, 9.779e-01, 1.131e+00, 1.285e+00,
    1.441e+00, 1.751e+00, 2.065e+00, 2.541e+00, 3.023e+00,
    3.820e+00, 4.625e+00, 5.437e+00, 6.254e+00, 7.075e+00,
    7.901e+00, 8.718e+00, 9.538e+00, 1.118e+01, 1.284e+01,
    1.450e+01, 1.616e+01,
        ]),
        "PAIR": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 1.438e-05, 8.199e-05, 1.588e-04,
    2.887e-04, 4.336e-04, 7.028e-04, 1.001e-03, 1.323e-03,
    1.666e-03, 2.028e-03, 2.408e-03, 2.814e-03, 3.234e-03,
    4.112e-03, 5.037e-03, 6.001e-03, 7.004e-03, 9.152e-03,
    1.142e-02, 1.499e-02, 1.876e-02, 2.550e-02, 3.262e-02,
    4.005e-02, 4.775e-02, 5.568e-02, 6.382e-02, 7.200e-02,
    8.033e-02, 9.742e-02, 1.150e-01, 1.330e-01, 1.514e-01,
    1.886e-01, 2.268e-01, 2.857e-01, 3.463e-01, 4.461e-01,
    5.483e-01, 6.524e-01, 7.583e-01, 8.656e-01, 9.742e-01,
    1.082e+00, 1.190e+00, 1.408e+00, 1.630e+00, 1.853e+00,
    2.079e+00, 2.523e+00, 2.972e+00, 3.652e+00, 4.338e+00,
    5.470e+00, 6.611e+00, 7.759e+00, 8.913e+00, 1.007e+01,
    1.124e+01, 1.239e+01, 1.354e+01, 1.586e+01, 1.818e+01,
    2.051e+01, 2.285e+01,
        ]),
        "NUCL": np.array([
    5.173e-05, 5.263e-05, 5.352e-05, 5.487e-05, 5.621e-05,
    5.844e-05, 6.068e-05, 6.292e-05, 6.515e-05, 6.739e-05,
    6.963e-05, 7.186e-05, 7.410e-05, 7.857e-05, 8.304e-05,
    8.752e-05, 9.199e-05, 1.009e-04, 1.099e-04, 1.233e-04,
    1.367e-04, 1.591e-04, 1.815e-04, 2.038e-04, 2.262e-04,
    2.485e-04, 2.709e-04, 2.933e-04, 3.156e-04, 3.604e-04,
    4.051e-04, 4.498e-04, 4.946e-04, 5.840e-04, 6.735e-04,
    8.077e-04, 9.451e-04, 1.186e-03, 1.431e-03, 1.677e-03,
    1.926e-03, 2.177e-03, 2.424e-03, 2.652e-03, 2.879e-03,
    3.331e-03, 3.780e-03, 4.227e-03, 4.671e-03, 5.532e-03,
    6.384e-03, 7.646e-03, 8.894e-03, 1.097e-02, 1.303e-02,
    1.506e-02, 1.708e-02, 1.908e-02, 2.107e-02, 2.309e-02,
    2.511e-02, 2.912e-02, 3.310e-02, 3.707e-02, 4.102e-02,
    4.905e-02, 5.705e-02, 6.902e-02, 8.094e-02, 1.011e-01,
    1.213e-01, 1.415e-01, 1.616e-01, 1.818e-01, 2.019e-01,
    2.226e-01, 2.434e-01, 2.849e-01, 3.267e-01, 3.685e-01,
    4.104e-01, 4.959e-01, 5.819e-01, 7.118e-01, 8.424e-01,
    1.065e+00, 1.290e+00, 1.517e+00, 1.745e+00, 1.975e+00,
    2.206e+00, 2.442e+00, 2.679e+00, 3.157e+00, 3.639e+00,
    4.124e+00, 4.613e+00,
        ]),
        "R10": 0.8400,   # CSDA range at 10 MeV [g/cm2]
    },
    "water": {
        "ION": np.array([
    7.902e+00, 6.897e+00, 6.166e+00, 5.378e+00, 4.817e+00,
    4.172e+00, 3.738e+00, 3.426e+00, 3.192e+00, 3.011e+00,
    2.867e+00, 2.750e+00, 2.654e+00, 2.506e+00, 2.398e+00,
    2.317e+00, 2.256e+00, 2.165e+00, 2.103e+00, 2.047e+00,
    2.015e+00, 1.989e+00, 1.981e+00, 1.982e+00, 1.987e+00,
    1.995e+00, 2.005e+00, 2.015e+00, 2.025e+00, 2.045e+00,
    2.064e+00, 2.081e+00, 2.098e+00, 2.128e+00, 2.154e+00,
    2.188e+00, 2.216e+00, 2.255e+00, 2.287e+00, 2.313e+00,
    2.336e+00, 2.356e+00, 2.373e+00, 2.389e+00, 2.403e+00,
    2.427e+00, 2.448e+00, 2.466e+00, 2.482e+00, 2.509e+00,
    2.531e+00, 2.559e+00, 2.581e+00, 2.611e+00, 2.635e+00,
    2.655e+00, 2.671e+00, 2.686e+00, 2.699e+00, 2.711e+00,
    2.721e+00, 2.739e+00, 2.755e+00, 2.769e+00, 2.781e+00,
    2.802e+00, 2.820e+00, 2.842e+00, 2.861e+00, 2.886e+00,
    2.907e+00, 2.925e+00, 2.940e+00, 2.953e+00, 2.965e+00,
    2.976e+00, 2.986e+00, 3.004e+00, 3.019e+00, 3.033e+00,
    3.045e+00, 3.066e+00, 3.084e+00, 3.107e+00, 3.127e+00,
    3.153e+00, 3.175e+00, 3.194e+00, 3.210e+00, 3.225e+00,
    3.238e+00, 3.249e+00, 3.260e+00, 3.280e+00, 3.296e+00,
    3.311e+00, 3.325e+00,
        ]),
        "BREMS": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    2.570e-07, 8.393e-07, 1.439e-06, 2.690e-06, 4.004e-06,
    5.380e-06, 6.813e-06, 9.843e-06, 1.307e-05, 1.826e-05,
    2.381e-05, 3.380e-05, 4.458e-05, 5.606e-05, 6.816e-05,
    8.082e-05, 9.399e-05, 1.076e-04, 1.217e-04, 1.510e-04,
    1.817e-04, 2.137e-04, 2.468e-04, 3.160e-04, 3.886e-04,
    5.031e-04, 6.233e-04, 8.340e-04, 1.056e-03, 1.286e-03,
    1.525e-03, 1.771e-03, 2.023e-03, 2.286e-03, 2.555e-03,
    3.106e-03, 3.675e-03, 4.259e-03, 4.857e-03, 6.103e-03,
    7.392e-03, 9.393e-03, 1.146e-02, 1.506e-02, 1.879e-02,
    2.263e-02, 2.657e-02, 3.059e-02, 3.469e-02, 3.884e-02,
    4.305e-02, 5.163e-02, 6.039e-02, 6.932e-02, 7.840e-02,
    9.681e-02, 1.156e-01, 1.446e-01, 1.741e-01, 2.239e-01,
    2.748e-01, 3.266e-01, 3.792e-01, 4.325e-01, 4.865e-01,
    5.401e-01, 5.941e-01, 7.034e-01, 8.139e-01, 9.257e-01,
    1.039e+00, 1.263e+00, 1.491e+00, 1.836e+00, 2.186e+00,
    2.766e+00, 3.352e+00, 3.942e+00, 4.537e+00, 5.136e+00,
    5.738e+00, 6.334e+00, 6.932e+00, 8.132e+00, 9.339e+00,
    1.055e+01, 1.177e+01,
        ]),
        "PAIR": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 6.995e-06, 5.275e-05, 1.048e-04,
    1.929e-04, 2.913e-04, 4.743e-04, 6.770e-04, 8.963e-04,
    1.130e-03, 1.376e-03, 1.636e-03, 1.918e-03, 2.210e-03,
    2.822e-03, 3.466e-03, 4.140e-03, 4.841e-03, 6.350e-03,
    7.943e-03, 1.046e-02, 1.312e-02, 1.788e-02, 2.291e-02,
    2.817e-02, 3.361e-02, 3.923e-02, 4.499e-02, 5.079e-02,
    5.671e-02, 6.885e-02, 8.135e-02, 9.416e-02, 1.073e-01,
    1.337e-01, 1.610e-01, 2.031e-01, 2.464e-01, 3.179e-01,
    3.913e-01, 4.662e-01, 5.423e-01, 6.196e-01, 6.978e-01,
    7.752e-01, 8.532e-01, 1.011e+00, 1.171e+00, 1.333e+00,
    1.496e+00, 1.817e+00, 2.142e+00, 2.634e+00, 3.132e+00,
    3.953e+00, 4.781e+00, 5.614e+00, 6.453e+00, 7.295e+00,
    8.142e+00, 8.978e+00, 9.816e+00, 1.150e+01, 1.319e+01,
    1.488e+01, 1.658e+01,
        ]),
        "NUCL": np.array([
    5.432e-05, 5.526e-05, 5.620e-05, 5.760e-05, 5.901e-05,
    6.136e-05, 6.371e-05, 6.606e-05, 6.841e-05, 7.075e-05,
    7.310e-05, 7.545e-05, 7.780e-05, 8.250e-05, 8.719e-05,
    9.189e-05, 9.658e-05, 1.060e-04, 1.154e-04, 1.295e-04,
    1.435e-04, 1.670e-04, 1.905e-04, 2.140e-04, 2.375e-04,
    2.610e-04, 2.844e-04, 3.079e-04, 3.314e-04, 3.784e-04,
    4.253e-04, 4.723e-04, 5.193e-04, 6.132e-04, 7.071e-04,
    8.480e-04, 9.922e-04, 1.245e-03, 1.500e-03, 1.757e-03,
    2.017e-03, 2.279e-03, 2.537e-03, 2.774e-03, 3.010e-03,
    3.480e-03, 3.947e-03, 4.411e-03, 4.871e-03, 5.765e-03,
    6.649e-03, 7.957e-03, 9.249e-03, 1.140e-02, 1.352e-02,
    1.563e-02, 1.771e-02, 1.978e-02, 2.183e-02, 2.392e-02,
    2.600e-02, 3.014e-02, 3.426e-02, 3.835e-02, 4.243e-02,
    5.071e-02, 5.897e-02, 7.130e-02, 8.359e-02, 1.044e-01,
    1.252e-01, 1.461e-01, 1.669e-01, 1.877e-01, 2.085e-01,
    2.298e-01, 2.512e-01, 2.941e-01, 3.371e-01, 3.803e-01,
    4.235e-01, 5.118e-01, 6.006e-01, 7.346e-01, 8.696e-01,
    1.100e+00, 1.332e+00, 1.567e+00, 1.803e+00, 2.041e+00,
    2.280e+00, 2.524e+00, 2.769e+00, 3.264e+00, 3.764e+00,
    4.267e+00, 4.773e+00,
        ]),
        "R10": 0.6998,   # CSDA range at 10 MeV [g/cm2]
    },
    "ice": {
        "ION": np.array([
    7.902e+00, 6.897e+00, 6.166e+00, 5.378e+00, 4.817e+00,
    4.172e+00, 3.738e+00, 3.426e+00, 3.192e+00, 3.011e+00,
    2.867e+00, 2.750e+00, 2.654e+00, 2.506e+00, 2.398e+00,
    2.317e+00, 2.256e+00, 2.167e+00, 2.106e+00, 2.050e+00,
    2.018e+00, 1.992e+00, 1.985e+00, 1.986e+00, 1.992e+00,
    2.000e+00, 2.009e+00, 2.019e+00, 2.029e+00, 2.049e+00,
    2.069e+00, 2.087e+00, 2.103e+00, 2.134e+00, 2.160e+00,
    2.194e+00, 2.222e+00, 2.261e+00, 2.293e+00, 2.320e+00,
    2.342e+00, 2.362e+00, 2.380e+00, 2.395e+00, 2.409e+00,
    2.434e+00, 2.455e+00, 2.473e+00, 2.489e+00, 2.516e+00,
    2.539e+00, 2.566e+00, 2.588e+00, 2.618e+00, 2.642e+00,
    2.662e+00, 2.679e+00, 2.693e+00, 2.706e+00, 2.718e+00,
    2.728e+00, 2.747e+00, 2.763e+00, 2.776e+00, 2.789e+00,
    2.810e+00, 2.828e+00, 2.850e+00, 2.868e+00, 2.894e+00,
    2.914e+00, 2.932e+00, 2.947e+00, 2.960e+00, 2.972e+00,
    2.983e+00, 2.993e+00, 3.011e+00, 3.026e+00, 3.040e+00,
    3.052e+00, 3.074e+00, 3.092e+00, 3.115e+00, 3.134e+00,
    3.160e+00, 3.182e+00, 3.201e+00, 3.217e+00, 3.232e+00,
    3.245e+00, 3.257e+00, 3.268e+00, 3.287e+00, 3.304e+00,
    3.319e+00, 3.332e+00,
        ]),
        "BREMS": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    2.570e-07, 8.393e-07, 1.439e-06, 2.690e-06, 4.004e-06,
    5.380e-06, 6.813e-06, 9.843e-06, 1.307e-05, 1.826e-05,
    2.381e-05, 3.380e-05, 4.458e-05, 5.606e-05, 6.816e-05,
    8.082e-05, 9.399e-05, 1.076e-04, 1.217e-04, 1.510e-04,
    1.817e-04, 2.137e-04, 2.468e-04, 3.160e-04, 3.886e-04,
    5.031e-04, 6.233e-04, 8.340e-04, 1.056e-03, 1.286e-03,
    1.525e-03, 1.771e-03, 2.023e-03, 2.286e-03, 2.555e-03,
    3.106e-03, 3.675e-03, 4.259e-03, 4.857e-03, 6.103e-03,
    7.392e-03, 9.393e-03, 1.146e-02, 1.506e-02, 1.879e-02,
    2.263e-02, 2.657e-02, 3.059e-02, 3.469e-02, 3.884e-02,
    4.305e-02, 5.163e-02, 6.039e-02, 6.932e-02, 7.840e-02,
    9.681e-02, 1.156e-01, 1.446e-01, 1.741e-01, 2.239e-01,
    2.748e-01, 3.266e-01, 3.792e-01, 4.325e-01, 4.865e-01,
    5.401e-01, 5.941e-01, 7.034e-01, 8.139e-01, 9.257e-01,
    1.039e+00, 1.263e+00, 1.491e+00, 1.836e+00, 2.186e+00,
    2.766e+00, 3.352e+00, 3.942e+00, 4.537e+00, 5.136e+00,
    5.738e+00, 6.334e+00, 6.932e+00, 8.132e+00, 9.339e+00,
    1.055e+01, 1.177e+01,
        ]),
        "PAIR": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 6.995e-06, 5.275e-05, 1.048e-04,
    1.929e-04, 2.913e-04, 4.743e-04, 6.770e-04, 8.963e-04,
    1.130e-03, 1.376e-03, 1.636e-03, 1.918e-03, 2.210e-03,
    2.822e-03, 3.466e-03, 4.140e-03, 4.841e-03, 6.350e-03,
    7.943e-03, 1.046e-02, 1.312e-02, 1.788e-02, 2.291e-02,
    2.817e-02, 3.361e-02, 3.923e-02, 4.499e-02, 5.079e-02,
    5.671e-02, 6.885e-02, 8.135e-02, 9.416e-02, 1.073e-01,
    1.337e-01, 1.610e-01, 2.031e-01, 2.464e-01, 3.179e-01,
    3.913e-01, 4.662e-01, 5.423e-01, 6.196e-01, 6.978e-01,
    7.752e-01, 8.532e-01, 1.011e+00, 1.171e+00, 1.333e+00,
    1.496e+00, 1.817e+00, 2.142e+00, 2.634e+00, 3.132e+00,
    3.953e+00, 4.781e+00, 5.614e+00, 6.453e+00, 7.295e+00,
    8.142e+00, 8.978e+00, 9.816e+00, 1.150e+01, 1.319e+01,
    1.488e+01, 1.658e+01,
        ]),
        "NUCL": np.array([
    5.432e-05, 5.526e-05, 5.620e-05, 5.760e-05, 5.901e-05,
    6.136e-05, 6.371e-05, 6.606e-05, 6.841e-05, 7.075e-05,
    7.310e-05, 7.545e-05, 7.780e-05, 8.250e-05, 8.719e-05,
    9.189e-05, 9.658e-05, 1.060e-04, 1.154e-04, 1.295e-04,
    1.435e-04, 1.670e-04, 1.905e-04, 2.140e-04, 2.375e-04,
    2.610e-04, 2.844e-04, 3.079e-04, 3.314e-04, 3.784e-04,
    4.253e-04, 4.723e-04, 5.193e-04, 6.132e-04, 7.071e-04,
    8.480e-04, 9.922e-04, 1.245e-03, 1.500e-03, 1.757e-03,
    2.017e-03, 2.279e-03, 2.537e-03, 2.774e-03, 3.010e-03,
    3.480e-03, 3.947e-03, 4.411e-03, 4.871e-03, 5.765e-03,
    6.649e-03, 7.957e-03, 9.249e-03, 1.140e-02, 1.352e-02,
    1.563e-02, 1.771e-02, 1.978e-02, 2.183e-02, 2.392e-02,
    2.600e-02, 3.014e-02, 3.426e-02, 3.835e-02, 4.243e-02,
    5.071e-02, 5.897e-02, 7.130e-02, 8.359e-02, 1.044e-01,
    1.252e-01, 1.461e-01, 1.669e-01, 1.877e-01, 2.085e-01,
    2.298e-01, 2.512e-01, 2.941e-01, 3.371e-01, 3.803e-01,
    4.235e-01, 5.118e-01, 6.006e-01, 7.346e-01, 8.696e-01,
    1.100e+00, 1.332e+00, 1.567e+00, 1.803e+00, 2.041e+00,
    2.280e+00, 2.524e+00, 2.769e+00, 3.264e+00, 3.764e+00,
    4.267e+00, 4.773e+00,
        ]),
        "R10": 0.6998,   # CSDA range at 10 MeV [g/cm2]
    },
    "iron": {
        "ION": np.array([
    5.494e+00, 4.817e+00, 4.321e+00, 3.783e+00, 3.399e+00,
    2.955e+00, 2.654e+00, 2.437e+00, 2.274e+00, 2.147e+00,
    2.047e+00, 1.965e+00, 1.897e+00, 1.793e+00, 1.717e+00,
    1.660e+00, 1.616e+00, 1.555e+00, 1.516e+00, 1.481e+00,
    1.463e+00, 1.452e+00, 1.452e+00, 1.458e+00, 1.467e+00,
    1.477e+00, 1.487e+00, 1.498e+00, 1.509e+00, 1.529e+00,
    1.548e+00, 1.565e+00, 1.581e+00, 1.610e+00, 1.635e+00,
    1.667e+00, 1.694e+00, 1.730e+00, 1.760e+00, 1.785e+00,
    1.806e+00, 1.824e+00, 1.841e+00, 1.855e+00, 1.868e+00,
    1.891e+00, 1.911e+00, 1.927e+00, 1.942e+00, 1.967e+00,
    1.987e+00, 2.012e+00, 2.032e+00, 2.059e+00, 2.080e+00,
    2.098e+00, 2.112e+00, 2.125e+00, 2.136e+00, 2.146e+00,
    2.155e+00, 2.171e+00, 2.184e+00, 2.196e+00, 2.207e+00,
    2.224e+00, 2.239e+00, 2.258e+00, 2.273e+00, 2.295e+00,
    2.312e+00, 2.327e+00, 2.339e+00, 2.351e+00, 2.361e+00,
    2.370e+00, 2.378e+00, 2.393e+00, 2.406e+00, 2.417e+00,
    2.428e+00, 2.446e+00, 2.461e+00, 2.480e+00, 2.496e+00,
    2.518e+00, 2.537e+00, 2.553e+00, 2.566e+00, 2.578e+00,
    2.589e+00, 2.599e+00, 2.608e+00, 2.624e+00, 2.638e+00,
    2.651e+00, 2.662e+00,
        ]),
        "BREMS": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 7.156e-07, 4.177e-06, 7.827e-06,
    1.166e-05, 1.566e-05, 2.413e-05, 3.320e-05, 4.779e-05,
    6.348e-05, 9.175e-05, 1.224e-04, 1.550e-04, 1.895e-04,
    2.257e-04, 2.633e-04, 3.023e-04, 3.426e-04, 4.265e-04,
    5.147e-04, 6.065e-04, 7.016e-04, 9.005e-04, 1.110e-03,
    1.440e-03, 1.786e-03, 2.395e-03, 3.035e-03, 3.703e-03,
    4.393e-03, 5.104e-03, 5.835e-03, 6.593e-03, 7.367e-03,
    8.958e-03, 1.060e-02, 1.228e-02, 1.401e-02, 1.758e-02,
    2.127e-02, 2.701e-02, 3.293e-02, 4.317e-02, 5.377e-02,
    6.467e-02, 7.584e-02, 8.723e-02, 9.883e-02, 1.105e-01,
    1.224e-01, 1.465e-01, 1.711e-01, 1.962e-01, 2.216e-01,
    2.730e-01, 3.254e-01, 4.058e-01, 4.879e-01, 6.255e-01,
    7.660e-01, 9.088e-01, 1.054e+00, 1.200e+00, 1.348e+00,
    1.495e+00, 1.643e+00, 1.942e+00, 2.244e+00, 2.549e+00,
    2.856e+00, 3.468e+00, 4.085e+00, 5.021e+00, 5.967e+00,
    7.533e+00, 9.112e+00, 1.070e+01, 1.230e+01, 1.391e+01,
    1.553e+01, 1.713e+01, 1.873e+01, 2.195e+01, 2.518e+01,
    2.842e+01, 3.168e+01,
        ]),
        "PAIR": np.array([
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,
    0.000e+00, 0.000e+00, 1.185e-05, 1.533e-04, 3.144e-04,
    5.876e-04, 8.929e-04, 1.461e-03, 2.091e-03, 2.773e-03,
    3.500e-03, 4.266e-03, 5.069e-03, 5.914e-03, 6.789e-03,
    8.618e-03, 1.054e-02, 1.255e-02, 1.463e-02, 1.900e-02,
    2.360e-02, 3.086e-02, 3.848e-02, 5.223e-02, 6.673e-02,
    8.185e-02, 9.752e-02, 1.136e-01, 1.302e-01, 1.467e-01,
    1.635e-01, 1.980e-01, 2.335e-01, 2.698e-01, 3.068e-01,
    3.815e-01, 4.583e-01, 5.766e-01, 6.983e-01, 8.972e-01,
    1.101e+00, 1.308e+00, 1.518e+00, 1.731e+00, 1.946e+00,
    2.158e+00, 2.370e+00, 2.799e+00, 3.233e+00, 3.670e+00,
    4.111e+00, 4.984e+00, 5.865e+00, 7.198e+00, 8.544e+00,
    1.076e+01, 1.300e+01, 1.524e+01, 1.750e+01, 1.977e+01,
    2.204e+01, 2.429e+01, 2.655e+01, 3.107e+01, 3.561e+01,
    4.016e+01, 4.472e+01,
        ]),
        "NUCL": np.array([
    4.757e-05, 4.839e-05, 4.922e-05, 5.045e-05, 5.168e-05,
    5.374e-05, 5.580e-05, 5.785e-05, 5.991e-05, 6.197e-05,
    6.402e-05, 6.608e-05, 6.814e-05, 7.225e-05, 7.636e-05,
    8.047e-05, 8.459e-05, 9.281e-05, 1.010e-04, 1.134e-04,
    1.257e-04, 1.463e-04, 1.668e-04, 1.874e-04, 2.080e-04,
    2.285e-04, 2.491e-04, 2.697e-04, 2.902e-04, 3.314e-04,
    3.725e-04, 4.136e-04, 4.548e-04, 5.370e-04, 6.193e-04,
    7.427e-04, 8.693e-04, 1.093e-03, 1.319e-03, 1.547e-03,
    1.778e-03, 2.010e-03, 2.240e-03, 2.452e-03, 2.663e-03,
    3.084e-03, 3.503e-03, 3.920e-03, 4.335e-03, 5.140e-03,
    5.937e-03, 7.119e-03, 8.289e-03, 1.024e-02, 1.216e-02,
    1.407e-02, 1.597e-02, 1.785e-02, 1.972e-02, 2.162e-02,
    2.351e-02, 2.728e-02, 3.103e-02, 3.477e-02, 3.849e-02,
    4.604e-02, 5.357e-02, 6.483e-02, 7.606e-02, 9.505e-02,
    1.140e-01, 1.330e-01, 1.520e-01, 1.710e-01, 1.900e-01,
    2.094e-01, 2.289e-01, 2.680e-01, 3.072e-01, 3.465e-01,
    3.859e-01, 4.662e-01, 5.470e-01, 6.688e-01, 7.914e-01,
    1.000e+00, 1.211e+00, 1.423e+00, 1.637e+00, 1.852e+00,
    2.068e+00, 2.288e+00, 2.510e+00, 2.956e+00, 3.406e+00,
    3.860e+00, 4.316e+00,
        ]),
        "R10": 1.0250,   # CSDA range at 10 MeV [g/cm2]
    },
}

_PROC_LOGT      = np.log(_PROC_T)
_B_RAD_ROCK_REF = 3.475e-6      # legacy constant b for rock; rescales custom materials


def _select_proc(mat):
    """
    Resolve (table, a_scale_factor, rad_scale_factor) for a material dict.
    Native materials (proc key in _PROC_TABLES) use their own table unscaled
    (or with explicit a_rescale/rad_rescale, e.g. seawater from water);
    custom materials rescale the Standard Rock table.
    """
    key = mat.get("proc")
    if key in _PROC_TABLES:
        return (_PROC_TABLES[key],
                mat.get("a_rescale", 1.0),
                mat.get("rad_rescale", 1.0))
    return (_PROC_TABLES["rock"],
            mat["a_scale"],
            mat["b_rad"] / _B_RAD_ROCK_REF)


def _loss_components(E_MeV, mat):
    """
    Per-process mean losses at kinetic energy E_MeV [MeV cm²/g]:
        returns (a_ion, L_brems, L_pair, L_nucl).
    """
    tab, s_a, s_r = _select_proc(mat)
    logE  = np.log(np.clip(np.asarray(E_MeV, float), _PROC_T[0], _PROC_T[-1]))
    a_ion = np.interp(logE, _PROC_LOGT, tab["ION"])   * s_a
    Lb    = np.interp(logE, _PROC_LOGT, tab["BREMS"]) * s_r
    Lp    = np.interp(logE, _PROC_LOGT, tab["PAIR"])  * s_r
    Ln    = np.interp(logE, _PROC_LOGT, tab["NUCL"])  * s_r
    return a_ion, Lb, Lp, Ln


# ─────────────────────────────────────────────────────────────────────────────
# Hard-event spectral shapes above v_cut (v = ΔE/E_tot)
#   brems        : Bethe-Heitler  φ(v) ∝ (1-v)/v
#   pair         : steeply falling φ(v) ∝ 1/v³   (asymptotic pair spectrum)
#   photonuclear : φ(v) ∝ 1/v
# Each process is sampled with rate λ_i = (1-v_cut)·L_i/(E_tot·<v>_i), which
# conserves the tabulated mean loss per process by construction.
# ─────────────────────────────────────────────────────────────────────────────

def _build_shape_icdf(shape_fn, v_cut, n=2000):
    """Return (v_grid, cdf, v_mean) for inverse-CDF sampling of shape_fn on [v_cut, 1]."""
    v   = np.geomspace(v_cut, 1.0, n)
    pdf = shape_fn(v)
    cdf = np.concatenate(([0.0], np.cumsum(0.5*(pdf[1:]+pdf[:-1])*np.diff(v))))
    m1  = np.concatenate(([0.0], np.cumsum(0.5*(v[1:]*pdf[1:]+v[:-1]*pdf[:-1])*np.diff(v))))
    v_mean = m1[-1] / cdf[-1]
    return v, cdf / cdf[-1], v_mean


_SHAPE_FNS = {
    "brems": lambda v: (1.0 - v) / v,
    "pair":  lambda v: v ** -3,
    "nucl":  lambda v: 1.0 / v,
}


def _proc_sampling_tables(v_cut):
    """Per-process (v_grid, cdf, v_mean) tuples for the given v_cut."""
    return {k: _build_shape_icdf(f, v_cut) for k, f in _SHAPE_FNS.items()}


# ─────────────────────────────────────────────────────────────────────────────
# δ-ray (knock-on electron) sampling — ionisation straggling
# Rutherford close-collision spectrum d²N/dT dX = k_δ/T², k_δ = 0.0767·a_scale/β²
# [MeV cm²/g].  Deltas with T > T_CUT_DELTA are sampled explicitly; their mean
# k_δ·ln(T_max/T_cut) is subtracted from the continuous ionisation term so the
# tabulated mean dE/dx is preserved exactly.  (Spin term (1-β²T/T_max) and
# Landau shape below T_cut neglected; for thick absorbers the explicit deltas
# carry the dominant variance.)
# ─────────────────────────────────────────────────────────────────────────────
T_CUT_DELTA = 10.0          # δ-ray explicit-sampling threshold [MeV]
M_E         = 0.510999      # electron mass [MeV]


def _delta_kinematics(E_kin_MeV, ZoverA):
    """Return (k_delta [MeV cm²/g], T_max [MeV], beta²) for muons at E_kin."""
    E_tot = np.asarray(E_kin_MeV, float) + M_MU
    gamma = E_tot / M_MU
    beta2 = np.clip(1.0 - 1.0/gamma**2, 1e-6, 1.0)
    bg2   = beta2 * gamma**2
    T_max = 2.0*M_E*bg2 / (1.0 + 2.0*gamma*M_E/M_MU + (M_E/M_MU)**2)
    k_del = 0.5 * 0.307075 * ZoverA / beta2
    return k_del, T_max, beta2

# ─────────────────────────────────────────────────────────────────────────────
# Bethe-Heitler (BH) hard-event inverse CDF
# dσ/dv ∝ (1-v)/v on [v_cut, 1]
# CDF: F(v) = [ln(v/v_cut) − (v − v_cut)] / C  where C = ln(1/v_cut)−(1−v_cut)
# No closed-form inverse; precomputed as a 2000-point lookup table.
# Rate: λ_BH = 2·b_total·C/(1−v_cut)  — conserves same b_hard×E as 1/v.
# ─────────────────────────────────────────────────────────────────────────────
def _build_bh_icdf(v_cut, n=2000):
    """Return (v_grid, cdf) for inverse-CDF sampling of (1-v)/v on [v_cut, 1]."""
    C      = np.log(1.0 / v_cut) - (1.0 - v_cut)
    v_grid = np.linspace(v_cut, 1.0, n)
    cdf    = (np.log(v_grid / v_cut) - (v_grid - v_cut)) / C
    cdf[0] = 0.0; cdf[-1] = 1.0
    return v_grid, cdf

_BH_DEFAULT_VCUT = 0.05
_BH_V, _BH_CDF   = _build_bh_icdf(_BH_DEFAULT_VCUT)


def _dedx(E_MeV, a_scale=1.0):
    """
    Total dE/dx [MeV cm2/g] at kinetic energy E_MeV.
    a_scale: ionisation multiplier for non-standard-rock materials
             (= (Z/A)_material / (Z/A)_rock, approximately).
    Derivation: 1/(dR/dT) from Groom (2001) range table.
    Validation:
      100 MeV  → 1.92  MeV cm2/g  (near minimum ionisation)
      1 GeV    → 1.81  MeV cm2/g  (minimum ionising)
      10 GeV   → 2.19  MeV cm2/g  (slight relativistic rise + radiative onset)
      100 GeV  → 2.75  MeV cm2/g  (radiative dominant)
      1 TeV    → 6.67  MeV cm2/g
    """
    E = np.clip(np.asarray(E_MeV, float), _T_FINE[0], _T_FINE[-1])
    return np.interp(E, _T_FINE, _DEDX_STD) * a_scale


def _csda_range(E_MeV, a_scale=1.0, t_fine=None, r_fine=None):
    """CSDA range [g/cm²] at kinetic energy E_MeV [MeV]."""
    if t_fine is None: t_fine = _T_FINE
    if r_fine is None: r_fine = _R_FINE
    E = np.clip(np.asarray(E_MeV, float), t_fine[0], t_fine[-1])
    return np.interp(E, t_fine, r_fine) * a_scale


def _det_range(E_MeV, mat, v_cut, delta_rays=False):
    """
    Deterministic-loss range [g/cm²]: range integrated with only the
    continuous components a(E) + v_cut·L_rad(E) from the PDG 2024
    per-process table (the same decomposition the transport loop uses).
    With delta_rays the continuous ionisation is the *restricted* term
    a_res (δ-rays are stochastic), keeping the bound strict.

    A muon can never lose less than the deterministic component, so this is
    a strict upper bound on penetration depth.  The pre-filter uses it so
    that only muons that cannot survive even with zero stochastic losses
    are killed outright; the mean-loss CSDA range (which includes the hard
    component) also killed muons that could have survived by avoiding hard
    events, biasing survival low near threshold.
    """
    a_ion, Lb, Lp, Ln = _loss_components(_PROC_T, mat)
    if delta_rays:
        k_del, T_max, _ = _delta_kinematics(
            _PROC_T, mat.get("ZoverA", 0.5 * mat["a_scale"]))
        m_del = np.where(T_max > 1.5 * T_CUT_DELTA,
                         k_del * np.log(np.maximum(T_max / T_CUT_DELTA, 1.0)),
                         0.0)
        a_ion = np.maximum(a_ion - m_del, 0.1 * a_ion)
    dedx_det = a_ion + v_cut * (Lb + Lp + Ln)
    inv      = 1.0 / np.maximum(dedx_det, 1e-12)
    seg      = 0.5 * (inv[1:] + inv[:-1]) * np.diff(_PROC_T)
    # Offset by the material's PDG CSDA range at the grid floor (radiative
    # losses are negligible at 10 MeV, so table ≈ deterministic there).
    tab, s_a, _ = _select_proc(mat)
    R_grid   = tab["R10"] / s_a + np.concatenate(([0.0], np.cumsum(seg)))
    E = np.clip(np.asarray(E_MeV, float), _PROC_T[0], _PROC_T[-1])
    return np.interp(E, _PROC_T, R_grid)


# ─────────────────────────────────────────────────────────────────────────────
# Material database
# b_rad [cm2/g]: total radiative coefficient (Groom 2001 Table V)
# a_scale     : ionisation scale = (Z/A)_mat / 0.500 where 0.500 = (Z/A)_rock
# X0_cm [cm]  : radiation length (PDG 2022)
# ─────────────────────────────────────────────────────────────────────────────
# "proc"   : key into _PROC_TABLES (native evaluated per-process tables)
# "ZoverA" : true <Z/A> from the PDG file header — used for δ-ray kinematics
# Seawater has no PDG table: it uses the water table with composition
# rescaling (a: Z/A ratio; radiative: b ratio) — both sub-percent.
_MAT_DB = {
    1: {"name": "Standard Rock", "Z": 11.0,  "A": 22.0,  "I_eV": 136.4,
        "X0_cm": 10.015, "b_rad": 3.475e-6, "a_scale": 1.000,
        "proc": "rock",  "ZoverA": 0.50000},  # 26.54 g/cm² / 2.65 g/cm³
    2: {"name": "Water",         "Z": 7.42,  "A": 14.18, "I_eV": 79.7,
        "X0_cm": 36.08, "b_rad": 3.20e-6,  "a_scale": 1.046,
        "proc": "water", "ZoverA": 0.55509},
    3: {"name": "Seawater",      "Z": 7.57,  "A": 14.75, "I_eV": 79.7,
        "X0_cm": 35.75, "b_rad": 3.22e-6,  "a_scale": 1.028,
        "proc": "water", "ZoverA": 0.55250,
        "a_rescale": 0.55250 / 0.55509, "rad_rescale": 3.22e-6 / 3.20e-6},
    4: {"name": "Iron",          "Z": 26.0,  "A": 55.85, "I_eV": 286.0,
        "X0_cm": 1.757, "b_rad": 4.06e-6,  "a_scale": 0.930,
        "proc": "iron",  "ZoverA": 0.46557},
    6: {"name": "Ice",           "Z": 7.42,  "A": 14.18, "I_eV": 79.7,
        "X0_cm": 39.31, "b_rad": 3.20e-6,  "a_scale": 1.046,
        "proc": "ice",   "ZoverA": 0.55509},  # 36.08 g/cm² / 0.918 g/cm³
}


def _split_dedx(E_MeV, mat, t_fine=None, dedx_fine=None):
    """
    Split total dE/dx into ionisation (a_ion) + radiative (b_total) components.
      Total  = a_ion × dx + b_soft × E_tot × dx + <hard events>
             = dE/dx × dx               (energy conserved)
    Returns (a_ion [MeV cm2/g], b_total [cm2/g]).
    t_fine / dedx_fine: fine-grid arrays to use (default: Groom 2001).
    """
    if t_fine is None:    t_fine    = _T_FINE
    if dedx_fine is None: dedx_fine = _DEDX_STD
    b_total = mat["b_rad"]
    E_tot   = np.asarray(E_MeV, float) + M_MU
    E_c     = np.clip(np.asarray(E_MeV, float), t_fine[0], t_fine[-1])
    dedx_t  = np.interp(E_c, t_fine, dedx_fine) * mat["a_scale"]
    a_ion   = dedx_t - b_total * E_tot
    a_ion   = np.maximum(a_ion, 1.5)        # floor at minimum ionising ~1.5
    return a_ion, b_total


# ─────────────────────────────────────────────────────────────────────────────
# Input file reader
# ─────────────────────────────────────────────────────────────────────────────

def _read_input(fpath, transport_all):
    """
    Parse a 13-col or 14-col CosmoALEPH surface file.
    13-col: EventID x y z p px py pz theta phi E charge det_mask
    14-col: EventID x y z p px py pz theta phi E charge hit_flag det_mask
    Returns dict of numpy arrays.  Positions in [cm], angles in [rad].
    Column E is the TOTAL energy [GeV] (generator writes E = sqrt(p²+m²));
    the returned Ekin_GeV/Ekin_MeV are the kinetic energy (E − m_mu) that
    the dE/dx and range tables are indexed by, E_tot_GeV is the raw column.

    Source-plane auto-detection
    ---------------------------
    The generator supports three source planes; the transport depth direction
    differs accordingly:
      XY plane (z = const): muons travel in ±Z → depth cosine = |pz / p|
      XZ plane (y = 0):     muons travel in ±Y → depth cosine = |py / p|
      YZ plane (x = const): muons travel in ±X → depth cosine = |px / p|
    The constant coordinate is detected from position variance (< 1 cm std).
    """
    rows = []
    with open(fpath) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append(s.split())
    if not rows:
        raise RuntimeError(f"Empty file: {fpath}")

    data = np.array(rows, dtype=float)
    nc   = data.shape[1]
    if nc not in (13, 14):
        raise RuntimeError(f"Expected 13 or 14 cols, got {nc}")

    EventID  = data[:, 0].astype(int)
    x, y, z  = data[:, 1], data[:, 2], data[:, 3]
    p_srf    = data[:, 4]
    px, py, pz = data[:, 5], data[:, 6], data[:, 7]
    theta    = data[:, 8]
    phi      = data[:, 9]
    E_tot_GeV = data[:, 10]                 # column 10 = TOTAL energy [GeV]
    charge   = data[:, 11]
    hit_flag = data[:, 12].astype(int) if nc == 14 else np.ones(len(data), dtype=int)

    if not transport_all:
        m        = (hit_flag == 1)
        EventID  = EventID[m]; x = x[m]; y = y[m]; z = z[m]
        p_srf    = p_srf[m]; px = px[m]; py = py[m]; pz = pz[m]
        theta    = theta[m];   phi = phi[m]
        E_tot_GeV = E_tot_GeV[m]; charge = charge[m]

    # Detect source plane: the coordinate with smallest spread is the constant one.
    coord_std = np.array([np.std(x), np.std(y), np.std(z)])
    if coord_std.max() < 0.1:
        # Degenerate point source (e.g. pencil beam at the origin): position
        # spread cannot identify the plane — use the dominant mean momentum
        # component as the depth axis instead.
        depth_axis = int(np.argmax([np.mean(np.abs(px)), np.mean(np.abs(py)),
                                    np.mean(np.abs(pz))]))
        print("  WARNING: point-like source (no position spread) — depth axis "
              f"taken from dominant momentum component ({'XYZ'[depth_axis]}).",
              flush=True)
    else:
        depth_axis = int(np.argmin(coord_std))    # 0=X(YZ), 1=Y(XZ), 2=Z(XY)

    p_mag = np.maximum(np.abs(p_srf), 1e-9)
    if depth_axis == 1:                            # XZ source plane, depth in Y
        depth_cos = np.abs(py) / p_mag
        _plane_name = "XZ (depth=Y)"
    elif depth_axis == 0:                          # YZ source plane, depth in X
        depth_cos = np.abs(px) / p_mag
        _plane_name = "YZ (depth=X)"
    else:                                          # XY source plane, depth in Z (default)
        depth_cos = np.abs(pz) / p_mag
        _plane_name = "XY (depth=Z)"

    print(f"  Source plane detected: {_plane_name}  (coord_std: "
          f"x={coord_std[0]:.1f}  y={coord_std[1]:.1f}  z={coord_std[2]:.1f} cm)",
          flush=True)

    # Direction cosines from momentum (preferred over theta/phi for accuracy).
    cx = px / p_mag
    cy = py / p_mag
    cz = pz / p_mag

    Ekin_GeV = np.maximum(E_tot_GeV - M_MU_GEV, 0.0)
    return dict(EventID=EventID, x=x, y=y, z=z,
                theta=theta, phi=phi,
                E_tot_GeV=E_tot_GeV,
                Ekin_GeV=Ekin_GeV, Ekin_MeV=Ekin_GeV * 1000.0,
                charge=charge, cx=cx, cy=cy, cz=cz,
                depth_cos=depth_cos,       # correct depth-direction cosine
                depth_axis=depth_axis)


# ─────────────────────────────────────────────────────────────────────────────
# Core transport
# ─────────────────────────────────────────────────────────────────────────────

def transport(muons, depth_m, rho, mat, n_steps=0, v_cut=0.05,
              ms_enable=True, rng=None, progress_cb=None,
              range_table='pdg2024', hard_spectrum='proc', delta_rays=True):
    """
    Forward stochastic muon transport through a uniform slab (v2 physics).

    Energy-loss decomposition is anchored to the PDG 2024 evaluated table
    (per-process columns), so the mean dE/dx is preserved exactly by
    construction at every energy:

        dE = a_res(E)·dx                       restricted ionisation
           + ΔE_δ                              δ-rays  (Poisson, T > 10 MeV)
           + v_cut·L_rad(E)·dx                 soft radiative (deterministic)
           + ΔE_hard                           catastrophic events (Poisson)

    with L_rad = L_brems + L_pair + L_nucl from the table; hard events are
    sampled per process with shape-specific spectra (hard_spectrum='proc').

    Parameters
    ----------
    muons      : dict from _read_input
    depth_m    : vertical overburden [m]
    rho        : density [g/cm³]
    mat        : material dict from _MAT_DB
    n_steps    : integration steps per muon; 0 = auto (per-muon adaptive:
                 dx ≈ 5 g/cm², clipped to [300, 20000] steps per muon)
    v_cut      : catastrophic event threshold (fraction of E_tot)
    ms_enable    : Highland multiple scattering deflections
    rng          : np.random.default_rng
    progress_cb  : callable(n_transported, n_survived, n_total)
    range_table  : retained for CLI compatibility; since v2 the loss model is
                   always anchored to the PDG 2024 per-process table
    hard_spectrum: 'proc' (per-process: brems (1-v)/v, pair 1/v³, nucl 1/v;
                   recommended), 'bh' (single (1-v)/v), 'groom' (single 1/v)
    delta_rays   : explicit δ-ray sampling above T_CUT_DELTA with restricted
                   continuous ionisation (ionisation straggling)

    Returns
    -------
    dict with alive, E_kin_f_MeV, cx_f, cy_f, cz_f, x_f, y_f, z_f, theta_f, phi_f
    """
    if rng is None:
        rng = np.random.default_rng()

    a_scale   = mat["a_scale"]
    ZoverA    = mat.get("ZoverA", 0.5 * a_scale)   # true <Z/A> for δ-ray kinematics
    X0_gcm2   = mat["X0_cm"] * rho             # radiation length [g/cm²]
    ln_ivc    = np.log(1.0 / v_cut)

    # Hard-event sampling tables.  All shapes conserve the tabulated mean:
    # λ_i·<v>_i·E_tot·dx = (1-v_cut)·L_i·dx  per process (or for the total).
    if hard_spectrum == 'proc':
        _shp = _proc_sampling_tables(v_cut)    # {name: (v_grid, cdf, v_mean)}
    elif hard_spectrum == 'bh':
        if abs(v_cut - _BH_DEFAULT_VCUT) < 1e-10:
            _bh_v, _bh_cdf = _BH_V, _BH_CDF
        else:
            _bh_v, _bh_cdf = _build_bh_icdf(v_cut)
        _C_bh    = ln_ivc - (1.0 - v_cut)
        _vmean_bh = (1.0 - v_cut)**2 / (2.0 * _C_bh)
    else:                                       # 'groom' 1/v
        _vmean_1v = (1.0 - v_cut) / ln_ivc

    N         = len(muons["Ekin_MeV"])
    d_cm      = depth_m * 100.0
    alive     = np.ones(N, dtype=bool)

    # Per-muon slant path [cm] and opacity [g/cm²].
    # Use the depth-direction cosine auto-detected from the source plane.
    depth_cos  = np.maximum(muons.get("depth_cos", np.abs(muons["cz"])), 0.02)
    slant_cm   = d_cm / depth_cos
    slant_gcm2 = slant_cm * rho

    # Pre-filter on the deterministic-loss range (strict upper bound on
    # penetration): only muons that cannot survive even with zero hard
    # radiative events are killed outright.  Muons between the mean-loss
    # CSDA range and this bound are transported stochastically.
    R_det   = _det_range(muons["Ekin_MeV"], mat, v_cut, delta_rays=delta_rays)
    prefilt = R_det < slant_gcm2
    alive[prefilt] = False

    # Per-muon adaptive step count: target dx ≈ 5 g/cm² for every muon (the
    # old global count, set from the median slant, gave near-horizontal muons
    # steps up to ~50× larger).  Work scales with Σ n_steps_i because only
    # muons with remaining steps are processed each iteration.
    DX_TARGET = 20.0                            # [g/cm²]  (5->20: ~2.5x faster,
    # survival/exit-KE unchanged to <0.02 pp at 25/100/200 m; step physics is
    # compound-Poisson and quadrature-MCS, invariant in the mean to step size;
    # forward-Euler dE/dx verified accurate to <0.1 GeV at this step)
    if n_steps > 0:
        n_steps_i = np.full(N, int(n_steps), dtype=np.int64)
    else:
        n_steps_i = np.clip(np.ceil(slant_gcm2 / DX_TARGET),
                            300, 20000).astype(np.int64)
    n_steps_max = int(n_steps_i.max()) if N else 0
    tot_steps   = float(n_steps_i.sum()) if N else 1.0

    dx    = slant_gcm2 / n_steps_i              # per-muon step size [g/cm²]
    E_cur = muons["Ekin_MeV"].copy()
    cx_c  = muons["cx"].copy()
    cy_c  = muons["cy"].copy()
    cz_c  = muons["cz"].copy()

    # Lateral position accumulators (step-by-step integration)
    # x_acc, y_acc track how far each muon has drifted laterally.
    # Updated BEFORE MCS so each step uses the direction at entry into that step.
    depth_axis = muons.get("depth_axis", 2)
    x_acc = np.zeros(N)
    y_acc = np.zeros(N)
    z_acc = np.zeros(N)   # z displacement — stopping depth + non-XY sources

    # Pre-filtered muons never enter the stepping loop; give them their CSDA
    # stopping displacement so z_stop reports a physical depth instead of 0.
    # (min with R_det: the deterministic-only range is an upper bound.)
    if prefilt.any():
        R_stop = np.minimum(
            _csda_range(muons["Ekin_MeV"], mat["a_scale"],
                        _PDG24_T_FINE, _PDG24_R_FINE),
            R_det)
        z_acc[prefilt] = muons["cz"][prefilt] * (R_stop[prefilt] / rho)

    E_stop = 1.0                                # stop if Ekin < 1 MeV
    prev_reported = 0

    for step in range(n_steps_max):
        act = alive & (step < n_steps_i)
        if not act.any():
            break

        idx   = np.where(act)[0]
        E_i   = E_cur[idx]
        E_tot = E_i + M_MU
        dx_i  = dx[idx]

        # ── Accumulate position BEFORE MCS (direction at entry of this step) ──
        # Convert dx [g/cm²] → geometric step [cm] along slant direction.
        # z is accumulated for every depth axis: survivors on an XY source are
        # snapped to the exact scoring plane below, but stopped muons need the
        # integrated z to report their true stopping depth.
        dx_cm_i = dx_i / rho
        x_acc[idx] += cx_c[idx] * dx_cm_i
        y_acc[idx] += cy_c[idx] * dx_cm_i
        z_acc[idx] += cz_c[idx] * dx_cm_i

        # Table-anchored loss decomposition at the current energy
        a_ion, Lb, Lp, Ln = _loss_components(E_i, mat)
        L_tot = Lb + Lp + Ln

        # 1. Ionisation: restricted continuous term + explicit δ-rays
        if delta_rays:
            k_del, T_max, beta2 = _delta_kinematics(E_i, ZoverA)
            has    = T_max > 1.5 * T_CUT_DELTA
            m_del  = np.where(has,
                              k_del * (np.log(np.maximum(T_max / T_CUT_DELTA, 1.0))
                                       - beta2 * (T_max - T_CUT_DELTA)
                                         / np.maximum(T_max, T_CUT_DELTA)),
                              0.0)
            a_res  = np.maximum(a_ion - m_del, 0.1 * a_ion)
            dE_ion = a_res * dx_i
            lam_d  = np.where(has,
                              k_del * (1.0 / T_CUT_DELTA
                                       - 1.0 / np.maximum(T_max, T_CUT_DELTA)),
                              0.0)
            n_d    = rng.poisson(lam_d * dx_i)
            md1 = (n_d == 1)
            if md1.any():
                U   = rng.uniform(size=md1.sum())
                Ts  = T_CUT_DELTA / (1.0 - U * (1.0 - T_CUT_DELTA / T_max[md1]))
                acc = rng.uniform(size=md1.sum()) < (1.0 - beta2[md1] * Ts / T_max[md1])
                dE_ion[md1] += np.where(acc, Ts, 0.0)
            for mi in np.where(n_d >= 2)[0]:
                U   = rng.uniform(size=int(n_d[mi]))
                Ts  = T_CUT_DELTA / (1.0 - U * (1.0 - T_CUT_DELTA / T_max[mi]))
                acc = rng.uniform(size=int(n_d[mi])) < (1.0 - beta2[mi] * Ts / T_max[mi])
                dE_ion[mi] += float((Ts * acc).sum())
        else:
            dE_ion = a_ion * dx_i

        # 2. Soft radiative (deterministic, v < v_cut)
        dE_soft = ((v_cut * L_tot * dx_i) if hard_spectrum != 'proc'
                   else (((2.0*v_cut - v_cut**2)*Lb + v_cut*Ln + Lp) * dx_i))

        # 3. Hard radiative (Poisson-sampled, v > v_cut)
        dE_hard = np.zeros(len(idx))
        if hard_spectrum == 'proc':
            # Per-process sampling: brems (1-v)/v, pair 1/v³, photonuclear 1/v
            for L_proc, key, hf in ((Lb, "brems", (1.0-v_cut)**2), (Ln, "nucl", (1.0-v_cut))):
                pv, pcdf, pvm = _shp[key]
                lam_i = hf * L_proc / (E_tot * pvm)  # events/(g/cm²)
                n_ev  = rng.poisson(lam_i * dx_i)
                m1 = (n_ev == 1)
                if m1.any():
                    v_samp = np.interp(rng.uniform(size=m1.sum()), pcdf, pv)
                    dE_hard[m1] += v_samp * E_tot[m1]
                for mi in np.where(n_ev >= 2)[0]:
                    U  = rng.uniform(size=int(n_ev[mi]))
                    vs = np.interp(U, pcdf, pv)
                    dE_hard[mi] += min(float(vs.sum()), 1.0) * E_tot[mi]
        else:
            # Single-shape sampling with energy-dependent total radiative loss
            if hard_spectrum == 'bh':
                lam = (1.0 - v_cut) * L_tot / (E_tot * _vmean_bh)
            else:
                lam = (1.0 - v_cut) * L_tot / (E_tot * _vmean_1v)
            n_ev = rng.poisson(lam * dx_i)

            m1 = (n_ev == 1)
            if m1.any():
                U = rng.uniform(size=m1.sum())
                if hard_spectrum == 'bh':
                    v_samp = np.interp(U, _bh_cdf, _bh_v)
                else:
                    v_samp = v_cut ** (1.0 - U)         # inverse CDF of 1/v
                dE_hard[m1] = np.minimum(v_samp, 1.0) * E_tot[m1]

            for mi in np.where(n_ev >= 2)[0]:
                Er = float(E_tot[mi])
                for _ in range(int(n_ev[mi])):
                    Uu = float(rng.uniform())
                    if hard_spectrum == 'bh':
                        v = float(np.interp(Uu, _bh_cdf, _bh_v))
                    else:
                        v = min(v_cut ** (1.0 - Uu), 1.0)
                    dE_hard[mi] += v * Er
                    Er = max(Er * (1.0 - v), M_MU)
                    if Er <= M_MU + E_stop:
                        dE_hard[mi] = E_tot[mi]; break

        # 4. Update energy
        E_new       = E_i - (dE_ion + dE_soft + dE_hard)
        E_cur[idx]  = np.maximum(E_new, 0.0)
        alive[idx[E_new < E_stop]] = False

        # 5. Muon decay  (p = √(E² − m²), not p ≈ E)
        la = alive[idx]
        if la.any():
            E_t_la  = E_tot[la]
            p_GeV   = np.sqrt(np.maximum(E_t_la**2 - M_MU**2, 1e-6)) / 1000.0
            dec_len = p_GeV * PCTAU_CM * rho    # [g/cm²]
            P_dec   = -np.expm1(-dx_i[la] / np.maximum(dec_len, 1.0))
            alive[idx[la][rng.uniform(size=la.sum()) < P_dec]] = False

        # 6. Highland multiple scattering
        if ms_enable:
            la2 = alive[idx]
            if la2.any():
                oi      = idx[la2]
                dx_ms   = dx_i[la2]
                p_GeV2  = E_tot[la2] / 1000.0
                beta    = p_GeV2 / np.sqrt(p_GeV2**2 + M_MU_GEV**2)
                t_X0    = dx_ms / X0_gcm2
                theta0  = (13.6e-3 / (beta * p_GeV2)
                           * np.sqrt(t_X0)
                           * (1.0 + 0.038 * np.log(np.maximum(t_X0, 1e-12))))

                phi_az  = rng.uniform(0.0, 2.0*np.pi, size=la2.sum())
                # Rayleigh(theta0) polar deflection: Highland theta0 is the RMS
                # projected angle; a Gaussian polar angle under-scatters by sqrt(2).
                dth     = rng.rayleigh(theta0)
                sin_dth = np.sin(dth); cos_dth = np.cos(dth)
                cos_phi = np.cos(phi_az); sin_phi = np.sin(phi_az)

                cx_o = cx_c[oi]; cy_o = cy_c[oi]; cz_o = cz_c[oi]

                # Build orthonormal basis (e1, e2) perpendicular to d so the
                # deflection azimuth phi_az is uniform around the track — same
                # construction as ucmuon_bb_driver.transport_bb.
                # Reference vector: y-axis if |cx| > 0.9, else x-axis
                near_x  = np.abs(cx_o) > 0.9
                wx = np.where(near_x, 0.0, 1.0)
                wy = np.where(near_x, 1.0, 0.0)
                wd = wx * cx_o + wy * cy_o          # w·d  (wz=0)
                e1x = wx - wd * cx_o
                e1y = wy - wd * cy_o
                e1z =    - wd * cz_o
                ne1 = np.maximum(np.sqrt(e1x**2 + e1y**2 + e1z**2), 1e-15)
                e1x /= ne1;  e1y /= ne1;  e1z /= ne1
                e2x = cy_o * e1z - cz_o * e1y      # e2 = d × e1
                e2y = cz_o * e1x - cx_o * e1z
                e2z = cx_o * e1y - cy_o * e1x

                npx = cos_phi * e1x + sin_phi * e2x
                npy = cos_phi * e1y + sin_phi * e2y
                npz = cos_phi * e1z + sin_phi * e2z

                cx_c[oi] = cx_o * cos_dth + npx * sin_dth
                cy_c[oi] = cy_o * cos_dth + npy * sin_dth
                cz_c[oi] = cz_o * cos_dth + npz * sin_dth

                norm = np.sqrt(cx_c[oi]**2 + cy_c[oi]**2 + cz_c[oi]**2)
                norm = np.maximum(norm, 1e-15)
                cx_c[oi] /= norm; cy_c[oi] /= norm; cz_c[oi] /= norm

        # 7. Progress  (fraction of total per-muon steps completed)
        if progress_cb:
            n_done = int(np.minimum(step + 1, n_steps_i).sum() / tot_steps * N)
            if n_done - prev_reported >= REPORT_EVERY:
                # Scale alive count by progress fraction so survived <= transported.
                # The final dedicated print (after the loop) overwrites with the
                # exact alive count, so the end-of-run metric is always correct.
                n_surv_est = int(n_done * alive.sum() / N) if N > 0 else 0
                progress_cb(n_done, n_surv_est, N)
                prev_reported = n_done

    # Final exit positions from step-by-step accumulated lateral drift
    x_f = muons["x"] + x_acc
    y_f = muons["y"] + y_acc
    if depth_axis != 2:
        z_f = muons["z"] + z_acc
    else:
        z_f = np.full(N, -d_cm)    # depth-axis=Z: z is the exact scoring plane
    theta_f = np.arccos(np.clip(-cz_c, -1.0, 1.0))
    phi_f   = np.arctan2(cy_c, cx_c)

    # For stopped muons z_f = -d_cm (plane depth) is wrong; use the actual
    # integrated position — same convention as the PROPOSAL/MUSIC drivers.
    alive_arr = alive.astype(int)
    z_stop = np.where(alive_arr == 0, muons["z"] + z_acc, z_f)

    return dict(
        alive=alive_arr,
        E_kin_f_MeV=np.maximum(E_cur, 0.0),
        cx_f=cx_c, cy_f=cy_c, cz_f=cz_c,
        x_f=x_f, y_f=y_f, z_f=z_f,
        z_stop=z_stop,                               # stopping depth for dead muons
        theta_f=theta_f, phi_f=phi_f,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Parallel transport (multiprocessing) — muons are embarrassingly parallel
# ─────────────────────────────────────────────────────────────────────────────

def _slice_muons(muons, sl):
    """Slice every per-muon array in the muons dict; pass scalars through."""
    return {k: (v[sl] if isinstance(v, np.ndarray) else v)
            for k, v in muons.items()}


def _transport_chunk(args):
    """Top-level worker: transport one chunk with an independent RNG stream."""
    (idx, mu_chunk, depth_m, rho, mat, n_steps, v_cut, ms_enable,
     child_seed, range_table, hard_spectrum, delta_rays) = args
    rng = np.random.default_rng(child_seed)
    res = transport(mu_chunk, depth_m, rho, mat,
                    n_steps=n_steps, v_cut=v_cut, ms_enable=ms_enable,
                    rng=rng, progress_cb=None,
                    range_table=range_table, hard_spectrum=hard_spectrum,
                    delta_rays=delta_rays)
    return idx, len(mu_chunk["Ekin_MeV"]), int(res["alive"].sum()), res


def transport_parallel(muons, depth_m, rho, mat, n_workers, seed=42,
                       n_steps=0, v_cut=0.05, ms_enable=True,
                       progress_cb=None, range_table='pdg2024',
                       hard_spectrum='proc', delta_rays=True):
    """
    Transport split across `n_workers` processes (same physics as transport()).

    Each chunk gets an independent SeedSequence-spawned RNG stream, so results
    are reproducible for a given (seed, n_workers).  Per-muon adaptive step
    counts depend only on each muon's own slant path, so chunking does not
    change the stepping.
    """
    import multiprocessing as mp

    N = len(muons["Ekin_MeV"])

    n_chunks = min(max(n_workers * 4, n_workers), N)
    bounds   = np.linspace(0, N, n_chunks + 1).astype(int)
    children = np.random.SeedSequence(seed).spawn(n_chunks)
    tasks    = [(i, _slice_muons(muons, slice(bounds[i], bounds[i + 1])),
                 depth_m, rho, mat, n_steps, v_cut, ms_enable,
                 children[i], range_table, hard_spectrum, delta_rays)
                for i in range(n_chunks)]

    results = [None] * n_chunks
    done_n = surv_n = 0
    with mp.Pool(processes=n_workers) as pool:
        for idx, n_c, n_s, res in pool.imap_unordered(_transport_chunk, tasks):
            results[idx] = res
            done_n += n_c
            surv_n += n_s
            if progress_cb:
                progress_cb(done_n, surv_n, N)

    return {k: np.concatenate([r[k] for r in results]) for k in results[0]}


# ─────────────────────────────────────────────────────────────────────────────
# Output writer  (18-col MUSIC underground format)
# ─────────────────────────────────────────────────────────────────────────────

def _write_output(muons, result, fpath):
    with open(fpath, "w") as fh:
        fh.write("# UCMuon-MC engine — UCLouvain Muography Group\n")
        fh.write("# E_convention: total_energy_GeV   (E = KE + 0.10566)\n")
        fh.write("# Cols: EventID xs ys zs Es[GeV] thetas phis charge alive"
                 " x y z E[GeV] cx cy cz theta phi\n")
        # Surface column Es is total energy; fall back to kinetic + m for
        # callers that build the muons dict without the raw file column.
        Es_arr = muons.get("E_tot_GeV", muons["Ekin_GeV"] + M_MU_GEV)
        for i in range(len(muons["Ekin_GeV"])):
            alive_i = int(result["alive"][i])
            Es_GeV  = Es_arr[i]
            if alive_i:
                x_i  = result["x_f"][i]
                y_i  = result["y_f"][i]
                z_i  = result["z_f"][i]
                E_i  = result["E_kin_f_MeV"][i] / 1000.0 + M_MU_GEV  # total energy
                cx_i = result["cx_f"][i]
                cy_i = result["cy_f"][i]
                cz_i = result["cz_f"][i]
                th_i = result["theta_f"][i]
                ph_i = result["phi_f"][i]
            else:
                # Spec: alive=0 → x=xs, y=ys, z=stop_depth, E=0,
                #                  cx=0, cy=0, cz=-1, theta=0, phi=0
                x_i  = muons["x"][i]
                y_i  = muons["y"][i]
                z_i  = result["z_stop"][i]
                E_i  = 0.0
                cx_i, cy_i, cz_i = 0.0, 0.0, -1.0
                th_i, ph_i = 0.0, 0.0
            fh.write(
                f"{muons['EventID'][i]:10d}"
                f" {muons['x'][i]:13.4f} {muons['y'][i]:13.4f} {muons['z'][i]:13.4f}"
                f" {Es_GeV:13.6f}"
                f" {muons['theta'][i]:13.6f} {muons['phi'][i]:13.6f}"
                f" {int(muons['charge'][i]):4d}"
                f" {alive_i:4d}"
                f" {x_i:13.4f} {y_i:13.4f} {z_i:13.4f}"
                f" {E_i:13.6f}"
                f" {cx_i:10.6f} {cy_i:10.6f} {cz_i:10.6f}"
                f" {th_i:13.6f} {ph_i:13.6f}\n"
            )


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    lines = [ln.strip() for ln in sys.stdin if ln.strip() and not ln.strip().startswith("#")]

    def _rd(i, default, typ=str):
        try:    return typ(lines[i]) if i < len(lines) else default
        except: return default

    infile        = _rd(0,  "muons_selected.dat")
    outfile       = _rd(1,  "muons_stochastic.dat")
    depth_m       = _rd(2,  500.0,    float)
    rho           = _rd(3,  2.65,     float)
    X0_cm         = _rd(4,  10.015,   float)
    mat_id        = _rd(5,  1,        int)
    transport_all = bool(_rd(6, 0,    int))
    _ncols        = _rd(7,  13,       int)
    n_steps       = _rd(8,  0,        int)
    v_cut         = _rd(9,  0.05,     float)
    ms_enable     = bool(_rd(10, 1,   int))
    Z_eff         = _rd(11, 11.0,     float)
    A_eff         = _rd(12, 22.0,     float)
    I_eV          = _rd(13, 136.4,    float)
    b_rad_cust    = _rd(14, 3.475e-6, float)
    range_table_id = _rd(15, 1,       int)   # 0=groom2001  1=pdg2024 (legacy; v2 loss model is PDG-anchored)
    hard_spec_id   = _rd(16, 2,       int)   # 0=groom (1/v)  1=bh ((1-v)/v)  2=per-process (default)
    range_table    = 'pdg2024'  if range_table_id == 1 else 'groom2001'
    hard_spectrum  = {0: 'groom', 1: 'bh'}.get(hard_spec_id, 'proc')
    seed           = _rd(17, 42, int)
    n_workers      = _rd(18, 1,  int)   # param 19: 1=serial  0=auto  >1 explicit
    delta_rays     = bool(_rd(19, 1, int))   # param 20: δ-ray straggling (default on)

    if mat_id == 5:
        mat = {"name": "Custom", "Z": Z_eff, "A": A_eff, "I_eV": I_eV,
               "X0_cm": X0_cm, "b_rad": b_rad_cust,
               "a_scale": (Z_eff / A_eff) / 0.5,
               "ZoverA": Z_eff / A_eff, "proc": None}
    else:
        mat = dict(_MAT_DB.get(mat_id, _MAT_DB[1]))
        mat["X0_cm"] = X0_cm

    print(f"  UCMuon-MC engine v2.0 — forward stochastic transport", flush=True)
    print(f"  {mat['name']}  rho={rho:.3f} g/cm3  X0={mat['X0_cm']:.2f} cm"
          f"  depth={depth_m:.1f} m  X={depth_m*100*rho:.1f} g/cm2", flush=True)
    print(f"  loss model: PDG 2024 per-process tables (ion/brems/pair/photonuc)"
          f"  v_cut={v_cut:.3f}  MS={'ON' if ms_enable else 'OFF'}", flush=True)
    print(f"  hard_spectrum={hard_spectrum}"
          f"  delta_rays={'ON (T>%.0f MeV)' % T_CUT_DELTA if delta_rays else 'OFF'}",
          flush=True)

    if not Path(infile).exists():
        print(f"  ERROR: input not found: {infile}", flush=True)
        sys.exit(1)

    muons = _read_input(infile, transport_all)
    N     = len(muons["Ekin_MeV"])
    if N == 0:
        print("  ERROR: no muons in input file.", flush=True); sys.exit(1)

    # Resolve worker count: auto = one worker per ~20k muons, up to all cores.
    if n_workers == 0:
        import os
        n_workers = min(os.cpu_count() or 1, max(1, N // 20000))
    n_workers = max(1, min(n_workers, N))

    print(f"  Read {N:,} muons from {infile}", flush=True)
    print(f"  Steps: {n_steps if n_steps > 0 else 'auto (>=300)'}", flush=True)
    print(f"  Workers: {n_workers}{' (serial)' if n_workers == 1 else ''}", flush=True)
    print(f"  Transported:     0  Survived:     0  Total: {N}", flush=True)

    def _cb(nd, ns, nt):
        print(f"  Transported: {nd}  Survived: {ns}  Total: {nt}", flush=True)

    # Seed 42 by default ensures reproducible output; override via stdin param 18.
    t_start = time.perf_counter()
    if n_workers > 1:
        result = transport_parallel(muons, depth_m, rho, mat, n_workers,
                                    seed=seed, n_steps=n_steps, v_cut=v_cut,
                                    ms_enable=ms_enable, progress_cb=_cb,
                                    range_table=range_table,
                                    hard_spectrum=hard_spectrum,
                                    delta_rays=delta_rays)
    else:
        rng    = np.random.default_rng(seed)
        result = transport(muons, depth_m, rho, mat,
                           n_steps=n_steps, v_cut=v_cut,
                           ms_enable=ms_enable, rng=rng, progress_cb=_cb,
                           range_table=range_table, hard_spectrum=hard_spectrum,
                           delta_rays=delta_rays)
    elapsed = time.perf_counter() - t_start

    n_surv = int(result["alive"].sum())
    print(f"  Transported: {N}  Survived: {n_surv}  Total: {N}", flush=True)
    print(f"  Survival rate: {100.0*n_surv/N:.4f} %", flush=True)
    print(f"  Elapsed: {elapsed:.1f} s", flush=True)

    _write_output(muons, result, outfile)
    print(f"  Output: {outfile}  ({N} rows, {n_surv} survived)", flush=True)

    # Timing file
    out_path    = Path(outfile)
    timing_file = str(out_path.with_suffix("")) + "_timing.txt"
    with open(timing_file, "w") as ft:
        ft.write(f"Elapsed : {elapsed:.1f}\n")
    print(f"  Timing: {timing_file}", flush=True)

    # Stopped muon file
    stopped_file = str(out_path.with_suffix("")) + "_stopped.dat"
    alive_arr    = result["alive"]
    z_stop       = result["z_stop"]
    mask_stopped = (alive_arr == 0)
    n_stopped    = int(mask_stopped.sum())
    with open(stopped_file, "w") as fs:
        fs.write(f"# UCMuon-MC stopped muons  depth={depth_m:.2f} m\n")
        fs.write("# EventID  InitKE_GeV  StopDepth_cm\n")
        evids   = muons["EventID"][mask_stopped].astype(int)
        initKEs = muons["Ekin_GeV"][mask_stopped]
        stop_zs = np.abs(z_stop[mask_stopped])
        lines   = [f"{evid:10d}  {ke:13.6f}  {sz:13.4f}\n"
                   for evid, ke, sz in zip(evids, initKEs, stop_zs)]
        fs.write("".join(lines))
    print(f"  Stopped muons: {stopped_file}  ({n_stopped} entries)", flush=True)


if __name__ == "__main__":
    main()
