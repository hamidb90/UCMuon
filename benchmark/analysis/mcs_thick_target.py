#!/usr/bin/env python3
"""Thick-target limit of multiple Coulomb scattering, the value UCMuon-MC and
the Bethe-Bloch engines use beyond T_SAT_X0 (ES_THICK_MEV in
gui/ucmuon_stochastic_driver.py).

Single scatterings are cut off at the angle hbar/(pR) by the size of the
nucleus. Once the multiple-scattering width is far above it (about 10
radiation lengths of rock), the angular distribution is Gaussian with the
second moment of the single-scattering cross section, which is finite and grows
linearly with the thickness: theta0 = E_s/(beta c p) sqrt(x/X0) with a constant
E_s, and Highland's log term, a fit to layers of 1e-3 to 100 X0, stops growing.

Per element, screened Rutherford (Moliere screening angle) times a Gaussian
nuclear form factor (rms radius sqrt(3/5) 1.2 A^(1/3) fm; the proton's 0.84 fm),
with Z(Z+1) for the atomic electrons:

    <theta_proj^2> per g/cm2 = sum_i w_i (N_A/A_i) 4 pi Z_i (Z_i + 1) r_e^2 (m_e c^2/p)^2 L_i
    L_i = ln( sqrt(3) a_TF / (R sqrt(1.13 + 3.76 (alpha Z)^2)) ) - (1 + gamma_E)/2

(the closed form of the integral for the Gaussian form factor; a numerical
integral with a uniform-sphere form factor differs by 0.2 %, a nuclear radius
15 % larger or smaller by 0.8 %). E_s = sqrt(X0 <theta_proj^2>) p.

    python3 benchmark/analysis/mcs_thick_target.py
"""
from __future__ import annotations

import numpy as np

ALPHA, R_E, ME_C2, N_A, A0_FM = 1 / 137.036, 2.8179403e-13, 0.51099895, 6.02214e23, 52917.72
GAMMA_E = 0.5772156649

MATERIALS = {   # (Z, A, mass fraction) and X0 [g/cm2] (PDG)
    "Standard Rock": ([(11, 22.0, 1.0)], 26.54),
    "water": ([(1, 1.008, 0.1119), (8, 15.999, 0.8881)], 36.08),
    "quartz (SiO2)": ([(14, 28.086, 0.4674), (8, 15.999, 0.5326)], 27.05),
    "limestone (CaCO3)": ([(20, 40.078, 0.4004), (6, 12.011, 0.1200), (8, 15.999, 0.4796)], 24.03),
    "iron": ([(26, 55.845, 1.0)], 13.84),
    "lead": ([(82, 207.2, 1.0)], 6.37),
}


def log_term(Z, A):
    a_tf = 0.88534 * A0_FM * Z ** (-1 / 3)                         # Thomas-Fermi radius [fm]
    r_rms = np.sqrt(3 / 5) * 1.2 * A ** (1 / 3) if A > 1.5 else 0.84
    return (np.log(np.sqrt(3) * a_tf / (r_rms * np.sqrt(1.13 + 3.76 * (ALPHA * Z) ** 2)))
            - (1 + GAMMA_E) / 2)


def es_thick_mev(elements, x0_gcm2):
    var = sum(w * N_A / A * 4 * np.pi * Z * (Z + 1) * R_E ** 2 * ME_C2 ** 2 * log_term(Z, A)
              for Z, A, w in elements)
    return float(np.sqrt(x0_gcm2 * var))


def main() -> int:
    for name, (elements, x0) in MATERIALS.items():
        es = es_thick_mev(elements, x0)
        print(f"{name:18s} E_s,thick = {es:5.2f} MeV  (Highland log saturating at "
              f"{np.exp((es / 13.6 - 1) / 0.038):4.1f} X0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
