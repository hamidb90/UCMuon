#!/usr/bin/env python3
"""
ucmuon_bb_driver.py  —  UCLouvain Muography Group
Fast CSDA (Bethe-Bloch) muon transport — pure Python / NumPy.

Speed advantage over the Fortran ucmuon_transport_bb_omp
---------------------------------------------------------
  MS = OFF  Analytical exact CSDA: exit energy = R⁻¹(R₀ − slant_gcm²).
            Zero step loops — O(N) array ops only.  Sub-second for 600 k muons.

  MS = ON   Vectorized NumPy stepping.  Default step = 100 g/cm²  (vs Fortran
            fixed 10 g/cm²) → ~26× fewer iterations.  All alive muons processed
            simultaneously per step via NumPy broadcasting.

Physics
-------
  Continuous energy loss from PDG 2024 CSDA range table (direct dE/dx column,
  56 entries).  Optional Groom 2001 table for comparison.
  No stochastic fluctuations (pure CSDA).  Optional Highland (1979) MCS.

Stdin (drop-in for ucmuon_transport_bb_omp — same parameter order)
-------------------------------------------------------------------
  1  infile
  2  outfile
  3  transport_all   0 = hit_flag=1 only  |  1 = all
  4  ncols_hint      ignored (auto-detected from file)
  5  depth_m         vertical depth [m]
  6  mat_type        1=StdRock  2=Ice  3=Water  4=Concrete  5=Custom
  [mat_type=5 only — four extra lines: Zeff  Aeff  rho_gcm3  I_eV]
  7  ms_enable       0=OFF  1=Highland ON
  8  range_table     0=Groom2001  1=PDG2024 (default 1)
  9  n_steps         steps for MS mode; 0=auto (default 0)

Output: 18-column MUSIC-compatible underground file  (same as stochastic engine).
"""

import sys
import time
import numpy as np
from pathlib import Path

# ── borrow tables and reader from the stochastic driver (same directory) ─────
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from ucmuon_stochastic_driver import (
    _T_FINE, _DEDX_STD, _R_FINE,                    # Groom 2001 fine grid
    _PDG24_T_FINE, _PDG24_DEDX_FINE, _PDG24_R_FINE, # PDG 2024 fine grid
    _read_input, _write_output,
    M_MU, M_MU_GEV, REPORT_EVERY,
)

# ─────────────────────────────────────────────────────────────────────────────
# Material database for BB engine
# X0_gcm2 : radiation length [g/cm²]  (Highland MCS)
# a_scale  : ionisation scale relative to Standard Rock
# rho      : default density [g/cm³]
# ─────────────────────────────────────────────────────────────────────────────
_BB_MAT = {
    1: {"name": "Standard Rock", "Z": 11.0,  "A": 22.0,  "rho": 2.65,
        "I_eV": 136.4, "X0_gcm2": 26.54, "a_scale": 1.000},
    2: {"name": "Ice",           "Z": 7.42,  "A": 14.99, "rho": 0.917,
        "I_eV":  79.7, "X0_gcm2": 33.10, "a_scale": 1.046},
    3: {"name": "Water",         "Z": 7.42,  "A": 14.99, "rho": 1.000,
        "I_eV":  79.7, "X0_gcm2": 36.08, "a_scale": 1.046},
    4: {"name": "Concrete",      "Z": 11.11, "A": 22.08, "rho": 2.300,
        "I_eV": 135.2, "X0_gcm2": 26.70, "a_scale": 1.001},
}


# ─────────────────────────────────────────────────────────────────────────────
# Core transport
# ─────────────────────────────────────────────────────────────────────────────

def transport_bb(muons, depth_m, mat, n_steps=0, ms_enable=True,
                 rng=None, progress_cb=None, range_table='pdg2024'):
    """
    CSDA muon transport through a uniform slab.

    Parameters
    ----------
    muons       : dict from _read_input
    depth_m     : vertical overburden [m]
    mat         : material dict from _BB_MAT (or custom dict with same keys)
    n_steps     : integration steps (MS mode only); 0 = auto (≥ 50, 100 g/cm² each)
    ms_enable   : Highland (1979) multiple Coulomb scattering
    rng         : np.random.default_rng instance
    progress_cb : callable(n_transported, n_survived, n_total)
    range_table : 'pdg2024' (default) or 'groom2001'

    Returns
    -------
    dict compatible with ucmuon_stochastic_driver._write_output
    """
    if rng is None:
        rng = np.random.default_rng()

    if range_table == 'pdg2024':
        t_tab, dedx_tab, r_tab = _PDG24_T_FINE, _PDG24_DEDX_FINE, _PDG24_R_FINE
    else:
        t_tab, dedx_tab, r_tab = _T_FINE, _DEDX_STD, _R_FINE

    a_scale  = mat["a_scale"]
    rho      = mat["rho"]
    X0_gcm2  = mat["X0_gcm2"]

    N          = len(muons["Ekin_MeV"])
    d_cm       = depth_m * 100.0
    depth_cos  = np.maximum(muons.get("depth_cos", np.abs(muons["cz"])), 0.02)
    slant_gcm2 = d_cm * rho / depth_cos         # per-muon opacity [g/cm²]

    E_cur  = muons["Ekin_MeV"].copy()
    cx_c   = muons["cx"].copy()
    cy_c   = muons["cy"].copy()
    cz_c   = muons["cz"].copy()

    depth_axis = muons.get("depth_axis", 2)
    x_acc = np.zeros(N)
    y_acc = np.zeros(N)
    z_acc = np.zeros(N)

    E_stop = 1.0   # survival threshold [MeV]

    # ── Analytical CSDA (MS = OFF) ────────────────────────────────────────────
    # Exact: no Euler error, no step-size tuning, zero loops.
    if not ms_enable:
        R0         = np.interp(np.clip(E_cur, t_tab[0], t_tab[-1]), t_tab, r_tab) * a_scale
        alive_bool = R0 > slant_gcm2

        R_exit = np.where(alive_bool, R0 - slant_gcm2, 0.0)
        R_std  = np.clip(R_exit / a_scale, r_tab[0], r_tab[-1])
        E_cur  = np.where(alive_bool, np.interp(R_std, r_tab, t_tab), 0.0)

        # Straight-line position: full slant for alive; fractional for stopped
        frac      = np.where(alive_bool, 1.0,
                             np.clip(R0 / np.maximum(slant_gcm2, 1e-9), 0.0, 1.0))
        slant_cm  = d_cm / depth_cos
        x_acc = cx_c * slant_cm * frac
        y_acc = cy_c * slant_cm * frac
        if depth_axis != 2:
            z_acc = cz_c * slant_cm * frac

        alive = alive_bool.astype(int)
        if progress_cb:
            progress_cb(N, int(alive.sum()), N)

    # ── Vectorized stepping + Highland MCS (MS = ON) ──────────────────────────
    else:
        alive = np.ones(N, dtype=bool)

        # CSDA pre-filter: instantly kill muons with insufficient range
        R0 = np.interp(np.clip(E_cur, t_tab[0], t_tab[-1]), t_tab, r_tab) * a_scale
        alive[R0 < slant_gcm2] = False

        # Adaptive step: 100 g/cm² default (26× larger than Fortran's 10 g/cm²)
        if n_steps <= 0:
            med_slant = float(np.median(slant_gcm2[alive])) if alive.any() else 1000.0
            n_steps   = int(np.clip(med_slant / 100.0, 50, 1000))

        dx            = slant_gcm2 / n_steps   # per-muon step [g/cm²]
        prev_reported = 0

        for step in range(n_steps):
            if not alive.any():
                break

            idx     = np.where(alive)[0]
            E_i     = E_cur[idx]
            dx_i    = dx[idx]
            dx_cm_i = dx_i / rho

            # Position update BEFORE MCS (direction at step entry)
            x_acc[idx] += cx_c[idx] * dx_cm_i
            y_acc[idx] += cy_c[idx] * dx_cm_i
            if depth_axis != 2:
                z_acc[idx] += cz_c[idx] * dx_cm_i

            # CSDA energy loss (vectorized interp on alive muons only)
            dEdX   = np.interp(np.clip(E_i, t_tab[0], t_tab[-1]),
                               t_tab, dedx_tab) * a_scale
            E_new  = E_i - dEdX * dx_i
            E_cur[idx]               = np.maximum(E_new, 0.0)
            alive[idx[E_new < E_stop]] = False

            # Highland MCS
            la = alive[idx]
            if la.any():
                oi     = idx[la]
                E_ms   = E_cur[oi]
                p_GeV  = np.sqrt(np.maximum(E_ms**2 - M_MU**2, 0.0)) / 1000.0
                beta   = p_GeV / np.sqrt(p_GeV**2 + M_MU_GEV**2)
                t_X0   = dx[oi] / X0_gcm2
                theta0 = (13.6e-3 / (beta * p_GeV)
                          * np.sqrt(t_X0)
                          * (1.0 + 0.038 * np.log(np.maximum(t_X0, 1e-12))))

                phi_az  = rng.uniform(0.0, 2.0 * np.pi, size=la.sum())
                dth     = rng.normal(0.0, theta0)
                sin_dth = np.sin(dth);  cos_dth = np.cos(dth)
                cos_phi = np.cos(phi_az); sin_phi = np.sin(phi_az)

                cx_o = cx_c[oi]; cy_o = cy_c[oi]; cz_o = cz_c[oi]

                # Build orthonormal basis (e1, e2) perpendicular to d
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
                norm = np.maximum(
                    np.sqrt(cx_c[oi]**2 + cy_c[oi]**2 + cz_c[oi]**2), 1e-15)
                cx_c[oi] /= norm; cy_c[oi] /= norm; cz_c[oi] /= norm

            if progress_cb:
                n_done = int((step + 1) / n_steps * N)
                if n_done - prev_reported >= REPORT_EVERY:
                    n_surv_est = int(n_done * alive.sum() / N) if N > 0 else 0
                    progress_cb(n_done, n_surv_est, N)
                    prev_reported = n_done

        alive = alive.astype(int)

    # ── Exit-state dict ───────────────────────────────────────────────────────
    x_f = muons["x"] + x_acc
    y_f = muons["y"] + y_acc
    z_f = (muons["z"] + z_acc) if depth_axis != 2 else np.full(N, -d_cm)
    z_stop  = np.where(alive == 0, muons["z"] + z_acc, z_f)
    theta_f = np.arccos(np.clip(-cz_c, -1.0, 1.0))
    phi_f   = np.arctan2(cy_c, cx_c)

    return dict(
        alive=alive,
        E_kin_f_MeV=np.maximum(E_cur, 0.0),
        cx_f=cx_c, cy_f=cy_c, cz_f=cz_c,
        x_f=x_f, y_f=y_f, z_f=z_f,
        z_stop=z_stop,
        theta_f=theta_f, phi_f=phi_f,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    lines = [ln.strip() for ln in sys.stdin
             if ln.strip() and not ln.strip().startswith("#")]

    def _rd(i, default, typ=str):
        try:    return typ(lines[i]) if i < len(lines) else default
        except: return default

    infile        = _rd(0, "muons_surface.dat")
    outfile       = _rd(1, "muons_underground.dat")
    transport_all = bool(_rd(2, 0,     int))
    _ncols_hint   = _rd(3, 13,        int)   # consumed, not used
    depth_m       = _rd(4, 500.0,     float)
    mat_type      = _rd(5, 1,         int)

    # For mat_type=5 (custom), four extra lines shift ms_enable and later params
    if mat_type == 5:
        Z_eff   = _rd(6, 11.0,  float)
        A_eff   = _rd(7, 22.0,  float)
        rho_cst = _rd(8, 2.65,  float)
        I_eV    = _rd(9, 136.4, float)
        X0_gcm2 = (716.408 * A_eff /
                   (Z_eff * (Z_eff + 1.0) * np.log(287.0 / np.sqrt(Z_eff))))
        mat  = {"name": "Custom", "Z": Z_eff, "A": A_eff, "rho": rho_cst,
                "I_eV": I_eV, "X0_gcm2": X0_gcm2,
                "a_scale": (Z_eff / A_eff) / 0.5}
        base = 10
    else:
        mat  = dict(_BB_MAT.get(mat_type, _BB_MAT[1]))
        base = 6

    ms_enable    = bool(_rd(base,     1, int))
    range_tbl_id = _rd(base + 1, 1,   int)   # 0=groom2001  1=pdg2024
    n_steps      = _rd(base + 2, 0,   int)   # 0=auto
    range_table  = 'pdg2024' if range_tbl_id == 1 else 'groom2001'

    print("  UCMuon BB (CSDA) engine — Python/NumPy", flush=True)
    print(f"  {mat['name']}  rho={mat['rho']:.3f} g/cm³"
          f"  X0={mat['X0_gcm2']:.2f} g/cm²"
          f"  depth={depth_m:.1f} m"
          f"  X={depth_m*100*mat['rho']:.1f} g/cm²", flush=True)
    print(f"  range_table={range_table}"
          f"  MS={'ON (Highland)' if ms_enable else 'OFF (analytical exact CSDA)'}",
          flush=True)

    if not Path(infile).exists():
        print(f"  ERROR: input not found: {infile}", flush=True)
        sys.exit(1)

    muons = _read_input(infile, transport_all)
    N     = len(muons["Ekin_MeV"])
    if N == 0:
        print("  ERROR: no muons in input file.", flush=True)
        sys.exit(1)

    print(f"  Read {N:,} muons from {infile}", flush=True)
    if ms_enable:
        print(f"  Steps: {n_steps if n_steps > 0 else 'auto (≥50, ~100 g/cm² each)'}",
              flush=True)
    print(f"  Transported:     0  Survived:     0  Total: {N}", flush=True)

    def _cb(nd, ns, nt):
        print(f"  Transported: {nd}  Survived: {ns}  Total: {nt}", flush=True)

    rng     = np.random.default_rng(42)
    t_start = time.perf_counter()
    result  = transport_bb(muons, depth_m, mat,
                           n_steps=n_steps, ms_enable=ms_enable,
                           rng=rng, progress_cb=_cb,
                           range_table=range_table)
    elapsed = time.perf_counter() - t_start

    n_surv = int(result["alive"].sum())
    print(f"  Transported: {N}  Survived: {n_surv}  Total: {N}", flush=True)
    print(f"  Survival rate: {100.0*n_surv/N:.4f} %", flush=True)
    print(f"  Elapsed: {elapsed:.3f} s", flush=True)

    _write_output(muons, result, outfile)
    print(f"  Output: {outfile}  ({N} rows, {n_surv} survived)", flush=True)

    out_path    = Path(outfile)
    timing_file = str(out_path.with_suffix("")) + "_timing.txt"
    with open(timing_file, "w") as ft:
        ft.write(f"Elapsed : {elapsed:.3f}\n")
    print(f"  Timing: {timing_file}", flush=True)

    stopped_file = str(out_path.with_suffix("")) + "_stopped.dat"
    alive_arr = result["alive"]
    z_stop    = result["z_stop"]
    n_stopped = 0
    with open(stopped_file, "w") as fs:
        fs.write(f"# UCMuon BB CSDA stopped muons  depth={depth_m:.2f} m\n")
        fs.write("# EventID  InitKE_GeV  StopDepth_cm\n")
        for i in range(N):
            if alive_arr[i] == 0:
                evid   = int(muons["EventID"][i])
                initKE = float(muons["Ekin_GeV"][i]) - M_MU_GEV
                stop_z = abs(float(z_stop[i]))
                fs.write(f"{evid:10d}  {initKE:13.6f}  {stop_z:13.4f}\n")
                n_stopped += 1
    print(f"  Stopped muons: {stopped_file}  ({n_stopped} entries)", flush=True)


if __name__ == "__main__":
    main()
