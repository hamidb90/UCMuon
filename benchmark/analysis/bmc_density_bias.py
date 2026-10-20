#!/usr/bin/env python3
"""Density bias of the Vesuvius example caused by the backward CSDA flux
integrator (Engine 5), propagated through the example's own inversion, for the
integrator without (mode 1) and with (mode 2) its forward calibration.

The shipped library T_2(rho) comes from mode 2; the mode-1 library is
T_1 = T_2 / C, C(X, theta) the calibration (exact: C is 1 at X = 0). The
forward engines give T_f = T_m / r_m, with r_m(X, theta) the measured ratio of
mode m to the forward engines (mean of UCMuon-MC, MUSIC and PROPOSAL):

    mode 1  benchmark/results/backward_validation_20261019.csv, at the
            calibration's nodes (T_backward / T_forward)
    mode 2  benchmark/results/backward_validation_offnode_20261020.csv, at
            zenith angles and opacities that are not nodes
            (T_backward_calibrated / T_forward)

interpolated in theta and ln X at each pixel's slant opacity
X(rho) = X(2.65) * rho / 2.65 and zenith angle. The script draws the example's
synthetic data (rho_true = 2.0, 10^5 events, seed 42) from T_f instead of T_m,
inverts it with T_m, and compares with the closure test (the same data drawn
from T_m).

    python3 benchmark/analysis/bmc_density_bias.py
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "gui"))
import ucmuon_backward_mc as bm  # noqa: E402
from ucmuon_density_analysis import (build_tsim_library, generate_synthetic_tdata,  # noqa: E402
                                     invert_density_map, load_transmission_map)

LIB = REPO / "examples" / "vesuvius" / "tsim_library"
RESULTS = REPO / "benchmark" / "results"
VALIDATION = {1: (RESULTS / "backward_validation_20261019.csv", "T_backward"),
              2: (RESULTS / "backward_validation_offnode_20261020.csv", "T_backward_calibrated")}
SPEC = 3                                       # Guan, as the library
RHOS = ("1.50", "2.00", "2.65", "3.00")
RHO_TRUE, N_EVENTS, SEED = 2.0, 100_000, 42


def ratio_table(path, column):
    """(thetas, ln X grid, r[theta, X]) from a validation, mean over engines."""
    acc = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            acc[(float(row["theta_deg"]), float(row["X_gcm2"]))].append(
                float(row[column]) / float(row["T_forward"]))
    th = sorted({k[0] for k in acc}); xs = sorted({k[1] for k in acc})
    r = np.array([[np.mean(acc[(a, x)]) for x in xs] for a in th])
    return np.array(th), np.log(np.array(xs)), r


def r_at(theta_deg, X, th, lnx, r):
    """Interpolate r; below the first X it goes linearly to 1 at X = 0, above
    the last it is held (counted separately)."""
    rows = np.array([np.interp(np.log(max(X, np.exp(lnx[0]))), lnx, rr) for rr in r])
    val = float(np.interp(theta_deg, th, rows))
    if X < np.exp(lnx[0]):
        val = 1.0 + (val - 1.0) * X / np.exp(lnx[0])
    return val


def overburden_265(az_c, el_c):
    """Overburden at rho = 2.65 [g/cm^2] on the library's (az, el) grid."""
    data = np.loadtxt(LIB / "terrain_overburden_ref2.65.dat", comments="#")
    ob = {}
    for az, ze, x, _ in data:
        ob[(round(az, 2), round(90.0 - ze, 2))] = x
    return np.array([[ob.get((round(a, 2), round(e, 2)), np.nan) for e in el_c] for a in az_c])


def pixel_map(T, x265, rho, el_c, fn):
    """fn(theta_deg, X) on every rock pixel of a library map, 1 elsewhere."""
    out = np.ones_like(T)
    for ia in range(T.shape[0]):
        for ie in range(T.shape[1]):
            X = x265[ia, ie] * rho / 2.65
            if np.isfinite(X) and X > 0:
                out[ia, ie] = fn(90.0 - el_c[ie], X)
    return out


def bias(lib, lib_f, x265, r_true, label):
    out = {}
    for name, truth in (("closure", lib), ("forward truth", lib_f)):
        T_data, sigma_T = generate_synthetic_tdata(truth, RHO_TRUE, n_events=N_EVENTS, seed=SEED)
        rho_map, sigma_rho, status = invert_density_map(T_data, lib, sigma_T)
        ok = status == 0
        out[name] = (rho_map, ok, sigma_rho)
        print(f"{label} {name:14s}: {ok.sum()} OK pixels, median rho {np.nanmedian(rho_map[ok]):.3f} "
              f"g/cm3, scatter {np.nanstd(rho_map[ok]):.3f}, median sigma_rho "
              f"{np.nanmedian(sigma_rho[ok]):.3f}")
    (rc, okc, _), (rf, okf, _) = out["closure"], out["forward truth"]
    both = okc & okf
    d = rf[both] - rc[both]
    X2 = (x265 * RHO_TRUE / 2.65)[both]
    print(f"{label} pixels OK in both: {both.sum()}; rho(forward truth) - rho(closure): median "
          f"{np.median(d):+.3f}, 16-84% {np.percentile(d, 16):+.3f} to {np.percentile(d, 84):+.3f} g/cm3")
    print(f"{label} slant opacity of those pixels at rho = 2.0: median {np.median(X2):.3g} g/cm2, "
          f"range {X2.min():.3g} to {X2.max():.3g}; T_{label[-1]} / T_forward there: median "
          f"{np.median(r_true[both]):.3f}, range {r_true[both].min():.3f} to {r_true[both].max():.3f}")


def main() -> int:
    lib2 = build_tsim_library([str(LIB / f"terrain_transmission_{r}.dat") for r in RHOS])
    az_c, el_c, _, _ = load_transmission_map(str(LIB / "terrain_transmission_2.65.dat"))
    x265 = overburden_265(az_c, el_c)
    cal = {rho: pixel_map(T, x265, rho, el_c,
                          lambda th, X: bm.forward_calibration(X, np.radians(th), SPEC))
           for rho, T in lib2.items()}
    libs = {2: lib2,
            1: {rho: np.where(np.isfinite(T), T / cal[rho], T) for rho, T in lib2.items()}}
    for mode in (1, 2):
        th, lnx, r = ratio_table(*VALIDATION[mode])
        lib = libs[mode]
        lib_f, rmaps = {}, {}
        for rho, T in lib.items():
            rmaps[rho] = pixel_map(T, x265, rho, el_c, lambda a, X: r_at(a, X, th, lnx, r))
            lib_f[rho] = np.where(np.isfinite(T), np.minimum(T / rmaps[rho], 1.0), T)
        X2 = x265 * RHO_TRUE / 2.65
        print(f"mode {mode}: r from {VALIDATION[mode][0].name}, theta {th.min():g}-{th.max():g} deg, "
              f"X {np.exp(lnx[0]):.3g}-{np.exp(lnx[-1]):.3g} g/cm2; rock pixels at rho = 2.0 beyond "
              f"its last X: {int(np.sum(np.isfinite(X2) & (X2 > np.exp(lnx[-1]))))}")
        bias(lib, lib_f, x265, rmaps[RHO_TRUE], f"mode {mode}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
