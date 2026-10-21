#!/usr/bin/env python3
"""Tests of the Vesuvius density inversion beyond the closure test.

The closure test of MURAVES_GUIDE.md draws synthetic data from the shipped
library at a library density (rho_true = 2.0) and inverts them with the same
library. This script measures what that test cannot show:

  1. pulls (rho_hat - rho_true) / sigma_rho, at and between library densities,
     for 10^5 and 10^4 open-sky events per pixel (20 seeds each);
  2. the interpolation bias between library densities: noise-free data
     computed by Engine 5 at the true density itself, inverted with the
     four-density library, interpolating ln T (the code) and T (up to 1.3.2);
  3. data from the forward engines instead of Engine 5: T_f = T / r, r the
     measured ratio of Engine 5 (mode 2) to UCMuon-MC, MUSIC and PROPOSAL
     (benchmark/results/backward_validation_offnode_20261020.csv);
  4. data from another flux model: Engine 5 with Frosin 2025 instead of Guan,
     inverted with the Guan library;
  5. Method 2 (direct opacity inversion with the analytic flux models of
     gui/fast_flux_estimator.py) on the closure-test data, for every model,
     at sea level and at the detector's 608 m;
  6. the exposure that n_events = 10^5 open-sky muons per 1 x 1 deg pixel
     implies for an ideal detector;
  7. the altitude: the library uses a sea-level spectrum, the detector is at
     608 m. PARMA at the site (parma_altitude.cc; skipped without PARMA or a
     C++ compiler) gives T at 608 m and at sea level in the CSDA limit, and
     the density a sea-level library would read from 608 m data.

    python3 examples/vesuvius/inversion_tests.py        # from the repo root
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "gui"))
sys.path.insert(0, str(REPO / "benchmark" / "analysis"))
import ucmuon_backward_mc as bm                       # noqa: E402
import ucmuon_density_analysis as da                  # noqa: E402
from bmc_density_bias import ratio_table, r_at        # noqa: E402

LIB = REPO / "examples" / "vesuvius" / "tsim_library"
RHOS = ("1.50", "2.00", "2.65", "3.00")
KW = dict(E_min_GeV=1.0, E_max_GeV=2500.0, n_E=40, mode=2)   # as make_tsim_library.py
N_EV, SEEDS = 100_000, range(20)
YEAR = 3.156e7


def load():
    lib = da.build_tsim_library([str(LIB / f"terrain_transmission_{r}.dat") for r in RHOS])
    az_c, el_c, _, _ = da.load_transmission_map(str(LIB / "terrain_transmission_2.65.dat"))
    d = np.loadtxt(LIB / "terrain_overburden_ref2.65.dat", comments="#")
    ob = {(round(a, 2), round(90.0 - z, 2)): x for a, z, x, _ in d}
    x265 = np.array([[ob[(round(a, 2), round(e, 2))] for e in el_c] for a in az_c])
    return lib, az_c, el_c, x265


def draw(T_true, n, seed):
    """Synthetic data as generate_synthetic_tdata: N_rock ~ Poisson(n T)."""
    rng = np.random.default_rng(seed)
    T = np.where(np.isfinite(T_true), T_true, 0.0)
    sky = T >= 0.999
    data = np.where(sky, 1.0, rng.poisson(n * np.clip(T, 0.0, 1.0)) / n)
    sig = np.where(sky, 1.0 / np.sqrt(n), np.sqrt(np.clip(T, 0.0, 1.0) / n))
    bad = ~np.isfinite(T_true)
    return np.where(bad, np.nan, data), np.where(bad, np.nan, sig)


def engine5_T(x265, el_c, rho, mask, spec=3):
    """Engine 5 transmission at density rho on the pixels of mask (NaN elsewhere)."""
    T = np.full(x265.shape, np.nan)
    sky = {}
    with bm.collect_range_warnings():      # under 5000 m w.e. mode 2 includes the cut
        for ia, ie in zip(*np.nonzero(mask)):
            th = np.radians(90.0 - el_c[ie])
            if ie not in sky:
                sky[ie] = bm.directional_flux(0.0, th, spec, **KW)
            T[ia, ie] = bm.directional_flux(x265[ia, ie] * rho / 2.65, th, spec, **KW) / sky[ie]
    return T


def invert_linear_T(T_data, lib):
    """The inversion of 1.3.2 and earlier: T interpolated linearly."""
    rv = np.array(sorted(lib)); Ts = np.stack([lib[r] for r in rv])
    out = np.full(T_data.shape, np.nan)
    for ia, ie in zip(*np.nonzero(np.isfinite(T_data))):
        Tv = Ts[:, ia, ie]
        if np.all(np.isfinite(Tv)) and Tv.min() <= T_data[ia, ie] <= Tv.max():
            out[ia, ie] = np.interp(T_data[ia, ie], Tv[::-1], rv[::-1])
    return out


def pulls(lib, T_true, rho_true, n):
    p, n_ok = [], []
    for seed in SEEDS:
        T, s = draw(T_true, n, seed)
        rho, sr, st = da.invert_density_map(T, lib, s)
        ok = st == 0
        n_ok.append(ok.sum())
        p.append((rho[ok] - rho_true) / sr[ok])
    p = np.concatenate(p)
    return np.mean(n_ok), p.mean(), np.sqrt(np.mean(p ** 2)), np.mean(np.abs(p) < 1.0)


def bias(lib, T_true, rho_true, n=N_EV, seeds=SEEDS):
    d = []
    for seed in seeds:
        T, s = draw(T_true, n, seed)
        rho, _, st = da.invert_density_map(T, lib, s)
        d.append(rho[st == 0] - rho_true)
    d = np.concatenate(d)
    return np.median(d), np.percentile(d, 16), np.percentile(d, 84)


def main() -> int:
    lib, az_c, el_c, x265 = load()
    rv = np.array(sorted(lib)); Ts = np.stack([lib[r] for r in rv])
    rock = np.all(np.isfinite(Ts), axis=0) & (Ts[0] < 0.999)

    # The closure test of the guide, to fix the pixel set used below
    T0, s0 = da.generate_synthetic_tdata(lib, 2.0, n_events=N_EV, seed=42)
    rho0, sr0, st0 = da.invert_density_map(T0, lib, s0)
    ok0 = st0 == 0
    print(f"closure test (rho_true 2.0, 1e5 events, seed 42): {ok0.sum()} OK, "
          f"{(st0 == 4).sum()} low sensitivity, mean {rho0[ok0].mean():.3f}, median "
          f"{np.median(rho0[ok0]):.3f}, scatter {rho0[ok0].std():.3f}, mean sigma {sr0[ok0].mean():.3f}")

    print("\n1. pulls, data from Engine 5 at the true density (20 seeds)")
    truth = {r: engine5_T(x265, el_c, r, rock) for r in (1.75, 2.0, 2.3, 2.8)}
    for r, T_true in truth.items():
        for n in (N_EV, 10_000):
            n_ok, mean, rms, inside = pulls(lib, T_true, r, n)
            print(f"   rho_true {r:4.2f}, {n:6d} events: {n_ok:4.0f} OK, pull mean {mean:+.3f}, "
                  f"RMS {rms:.3f}, |pull| < 1: {100 * inside:.0f} %")

    print("\n2. interpolation bias, noise-free data from Engine 5 at the true density")
    for r in (1.75, 2.25, 2.8):
        T_true = engine5_T(x265, el_c, r, ok0)
        log = da.invert_density_map(T_true, lib, np.full(T_true.shape, 1e-9))[0]
        lin = invert_linear_T(T_true, lib)
        m = ok0 & np.isfinite(log) & np.isfinite(lin)
        print(f"   rho_true {r:4.2f} on the {m.sum()} closure-test pixels: ln T {100 * np.median(log[m] / r - 1):+.2f} % "
              f"(median), T {100 * np.median(lin[m] / r - 1):+.2f} %")

    print("\n3. data from the forward engines (T / r, r from the test between calibration nodes)")
    th, lnx, rr = ratio_table(REPO / "benchmark/results/backward_validation_offnode_20261020.csv",
                              "T_backward_calibrated")
    for r in (2.0, 2.3):
        T_e5 = truth[r]
        rmap = np.ones_like(T_e5)
        for ia, ie in zip(*np.nonzero(rock)):
            rmap[ia, ie] = r_at(90.0 - el_c[ie], x265[ia, ie] * r / 2.65, th, lnx, rr)
        med, lo, hi = bias(lib, np.minimum(T_e5 / rmap, 1.0), r, seeds=range(5))
        c_med, _, _ = bias(lib, T_e5, r, seeds=range(5))
        print(f"   rho_true {r:4.2f}: rho_hat - rho_true median {med:+.3f} (16-84 % {lo:+.3f} to {hi:+.3f}); "
              f"same seeds from Engine 5: {c_med:+.3f} g/cm3")

    print("\n4. data from Frosin 2025 (Engine 5, spectrum 4), inverted with the Guan library")
    for r in (2.0, 2.3):
        med, lo, hi = bias(lib, engine5_T(x265, el_c, r, rock, spec=4), r, seeds=range(5))
        print(f"   rho_true {r:4.2f}: rho_hat - rho_true median {med:+.3f} (16-84 % {lo:+.3f} to {hi:+.3f}) g/cm3")

    print("\n5. Method 2 on the closure-test data (median over its OK pixels)")
    _, _, L, _ = da.load_overburden_as_L(str(LIB / "terrain_overburden_ref2.65.dat"))
    for model, label in da.flux_model_labels().items():
        meds = []
        for alt in (0.0, 608.0):
            with warnings.catch_warnings():    # plain Gaisser: shown on purpose, warned about
                warnings.simplefilter("ignore")
                op, _, st = da.invert_opacity_map(T0, el_c, model=model, altitude_m=alt)
            rho2, _ = da.opacity_to_density(op, L)
            m = ok0 & (st == 0) & np.isfinite(rho2)
            meds.append(np.median(rho2[m]))
        print(f"   {model:13s} ({label}): {meds[0]:.3f} at sea level, {meds[1]:.3f} at 608 m")

    print("\n6. exposure implied by 1e5 open-sky muons per 1 x 1 deg pixel (ideal detector)")
    el = np.broadcast_to(el_c, ok0.shape)[ok0]
    ze_lo, ze_hi = np.radians(90.0 - el - 0.5), np.radians(90.0 - el + 0.5)
    dOmega = np.radians(1.0) * (np.cos(ze_lo) - np.cos(ze_hi))
    phi = np.array([bm.directional_flux(0.0, np.radians(90.0 - e), 3, **KW) for e in el])
    At = N_EV / (phi * dOmega) / YEAR
    print(f"   {ok0.sum()} OK pixels at el {el.min():.1f}-{el.max():.1f} deg: area x time "
          f"{np.median(At):.1f} m2 yr (median), {At.min():.1f} to {At.max():.1f}")

    print("\n7. altitude: PARMA at the site, 608 m against sea level (CSDA)")
    altitude_effect()
    return 0


def altitude_effect():
    hdr = REPO / "ucmugen" / "include" / "UCMuGen_PARMA.h"
    data = REPO / "data" / "EXPACS" / "parma"
    cxx = shutil.which("c++")
    if not (hdr.exists() and data.exists() and cxx):
        print("   skipped: needs PARMA (python3 tools/install_parma.py) and c++")
        return
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / "parma_altitude"
        subprocess.run([cxx, "-O2", "-std=c++17", f"-I{REPO / 'ucmugen' / 'include'}",
                        str(Path(__file__).with_name("parma_altitude.cc")), "-o", str(exe)], check=True)
        d = np.loadtxt(subprocess.run([str(exe), str(data)], check=True, capture_output=True,
                                      text=True).stdout.splitlines())
    m_mu = 0.1056584

    def T_csda(alt, th, X):
        sel = np.isclose(d[:, 0], alt) & np.isclose(d[:, 1], th)
        p, I = d[sel, 2], d[sel, 3]
        e = np.sqrt(p * p + m_mu * m_mu)
        lnT, lnf = np.log(e - m_mu), np.log(I * p / e)          # per GeV of kinetic energy

        def integral(lo, hi):
            g = np.logspace(np.log10(lo), np.log10(hi), 2000)
            return np.trapezoid(np.exp(np.interp(np.log(g), lnT, lnf)), g)
        lo, hi = (bm._backward_energy(E, X)[0] / 1000.0 for E in (1.0e3, 2.5e6))
        return integral(lo, hi) / integral(1.0, 2500.0)

    for th in (75.0, 80.0, 85.0):
        out = []
        for X in (1.0e4, 4.3e4, 1.0e5, 2.0e5):
            t608 = T_csda(0.608, th, X)
            xs = X * np.linspace(0.8, 1.2, 81)                   # rho = 2.0 * xs / X
            lt = np.log([T_csda(0.0, th, x) for x in xs])
            rho_read = 2.0 * np.interp(np.log(t608), lt[::-1], xs[::-1]) / X
            out.append(f"{X / 100:.0f} m w.e. T608/T0 {t608 / T_csda(0.0, th, X):.3f}, reads {rho_read:.3f}")
        print(f"   zenith {th:.0f}: " + "; ".join(out))


if __name__ == "__main__":
    sys.exit(main())
