#!/usr/bin/env python3
"""The same spectrum must give the same flux in every implementation.

UCMuon carries each surface spectrum three times: in the Fortran generator
(src/generator/ucmuon_source_module.f90), in UCMuGen (ucmugen/include/UCMuGen.h)
and in Python (gui/fast_flux_estimator.py, which the backward MC and the
guaranteed-hit generator also use). Before the 2026-09 flux audit the three
disagreed: on the grid below, Reyna differed between UCMuGen and Python by
factors of 0.45 to 1.4e4, and the Gaisser-formula models in Python were 10x
low. This test is what keeps them together.

Compared, for spectra 1 (CosmoALEPH), 4 (Guan), 5 (Frosin), 6 (Gaisser) and
7 (Reyna):
  I   differential intensity I(p, cos theta) on an 8 x 5 grid   rel. tol 1e-6
  V   vertical integrated flux, E_min..1500 GeV                  rel. tol 1 %
      (Fortran: the "Integrated flux" line build_cosmoaleph_cdf prints,
      which is what the GUI rate estimate reads)
  Q   rate through a horizontal surface, theta <= 85 deg,
      by the same quadrature in all three languages             rel. tol 1e-6
  M   UCMuGen Generator::rateAndError against that quadrature    rel. tol 1 %
      (and within 4 sigma of its own Monte Carlo error)
  W   UCMuGen Generator::warnings() fires exactly when
      fast_flux_estimator.validity_warning does

    python3 tests/flux/test_flux_consistency.py        # exit 0 = all agree

Needs c++ (C++17) and gfortran. Exit 1 on any disagreement.
"""
from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
import fast_flux_estimator as ffe  # noqa: E402

M_MU = 0.10566          # GeV, the value all three implementations use
PI = 3.141592654        # kPi / PI in UCMuGen and the Fortran
EMAX, THETA_MAX_DEG = 1500.0, 85.0
PY_MODEL = {1: ffe._cosmoaleph, 4: ffe._guan_2015, 5: ffe._frosin_2025,
            6: ffe._bugaev, 7: ffe._reyna_bugaev}
NAME = {1: "CosmoALEPH", 4: "Guan", 5: "Frosin", 6: "Gaisser", 7: "Reyna"}

failures: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {what}")
    if not ok:
        failures.append(what)


def py_intensity(mode: int, p: float, c: float) -> float:
    """dN/dp [cm^-2 s^-1 sr^-1 (GeV/c)^-1] from the Python dΦ/dT."""
    E = math.sqrt(p * p + M_MU * M_MU)
    T = np.array([E - M_MU])
    theta = math.degrees(math.acos(min(1.0, c)))
    return float(PY_MODEL[mode](T, theta)[0]) * p / E


def py_horizontal_rate(mode: int, p_min: float, p_max: float) -> float:
    """Same quadrature as horizontalRate() in intensity_dump.cc."""
    np_, nc = 2000, 400
    c_lo = math.cos(THETA_MAX_DEG * PI / 180.0)
    p = p_min * np.exp(math.log(p_max / p_min) * np.arange(np_) / (np_ - 1))
    E = np.sqrt(p * p + M_MU * M_MU)
    total = 0.0
    for k in range(nc):
        c = c_lo + (k + 0.5) * (1.0 - c_lo) / nc
        theta = math.degrees(math.acos(c))
        f = PY_MODEL[mode](E - M_MU, theta) * p / E
        total += c * float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(p))) * (1.0 - c_lo) / nc
    return 2.0 * PI * total


def py_vertical(mode: int, p_min: float, p_max: float) -> float:
    p = np.logspace(math.log10(p_min), math.log10(p_max), 20000)
    E = np.sqrt(p * p + M_MU * M_MU)
    f = PY_MODEL[mode](E - M_MU, 0.0) * p / E
    return float(np.trapezoid(f, p)) if hasattr(np, "trapezoid") else float(np.trapz(f, p))


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def build_and_run(tmp: Path) -> tuple[str, str]:
    cxx = tmp / "dump_cc"
    f90 = tmp / "dump_f90"
    subprocess.run(["c++", "-std=c++17", "-O2", f"-I{ROOT/'ucmugen/include'}",
                    str(ROOT / "tests/flux/intensity_dump.cc"), "-o", str(cxx)],
                   check=True)
    subprocess.run(["gfortran", "-O2", "-fopenmp", f"-J{tmp}",
                    str(ROOT / "src/generator/rng_parallel.f90"),
                    str(ROOT / "src/generator/ucmuon_source_module.f90"),
                    str(ROOT / "tests/flux/intensity_dump.f90"), "-o", str(f90)],
                   check=True, cwd=tmp)
    return (subprocess.run([str(cxx)], check=True, capture_output=True, text=True).stdout,
            subprocess.run([str(f90)], check=True, capture_output=True, text=True).stdout)


def parse_cc(out: str):
    I, V, Q, M, W = {}, {}, {}, {}, {}
    for line in out.splitlines():
        t = line.split()
        if t[0] == "I":
            I[(int(t[1]), float(t[2]), float(t[3]))] = float(t[4])
        elif t[0] == "V":
            V[(int(t[1]), float(t[2]))] = float(t[3])
        elif t[0] == "Q":
            Q[(int(t[1]), float(t[2]))] = float(t[3])
        elif t[0] == "M":
            M[(int(t[1]), float(t[2]))] = (float(t[3]), float(t[4]))
        elif t[0] == "W":
            W[(int(t[1]), float(t[2]))] = int(t[3])
    return I, V, Q, M, W


def parse_f90(out: str):
    I, V, Q = {}, {}, {}
    key = None
    for line in out.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "I":
            I[(int(t[1]), float(t[2]), float(t[3]))] = float(t[4])
        elif t[0] == "DUMPCDF":
            key = (int(t[1]), float(t[2]))
        elif "Integrated flux" in line and key is not None:
            V[key] = float(line.split(":")[1].split()[0])
        elif t[0] == "Q":
            Q[(int(t[1]), float(t[2]))] = float(t[3])
    return I, V, Q


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        out_cc, out_f90 = build_and_run(Path(d))
    Icc, Vcc, Qcc, Mcc, Wcc = parse_cc(out_cc)
    Iff, Vff, Qff = parse_f90(out_f90)

    print("I(p, cos theta): C++ vs Fortran vs Python, 200 grid points")
    for mode in NAME:
        keys = [k for k in Icc if k[0] == mode]
        worst_f = max(rel(Icc[k], Iff[k]) for k in keys)
        worst_p = max(rel(Icc[k], py_intensity(mode, k[1], k[2])) for k in keys)
        check(worst_f < 1e-6 and worst_p < 1e-6,
              f"{NAME[mode]:10s} max rel diff  C++/F90 {worst_f:.1e}  C++/Py {worst_p:.1e}")

    print("\nVertical integrated flux, E_min..1500 GeV [cm^-2 s^-1 sr^-1]")
    for (mode, emin), vcc in sorted(Vcc.items()):
        p_min, p_max = math.sqrt(emin**2 - M_MU**2), math.sqrt(EMAX**2 - M_MU**2)
        vpy, vff = py_vertical(mode, p_min, p_max), Vff[(mode, emin)]
        check(rel(vcc, vff) < 0.01 and rel(vcc, vpy) < 0.01,
              f"{NAME[mode]:10s} E>{emin:5.0f}  C++ {vcc:.4e}  F90 {vff:.4e}  Py {vpy:.4e}")

    print(f"\nHorizontal-surface rate, theta <= {THETA_MAX_DEG:g} deg [cm^-2 s^-1]")
    for (mode, emin), qcc in sorted(Qcc.items()):
        p_min, p_max = math.sqrt(emin**2 - M_MU**2), math.sqrt(EMAX**2 - M_MU**2)
        qpy, qff = py_horizontal_rate(mode, p_min, p_max), Qff[(mode, emin)]
        r_mc, e_mc = Mcc[(mode, emin)]
        check(rel(qcc, qff) < 1e-6 and rel(qcc, qpy) < 1e-6,
              f"{NAME[mode]:10s} E>{emin:5.0f}  quadrature C++ {qcc:.5e}  F90 {qff:.5e}  Py {qpy:.5e}")
        check(rel(r_mc, qcc) < 0.01 and abs(r_mc - qcc) < 4.0 * e_mc + 1e-4 * qcc,
              f"{NAME[mode]:10s} E>{emin:5.0f}  Generator::rateAndError {r_mc:.5e} ± {e_mc:.1e}"
              f"  vs quadrature ({100*(r_mc/qcc-1):+.2f} %)")

    print("\nValidity warnings: UCMuGen vs fast_flux_estimator")
    py_key = {1: None, 4: "guan_2015", 5: "frosin_2025", 6: "bugaev", 7: "reyna_bugaev"}
    for (mode, emin), nw in sorted(Wcc.items()):
        if mode == 1:        # CosmoALEPH is not a fast-estimator model
            expect = emin < 99.0
        else:
            expect = ffe.validity_warning(py_key[mode], emin - M_MU) is not None
        check((nw > 0) == expect,
              f"{NAME[mode]:10s} E>{emin:5.0f}  UCMuGen warns: {nw > 0}  expected: {expect}")

    print()
    if failures:
        print(f"{len(failures)} FAILED")
        return 1
    print("ALL AGREE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
