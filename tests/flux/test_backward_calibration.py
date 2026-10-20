#!/usr/bin/env python3
"""The backward MC's forward calibration (mode 2, the default).

gui/ucmuon_backward_mc.py neglects energy-loss fluctuations, so its
transmission is up to 15 % below the forward engines' at 3000 m w.e. Mode 2
multiplies the mode-1 flux by C(spectrum, θ, X), measured with UCMuon-MC and
MUSIC (benchmark/results/backward_calibration_20261020.csv). This checks:

  - the table in the module is the one in the results file, every node;
  - mode 2 is mode 1 times C, exactly C at a node, and the open sky (X = 0) is
    unchanged, so T = Φ(X)/Φ(0) is mode 1's T times C;
  - C goes linearly to 1 below the first node and is held beyond the last,
    with one warning;
  - other materials get mode 1, with one warning;
  - mode 2 is the default of directional_flux, backward_mc_flux, the terrain
    driver's compute_flux_map and its command line.

    python3 tests/flux/test_backward_calibration.py      # exit 0 = all pass
"""
from __future__ import annotations

import csv
import inspect
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
import ucmuon_backward_mc as bmc  # noqa: E402
import ucmuon_terrain_driver as td  # noqa: E402

failures: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {what}")
    if not ok:
        failures.append(what)


def flux(X, theta_deg, mode, mat_id=1):
    return bmc.directional_flux(X, math.radians(theta_deg), 3, E_min_GeV=1.0, E_max_GeV=2500.0,
                                n_E=40, mode=mode, mat_id=mat_id)


def main() -> int:
    print("table")
    res = ROOT / "benchmark" / "results" / "backward_calibration_20261020.csv"
    with open(res) as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    bad = [r for r in rows
           if abs(bmc._CAL[int(r["spectrum_mode"])][
               list(bmc._CAL_THETA_DEG).index(float(r["theta_deg"])),
               list(bmc._CAL_X).index(float(r["X_gcm2"]))] - float(r["C"])) > 5e-5]
    check(len(rows) == 4 * len(bmc._CAL_THETA_DEG) * len(bmc._CAL_X) and not bad,
          f"{len(rows)} nodes in {res.name}, all equal to the module's table")
    check(all(bmc.has_calibration(s) for s in (1, 2, 3, 4)) and not bmc.has_calibration(5),
          "spectra 1-4 calibrated")

    print("\nmode 2 = mode 1 x C (Guan)")
    for th, X in ((60.0, 2.0e5), (70.0, 1.5e5), (82.0, 3.5e4)):
        c = bmc.forward_calibration(X, math.radians(th), 3)
        r = flux(X, th, 2) / flux(X, th, 1)
        check(abs(r / c - 1.0) < 1e-12, f"θ {th:g}°, {X / 100:.0f} m w.e.: ratio {r:.4f} = C")
    node = bmc._CAL[3][1, list(bmc._CAL_X).index(2.0e5)]
    check(abs(bmc.forward_calibration(2.0e5, math.radians(60.0), 3) - node) < 1e-12,
          f"at the node θ 60°, 2000 m w.e.: C = {node:.4f}")
    check(flux(0.0, 60.0, 2) == flux(0.0, 60.0, 1), "open sky unchanged")
    c1 = bmc.forward_calibration(1.0e4, 0.0, 3)
    check(abs(bmc.forward_calibration(5.0e3, 0.0, 3) - (1.0 + 0.5 * (c1 - 1.0))) < 1e-12,
          "half the first opacity: half of C - 1")

    print("\nwarnings, once per process")
    bmc._CAL_WARNED.clear()
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        c_last = bmc.forward_calibration(5.0e5, 0.0, 3)
        held = [bmc.forward_calibration(X, 0.0, 3) for X in (6.0e5, 7.0e5)]
        water = [flux(1.0e5, 0.0, m, mat_id=2) for m in (1, 2)]
        flux(1.0e5, 30.0, 2, mat_id=4)
    msgs = [str(w.message) for w in rec if "mode 2" in str(w.message)]
    check(held == [c_last, c_last] and sum("beyond" in m for m in msgs) == 1,
          f"beyond 5000 m w.e.: C held at {c_last:.4f}, one warning")
    check(water[0] == water[1] and sum("Standard Rock only" in m for m in msgs) == 1,
          "water and iron: mode 2 = mode 1, one warning")

    print("\ndefaults")
    dflt = {f.__name__: inspect.signature(f).parameters["mode"].default
            for f in (bmc.directional_flux, bmc.backward_mc_flux, td.compute_flux_map)}
    check(set(dflt.values()) == {2}, f"mode defaults: {dflt}")
    for path, line in (("gui/ucmuon_backward_mc.py", "_rd(10, 2,"),
                       ("gui/ucmuon_terrain_driver.py", "_rd(10,  2,")):
        check(line in (ROOT / path).read_text(), f"{path}: stdin mode defaults to 2")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        info = bmc.backward_mc_flux(depth_m=500.0, rho=2.65, mat_id=1, spectrum_mode=3,
                                    E_min_GeV=1.0, E_max_GeV=1000.0, theta_max_deg=30.0,
                                    n_E=20, n_theta=4)["info"]
    check("forward-calibrated" in info["mode"], f"backward_mc_flux info: {info['mode']}")

    print()
    if failures:
        print(f"{len(failures)} FAILED")
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
