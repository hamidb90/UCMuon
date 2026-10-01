#!/usr/bin/env python3
"""The GUI live time agrees with UCMuGen.

The GUI computes T = N_tried / R from what the Fortran generator prints
("Tried", "Surface rate R"), through gui/live_time.py, the functions
_compute_flux in gui/ucmuon_gui.py calls. UCMuGen's reference is
liveTime(N) = N / rate(), with rate() the rate through the surface without a
detector and into the detector with one. For each source (plane, disk,
hemisphere), with and without a detector filter:

  mode 6   T_GUI = N_tried / R    vs   T_UCMuGen = N_saved / rate()     1%
  legacy   Σw over the saved rows / T_GUI  vs  UCMuGen rate()          2%
           (angular modes 2 and 4, importance-weighted)
  R        Python quadrature vs the Fortran's printed R              0.1%

The 1% and 2% become 4σ when the run's statistical error is larger than a
quarter of them (σ from the hit count, the weights and UCMuGen's own rate
error): at the --quick size a filtered run has ±0.5%, and a fixed 1% would
fail a correct build one run in eight.

    python3 tests/flux/test_live_time.py [--quick]

Needs bin/ucmuon_gen_omp (make ucmuon_gen_omp) and a C++17 compiler.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _harness as H  # noqa: E402
import live_time as LT  # noqa: E402  (gui/, put on the path by _harness)

QUICK = "--quick" in sys.argv
N_HITS = 40_000 if QUICK else 120_000
failures: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {what}")
    if not ok:
        failures.append(what)


def within(dev: float, tol: float, sigma: float) -> bool:
    return abs(dev) < max(tol, 4.0 * sigma)


def gui_live_time(run: H.Run) -> tuple[float, float]:
    parsed = LT.parse_generator_output(run.stdout.splitlines())
    return LT.live_time(parsed["tried"], parsed["surface_rate"]), parsed["surface_rate"]


CYL = H.Detector("cyl", (0.0, 0.0, -2400.0, 0.0, 0.0, -2000.0, 150.0), margin=50.0)
BOX = H.Detector("box", (-200.0, 200.0, -200.0, 200.0, -2200.0, -2000.0), margin=20.0)

CASES = [
    ("plane",      dict(source_mode=2, half_m=30.0),   BOX),
    ("disk",       dict(source_mode=1, radius_m=30.0), CYL),
    ("hemisphere", dict(source_mode=3, radius_m=30.0), CYL),
]


def main() -> int:
    for name, geom, det in CASES:
        for with_det in (False, True):
            cfg = H.Config(spectrum=4, emin=15.0, emax=1500.0, theta_max_deg=85.0,
                           angular_mode=6, detector=det if with_det else None,
                           n=N_HITS, **geom)
            label = f"{name:10s} {'with' if with_det else 'no  '} filter"
            fr = H.run_fortran(cfg)
            uc = H.run_ucmugen(cfg, n=0)
            T_gui, R_f = gui_live_time(fr)
            T_uc = len(fr.p) / uc.rate
            err = math.hypot(1.0 / np.sqrt(len(fr.p)) if with_det else 0.0,
                             uc.rate_err / uc.rate)
            check(within(T_gui / T_uc - 1.0, 0.01, err),
                  f"{label}: T_GUI {T_gui:.5g} s vs UCMuGen {T_uc:.5g} s "
                  f"({100*(T_gui/T_uc-1):+.2f} %, stat ±{100*err:.2f} %)")
            if not with_det:
                R_py = LT.surface_rate(4, 15.0, 1500.0, 85.0, cfg.source())
                check(abs(R_py / R_f - 1.0) < 1e-3,
                      f"{label}: R Python {R_py:.6g} vs Fortran {R_f:.6g} s^-1")

    # Legacy modes: rows are not flux-distributed, the GUI weights them.
    for amode in (2, 4):
        cfg = H.Config(spectrum=4, emin=15.0, emax=1500.0, theta_max_deg=85.0,
                       angular_mode=amode, detector=CYL, n=N_HITS,
                       source_mode=1, radius_m=30.0)
        fr = H.run_fortran(cfg)
        uc = H.run_ucmugen(H.Config(**{**cfg.__dict__, "angular_mode": 6}), n=0)
        T_gui, R_f = gui_live_time(fr)
        w = LT.event_weights(4, amode, 15.0, 1500.0, 85.0, cfg.source(),
                             fr.p, fr.d, fr.x, R_f)
        rate_gui = w.sum() / T_gui
        err = math.hypot(np.sqrt((w ** 2).sum()) / w.sum(), uc.rate_err / uc.rate)
        check(within(rate_gui / uc.rate - 1.0, 0.02, err),
              f"disk with filter, legacy mode {amode}: weighted rate {rate_gui:.5g} vs "
              f"UCMuGen {uc.rate:.5g} s^-1 ({100*(rate_gui/uc.rate-1):+.2f} %, "
              f"stat ±{100*err:.2f} %; unweighted n/T would be "
              f"{100*(len(fr.p)/T_gui/uc.rate-1):+.1f} %)")

    print()
    if failures:
        print(f"{len(failures)} FAILED")
        return 1
    print("ALL AGREE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
