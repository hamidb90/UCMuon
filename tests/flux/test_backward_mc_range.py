#!/usr/bin/env python3
"""The backward MC says when its range table cuts a flux, instead of a silent zero.

gui/ucmuon_backward_mc.py maps a detector energy back to the surface through a
CSDA range table that ends at 2 TeV (3616 m w.e. of standard rock). Muons whose
surface energy would be beyond it are left out, so under thick rock the flux is
a lower bound, and zero once every detector energy is beyond the table. Up to
v1.1.2 that was silent. This checks, for Guan (spectrum 3), θ = 60°:

  - no RangeTableWarning where the cut is below 1 % (500 and 1000 m w.e.);
  - a warning with the estimated fraction at 2000 and 3000 m w.e., which
    matches the fraction measured by continuing the table (7 % and 37 %,
    ±5 points), and a flux that is unchanged by the warning;
  - a warning that says "0" where the flux is exactly zero (4000 m w.e.);
  - backward_mc_flux (the cone integral of the backward-MC tab) carries the
    note in info["range_note"] and warns; a shallow run carries none;
  - collect_range_warnings() gathers them without printing, passes other
    warnings on, and range_summary() counts them;
  - the terrain engine's compute_flux_map prints one summary line, and warns
    when a spectrum is used below its fitted range.

    python3 tests/flux/test_backward_mc_range.py      # exit 0 = all pass
"""
from __future__ import annotations

import contextlib
import io
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


def flux_and_warnings(X, theta_deg=60.0, E_max=2500.0):
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        f = bmc.directional_flux(X, math.radians(theta_deg), 3, E_min_GeV=1.0,
                                 E_max_GeV=E_max, n_E=60, mode=1)
    return f, [w for w in rec if w.category.__name__ == "RangeTableWarning"]


def main() -> int:
    print("directional_flux (Guan, θ = 60°, 1-2500 GeV)")
    for X in (5e4, 1e5):
        f, w = flux_and_warnings(X)
        check(f > 0 and not w, f"{X/100:.0f} m w.e.: flux {f:.3e}, no warning")

    for X, expect in ((2e5, 0.07), (3e5, 0.37)):
        f, w = flux_and_warnings(X)
        frac = w[0].message.fraction if w else float("nan")
        check(len(w) == 1 and abs(frac - expect) < 0.05 and "lower bound" in str(w[0].message),
              f"{X/100:.0f} m w.e.: one warning, {100*frac:.0f} % left out (expected ≈ {100*expect:.0f} %)")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f2 = bmc.directional_flux(X, math.radians(60.0), 3, E_min_GeV=1.0,
                                      E_max_GeV=2500.0, n_E=60, mode=1)
        check(f == f2, f"{X/100:.0f} m w.e.: the flux returned does not depend on the warning filter")

    f, w = flux_and_warnings(4e5)
    check(f == 0.0 and len(w) == 1 and w[0].message.zero and "is 0" in str(w[0].message),
          "4000 m w.e.: flux 0 and a warning that says so")

    print("\nbackward_mc_flux")
    for depth, want in ((100.0, False), (1600.0, True)):
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            res = bmc.backward_mc_flux(depth_m=depth, rho=2.65, mat_id=1, spectrum_mode=3,
                                       E_min_GeV=1.0, E_max_GeV=1000.0, theta_max_deg=30.0,
                                       n_E=40, n_theta=6, mode=1)
        got = [w for w in rec if w.category.__name__ == "RangeTableWarning"]
        note = res["info"].get("range_note")
        check(bool(note) == want and bool(got) == want,
              f"{depth:.0f} m of rock: range_note {'present' if note else 'absent'}, "
              f"{len(got)} warning(s) (expected {'one' if want else 'none'})")

    print("\ncollect_range_warnings / range_summary")
    with warnings.catch_warnings(record=True) as outer:
        warnings.simplefilter("always")
        with bmc.collect_range_warnings() as got:
            for X in (5e4, 3e5, 4e5):
                bmc.directional_flux(X, math.radians(60.0), 3, E_min_GeV=1.0,
                                     E_max_GeV=2500.0, n_E=60, mode=1)
            warnings.warn("an unrelated warning", UserWarning)
    leaked = [w for w in outer if w.category.__name__ == "RangeTableWarning"]
    passed = [w for w in outer if "unrelated" in str(w.message)]
    check(len(got) == 2 and not leaked, f"collected {len(got)} of 3 calls, none printed")
    check(len(passed) == 1, "other warnings are passed on")
    summ = bmc.range_summary(got, 3) or ""
    check("2 of 3 directions" in summ and "1 are 0" in summ, f"summary: {summ[:90]}…")
    check(bmc.range_summary([], 3) is None, "no summary when nothing was cut")

    print("\nterrain engine: compute_flux_map")
    ze = np.array([10.0, 30.0, 60.0])
    ob = np.array([[0.0, 1e5, 4e5]])               # open sky, 1000 m w.e., 4000 m w.e.
    sky = np.array([[True, False, False]])
    out = io.StringIO()
    with contextlib.redirect_stdout(out), warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        fm, osky = td.compute_flux_map(np.array([0.5]), ze, ob, sky, 2.65, spectrum_mode=3,
                                       mode=1, n_E=40, script_dir=ROOT / "gui")
    lines = [l for l in out.getvalue().splitlines() if "WARNING" in l]
    check(len(lines) == 1 and "1 of 3 directions" in lines[0] and fm[0, 2] == 0.0,
          f"one summary line for the map: {lines[0].strip()[:80] if lines else '(none)'}…")
    check(abs(fm[0, 0] / osky[0] - 1.0) < 1e-12, "open-sky bin keeps the open-sky flux")

    out = io.StringIO()
    with contextlib.redirect_stdout(out), warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        td.compute_flux_map(np.array([0.5]), ze[:1], ob[:, :1], sky[:, :1], 2.65,
                            spectrum_mode=3, mode=1, E_min_GeV=0.5, n_E=20,
                            script_dir=ROOT / "gui")
    check("fitted to data above 1 GeV" in out.getvalue(), "E_min = 0.5 GeV with Guan warns")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        td.compute_flux_map(np.array([0.5]), ze[:1], ob[:, :1], sky[:, :1], 2.65,
                            spectrum_mode=3, mode=1, n_E=20, script_dir=ROOT / "gui")
    check("WARNING" not in out.getvalue(), "the default E_min (1 GeV) does not warn")

    print()
    if failures:
        print(f"{len(failures)} FAILED")
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
