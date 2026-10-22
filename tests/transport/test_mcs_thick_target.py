#!/usr/bin/env python3
"""Multiple scattering in thick rock: the thick-target limit (mock-review item 7;
each check fails on 1.4.1).

UCMuon-MC and the Bethe-Bloch engines accumulate Highland's variance along the
path with the log term stopping at T_SAT_X0 = 8.4 X0, where 13.6 MeV F =
14.7 MeV: once theta0 is far above the nuclear-size angle the width grows like
sqrt(x) (benchmark/analysis/mcs_thick_target.py). Up to 1.4.1 the log stopped
at 1 X0 (13.6 MeV, 5-9 % below Geant4 in Standard Rock), and gui/mcs_margin.py
let it grow over the whole path (17.5 MeV at 200 m) with a 2 sigma_r margin.

  - the constants, and the asymptote of _ms_theta0 (exactly 14.7 MeV);
  - the thick-target value from first principles, 14.66-14.82 MeV from water
    to iron (mcs_thick_target.py);
  - UCMuon-MC, 100 GeV through 50 m of rock: core width within 3 % of
    14.7 MeV on the CSDA momentum profile;
  - mcs_margin: sigma_r at 200 m equals the thick-target integral, and the
    suggested margin is 2.5 sigma_r;
  - the Fortran Bethe-Bloch engine (bin/ucmuon_transport_bb_omp), the same
    thick layer; skipped without the binary.

    python3 tests/transport/test_mcs_thick_target.py      # exit 0 = all pass
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
sys.path.insert(0, str(ROOT / "benchmark" / "analysis"))
import fast_flux_estimator as ff                     # noqa: E402
import mcs_margin                                    # noqa: E402
import mcs_thick_target                              # noqa: E402
import ucmuon_stochastic_driver as U                 # noqa: E402

FAILS: list[str] = []
M_MU, X0, RHO = 0.1056584, 26.54, 2.65


def check(ok: bool, what: str) -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILS.append(what)


def csda_sum(e_tot, depth_m, n=4000):
    """sum dx/X0/(beta p)^2 [GeV^-2] along a vertical path, CSDA momentum."""
    r0 = ff._R_of_T(e_tot - M_MU)
    x = np.linspace(0.0, depth_m * 100.0 * RHO, n + 1)
    e = ff._T_of_R(np.maximum(r0 - 0.5 * (x[1:] + x[:-1]), 1e-3)) + M_MU
    return float(np.sum(np.diff(x) / X0 / (np.sqrt(e * e - M_MU * M_MU) ** 2 / e) ** 2))


def main() -> int:
    print("model")
    check(U.HIGHLAND_MEV == 13.6 and U.ES_THICK_MEV == 14.7 and abs(U.T_SAT_X0 - 8.402) < 1e-3,
          f"T_SAT_X0 = {U.T_SAT_X0:.3f} X0, E_s,thick = {U.ES_THICK_MEV} MeV")
    S = t = tot = 0.0
    for _ in range(20000):
        th, S, t = U._ms_theta0(S, t, 0.1, 100.0)
        tot += th ** 2
    es = np.sqrt(tot) * 100.0 / np.sqrt(t) * 1e3
    check(abs(es - 14.7) < 1e-6, f"_ms_theta0 over 2000 X0 at beta p = 100 GeV: E_s = {es:.4f} MeV")
    th1, _, _ = U._ms_theta0(0.0, 0.0, 1.0, 10.0)
    check(abs(th1 / (13.6e-3 / 10.0) - 1.0) < 1e-12, "one X0: Highland's 13.6 MeV, unchanged")
    vals = {k: mcs_thick_target.es_thick_mev(*v) for k, v in mcs_thick_target.MATERIALS.items()}
    check(all(14.6 < vals[k] < 14.9 for k in ("Standard Rock", "water", "iron")),
          "thick-target E_s from first principles: " + ", ".join(f"{k} {v:.2f}" for k, v in vals.items()))

    print("UCMuon-MC")
    n = 20000
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "b.dat"
        p = np.sqrt(100.0 ** 2 - M_MU ** 2)
        f.write_text("".join(f"{i + 1} 0 0 0 {p:.6f} 0 0 {-p:.6f} 0 0 100.0 1 1 0\n" for i in range(n)))
        mu = U._read_input(str(f), True)
        r = U.transport(mu, 50.0, RHO, dict(U._MAT_DB[1]), rng=np.random.default_rng(3))
    a = r["alive"] == 1
    th = np.arctan2(np.hypot(r["cx_f"][a], r["cy_f"][a]), -r["cz_f"][a])
    t0 = np.median(th) / np.sqrt(2.0 * np.log(2.0))
    es = t0 / np.sqrt(csda_sum(100.0, 50.0)) * 1e3
    check(abs(es / 14.7 - 1.0) < 0.03, f"100 GeV through 50 m of rock: core E_s = {es:.2f} MeV (14.7)")

    print("mcs_margin")
    sr = mcs_margin.sigma_r_cm(300.0, 20000.0)
    s = np.linspace(0.0, 20000.0, 4001)
    r0 = ff._R_of_T(300.0 - M_MU)
    e = ff._T_of_R(r0 - RHO * s) + M_MU
    bp = (e * e - M_MU * M_MU) / e
    f2 = (20000.0 - s) ** 2 * (14.7e-3 / bp) ** 2 / (X0 / RHO)
    ref = np.sqrt(2.0 * np.sum(0.5 * (f2[1:] + f2[:-1]) * np.diff(s)))
    check(abs(sr / ref - 1.0) < 0.01, f"sigma_r, 300 GeV through 200 m: {sr:.2f} cm, thick-target {ref:.2f} cm")
    rec = mcs_margin.suggested_margin(200.0, RHO, 300.0)
    check(rec["k"] == 2.5 and abs(rec["margin_cm"] / rec["sigma_r_cm"] - 2.5) < 1e-12,
          f"suggested margin = {rec['margin_cm'] / rec['sigma_r_cm']:.1f} sigma_r")

    print("Bethe-Bloch (Fortran, OpenMP)")
    bb = ROOT / "bin" / "ucmuon_transport_bb_omp"
    if not bb.exists():
        print("  skipped: build it with make ucmuon_transport_bb_omp")
    else:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "b.dat"
            p = np.sqrt(100.0 ** 2 - M_MU ** 2)
            # spread in x and y: the engine takes the source plane from it
            src.write_text("".join(f"{i + 1} {i % 100 * 10.0:.1f} {i // 100 % 100 * 10.0:.1f} 0 {p:.6f} "
                                   f"0 0 {-p:.6f} 0 0 100.0 1 1 0\n" for i in range(n)))
            out = Path(tmp) / "o.dat"
            subprocess.run([str(bb)], input=f"{src}\n{out}\n1\n14\n50\n1\n1\n", text=True,
                           capture_output=True, check=True, cwd=tmp,
                           env={"OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
            rows = np.loadtxt(out, comments="#")
        rows = rows[rows[:, 8] == 1]
        th = np.arctan2(np.hypot(rows[:, 13], rows[:, 14]), -rows[:, 15])
        t0 = np.median(th) / np.sqrt(2.0 * np.log(2.0))
        es = t0 / np.sqrt(csda_sum(100.0, 50.0)) * 1e3
        check(abs(es / 14.7 - 1.0) < 0.03, f"100 GeV through 50 m of rock: core E_s = {es:.2f} MeV (14.7)")

    print()
    print("ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
