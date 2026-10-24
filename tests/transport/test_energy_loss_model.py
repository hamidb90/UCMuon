#!/usr/bin/env python3
"""UCMuon-MC energy-loss model after 1.4.1 (mock-review item 8; the checks
marked * fail on 1.4.1).

  - * the defaults: hard_spectrum 'kkp' and v_cut = 0.01, in transport() and
    on the driver's standard input (1.4.1: 'proc' and 0.05);
  - * the tabulated spectra (gui/ucmuon_radiative_spectra.npz, written by
    tools/tabulate_radiative_spectra.py with PROPOSAL): present for rock,
    water and iron, shares of the mean loss above v_cut within [0, 1] and
    falling with v_cut, CDFs monotone, and at 100 GeV in rock the shares
    above v = 0.05 that PROPOSAL gave (brems 90 %, pair 10 %, photonuclear
    80 %);
  - the mean loss: the expected loss of a step, continuous part plus the
    mean of the sampled transfers, equals the PDG 2024 table within 0.02 %
    from 3 GeV to 30 TeV in the three materials, for v_cut 0.01 and 0.05;
  - * knock-on electrons of a spin-1/2 projectile: with radiative losses off,
    100 GeV muons through 6000 g/cm^2 lose the table's mean, and the variance
    of their loss is the spin-1/2 one, 27 % above the spin-0 form 1.4.1 used;
  - the range cut: the deterministic-loss range is an upper bound on the
    mean-loss range at every energy (it must be, or the cut would kill muons
    that can survive).

    python3 tests/transport/test_energy_loss_model.py     # exit 0 = all pass
"""
from __future__ import annotations

import inspect
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
import ucmuon_stochastic_driver as U                 # noqa: E402

FAILS: list[str] = []
M_MU = U.M_MU                                        # MeV


def check(ok: bool, what: str) -> None:
    print(f"  [{'ok  ' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILS.append(what)


def beam_file(path: Path, e_tot_gev: float, n: int) -> None:
    p = np.sqrt(e_tot_gev ** 2 - (M_MU / 1000.0) ** 2)
    with open(path, "w") as fh:
        fh.write("# EventID x_cm y_cm z_cm p_GeV px_GeV py_GeV pz_GeV theta_rad phi_rad "
                 "E_GeV charge hit_flag det_mask\n")
        for i in range(n):
            fh.write(f"{i + 1} {i % 300 * 10.0:.1f} {i // 300 * 10.0:.1f} 0.0 {p:.6f} 0.0 0.0 "
                     f"{-p:.6f} 0.0 0.0 {e_tot_gev:.6f} 1 1 0\n")


def defaults(tmp: Path) -> None:
    print("defaults")
    sig = inspect.signature(U.transport).parameters
    check(sig["hard_spectrum"].default == "kkp" and sig["v_cut"].default == 0.01,
          f"transport(): hard_spectrum={sig['hard_spectrum'].default!r}, "
          f"v_cut={sig['v_cut'].default}")
    src, out = tmp / "beam.dat", tmp / "out.dat"
    beam_file(src, 10.0, 20)
    r = subprocess.run([sys.executable, str(ROOT / "gui" / "ucmuon_stochastic_driver.py")],
                       input=f"{src}\n{out}\n1.0\n2.65\n", text=True, capture_output=True)
    log = r.stdout + r.stderr
    check(r.returncode == 0 and "v_cut=0.010" in log and "hard_spectrum=kkp" in log,
          "driver standard input with the options left out: v_cut 0.01, kkp")


def tables() -> None:
    print("tabulated spectra")
    check(U._RAD_FILE.exists(), f"{U._RAD_FILE.name} present")
    z = np.load(U._RAD_FILE)
    check("PROPOSAL" in str(z["provenance"]) and "KelnerKokoulinPetrukhin" in str(z["provenance"]),
          "provenance: PROPOSAL, Kelner-Kokoulin-Petrukhin")
    for mat in ("rock", "water", "iron"):
        t1, t5 = U._kkp_tables(mat, 0.01), U._kkp_tables(mat, 0.05)
        ok = True
        for proc in ("brems", "pair", "nucl"):
            ln_e, share, n_hard, v_grid, cdf = t1[proc]
            ok &= bool(np.all((share >= 0) & (share <= 1)) and np.all(share >= t5[proc][1] - 1e-6))
            ok &= bool(np.all(np.diff(cdf, axis=1) >= 0) and np.allclose(cdf[:, -1], 1.0))
            ok &= bool(np.all(n_hard >= t5[proc][2]) and np.isclose(v_grid[0], 0.01))
        check(ok, f"{mat}: shares in [0, 1] and falling with v_cut, CDFs monotone to 1")
    t5 = U._kkp_tables("rock", 0.05)
    s = {p: float(np.interp(np.log(100e3 + M_MU), t5[p][0], t5[p][1])) for p in t5}
    check(0.85 < s["brems"] < 0.95 and 0.07 < s["pair"] < 0.13 and 0.75 < s["nucl"] < 0.88,
          "rock, 100 GeV, above v = 0.05: brems {brems:.3f}, pair {pair:.3f}, "
          "photonuclear {nucl:.3f}".format(**s))


def mean_loss() -> None:
    print("mean loss of a step against the PDG 2024 table")
    u = (np.arange(100_000) + 0.5) / 100_000
    for key, mat in (("rock", U._MAT_DB[1]), ("water", U._MAT_DB[2]), ("iron", U._MAT_DB[4])):
        for v_cut in (0.01, 0.05):
            tab = U._kkp_tables(mat["proc"], v_cut)
            worst = 0.0
            for t in np.geomspace(3e3, 3e7, 37):
                e_tot = t + M_MU
                a, lb, lp, ln = (float(np.ravel(x)[0]) for x in U._loss_components(np.array([t]), mat))
                model = a
                for proc, L in (("brems", lb), ("pair", lp), ("nucl", ln)):
                    ln_e, share, n_hard, v_grid, cdf = tab[proc]
                    kb = int(U._kkp_bin(ln_e, np.array([np.log(e_tot)]))[0])
                    s, n = share[kb], n_hard[kb]          # as the transport takes them
                    v_mean = U._kkp_sample(tab[proc], np.full(len(u), np.log(e_tot)), u).mean()
                    model += (1.0 - s) * L + n * L * v_mean
                worst = max(worst, abs(model / (a + lb + lp + ln) - 1.0))
            check(worst < 2e-4, f"{key}, v_cut {v_cut}: within {100 * worst:.4f} % from 3 GeV to 30 TeV")


def knock_on(tmp: Path) -> None:
    print("knock-on electrons, spin 1/2 (radiative losses off, no scattering)")
    n, e_tot, x = 90_000, 100.0, 6000.0               # GeV, g/cm^2
    mat = dict(U._MAT_DB[1])
    mat.pop("proc")                                   # a custom material: rock's table
    mat["b_rad"] = 0.0                                # scaled by b_rad -> no radiative loss
    src = tmp / "beam100.dat"
    beam_file(src, e_tot, n)
    mu = U._read_input(str(src), True)
    r = U.transport(mu, x / 265.0, 2.65, mat, ms_enable=False, rng=np.random.default_rng(3))
    loss = (e_tot * 1000.0 - M_MU) - r["E_kin_f_MeV"]
    # prediction along the mean-loss path: mean = int a dX; variance from the
    # close collisions above T_cut, with and without the T^2/(2E^2) term
    h = x / 600
    t, mean, var1, var0 = e_tot * 1000.0 - M_MU, 0.0, 0.0, 0.0
    for _ in range(600):
        a = float(np.ravel(U._loss_components(np.array([t]), mat)[0])[0])
        k, tm, b2 = (float(np.asarray(q).ravel()[0]) for q in U._delta_kinematics(np.array([t]), mat["ZoverA"]))
        tc, e2 = U.T_CUT_DELTA, (t + M_MU) ** 2
        v0 = k * ((tm - tc) - b2 * (tm ** 2 - tc ** 2) / (2.0 * tm))
        var0 += v0 * h
        var1 += (v0 + k * (tm ** 3 - tc ** 3) / (6.0 * e2)) * h
        mean += a * h
        t -= a * h
    m_obs, v_obs = loss.mean(), loss.var()
    check(abs(m_obs / mean - 1.0) < 3e-3,
          f"mean loss {m_obs / 1000:.4f} GeV, table {mean / 1000:.4f} GeV ({100 * (m_obs / mean - 1):+.2f} %)")
    check(abs(v_obs / var1 - 1.0) < 0.10 and v_obs / var0 > 1.12,
          f"variance {v_obs / var1:.3f} x the spin-1/2 prediction, {v_obs / var0:.3f} x spin 0 "
          f"(spin 1/2 / spin 0 = {var1 / var0:.3f})")


def range_bound() -> None:
    print("range cut")
    t = np.geomspace(1e3, 1e7, 200)
    for key in (1, 2, 4):
        mat = U._MAT_DB[key]
        r_mean = U._mean_range(t, mat)
        ok = all(bool(np.all(U._det_range(t, mat, vc, delta_rays=dr) >= r_mean * (1 - 1e-9)))
                 for vc in (0.01, 0.05) for dr in (False, True))
        check(ok, f"{mat['name']}: deterministic-loss range >= mean-loss range, 1 GeV to 10 TeV")


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        defaults(tmp)
        tables()
        mean_loss()
        knock_on(tmp)
        range_bound()
    print(f"\n{'ALL PASS' if not FAILS else f'{len(FAILS)} FAILED'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
