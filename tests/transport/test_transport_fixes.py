#!/usr/bin/env python3
"""Fortran transport engine fixes after 1.3.2 (each check fails on 1.3.2).

  - Bethe-Bloch (bin/ucmuon_transport_bb_omp): ice and water had Z/A 0.495
    (Z=7.42, A=14.99; water's is 0.555) and b = 3.40e-6 (PDG 2.28e-6), so
    10 GeV muons lost 10 % too little energy in water. Checked against the
    PDG 2024 water table, within the engine's 2 % analytic-ionisation
    residual (the same as in rock).
  - Multiple scattering (Bethe-Bloch, OpenMP and MPI): Highland was applied
    per 10 g/cm2 step and the steps added in quadrature, which depends on
    the step size; one radiation length of rock gave 11 % less than Highland
    for that layer. The variance is now accumulated along the path.
  - Point-like sources (a pencil beam, all muons at one point): the engines
    guessed the source plane from the spread of the start positions, and a
    tie went to Y, so vertical muons were transported sideways and none
    survived. They now take the depth axis from the momentum. Checked for
    Bethe-Bloch, MUSIC (bin/ucmuon_transport_music_omp, needs the MUSIC
    tables) and the MPI Bethe-Bloch (bin/ucmuon_transport_bb, needs mpirun);
    each part is skipped when its binary is absent.

    python3 tests/transport/test_transport_fixes.py      # exit 0 = all pass
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
BIN = ROOT / "bin"
FAILS: list[str] = []
M_MU = 0.10566


def check(ok: bool, what: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILS.append(what)


def beam(path: Path, n: int, e_gev: float, spread_cm: float) -> None:
    """n vertical muons of total energy e_gev, start positions spread over
    +-spread_cm in x and y (0 = a pencil beam at the origin)."""
    p = np.sqrt(e_gev ** 2 - M_MU ** 2)
    xs = np.linspace(-spread_cm, spread_cm, n)
    with open(path, "w") as f:
        f.write("# EventID x y z p px py pz theta phi E charge det_mask\n")
        for i in range(n):
            f.write(f"{i + 1:10d} {xs[i]:.3f} {xs[::-1][i]:.3f} 0.0 {p:.6f} 0.0 0.0 {-p:.6f} "
                    f"0.0 0.0 {e_gev:.6f} {1 if i % 2 else -1} 1\n")


def survivors(files) -> np.ndarray:
    d = np.concatenate([np.loadtxt(f, comments="#", ndmin=2) for f in files])
    return d[d[:, 8] == 1]


def pdg_exit_total(key: str, e_gev: float, x_gcm2: float) -> float:
    """Total energy [GeV] after x_gcm2 (PDG 2024 table, CSDA, step by step)."""
    import ucmuon_stochastic_driver as U
    t, T = U._PROC_TABLES[key], U._PROC_T
    tot = t["ION"] + t["BREMS"] + t["PAIR"] + t["NUCL"]
    tk, h = e_gev * 1e3 - 105.6584, 0.5
    for _ in range(int(x_gcm2 / h)):
        tk -= np.interp(np.log(tk), np.log(T), tot) * h
    return (tk + 105.6584) / 1e3


def highland_check(rows: np.ndarray, what: str) -> None:
    """Projected RMS angle of 10 GeV survivors after one X0 vs Highland."""
    p = np.sqrt(10.0 ** 2 - M_MU ** 2)
    theta0 = 13.6e-3 / (p * p / 10.0)
    rms = np.sqrt(np.mean(np.arctan2(rows[:, 13], -rows[:, 15]) ** 2))
    check(abs(rms / theta0 - 1.0) < 0.03,
          f"{what}: 10 GeV through one X0 of rock, RMS angle {rms * 1e3:.4f} mrad vs Highland "
          f"{theta0 * 1e3:.4f}")


def run(cmd, stdin, cwd):
    p = subprocess.run(cmd, input=stdin, capture_output=True, text=True, cwd=cwd,
                       env={**os.environ, "OMP_NUM_THREADS": "1"}, timeout=600)
    if p.returncode != 0:
        raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])
    return p.stdout


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        beam(tmp / "v10.dat", 200, 10.0, 100.0)
        beam(tmp / "point.dat", 200, 100.0, 0.0)
        beam(tmp / "ms.dat", 20000, 10.0, 100.0)
        x0_m = 26.54 / 2.65 / 100.0                     # one X0 of Standard Rock

        bb = BIN / "ucmuon_transport_bb_omp"
        print("Bethe-Bloch (OpenMP)")
        if not bb.exists():
            print("  skipped: build it with make ucmuon_transport_bb_omp")
        else:
            run([str(bb)], f"{tmp}/v10.dat\n{tmp}/w.dat\n0\n13\n10\n3 1.0\n0\n\n", tmp)
            e = survivors([tmp / "w.dat"])[:, 12].mean()
            ref = pdg_exit_total("water", 10.0, 1000.0)
            loss, ref_loss = 10.0 - e, 10.0 - ref
            check(abs(loss / ref_loss - 1.0) < 0.02,
                  f"water, 10 GeV through 10 m: energy loss {loss:.4f} GeV vs PDG {ref_loss:.4f}")
            run([str(bb)], f"{tmp}/point.dat\n{tmp}/p.dat\n0\n13\n100\n1 2.65\n0\n\n", tmp)
            s = survivors([tmp / "p.dat"])
            check(len(s) == 200, f"pencil beam, 100 GeV through 100 m rock: {len(s)}/200 survive")
            run([str(bb)], f"{tmp}/ms.dat\n{tmp}/s.dat\n0\n13\n{x0_m}\n1 2.65\n1\n\n", tmp)
            highland_check(survivors([tmp / "s.dat"]), "multiple scattering")

        mu = BIN / "ucmuon_transport_music_omp"
        tables = sorted(glob.glob(str(BIN / "music-*rock*.dat")))
        print("MUSIC (OpenMP)")
        if not mu.exists() or len(tables) < 3:
            print("  skipped: needs bin/ucmuon_transport_music_omp and the MUSIC rock tables")
        else:
            for t in glob.glob(str(BIN / "music-*.dat")):
                os.symlink(os.path.realpath(t), tmp / os.path.basename(t))
            run([str(mu)], f"{tmp}/point.dat\n{tmp}/m.dat\n2.65\n26.48\n100\n1\n1\n-30\n1\n1\n0\n\n", tmp)
            s = survivors([tmp / "m.dat"])
            check(len(s) > 180, f"pencil beam, 100 GeV through 100 m rock: {len(s)}/200 survive")

        mpi = BIN / "ucmuon_transport_bb"
        print("Bethe-Bloch (MPI)")
        if not mpi.exists() or shutil.which("mpirun") is None:
            print("  skipped: needs bin/ucmuon_transport_bb (make ucmuon_transport_bb) and mpirun")
        else:
            run(["mpirun", "--oversubscribe", "-np", "2", str(mpi)],
                f"{tmp}/point.dat\n{tmp}/mp\n1\n100\n0\n0\n", tmp)
            s = survivors(sorted(glob.glob(str(tmp / "mp_*.dat"))))
            check(len(s) == 200, f"pencil beam, 100 GeV through 100 m rock: {len(s)}/200 survive")
            run(["mpirun", "--oversubscribe", "-np", "2", str(mpi)],
                f"{tmp}/ms.dat\n{tmp}/ms\n1\n{x0_m}\n1\n0\n", tmp)
            highland_check(survivors(sorted(glob.glob(str(tmp / "ms_*.dat")))),
                           "multiple scattering")

    print("ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
