#!/usr/bin/env python3
"""Command-line converters after 1.3.2 (the electron checks fail on 1.3.2).

  - src/converters/ucmuon_to_geant4.py, ucmuon_to_phits.py and
    ucmuon_to_phits.f90 wrote every row as a muon (PDG +-13, E - m_mu). A
    spectrum-8 file (cosmic e+-) has no particle column, so its electrons
    became muons, 99 % of them with zero kinetic energy (E < 105.7 MeV). The
    converters now take --electrons ("electrons" for the Fortran one) and
    refuse a row whose own E and p contradict the particle being written, in
    either direction.
  - Muon files convert as before: PDG +-13 (mu+ = -13), Ekin = E - m_mu.

Needs bin/ucmuon_gen_omp (make ucmuon_gen_omp); the Fortran converter part is
skipped when bin/ucmuon_to_phits is absent (make ucmuon_to_phits).

    python3 tests/converters/test_converters.py      # exit 0 = all pass
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "flux"))
from _harness import Config  # noqa: E402

CONV = ROOT / "src" / "converters"
F90 = ROOT / "bin" / "ucmuon_to_phits"
M_E, M_MU = 0.51099895e-3, 0.105658370
FAILS: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILS.append(what)


def generate(tmp: Path, spectrum: int, emin: float, emax: float, name: str) -> Path:
    cfg = Config(spectrum=spectrum, emin=emin, emax=emax, theta_max_deg=60.0, radius_m=5.0,
                 n=1000, seed=11)
    out = tmp / name
    env = {**os.environ, "UCMUON_SEED": "11", "OMP_NUM_THREADS": "1"}
    proc = subprocess.run([str(ROOT / "bin" / "ucmuon_gen_omp")], input=cfg.stdin(str(out), str(out) + ".sel"),
                          capture_output=True, text=True, env=env, cwd=tmp)
    if proc.returncode != 0:
        raise RuntimeError(proc.stdout[-2000:])
    return out


def load(path: Path) -> np.ndarray:
    return np.array([[float(x.replace("D", "E")) for x in line.split()]
                     for line in open(path) if line.strip() and not line.startswith("#")])


def run(kind: str, src: Path, dst: Path, electrons: bool) -> subprocess.CompletedProcess:
    if kind == "f90":
        with open(src) as fin, open(dst, "w") as fout:
            return subprocess.run([str(F90), "gen"] + (["electrons"] if electrons else []),
                                  stdin=fin, stdout=fout, stderr=subprocess.PIPE, text=True)
    return subprocess.run([sys.executable, str(CONV / kind), str(src), str(dst)]
                          + (["--electrons"] if electrons else []), capture_output=True, text=True)


def main() -> int:
    if not (ROOT / "bin" / "ucmuon_gen_omp").exists():
        print("build the generator first: make ucmuon_gen_omp")
        return 1
    kinds = ["ucmuon_to_geant4.py", "ucmuon_to_phits.py"] + (["f90"] if F90.exists() else [])
    if "f90" not in kinds:
        print("  (bin/ucmuon_to_phits absent: Fortran converter skipped)")
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        el = generate(tmp, 8, 0.01, 1.0, "electrons.dat")
        mu = generate(tmp, 4, 1.0, 1000.0, "muons.dat")
        for src, electrons, mass in ((el, True, M_E), (mu, False, M_MU)):
            raw = np.loadtxt(src, comments="#", ndmin=2)
            base = 11 if electrons else 13
            for kind in kinds:
                dst = tmp / f"{kind}_{src.stem}.out"
                proc = run(kind, src, dst, electrons)
                t_ = load(dst) if proc.returncode == 0 else np.empty((0, 8))
                ids = t_[:, 0].astype(int)
                dev = np.max(np.abs(t_[:, 7] - (raw[:, 10] - mass) * 1e3)) if len(t_) == len(raw) else np.inf
                check(proc.returncode == 0 and len(t_) == len(raw)
                      and np.array_equal(ids, np.where(raw[:, 11] == 1, -base, base)) and dev < 1e-3,
                      f"{kind}, {src.stem}: PDG +-{base} with the right sign, "
                      f"max |Ekin - (E - m)| = {dev:.1e} MeV")
                # the same file as the other particle is refused
                bad = tmp / f"{kind}_{src.stem}_wrong.out"
                proc = run(kind, src, bad, not electrons)
                gone = kind == "f90" or not bad.exists()
                check(proc.returncode == 2 and gone,
                      f"{kind}, {src.stem} {'without' if electrons else 'with'} the electron option: "
                      f"refused (exit {proc.returncode})")
    print("ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
