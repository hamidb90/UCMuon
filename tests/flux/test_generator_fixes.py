#!/usr/bin/env python3
"""Generator fixes after 1.3.2 (each check fails on 1.3.2).

  - PARMA "W from the date" (W <= -999): a date after the last year of
    FFPtable.day got W = 0 (solar minimum) silently, because getHP took the
    0.0 its table is initialised with for later years as neutron-monitor
    data. It now falls back to the latest day with data and says so, as it
    already did for a missing day inside the table.
  - Random deviates (src/generator/rng_parallel.f90 and its port in
    UCMuGen.h): only the top 24 bits of the LCG were used, as a real(4), so
    every deviate was a multiple of 2^-24 and the far momentum tail was cut
    off and quantised (Guan 1-10000 GeV/c never above 5923 GeV/c). Both now
    return 53 bits; ucmugen/validation/compare_legacy.py checks that they
    still agree bit for bit. Needs gfortran and a C++17 compiler.
  - MPI generator (bin/ucmuon_gen; skipped when it or mpirun is absent): it
    raised an E_min below 1 GeV to 1 GeV and turned spectrum 8 (electrons)
    into CosmoALEPH muons, both silently, and without a detector its files
    held the few extra muons the threads drew past the target. It now does
    what the OpenMP generator does. Every rank numbered its muons from 1, so
    the merged files (hpc/run_ucmuon_gen.sh concatenates them) repeated
    EventIDs; they are now unique, and the selected file runs 1..N.

    python3 tests/flux/test_generator_fixes.py      # exit 0 = all pass
"""
from __future__ import annotations

import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from _harness import ROOT, Config, Detector, run_fortran

MPI_BIN = ROOT / "bin" / "ucmuon_gen"
FAILS: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILS.append(what)


def parma_w(year: int, month: int, day: int) -> tuple[float, str]:
    cfg = Config(spectrum=3, emin=1.0, emax=100.0, theta_max_deg=60.0, radius_m=10.0, n=300,
                 seed=3)
    cfg.parma.update(year=year, month=month, day=day, sw=-999.0)
    out = run_fortran(cfg, threads=1).stdout
    w = float(re.search(r"PARMA W index:\s*([-0-9.]+)", out).group(1))
    return w, out


def run_mpi(cfg: Config, nranks: int = 2, threads: int = 2):
    """Run bin/ucmuon_gen; returns (rows of all rank files, stdout). With a
    detector it also writes the surface file and returns (surface rows,
    selected rows, stdout). Rank files are concatenated in rank order, as
    hpc/run_ucmuon_gen.sh does."""
    with tempfile.TemporaryDirectory() as tmp:
        out_all = os.path.join(tmp, "all.dat")
        out_sel = os.path.join(tmp, "sel.dat")
        lines = cfg.stdin(out_all, out_sel).split("\n")
        del lines[lines.index(out_all) - 1]          # no PHITS question in the MPI binary
        if cfg.detector is not None:
            lines[lines.index(out_all) - 1] = "1"    # save_all: the surface file too
        env = {**os.environ, "UCMUON_SEED": str(cfg.seed), "OMP_NUM_THREADS": str(threads)}
        proc = subprocess.run(["mpirun", "--oversubscribe", "-np", str(nranks), str(MPI_BIN)],
                              input="\n".join(lines).rstrip("\n") + "\n", capture_output=True,
                              text=True, env=env, cwd=tmp, timeout=600)
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout[-3000:] + proc.stderr[-2000:])
        def merged(path):
            return np.concatenate([np.loadtxt(f, comments="#", ndmin=2)
                                   for f in sorted(glob.glob(path[:-4] + "_*.dat"))])
        if cfg.detector is not None:
            return merged(out_all), merged(out_sel), proc.stdout
        return merged(out_all), proc.stdout


def deviates_on_24_bit_grid() -> tuple[int, int, int]:
    """(draws, Fortran draws on the 2^-24 grid, UCMuGen draws on it) for
    20000 draws of par_ranlux and of UCMuonLcg::flat() with the same seed."""
    n = 20000
    with tempfile.TemporaryDirectory() as tmp:
        f90 = os.path.join(tmp, "draw.f90")
        for kind in (8, 4):                    # the argument was real(4) up to 1.3.2
            with open(f90, "w") as f:
                f.write("program draw\n  use rng_parallel\n  implicit none\n"
                        f"  real({kind}) :: r\n  integer :: i, k\n  call par_init_rng(4321)\n"
                        f"  k = 0\n  do i = 1, {n}\n    call par_ranlux(r)\n"
                        "    if (r <= 0 .or. r >= 1) stop 2\n"
                        "    if (dble(r) * 16777216d0 == aint(dble(r) * 16777216d0)) k = k + 1\n"
                        "  end do\n  print *, k\nend program draw\n")
            if subprocess.run(["gfortran", "-O2", "-fopenmp", "-J", tmp,
                               str(ROOT / "src/generator/rng_parallel.f90"), f90,
                               "-o", os.path.join(tmp, "draw")], capture_output=True).returncode == 0:
                break
        k_f = int(subprocess.run([os.path.join(tmp, "draw")], check=True, capture_output=True,
                                 text=True, env={**os.environ, "OMP_NUM_THREADS": "1"}
                                 ).stdout.split()[-1])
        cc = os.path.join(tmp, "draw.cc")
        with open(cc, "w") as f:
            f.write('#include "UCMuGen.h"\n#include <cmath>\n#include <cstdio>\nint main() {\n'
                    "  ucmugen::UCMuonLcg g(4321);\n  int k = 0;\n"
                    f"  for (int i = 0; i < {n}; ++i) {{\n    const double r = g.flat();\n"
                    "    if (r <= 0.0 || r >= 1.0) return 2;\n"
                    "    if (std::floor(r * 16777216.0) == r * 16777216.0) ++k;\n  }\n"
                    '  std::printf("%d\\n", k);\n}\n')
        subprocess.run(["c++", "-std=c++17", "-O2", f"-I{ROOT / 'ucmugen/include'}", cc,
                        "-o", os.path.join(tmp, "drawcc")], check=True, capture_output=True)
        k_c = int(subprocess.run([os.path.join(tmp, "drawcc")], check=True, capture_output=True,
                                 text=True).stdout)
    return n, k_f, k_c


def main() -> int:
    print("Random deviates")
    n, k_f, k_c = deviates_on_24_bit_grid()
    check(k_f < 10, f"Fortran par_ranlux: {k_f} of {n} draws on the 2^-24 grid (all of them up to 1.3.2)")
    check(k_c < 10, f"UCMuGen UCMuonLcg: {k_c} of {n} draws on the 2^-24 grid (all of them up to 1.3.2)")

    print("PARMA W from the date")
    w_in, out_in = parma_w(2020, 1, 1)
    check(w_in > 0 and "WARNING" not in out_in, f"2020-01-01 inside the table: W = {w_in}, no warning")
    w_late, out_late = parma_w(2027, 3, 1)
    check(w_late > 0 and "no PARMA W index for this date; using the latest one" in out_late,
          f"2027-03-01 after the table: W = {w_late} (latest available), with a warning")
    w_far, out_far = parma_w(2050, 12, 31)
    check(w_far == w_late and "using the latest one" in out_far,
          f"2050-12-31 (latest GUI date): W = {w_far}, same fallback")

    print("MPI generator")
    if not MPI_BIN.exists() or shutil.which("mpirun") is None:
        print("  skipped: needs bin/ucmuon_gen (make ucmuon_gen) and mpirun")
    else:
        cfg = Config(spectrum=4, emin=0.5, emax=100.0, theta_max_deg=60.0, radius_m=10.0,
                     n=3000, seed=7)
        omp = run_fortran(cfg, threads=1)
        rows, out = run_mpi(cfg)
        r_mpi = float(re.search(r"Surface rate R\s*:\s*([0-9.Ee+\-]+)", out).group(1))
        tried = int(re.search(r"Total tried \(all ranks\):\s*(\d+)", out).group(1))
        check("Min Energy:" in out and re.search(r"Min Energy:\s*0\.50 GeV", out) is not None
              and rows[:, 10].min() < 1.0, f"E_min 0.5 GeV is kept (lowest E {rows[:, 10].min():.3f} GeV)")
        check(abs(r_mpi / omp.rate - 1.0) < 1e-6, f"surface rate R equals the OpenMP one ({r_mpi:.6e})")
        check(len(rows) == cfg.n and tried == cfg.n,
              f"no detector: {len(rows)} rows and Tried = {tried} for {cfg.n} requested")
        check(len(np.unique(rows[:, 0])) == len(rows),
              f"no detector: {len(np.unique(rows[:, 0]))} distinct EventIDs in {len(rows)} rows")
        cfg_d = Config(spectrum=4, emin=1.0, emax=100.0, theta_max_deg=60.0, radius_m=30.0,
                       n=200, seed=7, detector=Detector("cyl", (0.0, 0.0, -2400.0, 0.0, 0.0,
                                                                -2000.0, 150.0), margin=50.0))
        surf, sel, _ = run_mpi(cfg_d, nranks=3)              # 200 = 68 + 66 + 66
        check(np.array_equal(sel[:, 0], np.arange(1, cfg_d.n + 1)),
              f"detector, 3 ranks: selected EventIDs run 1..{cfg_d.n} "
              f"({len(np.unique(sel[:, 0]))} distinct in {len(sel)} rows)")
        check(len(np.unique(surf[:, 0])) == len(surf),
              f"detector, 3 ranks: {len(np.unique(surf[:, 0]))} distinct surface EventIDs "
              f"in {len(surf)} rows")
        cfg_e = Config(spectrum=8, emin=0.5, emax=10.0, theta_max_deg=60.0, radius_m=10.0,
                       n=2000, seed=7)
        rows_e, out_e = run_mpi(cfg_e)
        mass = np.sqrt(np.maximum(rows_e[:, 10] ** 2 - rows_e[:, 4] ** 2, 0.0))
        check("Cosmic electrons" in out_e and mass.max() < 0.01,
              f"spectrum 8 gives electrons (largest sqrt(E^2 - p^2) = {mass.max():.4f} GeV)")

    print("ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
