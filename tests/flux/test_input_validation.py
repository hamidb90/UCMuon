#!/usr/bin/env python3
"""The surface generator rejects invalid input (src/generator/input_module.f90).

Starting from the test run's input (test_run/input_gen.dat), each case changes
one value or cuts the input short; the generator must stop with exit status 2
and a message naming the parameter, and the unchanged input must still run.
Up to UCMuon 1.4.1 every one of these cases ran with status 0: angular mode 9
became 2, e_min > e_max wrote every muon at e_min, nmuons = -5 wrote an empty
file and a truncated input ended in a gfortran runtime error.

    python3 tests/flux/test_input_validation.py [path/to/ucmuon_gen_omp]
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "bin" / "ucmuon_gen_omp"


def base_input(out: Path) -> list[str]:
    lines = []
    for line in (ROOT / "test_run" / "input_gen.dat").read_text().splitlines():
        v = line.split("#")[0].strip()
        if v:
            lines.append(v.replace("test_run/output/", f"{out}/"))
    return lines


# (description, value index in the stripped input, new value or None to cut
#  the input there, text the message must contain)
CASES = [
    ("angular mode 9",          12, "9",     "angular_mode must be 1 to 6"),
    ("e_max below e_min",        2, "4.0",   "e_max (4.00000 GeV) must be above e_min"),
    ("negative muon count",     14, "-5",    "nmuons must be at least 1"),
    ("input cut before use_detector", 15, None, "the input ends before use_detector"),
    ("spectrum mode 11",         3, "11",    "spectrum_mode must be 1 to 8"),
    ("negative disk radius",     8, "-3",    "the disk radius must be positive"),
    ("theta_max not a number",  13, "abc",   "theta_max is not a number"),
    ("theta_max above 90",      13, "95",    "theta_max must be above 0 and at most 90"),
    ("negative safety margin",  18, "-1",    "the safety margin cannot be negative"),
    ("zero cylinder radius",    21, "0",     "the cylinder radius must be positive"),
    ("nothing to write",        15, "0",     "nothing to write"),
]


def run(lines: list[str], tmp: Path) -> tuple[int, str]:
    env = {**os.environ, "UCMUON_SEED": "1", "OMP_NUM_THREADS": "1"}
    p = subprocess.run([str(GEN)], input="\n".join(lines) + "\n", text=True, env=env,
                       capture_output=True, cwd=tmp)
    return p.returncode, p.stderr


def main() -> int:
    if not GEN.exists():
        print(f"[skip] {GEN} not built (make ucmuon_gen_omp)")
        return 0
    failed = []
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        base = base_input(tmp)
        for name, i, value, expect in CASES:
            lines = base[:i] if value is None else base[:i] + [value] + base[i + 1:]
            if name == "nothing to write":          # no detector and save_all = 0
                lines = base[:15] + ["0", "0", "0", base[-2]]
            code, err = run(lines, tmp)
            ok = code == 2 and "INPUT ERROR" in err and expect in err
            print(f"  [{'ok  ' if ok else 'FAIL'}] {name:30s} exit {code}: {err.strip().splitlines()[0] if err.strip() else '(no message)'}")
            if not ok:
                failed.append(name)
        base_small = list(base)
        base_small[14] = "200"                      # the valid input, smaller
        code, err = run(base_small, tmp)
        n = sum(1 for line in open(tmp / "muons_selected.dat") if not line.startswith("#")) if code == 0 else 0
        ok = code == 0 and n == 200
        print(f"  [{'ok  ' if ok else 'FAIL'}] {'valid input':30s} exit {code}, {n} muons selected")
        if not ok:
            failed.append("valid input")
    print("\nALL PASS" if not failed else f"\n{len(failed)} FAILED: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
