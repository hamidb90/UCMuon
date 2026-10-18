#!/usr/bin/env python3
"""Random number generators: src/common/ucmuon_rng.f90 against CERNLIB.

src/common/ucmuon_rng.f90 is an independent (MIT) implementation of RANLUX,
RANMAR and RNORML that replaces the GPL CERNLIB routines UCMuon used up to
1.4.1. It must give the same numbers, bit for bit, so that seeded runs of
the transport engines are unchanged. tests/rng/rng_dump.f90 prints the
streams (unseeded first calls, every luxury level, several seeds, skips,
vector lengths, RNORML in pairs, singles and sevens, four OpenMP threads)
and tests/rng/reference_vectors.txt holds its output with the CERNLIB
routines. Needs gfortran with OpenMP; skipped when gfortran is absent.

    python3 -m pytest tests/rng -q        (or python3 tests/rng/test_rng.py)
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "common" / "ucmuon_rng.f90"
DUMP = ROOT / "tests" / "rng" / "rng_dump.f90"
REFERENCE = ROOT / "tests" / "rng" / "reference_vectors.txt"
MODES = ("ldef", "mdef", "ndef", "main", "omp")     # the order of the reference file
DATA = re.compile(r"^(\S+) (\d+) (-?\d+)$")         # case, position, IEEE bit pattern


def data_lines(text: str) -> list[tuple[str, int, int]]:
    """The "case position bits" lines of a dump; comments and messages are dropped."""
    out = []
    for line in text.splitlines():
        m = DATA.match(line.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out


def dump_streams() -> list[tuple[str, int, int]]:
    """Build rng_dump with ucmuon_rng.f90 and run every mode, each in its own process."""
    with tempfile.TemporaryDirectory() as tmp:
        exe = os.path.join(tmp, "rng_dump")
        build = subprocess.run(["gfortran", "-O2", "-fopenmp", "-J", tmp, str(SRC), str(DUMP),
                                "-o", exe], capture_output=True, text=True, cwd=tmp)
        assert build.returncode == 0, f"gfortran failed:\n{build.stderr[-3000:]}"
        lines = []
        for mode in MODES:
            env = {**os.environ, "OMP_NUM_THREADS": "4" if mode == "omp" else "1"}
            run = subprocess.run([exe, mode], capture_output=True, text=True, env=env, cwd=tmp,
                                 timeout=300)
            assert run.returncode == 0, f"rng_dump {mode} failed:\n{run.stdout[-2000:]}{run.stderr[-2000:]}"
            lines += data_lines(run.stdout)
    return lines


def test_streams_match_cernlib_reference():
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not found")
    want = data_lines(REFERENCE.read_text())
    got = dump_streams()
    for i, (w, g) in enumerate(zip(want, got)):
        if w != g:
            pytest.fail(f"first mismatch at data line {i + 1}: case {w[0]}, position {w[1]}: "
                        f"expected bits {w[2]}, got {g[0]} {g[1]} {g[2]}")
    assert len(got) == len(want), f"{len(got)} data lines, reference has {len(want)}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
