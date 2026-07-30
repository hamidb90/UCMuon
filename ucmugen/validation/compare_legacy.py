#!/usr/bin/env python3
"""Compare UCMuGen's legacy driver against the Fortran generator, event by event.

Because UCMuGen reproduces the Fortran RNG bit-for-bit and consumes draws in
the same order, the two implementations must agree *exactly*, not merely
statistically. That is a far stronger check than any distributional test, so it
is the primary validation: a single differing event is a failure.

The statistical machinery in ucmuref.stats is still needed, but for the cases
where exact agreement is not expected by design (see README, "Known deliberate
divergence") and for the projection-correct API.

Usage:
    python3 compare_legacy.py [--n N] [--verbose]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from ucmuref import GenConfig, parse, run
from ucmuref.cases import build, label

HERE = Path(__file__).resolve().parent
DUMP = HERE / "dump_ucmugen"

# Columns that must agree exactly. event/hit_flag/det_mask are bookkeeping the
# legacy driver does not model.
COMPARED = ("x_cm", "y_cm", "z_cm", "p_GeV", "px_GeV", "py_GeV", "pz_GeV",
            "theta_rad", "phi_rad", "E_GeV", "charge")

# Configurations UCMuGen deliberately does not reproduce. The Fortran applies
# its tilt transform to positions only and permutes direction components for
# vertical planes, neither of which carries a projection weight; UCMuGen's
# surfaces do. See README, "Known deliberate divergence".
def is_divergent(cfg: GenConfig) -> str | None:
    if cfg.tilt_deg:
        return "tilted surface: Fortran applies no projection weighting"
    if cfg.source_mode in (1, 2) and cfg.source_plane != 1:
        return "vertical plane: Fortran rigidly rotates the sky"
    return None


def cpp_args(cfg: GenConfig) -> list[str]:
    """Map a GenConfig onto the dump_ucmugen command line."""
    if cfg.source_mode == 3:
        radius_cm = cfg.hemi_radius_m * 100.0
        source_z_cm = cfg.hemi_cz_m * 100.0
    else:
        radius_cm = cfg.disk_r_m * 100.0
        source_z_cm = cfg.w_m * 100.0
    half_lx_cm = (cfg.rect_u2_m - cfg.rect_u1_m) / 2.0 * 100.0
    half_ly_cm = (cfg.rect_v2_m - cfg.rect_v1_m) / 2.0 * 100.0
    return [str(DUMP),
            repr(cfg.emin), repr(cfg.emax),
            str(cfg.spectrum), str(cfg.angular_mode), repr(cfg.theta_max_deg),
            str(cfg.source_mode),
            repr(radius_cm), repr(half_lx_cm), repr(half_ly_cm),
            repr(source_z_cm),
            str(cfg.nmuons), str(cfg.seed)]


def run_cpp(cfg: GenConfig) -> np.ndarray:
    proc = subprocess.run(cpp_args(cfg), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"dump_ucmugen failed: {proc.stderr.strip()}")
    fd, path = tempfile.mkstemp(suffix=".dat")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(proc.stdout)
        return parse(Path(path))
    finally:
        os.unlink(path)


def compare_one(cfg: GenConfig, verbose: bool = False) -> tuple[bool, str]:
    fortran = run(cfg)
    cpp = run_cpp(cfg)
    if len(fortran) != len(cpp):
        return False, f"length {len(fortran)} vs {len(cpp)}"

    bad = []
    for col in COMPARED:
        if not np.array_equal(fortran[col], cpp[col]):
            d = np.abs(fortran[col].astype(float) - cpp[col].astype(float))
            bad.append(f"{col}: {int((d > 0).sum())}/{len(d)} differ, "
                       f"max|d|={d.max():.3e}")
    if bad:
        return False, "; ".join(bad)
    return True, f"{len(cpp)} events identical"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20_000,
                    help="events per configuration (default 20000)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not DUMP.exists():
        print(f"error: {DUMP.name} not built. Run:\n"
              f"  c++ -std=c++17 -O2 -I../include dump_ucmugen.cc -o dump_ucmugen",
              file=sys.stderr)
        return 2

    cases = build(n=args.n)
    n_pass = n_fail = n_skip = 0

    for i, cfg in enumerate(cases, 1):
        reason = is_divergent(cfg)
        if reason:
            n_skip += 1
            print(f"{i:3d}. SKIP  {label(cfg)[:66]:<66s} {reason}")
            continue
        try:
            ok, detail = compare_one(cfg, args.verbose)
        except Exception as exc:                      # noqa: BLE001
            ok, detail = False, f"{type(exc).__name__}: {exc}"
        if ok:
            n_pass += 1
            print(f"{i:3d}. ok    {label(cfg)[:66]:<66s} {detail}")
        else:
            n_fail += 1
            print(f"{i:3d}. FAIL  {label(cfg)[:66]:<66s} {detail}")

    print()
    print(f"exact match: {n_pass} passed, {n_fail} failed, "
          f"{n_skip} skipped by design ({args.n:,} events each)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
