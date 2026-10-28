#!/usr/bin/env python3
"""The paper's underground rate example: the muon rate in a cylindrical
detector (radius 10 cm, length 50 cm) with its top face 90 m below the
surface, in Standard Rock.

    python3 examples/underground_rate/ucmuon_underground_rate.py [--out DIR]

runs, from the repository root after `make`:
  1. bin/ucmuon_gen_omp: the Guan spectrum with the joint flux (angular
     mode 6), 45-2500 GeV, theta <= 70 deg, from a disk of the minimum radius
     (docs/GENERATOR_INPUT.md), keeping 4e5 muons aimed at the detector
     inflated by the suggested safety margin, 2.5 sigma_r (gui/mcs_margin.py);
  2. bin/ucmuon_transport_music_omp through 90 m of Standard Rock;
  3. the GUI's own ray test (apply_ug_detector_filter in gui/ucmuon_gui.py)
     with margin 0, to count the muons that cross the detector itself;
and prints R = hits / T, with T = N_tried / R_S the live time the generator
prints.

With one OpenMP thread (the default) and the default seed the run reproduces
the paper's numbers exactly (about 20 min on an Apple M4 core); with more
threads (--threads) it gives a statistically equivalent sample.
"""
from __future__ import annotations

import argparse
import ast
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "gui"))
import mcs_margin  # noqa: E402

DEPTH_M = 90.0                       # top face of the detector
DET = dict(shape=1, ax=0.0, ay=0.0, az=-9050.0, bx=0.0, by=0.0, bz=-9000.0, r=10.0)   # cm
RHO, X0 = 2.65, 26.48                # Standard Rock [g/cm3, g/cm2]
E_MIN, E_MAX, THETA_MAX, N_KEEP = 45.0, 2500.0, 70.0, 400_000


def gui_detector_filter():
    """apply_ug_detector_filter and its ray tests, taken from the GUI source so
    that Streamlit is not imported."""
    tree = ast.parse((REPO / "gui" / "ucmuon_gui.py").read_text())
    ns = {"np": np}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in (
                "_ray_hits_cylinder", "_ray_hits_box", "apply_ug_detector_filter"):
            exec(compile(ast.Module([node], []), "ucmuon_gui.py", "exec"), ns)
    return ns["apply_ug_detector_filter"]


def generator_input(out: Path, margin_cm: float, r_min_m: float) -> str:
    """Answers to bin/ucmuon_gen_omp's prompts (docs/GENERATOR_INPUT.md)."""
    a, b = DET["az"], DET["bz"]
    return "\n".join([
        "0",                                  # read the parameters below
        f"{E_MIN}", f"{E_MAX}",               # energy window [GeV]
        "4",                                  # Guan et al. (2015)
        "1", "1",                             # disk in the XY plane
        "0.0", "0.0", f"{r_min_m:.1f}", "0.0", "0.0", "0.0",   # centre, radius, w, no tilt
        "6", f"{THETA_MAX}",                  # joint J(p, theta) through the surface
        f"{N_KEEP}",
        "1", "1",                             # one detector
        "1", f"{margin_cm:.1f}",              # cylinder, safety margin [cm]
        f"0.0 0.0 {a}", f"0.0 0.0 {b}", f"{DET['r']}",
        "0", "0",                             # selected muons only, no PHITS file
        str(out / "muons_surface.dat"), str(out / "muons_selected.dat"), "", ""])


def music_input(out: Path) -> str:
    """Answers to bin/ucmuon_transport_music_omp's prompts."""
    return "\n".join([
        str(out / "muons_selected.dat"), str(out / "muons_underground.dat"),
        f"{RHO}", f"{X0}", f"{DEPTH_M}",
        "1", "1",                             # 3-D transport, other-process scattering
        "-30",                                # minv
        "1",                                  # cross-section tables on disk
        "1",                                  # rock
        "", ""])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=Path, default=REPO / "output" / "underground_rate")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--count-only", action="store_true",
                   help="only count the hits of an earlier run in --out")
    args = p.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    # The margin is evaluated at the bottom face, the deepest point of the detector.
    bottom_m = -DET["az"] / 100
    margin_cm = mcs_margin.suggested_margin(bottom_m, RHO, E_MIN)["margin_cm"]
    m = margin_cm / 100
    r_min_m = DET["r"] / 100 + m + (bottom_m + m) * math.tan(math.radians(THETA_MAX))
    print(f"safety margin 2.5 sigma_r = {margin_cm:.1f} cm, source disk radius {r_min_m:.1f} m")

    env = dict(os.environ, UCMUON_SEED=str(args.seed), OMP_NUM_THREADS=str(args.threads))
    if not args.count_only:
        (out / "input_gen.txt").write_text(generator_input(out, margin_cm, r_min_m))
        with open(out / "input_gen.txt") as fi, open(out / "gen.log", "w") as fo:
            subprocess.run([str(REPO / "bin" / "ucmuon_gen_omp")], stdin=fi, stdout=fo,
                           env=env, check=True)
        with open(out / "music.log", "w") as fo:   # MUSIC reads its tables from bin/
            subprocess.run([str(REPO / "bin" / "ucmuon_transport_music_omp")],
                           input=music_input(out), text=True, stdout=fo, cwd=REPO / "bin",
                           env=env, check=True)

    log = (out / "gen.log").read_text()
    n_tried = int(re.search(r"Tried:\s+(\d+)", log).group(1))
    rate_s = float(re.search(r"Surface rate R:\s+(\S+)", log).group(1))
    live = float(re.search(r"Live time:\s+(\S+)", log).group(1))
    cols = ("EventID xs ys zs Es theta_s phi_s charge alive "
            "x y z E cx cy cz theta phi").split()
    ug = pd.read_csv(out / "muons_underground.dat", sep=r"\s+", comment="#", header=None,
                     names=cols)
    _, hits = gui_detector_filter()(ug, [dict(DET, margin=0.0)])
    print(f"tried {n_tried:.3e}, R_S = {rate_s:.3e} /s, T = {live:.0f} s")
    print(f"muons {len(ug)}, survivors {int((ug.alive == 1).sum())}, in the detector {hits}")
    print(f"R = {hits}/T = {60 * hits / live:.2f} +- {60 * math.sqrt(hits) / live:.2f} muons/min")


if __name__ == "__main__":
    main()
