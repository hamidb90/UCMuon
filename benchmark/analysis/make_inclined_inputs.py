#!/usr/bin/env python3
"""Inputs of the inclined-beam benchmark: the six monoenergetic beams of the
vertical benchmark (5 to 300 GeV, 10^5 muons each, the same muons) tilted to a
zenith angle theta, travelling towards +x, with the scoring planes at the
slant lengths L = 1, 10, 25, 50, 100 and 200 m, that is at the vertical depths
L cos(theta).

    python3 benchmark/analysis/make_inclined_inputs.py [--theta 60,80,85]

For every angle it writes, in OUT (default benchmark/geant4_muon_rock_v5/inclined):
  benchmark_surface_th<theta>.dat  the engines' source (run_mono_bench.sh with
                                   SOURCE_FILE and INCL_DEG = theta)
  benchmark_geant4_th<theta>.txt   Geant4's (run_parallel.py --incl theta --xy 100000)
  benchmark_phits_th<theta>.dat    PHITS's dump (run_phits_mpi.sh)
  muon_rock_th<theta>.inp          the PHITS deck, benchmark/codes/phits/muon_rock.inp
                                   with its planes at L cos(theta)
Every code then labels its outputs by the slant length L, so they compare
directly with the vertical benchmark at the same L. A muon goes at most
L sin(theta) < 200 m sideways, well inside the slabs of +-1 km.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SIM = Path.home() / "simulations" / "geant4" / "geant4_muon_rock_v5" / "sources"
ENGINE_SRC = REPO / "benchmark" / "geant4_muon_rock_v5" / "sources" / "benchmark_surface.dat"
G4_SRC = SIM / "benchmark_geant4_new.txt"
PHITS_SRC = SIM / "benchmark_phits.dat"
PHITS_DECK = REPO / "benchmark" / "codes" / "phits" / "muon_rock.inp"
PLANES_CM = (100.0, 1000.0, 2500.0, 5000.0, 10000.0, 20000.0)


def check_vertical(p: float, pz: float, where: str) -> None:
    if p <= 0 or abs(pz) / p < 1.0 - 1e-9 or pz >= 0:
        sys.exit(f"{where}: not a vertical downward muon (p = {p}, pz = {pz})")


def tilt(lines, theta: float, fmt_line, tag: str):
    """Yield the source lines with every muon tilted to theta (comments kept)."""
    yield from (f"# Inclined to theta = {math.degrees(theta):g} deg towards +x by "
                f"benchmark/analysis/make_inclined_inputs.py ({tag})\n",)
    for line in lines:
        yield line if line.startswith("#") or not line.strip() else fmt_line(line.split(), theta)


def engine_line(t: list[str], th: float) -> str:
    # EventID x y z p px py pz theta phi E charge hit_flag det_mask
    p, pz = float(t[4]), float(t[7])
    check_vertical(p, pz, f"engine source event {t[0]}")
    t[5], t[6], t[7] = f"{p * math.sin(th):.6f}", "0.000000", f"{-p * math.cos(th):.6f}"
    t[8], t[9] = f"{th:.9f}", "0.000000000"
    return "  ".join(t) + "\n"


def g4_line(t: list[str], th: float) -> str:
    # PDG x[mm] y[mm] z[mm] px py pz [MeV/c] Ekin[MeV]
    px, py, pz = (float(v) for v in t[4:7])
    p = math.sqrt(px * px + py * py + pz * pz)
    check_vertical(p, pz, "Geant4 source")
    return (f"{t[0]:>4s} {t[1]:>15s} {t[2]:>15s} {t[3]:>14s} {p * math.sin(th):15.6f} "
            f"{0.0:15.6f} {-p * math.cos(th):15.6f} {t[7]:>15s}\n")


def d(v: float) -> str:
    return f"{v:.15E}".replace("E", "D")


def phits_line(t: list[str], th: float) -> str:
    # kf x y z u v w ek wgt time (dump = -10); w < 0 is downward
    u, v, w = (float(x.replace("D", "E")) for x in t[4:7])
    check_vertical(1.0, w, "PHITS dump")
    t[4], t[5], t[6] = d(math.sin(th)), d(0.0), d(-math.cos(th))
    return "  ".join(t) + "\n"


def phits_deck(text: str, theta_deg: float, dump: str) -> str:
    c = math.cos(math.radians(theta_deg))
    out, n_planes, n_mesh = [], 0, 0
    for line in text.splitlines(keepends=True):
        m = re.match(r"^(\s+)(11|12|13|14|15|6)(\s+pz\s+)(-?\d+\.\d*)(.*)$", line.rstrip("\n"))
        if m:
            z = float(m.group(4))
            assert -z in PLANES_CM, line
            line = (f"{m.group(1)}{m.group(2)}{m.group(3)}{z * c:.4f}"
                    f"      $ slant L = {-z / 100:g} m (vertical {-z * c / 100:.4f} m)\n")
            n_planes += 1
        elif line.strip() == "-20000.  -10000.  -5000.  -2500.  -1000.  -100.  0.":
            line = "  " + "  ".join(f"{-v * c:.4f}" for v in reversed(PLANES_CM)) + "  0.\n"
            n_mesh += 1
        elif line.strip() in ("-20000.  -19999.", "-100.  -99."):    # 1 cm x-y map layers
            z = float(line.split()[0])
            line = f"  {z * c:.4f}  {z * c + 1.0:.4f}\n"
            n_mesh += 1
        elif re.match(r"^\s*file\s*=\s*benchmark_phits\.dat", line):
            line = line.replace("benchmark_phits.dat", dump)
        elif line.startswith(" Cosmic Muon in Standard Rock"):
            line = line.rstrip("\n") + f", beams at {theta_deg:g} deg\n"
        out.append(line)
    assert n_planes == 6 and n_mesh == 3, (n_planes, n_mesh)
    head = (f"$  INCLINED VARIANT (make_inclined_inputs.py): beams at theta = {theta_deg:g} deg,\n"
            f"$  scoring planes at L cos(theta) for L = 1 ... 200 m; every result is labelled\n"
            f"$  by the slant length L. Otherwise identical to muon_rock.inp.\n")
    return head + "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--theta", default="60,80,85", help="zenith angles in degrees")
    ap.add_argument("--out", type=Path, default=REPO / "benchmark" / "geant4_muon_rock_v5" / "inclined")
    ap.add_argument("--engine-source", type=Path, default=ENGINE_SRC)
    ap.add_argument("--g4-source", type=Path, default=G4_SRC)
    ap.add_argument("--phits-source", type=Path, default=PHITS_SRC)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    for th_deg in (float(s) for s in a.theta.split(",")):
        th, tag = math.radians(th_deg), f"th{th_deg:g}"
        for src, name, fmt in ((a.engine_source, f"benchmark_surface_{tag}.dat", engine_line),
                               (a.g4_source, f"benchmark_geant4_{tag}.txt", g4_line),
                               (a.phits_source, f"benchmark_phits_{tag}.dat", phits_line)):
            with open(src) as fin, open(a.out / name, "w") as fout:
                lines = tilt(fin, th, fmt, src.name)
                if fmt is phits_line:       # PHITS reads every line of a dump: no comments
                    next(lines)
                fout.writelines(lines)
        deck = a.out / f"muon_rock_{tag}.inp"
        deck.write_text(phits_deck(PHITS_DECK.read_text(), th_deg, f"benchmark_phits_{tag}.dat"))
        print(f"[ok] theta = {th_deg:g} deg: {a.out}/*_{tag}.*")
    return 0


if __name__ == "__main__":
    sys.exit(main())
