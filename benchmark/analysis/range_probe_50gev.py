#!/usr/bin/env python3
"""The end of the 50 GeV beam's range, vertical and at 85 degrees, in Geant4 and
PHITS: the fraction of the beam that reaches the slant lengths L = 70 ... 100 m
(planes at the depth L cos(theta)).

At 85 degrees Geant4 passes 4.6 % of the 50 GeV beam at L = 100 m and PHITS
about none (inclined_summary.py). This probe shows why: both codes spread the
end of the range at 85 degrees the same way, but PHITS's ranges are shorter.

    python3 benchmark/analysis/range_probe_50gev.py inputs RUN_DIR
        writes RUN_DIR/probe_th{0,85}.inp and src_th{0,85}.dat: the first
        20 000 muons of the 50 GeV beam and the PHITS deck with its six planes
        moved to L = 70, 76, 82, 88, 94, 100 m and a tally that dumps every
        muon crossing with its history number. Run each with
        PHITSPATH=... bash benchmark/codes/phits/run_phits_chunks.sh \\
            RUN_DIR/probe_th0.inp RUN_DIR/src_th0.dat RUN_DIR/run_th0 1 1
        (and the same for th85).
    python3 benchmark/analysis/range_probe_50gev.py summary RUN_DIR
        writes benchmark/results/range_probe_50GeV_20261027.csv from
        RUN_DIR/run_th{0,85}/ (PHITS: the muons that cross each plane, from
        the dump, and its crossings, which count a muon that crosses a plane
        twice twice, as PHITS's own tallies do) and Geant4's stopped muons
        (<G4_OUT>/run_stdrock_xy1km_20261022_mono_stopped.csv and
        ..._20261027_incl85_stopped.csv; a muon reaches L if it stops
        below the plane at L or is never stopped).
"""
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SIM = Path.home() / "simulations" / "geant4" / "geant4_muon_rock_v5" / "sources"
G4_OUT = Path.home() / "simulations" / "geant4" / "geant4_muon_rock_v5" / "outputs"
INCL = REPO / "benchmark" / "geant4_muon_rock_v5" / "inclined"
OUT = REPO / "benchmark" / "results" / "range_probe_50GeV_20261027.csv"
L_M = (70, 76, 82, 88, 94, 100)
N_PHITS = 20_000
N_G4 = 100_000
DUMP = """[ T-Cross ]
  title  = every muon crossing, dumped (muons that cross a plane more than once)
  mesh   = reg
  reg    = 6
  r-from  r-to   area
  1  2  25000000.
  2  3  25000000.
  3  4  25000000.
  4  5  25000000.
  5  6  25000000.
  6  7  25000000.
  part   = muon+  muon-
  e-type = 1
  ne     = 1
  0.  1000000.
  output = f-curr
  axis   = reg
  file   = cross_dump.out
  unit   = 1
  dump   = -6
           4  7  8  9  18  19

"""
CASES = {0: (REPO / "benchmark" / "codes" / "phits" / "muon_rock.inp", SIM / "benchmark_phits.dat",
             "run_stdrock_xy1km_20261022_mono"),
         85: (INCL / "muon_rock_th85.inp", INCL / "benchmark_phits_th85.dat",
              "run_stdrock_xy1km_20261027_incl85")}


def inputs(run: Path) -> None:
    run.mkdir(parents=True, exist_ok=True)
    for th, (deck, src, _) in CASES.items():
        c = math.cos(math.radians(th))
        recs = []
        with open(src) as fh:
            for line in fh:
                t = line.split()
                if len(t) == 10 and 45e3 < float(t[7].replace("D", "E")) < 55e3:   # ek [MeV]
                    recs.append(line)
                    if len(recs) == N_PHITS:
                        break
        assert len(recs) == N_PHITS, (src, len(recs))
        (run / f"src_th{th}.dat").write_text("".join(recs))
        out, k = [], 0
        for line in deck.read_text().splitlines(keepends=True):
            m = re.match(r"^(\s+)(11|12|13|14|15|6)(\s+pz\s+)(-?\d+\.\d*)(.*)$", line.rstrip("\n"))
            if m:
                line = f"{m.group(1)}{m.group(2)}{m.group(3)}{-L_M[k] * 100 * c:.4f}      $ slant L = {L_M[k]} m\n"
                k += 1
            out.append(line)
        assert k == len(L_M), k
        deck_text = "".join(out)
        i = deck_text.rindex("[ End ]")
        (run / f"probe_th{th}.inp").write_text(deck_text[:i] + DUMP + deck_text[i:])
        print(f"[ok] {run}/probe_th{th}.inp, src_th{th}.dat ({N_PHITS} muons)")


def summary(run: Path) -> None:
    rows = []
    for th, (_, _, g4) in CASES.items():
        s = pd.read_csv(G4_OUT / f"{g4}_stopped.csv", usecols=["InitKEGeV", "StopDepthCm"])
        d = s[(s.InitKEGeV > 45) & (s.InitKEGeV < 55)].StopDepthCm.to_numpy() / 100
        for L in L_M:
            k = N_G4 - len(d) + int((d >= L).sum())
            rows.append(("Geant4", th, L, N_G4, k, ""))
        ph = list(csv.DictReader(open(run / f"run_th{th}" / "phits_summary.csv")))
        assert len(ph) == len(L_M)
        dumps = sorted((run / f"run_th{th}").glob("c*/cross_dump_dmp.out"))
        assert dumps, f"no crossing dump in {run}/run_th{th}"
        z, key = [], []
        for c, f in enumerate(dumps):     # z w e wt nocas nobch; a muon is (chunk, batch, history)
            a = np.loadtxt(f, dtype=str, ndmin=2)
            a = np.char.replace(a, "D", "E").astype(float)
            z.append(a[:, 0]); key.append(c * 1e12 + a[:, 5] * 1e7 + a[:, 4])
        z, key = np.concatenate(z), np.concatenate(key)
        c = math.cos(math.radians(th))
        for L, r in zip(L_M, ph):
            m = np.abs(z + L * 100 * c) < 0.05
            assert m.sum() == round(float(r["N_transmitted"])), (th, L, m.sum(), r["N_transmitted"])
            rows.append(("PHITS", th, L, N_PHITS, len(np.unique(key[m])), int(m.sum())))
    with open(OUT, "w", newline="") as fh:
        fh.write("# 50 GeV beam: fraction reaching the slant length L, vertical and at 85 deg.\n"
                 "# n_pass counts muons (PHITS: distinct histories in its crossing dump),\n"
                 "# n_crossings PHITS's crossings (benchmark/analysis/range_probe_50gev.py)\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["code", "theta_deg", "L_m", "n_source", "n_pass", "n_crossings", "pass_pct",
                    "pass_err_pct"])
        for code, th, L, n, k, nc in rows:
            p = k / n
            w.writerow([code, th, L, n, k, nc, f"{100 * p:.2f}", f"{100 * math.sqrt(p * (1 - p) / n):.2f}"])
    for code in ("Geant4", "PHITS"):
        for th in CASES:
            print(f"{code:6s} {th:2d} deg: " + "  ".join(f"{L} m {100 * k / n:5.1f}"
                                                        for c, t, L, n, k, _ in rows if c == code and t == th))
    twice = [(nc - k) / k for c, t, L, n, k, nc in rows if c == "PHITS" and k]
    print(f"PHITS: at most {100 * max(twice):.2f} % more crossings than muons at any plane")
    print(f"[ok] {OUT}")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("inputs", "summary"):
        sys.exit(__doc__)
    (inputs if sys.argv[1] == "inputs" else summary)(Path(sys.argv[2]))
