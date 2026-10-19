#!/usr/bin/env python3
"""Generate the depth-resolved survival-fraction table for the UCMuon paper.

Reads the v2 benchmark outputs (six monoenergetic vertical beams at
5/10/20/50/100/300 GeV, 1e5 muons each, identical source population fed to
every engine) and emits a LaTeX table of survival fraction (%) per
(surface energy, depth) for the four UCMuon transport engines.

Authoritative data: benchmark/geant4_muon_rock_v5/<ENGINE>/<ENGINE>_bench_<DEPTH>m.dat
(the post-fix "v2" run of 2026-05-25 documented in BENCHMARK_FEEDBACK.md; the
flat ../output/*_bench_*.dat copies are an earlier v1 UCMuon run and must NOT
be used). The 10 m plane, excluded up to UCMuon 1.4.1 because the Geant4
reference rock was wrong there (see benchmark_analysis.py), is included.

Column 5 (1-based) = surface energy [GeV]; column 9 = alive flag (0/1).
The .dat files are git-ignored (too large to track); regenerate with
benchmark/analysis/run_mono_bench.sh if absent.

Usage:  python make_survival_table.py  >  tab_survival_matrix.tex
"""
import sys
from pathlib import Path
from collections import defaultdict

BENCH    = Path(__file__).resolve().parent.parent.parent / "benchmark" / "geant4_muon_rock_v5"
ENGINES  = ["BB", "UCMuon", "MUSIC", "PROPOSAL"]      # column order
LABELS   = {"BB": "Bethe--Bloch", "UCMuon": "UCMuon-MC",
            "MUSIC": r"\music{}", "PROPOSAL": r"\proposal{}"}
DEPTHS   = [1, 10, 25, 50, 100, 200]
# Geant4 per beam and depth (Standard Rock reference), shipped as a small CSV
# because the per-crossing Geant4 output (400 MB) is not.
G4_CSV   = Path(__file__).resolve().parent.parent.parent / "benchmark" / "results" / "geant4_survival_by_beam.csv"


def geant4_survival():
    """Return {(E_GeV, depth_m): survival_percent} from G4_CSV."""
    out = {}
    with open(G4_CSV) as fh:
        for line in fh:
            if line.startswith(("#", "beam")) or not line.strip():
                continue
            b, d, s = line.strip().split(",")
            out[(float(b), int(d))] = float(s)
    return out


def survival(engine, depth):
    """Return {E_GeV: survival_percent} for one engine/depth file, or None."""
    f = BENCH / engine / f"{engine}_bench_{depth}m.dat"
    if not f.exists():
        return None
    n = defaultdict(int)
    a = defaultdict(int)
    with open(f) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            e = float(c[4])           # surface energy
            n[e] += 1
            if c[8] == "1":           # alive
                a[e] += 1
    return {e: 100.0 * a[e] / n[e] for e in n}


def main():
    # cube[(E, depth)][engine] = survival%
    cube = defaultdict(dict)
    for eng in ENGINES:
        for d in DEPTHS:
            s = survival(eng, d)
            if s is None:
                print(f"% WARNING: missing {eng}_bench_{d}m.dat", file=sys.stderr)
                continue
            for e, pct in s.items():
                cube[(e, d)][eng] = pct

    # keep only informative cells (not all-0 and not all->=99.95); the
    # 99.95 guard keeps the near-threshold cells the paper discusses
    g4 = geant4_survival()
    rows = []
    for (e, d) in sorted(cube):
        vals = [cube[(e, d)].get(x, float("nan")) for x in ENGINES]
        g = g4[(round(e), d)]
        if all(v < 0.05 for v in vals + [g]) or all(v >= 99.95 for v in vals + [g]):
            continue
        rows.append((e, d, vals, g))

    # the numbers the caption states, from the rows themselves
    col = {x: ENGINES.index(x) for x in ENGINES}
    first = rows[0]                                   # 5 GeV at 10 m
    rest = rows[1:]
    dmp_rest = max(abs(v[col["MUSIC"]] - v[col["PROPOSAL"]]) for _, _, v, _ in rest)
    up_all = max(abs(v[col["UCMuon"]] - v[col["PROPOSAL"]]) for _, _, v, _ in rows)
    um_rest = max(abs(v[col["UCMuon"]] - v[col["MUSIC"]]) for _, _, v, _ in rest)
    v0, g0 = first[2], first[3]
    # how far the stochastic engines sit above Geant4 in the other cells
    above = [min(v[col[x]] for x in ("UCMuon", "MUSIC", "PROPOSAL")) - g for _, _, v, g in rest]
    print(r"\begin{table}[ht]")
    print(r"\centering")
    print(r"\caption{Depth-resolved survival fraction (\%) for \geantfour{} and "
          r"the four \ucmuon{} transport engines on identical monoenergetic vertical "
          r"beams (Standard Rock, $\rho=2.65$\,g\,cm$^{-3}$, $10^5$ muons per "
          r"cell). Only the transition cells ($0<\text{survival}<100$) are "
          r"shown; all omitted cells are below $0.05\%$ (below CSDA range) or "
          r"$\ge99.95\%$ (well above range) for every code. $\Delta_{\text{M--P}}$ is the "
          r"\music{}--\proposal{} difference, which stays below "
          f"${dmp_rest:.1f}$ percentage points in every cell but the first: "
          f"at {first[0]:g}\\,GeV and {first[1]}\\,m, near the end of the beam's "
          f"range, \\music{{}} passes ${v0[col['MUSIC']]:.1f}\\%$ against "
          f"${v0[col['PROPOSAL']]:.1f}\\%$ (\\proposal{{}}), ${v0[col['UCMuon']]:.1f}\\%$ "
          f"(UCMuon-MC) and ${g0:.1f}\\%$ (\\geantfour{{}}). The deterministic "
          r"Bethe--Bloch engine returns $0$ or $100\%$ (no range straggling); "
          f"UCMuon-MC agrees with \\proposal{{}} within ${up_all:.1f}$ percentage "
          f"points in every cell and with \\music{{}} within ${um_rest:.1f}$ in "
          r"every cell but the first. \geantfour{} (Standard Rock) passes fewer "
          f"muons than all three stochastic engines in every other cell, by "
          f"${min(above):.2f}$ to ${max(above):.2f}$ percentage points, which is where "
          r"the aggregate difference of Table~\ref{tab:survival}.}")
    print(r"\label{tab:survival_matrix}")
    print(r"\begin{tabular}{rrc" + "c" * len(ENGINES) + "c}")
    print(r"\toprule")
    print(r"$E$ [GeV] & $d$ [m] & \geantfour{} & " +
          " & ".join(LABELS[x] for x in ENGINES) +
          r" & $\Delta_{\text{M--P}}$ \\")
    print(r"\midrule")
    for e, d, vals, g in rows:
        dmp = vals[ENGINES.index("MUSIC")] - vals[ENGINES.index("PROPOSAL")]
        cells = " & ".join(f"{v:.2f}" for v in vals)
        print(f"{e:g} & {d} & {g:.2f} & {cells} & {dmp:+.2f} \\\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")


if __name__ == "__main__":
    main()
