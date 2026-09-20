#!/usr/bin/env python3
"""Generate the depth-resolved survival-fraction table for the UCMuon paper.

Reads the v2 benchmark outputs (six monoenergetic vertical beams at
5/10/20/50/100/300 GeV, 1e5 muons each, identical source population fed to
every engine) and emits a LaTeX table of survival fraction (%) per
(surface energy, depth) for the four UCMuon transport engines.

Authoritative data: benchmark/geant4_muon_rock_v5/<ENGINE>/<ENGINE>_bench_<DEPTH>m.dat
(the post-fix "v2" run of 2026-05-25 documented in BENCHMARK_FEEDBACK.md; the
flat ../output/*_bench_*.dat copies are an earlier v1 UCMuon run and must NOT
be used). The 10 m plane is excluded (5 GeV near-threshold boundary effect;
see Table~\ref{tab:survival}).

Column 5 (1-based) = surface energy [GeV]; column 9 = alive flag (0/1).
The .dat files are git-ignored (too large to track); regenerate with
compare_engines.py if absent.

Usage:  python make_survival_table.py  >  tab_survival_matrix.tex
"""
import sys
from pathlib import Path
from collections import defaultdict

BENCH    = Path(__file__).resolve().parent.parent.parent / "benchmark" / "geant4_muon_rock_v5"
ENGINES  = ["BB", "UCMuon", "MUSIC", "PROPOSAL"]      # column order
LABELS   = {"BB": "Bethe--Bloch", "UCMuon": "UCMuon-MC",
            "MUSIC": r"\music{}", "PROPOSAL": r"\proposal{}"}
DEPTHS   = [1, 25, 50, 100, 200]


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
    rows = []
    for (e, d) in sorted(cube):
        vals = [cube[(e, d)].get(x, float("nan")) for x in ENGINES]
        if all(v == 0 for v in vals) or all(v >= 99.95 for v in vals):
            continue
        rows.append((e, d, vals))

    print(r"\begin{table}[ht]")
    print(r"\centering")
    print(r"\caption{Depth-resolved survival fraction (\%) for the four "
          r"\ucmuon{} transport engines on identical monoenergetic vertical "
          r"beams (Standard Rock, $\rho=2.65$\,g\,cm$^{-3}$, $10^5$ muons per "
          r"cell). Only the transition cells ($0<\text{survival}<100$) are "
          r"shown; all omitted cells are $0\%$ (below CSDA range) or $\ge99.95\%$ "
          r"(well above range) for every engine. $\Delta_{\text{M--P}}$ is the "
          r"\music{}--\proposal{} difference, which stays below $0.7$ "
          r"percentage points everywhere. The deterministic Bethe--Bloch engine "
          r"returns $100\%$ in every transition cell (no range straggling); "
          r"UCMuon-MC agrees with \music{} and \proposal{} to within $1$ "
          r"percentage point in every cell.}")
    print(r"\label{tab:survival_matrix}")
    print(r"\begin{tabular}{rr" + "c" * len(ENGINES) + "c}")
    print(r"\toprule")
    print(r"$E$ [GeV] & $d$ [m] & " +
          " & ".join(LABELS[x] for x in ENGINES) +
          r" & $\Delta_{\text{M--P}}$ \\")
    print(r"\midrule")
    for e, d, vals in rows:
        dmp = vals[ENGINES.index("MUSIC")] - vals[ENGINES.index("PROPOSAL")]
        cells = " & ".join(f"{v:.2f}" for v in vals)
        print(f"{e:g} & {d} & {cells} & {dmp:+.2f} \\\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")


if __name__ == "__main__":
    main()
