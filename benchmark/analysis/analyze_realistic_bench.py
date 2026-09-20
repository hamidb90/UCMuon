#!/usr/bin/env python3
"""
4-engine aggregate analysis for the realistic-source benchmark.

Reads the per-engine, per-depth 18-column outputs produced by
run_realistic_bench.sh and produces:
  - benchmark_realistic_summary.csv  (per engine x depth aggregate stats)
  - fig_survival.png                 (survival fraction vs depth)
  - fig_exitKE.png                   (mean exit KE vs depth)
  - fig_lateral.png                  (RMS lateral displacement vs depth)
  - fig_exitKE_hist_<depth>m.png     (exit-KE distributions at each depth)

No Geant4/PHITS reference (the realistic source does not match the
monoenergetic grid those external runs used).

18-col layout (0-indexed):
  0 eventid 1 xs 2 ys 3 zs 4 E_srf 5 th_srf 6 ph_srf 7 charge
  8 alive 9 x_ug 10 y_ug 11 z_ug 12 E_ug 13 cx 14 cy 15 cz 16 th_ug 17 ph_ug
"""
import sys, glob, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = Path(__file__).resolve().parent
OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else \
         HERE.parent / "run_realistic_20260614"
FIGDIR = OUTDIR / "figures"
FIGDIR.mkdir(exist_ok=True)

# MUSIC first = reference column. BB_fortran included only if its dir exists.
ENGINES = ["MUSIC", "PROPOSAL", "BB", "UCMuon"]
if (Path(sys.argv[1]) if len(sys.argv) > 1 else
        HERE.parent / "run_realistic_20260614").joinpath("BB_fortran").is_dir():
    ENGINES.append("BB_fortran")
LABELS  = {"BB": "Bethe-Bloch (Py)", "BB_fortran": "Bethe-Bloch (Fortran)",
           "UCMuon": "UCMuon-MC", "MUSIC": "MUSIC", "PROPOSAL": "PROPOSAL"}
COLORS  = {"BB": "#2196F3", "BB_fortran": "#00BCD4", "MUSIC": "#FF5722",
           "UCMuon": "#4CAF50", "PROPOSAL": "#9C27B0"}
RHO = 2.65


def depths_for(engine):
    out = []
    for f in glob.glob(str(OUTDIR / engine / f"{engine}_bench_*m.dat")):
        tag = os.path.basename(f).split("_bench_")[1].replace("m.dat", "")
        try: out.append(int(tag))
        except ValueError: pass
    return sorted(out)


def load(engine, depth):
    f = OUTDIR / engine / f"{engine}_bench_{depth}m.dat"
    if not f.exists():
        return None
    rows = []
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 18:
                rows.append([float(x) for x in p[:18]])
    return np.array(rows) if rows else None


def stats(d):
    alive = d[:, 8].astype(int) == 1
    n     = len(d)
    ns    = int(alive.sum())
    E_srf = d[:, 4]
    E_ug  = d[:, 12]
    dx    = (d[:, 9]  - d[:, 1]) / 100.0      # cm -> m
    dy    = (d[:, 10] - d[:, 2]) / 100.0
    dr    = np.sqrt(dx**2 + dy**2)
    s     = dict(n=n, n_surv=ns, surv=100.0*ns/n if n else 0.0)
    if ns:
        s.update(
            E_srf_mean=float(E_srf[alive].mean()),
            E_ug_mean=float(E_ug[alive].mean()),
            E_ug_med=float(np.median(E_ug[alive])),
            dE_mean=float((E_srf[alive] - E_ug[alive]).mean()),
            dr_rms=float(np.sqrt((dr[alive]**2).mean())),
            dr_p95=float(np.percentile(dr[alive], 95)),
            E_ug_alive=E_ug[alive],
        )
    return s


# ── gather ──────────────────────────────────────────────────────────────────
all_depths = sorted(set().union(*[set(depths_for(e)) for e in ENGINES]))
if not all_depths:
    sys.exit(f"No engine outputs found under {OUTDIR}")

data = {e: {} for e in ENGINES}
for e in ENGINES:
    for dep in all_depths:
        d = load(e, dep)
        if d is not None:
            data[e][dep] = stats(d)

# ── CSV summary ───────────────────────────────────────────────────────────────
csv = OUTDIR / "benchmark_realistic_summary.csv"
with open(csv, "w") as fh:
    fh.write("engine,depth_m,mwe,n,n_surv,survival_pct,"
             "E_srf_mean_GeV,E_ug_mean_GeV,E_ug_median_GeV,dE_mean_GeV,"
             "lateral_rms_m,lateral_p95_m\n")
    for e in ENGINES:
        for dep in all_depths:
            s = data[e].get(dep)
            if not s:
                continue
            mwe = dep * RHO
            fh.write(f"{e},{dep},{mwe:.1f},{s['n']},{s['n_surv']},{s['surv']:.4f},"
                     f"{s.get('E_srf_mean',float('nan')):.3f},"
                     f"{s.get('E_ug_mean',float('nan')):.3f},"
                     f"{s.get('E_ug_med',float('nan')):.3f},"
                     f"{s.get('dE_mean',float('nan')):.3f},"
                     f"{s.get('dr_rms',float('nan')):.4f},"
                     f"{s.get('dr_p95',float('nan')):.4f}\n")
print(f"Wrote {csv}")

# ── console table: survival + exit KE, with delta vs MUSIC ────────────────────
ref = "MUSIC"
print("\nSurvival fraction (%)  [delta vs MUSIC in pp]")
hdr = "  depth/mwe   " + "".join(f"{LABELS[e]:>22}" for e in ENGINES)
print(hdr); print("  " + "-" * (len(hdr)-2))
for dep in all_depths:
    cells = []
    r = data[ref].get(dep, {}).get("surv")
    for e in ENGINES:
        v = data[e].get(dep, {}).get("surv")
        if v is None: cells.append(f"{'--':>22}"); continue
        d = "" if (e == ref or r is None) else f" ({v-r:+.2f})"
        cells.append(f"{v:>10.3f}{d:>12}")
    print(f"  {dep:>4}m/{dep*RHO:>5.0f}" + "".join(cells))

print("\nMean exit KE of survivors (GeV)")
for dep in all_depths:
    cells = [f"{data[e].get(dep,{}).get('E_ug_mean',float('nan')):>22.2f}" for e in ENGINES]
    print(f"  {dep:>4}m/{dep*RHO:>5.0f}" + "".join(cells))

# ── figures ───────────────────────────────────────────────────────────────────
def line_fig(key, ylabel, title, fname, logy=False):
    plt.figure(figsize=(7, 5))
    for e in ENGINES:
        xs = [d for d in all_depths if key in data[e].get(d, {})]
        ys = [data[e][d][key] for d in xs]
        if xs:
            plt.plot([x*RHO for x in xs], ys, "o-", color=COLORS[e],
                     label=LABELS[e], lw=2, ms=6)
    plt.xlabel("Depth (m.w.e.)"); plt.ylabel(ylabel); plt.title(title)
    if logy: plt.yscale("log")
    plt.grid(alpha=0.3); plt.legend()
    plt.tight_layout(); plt.savefig(FIGDIR / fname, dpi=150); plt.close()
    print(f"Wrote {FIGDIR/fname}")

line_fig("surv",    "Survival fraction (%)", "Muon survival vs depth (realistic source)", "fig_survival.png")
line_fig("E_ug_mean","Mean exit KE (GeV)",   "Mean exit KE of survivors vs depth",        "fig_exitKE.png")
line_fig("dr_rms",  "RMS lateral disp. (m)", "Lateral displacement vs depth",             "fig_lateral.png")

# exit-KE histograms per depth
for dep in all_depths:
    have = [e for e in ENGINES if "E_ug_alive" in data[e].get(dep, {})]
    if not have: continue
    plt.figure(figsize=(7, 5))
    for e in have:
        v = data[e][dep]["E_ug_alive"]
        plt.hist(v, bins=80, density=True, histtype="step", lw=2,
                 color=COLORS[e], label=f"{LABELS[e]} ({data[e][dep]['surv']:.1f}%)")
    plt.xlabel("Exit KE (GeV)"); plt.ylabel("Density")
    plt.title(f"Exit-KE distribution at {dep} m ({dep*RHO:.0f} m.w.e.)")
    plt.yscale("log"); plt.grid(alpha=0.3); plt.legend()
    plt.tight_layout(); plt.savefig(FIGDIR / f"fig_exitKE_hist_{dep}m.png", dpi=150); plt.close()
print(f"Wrote per-depth exit-KE histograms to {FIGDIR}/")
print("\nDone.")
