#!/usr/bin/env python3
"""
Six-code benchmark comparison on the realistic generated source.

Combines:
  - UCMuon engines (18-col .dat): MUSIC, PROPOSAL, BB, UCMuon-MC, BB_fortran
        from benchmark/run_realistic_20260614/<ENGINE>/<ENGINE>_bench_<d>m.dat
  - Geant4 (per-crossing muons CSV): latest run_file_*_muons.csv
  - PHITS  (tally summary):          phits/phits_summary.csv

Geant4 is the external reference column (deltas computed vs it).
Produces benchmark_sixcode_summary.csv + survival/exit-KE figures.

PHITS is only included if its summary looks like the NEW run (total
N_transmitted at 1 m ~ 100k, not the old 600k monoenergetic run).
"""
import sys, glob, os, csv, collections, statistics
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO    = Path(__file__).resolve().parents[2]
OUTDIR  = REPO / "benchmark" / "run_realistic_20260614"
G4OUT   = REPO / "benchmark" / "geant4_muon_rock_v5" / "outputs"
PHSUM   = REPO / "benchmark" / "geant4_muon_rock_v5" / "phits" / "phits_summary.csv"
FIGDIR  = OUTDIR / "figures"; FIGDIR.mkdir(exist_ok=True)
N_SRC   = 100000
RHO     = 2.65
DEPTHS_M = [1, 10, 25, 50, 100, 200]
CM2M    = {100:1, 1000:10, 2500:25, 5000:50, 10000:100, 20000:200}

# display order: Geant4 reference first
ORDER  = ["Geant4", "PHITS", "MUSIC", "PROPOSAL", "BB", "BB_fortran", "UCMuon"]
LABELS = {"BB": "Bethe-Bloch (Py)", "BB_fortran": "Bethe-Bloch (Fortran)",
          "UCMuon": "UCMuon-MC", "MUSIC": "MUSIC", "PROPOSAL": "PROPOSAL",
          "Geant4": "Geant4", "PHITS": "PHITS"}
COLORS = {"BB": "#2196F3", "BB_fortran": "#00BCD4", "MUSIC": "#FF5722",
          "UCMuon": "#4CAF50", "PROPOSAL": "#9C27B0",
          "Geant4": "#111111", "PHITS": "#795548"}
ENGINE_DAT = ["MUSIC", "PROPOSAL", "BB", "UCMuon", "BB_fortran"]

data = {}   # data[code][depth_m] = dict(surv=%, ke=mean exit KE GeV, ke_list=optional)


def load_engine(engine, depth_m):
    f = OUTDIR / engine / f"{engine}_bench_{depth_m}m.dat"
    if not f.exists():
        return None
    al, eug = [], []
    with open(f) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 18:
                continue
            al.append(p[8] == "1" or p[8] == "1.0" or float(p[8]) == 1)
            eug.append(float(p[12]))
    al = np.array(al); eug = np.array(eug)
    ns = int(al.sum())
    return dict(surv=100.0*ns/len(al),
                ke=float(eug[al].mean()) if ns else float("nan"),
                ke_list=eug[al])


# ── engines ───────────────────────────────────────────────────────────────────
for e in ENGINE_DAT:
    if not (OUTDIR / e).is_dir():
        continue
    data[e] = {}
    for dm in DEPTHS_M:
        s = load_engine(e, dm)
        if s:
            data[e][dm] = s

# ── Geant4 (latest realistic muons CSV) ────────────────────────────────────────
g4files = sorted(glob.glob(str(G4OUT / "run_file_*_muons.csv")), key=os.path.getmtime)
if g4files:
    surv = collections.Counter(); ke = collections.defaultdict(list)
    with open(g4files[-1]) as fh:
        for r in csv.DictReader(fh):
            d = int(float(r["DepthCm"])); surv[d] += 1
            ke[d].append(float(r["ExitKEGeV"]))
    data["Geant4"] = {}
    for dcm, dm in CM2M.items():
        if surv[dcm]:
            data["Geant4"][dm] = dict(surv=100.0*surv[dcm]/N_SRC,
                                      ke=statistics.mean(ke[dcm]),
                                      ke_list=np.array(ke[dcm]))
    print(f"Geant4: {os.path.basename(g4files[-1])}")

# ── PHITS (only if summary is the new ~100k run) ───────────────────────────────
if PHSUM.exists():
    rows = list(csv.DictReader(open(PHSUM)))
    n1 = next((float(r["N_transmitted"]) for r in rows
               if int(float(r["DepthCm"])) == 100), 0)
    if n1 and n1 < 2 * N_SRC:          # ~100k, not the old 600k run
        data["PHITS"] = {}
        for r in rows:
            dm = CM2M.get(int(float(r["DepthCm"])))
            if dm:
                data["PHITS"][dm] = dict(surv=float(r["Transmission_%"]),
                                         ke=float(r["MeanExitKE_GeV"]),
                                         ke_list=None)
        print(f"PHITS: {PHSUM.name} (N@1m={n1:.0f})")
    else:
        print(f"PHITS: skipping stale summary (N@1m={n1:.0f}, looks like old run)")

codes = [c for c in ORDER if c in data]
print("Codes included:", ", ".join(codes))

# ── CSV ─────────────────────────────────────────────────────────────────────
csvf = OUTDIR / "benchmark_sixcode_summary.csv"
with open(csvf, "w") as fh:
    fh.write("code,depth_m,mwe,survival_pct,mean_exitKE_GeV\n")
    for c in codes:
        for dm in DEPTHS_M:
            s = data[c].get(dm)
            if s:
                fh.write(f"{c},{dm},{dm*RHO:.1f},{s['surv']:.4f},{s['ke']:.3f}\n")
print(f"Wrote {csvf}")

# ── console tables (delta vs Geant4) ──────────────────────────────────────────
ref = "Geant4" if "Geant4" in data else codes[0]
def table(metric, title, fmt):
    print(f"\n{title}  [delta vs {ref} in brackets]")
    head = "  d(m)/mwe " + "".join(f"{LABELS[c]:>22}" for c in codes)
    print(head); print("  " + "-"*(len(head)-2))
    for dm in DEPTHS_M:
        rv = data[ref].get(dm, {}).get(metric)
        cells = []
        for c in codes:
            v = data[c].get(dm, {}).get(metric)
            if v is None: cells.append(f"{'--':>22}"); continue
            d = "" if (c == ref or rv is None or metric=='ke') else f" ({v-rv:+.2f})"
            cells.append(f"{v:>{ '10' }{fmt}}{d:>12}")
        print(f"  {dm:>4}/{dm*RHO:>4.0f}" + "".join(cells))
table("surv", "Survival fraction (%)", ".3f")
table("ke",   "Mean exit KE of survivors (GeV)", ".2f")

# ── figures ───────────────────────────────────────────────────────────────────
def line_fig(metric, ylabel, title, fname):
    plt.figure(figsize=(8, 5.5))
    for c in codes:
        xs = [dm for dm in DEPTHS_M if metric in data[c].get(dm, {})]
        ys = [data[c][dm][metric] for dm in xs]
        ext = c in ("Geant4", "PHITS")
        plt.plot([x*RHO for x in xs], ys, "o-" if not ext else "s--",
                 color=COLORS[c], label=LABELS[c],
                 lw=3 if ext else 1.8, ms=8 if ext else 6,
                 zorder=5 if ext else 3)
    plt.xlabel("Depth (m.w.e.)"); plt.ylabel(ylabel); plt.title(title)
    plt.grid(alpha=0.3); plt.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(FIGDIR / fname, dpi=150); plt.close()
    print(f"Wrote {FIGDIR/fname}")

line_fig("surv", "Survival fraction (%)",
         "Muon survival vs depth — 6-code (realistic source)", "fig_sixcode_survival.png")
line_fig("ke", "Mean exit KE (GeV)",
         "Mean exit KE of survivors vs depth — 6-code", "fig_sixcode_exitKE.png")

# exit-KE histograms (codes with per-muon lists) at 100 m and 200 m
for dm in (100, 200):
    have = [c for c in codes if data[c].get(dm, {}).get("ke_list") is not None]
    if not have: continue
    plt.figure(figsize=(8, 5.5))
    for c in have:
        v = data[c][dm]["ke_list"]
        plt.hist(v, bins=80, density=True, histtype="step", lw=2.5 if c=="Geant4" else 1.8,
                 color=COLORS[c], label=f"{LABELS[c]} ({data[c][dm]['surv']:.1f}%)")
    plt.xlabel("Exit KE (GeV)"); plt.ylabel("Density"); plt.yscale("log")
    plt.title(f"Exit-KE distribution at {dm} m ({dm*RHO:.0f} m.w.e.)")
    plt.grid(alpha=0.3); plt.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(FIGDIR / f"fig_sixcode_exitKE_hist_{dm}m.png", dpi=150); plt.close()
    print(f"Wrote {FIGDIR}/fig_sixcode_exitKE_hist_{dm}m.png")
print("\nDone.")
