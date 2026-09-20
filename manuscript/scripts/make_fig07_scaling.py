#!/usr/bin/env python3
"""
make_fig07_scaling.py — Measure OpenMP strong scaling of the UCMuon
generator and MUSIC transport driver and generate fig07_scaling.pdf
(speedup + parallel efficiency panels).

Benchmark configurations (matches Table tab:scaling of the paper),
chosen to be compute-bound so the timing reflects the OpenMP region and
not the serialised ASCII output (every saved muon is written inside an
OMP critical section, so a run that saves every generated muon measures
file I/O, not parallel sampling):

  generator — detector-filtered configuration of hpc/input_params.dat
    (two 4 cm cylinder detectors at 35-90 m depth, disk source R = 500 m,
    CosmoALEPH 19-1500 GeV, cos^2-theta, theta_max = 85 deg);
    2000 *selected* muons, i.e. ~4e6 sampled trajectories;
  MUSIC — 1e6 muons, CosmoALEPH 10-2500 GeV, Standard Rock, d = 150 m.

Thread counts k = 1, 2, 4, 8.

Measured timings are written to manuscript/scripts/scaling_laptop.csv
(columns: program,threads,seconds).  Existing CSVs are reused unless
--rerun is given, so the figure can be re-styled without re-measuring.

HPC results: run run_scaling_hpc.sh on a cluster node to produce
manuscript/scripts/scaling_hpc.csv; if that file exists its curves are
added (dashed) to the plot.  The script also prints the LaTeX rows for
Table tab:scaling.

Usage:  python3 manuscript/scripts/make_fig07_scaling.py [--rerun]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ucmuon_run_helpers as rh

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF = SCRIPT_DIR.parent / "figs" / "fig07_scaling.pdf"
WORKDIR = SCRIPT_DIR / "_fig07_work"
LAPTOP_CSV = SCRIPT_DIR / "scaling_laptop.csv"
HPC_CSV = SCRIPT_DIR / "scaling_hpc.csv"

GEN_NSEL = 2_000           # selected muons (detector-filtered config)
MUSIC_N = 1_000_000        # muons in the MUSIC input file
DEPTH_M = 150.0
THREADS = [1, 2, 4, 8]

GEN_STDIN = None  # built in measure()


def measure() -> None:
    rh.check_binaries()
    rh.prepare_workdir(WORKDIR)
    rows = []
    gen_stdin = rh.gen_stdin(GEN_NSEL, e_min=19.0, e_max=1500.0,
                             detector=True, disk_radius_m=500.0,
                             outfile="muons_selected.dat")
    # MUSIC input file: 1e6 muons, 10-2500 GeV (generated once, not timed)
    if not (WORKDIR / "muons_surface.dat").exists():
        print(f"  [PREP] generating MUSIC input ({MUSIC_N:,} muons)")
        rh.run_binary(rh.GEN_EXE,
                      rh.gen_stdin(MUSIC_N, e_min=10.0, e_max=2500.0),
                      WORKDIR)
    # Warm-up run (touches cross-section tables, page cache)
    print("  [WARMUP]")
    rh.run_binary(rh.GEN_EXE, gen_stdin, WORKDIR, threads=2)
    rh.run_binary(rh.MUSIC_EXE, rh.music_stdin(depth_m=DEPTH_M), WORKDIR,
                  threads=2)
    for k in THREADS:
        t_gen = rh.run_binary(rh.GEN_EXE, gen_stdin, WORKDIR, threads=k)
        print(f"  [GEN]   k={k:2d}  {t_gen:7.2f} s")
        rows.append(("generator", k, t_gen))
        t_mus = rh.run_binary(rh.MUSIC_EXE, rh.music_stdin(depth_m=DEPTH_M),
                              WORKDIR, threads=k)
        print(f"  [MUSIC] k={k:2d}  {t_mus:7.2f} s")
        rows.append(("music", k, t_mus))
    with open(LAPTOP_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["program", "threads", "seconds"])
        w.writerows(rows)
    print(f"  [OK] wrote {LAPTOP_CSV}")


def load_csv(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    data: dict[str, list[tuple[int, float]]] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            data.setdefault(row["program"], []).append(
                (int(row["threads"]), float(row["seconds"])))
    out = {}
    for prog, vals in data.items():
        vals.sort()
        k = np.array([v[0] for v in vals])
        t = np.array([v[1] for v in vals])
        out[prog] = (k, t)
    return out


def latex_rows(timings, platform: str) -> list[str]:
    rows = []
    gen_k, gen_t = timings["generator"]
    mus_k, mus_t = timings["music"]
    for i, k in enumerate(gen_k):
        eff_g = gen_t[0] * gen_k[0] / (k * gen_t[i])
        eff_m = mus_t[0] * mus_k[0] / (k * mus_t[i])
        rows.append(f"  {k:2d} & {gen_t[i]:5.1f} & {eff_g:4.2f} & {platform}"
                    f" & {mus_t[i]:5.1f} & {eff_m:4.2f} & {platform} \\\\")
    return rows


def main():
    if "--rerun" in sys.argv or not LAPTOP_CSV.exists():
        measure()
    else:
        print(f"  [SKIP] reusing {LAPTOP_CSV.name} (--rerun to re-measure)")

    datasets = [("Laptop", load_csv(LAPTOP_CSV), "-")]
    if HPC_CSV.exists():
        datasets.append(("HPC", load_csv(HPC_CSV), "--"))
    else:
        print(f"  [INFO] {HPC_CSV.name} not found — laptop curves only "
              f"(run run_scaling_hpc.sh on the cluster)")

    style = {"generator": dict(color="#1565C0", marker="o"),
             "music": dict(color="#C62828", marker="s")}
    label = {"generator": "Generator", "music": "MUSIC transport"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    kmax = 1
    for plat, timings, ls in datasets:
        for prog, (k, t) in timings.items():
            s = t[0] * k[0] / t
            eff = s / k
            kmax = max(kmax, int(k.max()))
            ax1.plot(k, s, ls=ls, **style[prog],
                     label=f"{label[prog]} ({plat})")
            ax2.plot(k, eff, ls=ls, **style[prog],
                     label=f"{label[prog]} ({plat})")

    ax1.plot([1, kmax], [1, kmax], "k:", lw=1.0, label="Ideal")
    ax1.set_xlabel("Thread count $k$")
    ax1.set_ylabel("Speedup $S_k = T_1/T_k$")
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=2)
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8)

    ax2.axhline(0.85, ls="--", color="gray", lw=1.0)
    ax2.text(1.05, 0.855, "85% target", fontsize=8, color="gray",
             va="bottom")
    ax2.set_xlabel("Thread count $k$")
    ax2.set_ylabel(r"Parallel efficiency $\eta_k = S_k/k$")
    ax2.set_xscale("log", base=2)
    ax2.set_ylim(0, 1.1)
    ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF)
    print(f"[OK] wrote {OUT_PDF}")

    print("\nLaTeX rows for tab:scaling:")
    for plat, timings, _ in datasets:
        for row in latex_rows(timings, plat):
            print(row)


if __name__ == "__main__":
    main()
