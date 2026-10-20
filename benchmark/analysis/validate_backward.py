#!/usr/bin/env python3
"""Validate the backward CSDA flux integrator (Engine 5, gui/ucmuon_backward_mc.py)
against the forward stochastic engines: UCMuon-MC, MUSIC and PROPOSAL.

The quantity is the one the terrain driver hands to the density inversion,
the transmission T(X, theta) = Phi(X, theta) / Phi(0, theta), with Phi the
directional flux above E_min at the detector, for Standard Rock, the Guan 2015
spectrum and E = 1-2500 GeV (the Vesuvius library's settings).

Forward estimate per (theta, X): surface kinetic energies are drawn from the
integrator's own surface spectrum _flux_surface(E, theta), so only the
transport differs, above E_cut = half the CSDA surface energy that reaches
E_min (below the ionisation-only minimum, so no survivor is cut). The muons
cross a slab of X g/cm^2 at normal incidence (a slab's transport does not
depend on direction; theta enters through the spectrum). Then

    T_fwd = [int_{E_cut} Phi dE / int_{E_min}^{E_max} Phi dE] x n(E_det in [E_min, E_max]) / N.

    python3 benchmark/analysis/validate_backward.py [WORKDIR] [--spectra 3]
        [--thetas 0,60,80,85] [--x 2e4,5e4,1e5,2e5,3e5]
        [--engines UCMuon-MC,MUSIC,PROPOSAL]

Writes WORKDIR/backward_validation.csv (one row per cell and engine) and
re-uses the transport outputs it finds there. The defaults are the validation
of 2026-10-19 (Guan). The calibration of mode 2 (2026-10-20) used all four
spectra, X = 1e4 ... 5e5 (8 values) and UCMuon-MC and MUSIC
(build_backward_calibration.py turns that run into the table); its test
between the nodes used --thetas 40,70,75,82 --x 3.5e4,7e4,1.5e5,3.5e5,4.5e5
and all three engines (column T_backward_calibrated).
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "gui"))
import ucmuon_backward_mc as bm  # noqa: E402

M_MU = 0.1056584                     # GeV
RHO = 2.65                           # g/cm^3, Standard Rock
SPEC = 3                             # Guan 2015 in the integrator's numbering (default)
E_MIN, E_MAX = 1.0, 2500.0           # GeV at the detector (Vesuvius library)
E_UP = 50_000.0                      # GeV, upper end of the surface sampling
N_MU = int(os.environ.get("NMU", 20_000))   # muons per cell and engine
THETAS = [0.0, 60.0, 80.0, 85.0]     # deg
X_GCM2 = [2.0e4, 5.0e4, 1.0e5, 2.0e5, 3.0e5]   # 200 to 3000 m w.e.
SEED = 20261019
BIN = REPO / "bin"
PROP_PY = Path.home() / "venvs" / "ucmuon" / "bin" / "python3"
PROP_TABLES = Path.home() / ".proposal_tables"


def backward_T(X: float, theta_deg: float, mode: int, spec: int = SPEC) -> float:
    th = np.radians(theta_deg)
    kw = dict(E_min_GeV=E_MIN, E_max_GeV=E_MAX, n_E=40, mode=mode)
    return (bm.directional_flux(X, th, spec, **kw) /
            bm.directional_flux(0.0, th, spec, **kw))


def surface_sample(theta_deg: float, X: float, rng, spec: int = SPEC) -> tuple[np.ndarray, float]:
    """N kinetic energies [GeV] above E_cut, and the flux fraction above E_cut."""
    th = np.radians(theta_deg)
    e_s0 = bm._backward_energy(E_MIN * 1000.0, X)[0] / 1000.0
    e_cut = max(E_MIN, 0.5 * e_s0)
    grid = np.logspace(np.log10(E_MIN), np.log10(E_UP), 4000)
    phi = np.array([float(bm._flux_surface(e, th, spec)) for e in grid])
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (phi[1:] + phi[:-1]) * np.diff(grid))])
    f_open = np.interp(E_MAX, grid, cum)                  # int_{E_MIN}^{E_MAX}
    c_cut = np.interp(e_cut, grid, cum)
    u = c_cut + rng.random(N_MU) * (cum[-1] - c_cut)
    return np.interp(u, cum, grid), (cum[-1] - c_cut) / f_open


def write_surface(path: Path, t_kin: np.ndarray, rng) -> None:
    e = t_kin + M_MU
    p = np.sqrt(e * e - M_MU * M_MU)
    q = np.where(rng.random(len(e)) < 0.55, 1, -1)       # charge only labels
    with open(path, "w") as fh:
        fh.write("# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV  theta_rad  "
                 "phi_rad  E_GeV  charge  hit_flag  det_mask\n")
        for i in range(len(e)):
            fh.write(f"{i + 1:10d} 0.0 0.0 0.0 {p[i]:.6f} 0.0 0.0 {-p[i]:.6f} 0.0 0.0 "
                     f"{e[i]:.6f} {q[i]:d} 1 0\n")


def run_engine(engine: str, surf: Path, out: Path, depth_m: float) -> None:
    if out.exists() and out.stat().st_size > 0:
        return
    log = open(str(out) + ".log", "w")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    if engine == "UCMuon-MC":
        stdin = "\n".join(map(str, [surf, out, depth_m, RHO, 10.015, 1, 1, 14, 0, 0.05, 1,
                                    11, 22, 136.4, 3.475e-6, 1, 2, 42, 1, 1])) + "\n"
        subprocess.run([sys.executable, str(REPO / "gui" / "ucmuon_stochastic_driver.py")],
                       input=stdin, text=True, stdout=log, stderr=log, check=True, env=env)
    elif engine == "MUSIC":
        stdin = "\n".join(map(str, [surf, out, RHO, 26.48, depth_m, 1, 1, -30, 1, 1, 1, ""])) + "\n"
        subprocess.run([str(BIN / "ucmuon_transport_music_omp")], cwd=BIN, input=stdin,
                       text=True, stdout=log, stderr=log, check=True, env=env)
    elif engine == "PROPOSAL":
        stdin = "\n".join(map(str, [surf, out, depth_m, 1, 1, 500, 0.001, 3, PROP_TABLES])) + "\n"
        subprocess.run([str(PROP_PY), str(REPO / "gui" / "proposal_driver.py")],
                       input=stdin, text=True, stdout=log, stderr=log, check=True, env=env)


def count(out: Path) -> int:
    a = np.loadtxt(out, comments="#", usecols=(8, 12), ndmin=2)
    t = a[:, 1] - M_MU
    return int(np.sum((a[:, 0] == 1) & (t >= E_MIN) & (t <= E_MAX)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("work", nargs="?", default=str(REPO / "benchmark" / "backward_validation_work"))
    ap.add_argument("--spectra", default=str(SPEC))
    ap.add_argument("--thetas", default=",".join(f"{x:g}" for x in THETAS))
    ap.add_argument("--x", default=",".join(f"{x:g}" for x in X_GCM2))
    ap.add_argument("--engines", default="UCMuon-MC,MUSIC,PROPOSAL")
    a = ap.parse_args()
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    spectra = [int(s) for s in a.spectra.split(",")]
    thetas = [float(s) for s in a.thetas.split(",")]
    xs = [float(s) for s in a.x.split(",")]
    engines = a.engines.split(",")
    cells = []
    for spec in spectra:
        rng = np.random.default_rng(SEED + spec)
        for th in thetas:
            for X in xs:
                tag = f"th{th:g}_X{X:g}" if spec == SPEC else f"s{spec}_th{th:g}_X{X:g}"
                surf = work / f"surface_{tag}.dat"
                t_kin, frac = surface_sample(th, X, rng, spec)
                if not surf.exists():
                    write_surface(surf, t_kin, rng)
                cells.append(dict(spec=spec, theta=th, X=X, tag=tag, surf=surf, frac=frac))
    jobs = [(c, e) for c in cells for e in engines]
    with ThreadPoolExecutor(max_workers=int(os.environ.get("WORKERS", "8"))) as pool:
        list(pool.map(lambda ce: run_engine(ce[1], ce[0]["surf"],
                                            work / f"{ce[1]}_{ce[0]['tag']}.dat",
                                            ce[0]["X"] / (RHO * 100.0)), jobs))
    rows = []
    for c in cells:
        tb1 = backward_T(c["X"], c["theta"], 1, c["spec"])
        tb0 = backward_T(c["X"], c["theta"], 0, c["spec"])
        tb2 = backward_T(c["X"], c["theta"], 2, c["spec"]) if bm.has_calibration(c["spec"]) else float("nan")
        for e in engines:
            n = count(work / f"{e}_{c['tag']}.dat")
            p = n / N_MU
            tf = c["frac"] * p
            sig = tf * np.sqrt((1 - p) / max(n, 1))
            rows.append(dict(spectrum_mode=c["spec"], theta_deg=c["theta"], X_gcm2=c["X"],
                             mwe=c["X"] / 100.0,
                             engine=e, T_forward=f"{tf:.6g}", T_forward_err=f"{sig:.3g}",
                             T_backward=f"{tb1:.6g}", T_backward_csda=f"{tb0:.6g}",
                             T_backward_calibrated=f"{tb2:.6g}",
                             ratio_backward_forward=f"{tb1 / tf:.4f}" if tf > 0 else "nan",
                             n_survivors=n, n_sampled=N_MU))
            print(f"spec {c['spec']}  theta {c['theta']:4.0f}  X {c['X']:8.0f}  {e:10s} T_fwd {tf:.4e} "
                  f"+- {sig / tf * 100 if tf else 0:4.1f}%   T_bwd {tb1:.4e} "
                  f"(ratio {tb1 / tf if tf else float('nan'):.3f})   CSDA-only {tb0:.4e}")
    with open(work / "backward_validation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"[ok] {work / 'backward_validation.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
