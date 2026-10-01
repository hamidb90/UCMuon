"""Shared plumbing for the joint-sampling and live-time tests.

Runs the Fortran generator (bin/ucmuon_gen_omp) through its stdin protocol,
including the detector filter that ucmugen/validation/ucmuref deliberately
leaves out, and ucmugen::Generator through tests/flux/ucmugen_sample.cc, and
returns events in one common layout: p, direction (3), position (3).
"""
from __future__ import annotations

import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))
FORTRAN = ROOT / "bin" / "ucmuon_gen_omp"
PARMA_DIR = ROOT / "data" / "EXPACS" / "parma"
M_MU = 0.10566


@dataclass
class Detector:
    shape: str                      # "cyl" or "box"
    geom: tuple                     # cyl: ax ay az bx by bz r ; box: x0 x1 y0 y1 z0 z1 (cm)
    margin: float = 0.0             # cm


@dataclass
class Config:
    spectrum: int = 4
    emin: float = 15.0
    emax: float = 1500.0
    theta_max_deg: float = 85.0
    angular_mode: int = 6
    source_mode: int = 1            # 1 disk, 2 rectangle, 3 hemisphere
    source_plane: int = 1
    radius_m: float = 30.0          # disk / hemisphere
    half_m: float = 30.0            # rectangle half-side (square)
    w_m: float = 0.0                # fixed coordinate / hemisphere centre z
    tilt_deg: float = 0.0
    tilt_az_deg: float = 0.0
    detector: Detector | None = None
    n: int = 200_000
    seed: int = 12345
    parma: dict = field(default_factory=lambda: {"lat": 50.7, "lon": 4.4, "alt_km": 0.0,
                                                 "year": 2020, "month": 1, "day": 1,
                                                 "charge_mode": 0, "sw": 0.0})

    # ------------------------------------------------------------------
    def stdin(self, out_all: str, out_sel: str) -> str:
        L = ["0", repr(self.emin), repr(self.emax), str(self.spectrum)]
        if self.spectrum == 3:
            P = self.parma
            L += [repr(P["lat"]), repr(P["lon"]), repr(P["alt_km"]), str(P["year"]),
                  str(P["month"]), str(P["day"]), str(P["charge_mode"]), str(PARMA_DIR),
                  repr(P["sw"])]
        L.append(str(self.source_mode))
        if self.source_mode in (1, 2):
            L.append(str(self.source_plane))
            if self.source_mode == 1:
                L += ["0.0", "0.0", repr(self.radius_m), repr(self.w_m)]
            else:
                L += [repr(-self.half_m), repr(self.half_m), repr(-self.half_m),
                      repr(self.half_m), repr(self.w_m)]
            L += [repr(self.tilt_deg), repr(self.tilt_az_deg)]
        else:
            L += [repr(self.radius_m), repr(self.w_m)]
        L.append(str(self.angular_mode))
        if self.angular_mode >= 2:
            L.append(repr(self.theta_max_deg))
        L.append(str(self.n))
        if self.detector is None:
            L += ["0", "1", "0", out_all]
        else:
            d = self.detector
            L += ["1", "1", "1" if d.shape == "cyl" else "2", repr(d.margin)]
            g = d.geom
            if d.shape == "cyl":
                L += [f"{g[0]} {g[1]} {g[2]}", f"{g[3]} {g[4]} {g[5]}", repr(g[6])]
            else:
                L += [f"{g[0]} {g[1]}", f"{g[2]} {g[3]}", f"{g[4]} {g[5]}"]
            L += ["0", "0", out_all, out_sel]
        L.append("")
        return "\n".join(L) + "\n"

    # ------------------------------------------------------------------
    def source(self):
        import live_time as LT
        return LT.Source(mode=self.source_mode, plane=self.source_plane,
                         radius_cm=self.radius_m * 100.0, half_lx_cm=self.half_m * 100.0,
                         half_ly_cm=self.half_m * 100.0, tilt_deg=self.tilt_deg,
                         tilt_az_deg=self.tilt_az_deg,
                         centre_z_cm=self.w_m * 100.0 if self.source_mode == 3 else 0.0)


@dataclass
class Run:
    p: np.ndarray
    d: np.ndarray          # (N, 3) world direction
    x: np.ndarray          # (N, 3) world position, cm
    tried: int | None = None
    rate: float | None = None
    rate_err: float | None = None
    stdout: str = ""


def run_fortran(cfg: Config, threads: int | None = None) -> Run:
    if not FORTRAN.exists():
        raise RuntimeError("build the generator first: make ucmuon_gen_omp")
    with tempfile.TemporaryDirectory() as tmp:
        out_all, out_sel = os.path.join(tmp, "all.dat"), os.path.join(tmp, "sel.dat")
        env = {**os.environ, "UCMUON_SEED": str(cfg.seed)}
        if threads:
            env["OMP_NUM_THREADS"] = str(threads)
        proc = subprocess.run([str(FORTRAN)], input=cfg.stdin(out_all, out_sel),
                              capture_output=True, text=True, env=env, cwd=tmp, timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout[-3000:] + proc.stderr[-2000:])
        path = out_sel if cfg.detector is not None else out_all
        raw = np.loadtxt(path, comments="#", ndmin=2)
    p = raw[:, 4]
    d = raw[:, 5:8] / p[:, None]
    x = raw[:, 1:4]
    m_tried = re.search(r"^\s*Tried\s*:\s*(\d+)", proc.stdout, re.M)
    m_rate = re.search(r"Surface rate R\s*:\s*([0-9.Ee+\-]+)", proc.stdout)
    return Run(p=p, d=d, x=x, tried=int(m_tried.group(1)) if m_tried else None,
               rate=float(m_rate.group(1)) if m_rate else None, stdout=proc.stdout)


_BUILT: dict[bool, Path] = {}


def _ucmugen_exe(parma: bool) -> Path:
    if parma not in _BUILT:
        exe = Path(tempfile.mkdtemp()) / ("ucmugen_sample_parma" if parma else "ucmugen_sample")
        cmd = ["c++", "-std=c++17", "-O2", f"-I{ROOT/'ucmugen/include'}",
               str(ROOT / "tests/flux/ucmugen_sample.cc"), "-o", str(exe)]
        if parma:
            cmd.insert(3, "-DUCMUGEN_TEST_PARMA")
        subprocess.run(cmd, check=True)
        _BUILT[parma] = exe
    return _BUILT[parma]


def run_ucmugen(cfg: Config, n: int | None = None, parma_site=None) -> Run:
    """ucmugen::Generator on the same surface and detector as the Fortran run.

    parma_site: (w, cutoff_GV, depth_gcm2, local_g) for spectrum 3.
    """
    n = cfg.n if n is None else n
    args = []
    if cfg.spectrum == 3:
        args += [repr(v) for v in parma_site]
    args += [str(cfg.spectrum), repr(cfg.emin), repr(cfg.emax), repr(cfg.theta_max_deg)]
    src = cfg.source()
    if cfg.source_mode == 3:
        args += ["hsphere", repr(src.radius_cm), repr(src.centre_z_cm)]
    else:
        nrm = src.flat_normal()
        # The Fortran places a flat source's centre at (0, 0, w) in the
        # canonical frame and permutes it with the plane.
        c = [0.0, 0.0, cfg.w_m * 100.0]
        if cfg.source_plane == 2:
            c = [c[0], c[2], c[1]]
        elif cfg.source_plane == 3:
            c = [c[2], c[0], c[1]]
        if cfg.source_mode == 1:
            args += ["disk", repr(src.radius_cm)]
        else:
            args += ["plane", repr(src.half_lx_cm), repr(src.half_ly_cm)]
        args += [repr(v) for v in c] + [repr(float(v)) for v in nrm]
    if cfg.detector is None:
        args.append("none")
    else:
        args += [cfg.detector.shape] + [repr(float(v)) for v in cfg.detector.geom] \
                + [repr(cfg.detector.margin)]
    args += [str(n), str(cfg.seed)]
    proc = subprocess.run([str(_ucmugen_exe(cfg.spectrum == 3))] + args,
                          capture_output=True, text=True, check=True, timeout=7200)
    lines = proc.stdout.splitlines()
    _, _, r, e = lines[0].split()
    raw = np.loadtxt(lines[1:], ndmin=2) if n > 0 else np.zeros((0, 7))
    return Run(p=raw[:, 0], d=raw[:, 1:4], x=raw[:, 4:7], rate=float(r), rate_err=float(e))


# ---------------------------------------------------------------------------
E_BINS = (15.0, 30.0, 100.0, 1500.0)
THETA_BINS = (0.0, 20.0, 40.0, 60.0, 85.0)


def binned(run: Run, e_bins=E_BINS, theta_bins=THETA_BINS) -> np.ndarray:
    E = np.sqrt(run.p ** 2 + M_MU ** 2)
    th = np.degrees(np.arccos(np.clip(-run.d[:, 2], -1.0, 1.0)))
    H, _, _ = np.histogram2d(E, th, bins=[e_bins, theta_bins])
    return H


def compare_bins(Ha: np.ndarray, Hb: np.ndarray):
    """Fraction ratio a/b per bin, its 1-sigma error, and the pull."""
    na, nb = Ha.sum(), Hb.sum()
    fa, fb = Ha / na, Hb / nb
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = fa / fb
        err = ratio * np.sqrt(1.0 / np.maximum(Ha, 1) + 1.0 / np.maximum(Hb, 1))
        pull = (fa - fb) / np.sqrt(fa * (1 - fa) / na + fb * (1 - fb) / nb)
    return ratio, err, np.nan_to_num(pull)


def run_fortran_mpi(cfg: Config, binary: Path, nranks: int = 3, threads: int = 2) -> Run:
    """The MPI+OpenMP generator (bin/ucmuon_gen): its stdin has no PHITS
    question and no "Press Enter" line, and each rank writes its own file."""
    import glob
    with tempfile.TemporaryDirectory() as tmp:
        out_all, out_sel = os.path.join(tmp, "all.dat"), os.path.join(tmp, "sel.dat")
        lines = cfg.stdin(out_all, out_sel).split("\n")
        # drop the PHITS answer (the "0" right after save_all) and the final Enter
        k = lines.index(out_all)
        del lines[k - 1]
        text = "\n".join(lines).rstrip("\n") + "\n"
        env = {**os.environ, "UCMUON_SEED": str(cfg.seed), "OMP_NUM_THREADS": str(threads)}
        proc = subprocess.run(["mpirun", "--oversubscribe", "-np", str(nranks), str(binary)],
                              input=text, capture_output=True, text=True, env=env, cwd=tmp,
                              timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout[-3000:] + proc.stderr[-2000:])
        stem = (out_sel if cfg.detector is not None else out_all)[:-4]
        files = sorted(glob.glob(stem + "_*.dat"))
        raw = np.concatenate([np.loadtxt(f, comments="#", ndmin=2) for f in files])
    p = raw[:, 4]
    m_tried = re.search(r"Total tried \(all ranks\):\s*(\d+)", proc.stdout)
    m_rate = re.search(r"Surface rate R\s*:\s*([0-9.Ee+\-]+)", proc.stdout)
    return Run(p=p, d=raw[:, 5:8] / p[:, None], x=raw[:, 1:4],
               tried=int(m_tried.group(1)) if m_tried else None,
               rate=float(m_rate.group(1)) if m_rate else None, stdout=proc.stdout)
