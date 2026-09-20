#!/usr/bin/env python3
"""
ucmuon_run_helpers.py — drive the compiled UCMuon binaries (bin/) from the
manuscript figure scripts.  Used by make_fig05_energy_spectra.py and
make_fig07_scaling.py.

The binaries read their parameters from stdin (one value per line, same
order as the interactive prompts); see hpc/input_params.dat and
hpc/input_transport_music.dat for the documented sequences.
"""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BIN_DIR = REPO_ROOT / "bin"
GEN_EXE = BIN_DIR / "ucmuon_gen_omp"
MUSIC_EXE = BIN_DIR / "ucmuon_transport_music_omp"
MUSIC_TABLES = ["music-eloss-rock.dat", "music-cross-sections-rock.dat",
                "music-double-diff-rock.dat"]


def check_binaries() -> None:
    for exe in (GEN_EXE, MUSIC_EXE):
        if not exe.exists():
            raise SystemExit(
                f"[ERROR] {exe} not found — build the Fortran binaries first "
                f"(bash setup.sh or make).")


def prepare_workdir(workdir: Path) -> None:
    """Create workdir and link the MUSIC cross-section tables into it."""
    workdir.mkdir(parents=True, exist_ok=True)
    for name in MUSIC_TABLES:
        src, dst = BIN_DIR / name, workdir / name
        if src.exists() and not dst.exists():
            dst.symlink_to(src)


def gen_stdin(nmuons: int, e_min: float = 1.0, e_max: float = 2500.0,
              spectrum_mode: int = 1, angular_mode: int = 2,
              theta_max_deg: float = 85.0, disk_radius_m: float = 100.0,
              outfile: str = "muons_surface.dat",
              detector: bool = False) -> str:
    """stdin for ucmuon_gen_omp: disk source.

    detector=False — no filter, every generated muon is saved (save_all).
    detector=True  — the two-cylinder detector filter of
    hpc/input_params.dat is enabled and only detector-aimed muons are
    saved; `nmuons` is then the number of *selected* muons, with many
    sampled trajectories per saved muon (compute-bound configuration,
    used by the scaling benchmark)."""
    lines = [
        0,                # use_defaults
        e_min, e_max,
        spectrum_mode,
        1,                # source_mode  = disk
        1,                # source_plane = XY
        0.0, 0.0,         # disk centre U, V [m]
        disk_radius_m,
        0.0,              # W fixed [m]
        0.0, 0.0,         # tilt angle, tilt azimuth [deg]
        angular_mode,
        theta_max_deg,
        nmuons,
    ]
    if detector:
        lines += [
            1,                          # use_detector
            2,                          # ndet
            1, 3.0,                     # det 1: cylinder, 3 cm margin
            "990.0 0.0 -9000.0",        #   bottom cap centre [cm]
            "990.0 0.0 -3500.0",        #   top cap centre [cm]
            4.0,                        #   radius [cm]
            1, 3.0,                     # det 2: cylinder, 3 cm margin
            "0.0 -1650.0 -9000.0",
            "0.0 -1650.0 -3500.0",
            4.0,
            0,                          # save_all = 0 (selected only)
            0,                          # save_phits
            "muons_surface_unused.dat", # surface filename (not written)
            outfile,                    # selected muons filename
        ]
    else:
        lines += [
            0,                # use_detector
            1,                # save_all
            0,                # save_phits
            outfile,
        ]
    lines.append("")          # "Press Enter to start..."
    return "\n".join(str(v) for v in lines) + "\n"


def music_stdin(infile: str = "muons_surface.dat",
                outfile: str = "muons_underground.dat",
                depth_m: float = 100.0, rho: float = 2.65,
                rad_len: float = 26.48, init_tables: int = 1) -> str:
    """stdin for ucmuon_transport_music_omp (Standard Rock defaults)."""
    lines = [infile, outfile, rho, rad_len, depth_m,
             1,           # idim  — 3D lateral transport
             1,           # idim1 — other-process scattering
             -30,         # minv
             init_tables, # 1 = read cross-section tables from disk
             1,           # mat_type = rock
             1,           # transport ALL muons (asked for 14-col input)
             ""]          # "Press Enter to start..."
    return "\n".join(str(v) for v in lines) + "\n"


def run_binary(exe: Path, stdin: str, workdir: Path,
               threads: int | None = None) -> float:
    """Run a UCMuon binary; returns wall-clock seconds."""
    env = os.environ.copy()
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
    t0 = time.perf_counter()
    r = subprocess.run([str(exe)], input=stdin, capture_output=True,
                       text=True, cwd=workdir, env=env)
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        raise SystemExit(f"[ERROR] {exe.name} failed "
                         f"(exit {r.returncode}):\n{r.stdout[-2000:]}\n"
                         f"{r.stderr[-2000:]}")
    return dt
