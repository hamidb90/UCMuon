"""Drive the UCMuon Fortran generator headlessly and parse its output.

The Fortran generator (``bin/ucmuon_gen_omp``) is interactive: it reads its
whole configuration from stdin, one value per line.  This module reproduces
that protocol so the generator can be run reproducibly from a test suite.

The protocol is defined by the ``read(*,*)`` sequence in
``src/generator/ucmuon_gen_omp.f90`` (lines 159-400).  It is mirrored, for the
GUI, by ``build_ucmuon_input`` in ``gui/ucmuon_gui.py``.  If any of those three
drift apart, :func:`GenConfig.stdin_text` is the one that breaks tests first,
which is deliberate.

Reproducibility relies on the ``UCMUON_SEED`` environment variable and on
``OMP_NUM_THREADS=1``.  Both are required: the generator's RNG is seeded per
thread, so a multi-threaded run is only reproducible for a fixed thread count,
and the output rows are written under an OpenMP CRITICAL whose ordering is not
deterministic.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BINARY = PROJECT_ROOT / "bin" / "ucmuon_gen_omp"
DEFAULT_PARMA_DIR = PROJECT_ROOT / "data" / "EXPACS" / "parma"

# Column layout of muons_surface.dat, from the write statement at
# src/generator/ucmuon_gen_omp.f90:944.
COLUMNS = (
    "event", "x_cm", "y_cm", "z_cm",
    "p_GeV", "px_GeV", "py_GeV", "pz_GeV",
    "theta_rad", "phi_rad", "E_GeV",
    "charge", "hit_flag", "det_mask",
)

SPECTRUM_NAMES = {
    1: "CosmoALEPH", 2: "PowerLaw", 3: "PARMA", 4: "Guan",
    5: "Frosin", 6: "GaisserBugaev", 7: "ReynaBugaev", 8: "Electron",
}
ANGULAR_NAMES = {
    1: "vertical", 2: "cos2", 3: "uniform-cone", 4: "guan-selfconsistent",
    5: "cos3",
}
SURFACE_NAMES = {1: "disk", 2: "rect", 3: "hemisphere"}


class GeneratorError(RuntimeError):
    """The Fortran generator exited non-zero or produced no usable output."""


@dataclass(frozen=True)
class GenConfig:
    """One generator configuration, in the generator's own units.

    Lengths are metres (the generator converts to cm internally and writes cm),
    momenta are GeV/c, angles are degrees.
    """

    # [1/7] energy range, GeV
    emin: float = 1.0
    emax: float = 1000.0

    # [2/7] spectrum model, 1-8 (see SPECTRUM_NAMES)
    spectrum: int = 4

    # [2b/7] PARMA site and date, used only when spectrum == 3
    parma_lat: float = 50.7
    parma_lon: float = 4.4
    parma_alt_km: float = 0.0
    parma_year: int = 2026
    parma_month: int = 1
    parma_day: int = 20
    parma_charge_mode: int = 0
    parma_sw: float = 0.0
    parma_dir: Path | None = None

    # [3/7] generation surface
    source_mode: int = 1        # 1 disk, 2 rectangle, 3 hemisphere
    source_plane: int = 1       # 1 XY, 2 XZ, 3 YZ (disk and rectangle only)
    disk_cx_m: float = 0.0
    disk_cy_m: float = 0.0
    disk_r_m: float = 2.0
    rect_u1_m: float = -2.0
    rect_u2_m: float = 2.0
    rect_v1_m: float = -2.0
    rect_v2_m: float = 2.0
    w_m: float = 0.0            # fixed offset along the plane normal
    tilt_deg: float = 0.0
    tilt_az_deg: float = 0.0
    hemi_radius_m: float = 2.0
    hemi_cz_m: float = 0.0

    # [4/7] angular distribution
    angular_mode: int = 2       # 1-5, see ANGULAR_NAMES
    theta_max_deg: float = 70.0

    # [5/7] statistics
    nmuons: int = 100_000

    # [6/7] detector filter is deliberately unsupported here: it changes which
    # muons reach the output file and would confound spectrum validation.

    # RNG
    seed: int = 20260729

    def __post_init__(self) -> None:
        if self.spectrum not in SPECTRUM_NAMES:
            raise ValueError(f"spectrum must be 1-8, got {self.spectrum}")
        if self.angular_mode not in ANGULAR_NAMES:
            raise ValueError(f"angular_mode must be 1-5, got {self.angular_mode}")
        if self.source_mode not in SURFACE_NAMES:
            raise ValueError(f"source_mode must be 1-3, got {self.source_mode}")
        if self.emax < self.emin:
            raise ValueError("emax must be >= emin")

    # ------------------------------------------------------------------
    @property
    def name(self) -> str:
        """Stable, filesystem-safe identifier for this configuration."""
        parts = [
            f"spec{self.spectrum}-{SPECTRUM_NAMES[self.spectrum]}",
            f"ang{self.angular_mode}-{ANGULAR_NAMES[self.angular_mode]}",
            f"surf{self.source_mode}-{SURFACE_NAMES[self.source_mode]}",
        ]
        if self.source_mode in (1, 2) and self.source_plane != 1:
            parts.append(f"plane{self.source_plane}")
        if self.tilt_deg:
            parts.append(f"tilt{self.tilt_deg:g}")
        return "_".join(parts)

    @property
    def digest(self) -> str:
        """Hash of the exact stdin, so a changed protocol invalidates caches."""
        return hashlib.sha256(self.stdin_text().encode()).hexdigest()[:16]

    def replace(self, **kw) -> "GenConfig":
        return replace(self, **kw)

    # ------------------------------------------------------------------
    def stdin_text(self, output_all: str = "OUTPUT") -> str:
        """Build the stdin stream, matching ucmuon_gen_omp.f90's read order."""
        L: list[str] = []

        L.append("0")                       # [0] use built-in defaults? no
        L.append(repr(float(self.emin)))    # [1/7] energy range
        L.append(repr(float(self.emax)))
        L.append(str(self.spectrum))        # [2/7] spectrum model

        if self.spectrum == 3:              # [2b/7] PARMA site and date
            parma_dir = self.parma_dir or DEFAULT_PARMA_DIR
            L += [
                repr(float(self.parma_lat)),
                repr(float(self.parma_lon)),
                repr(float(self.parma_alt_km)),
                str(self.parma_year),
                str(self.parma_month),
                str(self.parma_day),
                str(self.parma_charge_mode),
                str(parma_dir),
                repr(float(self.parma_sw)),
            ]

        L.append(str(self.source_mode))     # [3/7] generation surface
        if self.source_mode in (1, 2):
            L.append(str(self.source_plane))
            if self.source_mode == 1:
                L += [
                    repr(float(self.disk_cx_m)),
                    repr(float(self.disk_cy_m)),
                    repr(float(self.disk_r_m)),
                    repr(float(self.w_m)),
                ]
            else:
                L += [
                    repr(float(self.rect_u1_m)),
                    repr(float(self.rect_u2_m)),
                    repr(float(self.rect_v1_m)),
                    repr(float(self.rect_v2_m)),
                    repr(float(self.w_m)),
                ]
            L.append(repr(float(self.tilt_deg)))
            L.append(repr(float(self.tilt_az_deg)))
        else:
            L.append(repr(float(self.hemi_radius_m)))
            L.append(repr(float(self.hemi_cz_m)))

        L.append(str(self.angular_mode))    # [4/7] angular distribution
        if self.angular_mode in (2, 3, 4, 5):
            L.append(repr(float(self.theta_max_deg)))

        L.append(str(int(self.nmuons)))     # [5/7] statistics
        L.append("0")                       # [6/7] detector filter: off
        L.append("1")                       # [7/7] save surface file
        L.append("0")                       #       save PHITS source: no
        L.append(output_all)
        L.append("")                        # "Press Enter to start"
        return "\n".join(L) + "\n"


# ----------------------------------------------------------------------
def run(
    cfg: GenConfig,
    binary: Path | None = None,
    keep_dir: Path | None = None,
    timeout: float = 1800.0,
) -> np.ndarray:
    """Run the Fortran generator and return its events as a structured array.

    Raises :class:`GeneratorError` if the binary is missing, exits non-zero, or
    writes fewer rows than requested.
    """
    binary = Path(binary or DEFAULT_BINARY)
    if not binary.exists():
        raise GeneratorError(
            f"generator binary not found: {binary}\n"
            "Build it with `make ucmuon_gen_omp` from the project root."
        )

    workdir = Path(keep_dir) if keep_dir else Path(tempfile.mkdtemp(prefix="ucmuref-"))
    workdir.mkdir(parents=True, exist_ok=True)
    out_path = workdir / "muons_surface.dat"

    env = {
        **os.environ,
        "UCMUON_SEED": str(int(cfg.seed)),
        "OMP_NUM_THREADS": "1",   # required for reproducibility, see module docstring
    }

    try:
        proc = subprocess.run(
            [str(binary)],
            input=cfg.stdin_text(str(out_path)),
            capture_output=True,
            text=True,
            env=env,
            cwd=workdir,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise GeneratorError(
                f"{binary.name} exited {proc.returncode} for {cfg.name}\n"
                f"--- stdout tail ---\n{proc.stdout[-2000:]}\n"
                f"--- stderr tail ---\n{proc.stderr[-2000:]}"
            )
        if "UCMUON_SEED" not in proc.stdout:
            raise GeneratorError(
                "generator did not report a UCMUON_SEED-derived seed; the binary "
                "predates the seed fix. Rebuild with `make ucmuon_gen_omp`."
            )
        if not out_path.exists():
            raise GeneratorError(
                f"no output file for {cfg.name}\n--- stdout tail ---\n{proc.stdout[-2000:]}"
            )
        events = parse(out_path)
    finally:
        if keep_dir is None:
            shutil.rmtree(workdir, ignore_errors=True)

    if len(events) != cfg.nmuons:
        raise GeneratorError(
            f"{cfg.name}: expected {cfg.nmuons} events, got {len(events)}"
        )
    return events


def parse(path: Path) -> np.ndarray:
    """Parse muons_surface.dat into a structured array with COLUMNS fields."""
    raw = np.loadtxt(path, comments="#", ndmin=2)
    if raw.shape[1] != len(COLUMNS):
        raise GeneratorError(
            f"{path}: expected {len(COLUMNS)} columns, found {raw.shape[1]}. "
            "The generator's output format changed; update COLUMNS."
        )
    dtype = [(c, "i8" if c in ("event", "charge", "hit_flag", "det_mask") else "f8")
             for c in COLUMNS]
    out = np.empty(raw.shape[0], dtype=dtype)
    for j, c in enumerate(COLUMNS):
        out[c] = raw[:, j]
    return out


# ----------------------------------------------------------------------
def derived(events: np.ndarray) -> dict[str, np.ndarray]:
    """Quantities the comparison tests operate on.

    ``cos_theta`` is taken from the direction cosines rather than the stored
    ``theta_rad`` column, so that a bug in either one shows up as a
    disagreement instead of cancelling.
    """
    p = events["p_GeV"]
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_theta = -events["pz_GeV"] / p     # -pz because muons travel in -Z
    return {
        "p": p,
        "log10_p": np.log10(p),
        "cos_theta": cos_theta,
        "theta": events["theta_rad"],
        "phi": events["phi_rad"],
        "x": events["x_cm"],
        "y": events["y_cm"],
        "z": events["z_cm"],
        "E": events["E_GeV"],
        "charge": events["charge"].astype(float),
    }
