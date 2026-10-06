"""Fixtures for the GUI tests (Streamlit AppTest).

The GUI chdirs to the project root of the copy it runs from and writes
ucmuon_autosave.json there, so the tests never run it in the checkout: the
`repo_copy` fixture copies what the GUI needs to a temporary directory once
per session, and every AppTest runs from there.
"""
import importlib.util
import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "gui"))      # basic_mode etc. for the unit tests

# Not needed by the GUI (large, or source/docs only).
_SKIP = {".git", "output", "benchmark", "manuscript", "references", "external",
         "build", "src", "ucmugen", "hpc", "misc", "docs", "tests", "tools",
         "__pycache__"}

HAVE_RASTERIO = importlib.util.find_spec("rasterio") is not None
HAVE_GENERATOR = (ROOT / "bin" / "ucmuon_gen_omp").exists()


def _ignore(path, names):
    out = {n for n in names if n in _SKIP or n.startswith("ucmuon_autosave.json")}
    if Path(path) == ROOT:
        out |= {n for n in names if n.endswith(".dat") or n.endswith(".tif")}
    return out


@pytest.fixture(scope="session")
def repo_copy(tmp_path_factory):
    dst = tmp_path_factory.mktemp("ucmuon") / "repo"
    shutil.copytree(ROOT, dst, ignore=_ignore, symlinks=True)
    (dst / "output").mkdir(exist_ok=True)
    cwd = os.getcwd()
    yield dst
    os.chdir(cwd)                          # the GUI chdirs into the copy


@pytest.fixture
def new_app(repo_copy):
    """new_app(mode=None, state=None, keep_autosave=False) -> a run AppTest."""
    from streamlit.testing.v1 import AppTest

    def make(mode=None, state=None, keep_autosave=False):
        if not keep_autosave:
            for f in repo_copy.glob("ucmuon_autosave.json*"):
                f.unlink()
        at = AppTest.from_file(str(repo_copy / "gui" / "ucmuon_gui.py"), default_timeout=600)
        if mode:
            at.session_state["ui_mode"] = mode
        for k, v in (state or {}).items():
            at.session_state[k] = v
        at.run()
        return at

    return make


# ── small helpers shared by the test modules ─────────────────────────────────
def widget(at, kind, key):
    found = [x for x in getattr(at, kind) if x.key == key]
    return found[0] if found else None


_TAB_LABELS = {"Generator": "🌌  Generator", "Transport": "🪨  Transport",
               "Terrain": "🗺  Terrain", "Results": "📊  Results",
               "Density": "🔬  Density", "Config": "📋  Config"}


def top_tab(at, name):
    """One of the six main tabs, by exact label (sub-tabs, such as the
    Density tab's, can share icons and words with them)."""
    return [t for t in at.tabs if t.label == _TAB_LABELS[name]][0]


def exceptions(at):
    return [str(e.value)[:500] for e in at.exception]


@pytest.fixture(scope="session")
def surface_run(repo_copy):
    """A small surface muon file plus the run record the GUI keeps for it, so
    the Results/Transport tests need no Fortran generator. The muons are
    drawn roughly like cos²θ, 1-100 GeV; the rate numbers are those of a
    real 20 000-muon Guan run (R = 31 395 s⁻¹ on a 10 m disk)."""
    import time
    import numpy as np
    rng = np.random.default_rng(7)
    n = 5000
    th = np.arccos(rng.uniform(0.2, 1.0, n) ** (1.0 / 3.0))
    ph = rng.uniform(0.0, 2.0 * np.pi, n)
    E = 1.0 * (100.0 ** rng.uniform(0.0, 1.0, n))
    p = np.sqrt(E ** 2 - 0.10566 ** 2)
    r = 1000.0 * np.sqrt(rng.uniform(0.0, 1.0, n))
    a = rng.uniform(0.0, 2.0 * np.pi, n)
    rows = np.column_stack([np.arange(1, n + 1), r * np.cos(a), r * np.sin(a), np.zeros(n), p,
                            p * np.sin(th) * np.cos(ph), p * np.sin(th) * np.sin(ph),
                            -p * np.cos(th), th, ph, E, rng.choice([-1, 1], n),
                            np.ones(n), np.zeros(n)])
    t0 = time.time()
    f = repo_copy / "output" / "muons_surface.dat"
    np.savetxt(f, rows, fmt="%d %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.9f %.9f %.6f %d %d %d",
               header="EventID x_cm y_cm z_cm p_GeV px_GeV py_GeV pz_GeV theta_rad phi_rad "
                      "E_GeV charge hit_flag det_mask")
    state = dict(gen_run_files=[str(f.resolve())], gen_run_started=t0 - 5.0,
                 gen_ntry=20000, gen_surface_rate=31395.03, gen_spectrum_run=4,
                 gen_angular_mode=6, gen_emin=1.0, gen_emax=100.0, gen_theta_max=78.5,
                 gen_source_mode=1, gen_radius=10.0, surface_file=str(f), gen_use_detector=False)
    return f, state


@pytest.fixture(scope="session")
def underground_run(repo_copy, surface_run):
    """The surface file transported through 5 m of standard rock with the
    (pure Python) Bethe-Bloch driver, recorded as a transport of that run."""
    import subprocess
    f, state = surface_run
    ug = repo_copy / "output" / "muons_underground.dat"
    stdin = "\n".join([str(f), str(ug), "1", "14", "5", "1 2.65", "0"]) + "\n"
    subprocess.run([sys.executable, str(repo_copy / "gui" / "ucmuon_bb_driver.py")],
                   input=stdin, text=True, capture_output=True, check=True)
    os.utime(ug)
    return ug, dict(state, ug_file=str(ug), ug_depth_m=5.0, ug_rho=2.65,
                    ug_sources={str(ug.resolve()): str(f.resolve())})
