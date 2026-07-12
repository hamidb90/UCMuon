# ucmuon_gui.py  —  streamlit run ucmuon_gui.py
# UCMuon — UCLouvain Muography Group
# Author : Hamid Basiri <hamid.basiri@uclouvain.be>
# License: MIT
__version__ = "1.0.0"          # app version — keep in sync with CITATION.cff
import streamlit as st
import sys
from pathlib import Path as _PathSetup
_GUI_DIR_SETUP = str(_PathSetup(__file__).resolve().parent)
if _GUI_DIR_SETUP not in sys.path:   # sibling gui modules importable regardless
    sys.path.insert(0, _GUI_DIR_SETUP)  # of how Streamlit was launched

# CSDA engine removed — use MUSIC, Bethe-Bloch, PROPOSAL, or UCMuon engines for transport

# ── UCMuon-MC / Backward MC (Engines 1 & 5 in Tab 2) ─────────────────────────
try:
    from gui_stochastic_engine import (
        stochastic_available, render_stochastic_settings,
        render_backward_mc_tab, build_stochastic_stdin, _STOCHASTIC_MAT_ID,
    )
    _STOCHASTIC_GUI_OK = True
except ImportError:
    _STOCHASTIC_GUI_OK = False
    render_stochastic_settings = render_backward_mc_tab = build_stochastic_stdin = None
    _STOCHASTIC_MAT_ID = {}
    def stochastic_available(): return False, "gui_stochastic_engine.py not found in gui/"

# ── UCMuon Terrain (Tab 5) ─────────────────────────────────────────────────────
_TERRAIN_IMPORT_ERROR = None
try:
    from gui_terrain_engine import terrain_available, render_terrain_tab
    _TERRAIN_GUI_OK = True
except Exception as _terr_exc:          # catch ALL errors, not just ImportError
    _TERRAIN_GUI_OK = False
    _TERRAIN_IMPORT_ERROR = str(_terr_exc)
    render_terrain_tab = None
    def terrain_available(): return False, "gui_terrain_engine.py not found in gui/"

# ── Density Analysis (Tab 6) ──────────────────────────────────────────────────
try:
    from gui_density_analysis import render_density_analysis_tab
    _DENSITY_GUI_OK = True
except Exception as _dens_exc:
    _DENSITY_GUI_OK = False
    render_density_analysis_tab = None

# ── MCS Acceptance Estimator (Generator tab) ──────────────────────────────────────────
try:
    from gui_mcs_estimator import render_mcs_panel
    _MCS_GUI_OK = True
except ImportError:
    _MCS_GUI_OK = False
    def render_mcs_panel(_fn): pass

# ── Source Size + Optimisation panel (Generator tab) ──────────────────────────────────
try:
    from gui_source_optimizer import render_combined_source_panel
    _SRCOPT_GUI_OK = True
except ImportError:
    _SRCOPT_GUI_OK = False
    def render_combined_source_panel(*_a, **_kw): pass

import sys
import subprocess, threading, time, re, os, json, queue
import math
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pathlib import Path
from scipy.stats import norm as _spnorm, gaussian_kde as _spkde

# All binaries and data files live next to this .py file.
# Anchoring CWD here makes all relative paths (./binary, music-eloss.dat,
# ucmuon_autosave.json, etc.) work regardless of how Streamlit was launched.
_SCRIPT_DIR  = Path(__file__).resolve().parent   # gui/
_PROJECT_DIR = _SCRIPT_DIR.parent                  # project root
_BIN_DIR     = _PROJECT_DIR / "bin"                # compiled binaries + MUSIC symlinks

def _abspath(p):
    """Return absolute path string. File paths passed via stdin to Fortran
    binaries must be absolute because CWD is bin/, not project root."""
    if not p:
        return p
    from pathlib import Path as _P
    pp = _P(p)
    return str(pp) if pp.is_absolute() else str(_PROJECT_DIR / pp)
os.chdir(_PROJECT_DIR)  # working dir = project root so dat files resolve correctly
os.makedirs("output", exist_ok=True)  # simulation outputs go here, not in root
from fast_flux_estimator import (
    integrated_flux, flux_vs_depth, exposure_time,
    emin_from_opacity, RHO_STANDARD_ROCK,
    differential_flux, angular_profile, MODEL_LABELS as _FFE_MODEL_LABELS,
)

AUTOSAVE_FILE = "ucmuon_autosave.json"

def save_settings():
    """Collect all widget-keyed and run-result values and persist to disk.

    Uses an atomic write (write to .tmp then os.replace) to prevent a
    mid-crash from leaving a truncated / corrupt autosave file.
    """
    keys_to_save = [
        # ── Widget keys — must exactly match the key= parameter on each widget ─
        "gen_spectrum_mode", "source_mode",           # spectrum + source shape
        "parma_lat", "parma_lon", "parma_alt",        # PARMA location
        "parma_year", "parma_month", "parma_day",     # PARMA date
        "parma_charge", "parma_sw",  # PARMA charge / solar
        "transport_engine", "density_mode",            # transport / material
        "music_rho", "music_rad", "music_rho_sigma",  # material params
        "emin", "emax",                               # energy range
        "gen_mono", "gen_mono_e",                     # mono-energetic beam
        "radius", "sourcezm", "planelx", "planely",   # legacy (hemisphere)
        "source_plane",                               # plane orientation
        "disk_cx", "disk_cy", "disk_r",               # disk: center + radius
        "disk_tilt", "disk_tilt_az",                  # disk: tilt angles
        "src_u1", "src_u2", "src_v1", "src_v2",      # rectangle bounding box
        "src_w",                                      # fixed coordinate
        "nmuonsgen",                                  # N muons to generate
        "angularmode", "thetamax",                    # angular distribution
        "outputall", "saveall",                       # output file settings
        "savephits", "outputphits",                   # PHITS output
        "savegeant4",                                 # Geant4 output
        "usedetector", "outputsel", "ndet",           # detector filter
        "gen_workflow",                                # Standard vs DAS-REM mode
        "musicpreset", "music_depth_m", "nthreads",   # material, depth, threads
        "proposal_e_cut", "proposal_v_cut",           # PROPOSAL stochastic cuts
        "proposal_tables_dir", "proposal_med_choice", # PROPOSAL medium & tables
        "proposal_scatter_choice",                    # PROPOSAL scattering model
        "stochastic_vcut", "stochastic_nsteps", "stochastic_ms_enable",
        "stochastic_n_workers", "stochastic_delta_rays",  # UCMuon-MC v2 settings
        "stochastic_hard_spectrum", "stochastic_range_table",
        "bmc_depth", "bmc_rho", "bmc_mat_id", "bmc_spec",
        "bmc_theta_max", "bmc_mode", "bmc_nE", "bmc_nth",
        "bmc_vcut", "bmc_Emin", "bmc_Emax", "bmc_outfile",
        "pumas_mode", "pumas_energy_loss", "pumas_scattering",  # PUMAS engine
        "pumas_E_min", "pumas_E_max", "pumas_theta_max",
        "pumas_n_events", "pumas_spectrum_id", "pumas_seed", "pumas_outfile",
        "terrain_lat", "terrain_lon", "terrain_alt", "terrain_rho",
        "terrain_engine_choice", "terrain_naz", "terrain_nze",
        "terrain_zemax", "terrain_step", "terrain_dem_path",
        "terrain_dem_mode",       # DEM source radio (Upload vs Download)
        "terrain_outfile",        # terrain_infile intentionally omitted:
                                  # its selectbox key is NOT autosaved because a
                                  # stale path in session_state causes
                                  # StreamlitAPIException ("value not in options").
        # ── Run-result keys (set manually after a run) ────────────────────────
        "gen_radius", "gen_source_mode", "gen_source_z_m", "gen_plane_lx",
        "gen_plane_ly", "gen_nmuons_done", "gen_use_detector", "gen_emin",
        "gen_emax", "gen_angular_mode", "gen_thetamax",
        "surface_file", "selected_file", "ug_file",
    ]
    data = {k: st.session_state[k] for k in keys_to_save if k in st.session_state}
    # ── Dynamic detector widget keys (up to 10 detectors) ────────────────────
    _det_suffixes = ["sh", "mg", "ax", "ay", "az", "bx", "by", "bz", "rr",
                     "xn", "yn", "zn", "xx", "yx", "zx"]
    for _i in range(10):
        for _sfx in _det_suffixes:
            _dk = f"{_sfx}{_i}"
            if _dk in st.session_state:
                data[_dk] = st.session_state[_dk]
    # Stamp a schema version so load_settings can detect stale files
    data["_ucmuon_schema_version"] = 3

    # Atomic write: write to .tmp then rename — prevents corrupt file on crash
    _tmp = Path(AUTOSAVE_FILE + ".tmp")
    try:
        with open(_tmp, "w") as f:
            json.dump(data, f, indent=2, default=str)
        os.replace(str(_tmp), AUTOSAVE_FILE)   # atomic on POSIX and Win32
    except Exception as _e:
        # Non-fatal — next save will retry
        try:
            _tmp.unlink(missing_ok=True)
        except Exception:
            pass


def load_settings():
    """Load persisted settings into session_state BEFORE widgets render.

    Handles:
      - Missing file  → silent no-op
      - Empty file    → warns user, renames to .corrupt, starts fresh
      - Invalid JSON  → warns user, renames to .corrupt, starts fresh
      - Wrong schema  → silently ignores (incompatible old format)
    """
    p = Path(AUTOSAVE_FILE)
    if not p.exists():
        return

    # Empty file guard (left by a mid-crash truncation before our atomic fix)
    if p.stat().st_size == 0:
        _corrupt = Path(AUTOSAVE_FILE + ".corrupt")
        try:
            p.rename(_corrupt)
        except Exception:
            pass
        st.session_state.setdefault("_autosave_warn",
            f"⚠️ Autosave file was empty (corrupted by a previous crash). "
            f"Renamed to `{_corrupt.name}` — starting with default settings.")
        return

    try:
        with open(p) as f:
            raw = f.read()
        data = json.loads(raw)
    except (json.JSONDecodeError, ValueError) as _je:
        _corrupt = Path(AUTOSAVE_FILE + ".corrupt")
        try:
            p.rename(_corrupt)
        except Exception:
            pass
        st.session_state.setdefault("_autosave_warn",
            f"⚠️ Autosave JSON was corrupt ({_je}). "
            f"Renamed to `{_corrupt.name}` — starting with default settings. "
            f"You can restore it manually or delete `{_corrupt.name}` to clear the warning.")
        return
    except Exception:
        return  # other IO error — silent

    for k, v in data.items():
        if k not in st.session_state:   # only set if widget hasn't claimed it yet
            st.session_state[k] = v
load_settings()


def _apply_pending_restore():
    """Apply a Config-tab JSON restore stashed on the previous run.

    Widget-keyed session-state entries can only be written BEFORE their
    widget is instantiated in the current run (Streamlit raises
    StreamlitAPIException otherwise), so the Config tab stores the parsed
    JSON under _cfg_restore_pending and reruns; the values are applied
    here, ahead of every widget.
    """
    pending = st.session_state.pop("_cfg_restore_pending", None)
    if not pending:
        return
    n_applied = 0
    for k, v in pending.items():
        if k.startswith("_"):           # skip schema stamp / internal keys
            continue
        try:
            st.session_state[k] = v
            n_applied += 1
        except Exception:
            pass                        # unknown/incompatible key — skip
    # Keep the material-preset change detector in sync with the restored
    # preset, or it would fire and overwrite the restored ρ/X₀ with the
    # preset defaults on this very run.
    if "musicpreset" in pending:
        st.session_state["music_preset_prev"] = pending["musicpreset"]
    st.session_state["_cfg_restore_msg"] = n_applied
_apply_pending_restore()

# ── Show autosave corruption warning (set by load_settings on bad file) ──────
def _maybe_warn_autosave():
    msg = st.session_state.pop("_autosave_warn", None)
    if msg:
        st.warning(
            msg + "  Use **⚙️ Config → Reset autosave** to clear or restore from backup.",
            icon="⚠️"
        )
_maybe_warn_autosave()


# ══════════════════════════════════════════════════════════════════════════════
# MUSIC MATERIAL PRESETS
# ══════════════════════════════════════════════════════════════════════════════
# "rad" is the radiation length X₀ in g/cm² — MUSIC's native unit (the
# Kudryavtsev driver divides path lengths in g/cm² by it).  Engines that
# want X₀ in cm (UCMuon-MC) divide by the density.
MUSIC_MATERIALS = {
    "Standard Rock": {"rho": 2.65,  "rad": 26.48, "mat_id": 1, "mat_suffix": "rock",
                      "desc": "ρ=2.65 g/cm³, X₀=26.48 g/cm² — standard rock (Z=11, A=22)"},
    "Limestone":     {"rho": 2.71,  "rad": 26.10, "mat_id": 1, "mat_suffix": "rock",
                      "desc": "ρ=2.71 g/cm³, X₀=26.10 g/cm² — rock elemental XS (good approx. for CaCO₃)"},
    "Rock Salt":     {"rho": 2.17,  "rad": 29.00, "mat_id": 1, "mat_suffix": "rock",
                      "desc": "ρ=2.17 g/cm³, X₀=29.00 g/cm² — rock elemental XS (approx. for NaCl)"},
    "Water":         {"rho": 1.00,  "rad": 36.08, "mat_id": 2, "mat_suffix": "water",
                      "desc": "ρ=1.00 g/cm³, X₀=36.08 g/cm² — water/ice composition (H₂O)"},
    "Ice":           {"rho": 0.917, "rad": 36.08, "mat_id": 2, "mat_suffix": "water",
                      "desc": "ρ=0.917 g/cm³, X₀=36.08 g/cm² — water XS (same H₂O composition as Water)"},
    "Seawater":      {"rho": 1.025, "rad": 35.75, "mat_id": 3, "mat_suffix": "seawater",
                      "desc": "ρ=1.025 g/cm³, X₀=35.75 g/cm² — seawater (H, O, Na, Cl)"},
    "Iron":          {"rho": 7.874, "rad": 13.84, "mat_id": 1, "mat_suffix": "rock",
                      "desc": "ρ=7.874 g/cm³, X₀=13.84 g/cm² — rock elemental XS (APPROXIMATE for Fe)"},
    "Custom":        {"rho": None,  "rad": None,   "mat_id": 1, "mat_suffix": "rock",
                      "desc": "Uses rock elemental tables — accurate for rock-like media"},
}

# ---------------------------------------------------------------------------
# Helper: derive material-specific MUSIC table filenames from mat_suffix.
# rock     → legacy names (backward compatible with MUSIC distribution)
# water    → -water suffix
# seawater → -seawater suffix
# ---------------------------------------------------------------------------
# The MUSIC driver looks for material-specific filenames first.
# For rock group: music-eloss-rock.dat (driver primary), music-eloss.dat (fallback / as distributed by Kudryavtsev).
# For water:      music-eloss-water.dat
# For seawater:   music-eloss-seawater.dat
_ELOSS_FILES = {
    "rock":     "music-eloss-rock.dat",
    "water":    "music-eloss-water.dat",
    "seawater": "music-eloss-seawater.dat",
}
# Fallback filenames (as distributed by Kudryavtsev)
_ELOSS_FALLBACK = {
    "rock":     "music-eloss-rock.dat",
    "water":    "music-eloss-water.dat",
    "seawater": "music-eloss-seawater.dat",
}
_XSEC_FILES = {
    "rock":     "music-cross-sections-rock.dat",
    "water":    "music-cross-sections-water.dat",
    "seawater": "music-cross-sections-seawater.dat",
}
_XS_QUALITY = {
    "Standard Rock": "exact",
    "Limestone":     "good",
    "Rock Salt":     "good",
    "Water":         "exact",
    "Ice":           "exact",
    "Seawater":      "exact",
    "Iron":          "approx",
    "Custom":        "good",
}

def mat_files(mat_name):
    """Return (eloss_file, xsec_file, xs_quality) for a MUSIC_MATERIALS key."""
    mat    = MUSIC_MATERIALS[mat_name]
    suffix = mat.get("mat_suffix", "rock")
    eloss  = mat.get("eloss_file", _ELOSS_FILES.get(suffix, "music-eloss-rock.dat"))
    xsec   = mat.get("xsec_file",  _XSEC_FILES.get(suffix,  "music-cross-sections-rock.dat"))
    qual   = mat.get("xs_quality", _XS_QUALITY.get(mat_name, "good"))
    return eloss, xsec, qual


# ══════════════════════════════════════════════════════════════════════════════
# GROOM 2001 — MUON CSDA RANGE IN STANDARD ROCK
# Source: Groom, Mokhov & Striganov, ADNDT 78, 183–356 (2001), Table IV-6
# Standard Rock: <Z/A>=0.50000, ρ=2.65 g/cm³, I=136.4 eV
# T in MeV (kinetic energy),  CSDA range in g/cm²
# ══════════════════════════════════════════════════════════════════════════════
_GROOM_T_MEV = np.array([
    1.0e+01, 1.4e+01, 2.0e+01, 3.0e+01, 4.0e+01,
    8.0e+01, 1.0e+02, 1.4e+02, 2.0e+02, 3.0e+02,
    4.0e+02, 8.0e+02, 1.0e+03, 1.4e+03, 2.0e+03,
    3.0e+03, 4.0e+03, 8.0e+03, 1.0e+04, 1.4e+04,
    2.0e+04, 3.0e+04, 4.0e+04, 8.0e+04, 1.0e+05,
    1.4e+05, 2.0e+05, 3.0e+05, 4.0e+05, 8.0e+05,
    1.0e+06, 1.4e+06, 2.0e+06,
])
_GROOM_R_GCM2 = np.array([
    8.516e-01, 1.542e+00, 2.866e+00, 5.698e+00, 9.145e+00,
    2.676e+01, 3.696e+01, 5.879e+01, 9.332e+01, 1.524e+02,
    2.115e+02, 4.418e+02, 5.534e+02, 7.712e+02, 1.088e+03,
    1.599e+03, 2.095e+03, 3.998e+03, 4.920e+03, 6.724e+03,
    9.360e+03, 1.362e+04, 1.776e+04, 3.343e+04, 4.084e+04,
    5.495e+04, 7.459e+04, 1.040e+05, 1.302e+05, 2.129e+05,
    2.453e+05, 2.990e+05, 3.616e+05,
])


def _groom_threshold_energy(opacity_gcm2):
    """
    Log-log interpolation in the Groom (2001) Standard Rock CSDA table.
    Returns (E_GeV, E_MeV) — minimum muon kinetic energy to traverse opacity_gcm2.
    """
    if opacity_gcm2 <= _GROOM_R_GCM2[0]:
        return _GROOM_T_MEV[0] / 1000.0, _GROOM_T_MEV[0]
    if opacity_gcm2 >= _GROOM_R_GCM2[-1]:
        return _GROOM_T_MEV[-1] / 1000.0, _GROOM_T_MEV[-1]
    log_T = np.interp(
        np.log(opacity_gcm2),
        np.log(_GROOM_R_GCM2),
        np.log(_GROOM_T_MEV),
    )
    E_MeV = float(np.exp(log_T))
    return E_MeV / 1000.0, E_MeV


st.set_page_config(page_title="UCMuon", page_icon="🌌", layout="wide")



# ══════════════════════════════════════════════════════════════════════════════
# PERSISTENT THREAD-SAFE STATE
# ══════════════════════════════════════════════════════════════════════════════
if "_state" not in st.session_state:
    st.session_state["_state"] = {
        "gen_lines":       [],    "music_lines":       [],   
        "gen_running":     False, "music_running":     False,
        "gen_success":     None,  "music_success":     None, 
        "gen_stop_req":    False, "music_stop_req":    False,
        "gen_nmuons":      0,     "music_nmuons":      0,    
        "gen_proc":        None,  "music_proc":        None,  
        "gen_start_time":  None,  "music_start_time":  None,  
        "gen_end_time":    None,  "music_end_time":    None,  
    }
if "_lock" not in st.session_state:
    st.session_state["_lock"] = threading.Lock()


_STATE: dict           = st.session_state["_state"]
_LOCK:  threading.Lock = st.session_state["_lock"]


def _gs(k, v):
    with _LOCK:
        _STATE[k] = v


def _gg(k, default=None):
    with _LOCK:
        return _STATE.get(k, default)


def _ga(k, v):
    with _LOCK:
        _STATE[k].append(v)


def _gl(k):
    with _LOCK:
        return list(_STATE[k])


for _k, _v in [("surface_file",  "output/muons_surface.dat"),
               ("selected_file", "output/muons_selected.dat"),
               ("ug_file",       "output/muons_underground.dat")]:
    if _k not in st.session_state:
        st.session_state[_k] = _v


# ══════════════════════════════════════════════════════════════════════════════
# BACKGROUND WORKER
# ══════════════════════════════════════════════════════════════════════════════
def _worker(cmd, stdin_str, prefix, state, lock, env=None):

    def ws(k, v):
        with lock:
            state[k] = v

    def wg(k):
        with lock:
            return state[k]

    def wa(k, v):
        with lock:
            state[k].append(v)

    try:
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,   # merge stderr into stdout
            bufsize=0,                  # unbuffered — get output as it arrives
            cwd=str(_BIN_DIR),
            env=env,
        )
        ws(f"{prefix}_proc", proc)

        proc.stdin.write(stdin_str.encode())
        proc.stdin.flush()
        proc.stdin.close()

        # A dedicated reader thread puts lines into a queue so the main loop
        # can poll non-blocking and honour the stop button.  This replaces
        # pty + select, which are POSIX-only and unavailable on Windows.
        line_q = queue.Queue()

        def _reader(stream, q):
            try:
                for raw in iter(stream.readline, b""):
                    q.put(raw)
            finally:
                q.put(None)  # EOF sentinel

        threading.Thread(target=_reader, args=(proc.stdout, line_q),
                         daemon=True).start()

        while True:
            if wg(f"{prefix}_stop_req"):
                proc.kill()
                proc.wait()
                wa(f"{prefix}_lines", "⛔  Stopped by user.")
                ws(f"{prefix}_running", False)
                return

            try:
                raw = line_q.get(timeout=0.25)
            except queue.Empty:
                continue

            if raw is None:     # EOF — process finished
                break

            txt = raw.decode("utf-8", errors="replace").rstrip("\r\n")
            if txt.strip() and not _is_noise_line(txt):
                wa(f"{prefix}_lines", txt)

        # Drain any lines that arrived after the EOF sentinel
        while True:
            try:
                raw = line_q.get_nowait()
            except queue.Empty:
                break
            if raw is None:
                continue
            txt = raw.decode("utf-8", errors="replace").rstrip("\r\n")
            if txt.strip() and not _is_noise_line(txt):
                wa(f"{prefix}_lines", txt)

        proc.wait()
        ws(f"{prefix}_end_time", time.time())
        ws(f"{prefix}_running", False)
        ws(f"{prefix}_success", proc.returncode == 0)

    except Exception as exc:
        wa(f"{prefix}_lines", f"❌  {exc}")
        ws(f"{prefix}_running", False)
        ws(f"{prefix}_success", False)


def _dasrem_worker(cfg, state, lock):
    """Background thread for the DAS-REM generator (pure Python, no subprocess)."""
    def ws(k, v):
        with lock: state[k] = v
    def wa(k, v):
        with lock: state[k].append(v)
    def wg(k):
        with lock: return state.get(k)

    try:
        from ucmuon_dasrem_driver import generate_dasrem

        def _progress(n_done, N, line):
            if wg("gen_stop_req"):
                raise KeyboardInterrupt
            wa("gen_lines", line)

        n_written, elapsed = generate_dasrem(
            N             = cfg["nmuons"],
            spectrum_mode = cfg["spectrum_mode"],
            e_min         = cfg["e_min"],
            e_max         = cfg["e_max"],
            angular_mode  = cfg["angular_mode"],
            theta_max     = cfg["theta_max"],
            detectors     = cfg["detectors"],
            output_file   = cfg["output_file"],
            source_z_cm   = cfg.get("source_z_cm", 0.0),
            progress_fn   = _progress,
        )
        wa("gen_lines", f"  Saved {n_written} ... tried {n_written}  (100.0%)")
        wa("gen_lines", f"✅  DAS-REM complete: {n_written:,} muons → {cfg['output_file']}"
                        f"  ({elapsed:.1f} s)")
        ws("gen_end_time", time.time())
        ws("gen_running",  False)
        ws("gen_success",  True)

    except KeyboardInterrupt:
        wa("gen_lines", "⛔  Stopped by user.")
        ws("gen_end_time", time.time())
        ws("gen_running",  False)
        ws("gen_success",  False)
    except Exception as exc:
        wa("gen_lines", f"❌  DAS-REM error: {exc}")
        ws("gen_end_time", time.time())
        ws("gen_running",  False)
        ws("gen_success",  False)


import stat
for _exe in ["ucmuon_gen_omp", "ucmuon_transport_music_omp", "ucmuon_transport_bb_omp"]:
    _p = _BIN_DIR / _exe
    if _p.exists() and not os.access(_p, os.X_OK):
        _p.chmod(_p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)


def start_run(cmd, stdin_str, prefix, nmuons=0, env=None):
    _gs(f"{prefix}_lines",     [])
    _gs(f"{prefix}_running",   True)
    _gs(f"{prefix}_success",   None)
    _gs(f"{prefix}_stop_req",  False)
    _gs(f"{prefix}_nmuons",    nmuons)
    _gs(f"{prefix}_start_time", time.time())
    _gs(f"{prefix}_end_time",   None)
    threading.Thread(
        target=_worker,
        args=(cmd, stdin_str, prefix, _STATE, _LOCK, env),
        daemon=True
    ).start()


def stop_run(prefix):
    _gs(f"{prefix}_stop_req", True)
    proc = _gg(f"{prefix}_proc")
    if proc and proc.poll() is None:
        proc.kill()
    _gs(f"{prefix}_running", False)


# ══════════════════════════════════════════════════════════════════════════════
# PROGRESS PARSER
# ══════════════════════════════════════════════════════════════════════════════
def parse_progress(lines, nmuons):
    """
    Parse Fortran driver stdout for progress metrics.
    Returns (frac, saved, tried) where:
      tried = muons transported so far (progress denominator)
      saved = muons survived so far
      frac  = tried/nmuons (progress bar fraction)

    Priority order (highest first):
      1. OMP TRANSPORT COMPLETE final summary block — looks for
         'Muons transported: N' and 'Survived: M' as SEPARATE lines.
         These are the true totals; must be preferred over per-thread subtotals.
         The OMP driver prints per-thread progress lines like
         '  Transported: 9700  Survived: 884  Total: 10003' — these are
         PARTIAL counts from each thread, NOT the cumulative total.
         Without this priority, the last thread's subtotal (e.g. 871/10000)
         is mistakenly shown as the final result instead of 3572/10003.

      2. OMP live progress line (per-thread subtotal) — used while running
         '  Transported: N  Survived: M  Total: T'

      3. Surface generator format:
         'Saved N ... tried M'  /  'Accepted N ... tried M'
    """
    saved, tried = 0, 0

    # ── Priority 0: generator COMPLETE summary block ──────────────────────────
    # Generator prints (after COMPLETE banner):
    #   '  Saved:    100000'   ← colon distinguishes from live '  Saved N / tried M'
    #   '  Tried:    100000'
    _gen_saved, _gen_tried = 0, 0
    for line in lines:
        m = re.search(r'^\s*Saved\s*:\s*(\d+)', line, re.IGNORECASE)
        if m:
            _gen_saved = int(m.group(1))
        m = re.search(r'^\s*Tried\s*:\s*(\d+)', line, re.IGNORECASE)
        if m:
            _gen_tried = int(m.group(1))
    if _gen_saved > 0 and _gen_tried > 0:
        frac = min(_gen_saved / nmuons, 1.0) if nmuons > 0 else 1.0
        return frac, _gen_saved, _gen_tried

    # ── Priority 1: final summary block (separate lines) ─────────────────────
    # Scan all lines for the authoritative summary printed by TRANSPORT COMPLETE.
    # Pattern: '  Muons transported:     10003'  (with optional 'Muons' prefix)
    # and:     '  Survived:               3572'  (standalone Survived: line)
    # These appear AFTER 'TRANSPORT COMPLETE' or 'Writing output...' in the log.
    # We collect ALL matches and take the LARGEST values — the final summary
    # always has the total, whereas per-thread subtotals are smaller partials.
    all_transported = []
    all_survived    = []
    for line in lines:
        # '  Muons transported: 10003'  or  '  Transported: 10003'
        # but NOT '  Transported: 9700  Survived: 884  Total: 10003' (has Survived on same line)
        m_t = re.search(r'Muons\s+transported[:\s]+(\d+)', line, re.IGNORECASE)
        if m_t:
            all_transported.append(int(m_t.group(1)))
            continue
        # Standalone 'Survived: N' line (no 'Transported' on the same line)
        if re.search(r'Transported', line, re.IGNORECASE):
            continue          # skip combined progress lines here
        m_s = re.search(r'^\s*Survived[:\s]+(\d+)', line, re.IGNORECASE)
        if m_s:
            all_survived.append(int(m_s.group(1)))

    if all_transported:
        tried = max(all_transported)    # largest = true total
    if all_survived:
        saved = max(all_survived)       # largest = true total

    if tried > 0 and saved > 0:
        frac = min(tried / nmuons, 1.0) if nmuons > 0 else 0.0
        return frac, saved, tried

    # ── Priority 2: OMP live progress line (per-thread subtotal) ─────────────
    # Used while the run is in progress — shows partial thread counts.
    # Format: '  Transported: N  Survived: M  Total: T'
    # N and M are THIS THREAD's subtotals; T is the TRUE total across all threads.
    # Use the LATEST line's N for the progress fraction (most recent thread update)
    # and M only as a floor estimate for the survival count during the run.
    latest_T = 0    # true total from the Total field
    for line in reversed(lines):
        m = re.search(
            r'Transported[:\s]+(\d+)\s+Survived[:\s]+(\d+)\s+Total[:\s]+(\d+)',
            line, re.IGNORECASE)
        if m:
            # Use Total field as denominator so frac is never >1 during run
            _n, _s, _T = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if tried == 0:
                tried = _n
            if saved == 0:
                saved = _s
            latest_T = _T
            break

        # Fallback: no Total field (older format)
        m = re.search(
            r'Transported[:\s]+(\d+)\s+Survived[:\s]+(\d+)',
            line, re.IGNORECASE)
        if m:
            if tried == 0:
                tried = int(m.group(1))
            if saved == 0:
                saved = int(m.group(2))
            break

        # ── Surface generator: Saved N ... tried M ────────────────────────────
        m = re.search(r'Saved\s+(\d+).*tried\s*(\d+)', line, re.IGNORECASE)
        if m:
            saved, tried = int(m.group(1)), int(m.group(2))
            break

        m = re.search(r'Accepted\s*(\d+).*[Tt]ried\s*(\d+)', line, re.IGNORECASE)
        if m:
            saved, tried = int(m.group(1)), int(m.group(2))
            break

    # ── Priority 3: fallback separate-line search ─────────────────────────────
    if tried == 0:
        for line in reversed(lines):
            m = re.search(r'(?:Muons\s+)?[Tt]ransported[:\s]+(\d+)', line)
            if m:
                tried = int(m.group(1))
                break
    if saved == 0:
        for line in reversed(lines):
            m = re.search(r'^\s*Survived[:\s]+(\d+)', line, re.IGNORECASE)
            if m:
                saved = int(m.group(1))
                break

    # For OMP live lines: use the True Total (from the Total field) as denominator
    # so the fraction reflects global progress, not one thread's local count.
    # If latest_T > 0 we had an OMP-format live line; use it.
    denom = latest_T if latest_T > 0 else nmuons
    frac = min(tried / denom, 1.0) if denom > 0 and tried > 0 else 0.0
    return frac, saved, tried




# ══════════════════════════════════════════════════════════════════════════════
# LIVE PANEL
# ══════════════════════════════════════════════════════════════════════════════
def _fmt_time(seconds):
    if seconds is None:
        return "—"
    if seconds < 1e-3:
        return f"{seconds*1e3:.2f} ms"
    if seconds < 1.0:
        return f"{seconds*1e3:.0f} ms"
    s = int(seconds)
    if s < 60:
        return f"{s} s"
    elif s < 3600:
        return f"{s//60}m {s%60:02d}s"
    elif s < 86400:
        return f"{s//3600}h {(s%3600)//60:02d}m"
    else:
        return f"{s//86400}d {(s%86400)//3600:02d}h {(s%3600)//60:02d}m"


def _check_proposal():
    """
    Return (ok, version_or_message) for PROPOSAL availability.

    IMPORTANT: `import proposal` segfaults on Anaconda/miniforge due to a pybind11
    ABI mismatch. We must probe it in a subprocess — never in-process — to avoid
    crashing the Streamlit server.

    If PROPOSAL segfaults: install it in a system Python venv (not conda/miniforge):
        /usr/bin/python3 -m venv ~/venvs/ucmuon
        source ~/venvs/ucmuon/bin/activate
        pip install -r requirements.txt && pip install proposal
        streamlit run gui/ucmuon_gui.py    ← run from that venv
    """
    import subprocess as _sp
    _driver = _SCRIPT_DIR / "proposal_driver.py"
    if not _driver.exists():
        return False, "`proposal_driver.py` not found in gui/ folder"
    # Probe in a subprocess — crash-safe even under miniforge
    try:
        result = _sp.run(
            [sys.executable, "-c",
             "import proposal as pp; print(pp.__version__)"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0 and result.stdout.strip():
            return True, result.stdout.strip()
        if result.returncode == -11 or "segmentation" in (result.stderr or "").lower():
            return False, (
                "PROPOSAL segfaults under Anaconda/miniforge (pybind11 ABI mismatch). "
                "Use a system Python venv instead — see README.md §PROPOSAL."
            )
        return False, "`proposal` not installed — run: pip install proposal"
    except FileNotFoundError:
        return False, "`proposal` not installed — run: pip install proposal"
    except _sp.TimeoutExpired:
        return False, "PROPOSAL import timed out (possible segfault/hang)"
    except Exception as exc:
        return False, f"Check failed: {exc}"


_NOISE_PATTERNS = (
    "CrossSection.DNDX",          # PROPOSAL: negative dNdx at kinematic boundary (harmless)
    "Negative dNdx value",
    "Setting dNdx to zero",
    "No stochastic interaction possible",  # muon near stopping range — harmless
)

def _is_noise_line(txt: str) -> bool:
    """Return True for verbose-but-harmless PROPOSAL/spdlog output lines."""
    return any(p in txt for p in _NOISE_PATTERNS)


def live_panel(prefix):
    # Guard: keys are only created by start_run(); skip unconditionally-called
    # panels whose prefix (e.g. "bbms") has never been started.
    if f"{prefix}_lines" not in _STATE:
        return
    lines   = _gl(f"{prefix}_lines")
    running = _gg(f"{prefix}_running")
    nmuons  = _gg(f"{prefix}_nmuons")
    success = _gg(f"{prefix}_success")
    t_start = _gg(f"{prefix}_start_time")
    t_end   = _gg(f"{prefix}_end_time")

    if not lines and not running and success is None:
        return

    if t_start is not None:
        elapsed = (t_end if t_end else time.time()) - t_start
    else:
        elapsed = None

    frac, saved, tried = parse_progress(lines, nmuons)

    # For the surface generator, the user's target is N *saved* muons.
    # Using tried/nmuons overshoots 1.0 immediately whenever acceptance << 1
    # (detector filter ON: tried >> saved).  Use saved/nmuons for gen runs.
    # For MUSIC transport, every muon is transported exactly once, so
    # tried/nmuons (already computed in parse_progress) is correct.
    if prefix != "music" and nmuons > 0 and saved > 0:
        frac = min(saved / nmuons, 1.0)

    # Hard cap: never display a completely full bar while the process is alive.
    if running:
        frac = min(frac, 0.99)

    # Snap to 100 % once the process has truly finished with no parsed progress.
    if not running and success is True and frac == 0.0 and len(lines) > 3:
        frac  = 1.0
        tried = nmuons if nmuons > 0 else tried

    dn   = nmuons if nmuons > 0 else max(tried, saved, 1)
    icon = "⏳" if running else ("✅" if success else "❌")

    # Progress bar text: show transported (tried) as the live counter,
    # since that's what drives the fraction. Survived shown in metrics below.
    if prefix == "music":
        if running and tried > 0:
            prog_text = (f"{icon}  {tried:,} / {dn:,} transported  "
                         f"({100*frac:.1f}%)  —  {saved:,} survived so far")
        elif not running and success is True:
            prog_text = f"✅  {dn:,} / {dn:,} transported  (100.0%)"
        else:
            prog_text = f"{icon}  {tried:,} / {dn:,} muons  ({100*frac:.1f}%)"
    else:
        prog_text = f"{icon}  {saved:,} / {dn:,} muons  ({100*frac:.1f}%)"
    st.progress(frac, text=prog_text)

    # Metrics row — always show once we have any numbers
    _show_tried = tried if tried > 0 else dn if (not running and success) else 0
    _show_saved = saved if saved > 0 else (dn if (not running and success and prefix != "music") else 0)
    if _show_tried > 0 or _show_saved > 0:
        m1, m2, m3, m4 = st.columns(4)
        if prefix == "music":
            m1.metric("Survived",    f"{_show_saved:,}" if _show_saved > 0 else "—")
            m2.metric("Transported", f"{_show_tried:,}" if _show_tried > 0 else "—")
            # Survival rate: only show after run completes.
            # During an OMP run, per-thread "Survived" values are local checkpoint
            # tallies (e.g. Survived:923 at Transported:10000), NOT cumulative totals.
            # Showing 923/10000=9% during run then 3590/10003=35% at end is misleading.
            # The correct final value only appears in the TRANSPORT COMPLETE summary.
            if not running and _show_tried > 0 and _show_saved > 0:
                m3.metric("Survival rate",
                           f"{100.0*_show_saved/_show_tried:.4f} %")
            elif running:
                m3.metric("Survival rate", "⏳ pending…",
                          help="OMP per-thread partial counts are not meaningful. "
                               "True rate shown after transport completes.")
            else:
                m3.metric("Survival rate", "—")
        else:
            m1.metric("Selected",   f"{_show_saved:,}")
            m2.metric("Total tried",f"{_show_tried:,}" if _show_tried > 0 else "—")
            if _show_tried > 0 and _show_saved > 0:
                m3.metric("Net rate",
                           f"{100.0*_show_saved/_show_tried:.5f} %")
            else:
                m3.metric("Net rate", "—")
        time_label = "Elapsed" if running else "Total time"
        m4.metric(time_label, _fmt_time(elapsed))


    st.markdown("**Console output** *(last 60 lines)*")
    st.code("\n".join(lines[-60:]) if lines else "— waiting for output —",
            language="text")

    if not running:
        if success is True:
            rate_str = (f"  |  Net rate: {100.0*saved/abs(tried):.5f}%"
                        if tried != 0 else "")
            st.success(f"✅  Run completed in **{_fmt_time(elapsed)}**"
                       f"  |  {saved:,} muons saved{rate_str}")
        elif success is False:
            err = next((l for l in reversed(lines)
                        if any(k in l for k in ("❌", "ERROR", "Error"))), "")
            st.error(f"❌  Run failed after {_fmt_time(elapsed)}.  {err}")

    if running:
        time.sleep(0.25)
        st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
# STDIN BUILDERS
# ══════════════════════════════════════════════════════════════════════════════
def build_ucmuon_input(cfg):
    L = []
    L.append("0")                                        # use defaults = No
    L.append(str(cfg["emin"]))
    L.append(str(cfg["emax"]))

    # ── [2b/7] PARMA block — only when spectrum_mode == 3 ─────────────────
    L.append(str(cfg["spectrum_mode"]))

    if cfg.get("spectrum_mode") == 3:
        L.append(str(cfg.get("parma_lat",       50.7)))
        L.append(str(cfg.get("parma_lon",        4.4)))
        L.append(str(cfg.get("parma_alt",        0.0)))
        L.append(str(cfg.get("parma_year",      2026)))
        L.append(str(cfg.get("parma_month",        1)))
        L.append(str(cfg.get("parma_day",         20)))
        L.append(str(cfg.get("parma_charge",       0)))
        L.append(_abspath(str(_PROJECT_DIR / "data" / "EXPACS" / "parma")))
        L.append(str(cfg.get("parma_sw",         0.0)))   # Wolf/sunspot W index

    L.append(str(cfg.get("source_mode", 1)))          # source shape

    sm = cfg.get("source_mode", 1)
    if sm in (1, 2):
        L.append(str(cfg.get("source_plane", 1)))      # XY / XZ / YZ
        if sm == 1:
            # Disk: center + radius + tilt
            L.append(str(cfg.get("disk_cx",      0.0)))
            L.append(str(cfg.get("disk_cy",      0.0)))
            L.append(str(cfg.get("disk_r",     200.0)))
            L.append(str(cfg.get("src_w_m",     0.0)))
            L.append(str(cfg.get("disk_tilt",   0.0)))
            L.append(str(cfg.get("disk_tilt_az", 0.0)))
        else:
            # Rectangle: bounding box + tilt
            L.append(str(cfg.get("src_u1_m", -200.0)))
            L.append(str(cfg.get("src_u2_m",  200.0)))
            L.append(str(cfg.get("src_v1_m", -200.0)))
            L.append(str(cfg.get("src_v2_m",  200.0)))
            L.append(str(cfg.get("src_w_m",   0.0)))
            L.append(str(cfg.get("disk_tilt",    0.0)))
            L.append(str(cfg.get("disk_tilt_az", 0.0)))
    else:
        # Hemisphere: radius + centre z
        L.append(str(cfg.get("hemi_radius", 200.0)))
        L.append(str(cfg.get("hemi_cz_m",   0.0)))
    L.append(str(cfg["angular_mode"]))
    if cfg["angular_mode"] in [2, 3, 4, 5]:
        L.append(str(cfg["theta_max"]))
    L.append(str(cfg["nmuons"]))
    L.append("1" if cfg["use_detector"] else "0")
    if cfg["use_detector"]:
        L.append(str(len(cfg["detectors"])))
        for d in cfg["detectors"]:
            L.append(str(d["shape"]))
            L.append(str(d["margin"]))
            if d["shape"] == 1:
                L.append(f"{d['ax']} {d['ay']} {d['az']}")
                L.append(f"{d['bx']} {d['by']} {d['bz']}")
                L.append(str(d["r"]))
            else:
                L.append(f"{d['xmin']} {d['xmax']}")
                L.append(f"{d['ymin']} {d['ymax']}")
                L.append(f"{d['zmin']} {d['zmax']}")
    L.append("1" if cfg["save_all"]   else "0")
    L.append("1" if cfg["save_phits"] else "0")
    L.append(_abspath(cfg["output_all"]))
    if cfg["use_detector"]:
        L.append(_abspath(cfg["output_sel"]))
    if cfg["save_phits"]:
        L.append(_abspath(cfg["output_phits"]))
    L.append("")                                         # "Press Enter" line
    return "\n".join(L) + "\n"



def build_phitsxs_input(cfg):
    """Build stdin for ucmuon_bb_driver.py (Python CSDA) and Fortran ucmuon_transport_bb_omp.

    Parameter order (matching both Python driver and Fortran binary):
      1  infile
      2  outfile
      3  transport_all   0|1
      4  ncols_hint      auto-detected; sent for Fortran compat
      5  depth_m
      6  mat_type        1=StdRock 2=Ice 3=Water 4=Concrete 5=Custom
      [mat_type=5 only: Zeff  Aeff  rho_gcm3  I_eV]
      7  ms_enable       0|1
    """
    transport_all = 1 if cfg.get("transport_all", False) else 0
    ms_enable     = 1 if cfg.get("phitsxs_ms_enable", True) else 0
    mat_type      = cfg.get("phitsxs_mat_type", 1)
    lines = [
        _abspath(cfg.get("infile", "")),
        _abspath(cfg.get("outfile", "")),
        str(transport_all),
        str(int(cfg.get("ncols", 13))),
        str(cfg["depth_m"]),
        str(mat_type),
    ]
    if mat_type == 5:
        lines += [
            str(cfg.get("phitsxs_Zeff", 11.0)),
            str(cfg.get("phitsxs_Aeff", 22.0)),
            str(cfg.get("phitsxs_rho",  2.65)),
            str(cfg.get("phitsxs_I_eV", 136.4)),
        ]
    lines.append(str(ms_enable))
    return "\n".join(lines) + "\n"



def build_proposal_input(cfg):
    """
    Build stdin string for proposal_driver.py.

    stdin protocol (one value per line):
      1  infile
      2  outfile
      3  depth_m
      4  medium_type  (1=StandardRock, 2=Water, 3=Ice, 4=Seawater, 5=Custom)
      [if 5: Z_eff  A_eff  rho  I_eV]
      5  transport_all  (0/1)
      6  e_cut_MeV  (stochastic absolute energy cut)
      7  v_cut      (stochastic relative cut)
      8  scattering (0=none, 1=Highland, 2=HighlandIntegral, 3=Moliere)
      9  tables_dir (path for interpolation tables, empty = ~/.proposal/tables)
    """
    lines = [
        _abspath(cfg["infile"]),
        _abspath(cfg["outfile"]),
        str(cfg["depth_m"]),
        str(cfg.get("proposal_medium_type", 1)),
    ]
    if cfg.get("proposal_medium_type", 1) == 5:
        lines += [
            str(cfg.get("proposal_Z",   11.0)),
            str(cfg.get("proposal_A",   22.0)),
            str(cfg.get("proposal_rho",  2.65)),
            str(cfg.get("proposal_I_eV", 136.4)),
        ]
    lines += [
        "1" if cfg.get("transport_all", False) else "0",
        str(cfg.get("proposal_e_cut",    500.0)),
        str(cfg.get("proposal_v_cut",    0.001)),
        str(cfg.get("proposal_scatter",  2)),
        cfg.get("proposal_tables_dir",   ""),
        "",
    ]
    return "\n".join(lines) + "\n"


def build_pumas_input(cfg):
    """Build stdin for ucmuon_pumas_driver.py (mode=forward or backward)."""
    mode = cfg.get("pumas_mode", "backward")
    infile = _abspath(cfg.get("infile", "-")) if mode == "forward" else "-"
    return "\n".join([
        mode,
        infile,
        _abspath(cfg.get("outfile", "output/pumas_flux.dat")),
        str(cfg.get("depth_m", 100.0)),
        str(cfg.get("pumas_mat_id", 1)),
        str(cfg.get("rho", 2.65)),
        str(cfg.get("pumas_energy_loss", 0)),
        str(cfg.get("pumas_scattering", 0)),
        str(1 if cfg.get("transport_all", False) else 0),
        str(cfg.get("pumas_E_min",      1.0)),
        str(cfg.get("pumas_E_max",   1000.0)),
        str(cfg.get("pumas_theta_max",  85.0)),
        str(cfg.get("pumas_n_events", 50000)),
        str(cfg.get("pumas_spectrum_id",   0)),
        str(cfg.get("pumas_seed",          0)),
    ]) + "\n"


def build_music_input(cfg):
    mat_id        = cfg.get("mat_id", 1)
    transport_all = 1 if cfg.get("transport_all", False) else 0
    return (
        f"{_abspath(cfg['infile'])}\n"   # 1
        f"{_abspath(cfg['outfile'])}\n"  # 2
        f"{cfg['rho']}\n"               # 3
        f"{cfg['rad']}\n"          # 4
        f"{cfg['depth_m']}\n"      # 5
        f"{cfg['idim']}\n"         # 6
        f"{cfg['idim1']}\n"        # 7
        f"{cfg['minv']}\n"         # 8
        f"{cfg['init']}\n"         # 9
        f"{mat_id}\n"              # 10  ← mat_type (NEW)
        f"{transport_all}\n"       # 11  ← only consumed if 14-col file
        f"\n"                      # 12  ← Press Enter to start
    )



# ══════════════════════════════════════════════════════════════════════════════
# MUSIC FILE PROBE
# ══════════════════════════════════════════════════════════════════════════════
@st.cache_data(show_spinner=False)
def probe_music_file(filepath, transport_all, mtime=None):


    try:
        ncols = 0
        with open(filepath) as fh:
            for raw in fh:
                if raw.strip() == "" or raw.startswith("#"):
                    continue
                ncols = len(raw.split())
                break
        if ncols == 0:
            return 0, 0
        if ncols == 14 and not transport_all:
            df = pd.read_csv(filepath, sep=r"\s+", comment="#",
                             header=None, usecols=[12])
            n = int((df[12] == 1).sum())
        else:
            n = sum(1 for ln in open(filepath)
                    if ln.strip() and not ln.startswith("#"))
        return ncols, n
    except Exception:
        return 0, 0


# ══════════════════════════════════════════════════════════════════════════════
# FILE LOADER
# ══════════════════════════════════════════════════════════════════════════════
@st.cache_data
def load_file(path, mtime=None):
    c14 = "EventID,x,y,z,p,px,py,pz,theta,phi,E,charge,hit_flag,det_mask".split(",")
    c13 = "EventID,x,y,z,p,px,py,pz,theta,phi,E,charge,det_mask".split(",")
    cug = ("EventID,xs,ys,zs,Es,theta_s,phi_s,charge,alive,"
           "x,y,z,E,cx,cy,cz,theta,phi").split(",")
    c10 = "kf,x,y,z,cx,cy,cz,Ekin_MeV,wt,td".split(",")
    try:
        # ── Streaming load — avoids 3× memory peak of the old list+join+StringIO
        # approach for large files (10M+ rows).
        #
        # Fortran I10 right-justifies EventID with leading spaces ("         1").
        # pandas regex sep=r'\s+' treats a leading space as a delimiter, creating
        # a spurious NaN first column.  We fix this by dropping that column when
        # it is all-NaN, which is O(1) work after the read.
        df = pd.read_csv(path, sep=r'\s+', comment='#', header=None,
                         engine='python', skip_blank_lines=True)

        # Drop the spurious empty first column produced by leading-space EventIDs.
        if df.shape[1] > 1 and df.iloc[:, 0].isna().all():
            df = df.iloc[:, 1:].reset_index(drop=True)

        nc = df.shape[1]
        if   nc == 18: df.columns = cug
        elif nc == 14: df.columns = c14
        elif nc == 13: df.columns = c13
        elif nc == 10:
            df.columns = c10
            df["charge"] = np.where(df["kf"] == 13, -1, 1)
            df["E"]      = df["Ekin_MeV"] / 1000.0 + 0.105658
            df["theta"]  = np.arccos(np.clip(-df["cz"], -1.0, 1.0))
            df["phi"]    = np.arctan2(df["cy"], df["cx"])
            st.info("ℹ️  PHITS s-type=17 dump loaded (10 col). "
                    "Derived: charge, E [GeV], theta, phi.")
        elif nc == 5:
            df.columns = ["E_det_GeV", "flux", "flux_err", "E_surf_mean_GeV", "n_events_in_bin"]
        elif nc == 6:
            df.columns = ["ev", "E_det_GeV", "cos_theta", "charge", "E_surf_GeV", "flux_contribution"]
        else:
            df.columns = [f"col{i}" for i in range(nc)]
            st.warning(f"⚠️  Unrecognised file format: {nc} columns in `{path}`. "
                       f"Expected 5, 6, 10, 13, 14 or 18. Showing raw data.")
        return df
    except pd.errors.EmptyDataError:
        # File exists but has a header comment and no data rows — the usual cause
        # is a selection/cut that left zero muons, not a corrupt file.
        st.warning(f"⚠️  `{path}` contains no data rows (header only). "
                   f"The run or selection produced zero muons.")
        return pd.DataFrame()
    except Exception as ex:
        st.error(f"Could not load `{path}`: {ex}")
        return pd.DataFrame()



# ══════════════════════════════════════════════════════════════════════════════
# GEANT4 CONVERTERS
# ══════════════════════════════════════════════════════════════════════════════
MUON_MASS_MEV = 105.658
MUON_MASS_GEV = 0.105658

def _phits_d(val):
    """Fortran 1p1d24.15 — matches PHITS dump-a format (30(1p1d24.15))"""
    return f"{float(val):24.15E}".replace("E", "D")


def _g4_momentum_MeV(row):
    p   = row["p"] * 1000.0
    th  = row["theta"]
    ph  = row["phi"]
    return (
        p * np.sin(th) * np.cos(ph),
        p * np.sin(th) * np.sin(ph),
        -p * np.cos(th),
    )


def write_geant4_ascii(df, outpath):
    ecol = "Es" if "Es" in df.columns else "E"
    xcol = "xs" if "xs" in df.columns else "x"
    ycol = "ys" if "ys" in df.columns else "y"
    zcol = "zs" if "zs" in df.columns else "z"
    lines = [
        "# UCMuon Geant4 source file — UCLouvain Muography Group",
        "# Format: PDG  x[mm]  y[mm]  z[mm]  px[MeV/c]  py[MeV/c]  pz[MeV/c]  Ekin[MeV]",
        "# PDG: 13 = mu-   -13 = mu+",
    ]
    for _, row in df.iterrows():
        pdg     = 13 if row["charge"] < 0 else -13
        x_mm    = row[xcol] * 10.0
        y_mm    = row[ycol] * 10.0
        z_mm    = row[zcol] * 10.0
        px, py, pz = _g4_momentum_MeV(row)
        Ekin    = max(row[ecol] * 1000.0 - MUON_MASS_MEV, 0.0)
        lines.append(
            f"{pdg:4d}  {x_mm:12.4f}  {y_mm:12.4f}  {z_mm:12.4f}"
            f"  {px:14.6f}  {py:14.6f}  {pz:14.6f}  {Ekin:14.6f}"
        )
    with open(outpath, "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(df)


def write_geant4_hepevt(df, outpath):
    """
    G4HEPEvtInterface ASCII format — per event:
        NHEP
        ISTHEP IDHEP JDAHEP1 JDAHEP2 PX PY PZ MASS
    Geant4 reads exactly these 8 values per particle (energy is derived
    from p and mass, NOT read); ISTHEP must be 1 for a tracked particle.
    Momenta and mass in GeV.
    """
    lines = []
    for _, row in df.iterrows():
        pdg    = 13 if row["charge"] < 0 else -13
        p      = row["p"]
        th, ph = row["theta"], row["phi"]
        px =  p * np.sin(th) * np.cos(ph)
        py =  p * np.sin(th) * np.sin(ph)
        pz = -p * np.cos(th)
        lines.append("1")
        lines.append(
            f"1  {pdg:d}  0  0"
            f"  {px:.6e}  {py:.6e}  {pz:.6e}"
            f"  {MUON_MASS_GEV:.6f}"
        )
    with open(outpath, "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(df)

# ══════════════════════════════════════════════════════════════════════════════
# PHITS DUMP WRITERS
# ══════════════════════════════════════════════════════════════════════════════
def write_phits_surface(df, outpath):
    """
    PHITS s-type=17 dump — format: (30(1p1d24.15))
    No headers. 10 values per line: kf x y z cx cy cz Ekin[MeV] wt td
    """
    xcol = "xs" if "xs" in df.columns else "x"
    ycol = "ys" if "ys" in df.columns else "y"
    zcol = "zs" if "zs" in df.columns else "z"
    ecol = "E"  if "E"  in df.columns else "Es"
    with open(outpath, "w") as f:
        for _, r in df.iterrows():
            kf   = 13 if r["charge"] < 0 else -13
            th, ph = r["theta"], r["phi"]
            cx   =  np.sin(th) * np.cos(ph)
            cy   =  np.sin(th) * np.sin(ph)
            cz   = -np.cos(th)
            Ekin = max(r[ecol] * 1000.0 - MUON_MASS_MEV, 0.0)
            f.write(
                _phits_d(kf)      + _phits_d(r[xcol]) + _phits_d(r[ycol]) +
                _phits_d(r[zcol]) + _phits_d(cx)      + _phits_d(cy)      +
                _phits_d(cz)      + _phits_d(Ekin)    + _phits_d(1.0)     +
                _phits_d(0.0)     + "\n"
            )
    return len(df)


def write_phits_underground(src_path_or_df, outpath, survived_only=True,
                             progress_cb=None):
    """
    PHITS s-type=17 dump — format: (30(1p1d24.15))
    No headers. 10 values per line: kf x y z cx cy cz Ekin[MeV] wt td
    Accepts a DataFrame (small files) or a file path (streaming, no RAM limit).
    """
    # ── DataFrame mode ────────────────────────────────────────────────────
    if isinstance(src_path_or_df, pd.DataFrame):
        df_out = src_path_or_df[src_path_or_df["alive"] == 1].copy() \
                 if (survived_only and "alive" in src_path_or_df.columns) \
                 else src_path_or_df.copy()
        with open(outpath, "w") as f:
            for _, r in df_out.iterrows():
                kf   = 13 if r["charge"] < 0 else -13
                Ekin = max(r["E"] * 1000.0 - MUON_MASS_MEV, 0.0)
                f.write(
                    _phits_d(kf)    + _phits_d(r["x"])  + _phits_d(r["y"])  +
                    _phits_d(r["z"])+ _phits_d(r["cx"]) + _phits_d(r["cy"]) +
                    _phits_d(r["cz"])+ _phits_d(Ekin)   + _phits_d(1.0)     +
                    _phits_d(0.0)   + "\n"
                )
        return len(df_out)

    # ── Streaming file mode (18-col underground file) ─────────────────────
    # 0:EventID 1:xs 2:ys 3:zs 4:Es 5:theta_s 6:phi_s 7:charge
    # 8:alive   9:x 10:y 11:z 12:E 13:cx 14:cy 15:cz 16:theta 17:phi
    n_read = n_written = 0
    with open(src_path_or_df, "r") as fin, open(outpath, "w") as fout:
        for raw in fin:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()
            if len(parts) < 18:
                continue
            n_read += 1
            try:
                alive  = int(parts[8])
                E      = float(parts[12])
                charge = int(parts[7])
            except ValueError:
                continue
            if survived_only and alive != 1:
                continue
            if E <= MUON_MASS_GEV:
                continue
            x  = float(parts[9]);  y  = float(parts[10]); z  = float(parts[11])
            cx = float(parts[13]); cy = float(parts[14]); cz = float(parts[15])
            Ekin = max(E * 1000.0 - MUON_MASS_MEV, 0.0)
            kf   = 13 if charge < 0 else -13
            n_written += 1
            fout.write(
                _phits_d(kf)  + _phits_d(x)    + _phits_d(y)    +
                _phits_d(z)   + _phits_d(cx)    + _phits_d(cy)   +
                _phits_d(cz)  + _phits_d(Ekin)  + _phits_d(1.0)  +
                _phits_d(0.0) + "\n"
            )
            if progress_cb and n_read % 100_000 == 0:
                progress_cb(n_read, n_written)
    return n_written



# ══════════════════════════════════════════════════════════════════════════════
# CHAINED TRANSPORT CONVERTER  (18-col underground → 13-col MUSIC input)
# ══════════════════════════════════════════════════════════════════════════════
def convert_ug_to_music_input(src_path, outpath, survived_only=True,
                               progress_cb=None):
    """
    Streaming 18-col → 13-col converter.
    Reads line by line — works on files of any size without RAM issues.
    Columns in 18-col underground file:
      0:EventID  1:xs  2:ys  3:zs  4:Es  5:theta_s  6:phi_s  7:charge
      8:alive    9:x  10:y  11:z  12:E  13:cx  14:cy  15:cz  16:theta  17:phi
    """
    n_read = n_written = 0
    with open(src_path, "r") as fin, open(outpath, "w") as fout:
        fout.write(
            "# EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV"
            "  theta_rad  phi_rad  E_GeV  charge  det_mask\n"
        )
        for raw in fin:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()
            if len(parts) < 18:
                continue
            n_read += 1

            try:
                alive  = int(parts[8])
                E      = float(parts[12])
                charge = int(parts[7])
            except ValueError:
                continue

            # Skip stopped muons if requested
            if survived_only and alive != 1:
                continue
            if E <= MUON_MASS_GEV:
                continue

            x  = float(parts[9]);  y  = float(parts[10]); z  = float(parts[11])
            cx = float(parts[13]); cy = float(parts[14]); cz = float(parts[15])
            th = float(parts[16]); ph = float(parts[17])

            p  = np.sqrt(max(E**2 - MUON_MASS_GEV**2, 0.0))
            px = p * cx; py = p * cy; pz = p * cz

            n_written += 1
            fout.write(
                f"{n_written:10d}  {x:13.4f}  {y:13.4f}  {z:13.4f}"
                f"  {p:13.6f}  {px:13.6f}  {py:13.6f}  {pz:13.6f}"
                f"  {th:13.9f}  {ph:13.9f}  {E:13.6f}  {charge:4d}  0\n"
            )

            if progress_cb and n_read % 100_000 == 0:
                progress_cb(n_read, n_written)

    return n_read, n_written



def _auto_geant4_convert():
    if not st.session_state.get("save_geant4", False):
        return
    if st.session_state.get("gen_geant4_done", False):
        return
    if _gg("gen_running") or _gg("gen_success") is not True:
        return

    dst = st.session_state.get("g4_filename", "output/muons_geant4.txt")
    fmt = st.session_state.get("g4_fmt",      "ascii")

    # Pick source file: prefer all-muons when user requested it
    g4_use_all = st.session_state.get("gen_g4_use_all", True)
    src_all = st.session_state.get("gen_output_all", "")
    src_sel = (st.session_state.get("selected_file") or
               st.session_state.get("surface_file", ""))

    if g4_use_all and src_all and Path(_abspath(src_all)).exists():
        src = _abspath(src_all)
        src_label = "all generated muons"
    else:
        src = src_sel
        src_label = "detector hits"
        if g4_use_all and src_all:
            st.caption("ℹ️  All-muons file not found — Geant4 file uses detector hits. "
                       "Enable **Save ALL muons** before running to get all muons.")

    if not src or not Path(src).exists():
        st.warning(f"⚠️  Geant4 conversion skipped — source file not found.")
        return
    try:
        df_g4 = load_file(src)
        n = write_geant4_ascii(df_g4, dst) if fmt == "ascii" \
            else write_geant4_hepevt(df_g4, dst)
        st.success(f"✅  Geant4 file written: `{dst}`  ({n:,} muons — {src_label})")
        st.session_state["gen_geant4_done"] = True
    except Exception as ex:
        st.error(f"❌  Geant4 conversion failed: {ex}")

# ── PHITS surface export ──────────────────────────────────────

def _auto_phits_convert():
    """After a successful run, write or confirm the PHITS file.

    DAS-REM: converts the UCMuon output to PHITS format (Fortran never runs).
    Fortran: file was already written by the generator; just confirm it exists.
    """
    if not st.session_state.get("savephits", True):
        return
    if st.session_state.get("gen_phits_done", False):
        return
    if _gg("gen_running") or _gg("gen_success") is not True:
        return

    phits_dst = st.session_state.get("outputphits", "output/muons_for_phits.dat")
    is_dasrem = st.session_state.get("gen_use_dasrem", False)

    if is_dasrem:
        # DAS-REM does not run Fortran — convert UCMuon output to PHITS now
        src = (st.session_state.get("selected_file") or
               st.session_state.get("surface_file", ""))
        if not src or not Path(src).exists():
            st.warning("⚠️  PHITS export skipped — source file not found.")
            return
        try:
            df_p = load_file(src)
            n = write_phits_surface(df_p, phits_dst)
            st.success(f"✅  PHITS file written: `{phits_dst}`  ({n:,} muons)")
            st.session_state["gen_phits_done"] = True
        except Exception as ex:
            st.error(f"❌  PHITS conversion failed: {ex}")
    else:
        # Fortran generator writes PHITS directly — just confirm the file exists
        phits_abs = Path(_abspath(phits_dst))
        if phits_abs.exists():
            st.success(f"✅  PHITS file written: `{phits_dst}`")
            st.session_state["gen_phits_done"] = True
        else:
            st.warning(f"⚠️  PHITS file not found at `{phits_dst}`")


def _store_gen_params():
    """Parse integrated flux from generator console output and store key params."""
    if _gg("gen_running") or _gg("gen_success") is not True:
        return
    if st.session_state.get("gen_params_stored", False):
        return
    lines = _gl("gen_lines")
    phi = None
    for line in lines:
        m = re.search(
             r"integrated\s+flux\s*[:\-]?\s*([0-9]+\.?[0-9]*[Ee][+\-]?[0-9]+|[0-9]+\.?[0-9]+)",
             line, re.IGNORECASE
              )
        if m:
            try:
                phi = float(m.group(1))
            except ValueError:
                pass
            break
    # Parse total tried from final summary — generator prints "Tried: N"
    ntry_val = None
    for line in lines:
        m = re.search(r'Tried[:\s]+([0-9]+)', line, re.IGNORECASE)
        if m:
            try: ntry_val = int(m.group(1))
            except ValueError: pass
            break

    if phi is not None:
        st.session_state["gen_integrated_flux"] = phi
        st.session_state["gen_params_stored"]   = True
    if ntry_val is not None:
        st.session_state["gen_ntry"] = ntry_val



def _detector_traces(det, di):
    """High-contrast amber wireframe + filled volume for one detector."""
    traces  = []
    COLWIRE = "rgba(255, 215, 50, 1.0)"
    COLFILL = "rgb(255, 195, 30)"
    n       = 60

    if det["shape"] == 1:   # ── Cylinder ──────────────────────────────────
        rm  = det["r"]  / 100.0
        ax  = det["ax"] / 100.0;  ay = det["ay"] / 100.0;  az = det["az"] / 100.0
        bx  = det["bx"] / 100.0;  by = det["by"] / 100.0;  bz = det["bz"] / 100.0

        th     = np.linspace(0, 2 * np.pi, n, endpoint=False)
        xr_bot = ax + rm * np.cos(th);  yr_bot = ay + rm * np.sin(th)
        xr_top = bx + rm * np.cos(th);  yr_top = by + rm * np.sin(th)

        # Side surface — two non-degenerate triangles per quad
        xs_s = np.concatenate([xr_bot, xr_top])
        ys_s = np.concatenate([yr_bot, yr_top])
        zs_s = np.concatenate([np.full(n, az), np.full(n, bz)])
        tis, tjs, tks = [], [], []
        for i in range(n):
            j = (i + 1) % n
            tis.append(i);  tjs.append(j);     tks.append(n + i)
            tis.append(j);  tjs.append(n + j); tks.append(n + i)
        traces.append(go.Mesh3d(
            x=xs_s, y=ys_s, z=zs_s, i=tis, j=tjs, k=tks,
            color=COLFILL, opacity=0.40,
            showlegend=False, hoverinfo="skip",
        ))

        # Caps — each at its own axis-end centre
        for idxz, (cx, cy, zz, xr, yr) in enumerate([
                (ax, ay, az, xr_bot, yr_bot),
                (bx, by, bz, xr_top, yr_top),
        ]):
            xcap = np.concatenate([[cx], xr])
            ycap = np.concatenate([[cy], yr])
            ic   = np.zeros(n, dtype=int)
            jc   = np.arange(1, n + 1)
            kc   = np.array([*range(2, n + 1), 1])
            traces.append(go.Mesh3d(
                x=xcap, y=ycap, z=np.full(n + 1, zz),
                i=ic, j=jc, k=kc,
                color=COLFILL, opacity=0.55,
                name=f"Detector {di + 1} (cylinder)" if idxz == 0 else None,
                showlegend=(idxz == 0), hoverinfo="skip",
            ))

        # Wireframe rim circles
        ang = np.linspace(0, 2 * np.pi, 121)
        for cx, cy, zz in [(ax, ay, az), (bx, by, bz)]:
            traces.append(go.Scatter3d(
                x=cx + rm * np.cos(ang), y=cy + rm * np.sin(ang),
                z=np.full(121, zz),
                mode="lines", line=dict(color=COLWIRE, width=3),
                showlegend=False, hoverinfo="skip",
            ))

        # 12 vertical staves
        for angle in np.linspace(0, 2 * np.pi, 13)[:-1]:
            traces.append(go.Scatter3d(
                x=[ax + rm * np.cos(angle), bx + rm * np.cos(angle)],
                y=[ay + rm * np.sin(angle), by + rm * np.sin(angle)],
                z=[az, bz],
                mode="lines", line=dict(color=COLWIRE, width=2),
                showlegend=False, hoverinfo="skip",
            ))

        # Centre spine A→B with endpoint markers
        traces.append(go.Scatter3d(
            x=[ax, bx], y=[ay, by], z=[az, bz],
            mode="lines+markers",
            line=dict(color=COLWIRE, width=4, dash="dot"),
            marker=dict(size=6, color=COLWIRE, symbol="circle"),
            showlegend=False, hoverinfo="skip",
        ))

    else:   # ── Box (AABB) ─────────────────────────────────────────────────
        xn, xx_ = det["xmin"] / 100, det["xmax"] / 100
        yn, yx_ = det["ymin"] / 100, det["ymax"] / 100
        zn, zx_ = det["zmin"] / 100, det["zmax"] / 100

        xv = [xn, xx_, xx_, xn,  xn, xx_, xx_, xn]
        yv = [yn, yn,  yx_, yx_, yn, yn,  yx_, yx_]
        zv = [zn, zn,  zn,  zn,  zx_, zx_, zx_, zx_]
        ti = [0, 0,  4, 4,  0, 0,  2, 2,  0, 0,  1, 1]
        tj = [1, 2,  5, 6,  1, 5,  3, 7,  3, 7,  2, 6]
        tk = [2, 3,  6, 7,  5, 4,  7, 6,  7, 4,  6, 5]
        traces.append(go.Mesh3d(
            x=xv, y=yv, z=zv, i=ti, j=tj, k=tk,
            color=COLFILL, opacity=0.40,
            name=f"Detector {di + 1} (box)", showlegend=True, hoverinfo="skip",
        ))

        # 12 wireframe edges
        for zz in [zn, zx_]:
            traces.append(go.Scatter3d(
                x=[xn, xx_, xx_, xn, xn], y=[yn, yn, yx_, yx_, yn], z=[zz] * 5,
                mode="lines", line=dict(color=COLWIRE, width=3),
                showlegend=False, hoverinfo="skip",
            ))
        for vx_, vy_ in [(xn, yn), (xx_, yn), (xx_, yx_), (xn, yx_)]:
            traces.append(go.Scatter3d(
                x=[vx_, vx_], y=[vy_, vy_], z=[zn, zx_],
                mode="lines", line=dict(color=COLWIRE, width=3),
                showlegend=False, hoverinfo="skip",
            ))

        # Centre diamond marker
        traces.append(go.Scatter3d(
            x=[(xn + xx_) / 2], y=[(yn + yx_) / 2], z=[(zn + zx_) / 2],
            mode="markers",
            marker=dict(size=7, color=COLWIRE, symbol="diamond"),
            showlegend=False, hoverinfo="skip",
        ))

    return traces



# ══════════════════════════════════════════════════════════════════════════════
# UNDERGROUND DETECTOR FILTER
# ══════════════════════════════════════════════════════════════════════════════
def _ray_hits_cylinder(ox, oy, oz, dx, dy, dz,
                        ax, ay, az, bx, by, bz, r):
    """
    Does ray (o + t*d, t>=0) intersect a finite capped cylinder (A→B, radius r)?

    Tests both the curved side surface AND the two flat end caps.
    This is essential for near-vertical muons entering through horizontal caps
    of a vertically-oriented cylinder — the most common underground detector geometry.
    """
    # ── Axis unit vector ──────────────────────────────────────────────────────
    abx, aby, abz = bx-ax, by-ay, bz-az
    ab2 = abx*abx + aby*aby + abz*abz
    if ab2 < 1e-12:
        return np.zeros(len(ox), dtype=bool)
    ab_len = ab2**0.5
    abx /= ab_len; aby /= ab_len; abz /= ab_len

    hit = np.zeros(len(ox), dtype=bool)

    # ── 1. Curved side surface ────────────────────────────────────────────────
    d_dot_ab = dx*abx + dy*aby + dz*abz
    rpx = dx - d_dot_ab*abx
    rpy = dy - d_dot_ab*aby
    rpz = dz - d_dot_ab*abz

    oax = ox - ax; oay = oy - ay; oaz = oz - az
    o_dot_ab = oax*abx + oay*aby + oaz*abz
    opx = oax - o_dot_ab*abx
    opy = oay - o_dot_ab*aby
    opz = oaz - o_dot_ab*abz

    A = rpx*rpx + rpy*rpy + rpz*rpz
    B = 2*(opx*rpx + opy*rpy + opz*rpz)
    C = opx*opx + opy*opy + opz*opz - r*r

    disc = B*B - 4*A*C
    ok   = (disc >= 0) & (A > 1e-12)
    sq   = np.where(ok, np.sqrt(np.maximum(disc, 0)), 0.0)

    for sign in [1, -1]:
        # A ≈ 0 for rays parallel to the cylinder axis (e.g. vertical muons in a
        # vertical detector); those elements are masked out by `ok`, but numpy
        # still evaluates the full division first, so silence the div-by-zero.
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.where(ok, (-B + sign*sq) / (2*A), -1.0)
        t_ok = ok & (t >= 0)
        hx = oax + t*dx; hy = oay + t*dy; hz = oaz + t*dz
        proj = hx*abx + hy*aby + hz*abz
        hit |= t_ok & (proj >= 0) & (proj <= ab_len)

    # ── 2. End caps (flat disks at A and B) ───────────────────────────────────
    # Ray-plane: t = (cap_centre - origin) · n  /  (d · n)
    # where n = axis unit vector (same for both caps)
    # Hit if: t >= 0  AND  distance from cap centre <= r
    for cap_cx, cap_cy, cap_cz in [(ax, ay, az), (bx, by, bz)]:
        d_dot_n = dx*abx + dy*aby + dz*abz
        non_par = np.abs(d_dot_n) > 1e-12
        # vector from ray origin to cap centre
        ocx = cap_cx - ox; ocy = cap_cy - oy; ocz = cap_cz - oz
        with np.errstate(divide='ignore', invalid='ignore'):
            t_cap = np.where(non_par,
                             (ocx*abx + ocy*aby + ocz*abz) / d_dot_n,
                             -1.0)
        t_ok = non_par & (t_cap >= 0)
        # Point of intersection on the cap plane
        hx = ox + t_cap*dx - cap_cx
        hy = oy + t_cap*dy - cap_cy
        hz = oz + t_cap*dz - cap_cz
        dist2 = hx*hx + hy*hy + hz*hz
        hit |= t_ok & (dist2 <= r*r)

    return hit


def _ray_hits_box(ox, oy, oz, dx, dy, dz,
                   xn, xx, yn, yx, zn, zx):
    """Slab method — ray vs axis-aligned box."""
    eps = 1e-12
    t_min = np.full(len(ox), -np.inf)
    t_max = np.full(len(ox),  np.inf)

    for o_c, d_c, lo, hi in [(ox,dx,xn,xx),(oy,dy,yn,yx),(oz,dz,zn,zx)]:
        non_par = np.abs(d_c) > eps
        with np.errstate(divide='ignore', invalid='ignore'):
            t1 = np.where(non_par, (lo - o_c) / d_c, np.where(o_c >= lo, -np.inf, np.inf))
            t2 = np.where(non_par, (hi - o_c) / d_c, np.where(o_c <= hi,  np.inf,-np.inf))
        t_min = np.maximum(t_min, np.minimum(t1, t2))
        t_max = np.minimum(t_max, np.maximum(t1, t2))

    return (t_max >= t_min) & (t_max >= 0)


def apply_ug_detector_filter(df_ug, detectors):
    """
    Ray-based filter: forward + backward ray so it works whether
    overburden > or < detector depth.
      alive=1 → forward ray (overburden ≤ det depth)
                backward ray (overburden > det depth, muon passed through)
      alive=0 → position-inside check (stopped inside detector volume)
    """
    if "alive" not in df_ug.columns:
        return df_ug, 0

    all_mu = df_ug.copy()
    if len(all_mu) == 0 or not detectors:
        return all_mu.iloc[0:0], 0

    x  = all_mu["x"].values.astype(float)
    y  = all_mu["y"].values.astype(float)
    z  = all_mu["z"].values.astype(float)
    cx = all_mu["cx"].values.astype(float)
    cy = all_mu["cy"].values.astype(float)
    cz = all_mu["cz"].values.astype(float)
    alive = all_mu["alive"].values
    hit   = np.zeros(len(all_mu), dtype=bool)

    for det in detectors:
        margin = float(det.get("margin", 0.0))
        if det["shape"] == 1:   # Cylinder
            hit |= _ray_hits_cylinder(
                x, y, z,  cx,  cy,  cz,
                det["ax"], det["ay"], det["az"],
                det["bx"], det["by"], det["bz"],
                det["r"] + margin)
            hit |= _ray_hits_cylinder(          # backward ray
                x, y, z, -cx, -cy, -cz,
                det["ax"], det["ay"], det["az"],
                det["bx"], det["by"], det["bz"],
                det["r"] + margin)
        else:                   # Box
            hit |= _ray_hits_box(
                x, y, z,  cx,  cy,  cz,
                det["xmin"] - margin, det["xmax"] + margin,
                det["ymin"] - margin, det["ymax"] + margin,
                det["zmin"] - margin, det["zmax"] + margin)
            hit |= _ray_hits_box(               # backward ray
                x, y, z, -cx, -cy, -cz,
                det["xmin"] - margin, det["xmax"] + margin,
                det["ymin"] - margin, det["ymax"] + margin,
                det["zmin"] - margin, det["zmax"] + margin)

    pos_inside = np.zeros(len(all_mu), dtype=bool)
    for det in detectors:
        margin = float(det.get("margin", 0.0))
        if det["shape"] == 1:
            r_eff = det["r"] + margin
            ax, ay, az = det["ax"], det["ay"], det["az"]
            bx, by, bz = det["bx"], det["by"], det["bz"]
            abx, aby, abz = bx-ax, by-ay, bz-az
            ab_len = (abx**2 + aby**2 + abz**2)**0.5
            if ab_len < 1e-12:
                continue
            t = ((x-ax)*abx + (y-ay)*aby + (z-az)*abz) / (ab_len**2)
            perp2 = (ax+t*abx-x)**2 + (ay+t*aby-y)**2 + (az+t*abz-z)**2
            pos_inside |= (t >= 0) & (t <= 1) & (perp2 <= r_eff*r_eff)
        else:
            pos_inside |= (
                (x >= det["xmin"] - margin) & (x <= det["xmax"] + margin) &
                (y >= det["ymin"] - margin) & (y <= det["ymax"] + margin) &
                (z >= det["zmin"] - margin) & (z <= det["zmax"] + margin))

    final_hit = ((alive == 1) & hit) | ((alive == 0) & pos_inside)
    result = all_mu[final_hit].copy()
    return result, int(final_hit.sum())



def _auto_ug_filter():
    """Called after a successful MUSIC run if underground filter was requested."""
    if not st.session_state.get("ug_use_filter", False):
        return
    if st.session_state.get("ug_filter_done", False):
        return
    if _gg("music_running") or _gg("music_success") is not True:
        return

    det_list = st.session_state.get("gen_detectors", []) \
               if st.session_state.get("gen_use_detector", False) else []
    if not det_list:
        st.warning("⚠️  Underground filter requested but no detector geometry found from the Generator tab.")
        return

    src = st.session_state.get("ug_file",        "output/muons_underground.dat")
    dst = st.session_state.get("ug_filter_file", "output/muons_ug_selected.dat")

    if not Path(src).exists():
        st.warning(f"⚠️  Underground filter skipped — `{src}` not found.")
        return

    # Pre-assign so Pylance/VSCode does not warn about potentially unbound names
    nhit          = 0
    n_survived_ug = 0

    try:
        mt            = Path(src).stat().st_mtime
        dfu           = load_file(src, mtime=mt)

        # Sanity check: warn if detector is shallower than transport depth
        if "alive" in dfu.columns and det_list:
            transport_depth_cm = float(
                st.session_state.get("music_depth_m", 0)) * 100.0
            d0 = det_list[0]
            det_max_depth_cm = abs(min(
                d0.get("az", 0), d0.get("bz", 0)
            )) if d0["shape"] == 1 else abs(min(
                d0.get("zmin", 0), d0.get("zmax", 0)
            ))
            if det_max_depth_cm < transport_depth_cm * 0.9:
                st.warning(
                    f"⚠️  Depth mismatch: detector deepest face is at "
                    f"{det_max_depth_cm/100:.1f} m but MUSIC transported to "
                    f"~{transport_depth_cm/100:.1f} m. "
                    f"Set MUSIC depth ≤ {det_max_depth_cm/100:.1f} m for correct results."
                )

        dfsel, nhit   = apply_ug_detector_filter(dfu, det_list)

        # Correct parentheses: int((series).sum()) not int(series).sum()
        n_survived_ug = int((dfu["alive"] == 1).sum()) \
                        if "alive" in dfu.columns else len(dfu)

        st.session_state["music_nmuons_survived"]    = n_survived_ug
        st.session_state["music_nmuons_transported"] = len(dfu)

        header = ("# EventID  xs_cm  ys_cm  zs_cm  Es_GeV  theta_s  phi_s  charge"
                  "  alive  x_cm  y_cm  z_cm  E_GeV  cx  cy  cz  theta  phi")
        np.savetxt(dst, dfsel.values, header=header, comments="", fmt="%s")

        rate = 100.0 * nhit / max(n_survived_ug, 1)
        st.success(
            f"✅  Underground detector filter applied: **{nhit:,} / {n_survived_ug:,}** "
            f"survived muons hit the detector ({rate:.2f}%)  →  `{dst}`"
        )
        st.session_state["ug_filtered_file"] = dst
        st.session_state["ug_filter_done"]   = True

    except Exception as ex:
        st.error(f"❌  Underground filter failed: {ex}")


# ══════════════════════════════════════════════════════════════════════════════
# 3D TRAJECTORY PLOT
# ══════════════════════════════════════════════════════════════════════════════
def plot_3d_trajectories(df, n_show, depth_m, detectors=None, radius_m=800.0):
    survived = df[df["alive"] == 1].head(n_show)
    stopped  = df[df["alive"] == 0].head(n_show)
    fig = go.Figure()

    # Compute detector depth from geometry (or fall back to MUSIC depth)
    det_depth_m = depth_m
    if detectors:
        d0 = detectors[0]
        if d0["shape"] == 1:
            det_depth_m = abs(min(d0["az"], d0["bz"])) / 100.0
        else:
            det_depth_m = abs(min(d0["zmin"], d0["zmax"])) / 100.0
    above = det_depth_m < depth_m  # detector shallower than MUSIC endpoint

    # ── Survived: rock path + direction ray ───────────────────────────────────
    for idx, (_, row) in enumerate(survived.iterrows()):
        fig.add_trace(go.Scatter3d(
            x=[row["xs"] / 100, row["x"] / 100],
            y=[row["ys"] / 100, row["y"] / 100],
            z=[0, row["z"] / 100],
            mode="lines", line=dict(color="rgba(0,180,216,0.45)", width=1),
            legendgroup="survived", showlegend=(idx == 0),
            name=f"Survived ({len(survived)})",
        ))
        sign    = -1.0 if above else 1.0
        cos_th  = max(abs(float(row["cz"])), 0.01)
        ray_len = min(depth_m * 0.25, depth_m * 0.4 / cos_th)
        x_end   = row["x"] / 100 + sign * ray_len * float(row["cx"])
        y_end   = row["y"] / 100 + sign * ray_len * float(row["cy"])
        z_end   = row["z"] / 100 + sign * ray_len * abs(float(row["cz"])) * (1 if above else -1)
        fig.add_trace(go.Scatter3d(
            x=[row["x"] / 100, x_end],
            y=[row["y"] / 100, y_end],
            z=[row["z"] / 100, z_end],
            mode="lines", line=dict(color="rgba(0,230,255,0.9)", width=2),
            legendgroup="survived_ray", showlegend=(idx == 0),
            name="Direction at depth (survived)",
        ))

    # ── Stopped: rock path ────────────────────────────────────────────────────
    for idx, (_, row) in enumerate(stopped.iterrows()):
        fig.add_trace(go.Scatter3d(
            x=[row["xs"] / 100, row["x"] / 100],
            y=[row["ys"] / 100, row["y"] / 100],
            z=[0,               row["z"] / 100],
            mode="lines", line=dict(color="rgba(255,107,107,0.35)", width=1),
            legendgroup="stopped", showlegend=(idx == 0),
            name=f"Stopped ({len(stopped)})",
        ))

    # ── Generation surface ring ───────────────────────────────────────────────
    ang = np.linspace(0, 2 * np.pi, 120)
    fig.add_trace(go.Scatter3d(
        x=radius_m * np.cos(ang), y=radius_m * np.sin(ang), z=np.zeros(120),
        mode="lines", line=dict(color="rgba(0,255,150,0.9)", width=3),
        name=f"Generation surface  r={radius_m:.0f} m", legendgroup="surface",
    ))

    # ── Detector plane — real footprint from geometry ─────────────────────────
    if detectors:
        d0 = detectors[0]
        if d0["shape"] == 1:                        # Cylinder → filled disk
            rm_d  = d0["r"]  / 100.0
            cx_d  = d0["ax"] / 100.0
            cy_d  = d0["ay"] / 100.0
            ang_d = np.linspace(0, 2 * np.pi, 61)
            xd    = cx_d + rm_d * np.cos(ang_d)
            yd    = cy_d + rm_d * np.sin(ang_d)
            nd    = len(ang_d)
            xcap  = np.concatenate([[cx_d], xd])
            ycap  = np.concatenate([[cy_d], yd])
            ic_   = np.zeros(nd, dtype=int)
            jc_   = np.arange(1, nd + 1)
            kc_   = np.array([*range(2, nd + 1), 1])
            fig.add_trace(go.Mesh3d(
                x=xcap, y=ycap, z=np.full(nd + 1, -det_depth_m),
                i=ic_, j=jc_, k=kc_,
                color="rgb(255,215,0)", opacity=0.45,
                name=f"Detector plane  {det_depth_m:.0f} m",
                legendgroup="plane", showlegend=True, hoverinfo="skip",
            ))
            fig.add_trace(go.Scatter3d(
                x=xd, y=yd, z=np.full(nd, -det_depth_m),
                mode="lines", line=dict(color="rgba(255,220,50,1.0)", width=4),
                legendgroup="plane", showlegend=False, hoverinfo="skip",
            ))
        else:                                        # Box → filled rectangle
            xn_d, xx_d = d0["xmin"] / 100, d0["xmax"] / 100
            yn_d, yx_d = d0["ymin"] / 100, d0["ymax"] / 100
            fig.add_trace(go.Scatter3d(
                x=[xn_d, xx_d, xx_d, xn_d, xn_d],
                y=[yn_d, yn_d, yx_d, yx_d, yn_d],
                z=[-det_depth_m] * 5,
                mode="lines", line=dict(color="rgba(255,220,50,1.0)", width=4),
                name=f"Detector plane  {det_depth_m:.0f} m",
                legendgroup="plane", showlegend=True,
            ))
            fig.add_trace(go.Mesh3d(
                x=[xn_d, xx_d, xx_d, xn_d],
                y=[yn_d, yn_d, yx_d, yx_d],
                z=[-det_depth_m] * 4,
                i=[0, 0], j=[1, 2], k=[2, 3],
                color="rgb(255,215,0)", opacity=0.40,
                legendgroup="plane", showlegend=False, hoverinfo="skip",
            ))
    else:
        # No detector — generic bounding square fallback
        xs_range = df["xs"].values / 100
        xr = max(abs(xs_range.min()), abs(xs_range.max())) * 0.6
        fig.add_trace(go.Scatter3d(
            x=[-xr, xr, xr, -xr, -xr], y=[-xr, -xr, xr, xr, -xr],
            z=np.full(5, -det_depth_m + 5),
            mode="lines", line=dict(color="rgba(255,220,50,0.9)", width=3),
            name=f"Detector plane  {det_depth_m:.0f} m", legendgroup="plane",
        ))

    # ── Detector volume ───────────────────────────────────────────────────────
    if detectors:
        for di, det in enumerate(detectors):
            for tr in _detector_traces(det, di):
                fig.add_trace(tr)

        # z range: surface (0) on top, underground (-depth) on bottom — no autorange needed
    _z_raw    = df["z"].min() / 100.0 if "z" in df.columns else 0.0
    _z_bottom = min(_z_raw, -det_depth_m) * 1.12   # 12 % headroom below deepest point
    _z_top    = abs(_z_bottom) * 0.05               # small margin above z = 0

    fig.update_layout(
        scene=dict(
            xaxis_title="X [m]", yaxis_title="Y [m]", zaxis_title="Depth [m]",
            bgcolor="rgb(15,17,23)",
            xaxis=dict(backgroundcolor="rgb(15,17,23)", gridcolor="#333",
                       zerolinecolor="#555"),
            yaxis=dict(backgroundcolor="rgb(15,17,23)", gridcolor="#333",
                       zerolinecolor="#555"),
            zaxis=dict(backgroundcolor="rgb(15,17,23)", gridcolor="#333",
                       zerolinecolor="#555",
                       range=[_z_bottom, _z_top]),
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
        ),
        paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(15,17,23)",
        font=dict(color="white"),
        legend=dict(bgcolor="rgba(0,0,0,0.5)", bordercolor="#555",
                    borderwidth=1, font=dict(color="white")),
        margin=dict(l=0, r=0, t=30, b=0), height=620,
        title=dict(text="3D Muon Trajectories — UG Selected",
                   font=dict(color="white")),
    )
    return fig


def plot_3d_surface(df, radius_m, detectors=None, source_mode=1,
                    disk_cx=0.0, disk_cy=0.0, src_w_m=0.0,
                    disk_tilt=0.0, disk_tilt_az=0.0, source_plane=1):
    """3D plot for surface/selected muon files — generation surface + trajectories."""
    n_show = min(1500, len(df))
    sample = df.sample(n=n_show, random_state=42) if len(df) > n_show else df.copy()
    xcol   = "x" if "x" in df.columns else "xs"
    ycol   = "y" if "y" in df.columns else "ys"
    x_m    = sample[xcol].values / 100.0
    y_m    = sample[ycol].values / 100.0
    th     = sample["theta"].values
    ph     = sample["phi"].values

    # ── Z origin: read from file for hemisphere; force 0 for flat sources ──────
    # For hemisphere source_mode==3 muons start on the sphere surface (z > 0).
    # Use the actual z column if it has non-zero variance (robust to stale
    # session_state: works even if source_mode defaults to 1 incorrectly).
    if "z" in sample.columns and sample["z"].abs().max() > 1.0:
        z_m = sample["z"].values / 100.0   # cm → m
    else:
        z_m = np.zeros(len(sample))

    # Derive detector depth from geometry if available
    det_depth_m    = None
    det_z_bottom_m = None
    det_z_top_m    = None
    if detectors:
        d0 = detectors[0]
        if d0["shape"] == 1:
            # Use signed min to handle detectors at any z (positive or negative)
            _z_min_det  = min(d0["az"], d0["bz"])
            _z_max_det  = max(d0["az"], d0["bz"])
            det_z_bottom_m = _z_min_det / 100.0
            det_z_top_m    = _z_max_det / 100.0
            det_depth_m    = abs(_z_min_det) / 100.0
        else:
            _z_min_det  = min(d0["zmin"], d0["zmax"])
            _z_max_det  = max(d0["zmin"], d0["zmax"])
            det_z_bottom_m = _z_min_det / 100.0
            det_z_top_m    = _z_max_det / 100.0
            det_depth_m    = abs(_z_min_det) / 100.0

    fig = go.Figure()

    # ── Generation surface ────────────────────────────────────────────────────
    if source_mode == 3:
        # Hemisphere wireframe: latitude rings + meridian lines + equator
        for cos_t in np.linspace(0.0, 1.0, 6):
            sin_t = np.sqrt(max(0.0, 1.0 - cos_t**2))
            ang   = np.linspace(0, 2 * np.pi, 120)
            is_first = bool(cos_t == 0.0)
            fig.add_trace(go.Scatter3d(
                x=radius_m * sin_t * np.cos(ang),
                y=radius_m * sin_t * np.sin(ang),
                z=np.full(120, radius_m * cos_t),
                mode="lines",
                line=dict(color="rgba(0,255,150,0.40)", width=1),
                legendgroup="hemiring",
                name=f"Generation hemisphere (R={radius_m:.0f} m)" if is_first else None,
                showlegend=is_first,
            ))
        for phi_m in np.linspace(0, 2 * np.pi, 12, endpoint=False):
            t_arr = np.linspace(0, np.pi / 2, 60)
            fig.add_trace(go.Scatter3d(
                x=radius_m * np.sin(t_arr) * np.cos(phi_m),
                y=radius_m * np.sin(t_arr) * np.sin(phi_m),
                z=radius_m * np.cos(t_arr),
                mode="lines",
                line=dict(color="rgba(0,255,150,0.25)", width=1),
                legendgroup="hemiring",
                name=None, showlegend=False,
            ))
        # Equator ring — brighter
        ang = np.linspace(0, 2 * np.pi, 120)
        fig.add_trace(go.Scatter3d(
            x=radius_m * np.cos(ang), y=radius_m * np.sin(ang), z=np.zeros(120),
            mode="lines",
            line=dict(color="rgba(0,255,150,0.80)", width=2),
            legendgroup="hemiring",
            name=None, showlegend=False,
        ))
    elif source_mode == 2:
        # Rectangle outline: 4 corners in local frame → tilt → source-plane rotation
        _a   = float(np.radians(disk_tilt))
        _ph  = float(np.radians(disk_tilt_az))
        # half-dimensions from session state (set after run); fall back to radius_m
        _hu  = float(st.session_state.get("gen_plane_lx", radius_m))
        _hv  = float(st.session_state.get("gen_plane_ly", radius_m))
        # corners in local (U,V) frame, closed loop
        _lu = np.array([-_hu,  _hu,  _hu, -_hu, -_hu])
        _lv = np.array([-_hv, -_hv,  _hv,  _hv, -_hv])
        _rx = _lu * np.cos(_a)*np.cos(_ph) - _lv * np.sin(_ph) + disk_cx
        _ry = _lu * np.cos(_a)*np.sin(_ph) + _lv * np.cos(_ph) + disk_cy
        _rz = -_lu * np.sin(_a) + src_w_m
        if source_plane == 2:
            _rx, _ry, _rz = _rx, _rz, _ry
        elif source_plane == 3:
            _rx, _ry, _rz = _rz, _rx, _ry
        _tilt_lbl = (f"  tilt {disk_tilt:.0f}° @ {disk_tilt_az:.0f}°"
                     if disk_tilt > 0.01 else "")
        fig.add_trace(go.Scatter3d(
            x=_rx, y=_ry, z=_rz,
            mode="lines", line=dict(color="rgba(0,255,150,0.9)", width=3),
            name=f"Generation rect ({2*_hu:.0f}×{2*_hv:.0f} m{_tilt_lbl})",
            legendgroup="disk",
        ))
    else:
        # Disk ring: apply center, W-offset, tilt, and source-plane rotation
        _ang = np.linspace(0, 2 * np.pi, 180)
        _lu  = radius_m * np.cos(_ang)   # local-frame U [m]
        _lv  = radius_m * np.sin(_ang)   # local-frame V [m]
        _a   = float(np.radians(disk_tilt))
        _ph  = float(np.radians(disk_tilt_az))
        # Tangent vectors (same as Fortran sample_position):
        #   t1 = (cos α cos φ, cos α sin φ, -sin α)
        #   t2 = (-sin φ,       cos φ,        0   )
        _rx = _lu * np.cos(_a)*np.cos(_ph) - _lv * np.sin(_ph) + disk_cx
        _ry = _lu * np.cos(_a)*np.sin(_ph) + _lv * np.cos(_ph) + disk_cy
        _rz = -_lu * np.sin(_a) + src_w_m
        # Apply source-plane permutation (matches Fortran plane rotation)
        if source_plane == 2:   # XZ plane: world(x,y,z) = canonical(x, w, y)
            _rx, _ry, _rz = _rx, _rz, _ry
        elif source_plane == 3: # YZ plane: world(x,y,z) = canonical(w, x, y)
            _rx, _ry, _rz = _rz, _rx, _ry
        _tilt_lbl = (f"  tilt {disk_tilt:.0f}° @ {disk_tilt_az:.0f}°"
                     if disk_tilt > 0.01 else "")
        fig.add_trace(go.Scatter3d(
            x=_rx, y=_ry, z=_rz,
            mode="lines", line=dict(color="rgba(0,255,150,0.9)", width=3),
            name=f"Generation disk (r={radius_m:.0f} m{_tilt_lbl})",
            legendgroup="disk",
        ))

    # ── Muon origins ──────────────────────────────────────────────────────────
    # z_m is non-zero for hemisphere (actual surface position), zero for flat sources
    fig.add_trace(go.Scatter3d(
        x=x_m, y=y_m, z=z_m,
        mode="markers", marker=dict(size=1.5, color="rgba(0,180,216,0.35)"),
        name=f"Muon origins ({n_show:,})",
        legendgroup="origins",
    ))

    # ── Trajectories ──────────────────────────────────────────────────────────
    if det_depth_m is not None:
        depth_cm = det_depth_m * 100.0

        if "hit_flag" in sample.columns:
            configs = [
                (sample["hit_flag"].values == 1, "rgba(0,220,100,0.45)",
                 "Hit — toward detector", "hit"),
                (sample["hit_flag"].values == 0, "rgba(255,80,80,0.18)",
                 "Miss — does not hit",   "miss"),
            ]
        else:
            configs = [(np.ones(len(sample), dtype=bool),
                        "rgba(0,180,216,0.35)",
                        "Selected trajectory", "selected")]

        for mask, col, label, grp in configs:
            idx_arr = np.where(mask)[0]
            step    = max(1, len(idx_arr) // 400)
            shown   = False
            for i in idx_arr[::step]:
                cos_th = max(np.cos(th[i]), 0.01)
                t_p    = depth_cm / cos_th
                x_end  = x_m[i] + (t_p / 100.0) * np.sin(th[i]) * np.cos(ph[i])
                y_end  = y_m[i] + (t_p / 100.0) * np.sin(th[i]) * np.sin(ph[i])
                fig.add_trace(go.Scatter3d(
                    x=[x_m[i], x_end],
                    y=[y_m[i], y_end],
                    z=[z_m[i], det_z_bottom_m],   # signed z: correct for +/- z detectors
                    mode="lines",
                    line=dict(color=col, width=1),
                    legendgroup=grp,
                    name=label if not shown else None,
                    showlegend=not shown,
                ))
                shown = True
    else:
        # No detector — short direction lines from actual origin position
        line_len = min(radius_m * 0.4, 200.0)
        step = max(1, len(sample) // 300)
        shown = False
        for i in range(0, len(sample), step):
            fig.add_trace(go.Scatter3d(
                x=[x_m[i], x_m[i] + line_len * np.sin(th[i]) * np.cos(ph[i])],
                y=[y_m[i], y_m[i] + line_len * np.sin(th[i]) * np.sin(ph[i])],
                z=[z_m[i], z_m[i] - line_len * np.cos(th[i])],  # ← start from actual muon z
                mode="lines",
                line=dict(color="rgba(0,180,216,0.22)", width=1),
                legendgroup="directions",
                name="Muon directions" if not shown else None,
                showlegend=not shown,
            ))
            shown = True

    # ── Detector geometry ─────────────────────────────────────────────────────
    if detectors:
        for di, det in enumerate(detectors):
            for tr in _detector_traces(det, di):
                fig.add_trace(tr)

    title = ("Surface Muon Trajectories — hits (green) vs misses (red)"
             if "hit_flag" in df.columns
             else "Selected Muon Trajectories → Detector")

    fig.update_layout(
        scene=dict(
            xaxis_title="X [m]", yaxis_title="Y [m]", zaxis_title="Z [m]",
            bgcolor="rgb(15,17,23)",
            xaxis=dict(backgroundcolor="rgb(15,17,23)",
                       gridcolor="#333", zerolinecolor="#555"),
            yaxis=dict(backgroundcolor="rgb(15,17,23)",
                       gridcolor="#333", zerolinecolor="#555"),
            zaxis=dict(backgroundcolor="rgb(15,17,23)",
                       gridcolor="#333", zerolinecolor="#555"),
            camera=dict(eye=dict(x=1.4, y=1.4, z=0.8)),
        ),
        paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(15,17,23)",
        font=dict(color="white"),
        legend=dict(bgcolor="rgba(0,0,0,0.5)", bordercolor="#555",
                    borderwidth=1, font=dict(color="white")),
        margin=dict(l=0, r=0, t=30, b=0), height=600,
        title=dict(text=title, font=dict(color="white")),
    )
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# SURVIVAL RATE vs DEPTH
# ══════════════════════════════════════════════════════════════════════════════
def plot_survival_vs_depth(df, depth_m, rho):
    depths_m    = np.linspace(0, depth_m * 1.5, 120)
    depths_gcm2 = depths_m * 100.0 * rho

    # Groom 2001 CSDA table (log-log interpolation) — replaces wrong linear approx.
    # depths_gcm2 already carries the material density (opacity equivalence);
    # the previous rho_stdrock/rho factor cancelled rho entirely, so the curve
    # was always evaluated at Standard Rock density.
    E_min_GeV = np.array([
        _groom_threshold_energy(x)[0] for x in depths_gcm2
    ])

    # Es is TOTAL energy; the Groom threshold is kinetic
    energies = df["Es"].values
    total    = max(len(df), 1)
    rates    = [100.0 * np.sum(energies > e + MUON_MASS_GEV) / total
                for e in E_min_GeV]

    # Bug 2 fix: prefer stored MUSIC result over re-computing from possibly
    # filtered df (which would give 100% if ug_selected file is loaded)
    n_survived    = st.session_state.get("music_nmuons_survived",    None)
    n_transported = st.session_state.get("music_nmuons_transported", None)
    if n_survived is not None and n_transported and n_transported > 0:
        actual_rate = 100.0 * n_survived / n_transported
    else:
        actual_rate = 100.0 * (df["alive"] == 1).sum() / total


    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=depths_m, y=rates, mode="lines",
        line=dict(color="#00b4d8", width=2.5),
        name="Estimated survival rate",
        fill="tozeroy", fillcolor="rgba(0,180,216,0.08)",
    ))
    fig.add_vline(
        x=depth_m, line=dict(color="#ffd700", width=2, dash="dash"),
        annotation_text=f"  MUSIC depth: {depth_m} m<br>  Rate: {actual_rate:.1f}%",
        annotation_font=dict(color="#ffd700", size=12),
        annotation_position="top right",
    )
    fig.add_trace(go.Scatter(
        x=[depth_m], y=[actual_rate], mode="markers",
        marker=dict(color="#ffd700", size=12, symbol="diamond",
                    line=dict(color="white", width=1.5)),
        name=f"MUSIC Monte Carlo result: {actual_rate:.1f}% at {depth_m} m",
    ))
    fig.update_layout(
        xaxis_title="Rock depth [m]",
        yaxis_title="Fraction of input muons surviving [%]",
        yaxis=dict(range=[0, 105]),
        paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
        font=dict(color="white"),
        xaxis=dict(gridcolor="#2a2a3a", zerolinecolor="#444"),
        yaxis_gridcolor="#2a2a3a",
        legend=dict(bgcolor="rgba(0,0,0,0.5)", bordercolor="#555", borderwidth=1,
                    font=dict(color="white")),
        margin=dict(l=60, r=30, t=80, b=50), height=440,
        title=dict(
            text=(
                "CSDA Analytical Estimate vs MUSIC Monte Carlo Result<br>"
                "<sup>Blue curve: Groom (2001) CSDA range table — "
                "fraction of input muons with E_surface > E_min(depth) | "
                "◆ Actual MUSIC survival at transport depth</sup>"
            ),
            font=dict(color="white", size=14),
        ),
        annotations=[dict(
            text=(
                "⚠️ CSDA neglects straggling & radiative fluctuations — "
                "gap between curve and ◆ is physically expected"
            ),
            xref="paper", yref="paper",
            x=0.01, y=0.1,
            showarrow=False,
            font=dict(color="rgba(200,200,200,0.7)", size=10),
            bgcolor="rgba(0,0,0,0.4)",
            bordercolor="#555", borderwidth=1,
        )],
    )
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
st.title(f"🌌 UCMuon  ·  v{__version__}")
st.caption(
    f"🌌 **UCMuon (You See Muon)!** v{__version__} — UCLouvain Muography Group | "
    "Hamid Basiri · [hamid.basiri@uclouvain.be](mailto:hamid.basiri@uclouvain.be) | "
    "MIT License 2026"
)
st.divider()

def compute_detector_solid_angle(detectors, origin_cm=(0.0, 0.0, 0.0),
                                  n_samples=600_000):
    """
    Monte Carlo solid angle of detector(s) as seen from origin_cm.
    Shoots uniform rays over the downward hemisphere (z < 0).

    Returns
    -------
    solid_angle_sr   : geometric solid angle [sr]
    cos2_acceptance  : cos²θ-weighted acceptance [sr]  ← used for muon rate
    solid_angle_msr  : geometric solid angle [msr]
    frac_hemisphere  : fraction of 2π hemisphere subtended
    """
    if not detectors:
        return 0.0, 0.0, 0.0, 0.0

    rng       = np.random.default_rng(42)
    cos_theta = rng.uniform(0.0, 1.0, n_samples)          # downward hemisphere
    phi       = rng.uniform(0.0, 2 * np.pi, n_samples)
    sin_theta = np.sqrt(1.0 - cos_theta ** 2)

    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = -cos_theta                                        # always downward

    ox = np.full(n_samples, float(origin_cm[0]))
    oy = np.full(n_samples, float(origin_cm[1]))
    oz = np.full(n_samples, float(origin_cm[2]))

    hit = np.zeros(n_samples, dtype=bool)

    for det in detectors:
        margin = float(det.get("margin", 0.0))
        if det["shape"] == 1:   # Cylinder
            hit |= _ray_hits_cylinder(
                ox, oy, oz, dx, dy, dz,
                det["ax"], det["ay"], det["az"],
                det["bx"], det["by"], det["bz"],
                det["r"] + margin,
            )
        else:                   # Box
            hit |= _ray_hits_box(
                ox, oy, oz, dx, dy, dz,
                det["xmin"] - margin, det["xmax"] + margin,
                det["ymin"] - margin, det["ymax"] + margin,
                det["zmin"] - margin, det["zmax"] + margin,
            )

    n_hit            = hit.sum()
    solid_angle_sr   = (n_hit / n_samples) * 2 * np.pi
    cos2_acceptance  = (cos_theta[hit] ** 2).sum() / n_samples * 2 * np.pi
    solid_angle_msr  = solid_angle_sr * 1e3
    frac_hemisphere  = n_hit / n_samples

    return solid_angle_sr, cos2_acceptance, solid_angle_msr, frac_hemisphere


def _render_mcs_margin_helper(detectors=None):
    """Standalone Highland MCS displacement / safety-margin calculator."""
    try:
        from gui_source_optimizer import _sr, _groom, _csda_range_gcm2
    except ImportError:
        st.info("gui_source_optimizer module not available.")
        return

    st.caption(
        "Lateral displacement of a muon after multiple Coulomb scattering in Standard Rock "
        "(ρ = 2.65 g/cm³, X₀ = 26.7 g/cm², Highland formula). Use the **suggested margin** "
        "for the detector *Safety margin* field."
    )

    # Pre-fill from the first detector, if one is defined
    _def_depth, _def_r = 100.0, 5.0
    if detectors:
        _d0 = detectors[0]
        if _d0["shape"] == 1:
            _def_depth = abs(min(float(_d0["az"]), float(_d0["bz"]))) / 100.0
            _def_r     = float(_d0["r"])
        else:
            _def_depth = abs(min(float(_d0["zmin"]), float(_d0["zmax"]))) / 100.0
            _def_r     = 0.5 * max(abs(float(_d0["xmax"]) - float(_d0["xmin"])),
                                   abs(float(_d0["ymax"]) - float(_d0["ymin"])))
        _def_depth = max(_def_depth, 0.1)

    _hc1, _hc2, _hc3 = st.columns(3)
    _depth = _hc1.number_input("Depth [m]", 0.1, 5000.0, _def_depth, 10.0,
                               key="hlp_mcs_depth",
                               help="Vertical rock thickness above the detector.")
    _zen   = _hc2.number_input("Zenith θ [°]", 0.0, 89.0,
                               float(st.session_state.get("thetamax", 85.0)), 5.0,
                               key="hlp_mcs_ze",
                               help="Slant path = depth / cos θ. Use θ_max for the worst case.")
    _r_det = _hc3.number_input("Detector radius [cm]", 0.1, 1e4, _def_r, 1.0,
                               key="hlp_mcs_rdet",
                               help="Lateral half-size of the detector.")

    _e_min   = float(st.session_state.get("emin", 1.0))
    _cos_ze  = max(math.cos(math.radians(_zen)), 0.01)
    _slant_m = _depth / _cos_ze
    _E_thr   = _groom(_slant_m * 100.0 * 2.65)   # min KE to traverse the slant path
    _E_eff   = max(_e_min, _E_thr)               # slowest muon that actually arrives

    _sig1 = _sr(_E_eff, _zen, _depth * 100.0)    # 1σ radial displacement [cm]
    _hm1, _hm2, _hm3 = st.columns(3)
    _hm1.metric("1σ displacement",  f"{_sig1:.1f} cm",
                help=f"σ_r at E = {_E_eff:.0f} GeV over the {_slant_m:.0f} m slant path.")
    _hm2.metric("Suggested margin", f"{1.5 * _sig1:.1f} cm",
                help="1.5 × σ_r — same convention as the source-size optimiser.")
    _hm3.metric("3σ displacement",  f"{3.0 * _sig1:.1f} cm",
                help="99.7 % coverage — compare with your detector size.")

    if _e_min < _E_thr:
        _stop_m = _csda_range_gcm2(_e_min) / 2.65 / 100.0
        st.warning(
            f"⚠️ E_min = {_e_min:.1f} GeV muons stop after ≈ {_stop_m:.1f} m of rock and never "
            f"reach {_depth:.0f} m at θ = {_zen:.0f}° (needs ≥ {_E_thr:.0f} GeV). "
            f"Displacement above is evaluated at {_E_thr:.0f} GeV — the slowest muon that arrives."
        )

    _eff = (_r_det / max(_r_det + _sig1, 1e-9)) ** 2 * 100.0
    if _eff < 5.0:
        st.warning(
            f"⚠️ 1σ scatter (**{_sig1:.0f} cm**) ≫ detector radius (**{_r_det:.1f} cm**) — the "
            f"detector pre-filter keeps < {_eff:.2f} % of useful muons. Disable the filter and "
            "generate over the full source area instead."
        )
    else:
        st.caption(f"Filter efficiency with this margin ≈ {_eff:.0f} %.")



def _compute_flux(nhits):
    """
    Returns (rate_per_s, None, valid).
    Uses cos²θ-weighted MC acceptance if detector geometry is available,
    otherwise falls back to analytical cone formula.
    """
    phi   = st.session_state.get("gen_integrated_flux", None)
    r_m   = st.session_state.get("gen_radius",          None)
    n_gen = st.session_state.get("gen_nmuons_done",     None)

    if any(v is None for v in [phi, r_m, n_gen]) or n_gen == 0:
        return None, None, False

    rm         = st.session_state.get("gen_radius", None)
    if rm is None:
        return None, None, False
    r_cm       = float(rm) * 100.0
    src_mode   = st.session_state.get("gen_source_mode", 1)
    if src_mode == 2:
        lx_cm  = float(st.session_state.get("gen_plane_lx", 0)) * 100.0
        ly_cm  = float(st.session_state.get("gen_plane_ly", 0)) * 100.0
        A_disk = 4.0 * lx_cm * ly_cm
    elif src_mode == 3:
        A_disk = 2.0 * np.pi * r_cm**2   # hemisphere surface area = 2πR²
    else:
        A_disk = np.pi * r_cm**2
    # Prefer MC acceptance if detector is defined
    det_list = st.session_state.get("gen_detectors", []) \
               if st.session_state.get("gen_use_detector", False) else []

    if det_list:
        _, cos2_acceptance, _, _ = compute_detector_solid_angle(det_list)
        Omega = cos2_acceptance if cos2_acceptance > 0 else 1e-12
    else:
        # Analytical fallback: full cone up to theta_max
        th_max = np.radians(float(st.session_state.get("gen_theta_max", 85.0)))
        Omega  = 2 * np.pi / 3.0 * (1 - np.cos(th_max) ** 3)

    weight     = phi * A_disk * Omega / float(n_gen)
    rate_per_s = nhits * weight
    return rate_per_s, None, True




tab_gen, tab_music, tab_terrain, tab_results, tab_density, tab_config = st.tabs([
    "🌌  Generator",
    "🪨  Transport",
    "🗺  Terrain",
    "📊  Results",
    "🔬  Density",
    "📋  Config",
])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — SURFACE GENERATOR
# ══════════════════════════════════════════════════════════════════════════════
with tab_gen:

    # ── Workflow selector ─────────────────────────────────────────────────────
    gen_workflow = st.radio(
        "Workflow",
        ["Standard", "DAS-REM"],
        format_func=lambda x: {
            "Standard": "🔵  Standard  —  forward generation",
            "DAS-REM":  "🔴  Guaranteed-hit mode  —  100 % detector hits",
        }[x],
        horizontal=True,
        key="gen_workflow",
        help=(
            "**Standard:** uniform sampling on the source surface, optional detector filter. "
            "N muons = total generated (acceptance < 100 %).\n\n"
            "**Guaranteed-hit mode:** samples the muon hit point uniformly on the detector face "
            "first, then back-projects to the source surface via straight-line geometry — every "
            "muon is guaranteed to reach the detector (no wasted trials). "
            "This reverse-sampling strategy follows the approach of "
            "Yao et al., *J. Appl. Phys.* **138**, 144901 (2025)."
        ),
    )
    use_dasrem = (gen_workflow == "DAS-REM")
    st.divider()

    # ══════════════════════════════════════════════════════════════════════════
    # ROW 1 — Physics (left) | Right column (conditional on workflow)
    # ══════════════════════════════════════════════════════════════════════════
    _col_l, _col_r = st.columns([4, 6])

    # ─────────────────────────── LEFT — Physics ───────────────────────────────
    with _col_l:
        st.markdown("#### Physics")

        # ── Mono-energetic beam toggle ────────────────────────────────────────
        mono_beam = st.checkbox(
            "⚡  Mono-energetic beam  (single energy E₀)",
            value=False, key="gen_mono",
            help="Every particle is generated at exactly E₀ (total energy). "
                 "The spectrum model is bypassed (delta-function energy); "
                 "the angular distribution below still applies. "
                 "Useful for detector-response and transport validation studies.",
        )

        # ── Spectrum model (compact selectbox) ────────────────────────────────
        spectrum_mode = st.selectbox(
            "Spectrum model",
            [1, 2, 3, 4, 5, 6, 7, 8],
            format_func=lambda x: {
                1: "① CosmoALEPH   p⁻³·¹⁹⁵²           ~100–2 500 GeV  (default — thick targets)",
                2: "② Power-law    E⁻³·⁷               legacy MUSIC cross-check",
                3: "③ PARMA/EXPACS                     0.1 GeV–100 TeV  (location & date-aware)",
                4: "④ Guan (2015)  a=3.64, b=1.29      > 10 GeV surface / any depth underground",
                5: "⑤ Frosin (2025) a=3.512, b=1.388   > 10 GeV surface / any depth underground",
                6: "⑥ Bugaev (1998) piecewise poly      1–1 000 GeV",
                7: "⑦ Reyna–Bugaev (2006) cos³θ        1–10 000 GeV  (best surface estimate)",
                8: "⑧ Cosmic electrons  E⁻³·⁰           10 MeV–1 GeV  (surface/shallow)",
            }[x],
            help=(
                "**① CosmoALEPH  —  ~100–2 500 GeV/c:**  "
                "Power-law fit dN/dp ∝ p⁻³·¹⁹⁵² anchored to the sea-level vertical "
                "muon spectrum measured by CosmoALEPH with the ALEPH detector at LEP "
                "(Schmelling et al. 2013). Reproduces the measurement to ~5% over "
                "112–2 239 GeV/c; being a pure power law it overestimates the flux "
                "at low momenta (×2.5 at 10 GeV/c). Use for thick targets where the "
                "detected flux is dominated by muons above ~50 GeV.\n\n"
                "**② Power-law  —  legacy MUSIC cross-check:**  "
                "Simple dN/dE ∝ E⁻³·⁷. A sampling shape, not an absolute flux model — "
                "provided for cross-checks with legacy MUSIC simulations and for "
                "underground studies where the spectrum above the transport threshold "
                "is approximately power-law.\n\n"
                "**③ PARMA/EXPACS  —  0.1 GeV–100 TeV:**  "
                "Full location- and date-aware spectrum with geomagnetic rigidity cutoff and "
                "correct μ⁺/μ⁻ ratio. Requires the PHITS PARMA data directory.\n\n"
                "**④ Guan (2015)  —  > 10 GeV surface, any depth underground:**  "
                "Modified Gaisser (arXiv:1509.06176) with E_eff = E + 3.64/cos^1.29(θ*) "
                "accounting for atmospheric energy loss. "
                "Use for underground detector acceptance. "
                "At the surface with E < 10 GeV the atmospheric correction suppresses "
                "the flux unrealistically — use Reyna-Bugaev for surface rate estimation.\n\n"
                "**⑤ Frosin (2025)  —  > 10 GeV surface, any depth underground:**  "
                "Same Guan formula re-fitted on 304 sea-level datasets "
                "(J.Phys.G 52, 035002). Same atmospheric-loss limitation at low energies.\n\n"
                "**⑥ Bugaev (1998)  —  1–1 000 GeV:**  "
                "Gaisser (1990) pion+kaon formula with Bugaev normalisation. "
                "Pair with cos²θ angular mode. Best range 1–1 000 GeV.\n\n"
                "**⑦ Reyna–Bugaev (2006)  —  1–10 000 GeV:**  "
                "Log-polynomial fit to p³·I_vert validated against PDG surface intensity "
                "(I_V > 1 GeV ≈ 70 m⁻²sr⁻¹s⁻¹, ~20% agreement). "
                "Best model for surface flux and measurement-time estimates. "
                "Pair with cos³θ angular mode.\n\n"
                "**⑧ Cosmic electrons  —  10 MeV–1 GeV:**  "
                "Primary atmospheric e⁺/e⁻ with dN/dE ∝ E⁻³·⁰. "
                "Generates electrons (PDG 11/−11) with equal e⁺/e⁻ ratio. "
                "Electrons are stopped by <1 m of rock — only relevant for surface "
                "or very shallow detectors (<1 m.w.e.). Pair with cos³θ angular mode."
            ),
            key="gen_spectrum_mode",
            disabled=mono_beam,
        )
        # reference line shown below the energy inputs — includes range guide
        _spec_refs = {
            1: ("📚 CosmoALEPH fit, Schmelling et al. (2013)  dN/dp ∝ p⁻³·¹⁹⁵²  "
                "| valid range: **~100–2 500 GeV/c**  (overestimates ×2.5 at 10 GeV/c)"),
            2: ("📚 Power-law  dN/dE ∝ E⁻³·⁷  "
                "| sampling shape for legacy MUSIC cross-checks — not an absolute flux model"),
            3: ("📚 Sato et al., PARMA/EXPACS  "
                "| valid range: **0.1 GeV–100 TeV**  (geomagnetic + solar modulation)"),
            4: ("📚 Guan et al. (2015), arXiv:1509.06176  a=3.64, b=1.29  "
                "| valid range: **> 10 GeV** surface  /  any depth underground"),
            5: ("📚 Frosin et al. (2025), J.Phys.G 52, 035002  a=3.512, b=1.388  "
                "| valid range: **> 10 GeV** surface  /  any depth underground"),
            6: ("📚 Bugaev et al. (1998) / Gaisser (1990)  pion+kaon terms  "
                "| valid range: **1–1 000 GeV**  — pair with angular mode ② cos²θ"),
            7: ("📚 Reyna (2006) / Bugaev (1998)  log-polynomial in p  "
                "| valid range: **1–10 000 GeV**  — pair with angular mode ⑤ cos³θ"),
            8: ("⚡ Cosmic electrons  dN/dE ∝ E⁻³·⁰  "
                "| valid range: **10 MeV–1 GeV**  — generates e⁺/e⁻, not muons  "
                "— pair with angular mode ⑤ cos³θ"),
        }

        # ── Energy range ──────────────────────────────────────────────────────
        if mono_beam:
            st.markdown("**Beam energy**")
            st.session_state.setdefault("gen_mono_e", 100.0)
            mono_E0 = st.number_input("E₀ [GeV]  (total energy)",
                                      min_value=0.001, max_value=100000.0,
                                      step=10.0, format="%.3f", key="gen_mono_e")
            e_min, e_max = float(mono_E0), float(mono_E0)
            # delta spectrum: keep electron mode if selected, otherwise use the
            # analytical power-law sampler (safe with a degenerate energy window)
            spectrum_mode = 8 if spectrum_mode == 8 else 2
            if spectrum_mode != 8 and mono_E0 <= 0.106:
                st.warning("⚠️  E₀ must exceed the muon rest mass (0.106 GeV).")
            st.caption(f"⚡ All particles at E₀ = {mono_E0:g} GeV — spectrum model bypassed.")
        else:
            st.markdown("**Energy range**")
            _ec1, _ec2 = st.columns(2)
            st.session_state.setdefault("emin", 1.0)
            st.session_state.setdefault("emax", 2500.0)
            e_min = _ec1.number_input("E min [GeV]", 0.01, 10000.0, step=0.1,
                                       format="%.2f", key="emin")
            e_max = _ec2.number_input("E max [GeV]", 0.1, 100000.0, step=100.0,
                                       format="%.0f", key="emax")
            if spectrum_mode in _spec_refs:
                st.caption(f"{_spec_refs[spectrum_mode]}  |  your range: E ∈ [{e_min:.2f}, {e_max:.0f}] GeV")
            st.caption("💡 CSDA floor & muography window for your detector depth → "
                       "[🧰 Helpers & calculators](#helpers)")
        # Out-of-range warnings per model
        _erange_warns = {
            1: (e_min < 100.0 or e_max > 2500,  "CosmoALEPH fit reproduces the measurement only in ~100–2 500 GeV/c; below that it overestimates the flux (×2.5 at 10 GeV/c) — for shallow targets prefer ④ Guan, ⑤ Frosin or ⑦ Reyna–Bugaev."),
            2: (e_min < 100.0,                   "Power-law E⁻³·⁷ is a sampling shape, not an absolute flux model — below ~100 GeV the real spectrum is much flatter."),
            4: (e_min < 10.0,                    "Guan (2015) suppresses surface flux by up to 50× below 10 GeV due to the atmospheric energy-loss correction."),
            5: (e_min < 10.0,                    "Frosin (2025) has the same atmospheric correction as Guan — surface flux unreliable below 10 GeV."),
            8: (e_max > 1.0,                     "Cosmic electrons are stopped by <1 m of rock above ~1 GeV. E max above 1 GeV is not physically meaningful for surface detectors."),
        }
        if not mono_beam and spectrum_mode in _erange_warns:
            _warn_cond, _warn_msg = _erange_warns[spectrum_mode]
            if _warn_cond:
                st.warning(f"⚠️  {_warn_msg}")
        if spectrum_mode == 8:
            st.info("⚡ **Electron mode:** output file `charge` column is ±1 (e⁺/e⁻). "
                    "PHITS dump uses PDG codes 11/−11. Angular mode ⑤ cos³θ recommended.")
        elif e_min < 1.0:
            st.caption("ℹ️  Sub-GeV muons are available via **guaranteed-hit mode** (workflow toggle at the top).")
        if use_dasrem and not mono_beam and spectrum_mode in (3, 8):
            _dr_unsup = "③ PARMA/EXPACS" if spectrum_mode == 3 else "⑧ Cosmic electrons"
            st.warning(
                f"⚠️  Spectrum {_dr_unsup} is **not available in guaranteed-hit mode** — "
                "the generator would fall back to ① CosmoALEPH. "
                "Pick another spectrum model or switch to the Standard workflow."
            )

        # ── PARMA (only when mode 3, Standard workflow) ───────────────────────
        if spectrum_mode == 3 and not use_dasrem:
            with st.expander("🌍  PARMA location & date", expanded=True):
                st.caption("Accounts for geomagnetic cutoff, atmospheric depth, solar modulation, and μ⁺/μ⁻ ratio.")
                _pc1, _pc2, _pc3 = st.columns(3)
                parma_lat = _pc1.number_input("Lat [°]", -90.0, 90.0, 50.7, 0.1, key="parma_lat",
                                              help="Positive = North. e.g. 50.7° for Brussels.")
                parma_lon = _pc2.number_input("Lon [°]", -180.0, 180.0, 4.4, 0.1, key="parma_lon",
                                              help="Positive = East. e.g. 4.4° for Brussels.")
                parma_alt = _pc3.number_input("Alt [km]", 0.0, 10.0, 0.0, 0.1, key="parma_alt")
                _pd1, _pd2, _pd3 = st.columns(3)
                import datetime as _dt
                _today = _dt.date.today()
                parma_year  = _pd1.number_input("Year",  2000, 2050, _today.year,  1, key="parma_year")
                parma_month = _pd2.number_input("Month", 1,    12,   _today.month, 1, key="parma_month")
                parma_day   = _pd3.number_input("Day",   1,    31,   _today.day,   1, key="parma_day")
                _pw1, _pw2 = st.columns(2)
                parma_charge = _pw1.radio("Charge", [0, 1, -1],
                                         format_func=lambda x: {0:"μ⁺+μ⁻", 1:"μ⁺", -1:"μ⁻"}[x],
                                         horizontal=True, key="parma_charge")
                parma_sw = _pw2.number_input("W (Wolf) index", -135.0, 300.0, 0.0, 5.0, key="parma_sw",
                                             help="0 = solar min, ~150 = solar max.")
                with st.expander("Citation requirement (EXPACS/PARMA)", expanded=False):
                    st.caption(
                        "Publishing results from this mode requires citing both papers below "
                        "and acknowledging http://phits.jaea.go.jp/expacs. "
                        "Send a copy of your publication to nsed-expacs@jaea.go.jp (EXPACS terms of use).\n\n"
                        "- T. Sato, *PLoS ONE* **10**(12): e0144679 (2015) — "
                        "[doi:10.1371/journal.pone.0144679](https://doi.org/10.1371/journal.pone.0144679)\n"
                        "- T. Sato, *PLoS ONE* **11**(8): e0160390 (2016) — "
                        "[doi:10.1371/journal.pone.0160390](https://doi.org/10.1371/journal.pone.0160390)"
                    )
                _plines = _gl("gen_lines")
                _diag   = [l for l in _plines if any(k in l for k in
                           ("atm. depth","cutoff rigid","FFP","W index","getHP","mu+ fraction"))]
                if _diag:
                    st.divider()
                    st.caption("**Last PARMA diagnostics:**")
                    st.code("\n".join(_diag), language="text")
        else:
            parma_lat = 50.7; parma_lon = 4.4; parma_alt = 0.0
            parma_year = 2026; parma_month = 1; parma_day = 20
            parma_charge = 0; parma_sw = 0.0

        st.divider()

        # ── Angular distribution ──────────────────────────────────────────────
        st.markdown("#### Angular distribution")
        _is_parma = (not mono_beam) and st.session_state.get("gen_spectrum_mode", 1) == 3
        st.session_state.setdefault("angularmode", 2)
        angular_mode = st.selectbox(
            "Mode",
            [1, 2, 3, 4, 5],
            format_func=lambda x: {
                1: "① Vertical only  (θ = 0)",
                2: ("② PARMA angular distribution  (recommended)"
                    if _is_parma else
                    "② cos²θ  — realistic (recommended)"),
                3: "③ Uniform cone",
                4: "④ Guan/Frosin  P(θ|E) — self-consistent",
                5: "⑤ cos³θ  — Reyna–Bugaev",
            }[x],
            help=(
                "**①** Pencil beam, θ=0. Quick geometry/acceptance checks.\n\n"
                "**② (non-PARMA spectra)** Φ ∝ cos²θ — standard empirical sea-level distribution.\n\n"
                "**② (PARMA spectrum)** Energy-averaged PARMA zenith angle distribution: "
                "P(cosθ) ∝ ∫ Φ(E)·F_ang(E,cosθ) dE. "
                "Accounts for location, altitude, and geomagnetic rigidity (Sato 2016).\n\n"
                "**③** Uniform in solid angle within [0°, θ_max]. Acceptance mapping.\n\n"
                "**④** Samples θ from P(θ|E) ∝ F(E,θ)·cosθ using the exact energy sampled. "
                "Physically exact for spectrum modes ④ and ⑤ — Guan/Frosin couple E and θ.\n\n"
                "**⑤** cos³θ analytical. Inverse CDF: cosθ=(1−u·(1−cos⁴θ_max))^(1/4). "
                "Designed for Reyna–Bugaev (spectrum ⑦)."
            ),
            key="angularmode",
        )
        theta_max = 85.0
        if angular_mode in [2, 3, 4, 5]:
            theta_max = st.slider("Max zenith angle θ_max [°]", 10.0, 89.0, 85.0, 1.0, key="thetamax")
            st.caption("💡 Recommended θ_max for your detector → "
                       "[🧰 Helpers & calculators](#helpers)")
        if _is_parma and angular_mode == 2:
            st.info(
                "**PARMA mode:** option ② samples from the physics-based PARMA zenith "
                "angle distribution (energy-averaged, Sato 2016), not cos²θ. "
                "This is the physically consistent choice for PARMA/EXPACS."
            )
        if _is_parma and angular_mode in [4, 5]:
            st.warning(
                "⚠️ With PARMA spectrum, modes ④ and ⑤ are not defined. "
                "Select ② for PARMA's own angular distribution, ① for vertical, or ③ for uniform cone."
            )
        _spec_now = 0 if mono_beam else st.session_state.get("gen_spectrum_mode", 1)
        if angular_mode == 4 and _spec_now not in [0, 4, 5]:
            if not _is_parma:
                st.warning("⚠️ Mode ④ is only physically meaningful with spectrum modes ④ or ⑤.")
        if angular_mode == 5 and _spec_now not in [0, 7, 8]:
            st.warning("⚠️ Mode ⑤ (cos³θ) is designed for spectrum ⑦ (Reyna–Bugaev) or ⑧ (cosmic electrons).")
        if _spec_now == 6 and angular_mode not in [1, 2]:
            st.info("ℹ️ Bugaev/Gaisser (⑥) uses a fixed vertical spectrum — pair with angular mode ② cos²θ for best consistency.")
        if _spec_now == 7 and angular_mode != 5:
            st.info("ℹ️ Reyna–Bugaev (⑦) is calibrated with cos³θ angular distribution — select angular mode ⑤.")
        if use_dasrem and angular_mode == 4:
            st.warning(
                "⚠️  Angular mode ④ P(θ|E) is **not available in guaranteed-hit mode** — "
                "a uniform cone would be used instead. Choose ②, ③ or ⑤."
            )

    # ────────────────────── RIGHT — conditional on workflow ──────────────────
    with _col_r:
        if use_dasrem:
            # ── DAS-REM: Detector is the starting point ────────────────────────
            st.markdown("#### Detector")
            st.caption(
                "Hit point sampled on the detector face first; source position computed "
                "by reverse geometry — every muon is guaranteed to arrive. "
                "The physics still propagates forward."
            )
            output_sel = "output/muons_selected.dat"
            _dr_f1, _dr_f2 = st.columns([3, 2])
            output_sel = _dr_f1.text_input("Output file", "output/muons_selected.dat", key="outputsel")
            _dr_f2.caption("ℹ️  Reverse mapping targets a **single detector** — "
                           "use the Standard workflow for multi-detector filtering.")
            ndet = 1
            detectors = []
            _det_cols_r = st.columns(1)
            for i in range(ndet):
                with _det_cols_r[i % len(_det_cols_r)]:
                    with st.expander(f"Detector {i+1}", expanded=(i == 0)):
                        shape  = st.selectbox("Shape", [1, 2], key=f"sh{i}",
                                     format_func=lambda x: "Cylinder" if x == 1 else "Box (AABB)")
                        margin = st.number_input(
                            "Safety margin [cm]", 0.0, 500.0, 0.0, key=f"mg{i}",
                            help="Acceptance halo added around the detector to catch "
                                 "muons scattered back in by MCS.")
                        st.caption("💡 Margin sizing → [🧰 Helpers & calculators](#helpers) "
                                   "→ 🎯 MCS margin")
                        d = {"shape": shape, "margin": margin}
                        if shape == 1:
                            _dc1, _dc2 = st.columns(2)
                            with _dc1:
                                st.markdown("**Bottom A**")
                                d["ax"] = st.number_input("Ax [cm]", value=0.0,     key=f"ax{i}")
                                d["ay"] = st.number_input("Ay [cm]", value=0.0,     key=f"ay{i}")
                                d["az"] = st.number_input("Az [cm]", value=-9000.0, key=f"az{i}")
                            with _dc2:
                                st.markdown("**Top B**")
                                d["bx"] = st.number_input("Bx [cm]", value=0.0, key=f"bx{i}")
                                d["by"] = st.number_input("By [cm]", value=0.0, key=f"by{i}")
                                d["bz"] = st.number_input("Bz [cm]", value=0.0, key=f"bz{i}")
                            d["r"] = st.number_input("Radius [cm]", 0.1, 1e4, 5.0, key=f"rr{i}")
                        else:
                            _dc1, _dc2 = st.columns(2)
                            with _dc1:
                                d["xmin"] = st.number_input("Xmin [cm]", value=-100.0,  key=f"xn{i}")
                                d["ymin"] = st.number_input("Ymin [cm]", value=-100.0,  key=f"yn{i}")
                                d["zmin"] = st.number_input("Zmin [cm]", value=-9000.0, key=f"zn{i}")
                            with _dc2:
                                d["xmax"] = st.number_input("Xmax [cm]", value=100.0, key=f"xx{i}")
                                d["ymax"] = st.number_input("Ymax [cm]", value=100.0, key=f"yx{i}")
                                d["zmax"] = st.number_input("Zmax [cm]", value=0.0,   key=f"zx{i}")
                        detectors.append(d)
            use_detector = True
        else:
            # ── Standard: Generation surface + Sampling ────────────────────────
            st.markdown("#### Generation surface")

            # ── Shape + plane in one row ──────────────────────────────────────
            _shape_labels = {1: "💿 Circular disk", 2: "▭  Rectangle", 3: "🌐 Hemisphere"}
            _sg_shape, _sg_plane = st.columns([3, 4])
            with _sg_shape:
                st.markdown("**Shape**")
                source_mode = st.radio("Shape", [1, 2, 3],
                                       format_func=lambda x: _shape_labels[x],
                                       horizontal=False,
                                       key="source_mode",
                                       label_visibility="collapsed")

            if source_mode == 3:
                source_plane = 1
                with _sg_plane:
                    st.markdown("**Parameters**")
                    hemi_radius = st.number_input("Radius [m]", 0.01, 10000.0, 200.0, 1.0,
                                                  key="radius",
                                                  help="Muons start on the upper hemisphere surface.")
                    hemi_cz_m   = st.number_input("Centre z [m]", -10000.0, 10000.0, 0.0, 1.0,
                                                  key="sourcezm",
                                                  help="z of the equator / sphere centre.")
                    st.caption(f"Area ≈ {2*np.pi*hemi_radius**2/1e6:.4f} km²  |  "
                               f"z ∈ [{hemi_cz_m:.1f}, {hemi_cz_m+hemi_radius:.1f}] m")
                src_u1_m = src_u2_m = src_v1_m = src_v2_m = 0.0
                src_w_m  = hemi_cz_m
                radius   = hemi_radius
                plane_lx = plane_ly = 0.0
                disk_cx = disk_cy = disk_tilt = disk_tilt_az = 0.0
                disk_r  = hemi_radius
            else:
                with _sg_plane:
                    st.markdown("**Plane**")
                    source_plane = st.radio(
                        "Source plane", [1, 2, 3], index=0, horizontal=False,
                        key="source_plane",
                        format_func=lambda x: {
                            1: "XY  horizontal  (muons → −Z)",
                            2: "XZ  vertical    (muons → −Y)",
                            3: "YZ  vertical    (muons → −X)",
                        }[x],
                        help="Orientation of the generation surface in world coordinates.",
                        label_visibility="collapsed")

                _uname, _vname, _wname = {1: ("X","Y","Z"), 2: ("X","Z","Y"), 3: ("Y","Z","X")}[source_plane]
                _w_help = {
                    1: "z = 0 → surface. Negative → underground.",
                    2: "Fixed Y coordinate of the vertical XZ source plane.",
                    3: "Fixed X coordinate of the vertical YZ source plane.",
                }[source_plane]

                if source_mode == 1:
                    # ── Disk ──────────────────────────────────────────────────
                    _dc1, _dc2, _dc3 = st.columns(3)
                    disk_cx = _dc1.number_input(f"Center {_uname} [m]", -1e5, 1e5, 0.0, 1.0, key="disk_cx")
                    disk_cy = _dc2.number_input(f"Center {_vname} [m]", -1e5, 1e5, 0.0, 1.0, key="disk_cy")
                    disk_r  = _dc3.number_input("Radius [m]", 0.01, 1e5, 200.0, 1.0, key="disk_r")
                    _dw, _dt1, _dt2 = st.columns(3)
                    src_w_m      = _dw.number_input(f"{_wname} fixed [m]", -1e5, 1e5, 0.0, 1.0,
                                                    key="src_w", help=_w_help)
                    disk_tilt    = _dt1.number_input("Tilt [°]", 0.0, 89.9, 0.0, 1.0,
                                                     key="disk_tilt",
                                                     help="0 = flat in chosen plane. 90 = edge-on.")
                    disk_tilt_az = _dt2.number_input(f"Tilt azimuth [°]", 0.0, 360.0, 0.0, 5.0,
                                                     key="disk_tilt_az",
                                                     help=f"Tilt direction toward {_uname} axis = 0°.")
                    st.caption(f"r = {disk_r:.1f} m  |  Area = {np.pi*disk_r**2/1e6:.4f} km²"
                               + (f"  |  Tilt {disk_tilt:.1f}° @ {disk_tilt_az:.0f}°"
                                  if disk_tilt > 0.01 else ""))
                    src_u1_m = disk_cx - disk_r;  src_u2_m = disk_cx + disk_r
                    src_v1_m = disk_cy - disk_r;  src_v2_m = disk_cy + disk_r
                    radius   = disk_r
                else:
                    # ── Rectangle ─────────────────────────────────────────────
                    _r1, _r2 = st.columns(2)
                    src_u1_m = _r1.number_input(f"{_uname} min [m]", -1e5, 1e5, -200.0, 0.01,
                                                format="%.3f", key="src_u1")
                    src_u2_m = _r2.number_input(f"{_uname} max [m]", -1e5, 1e5,  200.0, 0.01,
                                                format="%.3f", key="src_u2")
                    _r3, _r4 = st.columns(2)
                    src_v1_m = _r3.number_input(f"{_vname} min [m]", -1e5, 1e5, -200.0, 0.01,
                                                format="%.3f", key="src_v1")
                    src_v2_m = _r4.number_input(f"{_vname} max [m]", -1e5, 1e5,  200.0, 0.01,
                                                format="%.3f", key="src_v2")
                    _rw, _rt1, _rt2 = st.columns(3)
                    src_w_m      = _rw.number_input(f"{_wname} fixed [m]", -1e5, 1e5, 0.0, 0.01,
                                                    format="%.3f", key="src_w", help=_w_help)
                    disk_tilt    = _rt1.number_input("Tilt [°]", 0.0, 89.9, 0.0, 1.0,
                                                     key="disk_tilt",
                                                     help="0 = flat in chosen plane. 90 = edge-on.")
                    disk_tilt_az = _rt2.number_input(f"Tilt azimuth [°]", 0.0, 360.0, 0.0, 5.0,
                                                     key="disk_tilt_az",
                                                     help=f"Tilt direction toward {_uname} axis = 0°.")
                    _hu = abs(src_u2_m - src_u1_m);  _hv = abs(src_v2_m - src_v1_m)
                    _rect_area = _hu * _hv
                    _rect_area_str = (f"{_rect_area/1e6:.4f} km²" if _rect_area >= 1e4
                                      else f"{_rect_area:.4f} m²")
                    st.caption(f"{_hu:.3f} × {_hv:.3f} m  |  Area = {_rect_area_str}"
                               + (f"  |  Tilt {disk_tilt:.1f}° @ {disk_tilt_az:.0f}°"
                                  if disk_tilt > 0.01 else ""))
                    disk_cx = (src_u1_m + src_u2_m) / 2.0
                    disk_cy = (src_v1_m + src_v2_m) / 2.0
                    radius  = min(abs(src_u2_m - src_u1_m), abs(src_v2_m - src_v1_m)) / 2.0
                    disk_r  = radius

                plane_lx   = abs(src_u2_m - src_u1_m) / 2.0
                plane_ly   = abs(src_v2_m - src_v1_m) / 2.0
                source_z_m = src_w_m
                hemi_radius = radius
                hemi_cz_m   = src_w_m

            st.divider()

            # ── Sampling ──────────────────────────────────────────────────────
            st.markdown("#### Sampling")
            _samp_c1, _samp_c2 = st.columns(2)
            nmuons_gen = int(_samp_c1.number_input(
                "Muons to generate", min_value=100, max_value=None, value=100_000, step=1_000,
                key="nmuonsgen",
                help=("**Detector filter ON:** this many *accepted* muons (hits); "
                      "total tried = N / acceptance.\n\n"
                      "**Detector filter OFF:** this many muons generated and written directly."),
            ))
            n_threads = int(_samp_c2.slider(
                "⚡ OpenMP threads", 1, min(64, os.cpu_count() or 8),
                min(4, os.cpu_count() or 4), 1, key="nthreads",
                help="Set to number of physical cores. Use `nproc` to check."))
            if not (_BIN_DIR / "ucmuon_gen_omp").exists():
                _samp_c2.warning("`ucmuon_gen_omp` not found — run `make local`")
            else:
                _samp_c2.success("✅  Generator ready", icon="✅")

    # ══════════════════════════════════════════════════════════════════════════
    # ROW 2 — conditional on workflow
    # ══════════════════════════════════════════════════════════════════════════
    st.divider()

    if use_dasrem:
        # ── DAS-REM: source plane height + footprint + sampling ───────────────
        # The reverse-mapping engine only needs the height of the horizontal
        # source plane — x/y follow from back-projecting the detector hit point.
        from ucmuon_dasrem_driver import _det_depth_cm as _dr_depth

        st.markdown("#### Source surface  *(reverse-mapping target)*")
        st.caption(
            "Guaranteed-hit mode samples the hit point on the detector top face and "
            "back-projects along the muon direction to a horizontal plane — the x/y "
            "footprint follows from the geometry, so only the plane **height z** is needed."
        )
        _drs1, _drs2 = st.columns([2, 3])
        src_w_m = _drs1.number_input(
            "Source plane z [m]", -1e5, 1e5, 0.0, 1.0, key="src_w",
            help="Height of the horizontal generation plane. z = 0 → ground surface. "
                 "Must be above the detector top face.")

        _det0 = detectors[0]
        if _det0["shape"] == 1:
            _det_cx_m   = (_det0["ax"] + _det0["bx"]) / 200.0
            _det_cy_m   = (_det0["ay"] + _det0["by"]) / 200.0
            _det_half_m = _det0["r"] / 100.0
        else:
            _det_cx_m   = (_det0["xmin"] + _det0["xmax"]) / 200.0
            _det_cy_m   = (_det0["ymin"] + _det0["ymax"]) / 200.0
            _det_half_m = max(_det0["xmax"] - _det0["xmin"],
                              _det0["ymax"] - _det0["ymin"]) / 200.0
        _dr_depth_val = _dr_depth(_det0, src_w_m * 100.0)          # cm
        _th_eff   = 0.0 if angular_mode == 1 else float(theta_max)
        _foot_r_m = (_det_half_m
                     + max(_dr_depth_val, 0.0) / 100.0 * np.tan(np.radians(_th_eff)))
        with _drs2:
            if _dr_depth_val <= 0:
                st.warning(
                    "⚠️  Source plane must be **above** the detector top face.  \n"
                    "Raise *Source plane z* or lower the detector z-coordinates."
                )
            else:
                st.success(
                    f"✅  Separation: **{_dr_depth_val/100:.1f} m**  "
                    f"(source z = {src_w_m:.1f} m → detector top)  —  DAS-REM ready."
                )
                st.caption(
                    f"Source footprint: r ≈ **{_foot_r_m:.1f} m** around the detector axis  "
                    f"(detector half-size + separation × tan θ_max)."
                )

        # Internal geometry (footprint disk) — keeps helpers & saved state consistent
        source_mode  = 1
        source_plane = 1
        radius  = disk_r = _foot_r_m
        disk_cx, disk_cy = _det_cx_m, _det_cy_m
        disk_tilt = disk_tilt_az = 0.0
        src_u1_m = disk_cx - radius;  src_u2_m = disk_cx + radius
        src_v1_m = disk_cy - radius;  src_v2_m = disk_cy + radius
        plane_lx = plane_ly = radius
        source_z_m  = src_w_m
        hemi_radius = radius
        hemi_cz_m   = src_w_m

        st.divider()
        st.markdown("#### Sampling")
        _samp_dr1, _samp_dr2 = st.columns([3, 2])
        nmuons_gen = int(_samp_dr1.number_input(
            "Muons to generate", min_value=100, max_value=None, value=100_000, step=1_000,
            key="nmuonsgen",
            help="All N muons are guaranteed to hit the detector — no wasted trials.",
        ))
        _samp_dr2.info("ℹ️  Guaranteed-hit mode uses a Python engine — no OpenMP threads needed.")
        n_threads = 1

    else:
        # ── Standard: detector filter (optional) ──────────────────────────────
        use_detector = st.checkbox("🔍  Enable surface detector filter  (ray → geometry intersection)",
                                   key="usedetector")
        output_sel = "output/muons_selected.dat"
        detectors  = []
        if use_detector:
            _df1, _df2 = st.columns([3, 2])
            output_sel = _df1.text_input("Selected muons file", "output/muons_selected.dat", key="outputsel")
            ndet = int(_df2.number_input("Number of detectors", 1, 10, 1, key="ndet"))
            _det_cols = st.columns(min(ndet, 3))
            for i in range(ndet):
                with _det_cols[i % len(_det_cols)]:
                    with st.expander(f"Detector {i+1}", expanded=(i == 0)):
                        shape  = st.selectbox("Shape", [1, 2], key=f"sh{i}",
                                     format_func=lambda x: "Cylinder" if x == 1 else "Box (AABB)")
                        margin = st.number_input(
                            "Safety margin [cm]", 0.0, 500.0, 0.0, key=f"mg{i}",
                            help="Acceptance halo added around the detector to catch "
                                 "muons scattered back in by MCS.")
                        st.caption("💡 Margin sizing → [🧰 Helpers & calculators](#helpers) "
                                   "→ 🎯 MCS margin")
                        d = {"shape": shape, "margin": margin}
                        if shape == 1:
                            _dc1, _dc2 = st.columns(2)
                            with _dc1:
                                st.markdown("**Bottom A**")
                                d["ax"] = st.number_input("Ax [cm]", value=0.0,     key=f"ax{i}")
                                d["ay"] = st.number_input("Ay [cm]", value=0.0,     key=f"ay{i}")
                                d["az"] = st.number_input("Az [cm]", value=-9000.0, key=f"az{i}")
                            with _dc2:
                                st.markdown("**Top B**")
                                d["bx"] = st.number_input("Bx [cm]", value=0.0, key=f"bx{i}")
                                d["by"] = st.number_input("By [cm]", value=0.0, key=f"by{i}")
                                d["bz"] = st.number_input("Bz [cm]", value=0.0, key=f"bz{i}")
                            d["r"] = st.number_input("Radius [cm]", 0.1, 1e4, 5.0, key=f"rr{i}")
                        else:
                            _dc1, _dc2 = st.columns(2)
                            with _dc1:
                                d["xmin"] = st.number_input("Xmin [cm]", value=-100.0,  key=f"xn{i}")
                                d["ymin"] = st.number_input("Ymin [cm]", value=-100.0,  key=f"yn{i}")
                                d["zmin"] = st.number_input("Zmin [cm]", value=-9000.0, key=f"zn{i}")
                            with _dc2:
                                d["xmax"] = st.number_input("Xmax [cm]", value=100.0, key=f"xx{i}")
                                d["ymax"] = st.number_input("Ymax [cm]", value=100.0, key=f"yx{i}")
                                d["zmax"] = st.number_input("Zmax [cm]", value=0.0,   key=f"zx{i}")
                        detectors.append(d)


        if not use_detector:
            st.caption(
                "💡  Switch to **Guaranteed-hit mode** (top toggle) for 100 % detector hits — no wasted trials."
            )

    # ══════════════════════════════════════════════════════════════════════════
    # ROW 3 — Output & export
    # ══════════════════════════════════════════════════════════════════════════
    st.divider()
    with st.expander("📁  Output & export", expanded=False):
        if use_dasrem:
            output_all = "output/muons_surface.dat"
            st.caption("ℹ️  Guaranteed-hit mode writes a **single file** — its name is set "
                       "in the *Detector → Output file* field above. Every muon in it is a "
                       "detector hit, so there is no hits/misses split.")
        else:
            output_all = st.text_input("Surface muon file", "output/muons_surface.dat",
                                       key="outputall")
        st.markdown("**Optional exports**")
        _ex1, _ex2, _ex3 = st.columns(3)
        with _ex1:
            if use_dasrem:
                save_all = False
                st.caption("**Save ALL muons** — not applicable: every generated muon "
                           "hits the detector by construction.")
            else:
                save_all = st.checkbox("Save ALL muons (hits + misses)", value=False, key="saveall",
                                       help="Detector filter ON: also write the file with every "
                                            "generated muon, not just detector hits.")
        with _ex2:
            save_phits   = st.checkbox("PHITS source  (s-type=17)", value=True, key="savephits")
            output_phits = "output/muons_for_phits.dat"
            if save_phits:
                output_phits = st.text_input("PHITS filename", "output/muons_for_phits.dat",
                                             key="outputphits")
        with _ex3:
            save_geant4   = st.checkbox("Geant4 source file", value=False, key="savegeant4")
            output_geant4 = "output/muons_geant4.txt"
            g4_fmt_label  = "UCMuon ASCII (with positions)"
            g4_use_all    = True
            if save_geant4:
                g4_fmt_label  = st.radio("Format",
                                         ["UCMuon ASCII (with positions)", "HEPEvt (G4HEPEvtInterface)"],
                                         key="g4_fmt_radio")
                output_geant4 = st.text_input(
                    "Geant4 filename",
                    "output/muons_geant4.txt" if "ASCII" in g4_fmt_label else "output/muons_geant4.hepevt")
                if use_dasrem:
                    g4_use_all = False   # single output file — all muons are hits
                    st.caption("ℹ️  Converts the guaranteed-hit output file "
                               "(every muon is a detector hit).")
                elif use_detector:
                    _g4src = st.radio(
                        "Source muons",
                        ["All generated muons  (recommended)", "Detector hits only"],
                        key="g4_source_sel",
                        help="Geant4 applies its own geometry — use all muons so Geant4 "
                             "can track every particle. 'Detector hits only' gives a smaller "
                             "file but skips muons that miss the UCMuon detector shape.")
                    g4_use_all = "All" in _g4src
                    if g4_use_all and not save_all:
                        st.caption("ℹ️  Enable **Save ALL muons** so the all-muons file "
                                   "is written; otherwise Geant4 falls back to detector hits.")

    # defaults if expander never opened
    if "output_all"  not in dir(): output_all = "output/muons_surface.dat"  # noqa: E701
    if "save_all"    not in dir(): save_all   = False                        # noqa: E701
    if "save_phits"  not in dir(): save_phits = True                         # noqa: E701
    if "output_phits" not in dir(): output_phits = "output/muons_for_phits.dat"  # noqa: E701
    if "save_geant4"  not in dir(): save_geant4 = False                      # noqa: E701
    if "output_geant4" not in dir(): output_geant4 = "output/muons_geant4.txt"  # noqa: E701
    if "g4_fmt_label" not in dir(): g4_fmt_label = "UCMuon ASCII (with positions)"  # noqa: E701

    # ══════════════════════════════════════════════════════════════════════════
    # ROW 4 — Run bar
    # ══════════════════════════════════════════════════════════════════════════
    st.divider()
    write_surface = not (use_detector and not save_all)
    if use_dasrem:
        _dr_out_preview = output_sel if use_detector else output_all
        _dr_extras = []
        if save_phits:
            _dr_extras.append(f"PHITS → `{output_phits}`")
        if save_geant4:
            _dr_extras.append(f"Geant4 → `{output_geant4}`")
        _dr_extra_str = "  |  " + "  |  ".join(_dr_extras) if _dr_extras else ""
        st.info(f"📄 Will write: `{_dr_out_preview}`{_dr_extra_str}  "
                f"(all {nmuons_gen:,} muons guaranteed to hit the detector)")
    elif use_detector:
        st.info(f"📄 Will write: `{output_all}` + `{output_sel}`" if save_all
                else f"📄 Will write: `{output_sel}` only  (detector filter ON, save-all OFF)")
    else:
        st.info(f"📄 Will write: `{output_all}`")

    _rb1, _rb2 = st.columns([4, 1])
    _run_label = ("▶  Run Generator  —  guaranteed-hit mode  (Python)"
                  if use_dasrem else "▶  Run UCMuon Surface Generator")
    run_gen  = _rb1.button(_run_label, type="primary", width='stretch',
                           disabled=_gg("gen_running"))
    stop_gen = _rb2.button("⛔  Stop", key="stop_gen", width='stretch',
                           disabled=not _gg("gen_running"))

    # ── Run logic ─────────────────────────────────────────────────────────────
    if run_gen and not _gg("gen_running"):
        st.session_state["save_geant4"]       = save_geant4
        st.session_state["g4_filename"]       = output_geant4
        st.session_state["g4_fmt"]            = "ascii" if "ASCII" in g4_fmt_label else "hepevt"
        st.session_state["gen_g4_use_all"]    = g4_use_all
        st.session_state["gen_geant4_done"]   = False
        st.session_state["gen_phits_done"]    = False
        st.session_state["gen_params_stored"] = False
        # Drop the previous run's flux: runs that print no "Integrated flux"
        # line (PARMA, DAS-REM, power-law) must not silently reuse a stale
        # value from an earlier spectrum in the rate estimate.
        st.session_state.pop("gen_integrated_flux", None)
        st.session_state["gen_use_dasrem"]    = use_dasrem
        st.session_state["gen_save_all"]      = save_all
        st.session_state["gen_output_all"]    = output_all
        st.session_state["gen_ntry"]          = None
        st.session_state["gen_theta_max"]     = theta_max

        # Shared session-state update (same for both modes)
        st.session_state.update({
            "gen_radius": radius, "gen_source_mode": source_mode,
            "gen_source_plane": source_plane,
            "gen_disk_cx": disk_cx, "gen_disk_cy": disk_cy,
            "gen_disk_tilt": disk_tilt, "gen_disk_tilt_az": disk_tilt_az,
            "gen_src_u1_m": src_u1_m, "gen_src_u2_m": src_u2_m,
            "gen_src_v1_m": src_v1_m, "gen_src_v2_m": src_v2_m,
            "gen_src_w_m":  src_w_m,
            "gen_source_z_m": src_w_m, "gen_plane_lx": plane_lx,
            "gen_plane_ly": plane_ly, "gen_nmuons_done": nmuons_gen,
            "gen_use_detector": use_detector, "gen_detectors": detectors,
            "gen_emin": e_min, "gen_emax": e_max,
            "gen_angular_mode": angular_mode,
        })
        if write_surface:
            st.session_state["surface_file"] = output_all
        if use_detector:
            st.session_state["selected_file"] = output_sel
            if not save_all:
                st.session_state["surface_file"] = output_sel

        if use_dasrem:
            # ── DAS-REM mode: pure Python generator ──────────────────────────
            _dr_out = output_sel if use_detector else output_all
            _dr_cfg = {
                "nmuons":        nmuons_gen,
                "spectrum_mode": spectrum_mode,
                "e_min":         e_min,
                "e_max":         e_max,
                "angular_mode":  angular_mode,
                "theta_max":     theta_max,
                "detectors":     detectors,
                "output_file":   _dr_out,
                "source_z_cm":   src_w_m * 100.0,  # m → cm
            }
            _gs("gen_lines",      [])
            _gs("gen_running",    True)
            _gs("gen_success",    None)
            _gs("gen_stop_req",   False)
            _gs("gen_nmuons",     nmuons_gen)
            _gs("gen_start_time", time.time())
            _gs("gen_end_time",   None)
            threading.Thread(
                target=_dasrem_worker,
                args=(_dr_cfg, _STATE, _LOCK),
                daemon=True,
            ).start()
            st.session_state["surface_file"] = _dr_out
            st.session_state["selected_file"] = _dr_out
        else:
            # ── Standard mode: Fortran generator ─────────────────────────────
            cfg = {
                "emin": e_min, "emax": e_max, "spectrum_mode": spectrum_mode,
                "parma_lat": parma_lat, "parma_lon": parma_lon, "parma_alt": parma_alt,
                "parma_year": int(parma_year), "parma_month": int(parma_month),
                "parma_day": int(parma_day), "parma_charge": int(parma_charge),
                "parma_sw": float(parma_sw),
                "source_mode": source_mode, "source_plane": source_plane,
                "disk_cx": disk_cx, "disk_cy": disk_cy, "disk_r": disk_r,
                "disk_tilt": disk_tilt, "disk_tilt_az": disk_tilt_az,
                "src_u1_m": src_u1_m, "src_u2_m": src_u2_m,
                "src_v1_m": src_v1_m, "src_v2_m": src_v2_m,
                "src_w_m":  src_w_m,
                "hemi_radius": hemi_radius, "hemi_cz_m": hemi_cz_m,
                "radius": radius, "plane_lx": plane_lx, "plane_ly": plane_ly,
                "source_z_m": src_w_m,
                "nmuons": nmuons_gen, "angular_mode": angular_mode, "theta_max": theta_max,
                "use_detector": use_detector, "detectors": detectors,
                "save_all": save_all, "save_phits": save_phits,
                "output_all": output_all, "output_sel": output_sel,
                "output_phits": output_phits,
            }
            _omp_env = {**os.environ, "OMP_NUM_THREADS": str(n_threads)}
            start_run([str(_BIN_DIR / "ucmuon_gen_omp")], build_ucmuon_input(cfg), "gen",
                      nmuons_gen, env=_omp_env)

        save_settings()
        st.rerun()

    if stop_gen:
        stop_run("gen"); st.rerun()

    live_panel("gen")

    # ── Post-run: manual PHITS export — fallback only, shown when the automatic
    #    "PHITS source (s-type=17)" export in 📁 Output & export is disabled ────
    if (_gg("gen_success") is True and not _gg("gen_running")
            and not st.session_state.get("savephits", True)):
        st.divider()
        if st.checkbox("Export surface muons as PHITS source", key="phits_surf_enable",
                       help="Manual converter — available because the automatic PHITS "
                            "export (📁 Output & export) is switched off."):
            _ps1, _ps2 = st.columns(2)
            _surf_cands = [f for f in [
                st.session_state.get("selected_file",""),
                st.session_state.get("surface_file",""),
                "output/muons_surface.dat", "output/muons_selected.dat",
            ] if f and Path(f).exists()]
            if not _surf_cands:
                st.warning("⚠️  No surface/selected file found.")
            else:
                _phits_src = _ps1.selectbox("Surface file", _surf_cands, key="phits_surf_src")
                _phits_out = _ps2.text_input("PHITS output", "output/muons_surface_phits.dat", key="phits_surf_out")
                if st.button("Convert → PHITS surface", key="btn_phits_surf", width='stretch'):
                    try:
                        _df_ps = load_file(_phits_src, mtime=Path(_phits_src).stat().st_mtime)
                        _n_ps  = write_phits_surface(_df_ps, _phits_out)
                        st.success(f"✅  Written: `{_phits_out}` ({_n_ps:,} muons)")
                        st.session_state["phits_surf_file"] = _phits_out
                    except Exception as _ex:
                        st.error(f"❌  PHITS surface export failed: {_ex}")
                if st.session_state.get("phits_surf_file") and Path(st.session_state["phits_surf_file"]).exists():
                    with open(st.session_state["phits_surf_file"], "rb") as _fh:
                        st.download_button(f"⬇️  Download {st.session_state['phits_surf_file']}",
                                           data=_fh, file_name=st.session_state["phits_surf_file"],
                                           mime="text/plain", width='stretch', key="dl_phits_surf")
                st.code(f"[Source]\n  s-type = 17\n  file   = {st.session_state.get('phits_surf_out','muons_surface_phits.dat')}\n  dump   = -10\n  1 2 3 4 5 6 7 8 9 10", language="text")

    # ── Post-run: confirm Fortran save_all output ─────────────────────────────
    if (_gg("gen_success") is True and not _gg("gen_running")
            and st.session_state.get("gen_save_all", False)
            and not st.session_state.get("gen_use_dasrem", False)):
        _all_path = Path(_abspath(st.session_state.get("gen_output_all",
                                                        "output/muons_surface.dat")))
        if _all_path.exists():
            st.success(f"✅  All-muons file: `{_all_path.name}`  "
                       f"({_all_path.stat().st_size // 1024} KB)")

    _auto_geant4_convert()
    _auto_phits_convert()
    _store_gen_params()

    # ── MCS acceptance estimator ──────────────────────────────────────────────
    if (_MCS_GUI_OK
            and _gg("gen_success") is True
            and not _gg("gen_running")
            and st.session_state.get("gen_use_detector", False)
            and st.session_state.get("gen_detectors")
            and st.session_state.get("selected_file", "")):
        render_mcs_panel(load_file)

    if st.session_state.get("save_geant4", False):
        with st.expander("📋  Geant4 C++ snippets"):
            st.markdown("**Option A — UCMuon ASCII reader**")
            st.code("""// PrimaryGeneratorAction.cc
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    std::string line;
    while (std::getline(fFile, line)) {
        if (line.empty() || line[0] == '#') continue;
        G4int pdg; G4double x,y,z,px,py,pz,Ekin;
        std::istringstream(line) >> pdg >> x >> y >> z >> px >> py >> pz >> Ekin;
        fGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(pdg));
        fGun->SetParticlePosition(G4ThreeVector(x*mm, y*mm, z*mm));
        fGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz).unit());
        fGun->SetParticleEnergy(Ekin * MeV);
        fGun->GeneratePrimaryVertex(event); break;
    }
}""", language="cpp")
            st.markdown("**Option B — HEPEvt**")
            st.code('#include "G4HEPEvtInterface.hh"\n'
                    'fGenerator = new G4HEPEvtInterface("muons_geant4.hepevt");',
                    language="cpp")

    # ══════════════════════════════════════════════════════════════════════════
    # ROW 5 — Helpers & calculators (planning tools — do not affect generation)
    # ══════════════════════════════════════════════════════════════════════════
    st.divider()
    st.markdown('<div id="helpers"></div>', unsafe_allow_html=True)
    st.markdown("#### 🧰  Helpers & calculators")
    st.caption("Planning tools — they inform parameter choices and do **not** "
               "change the generation.")

    # ── 📏 Source & detector geometry helpers ───────────────────────────────
    with st.expander("📏  Source size, MCS margin & detector solid angle", expanded=False):
        _hlp_src, _hlp_mcs, _hlp_sa = st.tabs([
            "📏  Source size & energy window",
            "🎯  MCS margin",
            "📐  Detector solid angle",
        ])
        with _hlp_src:
            if detectors and _SRCOPT_GUI_OK:
                render_combined_source_panel(
                    detectors=detectors, e_min=e_min, e_max=e_max,
                    theta_max=theta_max, source_mode=source_mode,
                    radius=radius, plane_lx=plane_lx, plane_ly=plane_ly,
                )
            else:
                st.info("Define at least one detector (enable the detector filter or "
                        "guaranteed-hit mode) to get source size, energy-window and "
                        "θ_max recommendations.")
        with _hlp_mcs:
            _render_mcs_margin_helper(detectors)
        with _hlp_sa:
            if detectors:
                _sa_sr, _ca_sr, _sa_msr, _frac = compute_detector_solid_angle(detectors)
                _sg1, _sg2, _sg3, _sg4 = st.columns(4)
                _sg1.metric("Solid angle",      f"{_sa_sr:.4e} sr")
                _sg2.metric("Solid angle",      f"{_sa_msr:.4f} msr")
                _sg3.metric("Fraction of 2π",   f"{_frac*100:.4f} %")
                _sg4.metric("cos²θ acceptance", f"{_ca_sr:.4e} sr")
                st.caption("MC estimate (600k rays from source centre). "
                           "cos²θ acceptance is used in the flux formula.")
            else:
                st.info("Define at least one detector to compute its solid angle.")

    # ── Equivalent real measurement time ─────────────────────────────────────
    with st.expander("⏱️  Equivalent real measurement time", expanded=False):
        st.caption(
            "How long a real detector would need to collect an equivalent number of muons "
            "passing through the source surface area. "
            "Formula: **t = N / (I_vert × Ω_eff × A_src)**"
        )
        if mono_beam:
            st.warning(
                "⚡ Mono-energetic beam: the rate estimate below is **not meaningful** for a "
                "delta-function spectrum. Disable mono-energetic mode for realistic exposure-time estimates."
            )

        _mtime_c1, _mtime_c2 = st.columns([3, 2])
        with _mtime_c1:
            # ── Flux model ────────────────────────────────────────────────────
            _ffe_model = st.selectbox(
                "Flux model for rate estimate",
                list(_FFE_MODEL_LABELS.keys()),
                index=list(_FFE_MODEL_LABELS.keys()).index("reyna_bugaev"),
                format_func=lambda k: _FFE_MODEL_LABELS[k],
                key="mtime_ffe_model",
                help="**Reyna-Bugaev** (recommended): matches PDG sea-level intensity "
                     "(I_V > 1 GeV ≈ 70 m⁻²sr⁻¹s⁻¹) to within ~20%.\n\n"
                     "**Guan/Frosin**: include an atmospheric energy-loss correction "
                     "(E_eff = E + 3.64 GeV) designed for underground acceptance — "
                     "they under-estimate the surface rate by 50–100× below 20 GeV."
            )

        with _mtime_c2:
            st.markdown("")  # vertical alignment spacer
            if _ffe_model in ("guan_2015", "frosin_2025") and e_min < 20.0:
                _gf_supp = int(((e_min + 0.106 + 3.64) / (e_min + 0.106)) ** 2.7)
                st.warning(
                    f"⚠️  {_FFE_MODEL_LABELS[_ffe_model].split('[')[0].strip()}: "
                    f"surface rate suppressed ~{_gf_supp}× at E_min = {e_min:.1f} GeV. "
                    f"Use Reyna-Bugaev for surface estimates.",
                    icon="⚠️"
                )

        # ── Source area — computed from the geometry the user already set ──────
        if source_mode == 1:
            _area_m2 = np.pi * radius ** 2
            _a_label = f"disk  r = {radius:.1f} m"
        elif source_mode == 3:
            _area_m2 = np.pi * radius ** 2
            _a_label = f"hemisphere  R = {radius:.1f} m (projected)"
        else:
            _area_m2 = 4.0 * plane_lx * plane_ly
            _a_label = f"{2*plane_lx:.1f} × {2*plane_ly:.1f} m rectangle"
        A_src_cm2 = _area_m2 * 1e4

        # ── N muons ───────────────────────────────────────────────────────────
        _ntry = st.session_state.get("gen_ntry", None)
        if use_dasrem:
            _N_t         = int(nmuons_gen)
            _N_lbl       = "guaranteed detector hits"
            _N_note      = None
            _dasrem_note = True
        elif use_detector and _ntry and _ntry > 0:
            _N_t         = int(_ntry)
            _N_lbl       = f"total tried  ({_ntry:,} from last run)"
            _N_note      = None
            _dasrem_note = False
        elif use_detector:
            _N_t         = int(nmuons_gen)
            _N_lbl       = "target hits  (run once to get exact tried count)"
            _N_note      = ("Run the generator once to get the exact tried count. "
                            "Until then this uses N_hits as a lower bound — actual time will be **longer**.")
            _dasrem_note = False
        else:
            _N_t         = int(nmuons_gen)
            _N_lbl       = "total generated"
            _N_note      = None
            _dasrem_note = False

        # ── Flux integral ─────────────────────────────────────────────────────
        _emin_flux = max(e_min, 0.5)
        _emax_flux = min(e_max, 1.5e4)
        _clamped   = _emin_flux > e_min
        _flux_err  = None
        if _emin_flux >= _emax_flux:
            I_vert    = 0.0
            _I_full   = 0.0
            _flux_err = f"E_min ({_emin_flux:.2f} GeV) ≥ E_max — set E_max > 0.5 GeV"
        else:
            try:
                _trapz  = getattr(np, "trapezoid", None) or np.trapz
                _Egrid  = np.geomspace(_emin_flux, _emax_flux, 500)
                I_vert  = float(_trapz(differential_flux(_Egrid, theta_deg=0.0, model=_ffe_model), _Egrid))
                _Efull  = np.geomspace(0.5, 1.5e4, 500)
                _I_full = float(_trapz(differential_flux(_Efull, theta_deg=0.0, model=_ffe_model), _Efull))
            except Exception as _e:
                I_vert    = 0.0
                _I_full   = 0.0
                _flux_err = str(_e)

        # ── Angular factor: R = I_V × 2π × (1 − cos⁴θ_max) / 4 ─────────────
        _ang_rad  = np.radians(float(theta_max))
        Omega_eff = 2.0 * np.pi * (1.0 - np.cos(_ang_rad) ** 4) / 4.0

        # ── Rate & time ───────────────────────────────────────────────────────
        _rate_s   = I_vert * Omega_eff * A_src_cm2
        _rate_min = _rate_s * 60.0
        _t_s      = (_N_t / _rate_s) if _rate_s > 0 else float("inf")
        _t_str    = _fmt_time(_t_s) if _rate_s > 0 else "—"

        _band_frac   = I_vert / _I_full if _I_full > 0 else 0.0
        _rate_full_m = _I_full * (np.pi / 2.0) * A_src_cm2 * 60.0

        # ── Metrics ───────────────────────────────────────────────────────────
        _tm1, _tm2, _tm3 = st.columns(3)
        _tm1.metric("Equivalent exposure time", _t_str,
                    help="t = N / (I_vert × Ω_eff × A_src)  —  time a real detector "
                         "covering the source area would take to collect N muons.")
        _tm2.metric("Surface crossing rate",
                    f"{_rate_min:,.0f} /min" if _rate_min >= 1.0
                    else f"{_rate_s:.2g} /s" if _rate_s > 0
                    else "—",
                    help=f"I_vert × Ω_eff × A_src  in [{_emin_flux:.1f}, {_emax_flux:.0f}] GeV")
        _tm3.metric("N muons", f"{_N_t:,}", help=_N_lbl)

        # ── Detail caption ────────────────────────────────────────────────────
        _e_range_str = (f"{_emin_flux:.1f}–{_emax_flux:.0f} GeV"
                        + ("  *(E_min clamped — model unreliable below 0.5 GeV)*"
                           if _clamped else ""))
        st.caption(
            f"I_vert = {I_vert:.3g} cm⁻²sr⁻¹s⁻¹  |  "
            f"E: {_e_range_str}  ({100*_band_frac:.1f}% of full spectrum)  |  "
            f"Ω_eff = {Omega_eff:.4f} sr  (θ ≤ {theta_max:.0f}°)  |  "
            f"A_src = {_area_m2:.2g} m²  ({_a_label})"
        )

        # ── Contextual warnings ───────────────────────────────────────────────
        if _area_m2 <= 0:
            st.warning(
                "⚠️  **Source area = 0 m²** — adjust the source geometry above "
                "(" + ("rectangle: set X_min ≠ X_max and Y_min ≠ Y_max"
                        if source_mode == 2 else "set radius > 0") + ")."
            )
        if _flux_err:
            st.warning(f"⚠️  Flux computation error: {_flux_err}")

        if _N_note:
            st.warning(_N_note)
        if _dasrem_note:
            st.info(
                "ℹ️  **Guaranteed-hit mode:** this is the time for N muons to cross the source surface "
                "regardless of the detector. In a real measurement you'd see fewer detector hits "
                "(acceptance < 100%) — actual run time is longer by 1/acceptance.",
                icon="ℹ️"
            )
        if _t_s < 10.0 and _rate_s > 0:
            st.info(
                f"⚡ Very fast rate ({_rate_min:,.0f} /min). "
                f"Consider reducing N or the source area if a quick test is sufficient.",
                icon="⚡"
            )
        if _band_frac < 0.5 and _rate_s > 0:
            st.info(
                f"E band [{_emin_flux:.1f}–{_emax_flux:.0f} GeV] = **{100*_band_frac:.1f}%** of full spectrum.  "
                f"Full-spectrum rate through this source: ≈ {_rate_full_m:,.0f} /min  "
                f"(cf. ~10 000 m⁻²min⁻¹ rule of thumb)."
            )

    # ── Energy threshold estimator ────────────────────────────────────────────
    with st.expander("⚡  Energy threshold estimator  (Groom 2001)", expanded=False):
        st.caption("Minimum muon kinetic energy to traverse a given rock thickness via CSDA range.")
        _et1, _et2 = st.columns(2)
        _et_depth = _et1.number_input("Rock thickness [m]", 0.1, 10000.0,
                                      float(st.session_state.get("music_depth_m", 90.0)),
                                      5.0, key="et_depth")
        _et_rho   = _et2.number_input("Rock density [g/cm³]", 0.1, 20.0,
                                      float(st.session_state.get("music_rho", 2.65)),
                                      0.05, key="et_rho")
        _et_opacity          = _et_depth * 100.0 * _et_rho
        _et_E_GeV, _et_E_MeV = _groom_threshold_energy(_et_opacity)
        _thr_str = f"{_et_E_GeV:.3f} GeV" if _et_E_GeV >= 1.0 else f"{_et_E_MeV:.0f} MeV"
        _em1, _em2, _em3 = st.columns(3)
        _em1.metric("Opacity  ρ·L", f"{_et_opacity:,.0f} g/cm²")
        _em2.metric("Min. penetrating energy", _thr_str)
        _em3.metric("→ Set Emin ≥", _thr_str)
        _et_T_GeV = _GROOM_T_MEV / 1000.0
        _fig_et   = go.Figure()
        _fig_et.add_trace(go.Scatter(
            x=_et_T_GeV, y=_GROOM_R_GCM2, mode="lines",
            line=dict(color="#38bdf8", width=3), name="Groom (2001) CSDA",
            hovertemplate="T = %{x:.4g} GeV<br>Range = %{y:.3g} g/cm²<extra></extra>",
        ))
        _fig_et.add_trace(go.Scatter(
            x=[_et_E_GeV], y=[_et_opacity], mode="markers",
            marker=dict(color="#ffd700", size=14, symbol="star", line=dict(color="#000", width=1.5)),
            name=f"⚡ Threshold {_thr_str}",
            hovertemplate=f"Threshold = {_thr_str}<extra></extra>",
        ))
        _fig_et.update_layout(
            height=240, margin=dict(l=60, r=20, t=10, b=40),
            paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
            font=dict(color="#e0e0e0", size=10), showlegend=True,
            legend=dict(bgcolor="rgba(20,22,35,0.95)", bordercolor="#888", borderwidth=1,
                        font=dict(size=10, color="#fff"), x=0.02, y=0.02,
                        xanchor="left", yanchor="bottom"),
            xaxis=dict(type="log", title="Kinetic energy T [GeV]", gridcolor="#2a2a3a", zeroline=False),
            yaxis=dict(type="log", title="CSDA Range [g/cm²]",     gridcolor="#2a2a3a", zeroline=False),
        )
        st.plotly_chart(_fig_et, config={"displayModeBar": False})
        _approx = ("" if abs(_et_rho - 2.65) < 0.05 else
                   f"  ⚠️ Table is for Standard Rock (ρ=2.65 g/cm³); ρ={_et_rho:.2f} gives an approximation.")
        st.info(f"💡 For **{_et_depth:.0f} m** at **{_et_rho:.2f} g/cm³** "
                f"(opacity {_et_opacity:,.0f} g/cm²): set **Emin ≥ {_thr_str}**.{_approx}", icon="⚡")
        st.caption("📚 Groom, Mokhov & Striganov, ADNDT 78 (2001). CSDA neglects straggling — treat as lower bound.")

    # ── Directional flux ──────────────────────────────────────────────────────
    with st.expander("🧭  Directional flux  [cm⁻² sr⁻¹ s⁻¹]", expanded=False):
        st.caption(
            "Sea-level muon flux **per solid angle** as a function of zenith angle θ. "
            "All parametrizations are **azimuth-symmetric** — Φ does not depend on φ. "
            "For geomagnetic / East–West effects use PARMA (spectrum ③)."
        )
        _dfc1, _dfc2, _dfc3, _dfc4 = st.columns([2, 2, 1, 1])
        _df_model = _dfc1.selectbox(
            "Flux model", list(_FFE_MODEL_LABELS.keys()),
            format_func=lambda k: _FFE_MODEL_LABELS[k], index=0, key="df_model",
            help=(
                "**Reyna–Bugaev** (for total rate): calibrated to PDG ±20%. "
                "Uses energy-independent cos_th_star^1.85 factor — correct for "
                "total flux (dominated by low-E muons), wrong for angular shape at E > 46 GeV.\n\n"
                "**Guan 2015 / Frosin 2025** (for angular analysis): explicitly models "
                "pion/kaon angular enhancement — above ~46 GeV at 30°, I(θ)/I(0°) > 1 "
                "(more oblique muons than vertical). Confirmed by CMS and IceCube. "
                "Best choice for muography where E_min > 20 GeV. "
                "Absolute flux unreliable below 20 GeV E_min.\n\n"
                "**Bugaev / Gaisser–Tang**: Gaisser formula only valid above 10 GeV. "
                "10× too low below that."
            )
        )
        _df_theta = _dfc2.slider("Zenith angle θ [°]", 0, 89, 0, 1, key="df_theta",
                                 help="0° = vertical. Rock path = L/cosθ — doubles at 60°.")
        _df_emin  = _dfc3.number_input("E_min [GeV]", 0.1, 10000.0,
                                       max(0.1, float(st.session_state.get("emin", 1.0))),
                                       0.1, format="%.2f", key="df_emin",
                                       help="Lower integration cut-off.")
        _df_alt   = _dfc4.number_input("Altitude [m]", 0, 5000, 0, 100, key="df_alt",
                                       help="Applies exp(h/8500 m) correction.")
        if _df_model in ("guan_2015", "frosin_2025") and float(_df_emin) < 20.0:
            st.warning(
                f"⚠️  **{_FFE_MODEL_LABELS[_df_model].split('[')[0].strip()}** with "
                f"E_min = {_df_emin:.1f} GeV: absolute flux I is **~×100 underestimated** "
                f"(E_eff correction dominates below 20 GeV). "
                f"The **angular ratio I(θ)/I(0°) is unaffected** and remains reliable. "
                f"For absolute flux → use Reyna–Bugaev.", icon="⚠️")
        if _df_model in ("bugaev", "gaisser_tang") and float(_df_emin) < 10.0:
            st.warning(
                f"⚠️  Gaisser formula valid above 10 GeV only. "
                f"At E_min = {_df_emin:.1f} GeV the absolute flux is ~×10 too low. "
                f"→ Use Reyna–Bugaev for E_min < 10 GeV.", icon="⚠️")
        try:
            _df_I_theta, _df_I_theta_T = angular_profile(
                np.array([0.0, float(_df_theta)]), E_min_GeV=float(_df_emin),
                model=_df_model, altitude_m=float(_df_alt))
            _df_I_vert  = _df_I_theta[0]
            _df_I_at_th = _df_I_theta[1]
            _df_ratio   = _df_I_theta_T[1]
            _dfm1, _dfm2, _dfm3, _dfm4, _dfm5 = st.columns(5)
            def _fmt_flux(v):
                if v <= 0: return "—"
                exp = int(np.floor(np.log10(v)))
                man = v / 10**exp
                return f"{man:.2f}×10{exp:+d}".replace("+","⁺").replace("-","⁻").replace(
                    "0","⁰").replace("1","¹").replace("2","²").replace("3","³").replace(
                    "4","⁴").replace("5","⁵").replace("6","⁶").replace("7","⁷").replace(
                    "8","⁸").replace("9","⁹")
            _dfm1.metric("I(0°)",          f"{_df_I_vert:.3g}",  help="[cm⁻²sr⁻¹s⁻¹]")
            _dfm2.metric(f"I({_df_theta}°)", f"{_df_I_at_th:.3g}", help="[cm⁻²sr⁻¹s⁻¹]")
            _dfm3.metric("I(θ)/I(0°)", f"{_df_ratio:.4f}" if _df_I_vert > 0 else "—",
                         help="Angular ratio relative to vertical.")
            _cos2 = np.cos(np.radians(_df_theta))**2
            _dfm4.metric("cos²θ (naive)", f"{_cos2:.4f}", delta=f"Δ = {_df_ratio-_cos2:+.4f}",
                         help="Naive approximation — compare to model curve.")
            _dfm5.metric("Model E_min", f"{_df_emin:.0f} GeV", help="Lower cut-off for integration.")
            st.divider()
            _BG = "rgb(15,17,23)"
            _dfp_l, _dfp_r = st.columns(2)
            with _dfp_l:
                st.markdown("**dΦ/dT spectrum** (log–log)")
                _T_plot = np.logspace(np.log10(max(float(_df_emin), 0.5)), 4.0, 300)
                _phi_vert = differential_flux(_T_plot, theta_deg=0.0, model=_df_model, altitude_m=float(_df_alt))
                _phi_th   = differential_flux(_T_plot, theta_deg=float(_df_theta), model=_df_model, altitude_m=float(_df_alt))
                _fig_spec = go.Figure()
                _fig_spec.add_trace(go.Scatter(x=_T_plot, y=_phi_vert, name="θ = 0°",
                    mode="lines", line=dict(color="rgba(56,189,248,0.9)", width=2)))
                if _df_theta > 0:
                    _fig_spec.add_trace(go.Scatter(x=_T_plot, y=_phi_th, name=f"θ = {_df_theta}°",
                        mode="lines", line=dict(color="rgba(251,146,60,0.9)", width=2, dash="dash")))
                _fig_spec.update_layout(
                    height=300, margin=dict(l=60, r=10, t=10, b=50),
                    paper_bgcolor=_BG, plot_bgcolor="rgb(20,22,30)", font=dict(color="white", size=10),
                    xaxis=dict(title="T [GeV]", type="log", gridcolor="#2a2a3a", zeroline=False),
                    yaxis=dict(title="dΦ/dT  [cm⁻²s⁻¹sr⁻¹GeV⁻¹]", type="log", gridcolor="#2a2a3a", zeroline=False),
                    legend=dict(font=dict(color="white", size=10), bgcolor="rgba(0,0,0,0.4)"))
                st.plotly_chart(_fig_spec, config={"displayModeBar": False})
            with _dfp_r:
                st.markdown("**Angular profile** I(θ) / I(0°)")
                _th_arr  = np.arange(0, 90, 2, dtype=float)
                _I_prof, _T_prof = angular_profile(_th_arr, E_min_GeV=float(_df_emin),
                    model=_df_model, altitude_m=float(_df_alt))
                _cos2_arr = np.cos(np.radians(_th_arr))**2
                _fig_ang = go.Figure()
                _fig_ang.add_trace(go.Scatter(x=_th_arr, y=_T_prof,
                    name=_FFE_MODEL_LABELS[_df_model].split("←")[0].split("[")[0].strip(),
                    mode="lines", line=dict(color="rgba(56,189,248,0.9)", width=2.5)))
                _fig_ang.add_trace(go.Scatter(x=_th_arr, y=_cos2_arr, name="cos²θ",
                    mode="lines", line=dict(color="rgba(255,255,255,0.35)", width=1.5, dash="dot")))
                _fig_ang.add_trace(go.Scatter(x=[float(_df_theta)], y=[_df_ratio], mode="markers",
                    marker=dict(size=11, color="rgba(251,146,60,1)", symbol="diamond",
                                line=dict(color="#fff", width=1.5)),
                    name=f"θ = {_df_theta}°  ({_df_ratio:.3f})"))
                _y_max = max(1.2, float(np.nanmax(_T_prof)) * 1.1)
                _fig_ang.update_layout(
                    height=300, margin=dict(l=55, r=10, t=10, b=50),
                    paper_bgcolor=_BG, plot_bgcolor="rgb(20,22,30)", font=dict(color="white", size=10),
                    xaxis=dict(title="θ [°]", range=[0,89], gridcolor="#2a2a3a", zeroline=False),
                    yaxis=dict(title="I(θ) / I(0°)", range=[0,_y_max], gridcolor="#2a2a3a", zeroline=False),
                    legend=dict(font=dict(color="white", size=9), bgcolor="rgba(0,0,0,0.4)"))
                st.plotly_chart(_fig_ang, config={"displayModeBar": False})
            st.caption(
                f"Model: **{_FFE_MODEL_LABELS[_df_model].split('←')[0].strip()}**  |  "
                f"E_min = {_df_emin:.1f} GeV  |  altitude = {_df_alt} m a.s.l.  |  "
                "All models azimuth-symmetric. "
                "cosθ* corrects for Earth's curvature at large θ (Guan 2015, arXiv:1509.06176).")
        except Exception as _df_err:
            st.error(f"Flux computation error: {_df_err}")

    # ── Fast flux estimator ────────────────────────────────────────────────
    with st.expander("🪨  Fast Flux Estimator — flux through a rock slab (semi-analytical)", expanded=False):
        st.caption(
            "Semi-analytical flux I [cm⁻²sr⁻¹s⁻¹] after traversing a flat rock slab. "
            "No MC needed — useful for site planning and exposure-time estimates. "
            "All five models are azimuth-symmetric at sea level."
        )

        # ── Inline guidance ────────────────────────────────────────────────
        with st.expander("📖  How to read these numbers — model guide & limitations", expanded=False):
            st.markdown("""
**What I [cm⁻²sr⁻¹s⁻¹] means:**  
Muons crossing 1 cm² per second, per steradian of solid angle, from direction θ.  
This is **not** a count rate — multiply by your detector acceptance A [cm²·sr] to get R [s⁻¹].

**Detector acceptance A [cm²·sr]:**

| Geometry | Formula |
|---|---|
| Single upward panel, area S, all angles | A = S × π |
| Single panel, cone θ < θ_max | A = S × 2π(1 − cosθ_max) |
| Two-panel telescope, area S, separation d | A = S²/d² |
| UCMuon detector filter (Generator tab) | A = area × MC cos²θ acceptance |

GUI default A = 6 cm²·sr ≈ 100 cm² × 0.06 sr (narrow telescope-like).

---

**Which model to use:**

| Model | For absolute rate | For angular shape I(θ)/I(0°) | Min E_min |
|---|---|---|---|
| **Reyna–Bugaev** ← for rates | ✅ ±20% of PDG | ⚠️ energy-independent, wrong > 46 GeV | 1 GeV |
| Bugaev / Gaisser–Tang | ⚠️ ×10 too low < 10 GeV | ⚠️ energy-independent | **10 GeV** |
| **Guan 2015** ← for angular | ⚠️ wrong < 20 GeV | ✅ Best: pion/kaon + cosθ* | **20 GeV** |
| **Frosin 2025** ← for angular | ⚠️ wrong < 20 GeV | ✅ Best: pion/kaon + cosθ* | **20 GeV** |

**Why Guan/Frosin are better for muography angular analysis:**  
Above ~46 GeV at 30°, oblique muons are *more* abundant than vertical muons because pion/kaon decays are enhanced at large angles (pions re-interact less). Guan explicitly models this pion/kaon angular enhancement. Reyna-Bugaev applies a constant cos_th_star^1.85 factor regardless of energy — correct for total flux, wrong at muography energies. This is confirmed by CMS, IceCube, and AMANDA measurements.

**Why Guan/Frosin absolute flux is wrong at low E_min:**  
Their correction term E_eff = E·(1 + a/(E·cosθ*^b)) pushes the effective energy to 4.6× the real energy at 1 GeV, strongly suppressing the flux. Their **angular ratio** I(θ)/I(0°) is correct because this bias cancels.  
→ Use Guan/Frosin only if E_min shown in the metrics is > 20 GeV (i.e. significant overburden).

---

**Flat-slab limitations — when NOT to trust this tool:**

- **θ > 50° through real terrain**: the path is NOT L/cosθ for a volcano. Use the **Terrain** tab with a DEM.  
- **X > 100 000 g/cm²**: models diverge by factors of 2–4 at extreme opacity. Full MC (MUSIC/PROPOSAL) is needed.  
- **Azimuth φ**: all models ignore the ~2% East–West geomagnetic asymmetry. Use PARMA (generator spectrum ③) for φ-dependence.  
- **CSDA E_min is a lower bound**: stochastic losses let some muons below E_min,CSDA survive. MUSIC thresholds are lower.

**Transmission T definition:**  
T = I(rock, θ) / I(open sky, same θ) — the fraction of muons that survive the rock *relative to the open-sky flux from the same direction*. Not relative to the vertical flux.
""")

        # ── Controls ───────────────────────────────────────────────────────
        _ff1, _ff2, _ff3 = st.columns(3)

        ffe_thickness_m = _ff1.number_input(
            "Rock L [m]", 0.0, 5000.0, 100.0, 10.0, key="ffe_thickness",
            help="Vertical rock thickness. Slant path = L/cosθ.")
        ffe_rho         = _ff1.number_input(
            "ρ [g/cm³]", 1.0, 5.0, float(RHO_STANDARD_ROCK), 0.05, key="ffe_rho",
            help="Rock density. Standard Rock=2.65, limestone≈2.5, volcanic tuff≈1.7.")

        ffe_theta_deg   = _ff2.slider(
            "Zenith θ [°]", 0, 89,
            int(st.session_state.get("ffe_theta", 0)), 1,
            key="ffe_theta",
            help="Muon arrival zenith angle. 0°=vertical. Rock path = L/cosθ — doubles at 60°.")
        ffe_altitude_m  = _ff2.number_input(
            "Altitude [m a.s.l.]", 0, 5000, 0, 100, key="ffe_altitude",
            help="Surface altitude. Correction ≈ exp(h/8500 m). Valid below ~4 km.")

        ffe_model       = _ff3.selectbox(
            "Flux model",
            list(_FFE_MODEL_LABELS.keys()),
            format_func=lambda k: _FFE_MODEL_LABELS[k],
            key="ffe_model",
            help=(
                "**Reyna–Bugaev** (for total rate): calibrated to PDG ±20%. "
                "Uses energy-independent cos_th_star^1.85 factor — correct for "
                "total flux (dominated by low-E muons), wrong for angular shape at E > 46 GeV.\n\n"
                "**Guan 2015 / Frosin 2025** (for angular analysis): explicitly models "
                "pion/kaon angular enhancement — above ~46 GeV at 30°, I(θ)/I(0°) > 1 "
                "(more oblique muons than vertical). Confirmed by CMS and IceCube. "
                "Best choice for muography where E_min > 20 GeV. "
                "Absolute flux unreliable below 20 GeV E_min.\n\n"
                "**Bugaev / Gaisser–Tang**: Gaisser formula only valid above 10 GeV. "
                "10× too low below that."
            )
        )
        ffe_acceptance  = _ff3.number_input(
            "Acceptance A [cm²·sr]", 0.01, 1000.0, 6.0, 0.5, key="ffe_acceptance",
            help=(
                "Geometric acceptance = detector area × effective solid angle.\n\n"
                "Single 10×10 cm² panel (full sky): 314 cm²·sr\n"
                "Single panel, θ<30° cone: 84 cm²·sr\n"
                "Telescope 10×10 cm², d=50 cm: 4 cm²·sr\n"
                "MURAVES-style telescope: ~6 cm²·sr (this default)"
            ))
        ffe_n_threshold = _ff3.number_input(
            "Target N for t_exp", 1, 100000, 100, 10, key="ffe_n_thresh",
            help="Exposure time to collect this many muons: t = N / (I × A).")

        # ── Computed quantities ────────────────────────────────────────────
        # Path through flat slab is L/cosθ — correct for non-vertical muons
        _cos_ffe = float(np.cos(np.radians(ffe_theta_deg)))
        ffe_path_m  = ffe_thickness_m / _cos_ffe if _cos_ffe > 0.01 else ffe_thickness_m
        ffe_opacity = ffe_rho * ffe_path_m * 100.0   # [g/cm²]

        # ── Contextual warnings based on current inputs ───────────────────
        _E_min_warn = emin_from_opacity(ffe_opacity)  # None if too deep
        _guan_models = ("guan_2015", "frosin_2025")
        _gaisser_models = ("bugaev", "gaisser_tang")

        if ffe_model in _guan_models and (
                _E_min_warn is None or _E_min_warn < 20.0):
            st.warning(
                f"⚠️  **{_FFE_MODEL_LABELS[ffe_model].split('[')[0].strip()}**: "
                f"E_min = {f'{_E_min_warn:.0f} GeV' if _E_min_warn else '< 20 GeV'} "
                f"— absolute flux is unreliable below 20 GeV with this model. "
                f"Absolute I shown will be **~×100 underestimated**. "
                f"Angular ratio I(θ)/I(0°) remains correct. "
                f"→ Switch to **Reyna–Bugaev** for rate estimates.",
                icon="⚠️")

        if ffe_model in _gaisser_models and (
                _E_min_warn is None or _E_min_warn < 10.0):
            st.warning(
                f"⚠️  **{_FFE_MODEL_LABELS[ffe_model].split('[')[0].strip()}**: "
                f"Gaisser formula is valid only above 10 GeV. "
                f"Current E_min = {f'{_E_min_warn:.0f} GeV' if _E_min_warn else '< 10 GeV'} "
                f"— flux will be ~×10 underestimated. "
                f"→ Switch to **Reyna–Bugaev**.",
                icon="⚠️")

        if ffe_theta_deg > 50 and ffe_thickness_m > 0:
            st.info(
                f"ℹ️  θ = {ffe_theta_deg}° > 50°: for real terrain (volcano, glacier) "
                f"the flat-slab path L/cosθ = {ffe_path_m:.0f} m overestimates the "
                f"true slant path. Use the **Terrain** tab with a DEM for accurate results.",
                icon="🏔️")

        if ffe_opacity > 100_000:
            st.warning(
                f"⚠️  X = {ffe_opacity:,.0f} g/cm² is very high. "
                f"Models diverge by factors of 2–4 at this opacity. "
                f"The FFE transmission T is an order-of-magnitude estimate only. "
                f"Run MUSIC or PROPOSAL for reliable results.",
                icon="⚠️")

        try:
            I_flux, E_min = integrated_flux(ffe_opacity, ffe_theta_deg,
                                            model=ffe_model,
                                            altitude_m=float(ffe_altitude_m))
            I_open, _     = integrated_flux(0.0, 0.0,
                                            model=ffe_model,
                                            altitude_m=float(ffe_altitude_m))
            I_open_th, _  = integrated_flux(0.0, ffe_theta_deg,
                                            model=ffe_model,
                                            altitude_m=float(ffe_altitude_m))
            transmission  = I_flux / I_open_th if I_open_th > 0 else 0.0

            # Metrics row
            _ff_m1, _ff_m2, _ff_m3, _ff_m4, _ff_m5 = st.columns(5)
            _ff_m1.metric("E_min [GeV]",
                          f"{E_min:.2f}" if E_min else "∞",
                          help="Min. kinetic energy to traverse X = ρ·L/cosθ (Groom 2001 CSDA).")
            _ff_m2.metric("I(θ) [cm⁻²sr⁻¹s⁻¹]",
                          f"{I_flux:.3e}" if I_flux > 0 else "0",
                          help="Integrated flux per unit solid angle at the detector.")
            _ff_m3.metric("I(0°) [cm⁻²sr⁻¹s⁻¹]",
                          f"{I_open:.3e}" if I_open > 0 else "—",
                          help="Vertical open-sky flux (no rock) for reference.")
            _ff_m4.metric("Rock transmission T",
                          f"{transmission:.4f}" if transmission > 0 else "0",
                          delta=f"path {ffe_path_m:.0f} m  X={ffe_opacity:.0f} g/cm²",
                          help="I(θ,rock) / I(θ,open). Accounts for L/cosθ path length.")
            _rate_s  = I_flux * ffe_acceptance
            _t_exp_s = float(ffe_n_threshold) / _rate_s if _rate_s > 0 else float('inf')
            _t_str   = (f"{_t_exp_s:.1f} s"             if _t_exp_s < 3600         else
                        f"{_t_exp_s/3600:.2f} h"         if _t_exp_s < 86400        else
                        f"{_t_exp_s/86400:.1f} days"     if _t_exp_s < 86400*365    else
                        f"{_t_exp_s/86400/365.25:.2f} yr")
            _ff_m5.metric(f"t_exp ({ffe_n_threshold:,} μ)", _t_str,
                          delta=f"rate {_rate_s*86400:.3g} /day",
                          help="Time to accumulate target muon count: N / (I·A).")

            st.caption(
                f"X = ρ·L/cosθ = {ffe_opacity:.0f} g/cm²  |  "
                f"I₀(vertical) = {I_open:.3e} cm⁻²sr⁻¹s⁻¹  |  "
                f"Model: {_FFE_MODEL_LABELS[ffe_model].split('←')[0].strip()}"
            )

            # ── Plots: flux vs depth  +  angular profile side by side ────
            _BG = "rgb(15,17,23)"
            _ffe_pl, _ffe_pr = st.columns(2)

            with _ffe_pl:
                st.markdown("**I(L) and T vs depth** at θ = {}°".format(ffe_theta_deg))
                _depths = np.linspace(0.0, max(ffe_thickness_m * 1.5, 50.0), 120)
                try:
                    # Use L/cosθ path for each depth point
                    _paths  = _depths / _cos_ffe if _cos_ffe > 0.01 else _depths
                    _I_arr, _T_arr, _ = flux_vs_depth(
                        _paths, ffe_rho, ffe_theta_deg,
                        model=ffe_model, altitude_m=float(ffe_altitude_m))
                    _fig_ffe = go.Figure()
                    _fig_ffe.add_trace(go.Scatter(
                        x=_depths, y=_I_arr,
                        name="I(L, θ)  [cm⁻²sr⁻¹s⁻¹]",
                        mode="lines",
                        line=dict(color="rgba(0,230,255,0.9)", width=2),
                        yaxis="y1",
                    ))
                    _fig_ffe.add_trace(go.Scatter(
                        x=_depths, y=_T_arr,
                        name="Transmission T",
                        mode="lines",
                        line=dict(color="rgba(255,165,0,0.85)", width=2, dash="dash"),
                        yaxis="y2",
                    ))
                    if ffe_thickness_m > 0 and I_flux > 0:
                        _fig_ffe.add_trace(go.Scatter(
                            x=[ffe_thickness_m], y=[I_flux],
                            mode="markers",
                            marker=dict(size=10, color="rgba(255,80,80,1)",
                                        symbol="diamond",
                                        line=dict(color="#fff", width=1)),
                            yaxis="y1",
                            name=f"L = {ffe_thickness_m:.0f} m",
                        ))
                    _fig_ffe.update_layout(
                        height=300, margin=dict(l=60, r=70, t=10, b=50),
                        paper_bgcolor=_BG, plot_bgcolor=_BG,
                        xaxis=dict(title="Vertical depth L [m]", color="white",
                                   gridcolor="rgba(255,255,255,0.08)"),
                        yaxis=dict(title="Flux [cm⁻²sr⁻¹s⁻¹]", type="log",
                                   color="rgba(0,230,255,0.9)",
                                   gridcolor="rgba(255,255,255,0.08)"),
                        yaxis2=dict(title="T = I(rock)/I(open)",
                                    overlaying="y", side="right",
                                    range=[0, 1],
                                    color="rgba(255,165,0,0.85)"),
                        legend=dict(font=dict(color="white", size=9),
                                    bgcolor="rgba(0,0,0,0.4)"),
                    )
                    st.plotly_chart(_fig_ffe,                                         config={"displayModeBar": False})
                except Exception as _pe:
                    st.warning(f"Depth plot error: {_pe}")

            with _ffe_pr:
                st.markdown("**Angular profile** I(θ) / I(0°) — surface, L = 0")
                try:
                    _th_p = np.arange(0, 90, 2, dtype=float)
                    _I_ap, _T_ap = angular_profile(
                        _th_p, E_min_GeV=float(E_min) if E_min else 1.0,
                        model=ffe_model, altitude_m=float(ffe_altitude_m),
                    )
                    _cos2_p = np.cos(np.radians(_th_p))**2
                    _fig_ap = go.Figure()
                    _fig_ap.add_trace(go.Scatter(
                        x=_th_p, y=_T_ap,
                        name=_FFE_MODEL_LABELS[ffe_model].split("←")[0].split("[")[0].strip(),
                        mode="lines",
                        line=dict(color="rgba(56,189,248,0.9)", width=2.5),
                    ))
                    _fig_ap.add_trace(go.Scatter(
                        x=_th_p, y=_cos2_p,
                        name="cos²θ",
                        mode="lines",
                        line=dict(color="rgba(255,255,255,0.3)", width=1.5, dash="dot"),
                    ))
                    # Mark the currently chosen zenith angle
                    _I_at, _T_at = angular_profile(
                        np.array([float(ffe_theta_deg)]),
                        E_min_GeV=float(E_min) if E_min else 1.0,
                        model=ffe_model, altitude_m=float(ffe_altitude_m),
                    )
                    _fig_ap.add_trace(go.Scatter(
                        x=[float(ffe_theta_deg)], y=[float(_T_at[0])],
                        mode="markers",
                        marker=dict(size=11, color="rgba(251,146,60,1)",
                                    symbol="diamond",
                                    line=dict(color="#fff", width=1.5)),
                        name=f"θ = {ffe_theta_deg}°  ({_T_at[0]:.3f})",
                    ))
                    _fig_ap.update_layout(
                        height=300, margin=dict(l=55, r=10, t=10, b=50),
                        paper_bgcolor=_BG, plot_bgcolor=_BG,
                        xaxis=dict(title="θ [°]", range=[0, 89], color="white",
                                   gridcolor="rgba(255,255,255,0.08)"),
                        yaxis=dict(title="I(θ) / I(0°)", range=[0, 1.05],
                                   color="white",
                                   gridcolor="rgba(255,255,255,0.08)"),
                        legend=dict(font=dict(color="white", size=9),
                                    bgcolor="rgba(0,0,0,0.4)"),
                    )
                    st.plotly_chart(_fig_ap,                                         config={"displayModeBar": False})
                except Exception as _pe:
                    st.warning(f"Angular profile error: {_pe}")

        except Exception as _ffe_err:
            st.error(f"Fast flux estimator error: {_ffe_err}")
            I_flux = E_min = transmission = 0.0

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — TRANSPORT
# ══════════════════════════════════════════════════════════════════════════════
with tab_music:

    # ── Engine selector ───────────────────────────────────────────────────────
    # Internal option keys are kept stable (saved configs store them);
    # only the display labels and ordering changed when UCMuon-MC became Engine ①.
    _ENGINE_OPTIONS = [
        "UCMuon Stochastic (Python)",
        "MUSIC",
        "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS",
        "PROPOSAL",
        "Backward MC Flux Integrator",
        "PUMAS",
    ]
    _ENGINE_LABELS = {
        "UCMuon Stochastic (Python)":                                        "① ★ UCMuon-MC — native stochastic MC (flagship)",
        "MUSIC":                                                              "② MUSIC — Kudryavtsev stochastic (OMP)",
        "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS":          "③ UCMuon CSDA — Bethe-Bloch + MS (OMP)",
        "PROPOSAL":                                                           "④ PROPOSAL — full stochastic MC",
        "Backward MC Flux Integrator":                                        "⑤ Backward MC flux integrator",
        "PUMAS":                                                              "⑥ PUMAS — backward MC + forward (C)",
    }

    _eng1, _eng2 = st.columns([2, 3])
    with _eng1:
        transport_engine = st.selectbox(
            "Transport engine",
            _ENGINE_OPTIONS,
            format_func=lambda x: _ENGINE_LABELS.get(x, x),
            key="transport_engine",
            help=(
                "**① UCMuon-MC** (flagship): native stochastic MC — PDG per-process radiative "
                "sampling + δ-ray straggling + Highland MS + decay. Pure Python, no install needed.\n\n"
                "**② MUSIC**: Full stochastic MC — Kudryavtsev XS tables. OMP-parallel.\n\n"
                "**③ Bethe-Bloch**: Analytical a+bE loss + Highland MS. OMP-parallel. No external files.\n\n"
                "**④ PROPOSAL**: Full stochastic MC (IceCube/KM3NeT). Requires `pip install proposal`.\n\n"
                "**⑤ Backward MC**: Flux integrator — no muon file needed.\n\n"
                "**⑥ PUMAS**: True backward MC (Niess 2017). "
                "Forward or backward mode. No muon file needed in backward mode. Run `make pumas`.\n\n"
                "📌 For real terrain (volcano, glacier): use the **Terrain** tab."
            )
        )
    _ENGINE_DESC = {
        "MUSIC": (
            "**Full stochastic MC** (Fortran/OMP).  "
            "Kudryavtsev (2009) cross-section tables for ionisation, bremsstrahlung, "
            "pair production, and photonuclear interactions.  "
            "Industry standard for deep-underground muon flux."
        ),
        "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS": (
            "**Continuous slowing-down** (CSDA, Fortran/OMP).  "
            "PDG Bethe-Bloch for ionisation + Groom et al. (2001) parameterised radiative losses "
            "+ Highland formula for multiple Coulomb scattering.  "
            "No external files — fast and portable."
        ),
        "PROPOSAL": (
            "**Full stochastic MC** (external library).  "
            "PROPOSAL (Koehne et al. 2013) — used by IceCube and KM3NeT.  "
            "Models ionisation, bremsstrahlung, pair production, photonuclear, and decay "
            "with user-selectable cross-section parametrisations."
        ),
        "UCMuon Stochastic (Python)": (
            "**UCMuon-MC v2 — the native flagship engine** (pure Python/NumPy).  "
            "PDG 2024 per-process radiative losses (brems / pair / photonuclear) "
            "Poisson-sampled with process-specific spectra + explicit δ-ray "
            "straggling + Highland MS + muon decay; mean dE/dx anchored exactly "
            "to the evaluated table.  Multiprocess-parallel, no compilation, "
            "no external files, any platform — validated against MUSIC and "
            "PROPOSAL.  Default choice for slab and terrain transport."
        ),
        "Backward MC Flux Integrator": (
            "**Flux integrator** — no surface muon file needed.  "
            "Convolves the surface spectrum with pre-computed survival probabilities "
            "to return the underground flux directly.  "
            "Based on Kudryavtsev (2008) backward MC formalism."
        ),
        "PUMAS": (
            "**True backward MC** (PUMAS — Niess et al. 2017).  "
            "Forward mode: reads surface muon file, transports via PUMAS, writes 18-col output.  "
            "Backward mode: samples detector energy/angle, propagates backward through rock, "
            "returns weighted flux spectrum — no surface muon file needed.  "
            "Supports CSDA, mixed, and straggled energy-loss modes."
        ),
    }
    st.caption(_ENGINE_DESC.get(transport_engine, ""))

    with _eng2:
        if transport_engine == "MUSIC":
            _omp_ok_sel = (_BIN_DIR / "ucmuon_transport_music_omp").exists()
            _ser_ok_sel = (_BIN_DIR / "ucmuon_transport_music").exists()
            if _omp_ok_sel:
                st.success("✅ **MUSIC** ready — `ucmuon_transport_music_omp` found.")
            elif _ser_ok_sel:
                st.warning("⚠️ **MUSIC** serial binary found (no OMP). Run `make local` for parallel.")
            else:
                st.error("❌ **MUSIC** binary not found. Run `make local`.")
            st.caption("Kudryavtsev stochastic transport · Requires energy-loss + cross-section tables · Flat-slab geometry")
        elif transport_engine == "PROPOSAL":
            _prop_ok, _prop_ver = _check_proposal()
            _driver_ok = (_SCRIPT_DIR / "proposal_driver.py").exists()
            if _prop_ok and _driver_ok:
                st.success(f"✅ **PROPOSAL** v{_prop_ver} ready.")
            elif _driver_ok and not _prop_ok:
                st.error(f"❌ `proposal_driver.py` found but PROPOSAL not installed: {_prop_ver}")
                st.caption("`pip install proposal`")
            else:
                st.error("❌ `proposal_driver.py` not found in project directory.")
            st.caption("Full stochastic MC · ionization · brems · pair · photonuclear · LPM · Molière/Highland")
        elif transport_engine == "UCMuon Stochastic (Python)":
            _stochastic_ok, _stochastic_ver = stochastic_available()
            if _stochastic_ok:
                st.success(f"✅ **UCMuon-MC** v{_stochastic_ver} ready.")
            else:
                st.error(f"❌ UCMuon-MC not ready — {_stochastic_ver}")
            st.caption("PDG per-process radiative sampling + δ-ray straggling + Highland MS + decay · Pure Python · multiprocess")
        elif transport_engine == "Backward MC Flux Integrator":
            _stochastic_ok2, _ = stochastic_available()
            if _stochastic_ok2:
                st.success("✅ **Backward MC** ready — no input muon file needed.")
            else:
                st.warning("⚠️ **Backward MC** not ready — place `gui_stochastic_engine.py` in `gui/`.")
            st.caption("CSDA backward inversion · flux integrator · no generator step needed")
        elif transport_engine == "PUMAS":
            _pumas_bin_ok = (_BIN_DIR / "ucmuon_transport_pumas").exists()
            _pumas_drv_ok = (_SCRIPT_DIR / "ucmuon_pumas_driver.py").exists()
            _pumas_mdf_ok = (_PROJECT_DIR / "external" / "pumas-master" / "examples" / "data" / "materials.xml").exists()
            if _pumas_bin_ok and _pumas_drv_ok and _pumas_mdf_ok:
                st.success("✅ **PUMAS** ready — binary, driver, and MDF found.")
            elif not _pumas_bin_ok:
                st.error("❌ Binary not found — run `make pumas`.")
                st.code("make pumas", language="bash")
            elif not _pumas_mdf_ok:
                st.error("❌ `external/pumas-master/examples/data/materials.xml` not found.")
            else:
                st.error("❌ `ucmuon_pumas_driver.py` not found in `gui/`.")
            _pumas_dump_ok = (_BIN_DIR / "pumas_StandardRock.pumas").exists()
            if _pumas_dump_ok:
                st.caption("Physics dump cached · first run already complete · fast startup")
            else:
                st.caption("First run builds physics tables (~10 s) and caches them in `bin/`")
        else:  # BB+MS
            _bbms_ok_sel = (_BIN_DIR / "ucmuon_transport_bb_omp").exists()
            if _bbms_ok_sel:
                st.success("✅ **Bethe-Bloch + Highland MS** ready — `ucmuon_transport_bb_omp` found.")
            else:
                st.error("❌ `ucmuon_transport_bb_omp` not compiled. Run `make local`.")
            st.caption("Analytical a+bE energy loss (PDG/Groom) + Highland MS · OMP-parallel · no external files")

    with st.expander("⚛️  Engine comparison — physics & accuracy guide", expanded=False):
        st.markdown(r"""
| # | Engine | Physics model | Pros | Cons |
|---|---|---|---|---|
| ① | **UCMuon-MC** ★ | PDG per-process radiative MC + δ-ray straggling + Highland MS + decay | **Native flagship**; pure Python; no install; any platform; multiprocess-parallel; mean dE/dx exact by construction | No LPM (minor below 1 TeV) |
| ② | **MUSIC** | Full stochastic MC (Kudryavtsev 2009) | External reference; Landau fluctuations; OMP | Requires table files; init ~1 min on first run |
| ③ | **UCMuon CSDA** | Bethe-Bloch $a+bE$ + Highland MS | Fastest; OMP; no files | No Landau fluctuations → overestimates survival by ~5–20% |
| ④ | **PROPOSAL** | Full stochastic MC (Koehne/Alameddine 2013/2024) | Landau; LPM; 3D MS; independent check | Requires `pip install proposal`; first run ~60 s |
| ⑤ | **Backward MC** | Flux integrator (CSDA + stochastic) | No muon file needed; gives flux at depth | No individual muon tracking |
| ⑥ | **PUMAS** | True backward MC (Niess 2017) | No muon file in backward mode; 100% efficiency; CSDA/mixed/straggled | Requires `make pumas`; C binary |

**Survival fraction ordering** (same input, same geometry):

$$\text{MUSIC} \lesssim \text{PROPOSAL} \lesssim \text{UCMuon-MC} < \text{UCMuon CSDA}$$

| Overburden | Recommended engine | Notes |
|---|---|---|
| < 200 m.w.e. | **① UCMuon-MC** | All engines agree within ±5% |
| 200–1000 m.w.e. | **① UCMuon-MC** (cross-check ② MUSIC / ④ PROPOSAL) | Stochastic fluctuations matter; ③ CSDA +10% |
| > 1000 m.w.e. | ② MUSIC or ④ PROPOSAL (validate ①) | Hard radiative losses dominate; CSDA +20% |

UCMuon-MC agrees with MUSIC within 0.6 pp at the 500 m benchmark. MUSIC vs PROPOSAL differ by ~10% in bremsstrahlung parametrisation — both physically valid, known inter-code systematic.
""")
        st.caption("Kudryavtsev (2009) CPC 180, 339 · Koehne+ (2013) CPC 184, 2070 · Alameddine+ (2024) CPC 305, 109243 · Groom+ (2001) ADNDT 78, 183")

    st.divider()

    if transport_engine == "Backward MC Flux Integrator":
        if _STOCHASTIC_GUI_OK and render_backward_mc_tab is not None:
            render_backward_mc_tab(_SCRIPT_DIR)
        else:
            st.error("❌ Place `gui_stochastic_engine.py` and `ucmuon_backward_mc.py` in `gui/`.")
    else:  # transport engine block (shared material/depth/run widgets)
        proposal_medium_type = 1
        proposal_e_cut   = 500.0
        proposal_v_cut   = 0.001
        proposal_scatter = 2
        proposal_tables  = ""
        proposal_custom  = {}
        stochastic_v_cut         = 0.05
        stochastic_n_steps       = 0
        stochastic_ms_enable     = True
        stochastic_range_table   = 1
        stochastic_hard_spectrum = 2
        stochastic_n_workers     = 0
        stochastic_delta_rays    = True
        stochastic_mat_id    = 1
        stochastic_custom    = {}

        # ── Files  |  Material ────────────────────────────────────────────────────
        st.markdown("##### Input / Output")
        _tcol_l, _tcol_r = st.columns(2)

        with _tcol_l:
            # Prefer surface file over selected by default (selected requires detector filter)
            _infile_candidates = list(dict.fromkeys(f for f in [
                st.session_state.get("transport_infile", ""),   # persisted from any engine
                st.session_state.get("surface_file",  ""),
                st.session_state.get("selected_file", ""),
                "output/muons_surface.dat", "output/muons_selected.dat",
            ] if f))
            _infile_existing = [f for f in _infile_candidates if Path(f).exists()]

            _mi1, _mi2 = st.columns([3, 2])
            if _infile_existing:
                _saved_idx = 0
                _saved_infile = st.session_state.get("transport_infile", "")
                if _saved_infile in _infile_existing:
                    _saved_idx = _infile_existing.index(_saved_infile)
                _infile_choice = _mi1.selectbox("Input file", _infile_existing,
                                                index=_saved_idx,
                                                key="transport_infile_select",
                                                help="Files detected from the Generator tab output.")
            else:
                _infile_choice = None
                _mi1.warning("⚠️  No input files found — run the generator first.")
            _infile_custom = _mi2.text_input("Or custom filename", value="",
                                             key="transport_infile_custom",
                                             placeholder="e.g. muons_step2.dat",
                                             help="Overrides the dropdown when filled.")

            if _infile_custom:
                m_infile = _infile_custom
                if Path(_infile_custom).exists():
                    st.success(f"✅  Using: `{m_infile}`")
                else:
                    st.error(f"❌  Not found: `{_infile_custom}`")
            elif _infile_choice:
                m_infile = _infile_choice
            else:
                m_infile = "output/muons_surface.dat"

            # Persist for CSDA panel and other engines
            st.session_state["transport_infile"] = m_infile

            if Path(m_infile).exists():
                _mt = Path(m_infile).stat().st_mtime
                ncols_info, n_transport = probe_music_file(m_infile, False, mtime=_mt)
                if ncols_info == 14:
                    st.caption(f"📄 14-col · {n_transport:,} muons with hit_flag")
                elif ncols_info == 13:
                    st.caption(f"📄 13-col · {n_transport:,} muons")
                else:
                    st.warning("⚠️  Could not read input file.")
                    ncols_info, n_transport = 0, 0
            else:
                st.warning(f"⚠️  Not found: `{m_infile}`")
                ncols_info, n_transport = 0, 0

        with _tcol_r:
            m_outfile = st.text_input("Output underground file", "output/muons_underground.dat",
                                      key="transport_outfile")

            transport_all = st.checkbox(
                "Transport ALL muons  (ignore hit_flag)", value=False,
                key="transport_all_chk",
                help="Default: only hit_flag=1 muons (detector-intercepting). "
                     "Enable to transport all muons — useful without a detector filter.")

            if ncols_info == 14 and Path(m_infile).exists():
                _mt2 = Path(m_infile).stat().st_mtime
                ncols_info, n_transport = probe_music_file(m_infile, transport_all, mtime=_mt2)
                _label = "all" if transport_all else "hit_flag=1 only"
                st.info(f"**{n_transport:,}** muons will be transported  ({_label})")
    
        st.divider()
        st.markdown("##### Medium & Geometry")
        _mat_sect = st.container()
        with _mat_sect:
    
            mat_choice = st.selectbox("Material preset", list(MUSIC_MATERIALS.keys()), index=0,
                                       key="musicpreset",
                                      help="Sets default density and radiation length.")
            mat = MUSIC_MATERIALS[mat_choice]
            if mat["desc"]:
                st.caption(f"ℹ️ {mat['desc']}")
    
            _preset_prev = st.session_state.get("music_preset_prev")
            if _preset_prev != mat_choice:
                # Apply preset defaults only on a real preset CHANGE.
                # On the first render of a session _preset_prev is unset —
                # resetting there would clobber the ρ/X₀ values restored
                # from the autosave file.
                if _preset_prev is not None and mat_choice != "Custom":
                    st.session_state["music_rho"] = float(mat["rho"])
                    st.session_state["music_rad"] = float(mat["rad"])
                    st.session_state["music_rho_sigma"] = 0.0
                st.session_state["music_preset_prev"] = mat_choice
    
            default_rho = mat["rho"] if mat["rho"] is not None else 2.65
            default_rad = mat["rad"] if mat["rad"] is not None else 26.48
    
            # ── Density input ──────────────────────────────────────────────────────
            density_mode = st.radio("Density input", ["Fixed", "Gaussian prior", "KDE from samples"],
                                    horizontal=True, key="density_mode")
            if density_mode == "Fixed":
                m_rho = st.number_input("Density ρ [g/cm³]", 0.001, 25.0,
                                        value=st.session_state.get("music_rho", float(default_rho)),
                                        step=0.05, key="music_rho")
                m_rho_sigma = 0.0
                st.session_state["music_rho_sigma"] = 0.0
    
            elif density_mode == "Gaussian prior":
                _dc1, _dc2 = st.columns([3, 2])
                m_rho = _dc1.number_input("Mean ρ [g/cm³]", 0.001, 25.0,
                                          value=st.session_state.get("music_rho", float(default_rho)),
                                          step=0.05, key="music_rho")
                m_rho_sigma = _dc2.number_input("σ [g/cm³]", 0.001, 5.0,
                                                value=max(0.001, float(st.session_state.get("music_rho_sigma", 0.10))),
                                                step=0.01, key="music_rho_sigma")
                _px  = np.linspace(max(0.001, m_rho-4*m_rho_sigma), m_rho+4*m_rho_sigma, 300)
                _py  = _spnorm.pdf(_px, m_rho, m_rho_sigma)
                _p1s = np.linspace(max(0.001, m_rho-m_rho_sigma), m_rho+m_rho_sigma, 150)
                _fig_g = go.Figure()
                _fig_g.add_trace(go.Scatter(x=_px, y=_py, mode="lines", fill="tozeroy",
                    line=dict(color="#00b4d8", width=2), fillcolor="rgba(0,180,216,0.13)",
                    hovertemplate="ρ=%{x:.3f}<extra></extra>"))
                _fig_g.add_trace(go.Scatter(
                    x=np.concatenate([_p1s, _p1s[::-1]]),
                    y=np.concatenate([_spnorm.pdf(_p1s, m_rho, m_rho_sigma), np.zeros(150)]),
                    fill="toself", fillcolor="rgba(0,180,216,0.30)", line=dict(width=0),
                    showlegend=False, hoverinfo="skip"))
                _fig_g.add_vline(x=m_rho, line=dict(color="#ffd700", width=1.5, dash="dash"),
                    annotation_text=f"μ={m_rho:.3f}", annotation_font=dict(color="#ffd700", size=11))
                _fig_g.update_layout(height=130, margin=dict(l=42,r=10,t=4,b=28),
                    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                    font=dict(color="white", size=10), showlegend=False,
                    xaxis=dict(title="ρ [g/cm³]", gridcolor="#2a2a3a", zeroline=False),
                    yaxis=dict(visible=False))
                st.plotly_chart(_fig_g, config={"displayModeBar": False})
                st.caption(f"ρ ~ 𝒩({m_rho:.3f}, {m_rho_sigma:.3f}²) g/cm³")
    
            else:  # KDE
                _dens_file = st.file_uploader("Upload density samples CSV", type=["csv"],
                                              key="density_kde_file",
                                              help="One density [g/cm³] per row. Column: density/rho/ρ.")
                if _dens_file is not None:
                    try:
                        _df_k = pd.read_csv(_dens_file)
                        _cmap = {c.lower().strip(): c for c in _df_k.columns}
                        _kcol = next((_cmap[k] for k in ("density","rho","ρ","rho_gcm3","rho_gcm") if k in _cmap),
                                     _df_k.select_dtypes(include="number").columns[0])
                        _samp = _df_k[_kcol].dropna().values.astype(float)
                        _samp = _samp[(_samp > 0.1) & (_samp < 20.0)]
                        if len(_samp) < 3:
                            st.error("⛔  Need ≥ 3 valid density values.")
                            m_rho = st.session_state.get("music_rho", float(default_rho))
                            m_rho_sigma = float(st.session_state.get("music_rho_sigma", 0.0))
                        else:
                            _kde = _spkde(_samp, bw_method="scott")
                            m_rho = float(np.mean(_samp)); m_rho_sigma = float(np.std(_samp, ddof=1))
                            st.session_state["music_rho"] = m_rho
                            st.session_state["music_rho_sigma"] = m_rho_sigma
                            _lo, _hi = max(0.1, m_rho-4*m_rho_sigma), m_rho+4*m_rho_sigma
                            _rx = np.linspace(_lo, _hi, 400)
                            _fig_k = go.Figure()
                            _fig_k.add_trace(go.Scatter(x=_rx, y=_kde(_rx), mode="lines", fill="tozeroy",
                                line=dict(color="#00b4d8", width=2), fillcolor="rgba(0,180,216,0.10)", name="KDE"))
                            _fig_k.add_trace(go.Scatter(x=_rx, y=_spnorm.pdf(_rx, m_rho, m_rho_sigma), mode="lines",
                                line=dict(color="#ffd700", width=1.5, dash="dot"), name=f"𝒩({m_rho:.3f},{m_rho_sigma:.3f})"))
                            _fig_k.add_trace(go.Scatter(x=_samp, y=np.zeros(len(_samp)), mode="markers",
                                marker=dict(color="#ff6b6b", size=6, symbol="line-ns",
                                            line=dict(width=1.5, color="#ff6b6b")), name=f"n={len(_samp)}"))
                            _fig_k.update_layout(height=160, margin=dict(l=42,r=10,t=4,b=28),
                                paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                                font=dict(color="white", size=10),
                                xaxis=dict(title="ρ [g/cm³]", gridcolor="#2a2a3a", zeroline=False),
                                yaxis=dict(visible=False),
                                legend=dict(bgcolor="rgba(0,0,0,0.45)", bordercolor="#444",
                                            borderwidth=1, font=dict(size=9)))
                            st.plotly_chart(_fig_k, config={"displayModeBar": False})
                            st.caption(f"n={len(_samp)}  │  μ={m_rho:.4f}  │  σ={m_rho_sigma:.4f} g/cm³")
                    except Exception as _e:
                        st.error(f"❌  CSV error: {_e}")
                        m_rho = st.session_state.get("music_rho", float(default_rho))
                        m_rho_sigma = float(st.session_state.get("music_rho_sigma", 0.0))
                else:
                    st.info("⬆️  Upload CSV of lab density measurements.")
                    m_rho = st.session_state.get("music_rho", float(default_rho))
                    m_rho_sigma = float(st.session_state.get("music_rho_sigma", 0.0))
    
            # ── Depth + radiation length ───────────────────────────────────────────
            _geo1, _geo2 = st.columns(2)
            m_depth = _geo1.number_input("Overburden depth [m]", 1.0, 10000.0,
                                         st.session_state.get("music_depth_m", 90.0), 5.0,
                                         key="music_depth_m",
                                         help="Rock thickness the transport driver integrates through.")
            m_rad   = _geo2.number_input("Rad. length X₀ [g/cm²]", 0.1, 200.0,
                                         value=st.session_state.get("music_rad", float(default_rad)),
                                         step=0.5, key="music_rad",
                                         help="Radiation length in g/cm² (MUSIC "
                                              "convention). Standard rock: 26.48; "
                                              "water/ice: 36.08; iron: 13.84.")
    
            # ── Opacity metric ─────────────────────────────────────────────────────
            _X_mean  = m_depth * 100.0 * m_rho
            _X_sigma = m_depth * 100.0 * m_rho_sigma
            if m_rho_sigma > 0:
                st.metric("Overburden opacity  X = ρ·L", f"{_X_mean:.1f} g/cm²",
                          delta=f"±{_X_sigma:.1f} g/cm² (1σ)")
                # Opacity uncertainty chart
                _xr  = np.linspace(max(0.0, _X_mean-4*_X_sigma), _X_mean+4*_X_sigma, 300)
                _pX  = _spnorm.pdf(_xr, _X_mean, _X_sigma)
                _x1s = np.linspace(_X_mean-_X_sigma, _X_mean+_X_sigma, 150)
                _fig_x = go.Figure()
                _fig_x.add_trace(go.Scatter(x=_xr, y=_pX, mode="lines", fill="tozeroy",
                    line=dict(color="#90e0ef", width=1.5), fillcolor="rgba(144,224,239,0.10)"))
                _fig_x.add_trace(go.Scatter(
                    x=np.concatenate([_x1s, _x1s[::-1]]),
                    y=np.concatenate([_spnorm.pdf(_x1s, _X_mean, _X_sigma), np.zeros(150)]),
                    fill="toself", fillcolor="rgba(144,224,239,0.25)", line=dict(width=0),
                    showlegend=False, hoverinfo="skip"))
                _fig_x.add_vline(x=_X_mean, line=dict(color="#ffd700", width=1.5, dash="dash"),
                    annotation_text=f"X={_X_mean:.0f}",
                    annotation_font=dict(color="#ffd700", size=10))
                _fig_x.update_layout(height=120, margin=dict(l=42,r=10,t=4,b=24),
                    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                    font=dict(color="white", size=10), showlegend=False,
                    xaxis=dict(title="X [g/cm²]", gridcolor="#2a2a3a", zeroline=False),
                    yaxis=dict(visible=False))
                st.plotly_chart(_fig_x, config={"displayModeBar": False})
            else:
                st.metric("Overburden opacity  X = ρ·L", f"{_X_mean:.1f} g/cm²")
    
            # ── Depth vs detector mismatch warning ────────────────────────────────
            _det_chk = st.session_state.get("gen_detectors", []) if st.session_state.get("gen_use_detector", False) else []
            if _det_chk:
                _d0 = _det_chk[0]
                _det_depth_m = abs(min(_d0.get("az",0), _d0.get("bz",0))) / 100.0 \
                               if _d0["shape"] == 1 \
                               else abs(min(_d0.get("zmin",0), _d0.get("zmax",0))) / 100.0
                if abs(m_depth - _det_depth_m) > 1.0:
                    st.warning(f"⚠️  Overburden ({m_depth:.0f} m) ≠ detector depth ({_det_depth_m:.0f} m). "
                               f"Set overburden = **{_det_depth_m:.0f} m**.")
                else:
                    st.success(f"✅  Overburden matches detector depth ({_det_depth_m:.0f} m)")
    
        # ── BB depth warning ──────────────────────────────────────────────────────
        if (transport_engine == "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS"
                and m_depth > 15.0):
            st.warning(
                f"⚠️  **BB engine — no stochastic fluctuations.** "
                f"At {m_depth:.0f} m the energy-loss distribution is unreliable: "
                f"BB uses continuous CSDA with no Landau/Vavilov spread, so all muons of "
                f"the same initial energy exit with exactly the same final energy. "
                f"For depths > 15 m use **UCMuon-MC**, **MUSIC**, or **PROPOSAL** instead."
            )

        # ── Engine-specific settings ──────────────────────────────────────────────
        st.divider()
        st.markdown("##### Engine Settings")

        # Compute table file status here (needed for init default and expander)
        _eloss_f, _xsec_f, _xs_q = mat_files(mat_choice)
        _mat_suffix_tmp = MUSIC_MATERIALS[mat_choice].get("mat_suffix", "rock")
        _eloss_fallback = _ELOSS_FALLBACK.get(_mat_suffix_tmp, _eloss_f)
        # Driver accepts either the suffixed name (primary) or the plain Kudryavtsev name (fallback)
        # Absolute paths anchored to project root
        _eloss_abs    = _PROJECT_DIR / _eloss_f
        _eloss_fb_abs = _PROJECT_DIR / _eloss_fallback   # music-eloss.dat fallback
        _xsec_abs     = _PROJECT_DIR / _xsec_f
        _mat_suf_tmp  = MUSIC_MATERIALS[mat_choice].get("mat_suffix", "rock")

        # Differential cross-section table (FROM data/ — Kudryavtsev zip, required for init=0)
        _diff_table_names = {
            "rock":     ["music-double-diff-rock.dat"],
            "water":    ["music-double-diff-water.dat"],
            "seawater": ["music-double-diff-water.dat"],
        }
        _diff_ok = any(
            (_PROJECT_DIR / f).exists() or (_PROJECT_DIR / "data" / f).exists()
            for f in _diff_table_names.get(_mat_suf_tmp, [])
        )

        # Computed files (generated by init=0 run) — check root, data/, and bin/
        _eloss_ok = (_eloss_abs.exists() or _eloss_fb_abs.exists()
                     or (_PROJECT_DIR / "data" / _eloss_f).exists()
                     or (_BIN_DIR / _eloss_f).exists())
        _eloss_display_f = _eloss_f if _eloss_abs.exists() else (
            _eloss_fallback if _eloss_fb_abs.exists() else _eloss_f
        )
        _xsec_ok  = (_xsec_abs.exists()
                     or (_PROJECT_DIR / "data" / _xsec_f).exists()
                     or (_BIN_DIR / _xsec_f).exists())

        # files_ok: for init=1 need computed files; for init=0 just need diff table
        _files_ok_init1 = _eloss_ok and _xsec_ok
        _files_ok_init0 = _diff_ok
        _files_ok = _files_ok_init1   # used for init default selection
        _mat_suffix  = MUSIC_MATERIALS[mat_choice].get("mat_suffix", "rock")
        _group_label = {"rock":"rock", "water":"water / ice", "seawater":"seawater"}[_mat_suffix]
        _xs_q_note   = {"approx": (f"⚠️ Iron (Z=26) uses rock XS tables (Z≈11) — approximate. "
                                    f"Suitable for rough estimates only."),
                        "good": None, "exact": None}[_xs_q]
    
        if transport_engine == "MUSIC":
            _sc1, _sc2, _sc3 = st.columns(3)
            with _sc1:
                m_idim = st.radio("Lateral scattering", [1, 0],
                                  format_func=lambda x: {1:"3D — full MS (rec.)", 0:"1D — energy only"}[x])
            with _sc2:
                m_idim1 = st.radio("Other-process scattering", [1, 0],
                                   format_func=lambda x: {1:"ON (recommended)", 0:"OFF (cross-check only)"}[x])
            with _sc3:
                _init_default = 1 if _files_ok_init1 else 0
                m_init = st.radio("Cross-section tables", [0, 1], index=_init_default,
                                  format_func=lambda x: {
                                      0: "init=0 — recalculate (~1 min)",
                                      1: "init=1 — load from disk (fast)"}[x])
                if m_init == 1 and not _files_ok_init1:
                    st.error("⛔  init=1 requires both computed table files (run init=0 first).")
                if m_init == 0 and _files_ok:
                    st.caption("ℹ️  Files exist — will be overwritten.")
    
            # Table files expander
            with st.expander("📂  MUSIC table files", expanded=not _files_ok_init1):
                # Row 1: diff table (input — from Kudryavtsev zip, in data/)
                # Row 2+3: computed files (generated by init=0 run)
                _cf1, _cf2, _cf3 = st.columns(3)
                _diff_name = _diff_table_names.get(_mat_suf_tmp, ["?"])[0]
                _diff_exists = (_PROJECT_DIR / _diff_name).exists() or                                (_PROJECT_DIR / "data" / _diff_name).exists()
                with _cf1:
                    if _diff_exists:
                        st.success(f"✅ `{_diff_name}`  *(input)*")
                    else:
                        st.error(f"❌ `{_diff_name}`  *(input — copy from Kudryavtsev zip)*")
                with _cf2:
                    if _eloss_ok:
                        st.success(f"✅ `{_eloss_display_f}`  *(generated)*")
                    else:
                        st.warning(f"⏳ `{_eloss_f}`  *(generated by init=0)*")
                with _cf3:
                    if _xsec_ok:
                        st.success(f"✅ `{_xsec_f}`  *(generated)*")
                    else:
                        st.warning(f"⏳ `{_xsec_f}`  *(generated by init=0)*")

                st.divider()

                if not _diff_exists:
                    st.error(
                        f"**`{_diff_name}` not found.**\n\n"
                        "This file comes from the Kudryavtsev MUSIC zip and must be present "
                        "in either the project root or the `data/` folder.\n\n"
                        "```bash\n"
                        f"cp path/to/music-double-diff-rock.dat  data/\n"
                        "# or in project root:\n"
                        f"cp path/to/music-double-diff-rock.dat  .\n"
                        "```"
                    )
                elif not _files_ok_init1:
                    st.info(
                        "**First-run setup:**\n\n"
                        f"`{_diff_name}` ✅ is present — "
                        "the driver can compute the energy-loss and cross-section tables.\n\n"
                        "1. Select **init=0 — recalculate** above  \n"
                        "2. Click **▶ Run MUSIC Transport** — the driver will calculate and save "
                        f"`{_eloss_f}` and `{_xsec_f}` automatically (~1 min)  \n"
                        "3. After that, always use **init=1 — load from disk** (fast)"
                    )
                else:
                    st.success("All table files present. Use **init=1** for fast startup.")

                if _xs_q_note:
                    st.warning(_xs_q_note)
    
            with st.expander("⚙️  Expert parameters", expanded=False):
                st.warning("Do not change minv without expert advice.", icon="⚠️")
                m_minv = st.number_input("minv  (10^minv = stochastic/continuous loss cut)",
                                         min_value=-50, max_value=-1, value=-30, step=1)
    
            phitsxs_ms_enable = True   # unused in MUSIC mode
            phitsxs_mat_type  = 1
            phitsxs_custom    = {}
    
        elif transport_engine == "PROPOSAL":
            # ── PROPOSAL settings ─────────────────────────────────────────────────
            _pp1, _pp2, _pp3 = st.columns(3)
            _proposal_medium_map = {
                "Standard Rock  (ρ=2.65 g/cm³)": 1,
                "Water           (ρ=1.00 g/cm³)": 2,
                "Ice             (ρ=0.917 g/cm³)": 3,
                "Seawater        (ρ=1.025 g/cm³)": 4,
                "Custom": 5,
            }
            with _pp1:
                _pp_med_choice = st.selectbox(
                    "PROPOSAL medium", list(_proposal_medium_map.keys()),
                    index=0, key="proposal_med_choice",
                    help="Built-in PROPOSAL media. Custom: define Z, A, ρ, I manually.")
                proposal_medium_type = _proposal_medium_map[_pp_med_choice]
            with _pp2:
                _pp_scat_map = {
                    "HighlandIntegral (recommended)": 2,
                    "Highland (Gaussian approx.)":     1,
                    "Molière (exact)":                 3,
                    "None":                            0,
                }
                _pp_scat_choice = st.selectbox(
                    "Scattering model", list(_pp_scat_map.keys()),
                    index=0, key="proposal_scatter_choice",
                    help=("**HighlandIntegral**: integrated Highland — best accuracy/speed.\n\n"
                          "**Molière**: exact theory — most accurate, ~2× slower.\n\n"
                          "**None**: straight-line tracking."))
                proposal_scatter = _pp_scat_map[_pp_scat_choice]
            with _pp3:
                st.caption("⚛️ **PROPOSAL** — full stochastic MC\n"
                           "ionization · brems · pair · photonuclear · LPM\n"
                           "Ref: Koehne+ (2013) · Alameddine+ (2024)")

            with st.expander("⚙️  Stochastic energy cuts", expanded=False):
                st.caption("Losses > max(e_cut, v_cut·E) are sampled stochastically. "
                           "Defaults (e_cut = 500 MeV, v_cut = 0.001) are conservative; "
                           "the IceCube standard uses v_cut = 0.05 (faster, slightly coarser).")
                _ec1, _ec2 = st.columns(2)
                proposal_e_cut = _ec1.number_input(
                    "e_cut [MeV]", 1.0, 1e5, 500.0, 50.0, key="proposal_e_cut",
                    help="Absolute stochastic energy cut [MeV].")
                proposal_v_cut = _ec2.number_input(
                    "v_cut [fraction]", 1e-4, 1.0, 0.001, 0.0005,
                    format="%.4f", key="proposal_v_cut",
                    help="Fractional stochastic cut v=ΔE/E.")

            with st.expander("📂  Interpolation tables path", expanded=False):
                proposal_tables = st.text_input(
                    "Tables directory", value="",
                    key="proposal_tables_dir",
                    placeholder="~/.proposal/tables  (default)",
                    help="Tables built on first run (~60 s), then cached.")
                if proposal_tables:
                    _pt = Path(proposal_tables).expanduser()
                    st.success(f"✅  `{_pt}`") if _pt.exists() else st.info(f"ℹ️  Will create: `{_pt}`")
                else:
                    _dp = Path("~/.proposal/tables").expanduser()
                    st.caption(f"Default `~/.proposal/tables` — {'✅ exists' if _dp.exists() else '⏳ will be created on first run'}")

            proposal_custom = {}
            if proposal_medium_type == 5:
                st.markdown("**Custom medium parameters**")
                _cm1, _cm2 = st.columns(2)
                proposal_custom["proposal_Z"]    = _cm1.number_input("Z_eff", 1.0, 92.0, 11.0, 0.5, key="prop_Z")
                proposal_custom["proposal_A"]    = _cm2.number_input("A_eff [g/mol]", 1.0, 238.0, 22.0, 1.0, key="prop_A")
                proposal_custom["proposal_rho"]  = _cm1.number_input("ρ [g/cm³]", 0.01, 20.0, 2.65, 0.01, key="prop_rho")
                proposal_custom["proposal_I_eV"] = _cm2.number_input("I_mean [eV]", 10.0, 900.0, 136.4, 1.0, key="prop_Iev")

            _prop_ok, _prop_ver = _check_proposal()
            _driver_src = (_SCRIPT_DIR / "proposal_driver.py").exists()
            if _prop_ok and _driver_src:
                st.success(f"✅  PROPOSAL v{_prop_ver} + `proposal_driver.py` ready")
            elif _driver_src and not _prop_ok:
                st.error(f"❌  `proposal_driver.py` found but PROPOSAL not importable: {_prop_ver}")
                st.code("pip install proposal", language="bash")
            else:
                st.error("❌  `proposal_driver.py` not found in project directory")

            # Unused BB variables (needed by shared run block)
            phitsxs_mat_type  = 1
            phitsxs_ms_enable = True
            phitsxs_custom    = {}
            m_idim = 1; m_idim1 = 1; m_init = 0; m_minv = -30

        elif transport_engine == "UCMuon Stochastic (Python)":
            _stochastic_ok_s, _ = stochastic_available()
            if not _stochastic_ok_s:
                st.error("❌ Place `gui_stochastic_engine.py` and `ucmuon_stochastic_driver.py` in `gui/`.")
            else:
                stochastic_v_cut, stochastic_n_steps, stochastic_ms_enable, \
                    stochastic_range_table, stochastic_hard_spectrum, \
                    stochastic_n_workers, stochastic_delta_rays = render_stochastic_settings()
            stochastic_mat_id = _STOCHASTIC_MAT_ID.get(mat_choice, 1)
            stochastic_custom = {}
            if stochastic_mat_id == 5:
                _puc1, _puc2 = st.columns(2)
                stochastic_custom["stochastic_Z"]    = _puc1.number_input("Z_eff", 1.0, 92.0, 11.0, 0.5, key="stochastic_Zeff")
                stochastic_custom["stochastic_A"]    = _puc2.number_input("A_eff", 1.0,238.0, 22.0, 1.0, key="stochastic_Aeff")
                stochastic_custom["stochastic_I_eV"] = _puc1.number_input("I[eV]",10.0,900.0,136.4, 1.0, key="stochastic_Iev")
                stochastic_custom["stochastic_b_rad"]= _puc2.number_input("b_rad",1e-7, 1e-4,3.475e-6,1e-7,format="%.3e",key="stochastic_brad")
            phitsxs_mat_type=1; phitsxs_ms_enable=True; phitsxs_custom={}
            m_idim=1; m_idim1=1; m_init=0; m_minv=-30
        elif transport_engine == "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS":
            # ── Bethe-Bloch settings ──────────────────────────────────────────────
            _px1, _px2, _px3 = st.columns(3)
            with _px1:
                _pxs_mat_map = {
                    "Standard Rock  (Z=11, A=22, ρ=2.65)": 1,
                    "Limestone       (Z=15.6, A=31.2, ρ=2.71)": 2,
                    "Water / Ice     (Z=7.42, A=14.2, ρ=1.00)": 3,
                    "Iron            (Z=26, A=55.85, ρ=7.87)":  4,
                    "Custom": 5,
                }
                _pxs_choice = st.selectbox("Bethe-Bloch material", list(_pxs_mat_map.keys()),
                                           index=0, key="phitsxs_mat_choice")
                phitsxs_mat_type = _pxs_mat_map[_pxs_choice]
            with _px2:
                phitsxs_ms_enable = st.checkbox(
                    "🔀  Multiple scattering (Highland)", value=True,
                    key="phitsxs_ms_enable",
                    help="Highland formula: θ₀=(13.6/βp)·√(t/X₀)·[1+0.038·ln(t/X₀)]. "
                         "Box-Muller sampling of projected angles. ~10% CPU overhead.")
            with _px3:
                st.caption("💡 Bethe-Bloch (PDG) + Groom radiative b·E.\nNo external table files needed.")
    
            # ── Engine comparison info ─────────────────────────────────────────────
            with st.expander("📐  Bethe-Bloch engine — physics details", expanded=False):
                st.markdown(r"""
**② Bethe-Bloch + Highland MS (this engine)**

**Energy loss:** $dE/dx = a + b\cdot E_{\rm total}$ · CSDA range: $R = \frac{1}{b}\ln(1+\frac{b}{a}E_0)$

| Parameter | Standard Rock | Physical meaning |
|---|---|---|
| $a$ = 1.96 MeV cm²/g | ionisation (Bethe-Bloch plateau) | ≈ constant at muography energies |
| $b$ = 3.64×10⁻⁶ cm²/g | radiative total (brems + pair + photonuclear) | grows linearly with $E$ |

**Multiple scattering (Highland 1979):** applied per step.
$\theta_0 = \frac{13.6\,\text{MeV}}{\beta c p}\sqrt{x/X_0}\,\left[1 + 0.038\ln(x/X_0)\right]$

**Accuracy:** good for 50 GeV – 3 TeV in rock-like media.

**vs MUSIC / PROPOSAL:** BB+MS overestimates survival because it has no Landau/Vavilov
fluctuations — hard stochastic radiative events (which can stop a muon even when mean
$dE/dx$ says it survives) are absent. Typical bias: ~10% at 500 m.w.e., ~20% at 1500 m.w.e.

**Use for:** fast scans, parameter studies, upper-bound estimates.
""")
                st.caption("PDG 2022 §34 · Groom et al. (2001) ADNDT 78, 183 · Highland (1979) NIM 129, 497")
    
            phitsxs_custom = {}
            if phitsxs_mat_type == 5:
                _cx1, _cx2 = st.columns(2)
                phitsxs_custom["phitsxs_Zeff"] = _cx1.number_input("Z_eff", 1.0, 92.0, 11.0, 0.5, key="px_zeff")
                phitsxs_custom["phitsxs_Aeff"] = _cx2.number_input("A_eff", 1.0, 238.0, 22.0, 1.0, key="px_aeff")
                phitsxs_custom["phitsxs_rho"]  = _cx1.number_input("ρ [g/cm³]", 0.01, 20.0, 2.65, 0.01, key="px_rho")
                phitsxs_custom["phitsxs_I_eV"] = _cx2.number_input("I_mean [eV]", 10.0, 900.0, 136.4, 1.0, key="px_Iev")

            m_idim = 1; m_idim1 = 1; m_init = 0; m_minv = -30

        elif transport_engine == "PUMAS":
            # ── PUMAS sub-mode ────────────────────────────────────────────────────
            _pm1, _pm2 = st.columns([2, 3])
            with _pm1:
                pumas_mode = st.radio(
                    "Transport mode", ["backward", "forward"],
                    format_func=lambda x: {"backward": "⬅ Backward MC (no muon file)",
                                           "forward":  "➡ Forward transport"}[x],
                    key="pumas_mode",
                    help="**Backward**: start from detector, propagate backward — 100% efficiency, "
                         "no surface muon file needed.\n\n"
                         "**Forward**: read surface muon file, transport to depth — same as MUSIC/BB.")
            with _pm2:
                if pumas_mode == "backward":
                    st.info("No surface muon file needed — PUMAS samples at the detector "
                            "and propagates backward through rock.")
                else:
                    st.info("Reads the surface muon file selected above. "
                            "Output: standard 18-column underground file.")

            # ── Physics settings ──────────────────────────────────────────────────
            _pl1, _pl2 = st.columns(2)
            with _pl1:
                pumas_energy_loss = st.radio(
                    "Energy-loss mode",
                    [0, 1, 2],
                    format_func=lambda x: {0: "CSDA — deterministic",
                                           1: "Mixed — soft+hard",
                                           2: "Straggled — soft+hard+δe"}[x],
                    key="pumas_energy_loss",
                    help="**CSDA**: fastest, no fluctuations (same as BB engine).\n\n"
                         "**Mixed**: stochastic hard losses above v_cut (Geant4-like).\n\n"
                         "**Straggled**: as Mixed + electronic straggling.")
            with _pl2:
                pumas_scattering = st.radio(
                    "Multiple scattering",
                    [0, 1],
                    format_func=lambda x: {0: "Disabled", 1: "Mixed (Molière)"}[x],
                    key="pumas_scattering",
                    help="Enable Molière / Highland multiple Coulomb scattering. "
                         "Adds lateral displacement to the transported muon.")

            # ── Backward-specific settings ────────────────────────────────────────
            pumas_E_min = pumas_E_max = pumas_theta_max = 85.0
            pumas_n_events = 50000
            pumas_spectrum_id = 0
            pumas_seed = 0
            pumas_outfile = st.session_state.get("pumas_outfile", "output/pumas_flux.dat")

            if pumas_mode == "backward":
                st.markdown("**Backward MC sampling range**")
                _pb1, _pb2, _pb3 = st.columns(3)
                pumas_E_min = _pb1.number_input(
                    "E_det min [GeV]", 0.01, 1e5,
                    float(st.session_state.get("pumas_E_min", 1.0)),
                    step=0.5, key="pumas_E_min",
                    help="Minimum detector kinetic energy sampled.")
                pumas_E_max = _pb2.number_input(
                    "E_det max [GeV]", 0.1, 1e6,
                    float(st.session_state.get("pumas_E_max", 1000.0)),
                    step=10.0, key="pumas_E_max",
                    help="Maximum detector kinetic energy sampled.")
                pumas_theta_max = _pb3.number_input(
                    "θ_max [deg]", 1.0, 89.9,
                    float(st.session_state.get("pumas_theta_max", 85.0)),
                    step=1.0, key="pumas_theta_max",
                    help="Maximum zenith angle (0° = vertical, 90° = horizontal).")

                _pb4, _pb5, _pb6 = st.columns(3)
                pumas_n_events = int(_pb4.number_input(
                    "N events", 100, 10_000_000,
                    int(st.session_state.get("pumas_n_events", 50000)),
                    step=1000, key="pumas_n_events",
                    help="Number of backward MC events. 50 000 gives ~1% statistical error."))
                pumas_spectrum_id = _pb5.selectbox(
                    "Surface spectrum",
                    [0, 1],
                    format_func=lambda x: {0: "GCCLY (Guan et al. 2015)", 1: "Gaisser (PDG)"}[x],
                    key="pumas_spectrum_id",
                    help="Atmospheric muon flux model used to weight the backward-transported muons.")
                pumas_seed = int(_pb6.number_input(
                    "RNG seed (0=random)", 0, 2**31 - 1,
                    int(st.session_state.get("pumas_seed", 0)),
                    step=1, key="pumas_seed",
                    help="0 = time-based seed (different each run)."))

                pumas_outfile = st.text_input(
                    "Flux output file", pumas_outfile, key="pumas_outfile",
                    help="Binned flux spectrum written here. "
                         "Per-event file saved alongside as `<stem>_bwd_events.dat`.")

            # PUMAS material mapping (mat_name for PUMAS MDF)
            _PUMAS_MAT_ID = {
                "Standard Rock": 1, "Limestone": 1, "Rock Salt": 1, "Iron": 1, "Custom": 1,
                "Water": 2, "Ice": 4, "Seawater": 3,
            }
            pumas_mat_id = _PUMAS_MAT_ID.get(mat_choice, 1)

            # Unused variables needed by shared run block
            phitsxs_mat_type = 1; phitsxs_ms_enable = True; phitsxs_custom = {}
            m_idim = 1; m_idim1 = 1; m_init = 0; m_minv = -30

            _pumas_bin_ok2 = (_BIN_DIR / "ucmuon_transport_pumas").exists()
            _pumas_drv_ok2 = (_SCRIPT_DIR / "ucmuon_pumas_driver.py").exists()
            if _pumas_bin_ok2 and _pumas_drv_ok2:
                st.success("✅  `ucmuon_transport_pumas` + `ucmuon_pumas_driver.py` ready")
            elif not _pumas_bin_ok2:
                st.error("❌  Binary not found — run `make pumas`")
            else:
                st.error("❌  `ucmuon_pumas_driver.py` not found in `gui/`")
    
        # ── Underground detector filter ───────────────────────────────────────────
        st.divider()
        _det_ug = st.session_state.get("gen_detectors", []) \
                  if st.session_state.get("gen_use_detector", False) else []
        if _det_ug:
            ug_use_filter = st.checkbox(
                f"🔍  Filter survived muons by detector  ({len(_det_ug)} detector(s) from the Generator tab)",
                value=True)
            ug_filter_file = "output/muons_ug_selected.dat"
            if ug_use_filter:
                ug_filter_file = st.text_input("Underground selected file", "output/muons_ug_selected.dat")
                for _di, _det in enumerate(_det_ug):
                    _sn = "Cylinder" if _det["shape"] == 1 else "Box"
                    st.caption(f"  Detector {_di+1}: {_sn}, margin {_det.get('margin',0):.0f} cm")
        else:
            ug_use_filter  = False
            ug_filter_file = "output/muons_ug_selected.dat"
            st.info("ℹ️  No detector from the Generator tab — enable the detector filter there first.")
    
        # ── Run ───────────────────────────────────────────────────────────────────
        st.divider()
        st.markdown("##### Run")
        _mc1, _mc2 = st.columns([3, 2])
        with _mc1:
            n_threads_music = st.slider(
                "⚡  OpenMP threads", 1, min(64, os.cpu_count() or 8),
                min(4, os.cpu_count() or 4), 1, key="music_omp_threads",
                help="N muons transported in parallel. Each thread has private RANLUX/RANMAR stream.")
        with _mc2:
            _omp_bin     = (_BIN_DIR / "ucmuon_transport_music_omp").exists()
            _serial_bin  = (_BIN_DIR / "ucmuon_transport_music").exists()
            _bbms_bin    = (_BIN_DIR / "ucmuon_transport_bb_omp").exists()
            if transport_engine == "MUSIC":
                st.metric("Mode", "Serial" if n_threads_music == 1 else f"~{n_threads_music}× speedup")
                if _omp_bin:      st.success("✅  `ucmuon_transport_music_omp`")
                elif _serial_bin: st.warning("⚠️  Using serial driver — run `make local` for OMP")
                else:             st.error("❌  No MUSIC driver found — run `make local`")
            elif transport_engine == "PROPOSAL":
                st.metric("Mode", "Single-thread Python")
                _prop_ok2, _prop_ver2 = _check_proposal()
                _driver_py = (_SCRIPT_DIR / "proposal_driver.py").exists()
                if _prop_ok2 and _driver_py:
                    st.success(f"✅  PROPOSAL v{_prop_ver2}")
                elif _driver_py and not _prop_ok2:
                    st.error(f"❌  {_prop_ver2}")
                else:
                    st.error("❌  `proposal_driver.py` not found")
                st.caption("⏳ First run builds tables (~60 s).\nSubsequent runs: < 5 s.")
            elif transport_engine == "UCMuon Stochastic (Python)":
                st.metric("Mode", "Single-thread Python")
                _p_ok3, _p_ver3 = stochastic_available()
                if _p_ok3: st.success(f"✅  UCMuon-MC v{_p_ver3}")
                else:      st.error(f"❌  {_p_ver3}")
            elif transport_engine == "PUMAS":
                _pm = st.session_state.get("pumas_mode", "backward")
                st.metric("Mode", f"PUMAS {'backward MC' if _pm == 'backward' else 'forward'}")
                _pb_ok = (_BIN_DIR / "ucmuon_transport_pumas").exists()
                if _pb_ok:
                    _dump_cached = (_BIN_DIR / "pumas_StandardRock.pumas").exists()
                    st.success("✅  `ucmuon_transport_pumas`")
                    st.caption("⚡ Dump cached — fast startup" if _dump_cached
                               else "⏳ First run builds physics dump (~10 s)")
                else:
                    st.error("❌  Run `make pumas`")
            else:
                st.metric("Mode", "Python (fast analytical CSDA)")
                _bb_py_ok = (_SCRIPT_DIR / "ucmuon_bb_driver.py").exists()
                if _bb_py_ok:
                    st.success("✅  `ucmuon_bb_driver.py`  (MS=OFF: analytical; MS=ON: vectorised)")
                    if _bbms_bin:
                        st.caption("Fortran `ucmuon_transport_bb_omp` also present — Python driver used by default.")
                else:
                    if _bbms_bin:
                        st.warning("⚠️  Python driver missing — using Fortran `ucmuon_transport_bb_omp`")
                    else:
                        st.error("❌  Neither `ucmuon_bb_driver.py` nor `ucmuon_transport_bb_omp` found.")
    
        _rb1, _rb2 = st.columns([4, 1])
        _btn_short = {"MUSIC":"MUSIC","Bethe-Bloch (PDG) + Groom radiative losses + Highland MS":"Bethe-Bloch","PROPOSAL":"PROPOSAL","UCMuon Stochastic (Python)":"UCMuon-MC"}.get(transport_engine, transport_engine)
        _btn_label = f"▶  Run {_btn_short} Transport" 
        run_music  = _rb1.button(_btn_label, type="primary", width='stretch',
                                 disabled=_gg("music_running"))
        stop_music = _rb2.button("⛔  Stop", key="stop_music", width='stretch',
                                 disabled=not _gg("music_running"))
    
        # For PUMAS backward mode there is no surface muon file — use n_events
        _pumas_bwd = (transport_engine == "PUMAS"
                      and st.session_state.get("pumas_mode", "backward") == "backward")
        _n_run = int(st.session_state.get("pumas_n_events", 50000)) if _pumas_bwd else n_transport

        if run_music and not _gg("music_running"):
            if _n_run == 0:
                st.error("No transportable muons found — check the input file.")
            else:
                _source_mode_t2 = st.session_state.get("gen_source_mode", 1)
                _plane_lx_t2    = st.session_state.get("gen_plane_lx",    0.0)
                _plane_ly_t2    = st.session_state.get("gen_plane_ly",    0.0)
                cfg = dict(infile=m_infile, outfile=m_outfile, rho=m_rho, rad=m_rad,
                           depth_m=m_depth, idim=m_idim, idim1=m_idim1, minv=m_minv,
                           init=m_init, source_mode=_source_mode_t2,
                           plane_lx=_plane_lx_t2, plane_ly=_plane_ly_t2,
                           mat_id=MUSIC_MATERIALS[mat_choice].get("mat_id", 1),
                           ncols=ncols_info, transport_all=transport_all)
    
                if transport_engine == "MUSIC":
                    if (_BIN_DIR / "ucmuon_transport_music_omp").exists():
                        _bin = str(_BIN_DIR / "ucmuon_transport_music_omp")
                    elif (_BIN_DIR / "ucmuon_transport_music").exists():
                        _bin = str(_BIN_DIR / "ucmuon_transport_music")
                    else:
                        st.error("❌  No MUSIC driver found — run `make local`."); st.stop()
                    _env = {**os.environ, "OMP_NUM_THREADS": str(n_threads_music)}
                    start_run([_bin], build_music_input(cfg), "music", n_transport, env=_env)

                elif transport_engine == "PROPOSAL":
                    _prop_ok3, _prop_msg3 = _check_proposal()
                    if not _prop_ok3:
                        st.error(f"❌  PROPOSAL not available: {_prop_msg3}"); st.stop()
                    cfg["proposal_medium_type"] = proposal_medium_type
                    cfg["proposal_e_cut"]        = proposal_e_cut
                    cfg["proposal_v_cut"]        = proposal_v_cut
                    cfg["proposal_scatter"]      = proposal_scatter
                    cfg["proposal_tables_dir"]   = proposal_tables
                    cfg.update(proposal_custom)
                    _prop_env = {**os.environ, "PROPOSAL_LOG_LEVEL": "err", "SPDLOG_LEVEL": "err"}
                    start_run(
                        [sys.executable, str(_SCRIPT_DIR / "proposal_driver.py")],
                        build_proposal_input(cfg),
                        "music", n_transport,
                        env=_prop_env,
                    )
    
                elif transport_engine == "UCMuon Stochastic (Python)":
                    _p_ok4, _p_msg4 = stochastic_available()
                    if not _p_ok4:
                        st.error(f"❌ UCMuon-MC unavailable: {_p_msg4}"); st.stop()
                    cfg["stochastic_mat_id"]        = stochastic_mat_id
                    cfg["stochastic_v_cut"]         = stochastic_v_cut
                    cfg["stochastic_n_steps"]       = stochastic_n_steps
                    cfg["stochastic_ms_enable"]     = 1 if stochastic_ms_enable else 0
                    cfg["stochastic_range_table"]   = stochastic_range_table
                    cfg["stochastic_hard_spectrum"] = stochastic_hard_spectrum
                    cfg["stochastic_n_workers"]     = stochastic_n_workers
                    cfg["stochastic_delta_rays"]    = 1 if stochastic_delta_rays else 0
                    cfg.update(stochastic_custom)
                    _pcfg = {**cfg, "infile": _abspath(cfg.get("infile","")), "outfile": _abspath(cfg.get("outfile",""))}
                    start_run([sys.executable, str(_SCRIPT_DIR / "ucmuon_stochastic_driver.py")],
                              build_stochastic_stdin(_pcfg), "music", n_transport, env={**os.environ})

                elif transport_engine == "PUMAS":
                    _pb_bin = _BIN_DIR    / "ucmuon_transport_pumas"
                    _pb_drv = _SCRIPT_DIR / "ucmuon_pumas_driver.py"
                    if not _pb_bin.exists():
                        st.error("❌  `ucmuon_transport_pumas` not found — run `make pumas`."); st.stop()
                    if not _pb_drv.exists():
                        st.error("❌  `ucmuon_pumas_driver.py` not found in `gui/`."); st.stop()
                    _pm_now = st.session_state.get("pumas_mode", "backward")
                    _pumas_outfile_now = (pumas_outfile if _pm_now == "backward"
                                         else m_outfile)
                    cfg["pumas_mode"]        = _pm_now
                    cfg["pumas_mat_id"]      = pumas_mat_id
                    cfg["pumas_energy_loss"] = pumas_energy_loss
                    cfg["pumas_scattering"]  = pumas_scattering
                    cfg["pumas_E_min"]       = pumas_E_min
                    cfg["pumas_E_max"]       = pumas_E_max
                    cfg["pumas_theta_max"]   = pumas_theta_max
                    cfg["pumas_n_events"]    = pumas_n_events
                    cfg["pumas_spectrum_id"] = pumas_spectrum_id
                    cfg["pumas_seed"]        = pumas_seed
                    cfg["outfile"]           = _pumas_outfile_now
                    start_run(
                        [sys.executable, str(_pb_drv)],
                        build_pumas_input(cfg),
                        "music", _n_run,
                        env={**os.environ},
                    )

                else:  # Bethe-Bloch
                    _bb_py  = _SCRIPT_DIR / "ucmuon_bb_driver.py"
                    _bb_bin = _BIN_DIR    / "ucmuon_transport_bb_omp"
                    if not _bb_py.exists() and not _bb_bin.exists():
                        st.error("❌  Neither `ucmuon_bb_driver.py` nor "
                                 "`ucmuon_transport_bb_omp` found."); st.stop()
                    cfg["phitsxs_mat_type"]  = phitsxs_mat_type
                    cfg["phitsxs_ms_enable"] = 1 if phitsxs_ms_enable else 0
                    cfg.update(phitsxs_custom)
                    _bb_stdin = build_phitsxs_input(cfg)
                    # Prefer Python driver (faster, no compilation); fall back to Fortran
                    if _bb_py.exists():
                        _bb_cmd = [sys.executable, str(_bb_py)]
                        _env    = {**os.environ}
                    else:
                        _bb_cmd = [str(_bb_bin)]
                        _env    = {**os.environ, "OMP_NUM_THREADS": str(n_threads_music)}
                    start_run(_bb_cmd, _bb_stdin, "music", n_transport, env=_env)
    
                _ug_file_out = (pumas_outfile
                                if transport_engine == "PUMAS"
                                and st.session_state.get("pumas_mode","backward") == "backward"
                                else m_outfile)
                st.session_state.update({
                    "ug_file":           _ug_file_out,
                    "ug_use_filter":     ug_use_filter,
                    "ug_filter_file":    ug_filter_file,
                    "ug_filter_done":    False,
                    "music_last_mat_id": MUSIC_MATERIALS[mat_choice].get("mat_id", 1),
                })
                save_settings()
                st.rerun()
    
        if stop_music:
            stop_run("music"); st.rerun()
    
        live_panel("music")
    
        # ── Post-processing ────────────────────────────────────────────────────
        st.divider()
        st.markdown("##### 🔧  Post-processing")
        with st.expander("🔄  Convert output → input for chained transport", expanded=False):
            st.caption("Convert the underground output file to a 13-col input for a second transport stage.")
            _conv1, _conv2 = st.columns(2)
            _conv_src = _conv1.text_input("Underground output file", value=m_outfile, key="conv_src")
            _conv_dst = _conv2.text_input("New input filename", "muons_step2_input.dat", key="conv_dst")
            _conv_all = st.checkbox("Include all muons (not just survived)", value=False, key="conv_all")
            if st.button("Convert", key="btn_conv", width='stretch'):
                if not Path(_conv_src).exists():
                    st.error(f"❌  File not found: `{_conv_src}`")
                else:
                    try:
                        def _progress(n_read, n_written):
                            pass
                        convert_ug_to_music_input(_conv_src, _conv_dst,
                                                  survived_only=not _conv_all,
                                                  progress_cb=_progress)
                        st.success(f"✅  Written: `{_conv_dst}`")
                    except Exception as _e:
                        st.error(f"❌  Conversion failed: {_e}")
    
        # ── PHITS underground export ───────────────────────────────────────────────
        with st.expander("⚛️  Export underground muons as PHITS source", expanded=False):
            _pu1, _pu2 = st.columns(2)
            _phits_ug_src = _pu1.text_input("Underground file", value=m_outfile, key="phits_ug_src")
            _phits_ug_out = _pu2.text_input("PHITS output", "output/muons_underground_phits.dat", key="phits_ug_out")
            _phits_ug_all = st.checkbox("Include stopped muons", value=False, key="phits_ug_all")
            if st.button("Convert → PHITS underground", key="btn_phits_ug", width='stretch'):
                if not Path(_phits_ug_src).exists():
                    st.error(f"❌  File not found: `{_phits_ug_src}`")
                else:
                    try:
                        def _phits_progress(n_read, n_written): pass
                        write_phits_underground(_phits_ug_src, _phits_ug_out,
                                               survived_only=not _phits_ug_all,
                                               progress_cb=_phits_progress)
                        st.success(f"✅  Written: `{_phits_ug_out}`")
                        if Path(_phits_ug_out).exists():
                            with open(_phits_ug_out, "rb") as _fh:
                                st.download_button(f"⬇️  Download {_phits_ug_out}", data=_fh,
                                                   file_name=_phits_ug_out, mime="text/plain",
                                                   width='stretch', key="dl_phits_ug")
                        st.code(f"[Source]\n  s-type = 17\n  file   = {_phits_ug_out}\n  dump   = -10\n  1 2 3 4 5 6 7 8 9 10", language="text")
                    except Exception as _e:
                        st.error(f"❌  Failed: {_e}")
    
        _auto_ug_filter()
    
    
    # ══════════════════════════════════════════════════════════════════════════════
    # TAB 3 — RESULTS & VISUALIZATION
    # ══════════════════════════════════════════════════════════════════════════════
with tab_results:

    # ── File selector (compact single row) ────────────────────────────────────
    _pumas_out   = st.session_state.get("pumas_outfile", "output/pumas_flux.dat")
    _pumas_event = str(Path(_pumas_out).with_suffix("")) + "_bwd_events.dat"
    _candidates = [
        "output/muons_surface.dat", "output/muons_selected.dat",
        "output/muons_underground.dat", "output/muons_ug_selected.dat",
        _pumas_out, _pumas_event,
        st.session_state["surface_file"], st.session_state["selected_file"],
        st.session_state["ug_file"], st.session_state.get("ug_filtered_file", ""),
    ]
    available = list(dict.fromkeys(f for f in _candidates if Path(f).exists()))

    _prefer_selected = (
        st.session_state.get("gen_use_detector", False)
        and st.session_state.get("selected_file", "") in available
    )

    if _prefer_selected:
        _default_file = st.session_state["selected_file"]          # ✅ fix 1: was session_state["muons_surface.dat"]
    elif "output/muons_surface.dat" in available:                  # ✅ fix 2: explicit surface priority
        _default_file = "output/muons_surface.dat"
    elif available:
        _default_file = max(available, key=lambda f: Path(f).stat().st_mtime)
    else:
        _default_file = None

    _rfs1, _rfs2 = st.columns([3, 2])
    if available:
        chosen_auto = _rfs1.selectbox("Data file", available,
                                      index=available.index(_default_file) if _default_file in available else 0)
    else:
        chosen_auto = None
        _rfs1.info("ℹ️  No output files found yet. Run the generator →")
    manual_file = _rfs2.text_input("Load any file", value="", placeholder="e.g. muon_rok.dat",
                                   help="Overrides dropdown when filled.")

    if manual_file:
        if Path(manual_file).exists():
            chosen = manual_file
            st.success(f"✅  Using `{manual_file}`")
        else:
            st.error(f"❌  Not found: `{manual_file}`"); st.stop()
    elif chosen_auto:
        chosen = chosen_auto
    else:
        st.stop()

    try:
        _mt = Path(chosen).stat().st_mtime if Path(chosen).exists() else 0
        df  = load_file(chosen, mtime=_mt)
    except Exception as _ex:
        st.error(f"Could not load: {_ex}"); st.stop()

    # ── Classify the loaded file ───────────────────────────────────────────────
    ug_filter_file  = st.session_state.get("ug_filtered_file", "output/muons_ug_selected.dat")
    is_det_hits     = (chosen == ug_filter_file or "ug_selected" in chosen)
    is_underground  = ("alive" in df.columns and "xs" in df.columns and not is_det_hits)
    is_pumas_flux   = ("E_det_GeV" in df.columns and "flux" in df.columns)
    is_pumas_events = ("E_det_GeV" in df.columns and "flux_contribution" in df.columns)
    ecol            = "E" if "E" in df.columns else "Es"

    # Pre-compute survival / charge so they're available both in the quick bar
    # and inside the detailed expander without duplicating lookup logic.
    _ev          = df[ecol][df[ecol] > 0] if ecol in df.columns else pd.Series([], dtype=float)
    _s_survived  = int((df["alive"] == 1).sum()) if is_underground else None
    _n_surv_det  = st.session_state.get("music_nmuons_survived", None)
    _has_charge  = "charge" in df.columns
    if _has_charge:
        _np_c = int((df["charge"] == 1).sum())
        _nm_c = int((df["charge"] == -1).sum())
    rate_per_s, _, flux_ok = _compute_flux(len(df))
    _det_t3 = (st.session_state.get("gen_detectors", [])
               if st.session_state.get("gen_use_detector", False) else [])

    # ── Quick-glance bar (always visible — 4 most important metrics) ─────────
    _qb1, _qb2, _qb3, _qb4 = st.columns(4)
    if is_pumas_flux:
        _valid_bins = df[df["n_events_in_bin"] > 0]
        _qb1.metric("E_det bins (non-empty)", f"{len(_valid_bins):,}")
        _qb2.metric("E_det range [GeV]",
                    f"{df['E_det_GeV'].min():.2f} – {df['E_det_GeV'].max():.2f}")
        _pumas_peak = df.loc[df["flux"].idxmax(), "E_det_GeV"] if len(df) > 0 else 0
        _qb3.metric("Peak flux E_det [GeV]", f"{_pumas_peak:.2f}")
        _qb4.metric("Total events in file",  f"{int(df['n_events_in_bin'].sum()):,}")
    elif is_pumas_events:
        _qb1.metric("PUMAS backward events", f"{len(df):,}")
        _qb2.metric("Mean E_det [GeV]",      f"{df['E_det_GeV'].mean():.2f}")
        _qb3.metric("Mean E_surf [GeV]",     f"{df['E_surf_GeV'].mean():.2f}")
        _n_p = int((df["charge"] == 1).sum()); _n_m = int((df["charge"] == -1).sum())
        _qb4.metric("μ⁺/μ⁻", f"{_n_p/max(_n_m,1):.3f}")
    else:
        _qb1.metric("Muons" if not is_det_hits else "Detector hits", f"{len(df):,}")
        if len(_ev) > 0:
            _elabel = "Mean E at depth [GeV]" if (is_underground or is_det_hits) else "Mean E [GeV]"
            _qb2.metric(_elabel, f"{_ev.mean():.2f}")
        if is_underground and _s_survived is not None:
            _qb3.metric("Survived",      f"{_s_survived:,}")
            _qb4.metric("Survival rate", f"{100*_s_survived/max(len(df),1):.1f}%")
        elif is_det_hits and _n_surv_det:
            _qb3.metric("Survived (transport)", f"{_n_surv_det:,}")
            _qb4.metric("Detector hit rate",    f"{100*len(df)/max(_n_surv_det,1):.1f}%")
        elif flux_ok and rate_per_s is not None:
            _sigma = rate_per_s / np.sqrt(max(len(df), 1))
            _qb3.metric("Rate [/s]",   f"{rate_per_s:.4g}", delta=f"±{_sigma:.3g}")
            _qb4.metric("Rate [/day]", f"{rate_per_s*86400:.4g}")

    # ── Detailed statistics (collapsed — expand for full breakdown) ───────────
    with st.expander("📊  Detailed statistics", expanded=False):
        _ds_l, _ds_r = st.columns(2)

        with _ds_l:
            st.markdown("**Energy & survival**")
            _dsa1, _dsa2, _dsa3 = st.columns(3)
            _dsa1.metric("Total muons", f"{len(df):,}")
            if len(_ev) > 0:
                _dsa2.metric("Mean E [GeV]", f"{_ev.mean():.2f}")
                _dsa3.metric("Max E [GeV]",  f"{_ev.max():.2f}")
            if is_underground and _s_survived is not None:
                _dsb1, _dsb2, _dsb3 = st.columns(3)
                _dsb1.metric("Survived",      f"{_s_survived:,}")
                _dsb2.metric("Stopped",       f"{len(df)-_s_survived:,}")
                _dsb3.metric("Survival rate", f"{100*_s_survived/max(len(df),1):.1f}%")
            elif is_det_hits and _n_surv_det:
                _dsb1, _dsb2 = st.columns(2)
                _dsb1.metric("Hit / survived",  f"{len(df):,} / {_n_surv_det:,}")
                _dsb2.metric("Detector hit rate", f"{100*len(df)/max(_n_surv_det,1):.1f}%")

        with _ds_r:
            if _has_charge:
                st.markdown("**Charge composition**")
                _cr1, _cr2, _cr3 = st.columns(3)
                _cr1.metric("μ⁺  (+1)", f"{_np_c:,}")
                _cr2.metric("μ⁻  (−1)", f"{_nm_c:,}")
                _cr3.metric("μ⁺/μ⁻", f"{_np_c/max(_nm_c,1):.3f}",
                            help="Expected ~1.27 at sea level.")
            if _det_t3:
                st.markdown("**Detector solid angle**")
                _sa_sr, _ca_sr, _sa_msr, _frac = compute_detector_solid_angle(_det_t3)
                _dg1, _dg2, _dg3, _dg4 = st.columns(4)
                _dg1.metric("Solid angle",      f"{_sa_sr:.4e} sr")
                _dg2.metric("Solid angle",      f"{_sa_msr:.4f} msr")
                _dg3.metric("Fraction of 2π",   f"{_frac*100:.4f} %")
                _dg4.metric("cos²θ acceptance", f"{_ca_sr:.4e} sr")
                st.caption("MC estimate — 600k rays from disk centre.")

        if flux_ok and rate_per_s is not None:
            st.markdown("**Estimated detector rate**")
            _sigma = rate_per_s / np.sqrt(max(len(df), 1))
            _rc1, _rc2, _rc3, _rc4 = st.columns(4)
            _rc1.metric("Rate [/s]",   f"{rate_per_s:.4g}",   delta=f"±{_sigma:.3g}")
            _rc2.metric("Rate [/min]", f"{rate_per_s*60:.4g}")
            _rc3.metric("Rate [/h]",   f"{rate_per_s*3600:.4g}")
            _rc4.metric("Rate [/day]", f"{rate_per_s*86400:.4g}")
            st.caption(
                f"Φ = {st.session_state.get('gen_integrated_flux',0):.4e} cm⁻²s⁻¹sr⁻¹  |  "
                f"R = {st.session_state.get('gen_radius',0):.0f} m  |  "
                f"N_gen = {st.session_state.get('gen_nmuons_done',0):,}")

    st.divider()

    # ── PUMAS backward MC flux spectrum analysis ───────────────────────────────
    if is_pumas_flux:
        st.markdown("##### PUMAS Backward MC — Flux Spectrum")
        _pf_l, _pf_r = st.columns(2)

        with _pf_l:
            _valid = df[df["n_events_in_bin"] > 0].copy()
            _fig_flux = go.Figure()
            _fig_flux.add_trace(go.Scatter(
                x=_valid["E_det_GeV"].tolist(),
                y=_valid["flux"].tolist(),
                error_y=dict(type="data", array=_valid["flux_err"].tolist(),
                             visible=True, color="#aaaaaa", thickness=1),
                mode="lines+markers",
                line=dict(color="#00b4d8", width=2),
                marker=dict(size=4),
                name="dΦ/dE·dΩ",
            ))
            _fig_flux.update_layout(
                xaxis=dict(title="E_det [GeV]", type="log", gridcolor="#2a2a3a"),
                yaxis=dict(title="dΦ/dE·dΩ [m⁻²s⁻¹GeV⁻¹sr⁻¹]", type="log", gridcolor="#2a2a3a"),
                paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                font=dict(color="white"),
                margin=dict(l=70, r=20, t=30, b=50), height=320,
                title=dict(text="Underground Flux Spectrum", font=dict(color="white")),
            )
            st.plotly_chart(_fig_flux)

        with _pf_r:
            _fig_esurf = go.Figure()
            _fig_esurf.add_trace(go.Scatter(
                x=_valid["E_det_GeV"].tolist(),
                y=_valid["E_surf_mean_GeV"].tolist(),
                mode="lines+markers",
                line=dict(color="#ff6b6b", width=2),
                marker=dict(size=4),
                name="Mean E_surf",
            ))
            _fig_esurf.update_layout(
                xaxis=dict(title="E_det [GeV]", type="log", gridcolor="#2a2a3a"),
                yaxis=dict(title="Mean E_surface [GeV]", type="log", gridcolor="#2a2a3a"),
                paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                font=dict(color="white"),
                margin=dict(l=70, r=20, t=30, b=50), height=320,
                title=dict(text="Required Surface Energy per Bin", font=dict(color="white")),
            )
            st.plotly_chart(_fig_esurf)

        with st.expander("🗃️  Flux data table", expanded=False):
            st.dataframe(_valid.rename(columns={
                "E_det_GeV": "E_det [GeV]",
                "flux": "dΦ/dE·dΩ [m⁻²s⁻¹GeV⁻¹sr⁻¹]",
                "flux_err": "flux_err",
                "E_surf_mean_GeV": "E_surf [GeV]",
                "n_events_in_bin": "N events",
            }).reset_index(drop=True), width='stretch')
        st.divider()

    elif is_pumas_events:
        st.markdown("##### PUMAS Backward MC — Per-Event Data")
        _pe_l, _pe_r = st.columns(2)

        with _pe_l:
            _e_vals = df["E_det_GeV"].values
            _e_vals = _e_vals[(_e_vals > 0) & np.isfinite(_e_vals)]
            _fig_pe = go.Figure()
            if len(_e_vals):
                _eh, _eb = np.histogram(np.log10(_e_vals), bins=50)
                _ec = 0.5 * (_eb[:-1] + _eb[1:])
                _fig_pe.add_trace(go.Bar(x=(10**_ec).tolist(), y=_eh.tolist(),
                                         marker_color="#00b4d8", name="E_det"))
            _fig_pe.update_layout(
                xaxis=dict(title="E_det [GeV]", type="log", gridcolor="#2a2a3a"),
                yaxis=dict(title="Counts", gridcolor="#2a2a3a"),
                paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                font=dict(color="white"), margin=dict(l=60,r=20,t=30,b=50), height=280,
                title=dict(text="E_det distribution", font=dict(color="white")),
            )
            st.plotly_chart(_fig_pe)

        with _pe_r:
            _ct_vals = df["cos_theta"].values
            _fig_ct = go.Figure()
            _cth, _ctb = np.histogram(_ct_vals[np.isfinite(_ct_vals)], bins=40)
            _ctc = 0.5 * (_ctb[:-1] + _ctb[1:])
            _fig_ct.add_trace(go.Bar(x=_ctc.tolist(), y=_cth.tolist(),
                                     marker_color="#ff6b6b", name="cos_theta"))
            _fig_ct.update_layout(
                xaxis=dict(title="cos θ", gridcolor="#2a2a3a"),
                yaxis=dict(title="Counts", gridcolor="#2a2a3a"),
                paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                font=dict(color="white"), margin=dict(l=60,r=20,t=30,b=50), height=280,
                title=dict(text="cos θ distribution", font=dict(color="white")),
            )
            st.plotly_chart(_fig_ct)

        with st.expander("🗃️  Data preview (first 10 rows)", expanded=False):
            st.dataframe(df.head(10), width='stretch')
        st.divider()

    # ── Analysis sections in two parallel columns (skip for PUMAS-specific files) ─
    _show_generic = not (is_pumas_flux or is_pumas_events)
    if _show_generic:
        st.markdown("##### Analysis")
    _res_l, _res_r = st.columns(2)

    with _res_l:
        if _show_generic:
            # ── Distributions ─────────────────────────────────────────────────
            st.markdown("**Distributions**")
            ncols_list = [c for c in df.columns
                          if c not in ("EventID","charge","hit_flag","det_mask","alive")]
            if not ncols_list:
                st.warning("⚠️  No plottable columns.")
            else:
                pc = st.selectbox("Variable", ncols_list,
                                  index=ncols_list.index("theta") if "theta" in ncols_list else 0)
                if pc and pc in df.columns:
                    _has_ecol = ecol in df.columns
                    _ncols_fig = 2 if _has_ecol else 1
                    fig2, _axes = plt.subplots(1, _ncols_fig, figsize=(11 if _has_ecol else 5.5, 3.5))
                    axes = _axes if _has_ecol else [None, _axes]
                    fig2.patch.set_facecolor("#0e1117")
                    for ax in ([axes[0]] if _has_ecol else []) + [axes[1]]:
                        ax.set_facecolor("#1c1e26"); ax.tick_params(colors="white")
                        ax.xaxis.label.set_color("white"); ax.yaxis.label.set_color("white")
                        ax.title.set_color("white")
                        for sp in ax.spines.values(): sp.set_edgecolor("#444")
                    if _has_ecol:
                        _evals = df[ecol].values
                        _evals = _evals[(_evals > 0) & np.isfinite(_evals)]
                        if len(_evals):
                            _le = np.log10(_evals)
                            _eh, _eb = np.histogram(_le, bins=60)
                            axes[0].bar(_eb[:-1], _eh, width=np.diff(_eb),
                                        color="#00b4d8", edgecolor="none", alpha=0.85, align="edge")
                        axes[0].set_xlabel("log₁₀(E [GeV])", color="white")
                        axes[0].set_ylabel("Counts", color="white")
                        axes[0].set_title("Energy spectrum", color="white")
                    _raw = df[pc].dropna().values
                    _xlabel = pc
                    if pc in ("theta","theta_s","phi","phi_s"):
                        _raw = np.degrees(_raw); _xlabel = f"{pc} [deg]"
                    elif pc in ("x","y","z","xs","ys","zs"):
                        _raw = _raw / 100.0; _xlabel = f"{pc} [m]"
                    if len(_raw):
                        _dh, _db = np.histogram(_raw, bins=60)
                        axes[1].bar(_db[:-1], _dh, width=np.diff(_db),
                                    color="#ff6b6b", edgecolor="none", alpha=0.85, align="edge")
                    axes[1].set_xlabel(_xlabel, color="white")
                    axes[1].set_ylabel("Counts", color="white")
                    axes[1].set_title(f"Distribution: {pc}", color="white")
                    if pc in ("theta","theta_s") and len(_raw) > 0:
                        axes[1].set_xlim(left=0)
                        axes[1].axvline(35.26, color="#ffd700", lw=1.2, ls="--", alpha=0.7,
                                        label="cos²θ·sinθ peak (35.3°)")
                        axes[1].legend(fontsize=8, facecolor="#1c1e26", labelcolor="white", framealpha=0.7)
                    plt.tight_layout()
                    st.pyplot(fig2); plt.close(fig2)

        # ── Angular acceptance map ─────────────────────────────────────────────
        if _show_generic and "theta" in df.columns and "phi" in df.columns:
            with st.expander("🧭  Angular acceptance map  (θ vs φ)", expanded=False):
                _theta_deg = np.degrees(df["theta"].dropna().values)
                # Wrap phi into [-180, 180): the Fortran generator and MUSIC
                # write phi in [0, 2pi), the Python engines in [-pi, pi] —
                # without the wrap, half the muons fall outside the histogram.
                _phi_deg   = (np.degrees(df["phi"].dropna().values)
                              + 180.0) % 360.0 - 180.0
                # Pre-compute 2D histogram server-side: sends 36×18 bin counts
                # instead of 10M raw data points — avoids websocket size limit.
                _H, _phi_edges, _th_edges = np.histogram2d(
                    _phi_deg, _theta_deg,
                    bins=[36, 18],
                    range=[[-180, 180], [0, 90]],
                )
                _phi_centres = 0.5 * (_phi_edges[:-1] + _phi_edges[1:])
                _th_centres  = 0.5 * (_th_edges[:-1]  + _th_edges[1:])
                _fig_ang = go.Figure(go.Heatmap(
                    x=_phi_centres.tolist(), y=_th_centres.tolist(),
                    z=_H.T.tolist(),
                    colorscale="Plasma",
                    colorbar=dict(title=dict(text="Counts", font=dict(color="white")),
                                  tickfont=dict(color="white"))))
                _fig_ang.update_layout(
                    xaxis_title="φ [deg]", yaxis_title="θ [deg]",
                    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
                    font=dict(color="white"),
                    xaxis=dict(gridcolor="#2a2a3a"), yaxis=dict(gridcolor="#2a2a3a"),
                    margin=dict(l=60,r=30,t=30,b=50), height=320,
                    title=dict(text="Angular Distribution", font=dict(color="white")))
                st.plotly_chart(_fig_ang)

        # ── Data preview ───────────────────────────────────────────────────────
        if _show_generic:
            with st.expander("🗃️  Data preview  (first 10 rows)", expanded=False):
                st.dataframe(df.head(10), width='stretch')

    with _res_r:
        # ── Survival rate vs depth ─────────────────────────────────────────────
        if "alive" in df.columns and "Es" in df.columns:
            st.markdown("##### Survival rate vs depth")
            _sd1, _sd2 = st.columns(2)
            _depth_m_val = _sd1.number_input("Overburden depth [m]", 1.0, 10000.0,
                                             float(st.session_state.get("music_depth_m", 90.0)),
                                             5.0, key="res_depth")
            _rho_val     = _sd2.number_input("Rock density [g/cm³]", 0.1, 10.0,
                                             float(st.session_state.get("music_rho", 2.65)),
                                             0.05, key="res_rho")
            st.plotly_chart(plot_survival_vs_depth(df, _depth_m_val, _rho_val))
            st.caption(
                "**Blue curve** — analytical CSDA (Groom 2001). "
                "**◆ Diamond** — transport MC survival rate. "
                "CSDA is an upper bound; stochastic engines include straggling and radiative losses. "
                "Gap >5% may indicate E_min too low or (MUSIC) XS tables needing recalculation (init=0).")


    # ── 3D trajectories ───────────────────────────────────────────────────────
    st.divider()
    _det_3d    = st.session_state.get("gen_detectors",[]) if st.session_state.get("gen_use_detector",False) else []
    _radius_3d = float(st.session_state.get("gen_radius", 800.0))
    _src_3d    = int(st.session_state.get("gen_source_mode", 1))
    if _det_3d:
        _d0 = _det_3d[0]
        _z_min_det = min(_d0.get("az",0), _d0.get("bz",0)) if _d0["shape"]==1 \
                     else min(_d0.get("zmin",0), _d0.get("zmax",0))
        auto_depth = abs(_z_min_det) / 100.0
    else:
        auto_depth = float(st.session_state.get("music_depth_m", 90.0))

    with st.expander("🧊  3D Muon Trajectories", expanded=False):
        if "alive" in df.columns and "xs" in df.columns:
            _n_traj = st.slider("Trajectories per class", 10, 500, 100, 10,
                                help="N survived + N stopped. Higher = slower browser.")
            with st.expander("🔴  Overlay detector geometry", expanded=False):
                if len(_det_3d) > 0:
                    st.success(f"✅  {len(_det_3d)} detector(s) from the Generator tab.")
                else:
                    st.info("No detector from the Generator tab — define one here for overlay.")
                    _ov_shape = st.selectbox("Shape", [1,2], key="ov_sh",
                                             format_func=lambda x: "Cylinder" if x==1 else "Box")
                    _ov_det = {"shape": _ov_shape, "margin": 0.0}
                    if _ov_shape == 1:
                        _oc1, _oc2 = st.columns(2)
                        with _oc1:
                            _ov_det["ax"]=st.number_input("Ax",value=0.0,key="ov_ax")
                            _ov_det["ay"]=st.number_input("Ay",value=0.0,key="ov_ay")
                            _ov_det["az"]=st.number_input("Az",value=-float(auto_depth*100),key="ov_az")
                        with _oc2:
                            _ov_det["bx"]=st.number_input("Bx",value=0.0,key="ov_bx")
                            _ov_det["by"]=st.number_input("By",value=0.0,key="ov_by")
                            _ov_det["bz"]=st.number_input("Bz",value=0.0,key="ov_bz")
                        _ov_det["r"]=st.number_input("Radius [cm]",0.1,1e4,5.0,key="ov_r")
                    else:
                        _oc1, _oc2 = st.columns(2)
                        with _oc1:
                            _ov_det["xmin"]=st.number_input("Xmin",value=-100.0,key="ov_xn")
                            _ov_det["ymin"]=st.number_input("Ymin",value=-100.0,key="ov_yn")
                            _ov_det["zmin"]=st.number_input("Zmin",value=-float(auto_depth*100),key="ov_zn")
                        with _oc2:
                            _ov_det["xmax"]=st.number_input("Xmax",value=100.0,key="ov_xx")
                            _ov_det["ymax"]=st.number_input("Ymax",value=100.0,key="ov_yx")
                            _ov_det["zmax"]=st.number_input("Zmax",value=0.0,key="ov_zx")
                    if st.button("Apply overlay", key="ov_apply"):
                        _det_3d = [_ov_det]
                        st.session_state["gen_detectors"] = _det_3d
                        st.session_state["gen_use_detector"] = True
            st.plotly_chart(plot_3d_trajectories(df, _n_traj, auto_depth,
                                                  detectors=_det_3d, radius_m=_radius_3d))
            st.caption("🟢 Generation surface  🔵 Survived  🔴 Stopped  🟡 Detector")

        elif "theta" in df.columns and ("x" in df.columns or "xs" in df.columns):
            if _det_3d:
                st.success(f"✅  Detector from the Generator tab — extending to {auto_depth:.1f} m.")
            else:
                st.info("ℹ️  No detector from the Generator tab — showing direction lines only.")
            st.plotly_chart(plot_3d_surface(
                df, _radius_3d, detectors=_det_3d, source_mode=_src_3d,
                disk_cx=float(st.session_state.get("gen_disk_cx", 0.0)),
                disk_cy=float(st.session_state.get("gen_disk_cy", 0.0)),
                src_w_m=float(st.session_state.get("gen_src_w_m", 0.0)),
                disk_tilt=float(st.session_state.get("gen_disk_tilt", 0.0)),
                disk_tilt_az=float(st.session_state.get("gen_disk_tilt_az", 0.0)),
                source_plane=int(st.session_state.get("gen_source_plane", 1)),
            ))
            _lh = ("🟢 Hits detector  🔴 Misses" if "hit_flag" in df.columns
                   else "🔵 Selected trajectories")
            st.caption(f"🟢 Generation surface  {_lh}")
        else:
            st.info("3D trajectories require a surface, selected, or underground file.")

    # ── Export ────────────────────────────────────────────────────────────────
    with st.expander("⬇️  Export", expanded=False):
        _summary = {
            "timestamp":       time.strftime("%Y-%m-%d %H:%M"),
            "file_visualised": chosen, "N_muons": len(df),
            "N_generated":     st.session_state.get("gen_nmuons_done",   "—"),
            "gen_radius_m":    st.session_state.get("gen_radius",        "—"),
            "overburden_m":    st.session_state.get("music_depth_m",     "—"),
            "rock_density":    st.session_state.get("music_rho",         "—"),
            "N_survived":      st.session_state.get("music_nmuons_survived","—"),
            "rate_per_s":      f"{rate_per_s:.6g}" if flux_ok else "—",
            "solid_angle_sr":  "—",
        }
        _det_exp = st.session_state.get("gen_detectors",[]) if st.session_state.get("gen_use_detector",False) else []
        if _det_exp:
            _sa_e, _ca_e, _, _ = compute_detector_solid_angle(_det_exp)
            _summary["solid_angle_sr"] = f"{_sa_e:.6e}"
            _summary["cos2_acceptance_sr"] = f"{_ca_e:.6e}"
        _exp1, _exp2 = st.columns(2)
        _exp1.download_button("⬇️  Run summary (CSV)", data=pd.DataFrame([_summary]).to_csv(index=False),
                              file_name="ucmuon_summary.csv", mime="text/csv", width='stretch')
        with open(chosen, "rb") as _fh:
            _exp2.download_button(f"⬇️  {chosen}", data=_fh, file_name=Path(chosen).name,
                                  mime="text/plain", width='stretch')

    st.divider()
    st.caption("🌌 **UCMuon** — UCLouvain Muography Group | "
               "Hamid Basiri · [hamid.basiri@uclouvain.be](mailto:hamid.basiri@uclouvain.be) | "
               "MIT License 2026")


with tab_config:
    # ── Autosave status ───────────────────────────────────────────────────────
    import datetime as _dt
    if Path(AUTOSAVE_FILE).exists():
        _mtime = Path(AUTOSAVE_FILE).stat().st_mtime
        _saved_at = _dt.datetime.fromtimestamp(_mtime).strftime("%Y-%m-%d  %H:%M:%S")
        st.success(f"🕒  Last autosave: `{AUTOSAVE_FILE}`  —  {_saved_at}")
    else:
        st.info("ℹ️  No autosave file found yet — run the generator to create one.")
    st.caption(
        "💾  Settings are auto-saved to `ucmuon_autosave.json` on every run.  "
        "Upload a saved JSON below to restore a previous session."
    )

    st.divider()

    # ── Upload / Restore ──────────────────────────────────────────────────────
    with st.expander("🔄  Restore from JSON", expanded=False):
        _up = st.file_uploader(
            "Upload `ucmuon_config.json` or `ucmuon_autosave.json`",
            type=["json"], key="config_upload",
        )
        if _up is not None:
            # Widget keys cannot be written after their widgets rendered this
            # run (StreamlitAPIException), so stash the values and rerun —
            # _apply_pending_restore() applies them before any widget exists.
            _up_id = getattr(_up, "file_id", None) or f"{_up.name}:{getattr(_up, 'size', 0)}"
            if st.session_state.get("_cfg_restored_id") != _up_id:
                try:
                    _data = json.load(_up)
                    # Only restore scalar / list values; skip non-serialisable objects
                    _restored = {k: v for k, v in _data.items()
                                 if isinstance(v, (str, int, float, bool, list, dict, type(None)))}
                    st.session_state["_cfg_restore_pending"] = _restored
                    st.session_state["_cfg_restored_id"]     = _up_id
                    st.rerun()
                except Exception as _ex:
                    st.error(f"❌  Could not parse JSON: {_ex}")
            else:
                _n_restored = st.session_state.get("_cfg_restore_msg")
                if _n_restored is not None:
                    st.success(
                        f"✅  Restored **{_n_restored}** settings from `{_up.name}` — "
                        "widgets now reflect the restored values."
                    )

    st.divider()

    # ── Live settings snapshot ────────────────────────────────────────────────
    gen_detectors = (st.session_state.get("gen_detectors", [])
                     if st.session_state.get("gen_use_detector", False) else [])
    cfg_display = {
        "generator": {
            "e_min_GeV":       st.session_state.get("gen_emin",           "— run generator first"),
            "e_max_GeV":       st.session_state.get("gen_emax",           "— run generator first"),
            "spectrum_mode":   st.session_state.get("gen_spectrum_mode",  "— run generator first"),
            "radius_m":        st.session_state.get("gen_radius",         "— run generator first"),
            "n_muons_saved":   st.session_state.get("gen_nmuons_done",    "— run generator first"),
            "theta_max_deg":   st.session_state.get("gen_theta_max",      "— run generator first"),
            "angular_mode":    st.session_state.get("gen_angular_mode",   "— run generator first"),
            "use_detector":    st.session_state.get("gen_use_detector",   False),
            "surface_file":    st.session_state.get("surface_file",       "—"),
            "selected_file":   st.session_state.get("selected_file",      "—"),
            "integrated_flux": st.session_state.get("gen_integrated_flux","— run generator first"),
        },
        "transport": {
            "engine":           st.session_state.get("transport_engine",         "MUSIC"),
            "overburden_m":     st.session_state.get("music_depth_m",            "— run transport first"),
            "rock_density":     st.session_state.get("music_rho",                "— run transport first"),
            "radiation_length": st.session_state.get("music_rad",                "— run transport first"),
            "output_file":      st.session_state.get("ug_file",                  "—"),
            "n_transported":    st.session_state.get("music_nmuons_transported", "— run transport first"),
            "n_survived":       st.session_state.get("music_nmuons_survived",    "— run transport first"),
            "ug_filtered_file": st.session_state.get("ug_filtered_file",         "—"),
        },
        "detectors": gen_detectors,
    }

    _cc1, _cc2 = st.columns([2, 1])
    with _cc1:
        st.json(cfg_display)
    with _cc2:
        st.markdown("**Quick stats**")
        if st.session_state.get("gen_nmuons_done"):
            st.metric("Generated",   f"{st.session_state['gen_nmuons_done']:,}")
        if st.session_state.get("music_nmuons_transported"):
            st.metric("Transported", f"{st.session_state['music_nmuons_transported']:,}")
        if st.session_state.get("music_nmuons_survived"):
            _ntr = st.session_state.get("music_nmuons_transported", 1)
            _ns  = st.session_state["music_nmuons_survived"]
            st.metric("Survived",    f"{_ns:,}", delta=f"{100*_ns/max(_ntr,1):.1f}%")
        st.divider()
        config_str = json.dumps(cfg_display, indent=2, default=str)
        st.download_button("⬇️  Download config JSON", data=config_str,
                           file_name="ucmuon_config.json", mime="application/json",
                           width='stretch')
        if Path(AUTOSAVE_FILE).exists():
            with open(AUTOSAVE_FILE, "rb") as _fh:
                st.download_button("💾  Download autosave", data=_fh,
                                   file_name=AUTOSAVE_FILE, mime="application/json",
                                   width='stretch', key="dl_autosave")
        st.divider()
        if st.button("🗑️  Reset autosave", width='stretch',
                     help="Delete the autosave file so the app starts fresh next time."):
            try:
                Path(AUTOSAVE_FILE).unlink(missing_ok=True)
                st.success("Autosave deleted.")
            except Exception as _ex:
                st.error(f"Could not delete: {_ex}")
# ══════════════════════════════════════════════════════════════════════════════
# TAB 5 — TERRAIN (DEM-AWARE TRANSPORT)
# Note: save_settings() is called AFTER this block so terrain widget states
# are saved correctly (they render inside this block and claim their keys here).
# ══════════════════════════════════════════════════════════════════════════════
with tab_terrain:
    if _TERRAIN_GUI_OK and render_terrain_tab is not None:
        try:
            render_terrain_tab(
                script_dir              = _SCRIPT_DIR,
                project_dir             = _PROJECT_DIR,
                build_music_input_fn    = build_music_input,
                build_phitsxs_input_fn  = build_phitsxs_input,
                build_proposal_input_fn = build_proposal_input,
                music_materials         = MUSIC_MATERIALS,
                probe_music_file_fn     = probe_music_file,
                load_file_fn            = load_file,
            )
        except Exception as _terrain_exc:
            import traceback as _terrain_tb
            st.error(
                f"❌  Terrain tab error: `{_terrain_exc}`\n\n"
                "Full traceback shown below — please report this.",
                icon="❌"
            )
            st.code(_terrain_tb.format_exc())
    elif _TERRAIN_IMPORT_ERROR:
        st.error(
            f"❌  Failed to import gui_terrain_engine:\n\n`{_TERRAIN_IMPORT_ERROR}`",
            icon="❌"
        )
        st.code("pip install rasterio", language="bash")
    else:
        _terr_ok, _terr_msg = terrain_available()
        st.error(f"❌  UCMuon Terrain not ready: {_terr_msg}")
        st.code("pip install rasterio", language="bash")


# ══════════════════════════════════════════════════════════════════════════════
# TAB 6 — DENSITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
with tab_density:
    if _DENSITY_GUI_OK and render_density_analysis_tab is not None:
        try:
            render_density_analysis_tab()
        except Exception as _dens_exc:
            import traceback as _dens_tb
            st.error(f"❌  Density Analysis tab error: `{_dens_exc}`", icon="❌")
            st.code(_dens_tb.format_exc())
    else:
        st.error("❌  gui_density_analysis.py failed to import.", icon="❌")


# ── Autosave — runs AFTER all tabs so terrain widget keys are already claimed ──
save_settings()

