# gui/gui_terrain_engine.py  —  UCLouvain Muography Group
# Tab 5: UCMuon Terrain — DEM-aware transport using existing physics engines
#
# DESIGN CONCEPT
# ─────────────────────────────────────────────────────────────────────────────
# The DEM is NOT a transport engine. It is a geometry that provides
# per-direction overburden depth for a detector at known GPS coordinates.
# The actual muon transport is done by whichever engine the user selects
# (UCMuon-MC, MUSIC, Bethe-Bloch, or PROPOSAL).
#
# Workflow:
#   1. Load DEM (GeoTIFF) + set detector GPS position
#   2. Ray-trace through terrain → overburden [g/cm²] per (az, ze) bin
#   3. Load Tab 1 surface muon file (same file used in Tab 2)
#   4. For each direction bin:
#        - Filter muons in that angular bin from the surface file
#        - Compute depth = overburden / (rho × 100) [m]
#        - Transport with selected engine at that depth
#   5. Collect all survived muons → write 18-col underground file
#   6. Show polar flux maps, directional survival, muogram
#
# Public API consumed by cosmoaleph_gui.py:
#   terrain_available()         → (ok:bool, message:str)
#   render_terrain_tab(...)     → renders full Tab 5
#
# Author: Hamid Basiri <hamid.basiri@uclouvain.be>
# MIT License 2026

import sys
import os
import time
import subprocess
import importlib.util
import tempfile
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from pathlib import Path

# Synthetic geometry engine (optional — graceful fallback if missing)
try:
    from gui_geometry_engine import render_geometry_builder, SyntheticDEM
    _GEOM_ENGINE_OK = True
except ImportError:
    _GEOM_ENGINE_OK = False
    render_geometry_builder = None
    SyntheticDEM = None

# CSG / volumetric geometry engine (optional)
try:
    from gui_csg_engine import render_csg_builder, CSGGeometry
    _CSG_ENGINE_OK = True
except ImportError:
    _CSG_ENGINE_OK = False
    render_csg_builder = None
    CSGGeometry = None

# CSG per-muon transport through detector cell (optional — needs gui_csg_engine)
try:
    from gui_csg_transport import (
        transport_muons_through_csg,
        transport_muons_multi_detector,
        render_detector_cell_selector,
        plot_csg_transport_results,
        write_phits_dump_file,
    )
    _CSG_TRANSPORT_OK = True
except ImportError:
    _CSG_TRANSPORT_OK = False
    transport_muons_through_csg   = None
    transport_muons_multi_detector = None
    render_detector_cell_selector = None
    plot_csg_transport_results    = None
    write_phits_dump_file         = None

_TERRAIN_VERSION = "2.0.0"


def _sf(val, default):
    """Safely convert to float, returning default on failure."""
    try:
        return float(val)
    except (TypeError, ValueError):
        return float(default)


DARK = dict(
    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
    font=dict(color="white", size=11),
)

# ─────────────────────────────────────────────────────────────────────────────
# Availability check
# ─────────────────────────────────────────────────────────────────────────────

def terrain_available():
    """Return (ok, message). Never cached — always checks live state."""
    gui_dir  = Path(__file__).resolve().parent
    proj_dir = gui_dir.parent

    def _find(fname):
        for d in (gui_dir, proj_dir, Path(".")):
            if (d / fname).exists():
                return d / fname
        return None

    driver = _find("ucmuon_terrain_driver.py")
    bmc    = _find("ucmuon_backward_mc.py")
    missing = []
    if not driver: missing.append("ucmuon_terrain_driver.py")
    if not bmc:    missing.append("ucmuon_backward_mc.py")
    if missing:
        return False, f"Missing: {', '.join(missing)}  (checked gui/, project root, cwd)"
    try:
        import rasterio as _r
        return True, f"v{_TERRAIN_VERSION}  (rasterio {_r.__version__})"
    except ImportError:
        return False, "rasterio not installed — run:  pip install rasterio"


# ─────────────────────────────────────────────────────────────────────────────
# DEM loading and overburden computation (lazy-imports rasterio)
# ─────────────────────────────────────────────────────────────────────────────

def _load_terrain_driver(script_dir):
    path = Path(script_dir) / "ucmuon_terrain_driver.py"
    spec = importlib.util.spec_from_file_location("ucmuon_terrain_driver", str(path))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ─────────────────────────────────────────────────────────────────────────────
# Angular bin assignment (numpy-vectorised)
# ─────────────────────────────────────────────────────────────────────────────

def assign_direction_bins(surface_df, n_az, n_ze, ze_max_deg):
    """
    Assign each muon to an (az_bin, ze_bin) index pair.

    Azimuth:  geographic convention 0=N, 90=E, CW.  Derived from (cx, cy).
    Zenith:   angle from vertical.  Derived from cz.

    Returns az_idx, ze_idx arrays (int, shape N); invalid/out-of-range = -1.
    """
    df = surface_df
    # Direction cosines — surface muons go downward so cz < 0
    # theta = arccos(-cz) for downward; phi from cx, cy
    cz = df["cz"].values.astype(float)
    cx = df["cx"].values.astype(float)
    cy = df["cy"].values.astype(float)

    theta_rad = np.arccos(np.clip(-cz, -1.0, 1.0))    # zenith from vertical
    theta_deg = np.degrees(theta_rad)

    # Geographic azimuth: 0=N (−y in ENU where y=North), CW
    # In our coordinate system: x=East, y=North, z=Up
    # phi_geog = atan2(cx, cy) = atan2(East, North)
    phi_geog_deg = (np.degrees(np.arctan2(cx, cy)) + 360.0) % 360.0

    # Bin assignment
    az_step = 360.0 / n_az
    ze_step = ze_max_deg / n_ze

    az_idx = np.floor(phi_geog_deg / az_step).astype(int) % n_az
    ze_idx = np.floor(theta_deg    / ze_step).astype(int)

    # Mark out-of-range (theta > ze_max_deg) as invalid
    ze_idx[theta_deg > ze_max_deg] = -1
    ze_idx[theta_deg == 0.0]       = 0   # exactly vertical → bin 0

    return az_idx, ze_idx


# ─────────────────────────────────────────────────────────────────────────────
# Per-bin transport using UCMuon-MC (native per-muon depth)
# ─────────────────────────────────────────────────────────────────────────────

def transport_stochastic_terrain(surface_df, overburden_map, az_c, ze_c,
                                  open_sky_map, rho, n_az, n_ze, ze_max_deg,
                                  v_cut, ms_enable, script_dir,
                                  progress_container=None):
    """
    Transport all Tab-1 surface muons through terrain using UCMuon-MC.

    Each muon gets its own per-direction depth derived from the DEM overburden map.
    Open-sky muons are marked alive=1 with zero energy loss.

    Returns a DataFrame with 18-col underground format.
    """
    import pandas as pd

    gui_dir  = Path(script_dir)
    drv_path = gui_dir / "ucmuon_stochastic_driver.py"
    if not drv_path.exists():
        raise FileNotFoundError(f"ucmuon_stochastic_driver.py not found in {gui_dir}")

    spec = importlib.util.spec_from_file_location("ucmuon_stochastic_driver", str(drv_path))
    drv  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(drv)

    rng    = np.random.default_rng(42)
    mat    = drv._MAT_DB[1]   # Standard Rock — TODO: expose material choice

    az_idx, ze_idx = assign_direction_bins(surface_df, n_az, n_ze, ze_max_deg)

    n_total   = len(surface_df)
    out_rows  = []
    n_done    = 0

    for ia in range(n_az):
        for iz in range(n_ze):
            mask = (az_idx == ia) & (ze_idx == iz)
            if not mask.any():
                continue

            bin_df = surface_df[mask].copy()
            n_bin  = len(bin_df)

            ob_gcm2 = float(overburden_map[ia, iz])
            sky     = bool(open_sky_map[ia, iz])

            if sky or ob_gcm2 < 1.0:
                # Open sky: all muons survive unchanged — position stays at surface level
                for _, row in bin_df.iterrows():
                    out_rows.append({
                        "EventID": int(row["EventID"]),
                        "xs":      float(row["x"]),
                        "ys":      float(row["y"]),
                        "zs":      float(row.get("z", 0.0)),
                        "Es":      float(row["E"]),
                        "theta_s": float(row["theta"]),
                        "phi_s":   float(row["phi"]),
                        "charge":  int(row["charge"]),
                        "alive":   1,
                        "x":       float(row["x"]),
                        "y":       float(row["y"]),
                        "z":       float(row.get("z", 0.0)),
                        "E":       float(row["E"]),
                        "cx":      float(row["cx"]),
                        "cy":      float(row["cy"]),
                        "cz":      float(row["cz"]),
                        "theta":   float(row["theta"]),
                        "phi":     float(row["phi"]),
                    })
                n_done += n_bin
                continue

            # Vertical-equivalent depth from slant overburden
            ze_mid_rad = np.radians(ze_c[iz])
            cos_ze     = max(np.cos(ze_mid_rad), 0.02)
            depth_m    = ob_gcm2 / (rho * 100.0 * cos_ze)
            depth_m    = max(depth_m, 0.5)

            # Build per-muon dict matching transport() API.
            # Surface column E is TOTAL energy; transport() expects kinetic.
            _E_tot = bin_df["E"].values.astype(float)
            _E_kin = np.maximum(_E_tot - drv.M_MU_GEV, 0.0)
            mu = dict(
                EventID = bin_df["EventID"].values.astype(int),
                x       = bin_df["x"].values.astype(float),
                y       = bin_df["y"].values.astype(float),
                z       = bin_df["z"].values.astype(float),
                theta   = bin_df["theta"].values.astype(float),
                phi     = bin_df["phi"].values.astype(float),
                E_tot_GeV= _E_tot,
                Ekin_GeV= _E_kin,
                Ekin_MeV= _E_kin * 1000.0,
                charge  = bin_df["charge"].values.astype(int),
                cx      = bin_df["cx"].values.astype(float),
                cy      = bin_df["cy"].values.astype(float),
                cz      = bin_df["cz"].values.astype(float),
            )

            result = drv.transport(
                mu, depth_m, rho, mat,
                n_steps=0, v_cut=v_cut,
                ms_enable=ms_enable, rng=rng,
            )

            for k in range(n_bin):
                out_rows.append({
                    # Surface (source) columns — from input muon
                    "EventID": int(mu["EventID"][k]),
                    "xs":      float(mu["x"][k]),
                    "ys":      float(mu["y"][k]),
                    "zs":      float(mu["z"][k]),
                    "Es":      float(mu["E_tot_GeV"][k]),
                    "theta_s": float(mu["theta"][k]),
                    "phi_s":   float(mu["phi"][k]),
                    "charge":  int(mu["charge"][k]),
                    # Underground (result) columns — transport() returns *_f suffix keys
                    "alive":   int(result["alive"][k]),
                    "x":       float(result["x_f"][k]),
                    "y":       float(result["y_f"][k]),
                    "z":       float(result["z_f"][k]),
                    # 18-col convention: E = total energy for survivors, 0 stopped
                    "E":       (float(result["E_kin_f_MeV"][k]) / 1000.0 + drv.M_MU_GEV
                                if result["alive"][k] else 0.0),
                    "cx":      float(result["cx_f"][k]),
                    "cy":      float(result["cy_f"][k]),
                    "cz":      float(result["cz_f"][k]),
                    "theta":   float(result["theta_f"][k]),
                    "phi":     float(result["phi_f"][k]),
                })

            n_done += n_bin
            if progress_container:
                frac = min(n_done / max(n_total, 1), 0.99)
                progress_container.progress(
                    frac,
                    text=f"⏳  {n_done:,} / {n_total:,} muons transported  ({100*frac:.0f}%)"
                )

    return pd.DataFrame(out_rows)


# ─────────────────────────────────────────────────────────────────────────────
# Per-bin transport using subprocess engines (MUSIC / BB+MS / PROPOSAL)
# Each direction bin gets its own depth → separate subprocess call
# ─────────────────────────────────────────────────────────────────────────────

def transport_batched_terrain(surface_df, overburden_map, az_c, ze_c,
                               open_sky_map, rho, n_az, n_ze, ze_max_deg,
                               engine, engine_cfg, script_dir, project_dir,
                               build_music_input_fn, build_phitsxs_input_fn,
                               build_proposal_input_fn,
                               progress_container=None):
    """
    Transport Tab-1 surface muons per direction bin using subprocess engines.

    For each (az, ze) bin:
      1. Write a temporary input file with muons in that bin
      2. Call the chosen engine with depth = DEM overburden for that bin
      3. Read the output, append to result DataFrame

    Returns a DataFrame with 18-col underground format.
    """
    import pandas as pd

    az_idx, ze_idx = assign_direction_bins(surface_df, n_az, n_ze, ze_max_deg)
    out_dfs = []
    n_total = len(surface_df)
    n_done  = 0
    n_bins  = n_az * n_ze

    # Compiled binaries live in bin/, not the project root; MUSIC must also
    # RUN from bin/ because its table files (music-eloss-*.dat,
    # music-double-diff-*.dat) are opened relative to the CWD — same setup
    # as the Transport tab.
    bin_dir = Path(project_dir) / "bin"
    if engine == "MUSIC":
        if not ((bin_dir / "ucmuon_transport_music_omp").exists()
                or (bin_dir / "ucmuon_transport_music").exists()):
            raise FileNotFoundError(
                "MUSIC binary not found in bin/ "
                "(ucmuon_transport_music_omp) — run `make local`."
            )
    n_failed_bins = 0
    last_fail_msg = ""

    # Precompute column format — surface file is 13-col or 14-col
    ncols  = engine_cfg.get("ncols", 13)
    col13 = "EventID,x,y,z,p,px,py,pz,theta,phi,E,charge,det_mask".split(",")
    col14 = "EventID,x,y,z,p,px,py,pz,theta,phi,E,charge,hit_flag,det_mask".split(",")
    cols_s = col14 if ncols == 14 else col13

    for ia in range(n_az):
        for iz in range(n_ze):
            mask = (az_idx == ia) & (ze_idx == iz)
            if not mask.any():
                n_done += int(mask.sum())
                continue

            bin_df  = surface_df[mask].copy()
            n_bin   = len(bin_df)
            ob_gcm2 = float(overburden_map[ia, iz])
            sky     = bool(open_sky_map[ia, iz])

            if sky or ob_gcm2 < 1.0:
                # Open sky — mark all alive, copy through
                for _, row in bin_df.iterrows():
                    out_dfs.append(pd.DataFrame([{
                        "EventID": int(row["EventID"]),
                        "xs": row["x"], "ys": row["y"], "zs": row["z"],
                        "Es": row["E"], "theta_s": row["theta"], "phi_s": row["phi"],
                        "charge": int(row["charge"]), "alive": 1,
                        "x": row["x"], "y": row["y"], "z": row["z"],
                        "E": row["E"], "cx": row["cx"], "cy": row["cy"], "cz": row["cz"],
                        "theta": row["theta"], "phi": row["phi"],
                    }]))
                n_done += n_bin
                continue

            ze_mid_rad = np.radians(ze_c[iz])
            cos_ze     = max(np.cos(ze_mid_rad), 0.02)
            depth_m    = ob_gcm2 / (rho * 100.0 * cos_ze)
            depth_m    = max(depth_m, 0.5)

            # Write temp input file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".dat",
                                             delete=False, dir=str(project_dir)) as tf:
                tmp_in = tf.name
            with tempfile.NamedTemporaryFile(mode="w", suffix=".dat",
                                             delete=False, dir=str(project_dir)) as tf:
                tmp_out = tf.name

            # Write in correct column format
            _write_surface_bin(bin_df, tmp_in, ncols)

            # Build engine stdin
            cfg_bin = {**engine_cfg,
                       "infile": tmp_in, "outfile": tmp_out,
                       "depth_m": depth_m, "transport_all": True, "ncols": ncols}

            run_cwd = project_dir
            if engine == "MUSIC":
                stdin_str = build_music_input_fn(cfg_bin)
                if (bin_dir / "ucmuon_transport_music_omp").exists():
                    cmd = [str(bin_dir / "ucmuon_transport_music_omp")]
                else:
                    cmd = [str(bin_dir / "ucmuon_transport_music")]
                env_run = {**os.environ, "OMP_NUM_THREADS": "1"}
                run_cwd = bin_dir   # MUSIC table files are opened relative to CWD

            elif engine == "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS":
                stdin_str = build_phitsxs_input_fn(cfg_bin)
                _bb_py  = script_dir / "ucmuon_bb_driver.py"
                _bb_bin = bin_dir / "ucmuon_transport_bb_omp"
                if _bb_py.exists():
                    cmd     = [sys.executable, str(_bb_py)]
                    env_run = {**os.environ}
                else:
                    cmd     = [str(_bb_bin)]
                    env_run = {**os.environ, "OMP_NUM_THREADS": "1"}

            elif engine == "PROPOSAL":
                stdin_str = build_proposal_input_fn(cfg_bin)
                cmd       = [sys.executable, str(script_dir / "proposal_driver.py")]
                env_run   = {**os.environ, "PROPOSAL_LOG_LEVEL": "err",
                             "SPDLOG_LEVEL": "err"}

            else:
                # UCMuon-MC — should not reach here (handled separately)
                Path(tmp_in).unlink(missing_ok=True)
                Path(tmp_out).unlink(missing_ok=True)
                continue

            # Run subprocess
            proc_res = None
            try:
                proc_res = subprocess.run(
                    cmd, input=stdin_str, capture_output=True, text=True,
                    cwd=str(run_cwd), env=env_run, timeout=300,
                )
            except Exception as _sub_exc:
                last_fail_msg = f"{cmd[0]}: {_sub_exc}"

            # Read output
            _bin_ok = False
            if Path(tmp_out).exists() and Path(tmp_out).stat().st_size > 0:
                try:
                    import pandas as _pd
                    _df_bin = _pd.read_csv(tmp_out, sep=r"\s+", comment="#", header=None)
                    if _df_bin.shape[1] == 18:
                        _df_bin.columns = (
                            "EventID,xs,ys,zs,Es,theta_s,phi_s,charge,alive,"
                            "x,y,z,E,cx,cy,cz,theta,phi"
                        ).split(",")
                        out_dfs.append(_df_bin)
                        _bin_ok = True
                    else:
                        last_fail_msg = (f"unexpected column count "
                                         f"{_df_bin.shape[1]} in engine output")
                except Exception as _read_exc:
                    last_fail_msg = f"could not read engine output: {_read_exc}"
            elif proc_res is not None:
                _tail = ((proc_res.stderr or "") + (proc_res.stdout or ""))[-400:]
                last_fail_msg = (f"engine wrote no output "
                                 f"(exit {proc_res.returncode}): …{_tail}")
            if not _bin_ok:
                n_failed_bins += 1

            Path(tmp_in).unlink(missing_ok=True)
            Path(tmp_out).unlink(missing_ok=True)

            n_done += n_bin
            if progress_container:
                frac = min(n_done / max(n_total, 1), 0.99)
                progress_container.progress(
                    frac,
                    text=f"⏳  {n_done:,} / {n_total:,} muons  |  bin {ia*n_ze+iz+1}/{n_bins}  ({100*frac:.0f}%)"
                )

    if n_failed_bins and not out_dfs:
        raise RuntimeError(
            f"{engine} transport failed for all {n_failed_bins} occupied "
            f"direction bin(s). Last error: {last_fail_msg}"
        )
    if n_failed_bins:
        st.warning(
            f"⚠️  {engine}: {n_failed_bins} direction bin(s) produced no "
            f"transport output and were dropped. Last error: {last_fail_msg}"
        )

    if out_dfs:
        import pandas as pd
        return pd.concat(out_dfs, ignore_index=True)
    import pandas as pd
    return pd.DataFrame()


def _write_surface_bin(df, path, ncols):
    """Write a subset of surface muons to a temp file in 13-col or 14-col format."""
    import numpy as np
    with open(path, "w") as fh:
        for _, row in df.iterrows():
            E   = float(row["E"])
            p   = float(np.sqrt(max(E**2 - 0.105658**2, 0.0)))
            cx  = float(row["cx"]); cy = float(row["cy"]); cz = float(row["cz"])
            th  = float(row["theta"]); ph = float(row["phi"])
            if ncols == 14:
                fh.write(
                    f"{int(row['EventID']):10d}  {float(row['x']):13.4f}  {float(row['y']):13.4f}"
                    f"  {float(row['z']):13.4f}  {p:13.6f}  {p*cx:13.6f}  {p*cy:13.6f}"
                    f"  {p*cz:13.6f}  {th:13.9f}  {ph:13.9f}  {E:13.6f}"
                    f"  {int(row['charge']):4d}  1  0\n"
                )
            else:
                fh.write(
                    f"{int(row['EventID']):10d}  {float(row['x']):13.4f}  {float(row['y']):13.4f}"
                    f"  {float(row['z']):13.4f}  {p:13.6f}  {p*cx:13.6f}  {p*cy:13.6f}"
                    f"  {p*cz:13.6f}  {th:13.9f}  {ph:13.9f}  {E:13.6f}"
                    f"  {int(row['charge']):4d}  0\n"
                )


# ─────────────────────────────────────────────────────────────────────────────
# Results visualisation
# ─────────────────────────────────────────────────────────────────────────────

def _polar_heatmap(az_c, ze_c, data, title, unit, colorscale="Jet", mask=None,
                   log_scale=False):
    """
    Polar (az × ze) heatmap using scatterpolar.
    log_scale=True: display log10(data+1) for overburden maps with wide dynamic range.
    """
    theta_plot = (90.0 - az_c) % 360.0
    fig = go.Figure()

    # Prepare display values
    plot_data = data.copy().astype(float)
    if mask is not None:
        plot_data[mask] = np.nan

    if log_scale:
        # log10(x+1) so that zeros map to 0 cleanly
        valid = ~np.isnan(plot_data) & (plot_data > 0)
        log_data = np.full_like(plot_data, np.nan)
        log_data[valid] = np.log10(plot_data[valid])
        display_data = log_data
        disp_unit = f"log₁₀({unit})"
        # Build custom colorbar ticktext showing actual values
        vmax_log = float(np.nanmax(log_data)) if np.any(~np.isnan(log_data)) else 1.0
        cbar_ticks = np.linspace(0, vmax_log, 6)
        cbar_text  = [f"{10**v:,.0f}" for v in cbar_ticks]
    else:
        display_data = plot_data
        disp_unit = unit
        cbar_ticks = None
        cbar_text  = None

    d_min = float(np.nanmin(display_data)) if np.any(~np.isnan(display_data)) else 0.0
    d_max = float(np.nanmax(display_data)) if np.any(~np.isnan(display_data)) else 1.0
    if d_max <= d_min:
        d_max = d_min + 1.0

    colorbar_cfg = dict(
        title=dict(text=disp_unit, font=dict(color="white", size=11)),
        tickfont=dict(color="white"),
    )
    if cbar_ticks is not None:
        colorbar_cfg["tickvals"] = cbar_ticks.tolist()
        colorbar_cfg["ticktext"] = cbar_text

    for iz, ze in enumerate(ze_c):
        vals = display_data[:, iz].copy()
        # hover shows original (non-log) value
        orig_vals = data[:, iz].copy()
        if mask is not None:
            orig_vals[mask[:, iz]] = np.nan

        hover_unit = unit.replace("log₁₀(", "").rstrip(")")
        fig.add_trace(go.Scatterpolar(
            r      = np.full(len(az_c) + 1, ze),
            theta  = np.append(theta_plot, theta_plot[0]),
            mode   = "markers",
            marker = dict(
                size       = 11,
                color      = np.append(vals, vals[0]),
                colorscale = colorscale,
                showscale  = (iz == len(ze_c) - 1),
                colorbar   = colorbar_cfg,
                cmin       = d_min,
                cmax       = d_max,
            ),
            showlegend = False,
            customdata = np.append(orig_vals, orig_vals[0]).reshape(-1, 1),
            hovertemplate = (
                f"Az=%{{theta:.0f}}°  Ze={ze:.0f}°<br>"
                f"{hover_unit}=%{{customdata[0]:,.0f}}<extra></extra>"
            ),
        ))

    fig.update_layout(
        **DARK, height=480,
        title=dict(text=title, font=dict(size=12)),
        polar=dict(
            bgcolor="rgb(20,22,30)",
            radialaxis=dict(
                title="Zenith [°]", range=[0, float(ze_c[-1]) + 5],
                tickfont=dict(color="white"), gridcolor="#2a2a3a"
            ),
            angularaxis=dict(
                direction="clockwise", rotation=90,
                tickfont=dict(color="white"), gridcolor="#2a2a3a",
                tickvals=[0, 90, 180, 270], ticktext=["N", "E", "S", "W"],
            ),
        ),
        margin=dict(l=60, r=80, t=55, b=40),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# DEM download helper (OpenTopography — free, no account for SRTM)
# ─────────────────────────────────────────────────────────────────────────────

# Bundled sample DEM: Mt. Vesuvius, SRTM GL1 30 m (public domain, NASA/USGS).
# Matches the examples/vesuvius MURAVES example — see examples/vesuvius/DEM_SOURCE.md
_BUNDLED_VESUVIUS_DEM = (Path(__file__).resolve().parent.parent
                         / "examples" / "vesuvius" / "vesuvius_dem.tif")


def _download_dem(south, north, west, east, product, outpath, api_key="demoapikeyot2022"):
    try:
        import requests
    except ImportError:
        return False, "requests not installed: pip install requests"
    params = dict(demtype=product, south=south, north=north, west=west, east=east,
                  outputFormat="GTiff", API_Key=api_key)
    try:
        r = requests.get("https://portal.opentopography.org/API/globaldem",
                         params=params, timeout=120, stream=True)
        if r.status_code == 200:
            with open(outpath, "wb") as fh:
                for chunk in r.iter_content(chunk_size=65536): fh.write(chunk)
            return True, f"Downloaded {Path(outpath).stat().st_size/1e6:.1f} MB → {outpath}"
        return False, f"API error {r.status_code}: {r.text[:200]}"
    except Exception as e:
        return False, f"Download failed: {e}"


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PANEL
# ─────────────────────────────────────────────────────────────────────────────

_SPEC_LABELS = {
    1: "CosmoALEPH",
    2: "Power-law  (E⁻³·⁷)",
    3: "Guan et al. (2015)",
    4: "Frosin et al. (2025)",
}

_ENGINE_OPTIONS = [
    "★ UCMuon-MC (Python)  ← recommended for terrain",
    "MUSIC",
    "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS",
    "PROPOSAL",
]
# Internal token "UCMuon Stochastic" kept stable — used in run-logic comparisons.
_ENGINE_INTERNAL = {
    "★ UCMuon-MC (Python)  ← recommended for terrain": "UCMuon Stochastic",
    "MUSIC": "MUSIC",
    "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS": "Bethe-Bloch",
    "PROPOSAL": "PROPOSAL",
}



# ─────────────────────────────────────────────────────────────────────────────
# Rectangular muogram heatmap  (azimuth × elevation)
# ─────────────────────────────────────────────────────────────────────────────
def _muogram_heatmap(az_c, ze_c, data_2d, title, unit,
                     colorscale="Jet", zmin=None, zmax=None,
                     log_scale=False, mask=None):
    """
    2D rectangular heatmap:  X = azimuth [deg],  Y = elevation = 90 - zenith [deg]
    This is the standard muography display format (Carloganu 2013).
    data_2d shape: (n_az, n_ze)
    mask shape:    (n_az, n_ze) — True where data should be hidden (open sky)
    """
    import plotly.graph_objects as _go

    el_c = 90.0 - ze_c          # elevation centres
    az_centered = (az_c + 180.0) % 360.0 - 180.0  # shift so N=0 is centre

    # Sort by azimuth for display
    sort_idx = np.argsort(az_centered)
    az_s     = az_centered[sort_idx]
    data_s   = data_2d[sort_idx, :]
    mask_s   = mask[sort_idx, :] if mask is not None else None

    # Build Z array — NaN where masked (open sky)
    Z = data_s.astype(float).copy()
    if mask_s is not None:
        Z[mask_s] = np.nan

    if log_scale:
        with np.errstate(divide='ignore', invalid='ignore'):
            Z = np.where((Z > 0) & ~np.isnan(Z), np.log10(Z), np.nan)
        unit_label = f"log₁₀({unit})"
    else:
        unit_label = unit

    # Elevation on Y (ascending), azimuth on X
    # plotly heatmap: x=az, y=el, z shape=(n_el, n_az) i.e. (n_ze, n_az)
    z_plot = Z.T          # (n_ze, n_az)

    v_min = zmin if zmin is not None else float(np.nanmin(z_plot)) if np.any(~np.isnan(z_plot)) else 0
    v_max = zmax if zmax is not None else float(np.nanmax(z_plot)) if np.any(~np.isnan(z_plot)) else 1
    if v_max <= v_min: v_max = v_min + 1.0

    hover = (
        "Az: %{x:.1f}°  El: %{y:.1f}°<br>"
        + unit + ": %{z:.3g}<extra></extra>"
    )

    fig = _go.Figure(_go.Heatmap(
        x=az_s, y=el_c, z=z_plot,
        colorscale=colorscale,
        zmin=v_min, zmax=v_max,
        colorbar=dict(
            title=dict(text=unit_label, font=dict(color="white")),
            tickfont=dict(color="white"),
        ),
        hovertemplate=hover,
        hoverongaps=False,
        xgap=1, ygap=1,
    ))
    fig.update_layout(
        **DARK, height=380,
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(
            title="Azimuth [°]  (N=0°, E=90°, W=−90°)",
            gridcolor="#2a2a3a", zeroline=True,
            zerolinecolor="#ffd700", zerolinewidth=1.5,
        ),
        yaxis=dict(
            title="Elevation above horizon [°]",
            gridcolor="#2a2a3a",
        ),
        margin=dict(l=70, r=20, t=55, b=60),
    )
    # Add N/S/E/W annotations
    for az_ann, label in [(0,"N"), (90,"E"), (-90,"W"), (180,"S"), (-180,"S")]:
        if az_s.min() <= az_ann <= az_s.max():
            fig.add_vline(x=az_ann,
                line=dict(color="rgba(255,215,0,0.3)", width=1, dash="dot"),
                annotation_text=label,
                annotation_font=dict(color="rgba(255,215,0,0.6)", size=10))
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2D top-down terrain map with blocked direction footprints
# ─────────────────────────────────────────────────────────────────────────────
def _dem_2d_topdown(dem_path, det_lat, det_lon, det_alt,
                    az_c, ze_c, ob_map, sky_map, rho,
                    radius_deg=0.15,
                    elev_arr=None, tfm_arr=None):
    """
    2D top-down (plan-view) map of the terrain.

    Shows:
      • DEM elevation as a filled contour background (Earth colorscale)
      • Detector at the ENU origin (gold star)
      • Summit marker (red triangle)
      • For each blocked direction bin: a colored line radiating from the
        detector in the correct geographic azimuth, length = horizontal
        projection of the slant path (= slant_m × sin(ze)).  The endpoint
        is approximately where the ray exits the terrain.
      • Endpoints are colored by log₁₀(overburden) on a Jet scale.
      • A fake colorbar (via a hidden Scatter) reproduces the scale.

    ENU convention: X = East [m], Y = North [m], detector at (0, 0).
    """
    import plotly.graph_objects as _go
    try:
        import rasterio as _rio
    except ImportError:
        return None, "rasterio not installed"

    try:
        # Accept pre-loaded (elev_arr, tfm_arr) from synthetic DEM, or load from file.
        if elev_arr is not None and tfm_arr is not None:
            elev_full = elev_arr.astype(float)
            tfm = tfm_arr
        else:
            # Use _load_terrain_driver to call load_dem() — this converts XYZ/ASC
            # to WGS84 so clip arithmetic works with degree coordinates.
            _drv2 = _load_terrain_driver(Path(__file__).resolve().parent)
            elev_full, tfm = _drv2.load_dem(dem_path)
            elev_full = elev_full.astype(float)
        h, w = elev_full.shape

        # ── clip DEM to radius_deg around detector ────────────────────────
        lon_min_c = det_lon - radius_deg
        lon_max_c = det_lon + radius_deg
        lat_min_c = det_lat - radius_deg
        lat_max_c = det_lat + radius_deg
        col_min = max(0, int((lon_min_c - tfm.c) / tfm.a))
        col_max = min(w, int((lon_max_c - tfm.c) / tfm.a) + 1)
        row_min = max(0, int((lat_max_c - tfm.f) / tfm.e))
        row_max = min(h, int((lat_min_c - tfm.f) / tfm.e) + 1)
        elev = elev_full[row_min:row_max, col_min:col_max]
        n_rows, n_cols = elev.shape

        lons = tfm.c + (np.arange(col_min, col_max) + 0.5) * tfm.a
        lats = tfm.f + (np.arange(row_min, row_max) + 0.5) * tfm.e
        LON, LAT = np.meshgrid(lons, lats)
        x_m = (LON - det_lon) * 111320.0 * np.cos(np.radians(det_lat))
        y_m = (LAT - det_lat) * 111320.0

        # Downsample for browser performance
        step = max(1, max(n_rows, n_cols) // 250)
        x_d  = x_m[::step, ::step]
        y_d  = y_m[::step, ::step]
        z_d  = elev[::step, ::step]

        # ── helpers ───────────────────────────────────────────────────────
        def _jet(t):
            """Jet colorscale: 0=blue → 0.25=cyan → 0.5=green → 0.75=yellow → 1=red."""
            t = float(np.clip(t, 0, 1))
            if   t < 0.25: r, g, b = 0,          int(255*t*4),     255
            elif t < 0.50: r, g, b = 0,           255,              int(255*(1-(t-0.25)*4))
            elif t < 0.75: r, g, b = int(255*(t-0.5)*4), 255,       0
            else:          r, g, b = 255,          int(255*(1-(t-0.75)*4)), 0
            return f"rgb({r},{g},{b})"

        # ── figure ────────────────────────────────────────────────────────
        fig = _go.Figure()

        # DEM contour fill background
        # NOTE: 'line' is a top-level Contour argument, NOT inside contours=dict()
        fig.add_trace(_go.Contour(
            x=x_d[0, :], y=y_d[:, 0], z=z_d,
            colorscale="Earth",
            contours=dict(coloring="fill", showlines=True),
            line=dict(color="rgba(255,255,255,0.12)", width=0.5),
            colorbar=dict(
                title=dict(text="Elevation [m]", font=dict(color="white", size=11)),
                tickfont=dict(color="white", size=10),
                x=1.02, thickness=14, len=0.55,
                bgcolor="rgba(20,22,30,0.85)",
                bordercolor="#555", borderwidth=1,
            ),
            hovertemplate="E: %{x:.0f} m  N: %{y:.0f} m  Elev: %{z:.0f} m<extra></extra>",
        ))

        # ── blocked direction lines ───────────────────────────────────────
        if ob_map is not None and sky_map is not None:
            ob_blocked = ob_map[~sky_map]
            if ob_blocked.size > 0:
                ob_min_log = float(np.log10(max(ob_blocked.min(), 1.0)))
                ob_max_log = float(np.log10(max(ob_blocked.max(), 1.0)))

                for ia, az in enumerate(az_c):
                    for iz, ze in enumerate(ze_c):
                        if sky_map[ia, iz] or ob_map[ia, iz] < 10:
                            continue
                        ob_val  = float(ob_map[ia, iz])
                        slant_m = ob_val / (rho * 100.0)
                        az_r    = np.radians(float(az))
                        ze_r    = np.radians(float(ze))
                        horiz_m = slant_m * np.sin(ze_r)
                        ex      = horiz_m * np.sin(az_r)   # East
                        ny      = horiz_m * np.cos(az_r)   # North
                        t       = (np.log10(ob_val) - ob_min_log) / max(ob_max_log - ob_min_log, 0.1)
                        col     = _jet(t)

                        # Line: detector → terrain endpoint
                        fig.add_trace(_go.Scatter(
                            x=[0, ex], y=[0, ny],
                            mode="lines",
                            line=dict(color=col, width=2),
                            hoverinfo="skip",
                            showlegend=False,
                        ))
                        # Endpoint dot with hover info
                        fig.add_trace(_go.Scatter(
                            x=[ex], y=[ny],
                            mode="markers",
                            marker=dict(size=11, color=col,
                                        line=dict(color="white", width=1.5)),
                            hovertemplate=(
                                f"az = {az:.0f}°  ze = {ze:.0f}°<br>"
                                f"overburden = {ob_val:,.0f} g/cm²<br>"
                                f"slant path = {slant_m:.0f} m<br>"
                                f"horiz dist = {horiz_m/1000:.2f} km<extra></extra>"
                            ),
                            showlegend=False,
                        ))

                # Fake colorbar for overburden scale
                _cb_vals = np.linspace(ob_min_log, ob_max_log, 120)
                fig.add_trace(_go.Scatter(
                    x=[None], y=[None], mode="markers",
                    marker=dict(
                        colorscale="Jet",
                        cmin=ob_min_log, cmax=ob_max_log,
                        color=[0],
                        colorbar=dict(
                            title=dict(text="log₁₀(overburden g/cm²)",
                                       font=dict(color="white", size=11)),
                            tickfont=dict(color="white", size=10),
                            tickvals=np.linspace(ob_min_log, ob_max_log, 5).tolist(),
                            ticktext=[f"{10**v:,.0f}" for v in
                                      np.linspace(ob_min_log, ob_max_log, 5)],
                            x=1.16, thickness=14, len=0.55,
                            bgcolor="rgba(20,22,30,0.85)",
                            bordercolor="#555", borderwidth=1,
                        ),
                        showscale=True,
                    ),
                    hoverinfo="skip", showlegend=False,
                ))

        # ── detector (gold star at origin) ────────────────────────────────
        fig.add_trace(_go.Scatter(
            x=[0], y=[0],
            mode="markers+text",
            marker=dict(size=16, color="#ffd700", symbol="star",
                        line=dict(color="white", width=2)),
            text=[f"Detector<br>{det_alt:.0f} m"],
            textposition="top center",
            textfont=dict(color="#ffd700", size=10),
            name=f"Detector ({det_alt:.0f} m a.s.l.)",
            hovertemplate=(
                f"Detector<br>{det_lat:.5f}°N  {det_lon:.5f}°E<br>"
                f"Alt = {det_alt:.0f} m<extra></extra>"
            ),
        ))

        # ── summit marker (DEM maximum in clip region) ────────────────────
        imax = np.unravel_index(np.nanargmax(z_d), z_d.shape)
        sx, sy, sz = float(x_d[imax]), float(y_d[imax]), float(z_d[imax])
        summit_dist_km = np.sqrt(sx**2 + sy**2) / 1000.0
        summit_az_deg  = (np.degrees(np.arctan2(sx, sy)) + 360) % 360
        fig.add_trace(_go.Scatter(
            x=[sx], y=[sy],
            mode="markers+text",
            marker=dict(size=13, color="#ff6b6b", symbol="triangle-up",
                        line=dict(color="white", width=1.5)),
            text=[f"Summit<br>{sz:.0f} m"],
            textposition="top right",
            textfont=dict(color="#ff6b6b", size=10),
            name=f"Summit ({sz:.0f} m a.s.l.)",
            hovertemplate=(
                f"Summit<br>E={sx:.0f} m  N={sy:.0f} m<br>"
                f"Alt={sz:.0f} m  dist={summit_dist_km:.2f} km  "
                f"az={summit_az_deg:.0f}°<extra></extra>"
            ),
        ))

        # Dashed line detector → summit
        fig.add_trace(_go.Scatter(
            x=[0, sx], y=[0, sy],
            mode="lines",
            line=dict(color="rgba(255,107,107,0.4)", width=1.5, dash="dash"),
            hoverinfo="skip", showlegend=False,
        ))

        # Compass cardinal labels at map edge
        r_ann = radius_deg * 111320.0 * np.cos(np.radians(det_lat)) * 0.90
        r_ann_n = radius_deg * 111320.0 * 0.90
        for az_ann, lbl, xa, ya in [
            (0,   "N", 0,      r_ann_n),
            (90,  "E", r_ann,  0      ),
            (180, "S", 0,     -r_ann_n),
            (270, "W", -r_ann, 0      ),
        ]:
            fig.add_annotation(
                x=xa, y=ya, text=f"<b>{lbl}</b>",
                font=dict(color="rgba(255,255,255,0.55)", size=14),
                showarrow=False,
            )

        fig.update_layout(
            **DARK, height=560,
            title=dict(
                text=(
                    "Top-down terrain map — blocked muon directions "
                    "(lines = ray directions, dots = terrain endpoint, colour = overburden)"
                ),
                font=dict(color="white", size=12),
            ),
            xaxis=dict(
                title="East [m]", gridcolor="#2a2a3a",
                scaleanchor="y", scaleratio=1,
                zeroline=True, zerolinecolor="rgba(255,215,0,0.25)", zerolinewidth=1,
            ),
            yaxis=dict(
                title="North [m]", gridcolor="#2a2a3a",
                zeroline=True, zerolinecolor="rgba(255,215,0,0.25)", zerolinewidth=1,
            ),
            legend=dict(
                font=dict(color="white", size=10),
                bgcolor="rgba(15,17,23,0.85)",
                bordercolor="#555", borderwidth=1,
                x=0.01, y=0.99, xanchor="left", yanchor="top",
            ),
            margin=dict(l=70, r=160, t=60, b=60),
        )
        return fig, None

    except Exception as e:
        import traceback
        return None, traceback.format_exc()


# ─────────────────────────────────────────────────────────────────────────────
# 3D DEM surface plot with detector marker
# ─────────────────────────────────────────────────────────────────────────────
def _dem_3d_plot(dem_path, det_lat, det_lon, det_alt, script_dir,
                 radius_deg=0.08,
                 elev_arr=None, tfm_arr=None):
    """
    Interactive 3D surface of the DEM centred on the detector.
    Clips the DEM to ±radius_deg around the detector for performance.
    """
    import plotly.graph_objects as _go
    try:
        import rasterio as _rio
        import rasterio.transform as _rt
    except ImportError:
        return None

    try:
        # Accept pre-loaded (elev_arr, tfm_arr) from synthetic DEM, or load from file.
        if elev_arr is not None and tfm_arr is not None:
            elev_full = elev_arr.astype(float)
            tfm = tfm_arr
        else:
            # Use load_dem() so XYZ/ASC are in WGS84 before clip arithmetic.
            _drv3 = _load_terrain_driver(Path(__file__).resolve().parent)
            elev_full, tfm = _drv3.load_dem(dem_path)
            elev_full = elev_full.astype(float)
        h, w = elev_full.shape

        # Pixel coordinates for bounding box
        lon_min_clip = det_lon - radius_deg
        lon_max_clip = det_lon + radius_deg
        lat_min_clip = det_lat - radius_deg
        lat_max_clip = det_lat + radius_deg

        # Row/col bounds — clamp and validate
        col_min = max(0, int((lon_min_clip - tfm.c) / tfm.a))
        col_max = min(w, int((lon_max_clip - tfm.c) / tfm.a) + 1)
        # tfm.e is negative for north-up DEMs (latitude decreases with row index)
        if tfm.e < 0:
            row_min = max(0, int((lat_max_clip - tfm.f) / tfm.e))
            row_max = min(h, int((lat_min_clip - tfm.f) / tfm.e) + 1)
        else:
            row_min = max(0, int((lat_min_clip - tfm.f) / tfm.e))
            row_max = min(h, int((lat_max_clip - tfm.f) / tfm.e) + 1)

        # Guard against empty or invalid slice (detector outside DEM or radius too small)
        if col_min >= col_max or row_min >= row_max:
            return ("ERROR",
                    f"DEM clip is empty for this view radius ({radius_deg}°).  "
                    f"Try increasing the View radius slider, or check that the detector "
                    f"coordinates ({det_lat:.4f}°N, {det_lon:.4f}°E) are inside the DEM bounds.",
                    "")

        elev = elev_full[row_min:row_max, col_min:col_max]
        n_rows, n_cols = elev.shape

        # Geographic coordinates of each pixel
        lons = tfm.c + (np.arange(col_min, col_max) + 0.5) * tfm.a
        lats = tfm.f + (np.arange(row_min, row_max) + 0.5) * tfm.e

        LON, LAT = np.meshgrid(lons, lats)

        # Convert to local ENU (metres) centred on detector
        x_m = (LON - det_lon) * 111320.0 * np.cos(np.radians(det_lat))
        y_m = (LAT - det_lat) * 111320.0

        # Downsample for browser performance (max ~150×150)
        step = max(1, max(n_rows, n_cols) // 150)
        x_d = x_m[::step, ::step]
        y_d = y_m[::step, ::step]
        z_d = elev[::step, ::step]

        fig = _go.Figure(_go.Surface(
            x=x_d, y=y_d, z=z_d,
            colorscale="Earth",
            colorbar=None,
            showscale=False,
            lighting=dict(ambient=0.7, diffuse=0.8, specular=0.2),
            contours=dict(z=dict(show=True, usecolormap=True, highlightcolor="white",
                                 project_z=False, width=1)),
            hovertemplate="E: %{x:.0f} m  N: %{y:.0f} m  Alt: %{z:.0f} m<extra></extra>",
        ))

        # Detector marker — text in legend only to avoid label clipping
        det_ground = float(np.nanmin(z_d)) if not np.all(np.isnan(z_d)) else det_alt - 200
        fig.add_trace(_go.Scatter3d(
            x=[0, 0], y=[0, 0], z=[det_ground, det_alt],
            mode="lines",
            line=dict(color="#ffd700", width=2, dash="dot"),
            showlegend=False, hoverinfo="skip",
        ))
        fig.add_trace(_go.Scatter3d(
            x=[0], y=[0], z=[det_alt],
            mode="markers",
            marker=dict(size=12, color="#ffd700", symbol="diamond",
                        line=dict(color="white", width=2)),
            name=f"Detector ({det_alt:.0f} m a.s.l.)",
            hovertemplate=f"Detector<br>Alt={det_alt:.0f} m<extra></extra>",
        ))

        # Summit marker (DEM max in the clipped region)
        imax = np.unravel_index(np.nanargmax(z_d), z_d.shape)
        fig.add_trace(_go.Scatter3d(
            x=[float(x_d[imax])], y=[float(y_d[imax])], z=[float(z_d[imax])],
            mode="markers",
            marker=dict(size=10, color="#ff6b6b", symbol="circle",
                        line=dict(color="white", width=1)),
            name=f"Summit ({z_d[imax]:.0f} m a.s.l.)",
            hovertemplate=f"Summit<br>Alt={z_d[imax]:.0f} m<extra></extra>",
        ))

        # Line from detector to summit
        fig.add_trace(_go.Scatter3d(
            x=[0, float(x_d[imax])],
            y=[0, float(y_d[imax])],
            z=[det_alt, float(z_d[imax])],
            mode="lines",
            line=dict(color="#ffd700", width=3, dash="dash"),
            showlegend=False,
        ))

        fig.update_layout(
            **DARK, height=520,
            title=dict(text="DEM — 3D terrain view (ENU coordinates, detector at origin)",
                       font=dict(color="white", size=12)),
            scene=dict(
                xaxis=dict(title=dict(text="East [m]", font=dict(color="white", size=12)),
                           backgroundcolor="rgb(20,22,30)",
                           gridcolor="#2a2a3a", color="white"),
                yaxis=dict(title=dict(text="North [m]", font=dict(color="white", size=12)),
                           backgroundcolor="rgb(20,22,30)",
                           gridcolor="#2a2a3a", color="white"),
                zaxis=dict(title=dict(text="Elevation [m]", font=dict(color="white", size=12)),
                           backgroundcolor="rgb(20,22,30)",
                           gridcolor="#2a2a3a", color="white"),
                bgcolor="rgb(15,17,23)",
                camera=dict(eye=dict(x=0.5, y=-2.0, z=0.8)),
                aspectmode="manual",
                aspectratio=dict(x=1, y=1, z=0.4),
            ),
            legend=dict(
                font=dict(color="white", size=11),
                bgcolor="rgba(15,17,23,0.9)",
                bordercolor="#555", borderwidth=1,
                x=0.01, y=0.99, xanchor="left", yanchor="top",
            ),
            margin=dict(l=0, r=120, t=50, b=0),
        )
        return fig
    except Exception as e:
        import traceback
        return ("ERROR", str(e), traceback.format_exc())



# ─────────────────────────────────────────────────────────────────────────────
# Continuous skymap polar plot  (like MURAVES Fig. 1 / CCS monitoring papers)
# Fills the polar disk with a smooth interpolated heatmap
# ─────────────────────────────────────────────────────────────────────────────
def _skymap_polar(az_c, ze_c, data_2d, title, unit,
                  colorscale="Inferno", mask=None,
                  log_scale=False, vmin=None, vmax=None,
                  mark_az=None, mark_ze=None, mark_label=None):
    """
    Continuous filled polar skymap.
    az_c: azimuth centres [deg], geographic (N=0, E=90, clockwise)
    ze_c: zenith centres [deg]
    data_2d: (n_az, n_ze) array
    mask: (n_az, n_ze) bool — True = hide (gray)
    mark_az, mark_ze: optional marker direction (star symbol)
    """
    import plotly.graph_objects as _go

    n_az, n_ze = len(az_c), len(ze_c)

    Z = data_2d.astype(float).copy()
    if mask is not None:
        Z[mask] = np.nan

    if log_scale:
        with np.errstate(divide="ignore", invalid="ignore"):
            Z = np.where((Z > 0) & ~np.isnan(Z), np.log10(Z), np.nan)

    zmin = vmin if vmin is not None else float(np.nanmin(Z)) if np.any(~np.isnan(Z)) else 0
    zmax = vmax if vmax is not None else float(np.nanmax(Z)) if np.any(~np.isnan(Z)) else 1
    if zmax <= zmin:
        zmax = zmin + 1

    # Build polar Barpolar traces — one per zenith ring
    # Inner radius = ze[iz-1], outer radius = ze[iz]
    # Width = az_step, Base = inner_ze radius
    az_step = 360.0 / n_az
    ze_step = ze_c[-1] / n_ze if n_ze > 0 else 5.0

    fig = _go.Figure()

    # Background: dark filled circle
    fig.add_trace(_go.Barpolar(
        r=np.full(n_az, ze_c[-1] + ze_step/2),
        theta=az_c,
        width=np.full(n_az, az_step),
        base=0,
        marker_color="rgb(15,17,23)",
        showlegend=False,
        hoverinfo="skip",
    ))

    # One Barpolar trace per zenith ring, coloured by data
    norm = plt_colors_norm(zmin, zmax, Z, colorscale)

    for iz in range(n_ze):
        r_inner = 0 if iz == 0 else ze_c[iz-1] + ze_step/2
        r_outer = ze_c[iz] + ze_step/2
        colors = []
        hovers = []
        for ia in range(n_az):
            val = float(data_2d[ia, iz])
            z_val = float(Z[ia, iz]) if not np.isnan(Z[ia, iz]) else np.nan
            masked = (mask is not None and mask[ia, iz]) or np.isnan(z_val)
            if masked:
                colors.append("rgba(40,40,60,0.6)")
            else:
                colors.append(_interp_color(z_val, zmin, zmax, colorscale))
            u = unit if not log_scale else f"log10({unit})"
            hovers.append(
                f"az={az_c[ia]:.1f}°  ze={ze_c[iz]:.1f}°<br>"
                f"{u}={val:.3g}"
            )
        fig.add_trace(_go.Barpolar(
            r=np.full(n_az, r_outer - r_inner),
            theta=az_c,
            width=np.full(n_az, az_step),
            base=np.full(n_az, r_inner),
            marker_color=colors,
            hovertext=hovers,
            hoverinfo="text",
            showlegend=False,
        ))

    # Marker (star) for special direction
    if mark_az is not None and mark_ze is not None:
        fig.add_trace(_go.Scatterpolar(
            r=[mark_ze], theta=[mark_az],
            mode="markers+text",
            marker=dict(symbol="star", size=14, color="#00ffff",
                        line=dict(color="white", width=1)),
            text=[mark_label or f"az={mark_az:.0f}° ze={mark_ze:.0f}°"],
            textposition="top center",
            textfont=dict(color="#00ffff", size=10),
            name=mark_label or "Target direction",
        ))

    # Fake colorbar via scatter
    _cb_ze = np.linspace(zmin, zmax, 100)
    _cb_col = [_interp_color(v, zmin, zmax, colorscale) for v in _cb_ze]
    fig.add_trace(_go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(
            colorscale=colorscale,
            cmin=zmin, cmax=zmax,
            color=[0],
            colorbar=dict(
                title=dict(
                    text=("log₁₀(" + unit + ")" if log_scale else unit),
                    font=dict(color="white", size=11)
                ),
                tickfont=dict(color="white", size=10),
                x=1.05, thickness=14, len=0.7,
                bgcolor="rgba(20,22,30,0.8)",
                bordercolor="#444", borderwidth=1,
            ),
            showscale=True,
        ),
        hoverinfo="skip", showlegend=False,
    ))

    # Polar layout
    fig.update_layout(
        **DARK, height=480,
        title=dict(text=title, font=dict(color="white", size=13)),
        polar=dict(
            bgcolor="rgb(15,17,23)",
            angularaxis=dict(
                tickmode="array",
                tickvals=[0, 45, 90, 135, 180, 225, 270, 315],
                ticktext=["N", "45°", "E", "135°", "S", "225°", "W", "315°"],
                direction="clockwise",
                rotation=90,   # N at top
                gridcolor="#2a2a3a",
                tickfont=dict(color="white", size=11),
            ),
            radialaxis=dict(
                tickvals=list(ze_c[::3]),
                ticktext=[f"{v:.0f}°" for v in ze_c[::3]],
                tickfont=dict(color="rgba(255,255,255,0.7)", size=9),
                gridcolor="#2a2a3a",
                range=[0, float(ze_c[-1]) + float(ze_c[-1]/n_ze)],
                showgrid=True,
            ),
        ),
        margin=dict(l=60, r=80, t=60, b=40),
    )
    return fig


def _interp_color(val, vmin, vmax, colorscale_name):
    """Map a scalar value to an RGB string using a named Plotly colorscale."""
    try:
        import plotly.colors as _pc
        t = np.clip((val - vmin) / (vmax - vmin), 0, 1) if vmax > vmin else 0.0
        scale = _pc.get_colorscale(colorscale_name)
        return _pc.sample_colorscale(scale, [float(t)])[0]
    except Exception:
        # Fallback: simple blue-red
        t = np.clip((val - vmin) / (vmax - vmin), 0, 1) if vmax > vmin else 0.0
        r = int(255 * t)
        b = int(255 * (1 - t))
        return f"rgb({r},0,{b})"


def plt_colors_norm(vmin, vmax, Z, colorscale_name):
    return None  # placeholder used by colorscale logic above


def render_terrain_tab(script_dir, project_dir,
                       build_music_input_fn, build_phitsxs_input_fn,
                       build_proposal_input_fn,
                       music_materials, probe_music_file_fn, load_file_fn):
    """
    Full Tab 5: UCMuon Terrain panel.

    Parameters injected from cosmoaleph_gui.py to avoid re-importing
    everything here:
      script_dir, project_dir         — Path objects
      build_music_input_fn etc.       — stdin builder functions from main GUI
      music_materials                 — MUSIC_MATERIALS dict
      probe_music_file_fn             — probe_music_file() cache function
      load_file_fn                    — load_file() cache function
    """
    st.info(
        "**UCMuon Terrain** — Transport Generator-tab surface muons through real terrain.  \n"
        "The DEM provides per-direction rock overburden; the physics transport is done by "
        "whichever engine you select below.  No new muons are generated.",
        icon="🗺️",
    )

    ok, msg = terrain_available()
    if not ok:
        st.warning(
            f"⚠️ {msg}\n\n"
            "Install rasterio to enable DEM ray-tracing: "
            "`pip install rasterio`",
            icon="⚠️"
        )
        # Do NOT return — show the rest of the interface so the user
        # can see the settings and understand what is needed.

    st.caption(
        "**Workflow:** 📋 Setup → 🗺️ Overburden map → ▶️ Run & Results"
    )

    # ── Snapshot persistent state BEFORE creating sub-tabs ───────────────────
    # These reflect values from the PREVIOUS render (set via widget key= or
    # explicit session_state writes).  Sub-tabs cannot share local Python vars,
    # so all complex objects are stored in / read from session_state.
    _ss_lat        = st.session_state.get("terrain_lat")
    _ss_lon        = st.session_state.get("terrain_lon")
    _ss_alt        = st.session_state.get("terrain_alt")
    _ss_rho        = st.session_state.get("terrain_rho")
    _ss_naz        = st.session_state.get("terrain_naz")
    _ss_nze        = st.session_state.get("terrain_nze")
    _ss_zemax      = st.session_state.get("terrain_zemax")
    _ss_step       = st.session_state.get("terrain_step")
    _ss_dem_path   = st.session_state.get("terrain_dem_path", "")
    _ss_geom_mode  = st.session_state.get("terrain_geom_mode", "")
    _ss_synth_dem  = st.session_state.get("_terrain_synth_dem")
    _ss_csg_geom   = st.session_state.get("_terrain_csg_geom")

    # ── Create the 3 inner sub-tabs ───────────────────────────────────────────
    _ts_setup, _ts_ob, _ts_run = st.tabs([
        "📋  Setup",
        "🗺️  Overburden map",
        "▶️  Run & Results",
    ])

    # ══════════════════════════════════════════════════════════════════════════
    # TAB 1 — SETUP  (Sections 1-5)
    # ══════════════════════════════════════════════════════════════════════════
    with _ts_setup:
        st.divider()

        # ── SECTION 1: Surface muon file (from the Generator tab) ─────────────────────────
        st.markdown("### 1 — Surface muon file (Generator tab output)")
        st.caption(
            "Use the same file you would select in the Transport tab — either "
            "`muons_selected.dat` (detector filter ON) or `muons_surface.dat`."
        )

        _surf_cands = list(dict.fromkeys(f for f in [
            st.session_state.get("transport_infile", ""),
            st.session_state.get("selected_file",    ""),
            st.session_state.get("surface_file",     ""),
            "muons_selected.dat", "muons_surface.dat",
        ] if f and Path(f).exists()))

        _s1, _s2 = st.columns([3, 2])
        if _surf_cands:
            # No key= on this selectbox — avoids StreamlitAPIException when autosave
            # has a stale path that is no longer in the options list.
            t_infile = _s1.selectbox("Surface muon file", _surf_cands)
        else:
            t_infile = None
            _s1.info(
                "ℹ️  No surface muon file found yet.  "
                "Run the **Generator** tab first, then return here.",
                icon="ℹ️"
            )

        t_outfile = _s2.text_input("Output underground file", "muons_terrain_ug.dat",
                                    key="terrain_outfile")

        if t_infile and Path(t_infile).exists():
            _mt = Path(t_infile).stat().st_mtime
            ncols_t, n_transport_t = probe_music_file_fn(t_infile, False, mtime=_mt)
            st.info(f"📄 {ncols_t}-col file  |  **{n_transport_t:,}** surface muons")
        else:
            ncols_t, n_transport_t = 13, 0

        st.divider()

        # ── SECTION 2: Geometry source ────────────────────────────────────────
        st.markdown("### 2 — Geometry source")
        st.caption(
            "Select how rock overburden is defined.  "
            "**DEM** uses a real terrain elevation map.  "
            "**Synthetic** builds an analytical shape.  "
            "**PHITS** uses a volumetric CSG geometry (underground / CCS scenarios)."
        )

        _GEOM_MODE_OPTIONS = ["🗺️ DEM (GeoTIFF / Auto-download)"]
        if _GEOM_ENGINE_OK:
            _GEOM_MODE_OPTIONS.append("🔷 Synthetic Geometry")
        if _CSG_ENGINE_OK:
            _GEOM_MODE_OPTIONS.append("🔩 PHITS / STL Geometry")

        _geom_mode = st.radio(
            "Geometry mode", _GEOM_MODE_OPTIONS,
            horizontal=True, key="terrain_geom_mode",
            label_visibility="collapsed",
        )

        # ── Map _geom_mode onto the legacy _dem_mode string so the rest of the
        #    run handler (which checks _dem_mode) continues to work unchanged.
        if _geom_mode.startswith("🗺️"):
            _dem_mode = "DEM"
        elif _geom_mode.startswith("🔷"):
            _dem_mode = "🔷 Synthetic Geometry"
        else:
            _dem_mode = "🔩 CSG Geometry (STL / PHITS / MCNP)"  # run handler key

        # ════════════════════════════════════════════════════════════════════
        # DEM MODE
        # ════════════════════════════════════════════════════════════════════
        if _geom_mode.startswith("🗺️"):
            _dem_sub_opts = ["📁 Upload GeoTIFF", "⬇️ Auto-download (OpenTopography)"]
            _dem_sub = st.radio(
                "DEM source",
                _dem_sub_opts,
                horizontal=True, key="terrain_dem_mode",
            )

            if _dem_sub == "📁 Upload GeoTIFF":
                st.caption(
                    "Supported: **GeoTIFF** (.tif/.tiff), **XYZ** (.xyz), **ASC** (.asc).  "
                    "Free sources: [OpenTopography](https://opentopography.org) · "
                    "[Copernicus COP30](https://spacedata.copernicus.eu/) · "
                    "[USGS EarthExplorer](https://earthexplorer.usgs.gov/)"
                )
                def _save_dem_upload():
                    try:
                        _up = st.session_state.get("terrain_dem_upload")
                        if _up is None:
                            return
                        _up.seek(0); _bytes = _up.read()
                        if not _bytes:
                            return
                        _p = Path(tempfile.gettempdir()) / f"ucmuon_dem_{_up.name}"
                        _p.write_bytes(_bytes)
                        st.session_state["terrain_dem_path"] = str(_p)
                    except Exception:
                        pass

                st.file_uploader(
                    "DEM file (.tif  .tiff  .xyz  .asc  .txt)",
                    key="terrain_dem_upload", on_change=_save_dem_upload,
                )
                _up_saved = st.session_state.get("terrain_dem_path", "")
                if _up_saved and Path(_up_saved).exists():
                    st.success(
                        f"✅ DEM loaded: `{Path(_up_saved).name}`  "
                        f"({Path(_up_saved).stat().st_size / 1e6:.1f} MB)"
                    )

            else:  # Auto-download
                _dl1, _dl2 = st.columns(2)
                # Use terrain detector position as centre (falls back to PARMA position)
                _def_lat = _sf(
                    st.session_state.get("terrain_lat") or st.session_state.get("parma_lat", 50.668),
                    50.668,
                )
                _def_lon = _sf(
                    st.session_state.get("terrain_lon") or st.session_state.get("parma_lon", 4.615),
                    4.615,
                )
                dl_south = _dl1.number_input("South [°N]", -90.0, 90.0, _def_lat - 0.1, 0.01, key="terrain_dl_s")
                dl_north = _dl2.number_input("North [°N]", -90.0, 90.0, _def_lat + 0.1, 0.01, key="terrain_dl_n")
                dl_west  = _dl1.number_input("West  [°E]", -180.0, 180.0, _def_lon - 0.15, 0.01, key="terrain_dl_w")
                dl_east  = _dl2.number_input("East  [°E]", -180.0, 180.0, _def_lon + 0.15, 0.01, key="terrain_dl_e")
                _prods = {"SRTM GL1 (30m, global)": "SRTMGL1",
                          "SRTM GL3 (90m, global)": "SRTMGL3",
                          "COP30 (30m, Copernicus — needs personal key)": "COP30"}
                dl_prod = st.selectbox("Product", list(_prods.keys()), key="terrain_dl_prod")
                dl_out  = st.text_input("Save as", "dem_site.tif", key="terrain_dl_out")

                st.caption(
                    "SRTM GL1/GL3 are free — **paste your OpenTopography API key** below "
                    "(free account at [opentopography.org](https://opentopography.org/requestapi)).  "
                    "Leave blank to try the shared demo key (may hit rate limits)."
                )
                dl_apikey = st.text_input(
                    "OpenTopography API key", value="",
                    placeholder="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
                    type="password", key="terrain_otp_api_key",
                )

                if st.button("⬇️  Download DEM", key="terrain_dl_btn", width='stretch'):
                    _key_to_use = dl_apikey.strip() if dl_apikey.strip() else "demoapikeyot2022"
                    with st.spinner("Downloading from OpenTopography…"):
                        ok_dl, msg_dl = _download_dem(
                            dl_south, dl_north, dl_west, dl_east,
                            _prods[dl_prod], dl_out, api_key=_key_to_use)
                    if ok_dl:
                        st.success(f"✅ {msg_dl}")
                        st.session_state["terrain_dem_path"] = dl_out
                    else:
                        st.error(f"❌ {msg_dl}")
                        if "rate limit" in msg_dl.lower() or "401" in msg_dl:
                            st.info(
                                "💡 The shared demo key is rate-limited (50 calls/24 h shared "
                                "across all users). Register for a **free personal key** at "
                                "[opentopography.org/requestapi](https://opentopography.org/requestapi) "
                                "— takes ~1 minute, and your key has a much higher quota.",
                                icon="🔑",
                            )

        # ════════════════════════════════════════════════════════════════════
        # SYNTHETIC MODE
        # ════════════════════════════════════════════════════════════════════
        elif _geom_mode.startswith("🔷"):
            st.caption(
                "Build an analytical terrain shape (cone, Gaussian mound, hemisphere, slab…).  "
                "GPS coordinates are used to georeference the synthetic DEM."
            )

        # ════════════════════════════════════════════════════════════════════
        # PHITS / STL MODE
        # ════════════════════════════════════════════════════════════════════
        else:
            st.caption(
                "Load a PHITS input file (or STL mesh) describing the 3D subsurface geometry.  "
                "Each muon is ray-traced individually through the CSG cell structure.  "
                "GPS position is **not used** — the detector is defined by the cell you select below."
            )

        # ── Resolve geometry objects depending on mode ────────────────────────
        synth_dem = None
        csg_geom  = None

        if _geom_mode.startswith("🔷"):
            _sg_lat = _sf(st.session_state.get("terrain_lat",
                          st.session_state.get("parma_lat", 40.821)), 40.821)
            _sg_lon = _sf(st.session_state.get("terrain_lon",
                          st.session_state.get("parma_lon", 14.426)), 14.426)
            _sg_alt = _sf(st.session_state.get("terrain_alt",
                          st.session_state.get("parma_alt", 900.0)), 900.0)
            if _GEOM_ENGINE_OK:
                synth_dem = render_geometry_builder(_sg_lat, _sg_lon, _sg_alt)
            else:
                st.error(
                    "❌ `gui_geometry_engine.py` not found — place it alongside "
                    "`gui_terrain_engine.py` and restart."
                )
            dem_path = None

        elif _geom_mode.startswith("🔩"):
            if _CSG_ENGINE_OK:
                csg_geom = render_csg_builder()
            else:
                st.error(
                    "❌ `gui_csg_engine.py` not found — place it alongside "
                    "`gui_terrain_engine.py` and restart."
                )
            dem_path = None

        else:
            # File-based DEM path resolution (original logic unchanged)
            _dem_ss = st.session_state.get("terrain_dem_path", "")
            if _dem_ss and not Path(_dem_ss).exists():
                # Stale path — file was deleted (temp file from a previous session)
                st.warning(
                    f"⚠️  Previously loaded DEM `{Path(_dem_ss).name}` no longer exists on disk.  "
                    "Please upload or download the DEM again.",
                    icon="⚠️"
                )
                st.session_state.pop("terrain_dem_path", None)
                _dem_ss = ""

            # No DEM selected yet — fall back to the bundled Vesuvius sample.
            # Written to session_state so the Overburden-map and Run sub-tabs
            # (which read terrain_dem_path) see it too.
            if not _dem_ss and _BUNDLED_VESUVIUS_DEM.exists():
                _dem_ss = str(_BUNDLED_VESUVIUS_DEM)
                st.session_state["terrain_dem_path"] = _dem_ss

            if _dem_ss and Path(_dem_ss).exists():
                dem_path = _dem_ss
                if Path(dem_path) == _BUNDLED_VESUVIUS_DEM:
                    st.caption(
                        f"📍 Default DEM: bundled Mt. Vesuvius sample "
                        f"`{_BUNDLED_VESUVIUS_DEM.name}` (SRTM GL1 30 m, "
                        f"14.35–14.52 °E / 40.76–40.90 °N).  "
                        f"Press **🌋 muRAvES / Vesuvius example** below to set matching "
                        f"detector coordinates, or upload / download a DEM for your own site."
                    )
                else:
                    st.caption(f"📍 Current DEM: `{Path(dem_path).name}`  "
                               f"({Path(dem_path).stat().st_size/1e6:.1f} MB)")
            else:
                _manual_dem = st.text_input("Or enter path to existing GeoTIFF", value="",
                                            key="terrain_dem_manual", placeholder="/path/to/dem.tif")
                if _manual_dem and Path(_manual_dem).exists():
                    dem_path = _manual_dem
                    st.session_state["terrain_dem_path"] = _manual_dem
                elif _manual_dem:
                    st.error(f"❌ Not found: `{_manual_dem}`")
                    dem_path = None
                else:
                    dem_path = None

        # ── Persist complex objects to session_state so other tabs can read them
        if synth_dem is not None:
            st.session_state["_terrain_synth_dem"] = synth_dem
        if csg_geom is not None:
            st.session_state["_terrain_csg_geom"] = csg_geom

        st.divider()

        # ── SECTION 3: Detector position ──────────────────────────────────────
        _is_csg_mode = (csg_geom is not None)

        # ── PHITS mode: skip GPS entirely, go straight to cell selection ──────
        det_cell_id = None
        det_lat = _sf(st.session_state.get("parma_lat",  40.832), 40.832)
        det_lon = _sf(st.session_state.get("parma_lon",  14.412), 14.412)
        det_alt = _sf(st.session_state.get("parma_alt",  784.0),  784.0)
        underground = False

        if _is_csg_mode:
            # ── PHITS / STL: detector defined by cell selection only ──────────
            st.divider()
            st.markdown("### 3 — Detector Cell")
            st.caption(
                "Pick the geometry cell that acts as the detector volume.  "
                "Only muons whose ray physically enters this cell will be recorded.  "
                "GPS coordinates are **not used** — the geometry defines the rock structure."
            )
            if _CSG_TRANSPORT_OK:
                det_cell_id = render_detector_cell_selector(csg_geom)
            else:
                st.error(
                    "❌ `gui_csg_transport.py` not found — place it alongside "
                    "`gui_terrain_engine.py` and restart."
                )
            underground = True  # PHITS transport is always underground

        else:
            # ── DEM / Synthetic: show GPS + underground checkbox ───────────────
            st.divider()
            st.markdown("### 3 — Detector GPS position")
            st.caption(
                "Enter the detector location as GPS coordinates.  "
                "The DEM is sampled around this point to compute rock overburden "
                "in each direction."
            )

            # muRAvES / Vesuvius example preset button
            if st.button("🌋  muRAvES / Vesuvius example", key="terrain_muraves_preset"):
                # Detector = MURAVES site, SW flank (examples/vesuvius/MURAVES_GUIDE.md)
                st.session_state.update({
                    "terrain_lat": 40.8271, "terrain_lon": 14.4006,
                    "terrain_alt": 608.0, "terrain_rho": 2.65,
                    "terrain_naz": 36, "terrain_nze": 18,
                    "terrain_zemax": 85.0, "terrain_step": 50.0,
                })
                if _BUNDLED_VESUVIUS_DEM.exists():
                    st.session_state["terrain_dem_path"] = str(_BUNDLED_VESUVIUS_DEM)
                st.rerun()

            _gp1, _gp2, _gp3 = st.columns(3)
            det_lat = _gp1.number_input(
                "Latitude [°N]", -90.0, 90.0,
                _sf(st.session_state.get("parma_lat", 50.668), 50.668),
                0.0001, format="%.6f", key="terrain_lat",
            )
            det_lon = _gp2.number_input(
                "Longitude [°E]", -180.0, 180.0,
                _sf(st.session_state.get("parma_lon", 4.615), 4.615),
                0.0001, format="%.6f", key="terrain_lon",
            )
            det_alt = _gp3.number_input(
                "Altitude [m a.s.l.]", -500.0, 9000.0,
                _sf(st.session_state.get("parma_alt", 90.0), 90.0),
                1.0, key="terrain_alt",
            )

            underground = st.checkbox(
                "🔽  Underground detector",
                value=st.session_state.get("terrain_underground", False),
                key="terrain_underground",
                help=(
                    "Enable when the detector is BELOW the surface "
                    "(borehole, mine drift, tunnel cavern).  "
                    "Disables altitude clamping; the DEM defines the rock–air boundary above.  "
                    "Set Altitude to the true underground position "
                    "(e.g. −450 m for a borehole at 450 m depth below sea level)."
                ),
            )
            if underground:
                st.info(
                    "ℹ️  Underground mode: altitude clamping disabled.  "
                    "The altitude you entered is used as-is as the detector depth.",
                    icon="🔽",
                )

            # ── DEM quick-check (elevation profile + 3D view) ─────────────────
            if dem_path or synth_dem is not None:
                st.divider()
                st.markdown("### 4 — DEM check")
                import plotly.graph_objects as _go_chk
                try:
                    if synth_dem is not None:
                        _ec = synth_dem.elev.astype(float)
                        _tc = synth_dem.transform
                        st.caption(f"🔷 Synthetic DEM: {synth_dem.summary()}")
                    else:
                        _drv_c = _load_terrain_driver(script_dir)
                        _ec, _tc = _drv_c.load_dem(dem_path)
                        _ec = _ec.astype(float)

                    _hc, _wc = _ec.shape
                    _lon_min_c = float(_tc.c)
                    _lon_max_c = float(_tc.c + _wc * _tc.a)
                    _lat_max_c = float(_tc.f)
                    _lat_min_c = float(_tc.f + _hc * _tc.e)
                    if _lat_min_c > _lat_max_c:
                        _lat_min_c, _lat_max_c = _lat_max_c, _lat_min_c

                    _cm1, _cm2, _cm3, _cm4 = st.columns(4)
                    _cm1.metric("DEM lat", f"{_lat_min_c:.3f}° – {_lat_max_c:.3f}°")
                    _cm2.metric("DEM lon", f"{_lon_min_c:.3f}° – {_lon_max_c:.3f}°")
                    _cm3.metric("Min elev", f"{float(np.nanmin(_ec)):.0f} m")
                    _cm4.metric("Max elev", f"{float(np.nanmax(_ec)):.0f} m")

                    _inside_c = (_lat_min_c < det_lat < _lat_max_c and
                                 _lon_min_c < det_lon < _lon_max_c)
                    if _inside_c:
                        if synth_dem is None:
                            _det_elev_c = _drv_c.dem_elevation_at(
                                _ec.astype(np.float32), _tc, det_lat, det_lon)
                        else:
                            _det_elev_c = float(np.interp(
                                0.0,
                                [0.0], [float(_ec[_hc // 2, _wc // 2])],
                            ))
                        _diff_c = abs(det_alt - _det_elev_c)
                        if _diff_c < 50:
                            st.success(
                                f"✅ Detector inside DEM.  "
                                f"DEM elevation at detector = **{_det_elev_c:.0f} m**  "
                                f"(entered {det_alt:.0f} m,  Δ = {_diff_c:.0f} m)"
                            )
                        else:
                            st.warning(
                                f"⚠️ Altitude mismatch: DEM says **{_det_elev_c:.0f} m** "
                                f"at your position, you entered **{det_alt:.0f} m** "
                                f"(Δ = {_diff_c:.0f} m).  Consider updating the altitude.",
                                icon="⚠️",
                            )
                    else:
                        st.error(
                            f"❌ Detector ({det_lat:.4f}°N, {det_lon:.4f}°E) is **outside** "
                            f"the DEM (lat {_lat_min_c:.3f}–{_lat_max_c:.3f}, "
                            f"lon {_lon_min_c:.3f}–{_lon_max_c:.3f}).  "
                            f"Download a DEM centred on your detector coordinates."
                        )

                    # N–S elevation profile
                    if synth_dem is None:
                        _lats_c = np.linspace(_lat_min_c, _lat_max_c, 250)
                        _elvs_c = np.array([
                            _drv_c.dem_elevation_at(_ec.astype(np.float32), _tc, la, det_lon)
                            for la in _lats_c
                        ], dtype=float)
                        _fig_p = _go_chk.Figure()
                        _fig_p.add_trace(_go_chk.Scatter(
                            x=_lats_c, y=_elvs_c, mode="lines",
                            line=dict(color="#00b4d8", width=2),
                            name=f"Elevation at lon={det_lon:.3f}°",
                        ))
                        _fig_p.add_vline(x=det_lat,
                                         line=dict(color="#ffd700", width=2, dash="dash"),
                                         annotation_text=f"Detector {det_lat:.3f}°N",
                                         annotation_font=dict(color="#ffd700"))
                        _fig_p.add_hline(y=det_alt,
                                         line=dict(color="#ff6b6b", width=1.5, dash="dot"),
                                         annotation_text=f"Alt {det_alt:.0f} m",
                                         annotation_font=dict(color="#ff6b6b"))
                        _fig_p.update_layout(
                            **DARK, height=260,
                            xaxis=dict(title="Latitude [°N]", gridcolor="#2a2a3a"),
                            yaxis=dict(title="Elevation [m a.s.l.]", gridcolor="#2a2a3a"),
                            title=dict(
                                text=f"N–S elevation profile at lon={det_lon:.3f}°",
                                font=dict(size=12),
                            ),
                            margin=dict(l=60, r=20, t=45, b=45),
                        )
                        st.plotly_chart(_fig_p,                                         config={"displayModeBar": False},
                                        key="terrain_setup_profile")
                        st.caption(
                            "The summit should appear as a peak in the target direction.  "
                            "If the profile is flat the DEM does not cover the target."
                        )

                except Exception as _ce:
                    st.warning(f"DEM check failed: {_ce}")

                # ── 3D terrain view ────────────────────────────────────────────
                st.markdown("**3D terrain view**")
                _r3d_col1, _r3d_col2 = st.columns([4, 1])
                _r3d_km = _r3d_col2.slider(
                    "Radius [km]", 1, 30, 10, 1, key="terrain_setup_3d_radius"
                )
                _r3d_deg = _r3d_km / 111.0
                _ea3 = synth_dem.elev if synth_dem is not None else None
                _ta3 = synth_dem.transform if synth_dem is not None else None
                _fig3 = _dem_3d_plot(
                    dem_path, det_lat, det_lon, det_alt, script_dir,
                    radius_deg=_r3d_deg,
                    elev_arr=_ea3, tfm_arr=_ta3,
                )
                if isinstance(_fig3, tuple) and _fig3[0] == "ERROR":
                    st.warning(f"3D terrain: {_fig3[1]}")
                elif _fig3 is not None:
                    with _r3d_col1:
                        st.plotly_chart(_fig3,                                         config={"displayModeBar": True},
                                        key="terrain_setup_3d")

        st.divider()

        # ── SECTION 5: Physics — engine + material ────────────────────────────
        st.markdown("### 5 — Transport engine & material")
        st.caption(
            "The DEM provides per-direction overburden depth. "
            "The physics transport uses the engine you select here. "
            "**UCMuon-MC is recommended** — it processes each muon individually "
            "with its exact per-direction depth without needing separate subprocess calls."
        )

        _TERRAIN_ENGINE_DESC = {
            "★ UCMuon-MC (Python)  ← recommended for terrain": (
                "**UCMuon-MC — native stochastic MC** — Groom dE/dx + Poisson radiative losses "
                "+ Highland MS + decay.  "
                "Pure Python, processes each muon at its exact per-direction depth."
            ),
            "MUSIC": (
                "**Full stochastic MC** — Kudryavtsev (2009) XS tables "
                "(ionisation, bremsstrahlung, pair production, photonuclear).  "
                "Fortran/OMP; spawns one subprocess per direction bin."
            ),
            "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS": (
                "**CSDA** — PDG Bethe-Bloch ionisation + Groom radiative + Highland MS.  "
                "Fortran/OMP; no external files needed."
            ),
            "PROPOSAL": (
                "**Full stochastic MC** — PROPOSAL library (IceCube/KM3NeT).  "
                "Requires `pip install proposal`; slowest but most complete physics."
            ),
        }

        # Migrate sessions saved before the UCMuon-MC rename (old label string)
        if st.session_state.get("terrain_engine_choice") not in _ENGINE_OPTIONS:
            st.session_state.pop("terrain_engine_choice", None)

        _ph1, _ph2 = st.columns(2)
        terrain_engine_label = _ph1.radio(
            "Transport engine",
            options=_ENGINE_OPTIONS,
            key="terrain_engine_choice",
        )
        terrain_engine = _ENGINE_INTERNAL[terrain_engine_label]
        st.caption(_TERRAIN_ENGINE_DESC.get(terrain_engine_label, ""))

        if _is_csg_mode and csg_geom is not None:
            # ── CSG mode: density is defined per-cell in the .inp file ─────────
            # The terrain_rho widget has zero effect here — density_at() reads the
            # parsed cell densities directly.  Show a read-only summary instead so
            # the user knows exactly what densities will be used.
            _cell_mats = [
                c for c in csg_geom._cells
                if c.density > 0 and c.mat_id not in (-1, 0)
            ]
            # Fallback: read terrain_rho from session_state so downstream code that
            # still references it (overburden preview, Tab 3 plots) gets a sensible
            # value — use the mean density across non-void, non-detector cells.
            _non_det_cells = [
                c for c in _cell_mats
                if c.cell_id not in st.session_state.get("_csg_selected_det_cells", [])
            ]
            _mean_rho = (
                float(np.mean([c.density for c in _non_det_cells]))
                if _non_det_cells else 2.0
            )
            # Keep session_state in sync so Tab 3 / overburden display works
            st.session_state["terrain_rho"] = _mean_rho
            terrain_rho = _mean_rho

            with _ph2:
                st.markdown("**Material densities (from .inp)**")
                # Build colour-coded rows: one per distinct material
                _seen_mats = {}
                for c in _cell_mats:
                    if c.mat_id not in _seen_mats:
                        _seen_mats[c.mat_id] = (c.density, c.label)
                # Pick a colour ramp matching plotly_preview (blue→amber→red)
                def _rho_swatch(rho):
                    if rho < 1.2:   return "#4e9af1"
                    if rho < 1.6:   return "#5ab552"
                    if rho < 1.8:   return "#d4a017"
                    if rho < 2.0:   return "#e07b39"
                    if rho < 2.2:   return "#c94040"
                    return "#8b1a1a"
                _rows_html = "".join(
                    f"<div style='display:flex;align-items:center;gap:6px;"
                    f"margin:2px 0;font-size:0.82em;'>"
                    f"<span style='width:10px;height:10px;border-radius:2px;"
                    f"background:{_rho_swatch(rho)};display:inline-block;'></span>"
                    f"<span style='color:var(--text-color);'>"
                    f"mat {mid} — {rho:.3f} g/cm³"
                    f"</span></div>"
                    for mid, (rho, _lbl) in sorted(_seen_mats.items())
                )
                _mean_line = (
                    f"<div style='margin-top:4px;font-size:0.78em;"
                    f"color:var(--text-color);opacity:0.7;'>"
                    f"Mean (non-detector): {_mean_rho:.3f} g/cm³ — used for overburden display</div>"
                )
                st.markdown(
                    f"<div style='border:1px solid var(--border-color);"
                    f"border-radius:6px;padding:8px 10px;"
                    f"background:var(--secondary-background-color);'>"
                    f"{_rows_html}{_mean_line}</div>",
                    unsafe_allow_html=True,
                )
        else:
            # ── DEM / Synthetic mode: single uniform density for the whole column
            terrain_rho = _ph2.number_input(
                "Rock density ρ [g/cm³]", 0.1, 10.0,
                _sf(st.session_state.get("music_rho", 2.65), 2.65), 0.05,
                key="terrain_rho",
                help=(
                    "Applied uniformly along all rock paths. "
                    "Use the mean density of the overburden column.\n\n"
                    "Standard Rock: 2.65  |  Limestone: 2.70  |  "
                    "Basalt: 2.85  |  Ice: 0.917"
                ),
            )

        # UCMuon-MC extra settings (internal token "UCMuon Stochastic")
        if terrain_engine == "UCMuon Stochastic":
            with st.expander("⚙️  UCMuon-MC settings", expanded=False):
                _ps1, _ps2 = st.columns(2)
                t_vcut = _ps1.number_input("v_cut", 0.01, 0.5, 0.05, 0.01, key="terrain_vcut",
                                            help="Catastrophic event threshold (fraction of E). Lower = more accurate, slower.")
                t_ms   = _ps2.checkbox("Highland MS deflections", value=True, key="terrain_ms")
        else:
            t_vcut = 0.05
            t_ms   = True

        # Engine availability warnings
        if terrain_engine == "MUSIC":
            if not (project_dir / "ucmuon_transport_music_omp").exists():
                st.warning("⚠️  MUSIC driver not found — run `bash setup.sh`")
        elif terrain_engine == "Bethe-Bloch":
            if not (project_dir / "ucmuon_transport_bb_omp").exists():
                st.warning("⚠️  BB driver not found — run `bash setup.sh`")
        elif terrain_engine == "PROPOSAL":
            _prop_check = subprocess.run(
                [sys.executable, "-c", "import proposal; print(proposal.__version__)"],
                capture_output=True, text=True, timeout=10,
            )
            if _prop_check.returncode != 0:
                st.warning("⚠️  PROPOSAL not importable — run: `pip install proposal`")

        with st.expander("⚙️  Angular grid settings", expanded=False):
            _adv1, _adv2, _adv3 = st.columns(3)
            n_az    = int(_adv1.number_input("Azimuth bins",   8, 360, 36, 4,   key="terrain_naz"))
            n_ze    = int(_adv2.number_input("Zenith bins",    4, 90,  18, 2,   key="terrain_nze"))
            ze_max  = _adv3.number_input("Max zenith [°]", 20.0, 89.0, 85.0, 5.0, key="terrain_zemax")
            step_m  = _adv1.number_input("Ray-trace step [m]", 10.0, 500.0, 50.0, 10.0, key="terrain_step")
            _ze_last_c = ze_max - ze_max / n_ze / 2.0
            st.caption(
                f"Grid: {n_az} az x {n_ze} ze = {n_az*n_ze} direction bins.  "
                f"Last bin centre = {_ze_last_c:.1f} deg.  "
                f"Muons assigned by their direction; empty bins skipped."
            )
            if ze_max < 80.0:
                st.warning(
                    "Max zenith < 80 deg — terrain blocking near the summit "
                    "typically occurs only at ze > 70-80 deg. "
                    "Recommended: set Max zenith to 85 deg.",
                    icon="⚠️"
                )

    # end with _ts_setup

    # ══════════════════════════════════════════════════════════════════════════
    # TAB 2 — OVERBURDEN MAP  (Section 3b DEM validation + overburden preview)
    # ══════════════════════════════════════════════════════════════════════════
    with _ts_ob:
        # Read geometry state from session_state (set by Setup tab)
        _ob_dem_path  = st.session_state.get("terrain_dem_path", "")
        _ob_synth_dem = st.session_state.get("_terrain_synth_dem")
        _ob_csg_geom  = st.session_state.get("_terrain_csg_geom")
        _ob_lat       = _sf(st.session_state.get("terrain_lat",   50.668), 50.668)
        _ob_lon       = _sf(st.session_state.get("terrain_lon",    4.615),  4.615)
        _ob_alt       = _sf(st.session_state.get("terrain_alt",   90.0),   90.0)
        _ob_rho       = _sf(st.session_state.get("terrain_rho",   2.65),   2.65)
        _ob_naz       = int(st.session_state.get("terrain_naz",   36))
        _ob_nze       = int(st.session_state.get("terrain_nze",   18))
        _ob_zemax     = _sf(st.session_state.get("terrain_zemax", 85.0),   85.0)
        _ob_step      = _sf(st.session_state.get("terrain_step",  50.0),   50.0)
        _ob_underground = bool(st.session_state.get("terrain_underground", False))

        st.divider()
        st.info(
            "ℹ️  Elevation profile and 3D terrain view are in the **📋 Setup** tab "
            "(section 4 — DEM check), so you can validate and adjust before computing here.",
            icon="🔍",
        )

        # ── Overburden preview button ─────────────────────────────────────────
        if _ob_dem_path or (_ob_synth_dem is not None) or (_ob_csg_geom is not None):
            if _ob_csg_geom is not None:
                _prev_label = "👁️  Compute overburden map (CSG preview)"
            elif _ob_synth_dem is not None:
                _prev_label = "👁️  Compute overburden map (synthetic DEM preview, ~30s)"
            else:
                _prev_label = "👁️  Compute overburden map (DEM preview, ~30s)"
            if st.button(_prev_label, key="terrain_preview_btn", width='stretch'):
                with st.spinner("Ray tracing overburden map…"):
                    try:
                        if _ob_csg_geom is not None:
                            _csg_det_p = st.session_state.get("_csg_det_pos", np.array([0.,0.,0.]))
                            az_c, ze_c, ob_map, sky_map = _ob_csg_geom.compute_overburden_map(
                                det_pos_m  = _csg_det_p,
                                n_az       = _ob_naz,    n_ze       = _ob_nze,
                                ze_max_deg = _ob_zemax,  step_m     = float(st.session_state.get("_csg_step", 0.5)),
                                max_dist_m = float(st.session_state.get("_csg_maxdist", 5000.0)),
                            )
                        else:
                            drv = _load_terrain_driver(script_dir)
                            if _ob_synth_dem is not None:
                                elev, transform = _ob_synth_dem.elev, _ob_synth_dem.transform
                            else:
                                elev, transform = drv.load_dem(_ob_dem_path)
                            az_c, ze_c, ob_map, sky_map = drv.compute_overburden_map(
                                elev, transform, _ob_lat, _ob_lon, _ob_alt,
                                _ob_rho, _ob_naz, _ob_nze, _ob_zemax, _ob_step,
                                underground=_ob_underground,
                            )
                        st.session_state["terrain_preview"] = (az_c, ze_c, ob_map, sky_map)
                        st.success("✅  Overburden map computed.")
                    except Exception as _e:
                        st.error(f"❌  Preview failed: {_e}")
                        import traceback; st.code(traceback.format_exc())
        else:
            st.info(
                "ℹ️  Configure geometry in the **📋 Setup** tab first, "
                "then return here to compute the overburden map.",
                icon="ℹ️"
            )

        if "terrain_preview" in st.session_state:
            _prev_az_c, _prev_ze_c, ob_p, sky_p = st.session_state["terrain_preview"]
            ob_rock = ob_p[~sky_p]
            _v1, _v2, _v3, _v4 = st.columns(4)
            _v1.metric("Rock directions",  f"{(~sky_p).sum()}/{sky_p.size}")
            _v2.metric("Open sky",         f"{sky_p.sum()}/{sky_p.size}")
            _v3.metric("Max overburden",   f"{float(ob_p.max()):.0f} g/cm²")
            _v4.metric("Median overburden",f"{float(np.median(ob_rock)):.0f} g/cm²" if ob_rock.size else "—")
            if (~sky_p).sum() == 0:
                _ze_last = float(_prev_ze_c[-1]) if len(_prev_ze_c) else 0.0
                st.warning(
                    f"All {sky_p.size} direction bins show open sky — no terrain blocking detected.  "
                    f"The zenith grid only reaches **{_ze_last:.1f}°** (last bin centre).  "
                    f"**Fix:** increase Max zenith to **85°** and recompute.",
                    icon="⚠️"
                )
            if ob_rock.size > 0:
                _pv1b, _pv2b = st.columns(2)
                _pv1b.caption(
                    f"Blocked bin range:  "
                    f"{float(ob_rock.min()):.0f} – {float(ob_rock.max()):.0f} g/cm²  "
                    f"({float(ob_rock.min())/(_ob_rho*100):.0f} – "
                    f"{float(ob_rock.max())/(_ob_rho*100):.0f} m vertical equiv.)"
                )
                _pv2b.caption(
                    "How overburden is computed: for each (azimuth, zenith) bin the DEM "
                    "ray tracer shoots a ray from the detector in that direction and "
                    "integrates the column of rock it crosses until reaching open sky.  "
                    "Result: overburden in g/cm² = path length [cm] × density [g/cm³]."
                )
            st.plotly_chart(
                _polar_heatmap(_prev_az_c, _prev_ze_c, ob_p,
                               "Rock overburden — log scale [g/cm²]",
                               "g/cm²", colorscale="Jet", mask=sky_p, log_scale=True),
                config={"displayModeBar": False},
                key="terrain_ob_preview"
            )
            _ob_min_prev = float(ob_rock.min()) if ob_rock.size > 0 else 1.0
            _ob_max_prev = float(ob_p.max())
            _slant_prev  = _ob_max_prev / (_ob_rho * 100.0) if _ob_rho > 0 else 0
            _dec_prev    = np.log10(max(_ob_max_prev, 1.0) / max(_ob_min_prev, 1.0)) if _ob_min_prev > 0 else 0
            st.caption(
                f"Colour = log₁₀(overburden g/cm²). Hover for actual value.  "
                f"Gray = open sky. Red = max blocking "
                f"({_ob_max_prev:,.0f} g/cm² = {_slant_prev:.0f} m slant).  "
                f"Log scale spans {_dec_prev:.1f} orders of magnitude."
            )

    # end with _ts_ob

    # ══════════════════════════════════════════════════════════════════════════
    # TAB 3 — RUN & RESULTS  (readiness check + run button + all results)
    # ══════════════════════════════════════════════════════════════════════════
    with _ts_run:
        # Read all needed variables from session_state
        _run_dem_path   = st.session_state.get("terrain_dem_path", "")
        _run_synth_dem  = st.session_state.get("_terrain_synth_dem")
        _run_csg_geom   = st.session_state.get("_terrain_csg_geom")
        _run_lat        = _sf(st.session_state.get("terrain_lat",   50.668), 50.668)
        _run_lon        = _sf(st.session_state.get("terrain_lon",    4.615),  4.615)
        _run_alt        = _sf(st.session_state.get("terrain_alt",   90.0),   90.0)
        _run_rho        = _sf(st.session_state.get("terrain_rho",   2.65),   2.65)
        _run_naz        = int(st.session_state.get("terrain_naz",   36))
        _run_nze        = int(st.session_state.get("terrain_nze",   18))
        _run_zemax      = _sf(st.session_state.get("terrain_zemax", 85.0),   85.0)
        _run_step       = _sf(st.session_state.get("terrain_step",  50.0),   50.0)
        _run_underground = bool(st.session_state.get("terrain_underground", False))
        _run_engine_label = st.session_state.get("terrain_engine_choice", _ENGINE_OPTIONS[0])
        _run_engine     = _ENGINE_INTERNAL.get(_run_engine_label, "UCMuon Stochastic")
        _run_vcut       = _sf(st.session_state.get("terrain_vcut", 0.05), 0.05)
        _run_ms         = bool(st.session_state.get("terrain_ms", True))
        _run_outfile    = st.session_state.get("terrain_outfile", "muons_terrain_ug.dat")
        _run_is_csg     = (_run_csg_geom is not None)

        # Resolve surface muon file + probe from session_state
        _run_surf_cands = list(dict.fromkeys(f for f in [
            st.session_state.get("transport_infile", ""),
            st.session_state.get("selected_file",    ""),
            st.session_state.get("surface_file",     ""),
            "muons_selected.dat", "muons_surface.dat",
        ] if f and Path(f).exists()))
        _run_infile = _run_surf_cands[0] if _run_surf_cands else None
        if _run_infile and Path(_run_infile).exists():
            _mt_r = Path(_run_infile).stat().st_mtime
            _run_ncols, _run_n = probe_music_file_fn(_run_infile, False, mtime=_mt_r)
        else:
            _run_ncols, _run_n = 13, 0

        # Detector cell
        _run_det_cell_id = None
        _run_selected_det_cells = st.session_state.get("_csg_selected_det_cells", [])

        st.divider()

        # ── SECTION 6: Run ────────────────────────────────────────────────────
        st.markdown("### 6 — Run terrain transport")

        # Retrieve the full list of selected detector cells (multi-detector)
        if _run_is_csg and _run_det_cell_id is not None and not _run_selected_det_cells:
            _run_selected_det_cells = [_run_det_cell_id]

        _dem_ready    = bool(_run_dem_path) or (_run_synth_dem is not None) or (_run_csg_geom is not None)
        _cell_ready   = (_run_csg_geom is None) or (len(_run_selected_det_cells) > 0)
        _muons_ready  = (bool(_run_infile) and _run_n > 0) or (_run_csg_geom is not None)
        ready = _dem_ready and _cell_ready and _muons_ready
        if not ready:
            missing_items = []
            if not _dem_ready:   missing_items.append("geometry — upload a DEM, Synthetic, or PHITS (📋 Setup tab)")
            if not _cell_ready:  missing_items.append("detector cell — select one in 📋 Setup tab")
            if not _muons_ready and _run_csg_geom is None:
                missing_items.append("surface muon file (📋 Setup tab)")
            st.warning("⚠️  Not ready — missing: " + ", ".join(missing_items))

        if _run_csg_geom is not None and _run_selected_det_cells:
            _det_names = st.session_state.get("_csg_detector_names", {})
            _cells_str = " + ".join(
                f"cell {c}" + (f" ({_det_names[c]})" if c in _det_names else "")
                for c in _run_selected_det_cells
            )
            _run_label = f"▶  Run CSG transport  ({_run_n:,} muons → {_cells_str})"
        else:
            _run_label = f"▶  Run terrain transport  ({_run_n:,} muons × {_run_naz*_run_nze} direction bins)"

        run_btn = st.button(
            _run_label, type="primary", width='stretch',
            disabled=(not ready), key="terrain_run_btn",
        )

        # Convenience aliases used by the run handler below (same names as before)
        dem_path    = _run_dem_path
        synth_dem   = _run_synth_dem
        csg_geom    = _run_csg_geom
        det_lat     = _run_lat
        det_lon     = _run_lon
        det_alt     = _run_alt
        terrain_rho = _run_rho
        n_az        = _run_naz
        n_ze        = _run_nze
        ze_max      = _run_zemax
        step_m      = _run_step
        underground = _run_underground
        terrain_engine = _run_engine
        t_vcut      = _run_vcut
        t_ms        = _run_ms
        t_infile    = _run_infile
        t_outfile   = _run_outfile
        ncols_t     = _run_ncols
        n_transport_t = _run_n
        _is_csg_mode  = _run_is_csg
        _selected_det_cells = _run_selected_det_cells
        det_cell_id = _run_det_cell_id

        if run_btn and ready:
            with st.spinner(f"Running terrain transport ({terrain_engine})…"):
                t0 = time.time()
                prog = st.progress(0.0, text="⏳  Starting…")

                try:
                    # ── CSG mode: per-muon ray-trace transport ────────────────────
                    if csg_geom is not None:
                        if not _CSG_TRANSPORT_OK:
                            raise ImportError(
                                "gui_csg_transport.py not found.  "
                                "Place it alongside gui_terrain_engine.py and restart."
                            )
                        if not _selected_det_cells:
                            raise ValueError(
                                "No detector cells selected.  "
                                "Choose cells in Section 3 and re-run."
                            )

                        # Load surface muon file
                        if t_infile and Path(t_infile).exists():
                            df_surf = load_file_fn(t_infile,
                                                   mtime=Path(t_infile).stat().st_mtime)
                        else:
                            raise FileNotFoundError(
                                "Surface muon file not found — run the Generator tab first."
                            )

                        if "cx" not in df_surf.columns:
                            # Prefer momentum components (px/py/pz columns) for exact
                            # direction cosines — more accurate than re-deriving from
                            # the stored theta/phi angles.  Fall back to trig if absent.
                            if "px" in df_surf.columns and "p" in df_surf.columns:
                                _p = df_surf["p"].values.astype(float)
                                _p = np.where(_p > 0, _p, 1.0)   # guard /0
                                df_surf["cx"] = df_surf["px"].values / _p
                                df_surf["cy"] = df_surf["py"].values / _p
                                df_surf["cz"] = df_surf["pz"].values / _p
                            else:
                                df_surf["cx"] = np.sin(df_surf["theta"]) * np.cos(df_surf["phi"])
                                df_surf["cy"] = np.sin(df_surf["theta"]) * np.sin(df_surf["phi"])
                                df_surf["cz"] = -np.cos(df_surf["theta"])

                        _csg_step_r  = float(st.session_state.get("_csg_step_cm", 50.0))
                        _csg_maxd_r  = float(st.session_state.get("_csg_maxdist_cm", 500_000.0))
                        _csg_off     = st.session_state.get("_csg_offset_cm", np.zeros(3))

                        # Build detector name map for dump files
                        _det_name_map = {}
                        for cid in _selected_det_cells:
                            _det_name_map[cid] = st.session_state.get(
                                "_csg_detector_names", {}).get(cid, f"cell{cid}")

                        st.caption(
                            f"🔩 CSG: {csg_geom.summary()}  |  "
                            f"Detectors: {list(_selected_det_cells)}  |  "
                            f"Step: {_csg_step_r:.0f} cm  |  "
                            f"Offset: {np.asarray(_csg_off)} cm"
                        )

                        prog.progress(0.05, text="⏳  Per-muon CSG ray tracing (multi-detector)…")

                        # ── Multi-detector transport (single pass) ────────────────
                        _prog_container = st.empty()
                        result_dfs = transport_muons_multi_detector(
                            surface_df        = df_surf,
                            geom              = csg_geom,
                            detector_cell_ids = _selected_det_cells,
                            detector_names    = _det_name_map,
                            coord_offset_cm   = np.asarray(_csg_off, dtype=float),
                            step_cm           = _csg_step_r,
                            max_dist_cm       = _csg_maxd_r,
                            v_cut             = t_vcut,
                            ms_enable         = t_ms,
                            script_dir        = str(script_dir),
                            write_phits_dump  = True,
                            dump_dir          = str(script_dir),
                            progress_container = _prog_container,
                        )

                        # Store results per-detector
                        st.session_state["terrain_csg_multi_result"] = result_dfs
                        st.session_state["terrain_result_csg"]       = csg_geom
                        st.session_state["terrain_csg_det_cells"]    = _selected_det_cells
                        st.session_state["terrain_csg_det_names"]    = _det_name_map

                        # Primary detector result (backward compat)
                        df_ug = result_dfs[_selected_det_cells[0]]
                        det_cell_id = _selected_det_cells[0]

                        # Build placeholder overburden map
                        n_az_r = n_az; n_ze_r = n_ze
                        az_c  = np.linspace(0,360,n_az_r,endpoint=False)+180/n_az_r
                        ze_c  = np.linspace(0,ze_max,n_ze_r,endpoint=False)+ze_max/(2*n_ze_r)
                        ob_map  = np.zeros((n_az_r, n_ze_r))
                        sky_map = np.ones((n_az_r, n_ze_r), dtype=bool)
                        _surv = df_ug[df_ug["alive"] == 1]
                        if len(_surv) > 0:
                            _az_idx, _ze_idx = assign_direction_bins(
                                _surv, n_az_r, n_ze_r, ze_max)
                            for _ia, _iz in zip(_az_idx, _ze_idx):
                                if 0 <= _ia < n_az_r and 0 <= _iz < n_ze_r:
                                    sky_map[_ia, _iz] = False

                        st.session_state["terrain_preview"]         = (az_c, ze_c, ob_map, sky_map)
                        st.session_state["terrain_result_det_cell"] = det_cell_id

                    else:
                        drv = _load_terrain_driver(script_dir)
                        # DEM source: synthetic geometry or file-based
                        if synth_dem is not None:
                            elev      = synth_dem.elev
                            transform = synth_dem.transform
                            st.caption(f"🔷 Synthetic DEM: {synth_dem.summary()}")
                        else:
                            elev, transform = drv.load_dem(dem_path)
                        az_c, ze_c, ob_map, sky_map = drv.compute_overburden_map(
                            elev, transform, det_lat, det_lon, det_alt,
                            terrain_rho, n_az, n_ze, ze_max, step_m,
                            underground=underground,
                        )
                    st.session_state["terrain_preview"] = (az_c, ze_c, ob_map, sky_map)

                    prog.progress(0.1, text="⏳  Geometry done — loading muons…")

                    # ── Surface muon file — DEM path only (CSG already produced df_ug) ──
                    if csg_geom is None:
                        if t_infile and Path(t_infile).exists():
                            df_surf = load_file_fn(t_infile,
                                                   mtime=Path(t_infile).stat().st_mtime)
                        else:
                            raise FileNotFoundError(
                                "Surface muon file not found — run the Generator tab first."
                            )

                        # Ensure direction cosines are present
                        if "cx" not in df_surf.columns:
                            if "px" in df_surf.columns and "p" in df_surf.columns:
                                _p = df_surf["p"].values.astype(float)
                                _p = np.where(_p > 0, _p, 1.0)
                                df_surf["cx"] = df_surf["px"].values / _p
                                df_surf["cy"] = df_surf["py"].values / _p
                                df_surf["cz"] = df_surf["pz"].values / _p
                            else:
                                df_surf["cx"] = np.sin(df_surf["theta"]) * np.cos(df_surf["phi"])
                                df_surf["cy"] = np.sin(df_surf["theta"]) * np.sin(df_surf["phi"])
                                df_surf["cz"] = -np.cos(df_surf["theta"])

                        prog.progress(0.15, text="⏳  Transporting muons through terrain…")

                        # DEM transport engines — UCMuon-MC / MUSIC / BB / PROPOSAL
                        if terrain_engine == "UCMuon Stochastic":
                            stochastic_drv_path = Path(script_dir) / "ucmuon_stochastic_driver.py"
                            if not stochastic_drv_path.exists():
                                raise FileNotFoundError(
                                    "`ucmuon_stochastic_driver.py` not found in the gui/ directory.  "
                                    "Place it there and restart the app."
                                )
                            df_ug = transport_stochastic_terrain(
                                df_surf, ob_map, az_c, ze_c, sky_map, terrain_rho,
                                n_az, n_ze, ze_max, t_vcut, t_ms, script_dir,
                                progress_container=prog,
                            )
                        else:
                            # Map engine display name to internal name for batched transport
                            _eng_map = {
                                "MUSIC": "MUSIC",
                                "Bethe-Bloch": "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS",
                                "PROPOSAL": "PROPOSAL",
                            }
                            engine_cfg = {
                                "rho": terrain_rho, "rad": 26.48,
                                "idim": 1, "idim1": 1, "init": 1, "minv": -30,
                                "mat_id": 1, "source_mode": 1,
                                "plane_lx": 0.0, "plane_ly": 0.0,
                                "ncols": ncols_t,
                                "proposal_medium_type": 1, "proposal_e_cut": 500.0,
                                "proposal_v_cut": 0.001, "proposal_scatter": 2,
                                "proposal_tables_dir": "",
                                "phitsxs_mat_type": 1, "phitsxs_ms_enable": 1,
                            }
                            df_ug = transport_batched_terrain(
                                df_surf, ob_map, az_c, ze_c, sky_map, terrain_rho,
                                n_az, n_ze, ze_max, _eng_map.get(terrain_engine, terrain_engine),
                                engine_cfg, script_dir, project_dir,
                                build_music_input_fn, build_phitsxs_input_fn,
                                build_proposal_input_fn,
                                progress_container=prog,
                            )
                    # ── end if csg_geom is None (DEM transport) ──────────────────
    
                    elapsed = time.time() - t0
    
                    # Write output
                    if len(df_ug) > 0:
                        _cols18 = "EventID,xs,ys,zs,Es,theta_s,phi_s,charge,alive,x,y,z,E,cx,cy,cz,theta,phi".split(",")
                        _cols_out = [c for c in _cols18 if c in df_ug.columns]
                        df_ug[_cols_out].to_csv(t_outfile, sep=" ", index=False,
                                                float_format="%.6g", header=False)
    
                    n_surv  = int((df_ug["alive"] == 1).sum()) if "alive" in df_ug.columns else 0
                    n_total_t = len(df_ug)
                    prog.progress(1.0, text=f"✅  Done: {n_surv:,} survived / {n_total_t:,} transported")
    
                    st.session_state["terrain_result_df"]    = df_ug
                    st.session_state["terrain_result_ob"]    = ob_map
                    st.session_state["terrain_result_sky"]   = sky_map
                    st.session_state["terrain_result_az"]    = az_c
                    st.session_state["terrain_result_ze"]    = ze_c
                    st.session_state["terrain_result_outfile"]   = t_outfile
                    st.session_state["terrain_result_engine"]    = terrain_engine
                    st.session_state["terrain_result_underground"] = underground
                    st.session_state["terrain_result_det_cell"]  = det_cell_id  # None for DEM
                    st.session_state["ug_file"] = t_outfile
                    # Also expose to Tab 3 results
                    st.session_state["music_nmuons_transported"] = n_total_t
                    st.session_state["music_nmuons_survived"]    = n_surv
                    # Persist DEM/CSG arrays for results terrain viewer
                    if csg_geom is not None:
                        st.session_state["terrain_result_csg"]  = csg_geom
                        st.session_state.pop("terrain_result_elev", None)
                        st.session_state.pop("terrain_result_tfm",  None)
                    elif synth_dem is not None:
                        st.session_state["terrain_result_elev"]  = synth_dem.elev
                        st.session_state["terrain_result_tfm"]   = synth_dem.transform
                        st.session_state.pop("terrain_result_csg", None)
                        st.session_state.pop("terrain_dem_path", None)
                    else:
                        st.session_state.pop("terrain_result_elev", None)
                        st.session_state.pop("terrain_result_tfm",  None)
                        st.session_state.pop("terrain_result_csg",  None)

                    # Success message — CSG shows hit/survival breakdown
                    if csg_geom is not None and det_cell_id is not None:
                        n_reached = int((df_ug["x"] != df_ug["xs"]).sum())
                        st.success(
                            f"✅  CSG transport complete in **{elapsed:.1f} s**  |  "
                            f"**{n_reached:,}** reached detector cell {det_cell_id}  |  "
                            f"**{n_surv:,} / {n_total_t:,}** survived  "
                            f"({100*n_surv/max(n_total_t,1):.2f}%)  |  "
                            f"Output: `{t_outfile}`"
                        )
                    else:
                        st.success(
                            f"✅  Terrain transport complete in **{elapsed:.1f} s**  |  "
                            f"**{n_surv:,} / {n_total_t:,}** muons survived  "
                            f"({100*n_surv/max(n_total_t,1):.2f}%)  |  "
                            f"Output: `{t_outfile}`"
                        )
    
                except Exception as _e:
                    import traceback
                    st.error(f"❌  Transport failed: {_e}")
                    st.code(traceback.format_exc())
    
        # ── SECTION 7: Results ───────────────────────────────────────────────
        if "terrain_result_df" not in st.session_state:
            st.divider()
            st.info(
                "▶  Configure the 📋 Setup tab, then click **Run terrain transport** to "
                "compute per-direction survival rates and flux maps.  "
                "Results will appear here after the run completes.",
                icon="ℹ️"
            )
        else:
            df_ug  = st.session_state["terrain_result_df"]
            ob_map = st.session_state["terrain_result_ob"]
            sky_m  = st.session_state["terrain_result_sky"]
            az_c   = st.session_state["terrain_result_az"]
            ze_c   = st.session_state["terrain_result_ze"]

            if len(df_ug) == 0:
                st.warning("No output muons — check the input file and settings.")
            else:
                st.divider()
                st.markdown("### Results")

                n_surv  = int((df_ug["alive"] == 1).sum()) if "alive" in df_ug.columns else 0
                n_tot   = len(df_ug)
                _det_cell_res = st.session_state.get("terrain_result_det_cell")

                # ── CSG mode: show detector-cell specific metrics + plots ─────────────────
                if _det_cell_res is not None:
                    n_reached = int((df_ug["x"] != df_ug["xs"]).sum())
                    _m1, _m2, _m3, _m4 = st.columns(4)
                    _m1.metric("Surface muons",    f"{n_tot:,}")
                    _m2.metric("Reached cell",     f"{n_reached:,}",
                               delta=f"{100*n_reached/max(n_tot,1):.1f}%")
                    _m3.metric("Survived",         f"{n_surv:,}",
                               delta=f"{100*n_surv/max(n_reached,1):.1f}% of reached")
                    _m4.metric("Detector cell",    f"Cell {_det_cell_res}")
                    if _CSG_TRANSPORT_OK and plot_csg_transport_results is not None:
                        plot_csg_transport_results(df_ug, _det_cell_res)
                    st.divider()

                else:
                    # ── DEM mode: original metrics row ────────────────────────────────────
                    _m1, _m2, _m3, _m4 = st.columns(4)
                    _m1.metric("Transported",   f"{n_tot:,}")
                    _m2.metric("Survived",      f"{n_surv:,}")
                    _m3.metric("Survival rate", f"{100*n_surv/max(n_tot,1):.2f}%")
                    _m4.metric("Engine",        st.session_state.get("terrain_result_engine", "—"))
                    if n_surv == n_tot and n_tot > 0:
                        _ob_stored   = st.session_state.get("terrain_result_ob")
                        _ob_all_zero = (_ob_stored is not None and float(np.max(_ob_stored)) < 1.0)
                        if _ob_all_zero:
                            st.warning(
                                "Survival rate is 100% AND the overburden map is all zero — "
                                "no terrain blocking was detected. "
                                "Set Max zenith = 85° and recompute.",
                                icon="⚠️"
                            )
                        else:
                            st.info(
                                "ℹ️  99.9%+ survival is **physically correct** for this scenario. "
                                "The cos²θ distribution strongly concentrates muons near vertical (ze < 45°), "
                                "where the terrain is thin or absent. Only ~1–2% of muons travel in the "
                                "near-horizontal directions (ze > 74°) where blocking occurs.",
                                icon="ℹ️"
                            )
    
                # ── Compute per-direction survival and transmission grids (vectorised) ──────
                # Both surv_map/total_map and input_map/surv_map2 use the same bin assignment.
                # Call assign_direction_bins ONCE, use the result for both.
                # Replace pure-Python for-loops (O(N) CPython, ~60 s at 1 M muons) with
                # numpy np.add.at (vectorised, ~0.3 s).
                _n_az_v, _n_ze_v = len(az_c), len(ze_c)
                # Bin edges must match the overburden/transport grid exactly:
                # ze_c are bin CENTRES of linspace(0, ze_max, n_ze+1) edges, so
                # the grid edge is centre + half-step (85.0 for the defaults).
                _ze_step_v = (float(ze_c[1] - ze_c[0]) if _n_ze_v > 1
                              else 2.0 * float(ze_c[0]))
                _ze_max_v  = float(ze_c[-1]) + 0.5 * _ze_step_v
                # Bin by the SURFACE direction (theta_s/phi_s): it exists for
                # both alive and stopped rows — several engines reset a stopped
                # muon's exit direction to (0,0,-1), which would pile every
                # stopped muon into the vertical bin — and it matches the
                # binning used to assign per-direction depths during transport.
                if "theta_s" in df_ug.columns and "phi_s" in df_ug.columns:
                    import pandas as _pd_bins
                    _th_s = df_ug["theta_s"].values.astype(float)
                    _ph_s = df_ug["phi_s"].values.astype(float)
                    _df_bin_src = _pd_bins.DataFrame({
                        "cx":  np.sin(_th_s) * np.cos(_ph_s),
                        "cy":  np.sin(_th_s) * np.sin(_ph_s),
                        "cz": -np.cos(_th_s),
                    })
                else:
                    _df_bin_src = df_ug
                _az_idx, _ze_idx = assign_direction_bins(
                    _df_bin_src, _n_az_v, _n_ze_v, _ze_max_v)
    
                # Valid mask: bin indices are in range
                _valid_mask = ((_az_idx >= 0) & (_az_idx < _n_az_v) &
                               (_ze_idx >= 0) & (_ze_idx < _n_ze_v))
    
                # Flatten 2-D bin index → 1-D for np.bincount / np.add.at
                _flat_idx = _az_idx * _n_ze_v + _ze_idx   # shape (N,)
                _flat_idx_v = _flat_idx[_valid_mask]        # only valid bins
    
                total_map = np.bincount(_flat_idx_v, minlength=_n_az_v * _n_ze_v
                                        ).reshape(_n_az_v, _n_ze_v).astype(float)
    
                surv_map  = np.zeros(_n_az_v * _n_ze_v)
                if "alive" in df_ug.columns:
                    _alive_v  = df_ug["alive"].values[_valid_mask].astype(float)
                    np.add.at(surv_map, _flat_idx_v, _alive_v)
                surv_map  = surv_map.reshape(_n_az_v, _n_ze_v)
    
                # Transmission = surv_map / total_map × 100%  (NaN where no muons)
                # np.where evaluates the division for empty bins too — silence it.
                with np.errstate(divide='ignore', invalid='ignore'):
                    transmission_map = np.where(total_map > 0, surv_map / total_map * 100.0, np.nan)
                rate_map         = transmission_map   # same quantity — reuse variable name
                _trans_valid     = ~np.isnan(transmission_map) & (total_map > 0)
                input_map        = total_map           # alias — used by energy spectrum tab
                surv_map2        = surv_map            # alias — used by energy spectrum tab
                _az_idx_r        = _az_idx             # alias — used by energy spectrum tab
                _ze_idx_r        = _ze_idx             # alias

                # Pre-init radio key so the first click on "🏔️ 3D Terrain" tab
                # doesn't trigger a state-initializing re-run that resets the tab.
                if "terrain_3d_view_mode" not in st.session_state:
                    st.session_state["terrain_3d_view_mode"] = "🗺️ Top-down map"

                # Keyed segmented control instead of st.tabs: st.tabs keeps its
                # selection client-side only, so the rerun triggered by any
                # widget inside a view (3D radio, radius slider, …) snapped the
                # strip back to the first view.  The key persists the selection.
                _RESULT_VIEWS = [
                    "🗺️ Overburden (polar)",
                    "📐 Muogram (az×el)",
                    "📡 Survival rate",
                    "📈 Transmission",
                    "📊 Energy spectrum",
                    "🏔️ 3D Terrain",
                    "🔬 Literature cross-check",
                ]
                if st.session_state.get("terrain_results_view") not in _RESULT_VIEWS:
                    st.session_state["terrain_results_view"] = _RESULT_VIEWS[0]
                _res_view = st.segmented_control(
                    "Results view", _RESULT_VIEWS,
                    key="terrain_results_view",
                    label_visibility="collapsed",
                ) or _RESULT_VIEWS[0]   # clicking the active segment deselects → fall back
    
                if _res_view == _RESULT_VIEWS[0]:
                    st.caption(
                        "Skymap view: **radial axis = zenith angle** (0° = directly overhead, 90° = horizon).  "
                        "**Angular axis = geographic azimuth** (N=top, E=right, clockwise).  "
                        "Colour = log₁₀(overburden g/cm²). Gray = open sky.  "
                        "This format is the same as MURAVES / CCS monitoring skymaps (cf. screenshot)."
                    )
                    _ob_fig_sky = _skymap_polar(
                        az_c, ze_c, ob_map,
                        title="Rock overburden skymap — log₁₀ [g/cm²]",
                        unit="g/cm²", colorscale="Inferno",
                        mask=sky_m, log_scale=True,
                    )
                    st.plotly_chart(_ob_fig_sky,                                     config={"displayModeBar": True}, key="terrain_ob_result")
                    st.caption(
                        "How to read: the centre (ze=0°) is the direction directly overhead (sky above detector).  "
                        "Moving radially outward approaches the horizon (ze→85°).  "
                        "Coloured patches near the edge = rock directions at low elevation = the mountain mass.  "
                        "Compare directly to MURAVES Fig. 1 or CCS monitoring skymaps."
                    )
    
                if _res_view == _RESULT_VIEWS[1]:
                    st.caption(
                        "Standard muography display.  "
                        "X = azimuth (N=0, E=90), Y = elevation = 90 minus zenith.  "
                        "Log10(overburden) in colour.  Hover any bin for exact value."
                    )
                    _ob_max_val = float(np.nanmax(ob_map[ob_map > 0])) if np.any(ob_map > 0) else 1e5
                    _thresh = st.slider(
                        "Hide bins with overburden below [g/cm²]", 0, int(_ob_max_val),
                        100, 50, key="terrain_mug_thresh",
                        help=(
                            "Increase to ~1000 to highlight only significant terrain blocking "
                            "and suppress the diffuse background of small volcanic cones."
                        )
                    )
                    _ob_thresh = ob_map.copy().astype(float)
                    _ob_thresh[_ob_thresh < _thresh] = 0.0
                    _mug_mask = sky_m | (_ob_thresh <= 0)
                    _n_sig_bins = int(np.sum(_ob_thresh > 0))
                    _fig_mug = _muogram_heatmap(
                        az_c, ze_c, _ob_thresh,
                        title="Muography overburden — az x elevation  (log10 g/cm2)",
                        unit="g/cm2", colorscale="Jet",
                        log_scale=True,
                        mask=_mug_mask,
                    )
                    st.plotly_chart(_fig_mug,                                     config={"displayModeBar": True}, key="terrain_muogram_ob")
                    _mm1, _mm2 = st.columns(2)
                    _mm1.metric("Visible blocking bins", str(_n_sig_bins))
                    _mm2.metric("Max overburden", f"{_ob_max_val:.0f} g/cm2")
                    st.info(
                        "Red patch (az=0, el=7-16 deg): Puy de Dome lava dome - the main muographic target.  "
                        "Other coloured patches: secondary volcanic cones of the Chaine des Puys "
                        "(80+ cones in this volcanic chain - all visible when threshold is low).  "
                        "Increase the threshold slider to 1000-10000 g/cm2 to isolate only "
                        "the dominant geological signal from background terrain.",
                        icon="🌋"
                    )
    
                if _res_view == _RESULT_VIEWS[2]:
                    _valid = ~np.isnan(rate_map) & (rate_map > 0)
                    if _valid.any():
                        _c1, _c2 = st.columns(2)
                        with _c1:
                            st.caption("Skymap view — survival rate per direction bin.")
                            st.plotly_chart(
                                _skymap_polar(az_c, ze_c, np.nan_to_num(rate_map),
                                              "Muon survival rate skymap [%]", "%",
                                              colorscale="RdYlGn", mask=~_valid,
                                              log_scale=False),
                                config={"displayModeBar": False},
                                key="terrain_surv_result"
                            )
                        with _c2:
                            st.caption("Rectangular view — same data in az × elevation.")
                            st.plotly_chart(
                                _muogram_heatmap(az_c, ze_c, np.nan_to_num(rate_map),
                                                 "Survival rate [%]", "%",
                                                 colorscale="RdYlGn", mask=~_valid,
                                                 log_scale=False),
                                config={"displayModeBar": False},
                                key="terrain_surv_muogram"
                            )
                        st.caption(
                            "Green = high survival (open sky, thin overburden).  "
                            "Red = low survival (muons stopped by the summit).  "
                            "Gray = empty bins (no muons in that direction).  "
                            "The red patch in the north at low elevation IS the muographic signal."
                        )
                    else:
                        st.info("No directional survival data — run the transport first.")
    
                if _res_view == _RESULT_VIEWS[3]:
                    st.caption(
                        "**Transmission map** T(az, el) = survived/input per bin.  "
                        "This normalises out the cos²θ flux weighting and shows **only the terrain effect**.  "
                        "T=100% = open sky (no blocking). T→0% = fully blocked by rock."
                    )
                    if _trans_valid.any():
                        _fig_trans = _muogram_heatmap(
                            az_c, ze_c, transmission_map,
                            title="Transmission map T = survived / input  [%]",
                            unit="%", colorscale="RdYlGn",
                            zmin=0, zmax=100,
                            log_scale=False,
                            mask=~_trans_valid,
                        )
                        st.plotly_chart(_fig_trans,                                         config={"displayModeBar": True}, key="terrain_transmission")
                        # Stats on blocked bins
                        _t_blocked = transmission_map[_trans_valid & (transmission_map < 99.0)]
                        if len(_t_blocked) > 0:
                            _tm1, _tm2, _tm3 = st.columns(3)
                            _tm1.metric("Blocked bins (T<99%)", f"{len(_t_blocked)}")
                            _tm2.metric("Min transmission", f"{float(_t_blocked.min()):.2f}%")
                            _tm3.metric("Mean transmission (blocked)", f"{float(_t_blocked.mean()):.1f}%")
                        st.caption(
                            "T=0% = fully blocked by rock.  T=100% = open sky.  "
                            "Main signal (az near 0, el 7-16 deg): Puy de Dome lava dome.  "
                            "Secondary patches (az near 150 deg): other Chaine des Puys cones - "
                            "real geology, not simulation artefacts.  "
                            "The ratio T_blocked/T_open gives the muographic opacity per direction."
                        )
                    else:
                        st.info("Not enough muons per bin to compute transmission map. "
                                "Increase N muons to ≥500,000.")

                    # ── Save T_sim for density inversion ──────────────────────────────
                    if _trans_valid.any():
                        st.divider()
                        st.markdown("#### 💾 Save T_sim for density inversion")
                        st.caption(
                            "Run terrain transport at **several different densities** (e.g. 1.5, 2.0, 2.5, 3.0 g/cm³), "
                            "saving one file per density.  Load all files in the **Density Analysis** tab as the T_sim library."
                        )
                        _rho_curr_ts  = _sf(st.session_state.get("terrain_rho", 2.65), 2.65)
                        _tsim_c1, _tsim_c2 = st.columns([3, 1])
                        _tsim_fname = _tsim_c1.text_input(
                            "Output filename",
                            value=f"terrain_transmission_{_rho_curr_ts:.2f}.dat",
                            key="terrain_tsim_outfile",
                        )
                        if _tsim_c2.button("💾 Save", key="terrain_save_tsim", type="primary"):
                            try:
                                _T_norm = np.where(
                                    total_map > 0,
                                    surv_map / np.where(total_map > 0, total_map, 1.0),
                                    np.nan,
                                )
                                _lat_ts = _sf(st.session_state.get("terrain_lat", 0.0), 0.0)
                                _lon_ts = _sf(st.session_state.get("terrain_lon", 0.0), 0.0)
                                _alt_ts = _sf(st.session_state.get("terrain_alt", 0.0), 0.0)
                                _drv_ts = _load_terrain_driver(script_dir)
                                # Bins with no muons stay NaN ("no data"):
                                # writing 0 would make the density inversion
                                # read them as fully blocked rock (rho_max).
                                _drv_ts.write_transmission_map(
                                    az_c, ze_c, _T_norm,
                                    _tsim_fname,
                                    _lat_ts, _lon_ts, _alt_ts, _rho_curr_ts,
                                )
                                _pending = list(st.session_state.get("da_lib_paths_pending", []))
                                if _tsim_fname not in _pending:
                                    _pending.append(_tsim_fname)
                                st.session_state["da_lib_paths_pending"] = _pending
                                st.success(
                                    f"✅ Saved `{_tsim_fname}` (ρ = {_rho_curr_ts:.2f} g/cm³).  "
                                    "Now change density in Section 2 and run again.  "
                                    "Once you have ≥3 files, go to the **Density Analysis** tab."
                                )
                            except Exception as _tsim_err:
                                st.error(f"Failed to write T_sim file: {_tsim_err}")

                if _res_view == _RESULT_VIEWS[4]:
                    if "E" in df_ug.columns and "alive" in df_ug.columns:
                        _e_surv = df_ug.loc[df_ug["alive"] == 1, "E"].values
                        _e_surf = df_ug["Es"].values if "Es" in df_ug.columns else None
                        fig_e = go.Figure()
                        if _e_surf is not None:
                            fig_e.add_trace(go.Histogram(
                                x=np.log10(_e_surf[_e_surf > 0]), nbinsx=50, name="Surface (input)",
                                marker_color="rgba(100,200,255,0.5)", opacity=0.8,
                            ))
                        if len(_e_surv) > 0:
                            fig_e.add_trace(go.Histogram(
                                x=np.log10(_e_surv[_e_surv > 0]), nbinsx=50,
                                name="Survived (all directions)",
                                marker_color="rgba(255,120,80,0.7)", opacity=0.8,
                            ))
                        # Survived from blocked directions only — vectorised
                        _blocked_mask = np.zeros(len(df_ug), dtype=bool)
                        if "alive" in df_ug.columns:
                            # Use bin indices already computed; ob_map[ia,iz] > 1000 g/cm²
                            # flags significantly blocked bins (> ~6 m at std rock density)
                            _valid_en = ((_az_idx_r >= 0) & (_az_idx_r < len(az_c)) &
                                         (_ze_idx_r >= 0) & (_ze_idx_r < len(ze_c)))
                            _ia_en = _az_idx_r[_valid_en].astype(int)
                            _iz_en = _ze_idx_r[_valid_en].astype(int)
                            _blocked_mask[_valid_en] = ob_map[_ia_en, _iz_en] > 1000.0
                        _e_bl = df_ug.loc[_blocked_mask & (df_ug["alive"] == 1), "E"].values
                        if len(_e_bl) > 0:
                            fig_e.add_trace(go.Histogram(
                                x=np.log10(_e_bl[_e_bl > 0]), nbinsx=40,
                                name="Survived through terrain (blocked dirs)",
                                marker_color="rgba(80,255,80,0.8)", opacity=0.9,
                            ))
                        _log_y = st.checkbox("Log Y axis (recommended — reveals green bars)", value=True, key="terrain_energy_logy")
                        fig_e.update_layout(
                            **DARK, height=380, barmode="overlay",
                            xaxis=dict(title="log₁₀(E [GeV])", gridcolor="#2a2a3a"),
                            yaxis=dict(title="Counts", gridcolor="#2a2a3a",
                                       type="log" if _log_y else "linear"),
                            legend=dict(bgcolor="rgba(0,0,0,0.5)", font=dict(size=10)),
                            margin=dict(l=60, r=20, t=30, b=50),
                        )
                        st.plotly_chart(fig_e,                                         config={"displayModeBar": False}, key="terrain_energy_result")
                        _n_tr = int(len(_e_bl)) if len(_e_bl) > 0 else 0
                        st.caption(
                            f"Blue = surface spectrum.  Orange = all survived.  "
                            f"Green = survived through terrain ({_n_tr:,} muons, blocked dirs only).  "
                            "Green is nearly invisible on linear scale: only ~5% of CosmoALEPH muons "
                            "exceed the ~400 GeV CSDA threshold for 1293 m of rock.  "
                            "Enable log Y to reveal the hardened green spectrum."
                        )
    
                if _res_view == _RESULT_VIEWS[5]:
                    _dem_path_3d  = st.session_state.get("terrain_dem_path", "")
                    _synth_elev   = st.session_state.get("terrain_result_elev")
                    _synth_tfm    = st.session_state.get("terrain_result_tfm")
                    _csg_geom_3d  = st.session_state.get("terrain_result_csg")   # CSGGeometry or None
                    _det_lat_3d   = _sf(st.session_state.get("terrain_lat", 45.76), 45.76)
                    _det_lon_3d   = _sf(st.session_state.get("terrain_lon", 2.955), 2.955)
                    _det_alt_3d   = _sf(st.session_state.get("terrain_alt", 1094.0), 1094.0)
                    _rho_t3       = _sf(st.session_state.get("terrain_rho", 2.65), 2.65)
                    _ug_3d        = bool(st.session_state.get("terrain_result_underground", False))

                    # Four cases:
                    #   A) CSG geometry in session_state  → use csg_geom.plotly_preview()
                    #   B) synth arrays in session_state  → use DEM surface plot (synthetic)
                    #   C) file path exists on disk       → use file DEM plot
                    #   D) none                           → info message
                    _using_csg_3d   = (_csg_geom_3d is not None)
                    _using_synth_3d = (not _using_csg_3d and _synth_elev is not None and _synth_tfm is not None)
                    _dem_ok_3d      = _using_csg_3d or _using_synth_3d or \
                                      (bool(_dem_path_3d) and Path(_dem_path_3d).exists())

                    if _ug_3d:
                        st.caption("🔽 Underground detector mode was active for this run.")

                    if not _dem_ok_3d:
                        if _dem_path_3d:
                            st.warning(
                                f"⚠️  DEM file `{Path(_dem_path_3d).name}` no longer exists on disk.  "
                                "Re-upload the DEM in Section 2 to restore the terrain views.",
                                icon="⚠️"
                            )
                        else:
                            st.info(
                                "📂  Load a DEM file, **Synthetic Geometry**, or **CSG Geometry** "
                                "in Section 2 above to enable the terrain visualisation.  "
                                "The overburden polar map (first tab) and flux plots work without a DEM.",
                                icon="📂"
                            )

                    elif _using_csg_3d:
                        # ── CSG 3D preview ─────────────────────────────────────────────
                        st.caption("🔩 Showing CSG geometry (from Section 2)")
                        _csg_dp_3d = st.session_state.get("_csg_det_pos", np.array([0.,0.,0.]))
                        try:
                            _fig_csg = _csg_geom_3d.plotly_preview(_csg_dp_3d)
                            st.plotly_chart(_fig_csg,                                             config={"displayModeBar": True}, key="terrain_csg_3d")
                        except Exception as _csg_e:
                            st.warning(f"CSG 3D preview failed: {_csg_e}")

                    else:
                        if _using_synth_3d:
                            st.caption("🔷 Showing synthetic geometry DEM (from Section 2)")

                        # ── View mode radio ───────────────────────────────────────────
                        _view_mode = st.radio(
                            "Terrain view", ["🗺️ Top-down map", "🏔️ 3D terrain"],
                            horizontal=True, key="terrain_3d_view_mode", label_visibility="collapsed"
                        )

                        # ── Top-down 2D map ────────────────────────────────────────────
                        if _view_mode == "🗺️ Top-down map":
                            st.caption(
                                "Plan-view map centred on the detector (gold ★ at origin).  "
                                "Background = terrain elevation contour.  "
                                "Each coloured line = a blocked muon direction; length ∝ "
                                "horizontal distance to terrain endpoint; colour = overburden "
                                "(blue = shallow, red = deepest).  Hover any dot for exact values."
                            )
                            _r2d = st.slider("Map radius [°]", 0.03, 0.40, 0.15, 0.01,
                                             key="terrain_2d_radius")
                            _fig2d, _err2d = _dem_2d_topdown(
                                _dem_path_3d, _det_lat_3d, _det_lon_3d, _det_alt_3d,
                                az_c, ze_c, ob_map, sky_m, _rho_t3,
                                radius_deg=_r2d,
                                elev_arr=_synth_elev, tfm_arr=_synth_tfm,
                            )
                            if _err2d:
                                st.error(f"2D map error: {_err2d}", icon="❌")
                            elif _fig2d is not None:
                                st.plotly_chart(_fig2d,                                                 config={"displayModeBar": True},
                                                key="terrain_2d_map")
                                _n_bl_2d   = int((~sky_m).sum())
                                _ob_max_2d = float(ob_map[~sky_m].max()) if _n_bl_2d > 0 else 0
                                _slant_2d  = _ob_max_2d / (_rho_t3 * 100.0) if _rho_t3 > 0 else 0
                                st.caption(
                                    f"{_n_bl_2d} blocked direction bins.  "
                                    f"Longest path: {_slant_2d:.0f} m slant  "
                                    f"({_ob_max_2d:,.0f} g/cm²)."
                                )

                        # ── 3D terrain ────────────────────────────────────────────────
                        else:
                            st.caption(
                                "Interactive 3D terrain (ENU, detector at origin).  "
                                "🟡 Gold diamond = detector.  🔴 Red circle = DEM summit.  "
                                "Coloured lines = blocked ray directions (blue→red = shallow→deep).  "
                                "Rotate / zoom with mouse."
                            )
                            _r_deg = st.slider("View radius [°]", 0.02, 0.30, 0.08, 0.01,
                                               key="terrain_3d_radius")
                            _fig_3d = _dem_3d_plot(
                                _dem_path_3d, _det_lat_3d, _det_lon_3d,
                                _det_alt_3d, script_dir, radius_deg=_r_deg,
                                elev_arr=_synth_elev, tfm_arr=_synth_tfm,
                            )
                            if isinstance(_fig_3d, tuple) and _fig_3d[0] == "ERROR":
                                st.error(
                                    f"❌ 3D plot failed: {_fig_3d[1]}  \n"
                                    "Try increasing the View radius slider.",
                                    icon="❌"
                                )
                                if _fig_3d[2]:
                                    with st.expander("Technical details"):
                                        st.code(_fig_3d[2])
                                _fig_3d = None
    
                            if _fig_3d is not None:
                                import plotly.graph_objects as _go3
                                _ob_block = ob_map[~sky_m]
                                if _ob_block.size > 0:
                                    _ob_min_log3 = float(np.log10(max(_ob_block.min(), 1.0)))
                                    _ob_max_log3 = float(np.log10(max(_ob_block.max(), 1.0)))
    
                                    def _jet3(t):
                                        t = float(np.clip(t, 0, 1))
                                        if   t < 0.25: r,g,b = 0,                   int(255*t*4),            255
                                        elif t < 0.50: r,g,b = 0,                   255,                     int(255*(1-(t-0.25)*4))
                                        elif t < 0.75: r,g,b = int(255*(t-0.5)*4),  255,                     0
                                        else:          r,g,b = 255,                  int(255*(1-(t-0.75)*4)), 0
                                        return f"rgb({r},{g},{b})"
    
                                    _ex_pts, _ny_pts, _uz_pts = [], [], []
                                    _ob_log_pts, _ob_hover    = [], []
                                    _legend_shown = False
                                    for _ia3, _az3v in enumerate(az_c):
                                        for _iz3, _ze3v in enumerate(ze_c):
                                            if sky_m[_ia3, _iz3] or ob_map[_ia3, _iz3] < 10:
                                                continue
                                            _ob3v = float(ob_map[_ia3, _iz3])
                                            _sl3  = _ob3v / (_rho_t3 * 100.0)
                                            _ze3r = np.radians(float(_ze3v))
                                            _az3r = np.radians(float(_az3v))
                                            _ex3  = _sl3 * np.sin(_ze3r) * np.sin(_az3r)
                                            _ny3  = _sl3 * np.sin(_ze3r) * np.cos(_az3r)
                                            _uz3  = _sl3 * np.cos(_ze3r)
                                            _t3   = (np.log10(_ob3v) - _ob_min_log3) / max(
                                                        _ob_max_log3 - _ob_min_log3, 0.1)
                                            _col3 = _jet3(_t3)
                                            # One Scatter3d per ray (uniform colour — None not valid in 3d)
                                            _fig_3d.add_trace(_go3.Scatter3d(
                                                x=[0, _ex3], y=[0, _ny3],
                                                z=[_det_alt_3d, _det_alt_3d + _uz3],
                                                mode="lines",
                                                line=dict(color=_col3, width=4),
                                                hoverinfo="skip",
                                                showlegend=(not _legend_shown),
                                                name="Blocked ray directions" if not _legend_shown else "",
                                                legendgroup="rays",
                                            ))
                                            _legend_shown = True
                                            _ex_pts.append(_ex3)
                                            _ny_pts.append(_ny3)
                                            _uz_pts.append(_det_alt_3d + _uz3)
                                            _ob_log_pts.append(float(np.log10(_ob3v)))
                                            _ob_hover.append(
                                                f"az={_az3v:.0f}°  ze={_ze3v:.0f}°<br>"
                                                f"overburden = {_ob3v:,.0f} g/cm²<br>"
                                                f"slant path = {_sl3:.0f} m"
                                            )
                                    if _ex_pts:
                                        _fig_3d.add_trace(_go3.Scatter3d(
                                            x=_ex_pts, y=_ny_pts, z=_uz_pts,
                                            mode="markers",
                                            marker=dict(
                                                size=6, color=_ob_log_pts,
                                                colorscale="Jet",
                                                cmin=_ob_min_log3, cmax=_ob_max_log3,
                                                colorbar=dict(
                                                    title=dict(text="log₁₀(ob)<br>g/cm²",
                                                               font=dict(color="white", size=10)),
                                                    tickfont=dict(color="white", size=9),
                                                    x=1.02, thickness=12, len=0.50,
                                                    bgcolor="rgba(20,22,30,0.8)",
                                                    tickvals=np.linspace(
                                                        _ob_min_log3, _ob_max_log3, 4).tolist(),
                                                    ticktext=[
                                                        f"{10**v:,.0f}" for v in
                                                        np.linspace(_ob_min_log3, _ob_max_log3, 4)
                                                    ],
                                                ),
                                                line=dict(color="white", width=1),
                                            ),
                                            text=_ob_hover,
                                            hovertemplate="%{text}<extra></extra>",
                                            name="Terrain endpoints",
                                        ))
                                st.plotly_chart(_fig_3d,                                                 config={"displayModeBar": True},
                                                key="terrain_3d_dem")
                                st.caption(
                                    "Lines = ray directions from detector; endpoint = where ray "
                                    "exits terrain.  Blue = shallow path, Red = deepest path."
                                )
                            elif not isinstance(_fig_3d, tuple):
                                st.warning(
                                    "3D plot could not be rendered.  "
                                    "Check that rasterio is installed (`pip install rasterio`).",
                                    icon="⚠️"
                                )
    
                    st.divider()
                    with st.expander("📚  Physics validation — reference deployments",
                                     expanded=False):
                        st.markdown("""
            **How to validate your muogram against published results**
    
            | Site | Reference | Detector position | Notes |
            |------|-----------|-------------------|-------|
            | Puy de Dôme | Carloganu et al. 2013, GIDS **2**, 55 | 45.742°N 2.955°E 870 m | 2° bins, ~350,000 g/cm² max |
            | Mt. Vesuvius | Hong et al. 2025, J. Appl. Phys. **138**, doi:10.1063/5.0275078 | 40.827°N 14.401°E 608 m | MURAVES scintillator tracker, E_min ≈ 1 GeV |
    
            **Angular resolution:** to resolve a structure of angular width θ, use bin size ≤ θ/3.
            A 600 m dome at 2 km subtends ≈ 17° → use ≤ 5° bins (180 az × 45 ze) for shape imaging.
    
            **Minimum statistics for 10% precision at transmission T:**
            N ≥ 100/T muons per bin.  At T = 1%: N ≥ 10,000 per bin → ≥ 10 M muons total.
    
            **Low-threshold detectors (nuclear emulsion, RPC):**
            Use Guan et al. 2015 spectrum with E_min = 1 GeV.
            CosmoALEPH with E_min = 100 GeV cannot image paths with X < ~50,000 g/cm².
            """)
                        _ob_rock_v = ob_map[~sky_m]
                        if _ob_rock_v.size > 0:
                            _n_bl_v    = int((~sky_m).sum())
                            _ob_med_v  = float(np.median(_ob_rock_v))
                            _ob_max_v2 = float(np.max(_ob_rock_v))
                            _rho_vv    = _sf(st.session_state.get("terrain_rho", 2.65), 2.65)
                            st.info(
                                f"**Current result:** {_n_bl_v} blocked bins  |  "
                                f"median {_ob_med_v:,.0f} g/cm²  |  "
                                f"max {_ob_max_v2:,.0f} g/cm²  "
                                f"({_ob_max_v2/(_rho_vv*100):.0f} m slant at ρ={_rho_vv:.2f} g/cm³).",
                                icon="✅"
                            )
    
                # ── View 7: Literature cross-check ────────────────────────────────────────
                if _res_view == _RESULT_VIEWS[6]:
                    st.markdown(
                        "**Literature cross-check.** Compare UCMuon's terrain overburden and "
                        "integrated flux against published muography references for Mt. Vesuvius. "
                        "UCMuon curves are computed from your DEM and detector position; external "
                        "overlays are drawn only from peer-reviewed sources and labelled with their "
                        "citation."
                    )
                    st.caption(
                        "Convention: the published MURAVES analysis (Hong et al. 2025, "
                        "J. Appl. Phys. 138) uses a detector at ~40.81 °N, 14.41 °E, 598 m a.s.l. "
                        "and a detector-centric azimuth with the summit near 180°. If your detector "
                        "differs, the thickness-profile *shape* should still agree, but the absolute "
                        "azimuth of the blocked region shifts accordingly. UCMuon uses geographic "
                        "azimuth (N = 0°, E = 90°)."
                    )

                    _rho_comp  = _sf(st.session_state.get("terrain_rho", 2.65), 2.65)
                    _det_lat_c = _sf(st.session_state.get("terrain_lat", 40.827), 40.827)
                    _det_lon_c = _sf(st.session_state.get("terrain_lon", 14.401), 14.401)
                    _det_alt_c = _sf(st.session_state.get("terrain_alt", 608.0), 608.0)

                    # Load backward MC module once
                    try:
                        import importlib as _il
                        _bmc_spec = _il.util.spec_from_file_location(
                            "ucmuon_backward_mc",
                            str(Path(script_dir) / "ucmuon_backward_mc.py")
                        )
                        _bmc = _il.util.module_from_spec(_bmc_spec)
                        _bmc_spec.loader.exec_module(_bmc)
                        _bmc_ok = True
                    except Exception as _bmc_err:
                        _bmc_ok = False
                        st.warning(f"ucmuon_backward_mc.py not found — flux plots unavailable: {_bmc_err}")

                    # ── 1 — Select target azimuth ─────────────────────────────────────────
                    st.markdown("#### 1 — Select target azimuth")
                    st.caption(
                        "Pick the azimuth pointing toward the geological target (e.g. the crater). "
                        "It defaults to the direction of maximum rock overburden. The 2D map below "
                        "marks your choice; the profiles further down are taken along it."
                    )
                    _col_sl1, _col_sl2 = st.columns([2, 1])
                    _az_target = _col_sl1.number_input(
                        "Target azimuth [° geographic, N=0°, E=90°]", 0.0, 360.0,
                        float(az_c[np.argmax(ob_map.max(axis=1))]),
                        1.0, key="comp_az_target",
                        help="The bin nearest to this value will be selected.",
                    )
                    _az_slice_idx = int(np.argmin(np.abs(az_c - _az_target)))
                    _az_slice_deg = float(az_c[_az_slice_idx])
                    _col_sl2.metric("Nearest bin centre", f"{_az_slice_deg:.1f}°")
                    _el_c = 90.0 - ze_c   # elevation centres

                    # ── 2 — Full rock-thickness map L(az, el) ─────────────────────────────
                    st.markdown("#### 2 — Rock thickness map  L(az, el)")
                    st.caption(
                        "The 2D slant-rock-path map — the overburden input a TURTLE/Mulder-style "
                        "backward transport would use. Colour = slant path [m]. Use it to choose the "
                        "target azimuth above; the gold dotted line marks the current selection."
                    )
                    _thick_map = ob_map / (_rho_comp * 100.0)   # m, shape (n_az, n_ze)
                    _thick_map_masked = _thick_map.copy().astype(float)
                    _thick_map_masked[sky_m] = np.nan
                    _el_centres = 90.0 - ze_c
                    _az_cent_shifted = (az_c + 180.0) % 360.0 - 180.0   # centre on N
                    _sort_az = np.argsort(_az_cent_shifted)
                    _Z_thick = _thick_map_masked[_sort_az, :].T   # (n_ze, n_az) for plotly heatmap
                    _fig_thmap = go.Figure(go.Heatmap(
                        x=_az_cent_shifted[_sort_az],
                        y=_el_centres,
                        z=_Z_thick,
                        colorscale="Plasma",
                        colorbar=dict(
                            title=dict(text="Rock thickness [m]", font=dict(color="white")),
                            tickfont=dict(color="white"),
                        ),
                        hoverongaps=False,
                        xgap=1, ygap=1,
                        hovertemplate="az=%{x:.1f}°  el=%{y:.1f}°  L=%{z:.0f} m<extra></extra>",
                    ))
                    _fig_thmap.add_vline(
                        x=(_az_slice_deg + 180.0) % 360.0 - 180.0,
                        line=dict(color="#ffd700", width=2, dash="dot"),
                        annotation_text=f"selected az={_az_slice_deg:.0f}°",
                        annotation_font=dict(color="#ffd700", size=10),
                    )
                    _fig_thmap.update_layout(
                        **DARK, height=380,
                        title=dict(text="Rock thickness map  L(az, el)  [m]", font=dict(size=12)),
                        xaxis=dict(title="Azimuth [°]  (N=0°, E=90°, W=−90°)", gridcolor="#2a2a3a",
                                   zeroline=True, zerolinecolor="#ffd700", zerolinewidth=1),
                        yaxis=dict(title="Elevation above horizon [°]", gridcolor="#2a2a3a"),
                        margin=dict(l=65, r=20, t=50, b=55),
                    )
                    st.plotly_chart(_fig_thmap, config={"displayModeBar": True}, key="comp_thickness_map")
                    _th_max_v = float(np.nanmax(_thick_map_masked)) if np.any(~np.isnan(_thick_map_masked)) else 0
                    _imax_all = np.unravel_index(np.argmax(ob_map), ob_map.shape)
                    st.caption(
                        f"Max thickness: **{_th_max_v:.0f} m** at az={az_c[_imax_all[0]]:.0f}°, "
                        f"el={90-ze_c[_imax_all[1]]:.0f}°.  "
                        "Compare with the published MURAVES thickness map for Mt. Vesuvius in "
                        "Hong et al. (2025), J. Appl. Phys. 138.  "
                        "Gray bins = open sky (no terrain intersection)."
                    )

                    # ── 3 — Rock thickness vs elevation at selected azimuth ───────────────
                    st.markdown("#### 3 — Rock thickness vs elevation at selected azimuth")
                    st.caption(
                        "Slant rock path [m] = overburden [g/cm²] / (ρ × 100), along the azimuth "
                        "selected above. Open-sky bins are excluded. Elevation = 90° − zenith."
                    )
                    _thick_slice = ob_map[_az_slice_idx, :] / (_rho_comp * 100.0)  # m
                    _sky_slice   = sky_m[_az_slice_idx, :]
                    _fig_thick = go.Figure()
                    _fig_thick.add_trace(go.Scatter(
                        x=_el_c[~_sky_slice],
                        y=_thick_slice[~_sky_slice],
                        mode="lines+markers",
                        name=f"UCMuon  az = {_az_slice_deg:.1f}°  ρ = {_rho_comp:.2f} g/cm³",
                        line=dict(color="#00b4d8", width=2.5),
                        marker=dict(size=7, color="#00b4d8"),
                        hovertemplate="el = %{x:.1f}°  thickness = %{y:.0f} m<extra></extra>",
                    ))
                    _fig_thick.update_layout(
                        **DARK, height=360,
                        title=dict(text=f"Rock thickness vs elevation — az = {_az_slice_deg:.1f}°",
                                   font=dict(size=12)),
                        xaxis=dict(title="Elevation above horizon [°]", gridcolor="#2a2a3a"),
                        yaxis=dict(title="Rock thickness [m]", gridcolor="#2a2a3a"),
                        legend=dict(bgcolor="rgba(0,0,0,0.5)", font=dict(size=10)),
                        margin=dict(l=70, r=20, t=50, b=55),
                    )
                    st.plotly_chart(_fig_thick, config={"displayModeBar": True}, key="comp_thick_profile")
                    st.caption(
                        f"UCMuon slant path along geographic az = {_az_slice_deg:.1f}° at detector "
                        f"{_det_lat_c:.4f} °N, {_det_lon_c:.4f} °E, {_det_alt_c:.0f} m.  "
                        "For a published MURAVES thickness/transmission comparison at Mt. Vesuvius see "
                        "Hong et al. (2025), J. Appl. Phys. 138, doi:10.1063/5.0275078 — note their "
                        "profile uses a 5 m LIDAR DEM and their own detector position and azimuth "
                        "convention, so absolute values and azimuth will differ from a 30 m SRTM run."
                    )

                    # ── 4 — Integrated muon flux vs elevation (backward MC) ────────────────
                    st.markdown("#### 4 — Integrated muon flux vs elevation at selected azimuth")
                    st.caption(
                        "Backward MC (Guan et al. 2015 CSDA + stochastic correction) run for each "
                        "elevation bin at the selected azimuth, for both terrain-blocked and "
                        "open-sky directions."
                    )
                    if _bmc_ok:
                        _spec_mode_comp = st.selectbox(
                            "Spectrum model for flux calculation",
                            options=[3, 4, 1, 2],
                            format_func=lambda x: {1:"CosmoALEPH", 2:"Power-law", 3:"Guan et al. 2015",
                                                   4:"Frosin et al. 2025"}[x],
                            index=0, key="comp_spec_mode",
                        )
                        _emin_comp = st.number_input("E_min [GeV]", 0.1, 100.0, 1.0, 0.1, key="comp_emin")
                        _emax_comp = st.number_input("E_max [GeV]", 10.0, 5000.0, 1000.0, 10.0, key="comp_emax")

                        if st.button("▶  Compute flux vs elevation (backward MC)", key="comp_flux_btn",
                                     width='stretch', type="primary"):
                            _flux_terrain  = np.full(len(ze_c), np.nan)
                            _flux_opensky  = np.full(len(ze_c), np.nan)
                            _prog_comp = st.progress(0.0, text="⏳  Running backward MC per elevation bin…")
                            for _iz_c, _ze_c_v in enumerate(ze_c):
                                _el_v = 90.0 - _ze_c_v
                                if _el_v <= 0:
                                    continue
                                # Open-sky flux: depth→0 (thin air layer)
                                try:
                                    _res_sky = _bmc.backward_mc_flux(
                                        depth_m=0.01, rho=1.0, mat_id=1,
                                        spectrum_mode=_spec_mode_comp,
                                        E_min_GeV=_emin_comp, E_max_GeV=_emax_comp,
                                        theta_max_deg=_ze_c_v, n_E=40, n_theta=1,
                                        mode=1,
                                    )
                                    _dze  = ze_c[1] - ze_c[0] if len(ze_c) > 1 else 5.0
                                    _daz  = 360.0 / len(az_c)
                                    _dO   = np.radians(_daz) * np.radians(_dze) * np.cos(np.radians(_ze_c_v))
                                    _flux_opensky[_iz_c] = _res_sky["rate_m2_s"] / max(_dO, 1e-6)
                                except Exception:
                                    pass
                                # Terrain flux: use slant overburden for this bin
                                if not _sky_slice[_iz_c] and ob_map[_az_slice_idx, _iz_c] > 1.0:
                                    _ob_sl = float(ob_map[_az_slice_idx, _iz_c])
                                    _depth = _ob_sl / (_rho_comp * 100.0) * np.cos(np.radians(_ze_c_v))
                                    _depth = max(_depth, 0.1)
                                    try:
                                        _res_ter = _bmc.backward_mc_flux(
                                            depth_m=_depth, rho=_rho_comp, mat_id=1,
                                            spectrum_mode=_spec_mode_comp,
                                            E_min_GeV=_emin_comp, E_max_GeV=_emax_comp,
                                            theta_max_deg=_ze_c_v, n_E=40, n_theta=1,
                                            mode=1,
                                        )
                                        _flux_terrain[_iz_c] = _res_ter["rate_m2_s"] / max(_dO, 1e-6)
                                    except Exception:
                                        pass
                                _prog_comp.progress((_iz_c + 1) / len(ze_c),
                                                    text=f"⏳  el = {_el_v:.1f}°  ({_iz_c+1}/{len(ze_c)})")
                            st.session_state["comp_flux_terrain"] = _flux_terrain
                            st.session_state["comp_flux_opensky"] = _flux_opensky
                            st.session_state["comp_flux_el"]      = _el_c
                            _prog_comp.progress(1.0, text="✅  Done.")

                        if "comp_flux_terrain" in st.session_state:
                            _ft = st.session_state["comp_flux_terrain"]
                            _fo = st.session_state["comp_flux_opensky"]
                            _el_ax = st.session_state["comp_flux_el"]
                            _log_flux = st.checkbox("Log scale (recommended)", value=True, key="comp_logscale")
                            _fig_flux = go.Figure()
                            _ok_sky = ~np.isnan(_fo) & (_fo > 0)
                            if _ok_sky.any():
                                _fig_flux.add_trace(go.Scatter(
                                    x=_el_ax[_ok_sky], y=_fo[_ok_sky],
                                    mode="lines+markers",
                                    name=f"UCMuon open-sky (Guan, E={_emin_comp:.1f}–{_emax_comp:.0f} GeV)",
                                    line=dict(color="#ffd700", width=2, dash="dash"),
                                    marker=dict(size=6, symbol="circle", color="#ffd700"),
                                ))
                            _ok_ter = ~np.isnan(_ft) & (_ft > 0)
                            if _ok_ter.any():
                                _fig_flux.add_trace(go.Scatter(
                                    x=_el_ax[_ok_ter], y=_ft[_ok_ter],
                                    mode="lines+markers",
                                    name=f"UCMuon terrain az={_az_slice_deg:.0f}° ρ={_rho_comp:.2f} g/cm³",
                                    line=dict(color="#00b4d8", width=2.5),
                                    marker=dict(size=7, color="#00b4d8"),
                                ))
                            _fig_flux.update_layout(
                                **DARK, height=400,
                                title=dict(
                                    text=f"Integrated flux vs elevation — az = {_az_slice_deg:.1f}°",
                                    font=dict(size=12)
                                ),
                                xaxis=dict(title="Elevation above horizon [°]", gridcolor="#2a2a3a"),
                                yaxis=dict(
                                    title="Integrated flux [m⁻² s⁻¹ sr⁻¹]",
                                    gridcolor="#2a2a3a",
                                    type="log" if _log_flux else "linear",
                                ),
                                legend=dict(bgcolor="rgba(0,0,0,0.5)", font=dict(size=10)),
                                margin=dict(l=75, r=20, t=50, b=55),
                            )
                            st.plotly_chart(_fig_flux, config={"displayModeBar": True}, key="comp_flux_plot")
                            st.caption(
                                "UCMuon backward CSDA + stochastic MC: open-sky (gold dashed) vs "
                                "terrain-blocked (blue) integrated flux. A published MURAVES flux "
                                "reference can be overlaid here once digitised from a peer-reviewed source."
                            )
                    else:
                        st.info("Install ucmuon_backward_mc.py alongside gui_terrain_engine.py to enable flux plots.")

                    # ── 5 — Export thickness map ──────────────────────────────────────────
                    st.markdown("#### 5 — Export thickness map for external comparison")
                    st.caption(
                        "Download the thickness map as text (cols: az[°], el[°], thickness[m], "
                        "overburden[g/cm²], open_sky). The format matches TURTLE/Mulder-style "
                        "overburden input for external cross-checks."
                    )
                    _lines_export = ["# az[deg]  el[deg]  thickness[m]  overburden[g/cm2]  open_sky\n"]
                    for _ia_e, _az_e in enumerate(az_c):
                        for _iz_e, _ze_e in enumerate(ze_c):
                            _el_e = 90.0 - _ze_e
                            _th_e = float(ob_map[_ia_e, _iz_e]) / (_rho_comp * 100.0)
                            _ob_e = float(ob_map[_ia_e, _iz_e])
                            _sk_e = int(sky_m[_ia_e, _iz_e])
                            _lines_export.append(f"{_az_e:8.2f}  {_el_e:7.2f}  {_th_e:12.2f}  {_ob_e:14.2f}  {_sk_e}\n")
                    st.download_button(
                        "⬇️  Download thickness map (az, el, thickness [m], overburden [g/cm²])",
                        data="".join(_lines_export),
                        file_name="ucmuon_thickness_map.dat",
                        mime="text/plain",
                        width='stretch',
                        key="comp_export_thick",
                    )

                # Download
                out_path = st.session_state.get("terrain_result_outfile", t_outfile)
                if Path(out_path).exists():
                    st.markdown("**Output file** (18-col underground format, same as MUSIC / Transport-tab output)")
                    st.caption("This file can be loaded in the Results tab for full analysis, 3D plots, and export.")
                    with open(out_path, "rb") as fh:
                        st.download_button(
                            f"⬇️  Download `{Path(out_path).name}`  ({Path(out_path).stat().st_size/1024:.0f} kB)",
                            data=fh, file_name=Path(out_path).name, mime="text/plain",
                            width='stretch', key="terrain_dl_ug",
                        )
                    st.info(
                        f"💡 Load `{out_path}` in the **Results** tab to see energy spectra, "
                        "3D trajectories, survival rate curves, and PHITS/Geant4 export.",
                        icon="💡",
                    )

