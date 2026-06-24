# gui/gui_density_analysis.py  —  UCLouvain Muography Group
#
# Streamlit GUI tab for density inversion from muon transmission maps.
# Public API:
#   render_density_analysis_tab()  — renders the full Density Analysis panel
#
# Author: Hamid Basiri <hamid.basiri@uclouvain.be>
# MIT License 2026

import sys
import importlib.util
import io
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from pathlib import Path

_DA_VERSION = "1.0.0"

# ─────────────────────────────────────────────────────────────────────────────
# Dark theme — matches the rest of the GUI
# ─────────────────────────────────────────────────────────────────────────────

DARK = dict(
    paper_bgcolor="rgb(15,17,23)",
    plot_bgcolor="rgb(20,22,30)",
    font=dict(color="white", size=11),
)

# ─────────────────────────────────────────────────────────────────────────────
# Lazy-load the pure-Python analysis engine
# ─────────────────────────────────────────────────────────────────────────────

def _load_da_engine():
    gui_dir  = Path(__file__).resolve().parent
    da_path  = gui_dir / "ucmuon_density_analysis.py"
    if not da_path.exists():
        raise FileNotFoundError(
            f"ucmuon_density_analysis.py not found in {gui_dir}"
        )
    spec = importlib.util.spec_from_file_location("ucmuon_density_analysis",
                                                   str(da_path))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ─────────────────────────────────────────────────────────────────────────────
# Plotly helpers
# ─────────────────────────────────────────────────────────────────────────────

def _az_el_heatmap(az_c, el_c, data_2d, title, unit,
                   colorscale="plasma", zmin=None, zmax=None,
                   height=340, colorbar_x=1.01):
    z_plot = data_2d.T.copy()
    finite = z_plot[~np.isnan(z_plot)]
    if zmin is None:
        zmin = float(finite.min()) if finite.size > 0 else 0.0
    if zmax is None:
        zmax = float(finite.max()) if finite.size > 0 else 1.0
    if zmax <= zmin:
        zmax = zmin + 1.0

    fig = go.Figure(go.Heatmap(
        x=az_c, y=el_c, z=z_plot,
        colorscale=colorscale,
        zmin=zmin, zmax=zmax,
        colorbar=dict(
            title=dict(text=unit, font=dict(color="white", size=11)),
            tickfont=dict(color="white"),
            x=colorbar_x, thickness=14,
        ),
        hovertemplate="Az: %{x:.1f}°  El: %{y:.1f}°<br>" + unit + ": %{z:.3g}<extra></extra>",
        hoverongaps=False,
        xgap=1, ygap=1,
    ))
    fig.update_layout(
        **DARK,
        height=height,
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(title="Azimuth [°]  (N=0°, E=90°)", gridcolor="#2a2a3a"),
        yaxis=dict(title="Elevation [°]", gridcolor="#2a2a3a"),
        margin=dict(l=60, r=70, t=50, b=55),
    )
    return fig


def _status_heatmap(az_c, el_c, status_map, height=340):
    _STATUS_COLORS = {
        0: "rgb(0,200,80)",
        1: "rgb(80,80,100)",
        2: "rgb(255,200,0)",
        3: "rgb(220,50,50)",
        4: "rgb(50,120,255)",
    }
    _STATUS_LABELS = {0: "OK", 1: "Open sky", 2: "ρ < ρ_min", 3: "ρ > ρ_max", 4: "Low sensitivity"}
    z_plot = status_map.T.astype(float)
    colorscale = [
        [0.00, _STATUS_COLORS[0]], [0.20, _STATUS_COLORS[0]],
        [0.20, _STATUS_COLORS[1]], [0.40, _STATUS_COLORS[1]],
        [0.40, _STATUS_COLORS[2]], [0.60, _STATUS_COLORS[2]],
        [0.60, _STATUS_COLORS[3]], [0.80, _STATUS_COLORS[3]],
        [0.80, _STATUS_COLORS[4]], [1.00, _STATUS_COLORS[4]],
    ]
    fig = go.Figure(go.Heatmap(
        x=az_c, y=el_c, z=z_plot,
        colorscale=colorscale,
        zmin=0, zmax=4,
        colorbar=dict(
            title=dict(text="Status", font=dict(color="white", size=11)),
            tickfont=dict(color="white"),
            tickvals=[0.4, 1.2, 2.0, 2.8, 3.6],
            ticktext=[_STATUS_LABELS[k] for k in range(5)],
            x=1.01, thickness=14,
        ),
        hovertemplate="Az: %{x:.1f}°  El: %{y:.1f}°<br>Status: %{z:.0f}<extra></extra>",
        hoverongaps=False,
        xgap=1, ygap=1,
    ))
    fig.update_layout(
        **DARK,
        height=height,
        title=dict(text="Inversion status map", font=dict(size=12)),
        xaxis=dict(title="Azimuth [°]", gridcolor="#2a2a3a"),
        yaxis=dict(title="Elevation [°]", gridcolor="#2a2a3a"),
        margin=dict(l=60, r=100, t=50, b=55),
    )
    return fig


def _chi2_line_plot(rho_fine, chi2_curve, rho_hat, title="χ²(ρ) landscape"):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=rho_fine, y=chi2_curve,
        mode="lines",
        line=dict(color="#4fc3f7", width=2),
        name="χ²",
        hovertemplate="ρ = %{x:.3f} g/cm³<br>χ² = %{y:.2f}<extra></extra>",
    ))
    if not np.isnan(rho_hat):
        fig.add_vline(
            x=float(rho_hat),
            line=dict(color="#ffd700", width=1.5, dash="dash"),
            annotation_text=f"ρ̂ = {rho_hat:.3f}",
            annotation_font=dict(color="#ffd700", size=11),
        )
    fig.update_layout(
        **DARK,
        height=280,
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(title="Density [g/cm³]", gridcolor="#2a2a3a"),
        yaxis=dict(title="χ²", gridcolor="#2a2a3a"),
        margin=dict(l=65, r=20, t=50, b=55),
        showlegend=False,
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Tab 1 — T_sim library
# ─────────────────────────────────────────────────────────────────────────────

def _tab_library(da):
    st.caption(
        "Run the Terrain engine at **3–5 different densities** (same DEM and az/el grid each time).  "
        "After each run, open the **📈 Transmission** sub-tab and click **💾 Save T_sim** — "
        "the file is auto-queued here.  Then click **📥 Import** or paste paths manually."
    )

    # ── Auto-import from Terrain tab ─────────────────────────────────────────
    _pending = list(st.session_state.get("da_lib_paths_pending", []))
    if _pending:
        _pend_names = ", ".join(Path(p).name for p in _pending)
        if st.button(
            f"📥 Import {len(_pending)} file(s) queued from Terrain tab  ({_pend_names})",
            key="da_import_terrain_btn",
            type="primary",
        ):
            _cur_text = st.session_state.get("da_lib_paths_text", "")
            _cur_set  = {l.strip() for l in _cur_text.splitlines() if l.strip()}
            _new      = [p for p in _pending if p not in _cur_set]
            if _new:
                _merged = (_cur_text.rstrip("\n") + "\n" + "\n".join(_new)).lstrip("\n")
                st.session_state["da_lib_paths_text"] = _merged
            st.session_state["da_lib_paths_pending"] = []
            st.rerun()

    # ── File path text area ───────────────────────────────────────────────────
    default_paths = "\n".join(
        str(p) for p in sorted(Path(".").glob("terrain_transmission*.dat"))
    )
    file_text = st.text_area(
        "T_sim files — one path per line",
        value=st.session_state.get("da_lib_paths_text", default_paths),
        height=110,
        key="da_lib_paths_text",
        help="Full paths to terrain_transmission.dat files, one density per file.",
    )

    if st.button("Load library", key="da_load_lib_btn", type="primary"):
        paths = [p.strip() for p in file_text.splitlines() if p.strip()]
        if not paths:
            st.warning("Enter at least one file path above.")
        else:
            with st.spinner("Loading T_sim library…"):
                try:
                    lib    = da.build_tsim_library(paths)
                    az_c, el_c, _, _ = da.load_transmission_map(paths[0])
                    st.session_state["da_tsim_lib"]   = lib
                    st.session_state["da_lib_az_c"]   = az_c
                    st.session_state["da_lib_el_c"]   = el_c
                    st.session_state["da_lib_paths"]  = paths
                    st.success(
                        f"✅  {len(lib)} density maps loaded  |  "
                        f"ρ ∈ [{min(lib):.2f}, {max(lib):.2f}] g/cm³  |  "
                        f"grid {len(az_c)} az × {len(el_c)} el"
                    )
                except Exception as _e:
                    import traceback
                    st.error(f"Failed to load library: {_e}")
                    st.code(traceback.format_exc())

    # ── Library summary + preview ─────────────────────────────────────────────
    lib  = st.session_state.get("da_tsim_lib")
    az_c = st.session_state.get("da_lib_az_c")
    el_c = st.session_state.get("da_lib_el_c")
    if lib:
        import pandas as pd
        rows = [
            {
                "ρ [g/cm³]": f"{rho:.3f}",
                "n_az": T2d.shape[0],
                "n_el": T2d.shape[1],
                "T_min": f"{float(np.nanmin(T2d)):.4f}",
                "T_max": f"{float(np.nanmax(T2d)):.4f}",
                "T_mean": f"{float(np.nanmean(T2d)):.4f}",
            }
            for rho, T2d in lib.items()
        ]
        st.dataframe(pd.DataFrame(rows), width='stretch', hide_index=True)

        if az_c is not None and el_c is not None:
            rho_keys = sorted(lib.keys())
            sel_rho  = st.select_slider(
                "Preview density",
                options=rho_keys,
                format_func=lambda r: f"{r:.3f} g/cm³",
                key="da_preview_rho",
            )
            st.plotly_chart(
                _az_el_heatmap(
                    az_c, el_c, lib[sel_rho],
                    title=f"T_sim(az, el)  |  ρ = {sel_rho:.3f} g/cm³",
                    unit="T_sim", colorscale="viridis",
                    zmin=0.0, zmax=1.0,
                ),
            )

    # ── Detector geometry check (secondary — tucked away) ────────────────────
    with st.expander("🔭 Detector geometry — angular resolution check", expanded=False):
        _detector_geometry_check(az_c, el_c)


def _detector_geometry_check(az_c, el_c):
    st.caption(
        "Optional consistency check: verify that the simulation grid bin size "
        "matches your detector's native angular resolution."
    )
    import math

    _dc1, _dc2 = st.columns(2)
    _baseline = _dc1.number_input("Baseline L [m]",   0.01, 50.0, 1.0, 0.1, key="da_det_baseline",
                                   help="Distance between the two outermost tracking planes.")
    _pitch_x  = _dc1.number_input("Strip pitch X [cm]", 0.1, 50.0, 5.0, 0.5, key="da_det_pitch_x")
    _pitch_y  = _dc1.number_input("Strip pitch Y [cm]", 0.1, 50.0, 5.0, 0.5, key="da_det_pitch_y")
    _n_ch_x   = _dc2.number_input("Channels X", 1, 10000, 16, 1, key="da_det_nch_x")
    _n_ch_y   = _dc2.number_input("Channels Y", 1, 10000, 16, 1, key="da_det_nch_y")

    _daz = math.degrees(math.atan(_pitch_x * 1e-2 / _baseline))
    _del = math.degrees(math.atan(_pitch_y * 1e-2 / _baseline))
    _faz = math.degrees(math.atan(_n_ch_x * _pitch_x * 1e-2 / 2.0 / _baseline)) * 2
    _fel = math.degrees(math.atan(_n_ch_y * _pitch_y * 1e-2 / 2.0 / _baseline)) * 2

    _m1, _m2, _m3, _m4 = st.columns(4)
    _m1.metric("δaz / pixel", f"{_daz:.2f}°")
    _m2.metric("δel / pixel", f"{_del:.2f}°")
    _m3.metric("FoV az",      f"±{_faz/2:.1f}°")
    _m4.metric("FoV el",      f"±{_fel/2:.1f}°")

    if az_c is not None and el_c is not None and len(az_c) > 1 and len(el_c) > 1:
        _sdaz = float(az_c[1] - az_c[0])
        _sdel = float(el_c[1] - el_c[0])
        _ok   = (abs(_sdaz - _daz) / _daz < 0.5) and (abs(_sdel - _del) / _del < 0.5)
        if _ok:
            st.success(f"✅ Grid (Δaz={_sdaz:.2f}°, Δel={_sdel:.2f}°) matches detector resolution.")
        else:
            st.warning(
                f"Simulation grid Δaz={_sdaz:.2f}°, Δel={_sdel:.2f}° differs from "
                f"detector resolution δaz={_daz:.2f}°, δel={_del:.2f}°.  "
                "Consider re-running the terrain simulation with matching grid spacing."
            )
    else:
        st.info("Load the T_sim library to compare the grid spacing.")


# ─────────────────────────────────────────────────────────────────────────────
# Tab 2 — Measured / synthetic T_data
# ─────────────────────────────────────────────────────────────────────────────

def _tab_tdata(da):
    lib  = st.session_state.get("da_tsim_lib")
    az_c = st.session_state.get("da_lib_az_c")
    el_c = st.session_state.get("da_lib_el_c")

    if lib is None:
        st.info("⬅️  Load the T_sim library in **① T_sim Library** before loading T_data.", icon="ℹ️")

    mode = st.radio(
        "Data source",
        ["📂 Upload / file path", "🔬 Generate synthetic (test)"],
        key="da_tdata_mode",
        horizontal=True,
    )

    if mode == "📂 Upload / file path":
        st.caption(
            "Load a `terrain_transmission.dat` from a real measurement campaign.  "
            "The grid must match the T_sim library (same n_az, n_el, same bin edges)."
        )
        _r1, _r2 = st.columns([3, 1])
        tdata_path = _r1.text_input(
            "File path",
            value=st.session_state.get("da_tdata_path", ""),
            key="da_tdata_path",
        )
        uploaded = _r2.file_uploader("Or upload", type=["dat", "txt"],
                                     key="da_tdata_upload")

        if st.button("Load T_data", key="da_load_tdata_btn", type="primary"):
            try:
                if uploaded is not None:
                    import tempfile, os
                    suffix = Path(uploaded.name).suffix or ".dat"
                    with tempfile.NamedTemporaryFile(mode="wb", suffix=suffix, delete=False) as tmp:
                        tmp.write(uploaded.read())
                        tmp_path = tmp.name
                    az_c_d, el_c_d, T_2d, meta = da.load_transmission_map(tmp_path)
                    os.unlink(tmp_path)
                elif tdata_path:
                    az_c_d, el_c_d, T_2d, meta = da.load_transmission_map(tdata_path)
                else:
                    st.warning("Provide a file path or upload a file.")
                    return

                if az_c is not None:
                    if not (np.allclose(az_c, az_c_d, atol=1e-3) and
                            np.allclose(el_c, el_c_d, atol=1e-3)):
                        st.error(
                            "Grid mismatch: T_data grid does not match the T_sim library.  "
                            "Ensure all files were produced with the same az/el grid."
                        )
                        return

                st.session_state["da_T_data"]   = T_2d
                st.session_state["da_sigma_T"]  = None
                st.session_state["da_lib_az_c"] = az_c_d
                st.session_state["da_lib_el_c"] = el_c_d
                st.success(
                    f"✅  T_data loaded — "
                    f"grid {len(az_c_d)} az × {len(el_c_d)} el  |  "
                    f"T ∈ [{float(np.nanmin(T_2d)):.4f}, {float(np.nanmax(T_2d)):.4f}]  |  "
                    f"ρ_file = {meta['density']:.3f} g/cm³"
                )

            except Exception as _e:
                import traceback
                st.error(f"Failed to load T_data: {_e}")
                st.code(traceback.format_exc())

    else:
        st.caption(
            "Creates a mock dataset by evaluating T_sim at ρ_true and adding Poisson noise.  "
            "Use this to verify the inversion recovers the known density before applying to real data."
        )
        if lib is None:
            st.info("Load the T_sim library first (**① T_sim Library**).")
            return

        _s1, _s2 = st.columns(2)
        rho_true = _s1.slider(
            "True density ρ_true [g/cm³]",
            min_value=float(min(lib.keys())), max_value=float(max(lib.keys())),
            value=min(2.65, float(max(lib.keys()))),
            step=0.05, key="da_synth_rho",
        )
        n_events = _s2.number_input(
            "Muon statistics N",
            min_value=100, max_value=1_000_000, value=10_000, step=1_000,
            key="da_synth_nevents",
            help="Higher N → smaller Poisson noise → tighter σ_ρ.",
        )

        if st.button("Generate synthetic T_data", key="da_gen_synth_btn", type="primary"):
            with st.spinner("Generating with Poisson noise…"):
                try:
                    T_data, sigma_T = da.generate_synthetic_tdata(
                        lib, true_rho=rho_true, n_events=int(n_events), seed=42
                    )
                    st.session_state["da_T_data"]         = T_data
                    st.session_state["da_sigma_T"]        = sigma_T
                    st.session_state["da_synth_rho_used"] = rho_true
                    st.success(
                        f"✅  Synthetic T_data generated — ρ_true = {rho_true:.2f} g/cm³  |  "
                        f"N = {int(n_events):,}  |  σ_T mean = {float(np.nanmean(sigma_T)):.4f}"
                    )
                except Exception as _e:
                    import traceback
                    st.error(f"Synthetic generation failed: {_e}")
                    st.code(traceback.format_exc())

    # ── T_data preview ────────────────────────────────────────────────────────
    T_data = st.session_state.get("da_T_data")
    if T_data is not None and az_c is not None and el_c is not None:
        st.plotly_chart(
            _az_el_heatmap(
                az_c, el_c, T_data,
                title="T_data(az, el)",
                unit="T", colorscale="cividis",
                zmin=0.0, zmax=1.0,
            ),
        )


# ─────────────────────────────────────────────────────────────────────────────
# Tab 3 — Run inversion
# ─────────────────────────────────────────────────────────────────────────────

def _tab_inversion(da):
    lib     = st.session_state.get("da_tsim_lib")
    T_data  = st.session_state.get("da_T_data")
    sigma_T = st.session_state.get("da_sigma_T")
    az_c    = st.session_state.get("da_lib_az_c")
    el_c    = st.session_state.get("da_lib_el_c")

    # Prerequisite checklist
    _ready = True
    if lib is None:
        st.info("T_sim library not loaded — load it in **① T_sim Library** first.", icon="ℹ️")
        _ready = False
    else:
        st.success(
            f"T_sim library: {len(lib)} maps, ρ ∈ [{min(lib):.2f}, {max(lib):.2f}] g/cm³",
            icon="✅",
        )

    if T_data is None:
        st.info("T_data not loaded — load it in **② Measured Data** first.", icon="ℹ️")
        _ready = False
    else:
        _td_ok = int(np.sum(~np.isnan(T_data)))
        st.success(f"T_data: {_td_ok} valid pixels", icon="✅")

    if not _ready:
        return

    st.divider()

    rho_keys = sorted(lib.keys())
    _c1, _c2 = st.columns(2)
    ref_rho = _c1.selectbox(
        "Reference density for double-ratio map",
        options=rho_keys,
        index=len(rho_keys) // 2,
        format_func=lambda r: f"{r:.3f} g/cm³",
        key="da_ref_rho",
        help="D = T_data / T_sim(ρ_ref).  Choose a plausible background density.",
    )
    min_sens = _c2.number_input(
        "Min sensitivity |dT/dρ|",
        min_value=0.0001, max_value=0.1,
        value=0.005, step=0.001, format="%.4f",
        key="da_min_sens",
        help="Pixels with sensitivity below this are flagged as unreliable (status 4).",
    )

    if st.button("▶  Run Inversion", key="da_run_inversion_btn", type="primary",
                 width='stretch'):
        if T_data.shape != lib[rho_keys[0]].shape:
            st.error(
                f"Shape mismatch: T_data {T_data.shape} vs library {lib[rho_keys[0]].shape}.  "
                "Reload T_data after loading the library."
            )
            return

        with st.status("Running density inversion…", expanded=True) as status_box:
            st.write("Inverting density map pixel-by-pixel…")
            try:
                rho_map, sigma_rho, status_map = da.invert_density_map(
                    T_data, lib, sigma_T=sigma_T,
                    min_sensitivity=float(min_sens),
                )
                st.write("Computing double ratio…")
                D_map = da.compute_double_ratio(T_data, lib[ref_rho])

                st.session_state["da_rho_map"]      = rho_map
                st.session_state["da_sigma_rho"]    = sigma_rho
                st.session_state["da_status_map"]   = status_map
                st.session_state["da_D_map"]        = D_map
                st.session_state["da_ref_rho_used"] = ref_rho
                st.session_state["da_az_c_result"]  = az_c
                st.session_state["da_el_c_result"]  = el_c

                n_ok    = int(np.sum(status_map == 0))
                n_sky   = int(np.sum(status_map == 1))
                n_flag  = int(np.sum(status_map >= 2))
                rho_ok  = rho_map[status_map == 0]
                rho_med = float(np.nanmedian(rho_ok)) if rho_ok.size else float("nan")

                status_box.update(label="Inversion complete!", state="complete")
                st.write(
                    f"OK: {n_ok}  |  Open sky: {n_sky}  |  "
                    f"Flagged: {n_flag}  |  Median ρ̂ = {rho_med:.3f} g/cm³"
                )

            except Exception as _e:
                import traceback
                status_box.update(label="Inversion failed!", state="error")
                st.error(f"Inversion error: {_e}")
                st.code(traceback.format_exc())


# ─────────────────────────────────────────────────────────────────────────────
# Tab 4 — Results (sub-tabbed)
# ─────────────────────────────────────────────────────────────────────────────

def _tab_results(da):
    rho_map    = st.session_state.get("da_rho_map")
    sigma_rho  = st.session_state.get("da_sigma_rho")
    status_map = st.session_state.get("da_status_map")
    D_map      = st.session_state.get("da_D_map")
    az_c       = st.session_state.get("da_az_c_result")
    el_c       = st.session_state.get("da_el_c_result")
    T_data     = st.session_state.get("da_T_data")
    sigma_T    = st.session_state.get("da_sigma_T")
    lib        = st.session_state.get("da_tsim_lib")
    ref_rho    = st.session_state.get("da_ref_rho_used", 2.65)
    synth_rho  = st.session_state.get("da_synth_rho_used")

    if rho_map is None:
        st.info("⬅️  Run the inversion in **③ Run Inversion** to see results here.", icon="ℹ️")
        return

    lib_rho_keys = sorted(lib.keys()) if lib else []
    rho_lo = lib_rho_keys[0]  if lib_rho_keys else 0.5
    rho_hi = lib_rho_keys[-1] if lib_rho_keys else 3.5

    # Summary metrics row — always visible
    n_ok   = int(np.sum(status_map == 0))
    n_sky  = int(np.sum(status_map == 1))
    n_flag = int(np.sum(status_map >= 2))
    rho_ok_vals = rho_map[status_map == 0]
    rho_med = float(np.nanmedian(rho_ok_vals)) if rho_ok_vals.size > 0 else float("nan")
    rho_std = float(np.nanstd(rho_ok_vals))    if rho_ok_vals.size > 0 else float("nan")

    _m1, _m2, _m3, _m4 = st.columns(4)
    _m1.metric("Median ρ̂  (OK bins)",   f"{rho_med:.3f} g/cm³",
               delta=(f"{rho_med - ref_rho:+.3f}" if not np.isnan(rho_med) else None))
    _m2.metric("σ(ρ̂) across OK bins",  f"{rho_std:.3f} g/cm³")
    _m3.metric("OK bins",                f"{n_ok} / {rho_map.size}")
    _m4.metric("Open sky / flagged",     f"{n_sky} / {n_flag}")

    if synth_rho is not None and not np.isnan(rho_med):
        bias = rho_med - synth_rho
        st.info(
            f"Synthetic validation — ρ_true = {synth_rho:.3f} g/cm³  |  "
            f"ρ̂_median = {rho_med:.3f} g/cm³  |  "
            f"bias = {bias:+.4f} g/cm³  ({100 * bias / synth_rho:+.2f}%)"
        )

    # Sub-tabs for the four result views
    _rt1, _rt2, _rt3, _rt4 = st.tabs([
        "🗺️ Density map",
        "📊 Diagnostics",
        "📉 Chi² landscape",
        "💾 Download",
    ])

    # ── Sub-tab 1: density + uncertainty ─────────────────────────────────────
    with _rt1:
        _col_a, _col_b = st.columns(2)
        with _col_a:
            st.plotly_chart(
                _az_el_heatmap(
                    az_c, el_c, rho_map,
                    title="Inverted density  ρ̂(az, el)",
                    unit="ρ̂ [g/cm³]", colorscale="plasma",
                    zmin=rho_lo, zmax=rho_hi,
                ),
            )
        with _col_b:
            has_sigma = (sigma_rho is not None and not np.all(np.isnan(sigma_rho)))
            if has_sigma:
                st.plotly_chart(
                    _az_el_heatmap(
                        az_c, el_c, sigma_rho,
                        title="Density uncertainty  σ_ρ(az, el)",
                        unit="σ_ρ [g/cm³]", colorscale="viridis",
                    ),
                )
            else:
                st.info(
                    "Uncertainty map unavailable.  "
                    "Use the synthetic data mode (**② Measured Data**) or supply a σ_T file to enable error propagation."
                )

    # ── Sub-tab 2: double ratio + status map ─────────────────────────────────
    with _rt2:
        _col_c, _col_d = st.columns(2)
        with _col_c:
            D_plot = D_map.copy()
            D_plot[status_map == 1] = np.nan
            st.plotly_chart(
                _az_el_heatmap(
                    az_c, el_c, D_plot,
                    title=f"Double ratio  D = T_data / T_sim(ρ={ref_rho:.2f})",
                    unit="D", colorscale="RdBu",
                    zmin=0.5, zmax=2.0,
                ),
            )
        with _col_d:
            st.plotly_chart(
                _status_heatmap(az_c, el_c, status_map),
            )
        st.caption(
            "**Double ratio D > 1** → more muons observed than simulated at ρ_ref → "
            "actual density lower than ρ_ref (void / low-density zone).  "
            "**D < 1** → denser than ρ_ref.  "
            "**Status colours**: green=OK, grey=open sky, yellow=below ρ_min, "
            "red=above ρ_max, blue=low sensitivity."
        )

    # ── Sub-tab 3: chi² landscape ─────────────────────────────────────────────
    with _rt3:
        if sigma_T is None or lib is None:
            st.info(
                "Chi-squared landscape requires σ_T.  "
                "Use the synthetic mode (**② Measured Data**) to produce σ_T automatically, "
                "or supply a sigma file.",
                icon="ℹ️",
            )
        else:
            n_az, n_el = rho_map.shape
            _cp1, _cp2 = st.columns(2)
            az_pick_idx = _cp1.slider("Azimuth bin", 0, n_az - 1, n_az // 2, key="da_chi2_az_idx")
            el_pick_idx = _cp2.slider("Elevation bin", 0, n_el - 1, n_el // 2, key="da_chi2_el_idx")

            az_val  = float(az_c[az_pick_idx])
            el_val  = float(el_c[el_pick_idx])
            rho_hat = float(rho_map[az_pick_idx, el_pick_idx])
            sT_val  = float(sigma_T[az_pick_idx, el_pick_idx])
            td_val  = float(T_data[az_pick_idx, el_pick_idx])

            st.metric(
                f"ρ̂  at  az = {az_val:.1f}°,  el = {el_val:.1f}°",
                f"{rho_hat:.4f} g/cm³" if not np.isnan(rho_hat) else "NaN",
                delta=(f"{rho_hat - ref_rho:+.4f} vs ref" if not np.isnan(rho_hat) else None),
            )

            if np.isnan(td_val) or sT_val <= 0.0 or np.isnan(sT_val):
                st.caption(
                    f"Pixel ({az_pick_idx}, {el_pick_idx}): "
                    f"T_data={td_val:.4f}, σ_T={sT_val:.4g} — χ² not available."
                )
            else:
                try:
                    rho_fine, chi2_curve = da.compute_chi2_landscape(
                        T_data, lib, sigma_T,
                        az_idx=az_pick_idx, el_idx=el_pick_idx, n_rho_fine=200,
                    )
                    st.plotly_chart(
                        _chi2_line_plot(
                            rho_fine, chi2_curve, rho_hat,
                            title=f"χ²(ρ)  —  az={az_val:.1f}°  el={el_val:.1f}°",
                        ),
                    )
                except Exception as _e:
                    st.caption(f"Chi-squared plot unavailable: {_e}")

    # ── Sub-tab 4: download ───────────────────────────────────────────────────
    with _rt4:
        st.caption("Write the full density map (ρ̂, σ_ρ, status) to a .dat file for further analysis.")
        if st.button("Prepare density_map.dat", key="da_prep_download", type="primary"):
            with st.spinner("Writing…"):
                try:
                    import tempfile, os
                    with tempfile.NamedTemporaryFile(mode="w", suffix=".dat", delete=False) as tmp:
                        tmp_path = tmp.name
                    da.write_density_map(
                        az_c, el_c, rho_map, sigma_rho, status_map,
                        tmp_path,
                        metadata={"rho_ref": ref_rho,
                                  "n_events": st.session_state.get("da_synth_nevents")},
                    )
                    with open(tmp_path, "rb") as fh:
                        file_bytes = fh.read()
                    os.unlink(tmp_path)
                    st.session_state["da_download_bytes"] = file_bytes
                    st.success(f"Ready — {len(file_bytes)//1024 + 1} kB")
                except Exception as _e:
                    st.error(f"Failed to write density map: {_e}")

        dl_bytes = st.session_state.get("da_download_bytes")
        if dl_bytes is not None:
            st.download_button(
                label="⬇ Download density_map.dat",
                data=dl_bytes,
                file_name="density_map.dat",
                mime="text/plain",
                key="da_download_btn",
            )


# ─────────────────────────────────────────────────────────────────────────────
# Direct opacity inversion workflow  (methods 2 & 3 — single open-sky + target)
# ─────────────────────────────────────────────────────────────────────────────

def _grid_matches(az_a, el_a, az_b, el_b):
    return (az_a is not None and az_b is not None and
            len(az_a) == len(az_b) and len(el_a) == len(el_b) and
            np.allclose(az_a, az_b, atol=1e-3) and
            np.allclose(el_a, el_b, atol=1e-3))


def _di_load_input(da, method_key):
    """Step ① — load the measured transmission map (or two flux maps)."""
    st.markdown("**① Measured data**")

    conv = st.radio(
        "Angle column convention",
        ["Elevation (0°=horizon)", "Zenith (0°=vertical)"],
        key=f"di_conv_{method_key}", horizontal=True,
        help="Terrain transmission files use elevation; raw flux maps use zenith.",
    )
    angle_is_zenith = conv.startswith("Zenith")

    if method_key == "direct":
        st.caption(
            "Load one measured transmission map  T = Φ_target / Φ_open  "
            "(columns: az, angle, T).  No density header required."
        )
        _c1, _c2 = st.columns([3, 1])
        path = _c1.text_input("Transmission file path", key="di_tpath",
                              value=st.session_state.get("di_tpath", ""))
        up = _c2.file_uploader("Or upload", type=["dat", "txt"], key="di_tupload")
        if st.button("Load transmission", key="di_load_t", type="primary"):
            try:
                src = _stash_upload(up) if up is not None else path
                if not src:
                    st.warning("Provide a path or upload a file.")
                    return
                az_c, el_c, T_2d, _ = da.load_generic_map(
                    src, angle_is_zenith=angle_is_zenith, value_col=2)
                st.session_state["di_T"]    = T_2d
                st.session_state["di_az"]   = az_c
                st.session_state["di_el"]   = el_c
                st.session_state["di_sigma"] = None
                st.success(f"✅ T loaded — {len(az_c)} az × {len(el_c)} el  |  "
                           f"T ∈ [{np.nanmin(T_2d):.3f}, {np.nanmax(T_2d):.3f}]")
            except Exception as _e:
                import traceback
                st.error(f"Load failed: {_e}")
                st.code(traceback.format_exc())

    else:  # fluxratio
        st.caption(
            "Load two raw flux maps — open-sky reference and target — "
            "(columns: az, angle, flux).  T = Φ_target / Φ_open is formed here."
        )
        _c1, _c2 = st.columns(2)
        with _c1:
            st.markdown("Target (through rock)")
            tgt_path = st.text_input("Target flux path", key="di_tgt_path",
                                     value=st.session_state.get("di_tgt_path", ""))
            tgt_up = st.file_uploader("Or upload target", type=["dat", "txt"],
                                      key="di_tgt_up")
        with _c2:
            st.markdown("Open-sky (reference)")
            opn_path = st.text_input("Open-sky flux path", key="di_opn_path",
                                     value=st.session_state.get("di_opn_path", ""))
            opn_up = st.file_uploader("Or upload open-sky", type=["dat", "txt"],
                                      key="di_opn_up")
        if st.button("Form transmission ratio", key="di_load_ratio", type="primary"):
            try:
                tsrc = _stash_upload(tgt_up) if tgt_up is not None else tgt_path
                osrc = _stash_upload(opn_up) if opn_up is not None else opn_path
                if not tsrc or not osrc:
                    st.warning("Provide both target and open-sky maps.")
                    return
                az_t, el_t, F_t, _ = da.load_generic_map(
                    tsrc, angle_is_zenith=angle_is_zenith, value_col=2)
                az_o, el_o, F_o, _ = da.load_generic_map(
                    osrc, angle_is_zenith=angle_is_zenith, value_col=2)
                if not _grid_matches(az_t, el_t, az_o, el_o):
                    st.error("Target and open-sky maps are on different grids.")
                    return
                T_2d = da.transmission_from_flux_maps(F_t, F_o)
                st.session_state["di_T"]     = T_2d
                st.session_state["di_az"]    = az_t
                st.session_state["di_el"]    = el_t
                st.session_state["di_sigma"] = None
                st.success(f"✅ T = Φ_target/Φ_open formed — "
                           f"{len(az_t)} az × {len(el_t)} el  |  "
                           f"T ∈ [{np.nanmin(T_2d):.3f}, {np.nanmax(T_2d):.3f}]")
            except Exception as _e:
                import traceback
                st.error(f"Ratio failed: {_e}")
                st.code(traceback.format_exc())

    # Optional σ_T from counting statistics
    T_2d = st.session_state.get("di_T")
    if T_2d is not None:
        with st.expander("σ_T — uncertainty from counting statistics (optional)"):
            use_sig = st.checkbox("Estimate σ_T ≈ √[T(1−T)/N]", key="di_use_sig",
                                  help="Enables σ_ρ error maps. N = open-sky counts/bin.")
            if use_sig:
                N = st.number_input("Open-sky counts N per bin", 10, 10_000_000,
                                    10_000, 100, key="di_N")
                with np.errstate(invalid="ignore"):
                    st.session_state["di_sigma"] = np.sqrt(
                        np.clip(T_2d * (1.0 - T_2d), 0.0, None) / float(N))
            else:
                st.session_state["di_sigma"] = None

        az_c = st.session_state.get("di_az")
        el_c = st.session_state.get("di_el")
        if az_c is not None:
            st.plotly_chart(
                _az_el_heatmap(az_c, el_c, T_2d, "T_data(az, el)",
                               "T", colorscale="cividis", zmin=0.0, zmax=1.0,
                               height=300))


def _stash_upload(uploaded):
    """Write a Streamlit upload to a temp file and return its path."""
    import tempfile
    suffix = Path(uploaded.name).suffix or ".dat"
    with tempfile.NamedTemporaryFile(mode="wb", suffix=suffix, delete=False) as tmp:
        tmp.write(uploaded.read())
        return tmp.name


def _di_path_length(da):
    """Step ② — configure the line-of-sight path length L(az, el)."""
    st.markdown("**② Path length  L(az, el)**  —  needed to convert opacity → density")

    az_c = st.session_state.get("di_az")
    el_c = st.session_state.get("di_el")

    src = st.radio(
        "Path-length source",
        ["Derive from a T_sim map (known ρ_sim)",
         "From a terrain overburden map",
         "Opacity only — skip L (no density)"],
        key="di_L_source",
    )
    model = st.session_state.get("di_model", "reyna_bugaev")
    alt   = st.session_state.get("di_alt", 0.0)

    if src.startswith("Derive"):
        st.caption("Invert one simulated T_sim map at a known density to recover "
                   "L = ϱ_sim / (100·ρ_sim).  Reuses the analytic flux model.")
        _c1, _c2, _c3 = st.columns([3, 1, 1])
        tsim_path = _c1.text_input("T_sim file path", key="di_Lsim_path",
                                   value=st.session_state.get("di_Lsim_path", ""))
        tsim_up   = _c2.file_uploader("Upload", type=["dat", "txt"], key="di_Lsim_up")
        rho_sim   = _c3.number_input("ρ_sim", 0.3, 6.0, 2.65, 0.05, key="di_rho_sim")
        if st.button("Build L from T_sim", key="di_build_Lsim"):
            try:
                src_f = _stash_upload(tsim_up) if tsim_up is not None else tsim_path
                if not src_f:
                    st.warning("Provide a T_sim file.")
                    return
                conv_z = st.session_state.get("di_Lsim_zenith", False)
                az_s, el_s, T_sim, _ = da.load_generic_map(
                    src_f, angle_is_zenith=conv_z, value_col=2)
                if not _grid_matches(az_c, el_c, az_s, el_s):
                    st.error("T_sim grid differs from the measured-data grid.")
                    return
                L_map = da.path_length_from_tsim(T_sim, float(rho_sim), el_s,
                                                 model=model, altitude_m=alt)
                st.session_state["di_L"]      = L_map
                st.session_state["di_L_desc"] = f"T_sim @ ρ={rho_sim:.2f}"
                st.success(f"✅ L map built — L ∈ "
                           f"[{np.min(L_map[L_map>0]) if np.any(L_map>0) else 0:.0f}, "
                           f"{np.max(L_map):.0f}] m")
            except Exception as _e:
                import traceback
                st.error(f"Failed: {_e}")
                st.code(traceback.format_exc())
        st.checkbox("T_sim file uses zenith angle", key="di_Lsim_zenith")

    elif src.startswith("From a terrain"):
        st.caption("Load the terrain engine's overburden_map.dat "
                   "(az, zenith, overburden g/cm², open_sky).  L = ϱ / (100·ρ).")
        _c1, _c2 = st.columns([3, 1])
        ob_path = _c1.text_input("Overburden file path", key="di_ob_path",
                                 value=st.session_state.get("di_ob_path", ""))
        ob_up   = _c2.file_uploader("Upload", type=["dat", "txt"], key="di_ob_up")
        if st.button("Build L from overburden", key="di_build_ob"):
            try:
                src_f = _stash_upload(ob_up) if ob_up is not None else ob_path
                if not src_f:
                    st.warning("Provide an overburden file.")
                    return
                az_o, el_o, L_map, _ = da.load_overburden_as_L(src_f)
                if not _grid_matches(az_c, el_c, az_o, el_o):
                    st.error("Overburden grid differs from the measured-data grid.")
                    return
                st.session_state["di_L"]      = L_map
                st.session_state["di_L_desc"] = "overburden map"
                st.success(f"✅ L map built — L ∈ "
                           f"[{np.min(L_map[L_map>0]) if np.any(L_map>0) else 0:.0f}, "
                           f"{np.max(L_map):.0f}] m")
            except Exception as _e:
                import traceback
                st.error(f"Failed: {_e}")
                st.code(traceback.format_exc())

    else:
        st.info("Opacity-only mode — the result is column density ϱ [g/cm²]; "
                "no path length is used and no density map is produced.")
        st.session_state["di_L"]      = None
        st.session_state["di_L_desc"] = "none (opacity only)"

    L_map = st.session_state.get("di_L")
    if L_map is not None and az_c is not None and src and not src.startswith("Opacity"):
        st.plotly_chart(
            _az_el_heatmap(az_c, el_c, np.where(L_map > 0, L_map, np.nan),
                           "Path length L(az, el)", "L [m]",
                           colorscale="turbo", height=300))


def _di_run_and_results(da, method_key):
    """Step ③ — run the opacity inversion and show results."""
    st.markdown("**③ Run inversion**")

    T_2d    = st.session_state.get("di_T")
    az_c    = st.session_state.get("di_az")
    el_c    = st.session_state.get("di_el")
    sigma_T = st.session_state.get("di_sigma")
    L_map   = st.session_state.get("di_L")
    L_desc  = st.session_state.get("di_L_desc", "—")
    model   = st.session_state.get("di_model", "reyna_bugaev")
    alt     = st.session_state.get("di_alt", 0.0)

    if T_2d is None:
        st.error("Load the measured data in step ① first.", icon="❌")
        return

    opac_only = st.session_state.get("di_L_source", "").startswith("Opacity")
    if L_map is None and not opac_only:
        st.warning("Configure a path-length source in step ② "
                   "(or choose opacity-only mode).", icon="⚠️")

    if st.button("▶  Run direct inversion", key="di_run", type="primary",
                 width='stretch'):
        with st.status("Inverting opacity…", expanded=True) as box:
            try:
                st.write(f"Flux model: {model}  |  altitude {alt:.0f} m  |  "
                         f"L source: {L_desc}")
                opac, sig_opac, status = da.invert_opacity_map(
                    T_2d, el_c, model=model, altitude_m=alt, sigma_T=sigma_T)
                rho_map = sig_rho = None
                if L_map is not None:
                    rho_map, sig_rho = da.opacity_to_density(opac, L_map, sig_opac)
                st.session_state["di_opac"]   = opac
                st.session_state["di_sigopac"] = sig_opac
                st.session_state["di_status"] = status
                st.session_state["di_rho"]    = rho_map
                st.session_state["di_sigrho"] = sig_rho
                box.update(label="Inversion complete!", state="complete")
                n_ok  = int(np.sum(status == 0))
                n_sky = int(np.sum(status == 1))
                msg = f"OK: {n_ok}  |  Open sky: {n_sky}"
                if rho_map is not None:
                    rr = rho_map[np.isfinite(rho_map)]
                    if rr.size:
                        msg += f"  |  Median ρ̄ = {np.median(rr):.3f} g/cm³"
                st.write(msg)
            except Exception as _e:
                import traceback
                box.update(label="Inversion failed!", state="error")
                st.error(f"{_e}")
                st.code(traceback.format_exc())

    opac = st.session_state.get("di_opac")
    if opac is None:
        return

    status  = st.session_state.get("di_status")
    rho_map = st.session_state.get("di_rho")
    sig_rho = st.session_state.get("di_sigrho")
    sig_opac = st.session_state.get("di_sigopac")

    st.divider()
    # Metrics
    n_ok  = int(np.sum(status == 0))
    n_sky = int(np.sum(status == 1))
    _m1, _m2, _m3 = st.columns(3)
    op_ok = opac[status == 0]
    _m1.metric("Median opacity (OK)",
               f"{np.nanmedian(op_ok):.0f} g/cm²" if op_ok.size else "—")
    if rho_map is not None:
        rr = rho_map[np.isfinite(rho_map)]
        _m2.metric("Median ρ̄", f"{np.median(rr):.3f} g/cm³" if rr.size else "—")
    else:
        _m2.metric("Median ρ̄", "opacity-only")
    _m3.metric("OK / open-sky bins", f"{n_ok} / {n_sky}")

    # Heatmaps
    _t1, _t2, _t3 = st.tabs(["🗺️ Maps", "📊 Status", "💾 Download"])
    with _t1:
        _ca, _cb = st.columns(2)
        with _ca:
            st.plotly_chart(
                _az_el_heatmap(az_c, el_c, np.where(np.isfinite(opac), opac, np.nan),
                               "Opacity  ϱ̂(az, el)", "ϱ [g/cm²]",
                               colorscale="inferno"))
        with _cb:
            if rho_map is not None:
                st.plotly_chart(
                    _az_el_heatmap(az_c, el_c, rho_map, "Mean density  ρ̄(az, el)",
                                   "ρ̄ [g/cm³]", colorscale="plasma",
                                   zmin=1.0, zmax=3.2))
            else:
                st.info("Density map requires a path length (step ②).")
        # Uncertainty row
        if rho_map is not None and sig_rho is not None and np.any(np.isfinite(sig_rho)):
            st.plotly_chart(
                _az_el_heatmap(az_c, el_c, sig_rho, "Density uncertainty  σ_ρ̄",
                               "σ_ρ̄ [g/cm³]", colorscale="viridis", height=300))
        elif sig_opac is not None and np.any(np.isfinite(sig_opac)):
            st.plotly_chart(
                _az_el_heatmap(az_c, el_c, sig_opac, "Opacity uncertainty  σ_ϱ",
                               "σ_ϱ [g/cm²]", colorscale="viridis", height=300))

    with _t2:
        st.plotly_chart(_status_heatmap(az_c, el_c, status))
        st.caption("green=OK · grey=open sky · red=beyond CSDA range "
                   "(opacity clamped to table maximum).")

    with _t3:
        if st.button("Prepare opacity_map.dat", key="di_prep_dl", type="primary"):
            try:
                import tempfile, os
                with tempfile.NamedTemporaryFile(mode="w", suffix=".dat",
                                                 delete=False) as tmp:
                    tpath = tmp.name
                da.write_opacity_map(
                    az_c, el_c, opac, sig_opac, status, tpath,
                    rho_map=rho_map, sigma_rho=sig_rho,
                    metadata={"model": model, "altitude_m": alt,
                              "L_source": L_desc})
                with open(tpath, "rb") as fh:
                    st.session_state["di_dl"] = fh.read()
                os.unlink(tpath)
                st.success("Ready.")
            except Exception as _e:
                st.error(f"Write failed: {_e}")
        if st.session_state.get("di_dl") is not None:
            st.download_button("⬇ Download opacity_map.dat",
                               data=st.session_state["di_dl"],
                               file_name="opacity_map.dat", mime="text/plain",
                               key="di_dl_btn")


def _render_direct_workflow(da, method_key):
    """Methods 2 & 3 — analytical opacity inversion from a single measurement."""
    if method_key == "direct":
        st.markdown(
            "**Direct opacity inversion** — recover column density "
            "ϱ = ∫ρ·dl from one measured transmission map "
            "T = Φ_target/Φ_open using an analytic flux model + CSDA range, "
            "then ρ̄ = ϱ/L.  No T_sim density library needed.")
    else:
        st.markdown(
            "**Two-flux-map ratio** — same physics as direct inversion, but the "
            "transmission is formed here from two raw flux maps "
            "(open-sky + target).")

    # Flux-model config (shared, sticky across steps)
    _f1, _f2 = st.columns([3, 1])
    labels = da.flux_model_labels()
    keys   = list(labels.keys())
    sel    = _f1.selectbox("Sea-level flux model", keys,
                           index=keys.index(st.session_state.get("di_model", keys[0]))
                                 if st.session_state.get("di_model", keys[0]) in keys else 0,
                           format_func=lambda k: labels[k], key="di_model_sel")
    st.session_state["di_model"] = sel
    st.session_state["di_alt"] = _f2.number_input("Altitude [m]", 0.0, 6000.0,
                                                  float(st.session_state.get("di_alt", 0.0)),
                                                  10.0, key="di_alt_in")
    st.caption("Flux models are azimuth-symmetric → T(ϱ) depends only on elevation.")
    st.divider()

    _di_load_input(da, method_key)
    st.divider()
    _di_path_length(da)
    st.divider()
    _di_run_and_results(da, method_key)


# ─────────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────────

def render_density_analysis_tab():
    """Render the full Density Analysis tab."""

    # ── Header ────────────────────────────────────────────────────────────────
    st.markdown(
        "**Density Analysis** — Recover rock density ρ̂(az, el) from muon "
        "transmission.  Pick an inversion method below."
    )
    st.caption(f"ucmuon_density_analysis v{_DA_VERSION}")

    # ── Load engine ───────────────────────────────────────────────────────────
    try:
        da = _load_da_engine()
    except FileNotFoundError as _fe:
        st.error(f"Cannot find ucmuon_density_analysis.py: {_fe}")
        return
    except Exception as _e:
        import traceback
        st.error(f"Failed to import ucmuon_density_analysis: {_e}")
        st.code(traceback.format_exc())
        return

    # ── Method selector ───────────────────────────────────────────────────────
    _METHODS = {
        "library":   "📚 Library matching — fit a T_sim(ρ) density library",
        "direct":    "🎯 Direct opacity inversion — single open-sky + target",
        "fluxratio": "⚖️ Two-flux-map ratio — open-sky & target flux maps",
    }
    method = st.radio(
        "Inversion method",
        list(_METHODS.keys()),
        format_func=lambda k: _METHODS[k],
        key="da_method",
        horizontal=False,
    )

    with st.expander("ℹ️ Which method should I use?", expanded=False):
        st.markdown(
            "- **Library matching** — most robust; folds in your full forward "
            "transport sim, but needs 3–5 terrain runs at different densities to "
            "build the T_sim(ρ) lookup.\n"
            "- **Direct opacity inversion** — the classic muography estimate. "
            "From one measured transmission map T = Φ_target/Φ_open it analytically "
            "recovers the opacity ϱ = ∫ρ·dl (flux model + CSDA range), then "
            "ρ̄ = ϱ/L using a geometric path length L.  No density library needed.\n"
            "- **Two-flux-map ratio** — identical physics to the direct method, "
            "but you supply the open-sky and target flux maps separately and the "
            "ratio is formed for you.\n\n"
            "A single detector resolves only the *mean* density per direction. "
            "Resolving density *along* the line of sight (3-D tomography) needs "
            "several detector viewpoints."
        )

    st.divider()

    if method == "library":
        _render_library_workflow(da)
    else:
        _render_direct_workflow(da, method)


def _render_library_workflow(da):
    """Method 1 — the original T_sim(ρ) library matching workflow."""
    st.markdown(
        "Workflow: **① Build T_sim library** (Terrain tab, multiple densities) → "
        "**② Load T_data** (measured or synthetic) → "
        "**③ Run inversion** → **④ Inspect results**"
    )

    # ── Status strip ──────────────────────────────────────────────────────────
    _lib   = st.session_state.get("da_tsim_lib")
    _tdata = st.session_state.get("da_T_data")
    _rmap  = st.session_state.get("da_rho_map")

    _s1, _s2, _s3, _s4 = st.columns(4)
    _s1.metric(
        "① T_sim library",
        f"{len(_lib)} maps" if _lib else "—",
        delta=(f"ρ [{min(_lib):.2f}–{max(_lib):.2f}]" if _lib else "not loaded"),
        delta_color="normal" if _lib else "off",
    )
    _s2.metric(
        "② T_data",
        f"{int(np.sum(~np.isnan(_tdata)))} px" if _tdata is not None else "—",
        delta="loaded" if _tdata is not None else "not loaded",
        delta_color="normal" if _tdata is not None else "off",
    )
    _s3.metric(
        "③ Inversion",
        "done" if _rmap is not None else "—",
        delta="complete" if _rmap is not None else "pending",
        delta_color="normal" if _rmap is not None else "off",
    )
    _s4.metric(
        "④ Results",
        "ready" if _rmap is not None else "—",
        delta="view below" if _rmap is not None else "run inversion first",
        delta_color="normal" if _rmap is not None else "off",
    )

    st.divider()

    # ── Main tabs ─────────────────────────────────────────────────────────────
    _tab1, _tab2, _tab3, _tab4 = st.tabs([
        "① T_sim Library",
        "② Measured Data",
        "③ Run Inversion",
        "④ Results",
    ])

    with _tab1:
        _tab_library(da)

    with _tab2:
        _tab_tdata(da)

    with _tab3:
        _tab_inversion(da)

    with _tab4:
        _tab_results(da)
