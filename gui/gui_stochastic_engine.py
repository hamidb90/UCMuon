# gui/gui_stochastic_engine.py  —  UCLouvain Muography Group
#
# Provides four public symbols consumed by gui/ucmuon_gui.py:
#
#   stochastic_available()              → (ok:bool, version_str:str)
#   render_stochastic_settings()        → renders stochastic settings widgets,
#                                         returns (v_cut, n_steps, ms_enable)
#   render_backward_mc_tab(script_dir)  → full self-contained Backward MC panel
#   build_stochastic_stdin(cfg)         → stdin string for ucmuon_stochastic_driver.py
#
# Design notes
# ────────────
# The UCMuon-MC engine follows the same integration pattern as PROPOSAL:
#   • Shared widgets (material, depth, infile) are owned by ucmuon_gui.py
#   • This module only renders the engine-specific settings block
#   • start_run() is called from ucmuon_gui.py using build_stochastic_stdin()
#
# Author: Hamid Basiri <hamid.basiri@uclouvain.be>
# MIT License 2026

import sys
import importlib.util
import time
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from pathlib import Path

_STOCHASTIC_ENGINE_VERSION = "2.0.0"

DARK_LAYOUT = dict(
    paper_bgcolor="rgb(15,17,23)", plot_bgcolor="rgb(20,22,30)",
    font=dict(color="white", size=11),
)

# ─────────────────────────────────────────────────────────────────────────────
# Availability check
# ─────────────────────────────────────────────────────────────────────────────

def stochastic_available():
    """
    Return (ok, version_str).
    Checks that the two driver scripts exist alongside this file (gui/).
    """
    gui_dir  = Path(__file__).resolve().parent
    fwd_ok   = (gui_dir / "ucmuon_stochastic_driver.py").exists()
    back_ok  = (gui_dir / "ucmuon_backward_mc.py").exists()
    if fwd_ok and back_ok:
        return True, _STOCHASTIC_ENGINE_VERSION
    missing = []
    if not fwd_ok:  missing.append("ucmuon_stochastic_driver.py")
    if not back_ok: missing.append("ucmuon_backward_mc.py")
    return False, f"Missing in gui/: {', '.join(missing)}"


# ─────────────────────────────────────────────────────────────────────────────
# Stochastic engine material constants  (mirrors MUSIC_MATERIALS; mat_id matches driver _MAT_DB)
# ─────────────────────────────────────────────────────────────────────────────
_STOCHASTIC_MAT_ID = {
    # MUSIC_MATERIALS key  →  stochastic _MAT_DB mat_id
    # 1, 2, 4, 6 use native embedded PDG per-process tables in the driver;
    # 3 (seawater) rescales the water table, 5 (custom) the rock table.
    "Standard Rock": 1, "Limestone": 1, "Rock Salt": 1,
    "Water": 2, "Ice": 6,
    "Seawater": 3,
    "Iron": 4,
    "Custom": 5,
}

# Per-material b_rad [cm²/g] for the stochastic driver (Groom 2001 Table V);
# only used for custom-material rescaling — native tables ignore it.
_STOCHASTIC_BRAD = {1: 3.475e-6, 2: 3.20e-6, 3: 3.22e-6, 4: 4.06e-6, 6: 3.20e-6}


# ─────────────────────────────────────────────────────────────────────────────
# Engine-specific settings block (renders inside the existing if/elif chain)
# ─────────────────────────────────────────────────────────────────────────────

def render_stochastic_settings():
    """
    Render the UCMuon-MC engine settings widgets.

    Called from within the Tab 2 `if True:` block, after the shared
    material/depth widgets are already rendered.

    Returns
    -------
    v_cut        : float  — catastrophic energy-loss threshold
    n_steps      : int    — transport steps (0 = auto, per-muon adaptive)
    ms_enable    : bool   — Highland MS on/off
    range_table  : int    — 0=Groom 2001  1=PDG 2024 (legacy; v2 is PDG-anchored)
    hard_spectrum: int    — 2=per-process (v2)  1=(1-v)/v (BH)  0=1/v (Groom)
    n_workers    : int    — parallel worker processes (0 = auto, 1 = serial)
    delta_rays   : bool   — explicit δ-ray straggling (T > 10 MeV)
    """
    _ps1, _ps2, _ps3 = st.columns(3)

    with _ps1:
        v_cut = st.number_input(
            "v_cut — catastrophic threshold",
            min_value=0.005, max_value=0.5, value=0.05, step=0.005,
            key="stochastic_vcut",
            help=(
                "Minimum fractional energy loss v = ΔE/E for a radiative event "
                "to be sampled stochastically (brem + pair + photonuclear).\n\n"
                "Events below v_cut contribute to the deterministic mean loss.\n\n"
                "Recommended: **0.02–0.10**.  "
                "Lower = more accurate at the cost of slower runtime."
            ),
        )

    with _ps2:
        n_steps = int(st.number_input(
            "Integration steps (0 = auto)",
            min_value=0, max_value=5000, value=0, step=100,
            key="stochastic_nsteps",
            help=(
                "Number of transport steps per muon.  **0 = adaptive auto** "
                "(≥300, scales with opacity).  Increase for higher accuracy "
                "at deep overburden (>1000 m.w.e.)."
            ),
        ))

    with _ps3:
        ms_enable = st.checkbox(
            "🔀 Highland MS deflections", value=True,
            key="stochastic_ms_enable",
            help=(
                "Apply Highland (1975) multiple Coulomb scattering per step.\n\n"
                "θ₀ = (13.6 MeV / βp) × √(x/X₀) × [1 + 0.038 ln(x/X₀)]\n\n"
                "~5% runtime overhead.  Disable only for sensitivity studies."
            ),
        )

    _pa, _pb, _pc = st.columns(3)

    with _pa:
        range_table = st.selectbox(
            "Range / dE/dx table",
            options=[1, 0],
            format_func=lambda x: "PDG 2024 (recommended)" if x == 1 else "Groom 2001",
            index=0,
            key="stochastic_range_table",
            help=(
                "CSDA range and dE/dx table used for continuous energy loss and the "
                "CSDA pre-filter.\n\n"
                "**PDG 2024** — 56 entries, direct dE/dx column (no differentiation "
                "noise), same Standard Rock parametrisation.\n\n"
                "**Groom 2001** — original 33-entry table; kept for reproducibility."
            ),
        )

    with _pb:
        hard_spectrum = st.selectbox(
            "Hard-event spectrum",
            options=[2, 1, 0],
            format_func=lambda x: {
                2: "Per-process: brems+pair+photonuc (v2)",
                1: "Bethe-Heitler (1-v)/v",
                0: "1/v (Groom)",
            }[x],
            index=0,
            key="stochastic_hard_spectrum",
            help=(
                "Spectral shape for stochastic (catastrophic) radiative events "
                "above v_cut.\n\n"
                "**Per-process (v2, recommended)** — bremsstrahlung (1−v)/v, "
                "pair production 1/v³, photonuclear 1/v, each Poisson-sampled "
                "with its own energy-dependent rate from the PDG 2024 "
                "per-process loss tables.\n\n"
                "**Bethe-Heitler (1-v)/v** — single shape; v1 default.\n\n"
                "**1/v (Groom)** — flat log-v spectrum; original behaviour."
            ),
        )
        delta_rays = st.checkbox(
            "⚡ δ-ray straggling", value=True,
            key="stochastic_delta_rays",
            help=(
                "Explicit knock-on electrons (δ-rays) with T > 10 MeV sampled "
                "from the Rutherford 1/T² spectrum; the continuous ionisation "
                "term is restricted accordingly, so the mean dE/dx is "
                "unchanged.  Adds Landau-like ionisation straggling "
                "(the dominant fluctuation below ~100 GeV)."
            ),
        )

    with _pc:
        n_workers = int(st.number_input(
            "CPU workers (0 = auto)",
            min_value=0, max_value=64, value=0, step=1,
            key="stochastic_n_workers",
            help=(
                "Parallel transport across worker processes — muons are "
                "independent, so speedup is near-linear.\n\n"
                "**0 = auto**: one worker per ~20k muons, up to all cores.\n\n"
                "**1 = serial**: bit-identical to previous releases for a "
                "given seed.  Results are reproducible for a fixed "
                "(seed, workers) pair."
            ),
        ))

    return v_cut, n_steps, ms_enable, range_table, hard_spectrum, n_workers, delta_rays


# ─────────────────────────────────────────────────────────────────────────────
# stdin builder for ucmuon_stochastic_driver.py
# ─────────────────────────────────────────────────────────────────────────────

def build_stochastic_stdin(cfg):
    """
    Build stdin string for ucmuon_stochastic_driver.py.

    cfg keys consumed:
        infile, outfile, depth_m, rho, rad (=X0_cm), mat_id (_MAT_DB key),
        transport_all, ncols, stochastic_n_steps, stochastic_v_cut,
        stochastic_ms_enable, stochastic_range_table, stochastic_hard_spectrum,
        [custom only] stochastic_Z, stochastic_A, stochastic_I_eV, stochastic_b_rad
    """
    mat_id  = cfg.get("stochastic_mat_id", 1)
    b_rad   = _STOCHASTIC_BRAD.get(mat_id, 3.475e-6)

    rho_val = cfg.get("rho", 2.65)
    rad_val = cfg.get("rad", 26.48)
    # GUI stores X₀ in g/cm² for every material (MUSIC's native unit);
    # the driver expects X₀ in cm.
    x0_cm = rad_val / max(rho_val, 1e-3)   # g/cm² → cm

    # Params 12-15: custom material (always emitted so params 16/17 land correctly)
    Z_eff  = cfg.get("stochastic_Z",     11.0)
    A_eff  = cfg.get("stochastic_A",     22.0)
    I_eV   = cfg.get("stochastic_I_eV", 136.4)
    b_cust = cfg.get("stochastic_b_rad", b_rad)

    lines = [
        cfg["infile"],
        cfg["outfile"],
        str(cfg["depth_m"]),
        str(rho_val),
        f"{x0_cm:.6f}",                                     # X₀ [cm]
        str(mat_id),
        str(int(cfg.get("transport_all", 0))),
        str(int(cfg.get("ncols", 13))),
        str(cfg.get("stochastic_n_steps", 0)),
        str(cfg.get("stochastic_v_cut",   0.05)),
        str(int(cfg.get("stochastic_ms_enable", 1))),
        str(Z_eff),                                         # param 12
        str(A_eff),                                         # param 13
        str(I_eV),                                          # param 14
        str(b_cust),                                        # param 15
        str(int(cfg.get("stochastic_range_table",   1))),  # param 16: 0=groom2001 1=pdg2024
        str(int(cfg.get("stochastic_hard_spectrum", 1))),  # param 17: 0=1/v 1=bh
        str(int(cfg.get("stochastic_seed",         42))),  # param 18: RNG seed
        str(int(cfg.get("stochastic_n_workers",     0))),  # param 19: 0=auto 1=serial
        str(int(cfg.get("stochastic_delta_rays",    1))),  # param 20: δ-ray straggling
    ]

    return "\n".join(lines) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Backward MC flux integrator — full self-contained panel
# ─────────────────────────────────────────────────────────────────────────────

_SPEC_LABELS = {
    1: "CosmoALEPH  (dN/dp ∝ p⁻³·¹⁹⁵²)",
    2: "Power-law                  (dN/dE ∝ E⁻³·⁷)",
    3: "Guan et al. (2015)",
    4: "Frosin et al. (2025)",
}


def render_backward_mc_tab(script_dir):
    """
    Full self-contained Backward MC panel.  No input muon file required.
    Physics is imported directly from gui/ucmuon_backward_mc.py.

    Parameters
    ----------
    script_dir : Path   — gui/ directory (where ucmuon_backward_mc.py lives)
    """
    st.info(
        "🔁  **Backward MC Flux Integrator** — Computes the expected muon flux "
        "at a detector sited at depth *d* without a surface muon file.  "
        "Uses backward CSDA energy inversion (E_det → E_surf) weighted by the surface "
        "spectrum, with a stochastic correction for catastrophic radiative losses.",
        icon="🔁",
    )

    # ── Parameters ────────────────────────────────────────────────────────────
    st.markdown("##### Detector & overburden")
    _bc1, _bc2 = st.columns(2)

    # Sync depth from main transport panel if available
    _default_depth = float(st.session_state.get("music_depth_m", 500.0))
    _default_rho   = float(st.session_state.get("music_rho", 2.65))

    with _bc1:
        bmc_depth = st.number_input(
            "Overburden depth [m]", 10.0, 10000.0, _default_depth, 10.0, key="bmc_depth",
        )
        bmc_rho = st.number_input(
            "Density ρ [g/cm³]", 0.1, 20.0, _default_rho, 0.05, key="bmc_rho",
        )
        bmc_mat_id = st.selectbox(
            "Material",
            options=[1, 2, 3, 4],
            format_func=lambda x: {1:"Standard Rock",2:"Water/Ice",3:"Seawater",4:"Iron"}[x],
            key="bmc_mat_id",
        )
        X_vert = bmc_depth * 100.0 * bmc_rho
        st.metric("Vertical opacity X = ρ·d", f"{X_vert:.0f} g/cm²")

    with _bc2:
        bmc_spec = st.radio(
            "Surface spectrum model",
            options=list(_SPEC_LABELS.keys()),
            format_func=lambda k: _SPEC_LABELS[k],
            key="bmc_spec",
        )
        bmc_theta_max = st.slider("θ_max [°]", 10, 80, 70, 5, key="bmc_theta_max")
        bmc_mode = st.radio(
            "Survival probability",
            options=[1, 0],
            format_func=lambda x: {
                1: "CSDA + stochastic P_surv (recommended)",
                0: "CSDA only  (P_surv = 1)",
            }[x],
            key="bmc_mode",
            help=(
                "Stochastic correction: P_surv = exp(−λ_fatal × X_slant)  "
                "where λ_fatal accounts for catastrophic radiative events that "
                "stop a muon even when mean CSDA loss says it survives.  "
                "Reduces to 1 in the CSDA limit."
            ),
        )

    with st.expander("Advanced settings", expanded=False):
        _adv1, _adv2, _adv3 = st.columns(3)
        bmc_nE       = int(_adv1.number_input("Energy bins", 20, 300, 60, 10, key="bmc_nE"))
        bmc_nth      = int(_adv2.number_input("Zenith bins",  10, 100, 25, 5,  key="bmc_nth"))
        bmc_vcut     = _adv3.number_input("v_cut", 0.01, 0.5, 0.05, 0.01, key="bmc_vcut")
        _adv4, _adv5 = st.columns(2)
        bmc_Emin = _adv4.number_input("E_min [GeV]",   0.01, 1000.0, 0.5,    0.1,   key="bmc_Emin")
        bmc_Emax = _adv5.number_input("E_max [GeV]",  100.0, 1e5,    5000.0, 100.0, key="bmc_Emax")
        bmc_outfile  = st.text_input(
            "Results output file", "backward_mc_results.dat", key="bmc_outfile",
            help="Written to the project root directory (CWD).",
        )

    # ── Run ───────────────────────────────────────────────────────────────────
    if st.button("▶  Compute Backward MC Flux", type="primary",
                 width='stretch', key="bmc_run_btn"):

        bmc_path = script_dir / "ucmuon_backward_mc.py"
        if not bmc_path.exists():
            st.error(f"❌  `ucmuon_backward_mc.py` not found in `{script_dir}`.")
            return

        with st.spinner("Running backward MC integration…"):
            try:
                spec = importlib.util.spec_from_file_location(
                    "ucmuon_backward_mc", str(bmc_path)
                )
                bmc_mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(bmc_mod)

                t0  = time.time()
                res = bmc_mod.backward_mc_flux(
                    depth_m       = bmc_depth,
                    rho           = bmc_rho,
                    mat_id        = int(bmc_mat_id),
                    spectrum_mode = int(bmc_spec),
                    E_min_GeV     = float(bmc_Emin),
                    E_max_GeV     = float(bmc_Emax),
                    theta_max_deg = float(bmc_theta_max),
                    n_E           = bmc_nE,
                    n_theta       = bmc_nth,
                    mode          = int(bmc_mode),
                    v_cut         = float(bmc_vcut),
                )
                elapsed = time.time() - t0

                st.session_state["bmc_result"]  = res
                st.session_state["bmc_elapsed"] = elapsed
                st.session_state["bmc_outfile_path"] = bmc_outfile

                bmc_mod._write_results(res, bmc_outfile)
                st.success(
                    f"✅  Computed in {elapsed:.1f} s  |  "
                    f"Expected rate: **{res['rate_m2_s']:.4e} m⁻² s⁻¹**"
                )
            except Exception as _e:
                import traceback
                st.error(f"❌  Computation failed: {_e}")
                st.code(traceback.format_exc())

    # ── Results ───────────────────────────────────────────────────────────────
    res = st.session_state.get("bmc_result")
    if res is None:
        st.info("Press **Compute** to run the backward MC integration.")
        return

    info    = res["info"]
    elapsed = st.session_state.get("bmc_elapsed", 0)

    _m1, _m2, _m3, _m4 = st.columns(4)
    _m1.metric("Expected rate",    f"{res['rate_m2_s']:.3e} m⁻² s⁻¹")
    _m2.metric("Vertical opacity", f"{res['X_vert_gcm2']:.0f} g/cm²")
    imax = int(np.argmax(res["flux_det"]))
    _m3.metric("Peak E (detector)", f"{res['E_det_GeV'][imax]:.1f} GeV")
    _ps = res["P_survival"]
    _m4.metric("Mean P_survival",  f"{float(_ps[_ps > 1e-4].mean()):.4f}" if (_ps > 1e-4).any() else "—")

    # Spectrum plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=res["E_det_GeV"], y=res["flux_det"], mode="lines",
        name=f"Detector  (depth {info['depth_m']:.0f} m, {info['mat']})",
        line=dict(color="#00b4d8", width=2.5),
    ))
    mask = res["E_surf_GeV"] > 0
    if mask.any():
        fig.add_trace(go.Scatter(
            x=res["E_surf_GeV"][mask], y=res["flux_surf"][mask], mode="lines",
            name="Surface (reference)",
            line=dict(color="#ffd700", width=1.5, dash="dash"),
        ))
    fig.update_layout(
        **DARK_LAYOUT, height=380,
        xaxis=dict(type="log", title="Energy [GeV]",              gridcolor="#2a2a3a"),
        yaxis=dict(type="log", title="dΦ/dE  [m⁻² s⁻¹ GeV⁻¹]", gridcolor="#2a2a3a"),
        title=dict(
            text=(f"Backward MC — Muon flux at {info['depth_m']:.0f} m depth  |  "
                  f"Spectrum {info['spectrum_mode']}  |  {info['mode']}  |  "
                  f"computed in {elapsed:.1f} s"),
            font=dict(size=12),
        ),
        legend=dict(bgcolor="rgba(0,0,0,0.5)", bordercolor="#555", borderwidth=1),
        margin=dict(l=60, r=20, t=55, b=50),
    )
    st.plotly_chart(fig, config={"displayModeBar": False})

    # Survival probability plot
    if info["mode"] != "CSDA only" and mask.any():
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            x=res["E_det_GeV"][mask], y=_ps[mask],
            mode="lines+markers", name="P_survival",
            line=dict(color="#ff9f43", width=2),
            marker=dict(size=5, color="#ff9f43"),
        ))
        fig2.add_hline(y=1.0, line=dict(color="#aaa", width=1, dash="dot"),
                       annotation_text="CSDA limit (P=1)", annotation_position="top left")
        fig2.update_layout(
            **DARK_LAYOUT, height=260,
            xaxis=dict(type="log", title="Detector energy E_det [GeV]", gridcolor="#2a2a3a"),
            yaxis=dict(title="P_survival", range=[0, 1.05],              gridcolor="#2a2a3a"),
            title=dict(text="Stochastic survival probability", font=dict(size=12)),
            margin=dict(l=60, r=20, t=40, b=50),
        )
        st.plotly_chart(fig2, config={"displayModeBar": False})

    # Energy mapping expander
    with st.expander("📊 Backward energy mapping  E_det → E_surface", expanded=False):
        fig3 = go.Figure()
        if mask.any():
            fig3.add_trace(go.Scatter(
                x=res["E_det_GeV"][mask], y=res["E_surf_GeV"][mask],
                mode="lines", name="E_surf(E_det)",
                line=dict(color="#a29bfe", width=2),
            ))
            fig3.add_trace(go.Scatter(
                x=[res["E_det_GeV"][0], res["E_det_GeV"][-1]],
                y=[res["E_det_GeV"][0], res["E_det_GeV"][-1]],
                mode="lines", name="E_surf = E_det",
                line=dict(color="#dfe6e9", dash="dot", width=1),
            ))
        fig3.update_layout(
            **DARK_LAYOUT, height=280,
            xaxis=dict(type="log", title="E_det [GeV]",  gridcolor="#2a2a3a"),
            yaxis=dict(type="log", title="E_surf [GeV]", gridcolor="#2a2a3a"),
            title=dict(text="CSDA backward energy mapping (Jacobian included in flux)",
                       font=dict(size=12)),
            margin=dict(l=60, r=20, t=40, b=50),
        )
        st.plotly_chart(fig3, config={"displayModeBar": False})

    # Download
    outpath = st.session_state.get("bmc_outfile_path", "backward_mc_results.dat")
    if Path(outpath).exists():
        with open(outpath, "rb") as _fh:
            st.download_button(
                "⬇️  Download results (ASCII, 5-col)",
                data=_fh,
                file_name=Path(outpath).name,
                mime="text/plain",
                width='stretch',
                key="bmc_dl_btn",
            )
