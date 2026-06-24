# gui_mcs_estimator.py  —  Multiple Coulomb Scattering acceptance estimator
# UCMuon — UCLouvain Muography Group
# Author : Hamid Basiri <hamid.basiri@uclouvain.be>
# License: MIT
#
# Physics: Highland–Lynch–Dahl formula (PDG Rev. Part. Phys., §34.3)
#
#   θ₀ = (13.6 MeV / βcp) · √(x/X₀) · [1 + 0.038 ln(x/X₀)]
#
#   σ_⊥  = (L / √3) · θ₀          RMS lateral displacement, one projected plane
#   σ_r  = √2 · σ_⊥               RMS radial displacement in 2D transverse plane
#
# Geometry: nominal muon impact position is computed as the closest-approach
#   point to the detector axis, clipped to the cylinder z-range.
#   (Using mid-plane projection is incorrect for oblique tracks: a muon entering
#   the top face at 40° projects to >300 cm at the mid-plane even though it
#   passes through the centre — this gives spuriously low acceptance.)
#
# Usage: called from cosmoaleph_gui.py (Tab 1) via render_mcs_panel().

from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

# ── Physical constants ─────────────────────────────────────────────────────────
_M_MU_GEV  = 0.10566           # muon mass [GeV/c²]
_HIGHLAND_K = 13.6              # [MeV] — coefficient in Highland formula
_BG         = "rgb(15,17,23)"   # Streamlit dark background

# ── Material presets: (ρ [g/cm³], X₀ [g/cm²]) ────────────────────────────────
MCS_MATERIALS: dict[str, tuple[float, float] | None] = {
    "Standard Rock":  (2.65, 26.7),
    "Limestone":      (2.71, 27.0),
    "Basalt":         (2.80, 25.5),
    "Sandstone":      (2.32, 26.1),
    "Concrete":       (2.30, 26.7),
    "Salt (halite)":  (2.16, 29.5),
    "Ice":            (0.917, 36.1),
    "Custom":         None,
}


# ── Core physics ───────────────────────────────────────────────────────────────

def highland_theta0(p_GeV: np.ndarray,
                    x_gcm2: np.ndarray,
                    X0_gcm2: float = 26.7) -> np.ndarray:
    """Highland–Lynch–Dahl RMS projected scattering angle [rad]."""
    p_GeV  = np.asarray(p_GeV,  dtype=float)
    x_gcm2 = np.asarray(x_gcm2, dtype=float)
    E_GeV  = np.sqrt(p_GeV**2 + _M_MU_GEV**2)
    beta   = p_GeV / np.maximum(E_GeV, 1e-12)
    p_MeV  = p_GeV * 1000.0
    xr     = np.maximum(x_gcm2 / X0_gcm2, 1e-12)
    theta0 = (_HIGHLAND_K / (beta * p_MeV)) * np.sqrt(xr) * (1.0 + 0.038 * np.log(xr))
    return np.clip(theta0, 0.0, np.pi)


def lateral_sigma_cm(theta0_rad: np.ndarray,
                     path_cm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    RMS lateral displacement from Highland θ₀.
    Returns (sigma_perp [cm], sigma_r [cm]) where sigma_r = √2 · sigma_perp.
    """
    sigma_perp = (path_cm / np.sqrt(3.0)) * theta0_rad
    sigma_r    = np.sqrt(2.0) * sigma_perp
    return sigma_perp, sigma_r


def _cylinder_closest_approach(x0_cm, y0_cm, cx, cy, cz,
                                ax, ay, az, bx, by, bz):
    """
    For each muon (vectorised), compute the nominal transverse distance
    from the cylinder axis at the CLOSEST APPROACH point along the track,
    clipped to the z-range of the cylinder.

    This is the correct nominal impact position for a general cylinder
    regardless of muon zenith angle.  Mid-plane projection is wrong for
    oblique tracks: a 40° muon entering the top face at r=100 cm projects
    to r>300 cm at the mid-plane, even though it passes through the centre.

    Parameters (all arrays of length n_mu, or scalars for det geometry)
    ----------
    x0_cm, y0_cm : muon surface positions [cm]
    cx, cy, cz   : muon direction cosines (unit vector, cz < 0 for downward)
    ax,ay,az     : cylinder axis bottom point [cm]
    bx,by,bz     : cylinder axis top point [cm]

    Returns
    -------
    x_nom, y_nom : nominal transverse positions relative to lab origin [cm]
                   (distance from axis = √((x_nom-det_cx)²+(y_nom-det_cy)²))
    depth_cm     : representative vertical depth [cm] (for path-length calc)
    """
    n_mu = len(x0_cm)

    # Axis unit vector
    ab    = np.array([bx - ax, by - ay, bz - az], dtype=float)
    ab_len = float(np.linalg.norm(ab))
    if ab_len < 1e-6:
        return x0_cm * 0.0, y0_cm * 0.0, abs(az)

    ahat = ab / ab_len  # unit axis direction

    # Muon track: r(t) = r0 + t*d  (r0 at z=0)
    r0 = np.stack([x0_cm, y0_cm, np.zeros(n_mu)], axis=1)   # (n, 3)
    d  = np.stack([cx, cy, cz], axis=1)                       # (n, 3)
    A  = np.array([ax, ay, az], dtype=float)

    # Vector from axis origin A to muon surface position
    p = r0 - A[np.newaxis, :]                                 # (n, 3)

    # Projections along axis
    p_da = p @ ahat      # (n,)
    d_da = d @ ahat      # (n,)

    # Perpendicular components (to axis)
    p_perp = p - p_da[:, np.newaxis] * ahat[np.newaxis, :]   # (n, 3)
    d_perp = d - d_da[:, np.newaxis] * ahat[np.newaxis, :]   # (n, 3)

    dp_sq = np.maximum(np.sum(d_perp**2, axis=1), 1e-12)     # (n,)

    # t of minimum transverse distance to axis
    t_min_tr = -np.sum(p_perp * d_perp, axis=1) / dp_sq      # (n,)

    # Clip t to the extent of the cylinder along its axis
    # Cylinder extends from axis-projection 0 to ab_len
    # At parametric t: projection = p_da + t * d_da
    with np.errstate(divide="ignore", invalid="ignore"):
        d_da_safe = np.where(np.abs(d_da) > 1e-9, d_da, np.nan)
        t_at_0    = -p_da / d_da_safe           # t where projection = 0
        t_at_L    = (ab_len - p_da) / d_da_safe  # t where projection = ab_len

    t_lo = np.where(np.isfinite(t_at_0) & np.isfinite(t_at_L),
                    np.minimum(t_at_0, t_at_L), 0.0)
    t_hi = np.where(np.isfinite(t_at_0) & np.isfinite(t_at_L),
                    np.maximum(t_at_0, t_at_L), 1e8)

    t_nom = np.clip(t_min_tr, t_lo, t_hi)
    t_nom = np.maximum(t_nom, 0.0)   # only forward propagation

    # Nominal 3D position at t_nom
    r_nom = r0 + t_nom[:, np.newaxis] * d  # (n, 3)

    # Transverse components (displacement from axis at the closest point)
    r_from_A  = r_nom - A[np.newaxis, :]                        # (n, 3)
    along_ax  = (r_from_A @ ahat)[:, np.newaxis] * ahat        # (n, 3) along axis
    r_transv  = r_from_A - along_ax                             # (n, 3) perpendicular

    # Nominal impact in lab frame (add back axis origin)
    x_nom = r_transv[:, 0] + ax
    y_nom = r_transv[:, 1] + ay

    # Representative depth = midpoint of cylinder along z
    depth_cm = abs((az + bz) / 2.0)
    if depth_cm < 1.0:
        depth_cm = abs(az)

    return x_nom, y_nom, depth_cm


# ── Main estimator ─────────────────────────────────────────────────────────────

def estimate_mcs_acceptance(df: pd.DataFrame,
                             detectors: list[dict],
                             rho: float,
                             X0: float,
                             n_mc: int = 500,
                             seed: int = 42) -> dict | None:
    """
    Per-muon Monte Carlo MCS acceptance estimate using the Highland formula.

    For each selected muon the nominal trajectory (already confirmed to hit
    the detector by the Fortran ray-caster) is displaced by a 2-D Gaussian
    representing the MCS smearing:

        δx, δy  ∼  𝒩(0, σ_⊥²)  independent per projected plane

    The nominal impact position is computed as the CLOSEST-APPROACH point
    to the cylinder axis, clipped to the cylinder z-range (not the mid-plane
    projection, which is incorrect for oblique tracks).

    Path length per muon:  L = D_det / |cos θ_z|   (flat-slab)
    """
    if df is None or len(df) == 0 or not detectors:
        return None

    rng = np.random.default_rng(seed)

    # ── Kinematics ─────────────────────────────────────────────────────────────
    p_arr = np.asarray(df["p"].values, dtype=float)
    E_arr = np.sqrt(p_arr**2 + _M_MU_GEV**2)

    if {"px", "py", "pz"}.issubset(df.columns):
        p_safe = np.maximum(p_arr, 1e-10)
        cx = df["px"].values.astype(float) / p_safe
        cy = df["py"].values.astype(float) / p_safe
        cz = df["pz"].values.astype(float) / p_safe
    else:
        th = np.radians(df["theta"].values.astype(float))
        ph = np.radians(df["phi"].values.astype(float))
        cx = np.sin(th) * np.cos(ph)
        cy = np.sin(th) * np.sin(ph)
        cz = -np.cos(th)

    x0_cm = df["x"].values.astype(float)
    y0_cm = df["y"].values.astype(float)

    # ── Detector geometry ───────────────────────────────────────────────────────
    det    = detectors[0]
    margin = float(det.get("margin", 0.0))
    n_mu   = len(p_arr)

    if det["shape"] == 1:   # ── Cylinder ──────────────────────────────────────
        ax_, ay_, az_ = float(det["ax"]), float(det["ay"]), float(det["az"])
        bx_, by_, bz_ = float(det["bx"]), float(det["by"]), float(det["bz"])
        R_eff = float(det["r"]) + margin

        # Nominal impact: closest approach to cylinder axis (not mid-plane!)
        x_nom, y_nom, depth_cm = _cylinder_closest_approach(
            x0_cm, y0_cm, cx, cy, cz,
            ax_, ay_, az_, bx_, by_, bz_
        )

        # Axis origin for transverse distance check
        det_cx, det_cy = ax_, ay_

    else:   # ── AABB box ──────────────────────────────────────────────────────
        xn, xx = float(det["xmin"]), float(det["xmax"])
        yn, yx = float(det["ymin"]), float(det["ymax"])
        zn, zx = float(det["zmin"]), float(det["zmax"])
        det_cx = (xn + xx) / 2.0
        det_cy = (yn + yx) / 2.0
        depth_cm = abs((zn + zx) / 2.0) or 1000.0
        hx = (xx - xn) / 2.0 + margin
        hy = (yx - yn) / 2.0 + margin

        # Project to entry face (z = zmax for downward muons)
        with np.errstate(divide="ignore", invalid="ignore"):
            t_entry = np.where(np.abs(cz) > 1e-8, zx / cz, 0.0)
        t_entry = np.maximum(t_entry, 0.0)
        x_nom = x0_cm + t_entry * cx
        y_nom = y0_cm + t_entry * cy

    # ── Path length and overburden ─────────────────────────────────────────────
    cos_zen = np.maximum(np.abs(cz), 0.02)
    path_cm = depth_cm / cos_zen
    x_gcm2  = rho * path_cm

    # ── Highland scattering ────────────────────────────────────────────────────
    theta0              = highland_theta0(p_arr, x_gcm2, X0)
    sigma_perp, sigma_r = lateral_sigma_cm(theta0, path_cm)

    # ── Monte Carlo miss probability ───────────────────────────────────────────
    dx_mc = rng.normal(0.0, 1.0, (n_mc, n_mu)) * sigma_perp[np.newaxis, :]
    dy_mc = rng.normal(0.0, 1.0, (n_mc, n_mu)) * sigma_perp[np.newaxis, :]

    x_sc = x_nom[np.newaxis, :] + dx_mc  # (n_mc, n_mu)
    y_sc = y_nom[np.newaxis, :] + dy_mc

    if det["shape"] == 1:
        r_sc    = np.sqrt((x_sc - det_cx)**2 + (y_sc - det_cy)**2)
        miss_mc = r_sc > R_eff
    else:
        miss_mc = ((np.abs(x_sc - det_cx) > hx) |
                   (np.abs(y_sc - det_cy) > hy))

    miss_prob = miss_mc.mean(axis=0)
    hit_prob  = 1.0 - miss_prob

    n_selected  = n_mu
    n_after_mcs = float(hit_prob.sum())
    acceptance  = n_after_mcs / n_selected if n_selected > 0 else 0.0

    # ── Per-energy-bin statistics ──────────────────────────────────────────────
    n_bins    = 10
    e_lo      = max(E_arr.min() * 0.95, 0.1)
    e_hi      = E_arr.max() * 1.05
    log_edges = np.logspace(np.log10(e_lo), np.log10(e_hi), n_bins + 1)

    bin_acc = np.full(n_bins, np.nan)
    bin_n   = np.zeros(n_bins, dtype=int)
    bin_sig = np.full(n_bins, np.nan)
    bin_th0 = np.full(n_bins, np.nan)

    for i in range(n_bins):
        mask = (E_arr >= log_edges[i]) & (E_arr < log_edges[i + 1])
        cnt  = int(mask.sum())
        if cnt > 0:
            bin_acc[i] = float(hit_prob[mask].mean())
            bin_n[i]   = cnt
            bin_sig[i] = float(np.median(sigma_r[mask]))
            bin_th0[i] = float(np.median(theta0[mask]) * 1e3)   # mrad

    return {
        "p_arr":                  p_arr,
        "E_arr":                  E_arr,
        "sigma_r_arr":            sigma_r,
        "theta0_arr":             theta0,
        "miss_prob_arr":          miss_prob,
        "hit_prob_arr":           hit_prob,
        "n_selected":             n_selected,
        "n_after_mcs":            n_after_mcs,
        "acceptance_frac":        acceptance,
        "depth_cm":               depth_cm,
        "path_cm_arr":            path_cm,
        "x_gcm2_arr":             x_gcm2,
        "energy_bins":            log_edges,
        "bin_acceptance":         bin_acc,
        "bin_n":                  bin_n,
        "bin_sigma_r_median":     bin_sig,
        "bin_theta0_mrad_median": bin_th0,
    }


# ── Streamlit renderer ─────────────────────────────────────────────────────────

def render_mcs_panel(load_file_fn) -> None:
    """
    Render the MCS Acceptance Estimator panel inside Tab 1.
    All gating (gen_success, gen_running, gen_use_detector) is handled
    by the caller in the GUI file where _gg() is in scope.
    """
    sel_file  = st.session_state.get("selected_file", "")
    detectors = st.session_state.get("gen_detectors", [])
    if not sel_file or not Path(sel_file).exists() or not detectors:
        return

    st.divider()
    st.markdown("#### 🎯 MCS Acceptance Estimator")
    st.caption(
        "Estimates how many selected muons still hit the detector after "
        "multiple Coulomb scattering through the overburden — "
        "Highland–Lynch–Dahl formula (PDG §34.3), flat-slab geometry."
    )

    # ── Controls row ──────────────────────────────────────────────────────────
    c1, c2, c3, c4 = st.columns([2, 1, 1, 1])

    with c1:
        mat_choice = st.selectbox(
            "Overburden material",
            list(MCS_MATERIALS.keys()),
            index=0,
            key="mcs_material",
        )
        preset  = MCS_MATERIALS[mat_choice]
        rho_def = preset[0] if preset else 2.65
        X0_def  = preset[1] if preset else 26.7

        if mat_choice == "Custom":
            cc1, cc2 = st.columns(2)
            rho = cc1.number_input("ρ [g/cm³]", 0.1, 20.0, rho_def,
                                   step=0.01, format="%.3f", key="mcs_rho_custom")
            X0  = cc2.number_input("X₀ [g/cm²]", 1.0, 100.0, X0_def,
                                   step=0.1, format="%.1f", key="mcs_X0_custom")
        else:
            rho, X0 = rho_def, X0_def
            st.caption(f"ρ = {rho:.3f} g/cm³  ·  X₀ = {X0:.1f} g/cm²")

    with c2:
        n_mc = st.select_slider(
            "MC samples / muon",
            options=[100, 200, 500, 1000],
            value=500,
            key="mcs_nmc",
            help="MC draws per muon for miss-probability estimate. 500 is sufficient.",
        )

    with c3:
        det0 = detectors[0]
        if det0["shape"] == 1:
            auto_d = abs((float(det0["az"]) + float(det0["bz"])) / 2.0) / 100.0
        else:
            auto_d = abs((float(det0["zmin"]) + float(det0["zmax"])) / 2.0) / 100.0
        if auto_d < 0.01:
            auto_d = 100.0
        st.metric("Detector depth", f"{auto_d:.1f} m",
                  help="Vertical depth to detector mid-plane (from geometry).")

    with c4:
        st.markdown("&nbsp;")
        run_mcs = st.button("▶  Estimate", key="btn_mcs_run",
                            width='stretch', type="primary")

    # ── Formula & assumptions ─────────────────────────────────────────────────
    with st.expander("ℹ️  Formula & assumptions", expanded=False):
        st.markdown("**Highland–Lynch–Dahl formula** (PDG Rev. Part. Phys. §34.3):")
        st.latex(
            r"\theta_0 = \frac{13.6\,\text{MeV}}{\beta c p}"
            r"\sqrt{\frac{x}{X_0}}"
            r"\left[1 + 0.038\,\ln\frac{x}{X_0}\right]"
        )
        st.markdown(
            "- **Path length:** L = D_det / |cos θ_z| — flat-slab (same as MUSIC/BB engines)  \n"
            "- **Lateral displacement (1 plane):** σ_⊥ = (L/√3) · θ₀  \n"
            "- **Radial displacement (2D):** σ_r = √2 · σ_⊥  \n"
            "- **Nominal position:** closest-approach to cylinder axis, clipped to "
            "cylinder z-range (corrects mid-plane projection error for oblique tracks)  \n"
            "- **MC step:** N draws of (δx, δy) ~ 𝒩(0, σ_⊥²) tested against detector  \n"
            "  \n"
            "⚠️ Flat-slab only · first detector · no energy loss during propagation"
        )

    # ── Run ───────────────────────────────────────────────────────────────────
    if run_mcs:
        with st.spinner("Computing MCS acceptance…"):
            try:
                df_sel = load_file_fn(sel_file,
                                      mtime=Path(sel_file).stat().st_mtime)
                result = estimate_mcs_acceptance(
                    df_sel, detectors, rho=rho, X0=X0, n_mc=n_mc, seed=42
                )
                st.session_state["mcs_result"]     = result
                st.session_state["mcs_mat_label"]  = mat_choice
                st.session_state["mcs_rho_used"]   = rho
                st.session_state["mcs_X0_used"]    = X0
                st.session_state["mcs_nmc_used"]   = n_mc
            except Exception as exc:
                st.error(f"❌  MCS estimator failed: {exc}")
                st.session_state.pop("mcs_result", None)
                return

    result = st.session_state.get("mcs_result")
    if result is None:
        return

    # ── Summary metrics ────────────────────────────────────────────────────────
    n_sel   = result["n_selected"]
    n_aft   = result["n_after_mcs"]
    acc     = result["acceptance_frac"]
    depth_m = result["depth_cm"] / 100.0
    sig_med = float(np.median(result["sigma_r_arr"]))
    mat_lbl = st.session_state.get("mcs_mat_label", "—")
    rho_s   = st.session_state.get("mcs_rho_used", rho)
    X0_s    = st.session_state.get("mcs_X0_used", X0)
    nmc_s   = st.session_state.get("mcs_nmc_used", n_mc)

    st.markdown("---")
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Selected muons",     f"{n_sel:,}")
    m2.metric("Expected after MCS", f"{n_aft:,.0f}",
              delta=f"−{n_sel - n_aft:.0f}", delta_color="inverse")
    m3.metric("MCS acceptance",     f"{acc*100:.1f} %",
              delta=f"−{(1-acc)*100:.1f} %", delta_color="inverse")
    m4.metric("Depth / med. σ_r",   f"{depth_m:.0f} m / {sig_med:.1f} cm",
              help="Vertical depth to detector mid-plane | Median RMS radial MCS displacement")
    m5.metric("Material",           mat_lbl)

    # Severity banner
    if acc >= 0.95:
        st.success(f"✅  MCS effect negligible — {(1-acc)*100:.1f}% loss. "
                   "Detector acceptance is robust against scattering.")
    elif acc >= 0.80:
        st.warning(f"⚠️  Moderate MCS loss ({(1-acc)*100:.1f}%). "
                   "Low-energy muons are significantly affected — consider raising E_min.")
    else:
        st.error(f"🚨  Severe MCS loss ({(1-acc)*100:.1f}%). "
                 "E_min is likely too low for this depth/geometry. "
                 "Raise E_min or increase the detector acceptance margin.")

    # ── Plots ──────────────────────────────────────────────────────────────────
    E_bins  = result["energy_bins"]
    bin_E_c = np.sqrt(E_bins[:-1] * E_bins[1:])
    bin_acc = result["bin_acceptance"]
    bin_n   = result["bin_n"]
    bin_sig = result["bin_sigma_r_median"]
    valid   = bin_n > 0

    # Detector reference size
    det0 = detectors[0]
    if det0["shape"] == 1:
        r_ref   = float(det0["r"]) + float(det0.get("margin", 0.0))
        r_label = f"R_det + margin = {r_ref:.0f} cm"
    else:
        r_ref   = (min(float(det0["xmax"]) - float(det0["xmin"]),
                       float(det0["ymax"]) - float(det0["ymin"])) / 2.0
                   + float(det0.get("margin", 0.0)))
        r_label = f"Half-width + margin = {r_ref:.0f} cm"

    def _bar_colour(v):
        if np.isnan(v):  return "rgba(80,80,80,0.4)"
        if v >= 0.95:    return "rgba(52,211,153,0.85)"
        if v >= 0.80:    return "rgba(251,191,36,0.85)"
        return                  "rgba(239,68,68,0.85)"

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Acceptance per energy bin", "Median σ_r vs energy"),
        horizontal_spacing=0.14,
    )

    # ── Left: acceptance bar chart ─────────────────────────────────────────────
    vx = [f"{E_bins[i]:.0f}–{E_bins[i+1]:.0f}"
          for i in range(len(bin_E_c)) if valid[i]]
    vy = [float(bin_acc[i]) for i in range(len(bin_E_c)) if valid[i]]
    vn = [int(bin_n[i])      for i in range(len(bin_E_c)) if valid[i]]
    vc = [_bar_colour(v)     for v in vy]
    vt = [f"{v*100:.0f}%<br>N={n}" for v, n in zip(vy, vn)]

    fig.add_trace(go.Bar(
        x=vx, y=[v * 100 for v in vy],
        marker_color=vc, text=vt,
        textposition="outside", textfont=dict(size=9),
        hovertemplate="<b>%{x} GeV</b><br>Acceptance: %{y:.1f}%<extra></extra>",
        showlegend=False,
    ), row=1, col=1)

    for thresh, col in [(95, "rgba(52,211,153,0.35)"),
                        (80, "rgba(251,191,36,0.35)")]:
        fig.add_hline(y=thresh,
                      line=dict(color=col, width=1, dash="dash"),
                      row=1, col=1)

    fig.update_yaxes(title_text="Acceptance [%]", range=[0, 115],
                     row=1, col=1, gridcolor="rgba(255,255,255,0.06)")
    fig.update_xaxes(title_text="Energy bin [GeV]", tickangle=-35, row=1, col=1)

    # ── Right: σ_r vs energy (linear scale, explicit range) ────────────────────
    # Use linear scale — range is ~1 decade so linear is cleaner than log.
    # Explicit range avoids Plotly autoscale bug with add_hline on log axes.
    vs = valid & ~np.isnan(bin_sig) & np.isfinite(bin_sig)

    if vs.any():
        sig_vals = bin_sig[vs]
        y_max_plot = float(np.max(sig_vals)) * 1.25
        y_max_plot = max(y_max_plot, r_ref * 1.1)

        fig.add_trace(go.Scatter(
            x=bin_E_c[vs], y=sig_vals,
            mode="lines+markers",
            line=dict(color="rgba(56,189,248,0.9)", width=2),
            marker=dict(size=7),
            name="σ_r [cm]",
            hovertemplate="E = %{x:.0f} GeV<br>σ_r = %{y:.1f} cm<extra></extra>",
        ), row=1, col=2)

        # Detector reference line — draw as a scatter trace to avoid log-axis bug
        e_range = [float(bin_E_c[vs][0]) * 0.8, float(bin_E_c[vs][-1]) * 1.2]
        fig.add_trace(go.Scatter(
            x=e_range, y=[r_ref, r_ref],
            mode="lines",
            line=dict(color="rgba(239,68,68,0.6)", width=1.5, dash="dash"),
            name=r_label,
            hoverinfo="skip",
        ), row=1, col=2)

        fig.update_xaxes(title_text="Energy [GeV]", type="log",
                         row=1, col=2, gridcolor="rgba(255,255,255,0.06)")
        fig.update_yaxes(title_text="Median σ_r [cm]",
                         range=[0, y_max_plot],
                         row=1, col=2, gridcolor="rgba(255,255,255,0.06)")

    fig.update_layout(
        height=380,
        margin=dict(l=60, r=20, t=50, b=80),
        paper_bgcolor=_BG,
        plot_bgcolor="rgb(20,22,30)",
        font=dict(color="white", size=11),
        showlegend=True,
        legend=dict(bgcolor="rgba(0,0,0,0)", font_size=9,
                    x=0.99, xanchor="right", y=0.99, yanchor="top"),
    )
    fig.update_annotations(font_size=12, font_color="rgba(200,200,200,0.85)")
    st.plotly_chart(fig)

    # ── Details & export ───────────────────────────────────────────────────────
    with st.expander("📋  Details & export", expanded=False):

        sigma_r_arr = result["sigma_r_arr"]
        p50 = float(np.median(sigma_r_arr))
        p95 = float(np.percentile(sigma_r_arr, 95))

        d1, d2, d3, d4 = st.columns(4)
        d1.metric("Median σ_r",          f"{p50:.1f} cm")
        d2.metric("95th-pct σ_r",        f"{p95:.1f} cm")
        d3.metric("Detector ref",         f"{r_ref:.0f} cm", help=r_label)
        d4.metric("σ_r > R_det fraction",
                  f"{(sigma_r_arr > r_ref).mean()*100:.1f} %",
                  help="Fraction of muons whose 1σ scatter circle exceeds the detector radius.")

        st.markdown("&nbsp;")
        st.markdown("**Per-bin summary**")

        valid_idx = np.where(valid)[0]
        th0_arr   = result["bin_theta0_mrad_median"]
        rows = []
        for i in valid_idx:
            rows.append({
                "Energy [GeV]":    f"{E_bins[i]:.1f} – {E_bins[i+1]:.1f}",
                "N muons":         int(bin_n[i]),
                "Acceptance [%]":  f"{bin_acc[i]*100:.1f}",
                "Med. σ_r [cm]":   f"{bin_sig[i]:.1f}" if not np.isnan(bin_sig[i]) else "—",
                "Med. θ₀ [mrad]":  f"{th0_arr[i]:.3f}" if not np.isnan(th0_arr[i]) else "—",
            })
        if rows:
            st.dataframe(pd.DataFrame(rows).set_index("Energy [GeV]"),
                         width='stretch')

        st.markdown("&nbsp;")
        dl1, dl2 = st.columns(2)

        df_exp = pd.DataFrame({
            "E_GeV":           result["E_arr"],
            "p_GeV_c":         result["p_arr"],
            "path_cm":         result["path_cm_arr"],
            "overburden_gcm2": result["x_gcm2_arr"],
            "theta0_mrad":     result["theta0_arr"] * 1e3,
            "sigma_r_cm":      sigma_r_arr,
            "hit_probability": result["hit_prob_arr"],
        })
        dl1.download_button(
            "⬇️  Per-muon CSV",
            data=df_exp.to_csv(index=False).encode(),
            file_name="ucmuon_mcs_per_muon.csv",
            mime="text/csv",
            width='stretch',
            key="dl_mcs_csv",
        )

        summary_lines = [
            "# UCMuon MCS Acceptance Summary",
            f"# Material: {mat_lbl}  rho={rho_s:.3f} g/cm3  X0={X0_s:.1f} g/cm2",
            f"# Depth: {depth_m:.1f} m  |  MC samples/muon: {nmc_s}",
            f"# N_selected={n_sel}  N_after_MCS={n_aft:.0f}  Acceptance={acc*100:.2f}%",
            "#",
            "# E_low    E_high   N      Acc_%   sig_r_cm  th0_mrad",
        ]
        for i in np.where(valid)[0]:
            sig_v = bin_sig[i] if np.isfinite(bin_sig[i]) else 0.0
            th0_v = th0_arr[i] if np.isfinite(th0_arr[i]) else 0.0
            summary_lines.append(
                f"{E_bins[i]:8.2f}  {E_bins[i+1]:8.2f}  {bin_n[i]:5d}  "
                f"{bin_acc[i]*100:6.2f}  {sig_v:9.2f}  {th0_v:9.4f}"
            )
        dl2.download_button(
            "⬇️  Summary TXT",
            data="\n".join(summary_lines).encode(),
            file_name="ucmuon_mcs_summary.txt",
            mime="text/plain",
            width='stretch',
            key="dl_mcs_txt",
        )
