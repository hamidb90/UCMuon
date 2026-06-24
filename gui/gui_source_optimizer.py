# gui_source_optimizer.py  —  Source size + optimisation panel
# UCMuon — UCLouvain Muography Group
# Author : Hamid Basiri <hamid.basiri@uclouvain.be>
# License: MIT

from __future__ import annotations
import numpy as np
import pandas as pd
import streamlit as st

_GAMMA = 3.7   # CosmoALEPH spectral index

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


def _groom(opacity_gcm2: float) -> float:
    """Min muon KE [GeV] to traverse opacity_gcm2 [g/cm²] (Groom 2001)."""
    x = float(np.clip(opacity_gcm2, _GROOM_R_GCM2[0], _GROOM_R_GCM2[-1]))
    return float(np.exp(
        np.interp(np.log(x), np.log(_GROOM_R_GCM2), np.log(_GROOM_T_MEV))
    )) / 1000.0


def _csda_range_gcm2(E_GeV: float) -> float:
    """CSDA range [g/cm²] of a muon with KE E_GeV in Standard Rock (Groom 2001)."""
    T = float(np.clip(E_GeV * 1000.0, _GROOM_T_MEV[0], _GROOM_T_MEV[-1]))
    return float(np.exp(
        np.interp(np.log(T), np.log(_GROOM_T_MEV), np.log(_GROOM_R_GCM2))
    ))


def _sr(E: float, th: float, d_cm: float, rho=2.65, X0=26.7) -> float:
    """Highland RMS radial MCS displacement [cm]."""
    M = 0.10566; p = np.sqrt(max(E**2 - M**2, 1e-6))
    cz = max(np.cos(np.radians(th)), 0.02); L = d_cm / cz
    xr = max(rho * L / X0, 1e-12); beta = p / np.sqrt(p**2 + M**2)
    th0 = (13.6 / (beta * p * 1000)) * np.sqrt(xr) * (1 + 0.038 * np.log(xr))
    return float(np.sqrt(2) * (L / np.sqrt(3)) * th0)


def _mcs_emin(th: float, d_cm: float, R_cm: float, rho=2.65, X0=26.7) -> float:
    """Min E [GeV] such that σ_r(E, θ, d) ≤ R. Returns inf if impossible."""
    from scipy.optimize import brentq
    if _sr(1.0,    th, d_cm, rho, X0) <= R_cm: return 1.0
    if _sr(5000.0, th, d_cm, rho, X0) >  R_cm: return float("inf")
    return brentq(lambda E: _sr(E, th, d_cm, rho, X0) - R_cm, 1.0, 5000.0, xtol=0.5)


def _fmt_E(E: float) -> str:
    """Energy with sensible precision: 0.05, 2.3, 47."""
    if not np.isfinite(E): return "∞"
    if E >= 10: return f"{E:.0f}"
    if E >= 1:  return f"{E:.1f}"
    return f"{E:.2f}"


def _det_geom(det: dict) -> dict:
    """Unified geometry from one detector dict."""
    if det["shape"] == 1:
        cx = (float(det["ax"]) + float(det["bx"])) / 2.0
        cy = (float(det["ay"]) + float(det["by"])) / 2.0
        za = float(det["az"]); zb = float(det["bz"])
        d_top = abs(max(za, zb)); d_bot = abs(min(za, zb))
        R_eff = float(det["r"]) + float(det.get("margin", 0.0))
    else:
        cx = (float(det["xmin"]) + float(det["xmax"])) / 2.0
        cy = (float(det["ymin"]) + float(det["ymax"])) / 2.0
        zlo = float(det["zmin"]); zhi = float(det["zmax"])
        d_top = abs(max(zlo, zhi)); d_bot = abs(min(zlo, zhi))
        R_eff = (0.5 * max(abs(det["xmax"] - det["xmin"]), abs(det["ymax"] - det["ymin"]))
                 + float(det.get("margin", 0.0)))
    return dict(cx=cx, cy=cy,
                depth_top=d_top, depth_bot=d_bot,
                depth_mid=(d_top + d_bot) / 2.0,
                R_eff=R_eff, half_m=R_eff / 100.0,
                offset_m=np.sqrt(cx**2 + cy**2) / 100.0,
                margin=float(det.get("margin", 0.0)))


def _compute(detectors, e_min, theta_max_deg, rho, X0):
    geoms   = [_det_geom(d) for d in detectors]
    src_cx  = float(np.mean([g["cx"] for g in geoms]))
    src_cy  = float(np.mean([g["cy"] for g in geoms]))

    # CSDA energy floors
    depth_top_all    = max(g["depth_top"] for g in geoms)
    depth_bot_all    = max(g["depth_bot"] for g in geoms)
    E_floor_enter    = _groom(rho * depth_top_all) * 1.1
    E_floor_traverse = _groom(rho * depth_bot_all) * 1.1

    # Effective energy per detector: a muon contributing at depth_mid must
    # survive to depth_mid, so scattering is never evaluated below the CSDA
    # threshold even when the user sets a lower E_min.
    E_effs = [max(float(e_min), _groom(rho * g["depth_mid"])) for g in geoms]

    # MCS in-scattering margins (1.5 × σ_r at E_eff)
    mcs_margins = []
    for g, E_eff in zip(geoms, E_effs):
        dx = g["cx"] - src_cx; dy = g["cy"] - src_cy
        th = float(np.degrees(np.arctan2(np.sqrt(dx**2 + dy**2), g["depth_mid"])))
        mcs_margins.append(1.5 * _sr(E_eff, th, g["depth_mid"], rho, X0))

    # MCS E_min (for directional detectors)
    E_mcs_list = []
    for g in geoms:
        dx = g["cx"] - src_cx; dy = g["cy"] - src_cy
        th = float(np.degrees(np.arctan2(np.sqrt(dx**2 + dy**2), g["depth_mid"])))
        E_mcs_list.append(_mcs_emin(th, g["depth_mid"], g["R_eff"], rho, X0))
    E_mcs = float(max(E_mcs_list))

    # Muography energy window
    E_thr_shallowest     = _groom(rho * depth_top_all)
    E_thr_deepest_vert   = _groom(rho * depth_bot_all)
    cap_deg              = min(float(theta_max_deg), 60.0)
    E_thr_deepest_oblique = _groom(rho * depth_bot_all / max(np.cos(np.radians(cap_deg)), 0.02))
    E_max_muogr          = 5.0 * E_thr_deepest_oblique
    muogr_window_ok      = E_thr_shallowest < E_max_muogr

    # Optimised source (centred at detector centroid, no depth×tan term)
    src_radii = []
    for i, g in enumerate(geoms):
        dist = float(np.sqrt((g["cx"] - src_cx)**2 + (g["cy"] - src_cy)**2))
        src_radii.append(dist + g["R_eff"] + mcs_margins[i])
    src_r_cm = float(max(src_radii))

    if len(geoms) == 1:
        g = geoms[0]
        max_x = abs(g["cx"] - src_cx) + g["R_eff"] + mcs_margins[0]
        max_y = abs(g["cy"] - src_cy) + g["R_eff"] + mcs_margins[0]
    else:
        offsets_x = [abs(g["cx"] - src_cx) + g["R_eff"] + mcs_margins[i]
                     for i, g in enumerate(geoms)]
        offsets_y = [abs(g["cy"] - src_cy) + g["R_eff"] + mcs_margins[i]
                     for i, g in enumerate(geoms)]
        max_x = float(max(offsets_x)); max_y = float(max(offsets_y))

    # Geometric source radius from origin at current θ_max
    tan_th   = float(np.tan(np.radians(theta_max_deg)))
    r_geo_rows = []
    for i, g in enumerate(geoms):
        r_geo = g["offset_m"] + g["half_m"] + g["depth_bot"] / 100.0 * tan_th
        r_geo_rows.append(dict(det=i + 1, depth_m=g["depth_bot"] / 100.0,
                               offset_m=g["offset_m"], half_m=g["half_m"],
                               r_min_m=r_geo))
    r_geo_max = max(r["r_min_m"] for r in r_geo_rows)

    # θ recommendations
    theta_geom_list = []
    for i, g in enumerate(geoms):
        if g["depth_top"] < 100.0:
            # Entry face essentially at the source plane: the source area
            # covers it directly, so it imposes no zenith requirement.
            theta_geom_list.append(0.0)
            continue
        dist = float(np.sqrt((g["cx"] - src_cx)**2 + (g["cy"] - src_cy)**2))
        th_geom = float(np.degrees(
            np.arctan2(dist + g["R_eff"] + mcs_margins[i], g["depth_top"])))
        theta_geom_list.append(min(th_geom, 88.0))
    theta_geom_rec = float(max(theta_geom_list))

    theta_mcs_list = []
    for g, E_eff in zip(geoms, E_effs):
        f_lo = _sr(E_eff, 0.1,  g["depth_top"], rho, X0) - g["R_eff"]
        f_hi = _sr(E_eff, 88.0, g["depth_top"], rho, X0) - g["R_eff"]
        if f_lo <= 0:
            th_mcs = 88.0
        elif f_hi > 0:
            th_mcs = 0.1
        else:
            from scipy.optimize import brentq
            th_mcs = brentq(lambda th: _sr(E_eff, th, g["depth_top"], rho, X0) - g["R_eff"],
                            0.1, 88.0, xtol=0.1)
        theta_mcs_list.append(float(th_mcs))
    theta_mcs_rec = float(min(theta_mcs_list))
    theta_rec     = float(np.clip(theta_mcs_rec, theta_geom_rec, 88.0))
    theta_conflict      = theta_geom_rec > theta_mcs_rec + 2.0
    theta_unconstrained = theta_mcs_rec >= 87.9 and theta_geom_rec <= 0.0

    return dict(
        geoms=geoms,
        src_cx_m=src_cx / 100.0, src_cy_m=src_cy / 100.0,
        src_r_m=src_r_cm / 100.0,
        src_lx_m=max_x / 100.0, src_ly_m=max_y / 100.0,
        mcs_margins=mcs_margins,
        E_eff=float(max(E_effs)), per_det_E_eff=E_effs,
        depth_top_m=depth_top_all / 100.0, depth_bot_m=depth_bot_all / 100.0,
        E_floor_enter=E_floor_enter, E_floor_traverse=E_floor_traverse,
        E_mcs=E_mcs,
        E_thr_shallowest=E_thr_shallowest,
        E_thr_deepest_vert=E_thr_deepest_vert,
        E_max_muogr=E_max_muogr,
        muogr_window_ok=muogr_window_ok,
        muogr_window_E_lo=E_thr_shallowest,
        muogr_window_E_hi=E_max_muogr,
        theta_geom_rec=theta_geom_rec,
        theta_mcs_rec=theta_mcs_rec,
        theta_rec=theta_rec,
        theta_conflict=theta_conflict,
        theta_unconstrained=theta_unconstrained,
        per_det_theta_mcs=theta_mcs_list,
        per_det_theta_geom=theta_geom_list,
        per_det_mcs_margin=mcs_margins,
        per_det_r=src_radii,
        r_geo_max=r_geo_max,
        r_geo_rows=r_geo_rows,
    )


def render_combined_source_panel(
        detectors, e_min, e_max, theta_max, source_mode, radius, plane_lx, plane_ly
) -> None:
    """Source size, energy window and θ_max recommendations.

    Renders directly into the current container (no outer expander) so it can
    be embedded in the 🧰 Helpers & calculators section.
    """
    if not detectors:
        st.info("Define at least one detector to get recommendations.")
        return

    # ── Material ──────────────────────────────────────────────────────────
    _mc1, _mc2, _mc3 = st.columns([2, 1, 1])
    mat = _mc1.selectbox("Overburden material",
                         ["Standard Rock (ρ=2.65, X₀=26.7)", "Custom"],
                         key="srcopt_mat")
    if mat == "Custom":
        rho = _mc2.number_input("ρ [g/cm³]", 0.1, 20.0, 2.65,
                                step=0.01, format="%.3f", key="srcopt_rho")
        X0  = _mc3.number_input("X₀ [g/cm²]", 1.0, 100.0, 26.7,
                                step=0.1,  format="%.1f",  key="srcopt_X0")
    else:
        rho, X0 = 2.65, 26.7
        _mc2.metric("ρ", "2.650 g/cm³")
        _mc3.metric("X₀", "26.7 g/cm²")

    try:
        opt = _compute(detectors, e_min, theta_max, rho, X0)
    except Exception as exc:
        st.error(f"❌  Computation failed: {exc}"); return

    geoms  = opt["geoms"]
    n_det  = len(detectors)
    _cur_r = radius if source_mode in (1, 3) else max(plane_lx, plane_ly)

    # ══════════════════════════════════════════════════════════════════════
    # A — Source size
    # ══════════════════════════════════════════════════════════════════════
    st.markdown("#### 📐 Source size")

    if source_mode == 3:
        st.info("Hemisphere mode: source centring not applicable. "
                "Energy and θ_max recommendations below still apply.")
    else:
        # Metrics: optimised (centred at detector centroid) + geometric at current θ_max
        _a1, _a2, _a3, _a4, _a5 = st.columns(5)
        _a1.metric("Centre X", f"{opt['src_cx_m']:.2f} m")
        _a2.metric("Centre Y", f"{opt['src_cy_m']:.2f} m")
        if source_mode == 1:
            _opt_r = opt["src_r_m"]
            _a3.metric("Optimised radius",   f"{_opt_r:.2f} m",
                       help="Centred at detector centroid + 1.5 × MCS margin "
                            f"(σ_r at E = {_fmt_E(opt['E_eff'])} GeV — never below "
                            "the CSDA survival threshold)")
            _a4.metric("Geometric (θ_max)",  f"{opt['r_geo_max']:.0f} m",
                       help="Min radius from origin to cover all detectors at current θ_max")
        else:
            _opt_r = max(opt["src_lx_m"], opt["src_ly_m"])
            _a3.metric("Half-width Lx", f"{opt['src_lx_m']:.2f} m")
            _a4.metric("Half-width Ly", f"{opt['src_ly_m']:.2f} m")
        _area_gain = (_cur_r / max(_opt_r, 0.01)) ** 2
        _a5.metric("Area gain vs current",
                   f"×{_area_gain:.0f}" if _area_gain >= 2 else "—",
                   help="How many fewer wasted muons with the optimised source")

        st.caption(
            f"MCS margin evaluated at E = {_fmt_E(opt['E_eff'])} GeV "
            "(max of E_min and the CSDA threshold at detector depth).  "
            "ℹ️ The optimised source omits the depth × tan θ reach, so oblique "
            "trajectories are under-sampled — use it with the detector filter for "
            "hit-count studies; for unbiased angular spectra use the geometric radius."
        )

        # Alerts only when action is needed
        _shift = np.sqrt(opt["src_cx_m"]**2 + opt["src_cy_m"]**2)
        if _shift > 0.5 or (_cur_r - _opt_r > 0.5 * _opt_r):
            _parts = []
            if _shift > 0.5:
                _parts.append(f"shift centre to ({opt['src_cx_m']:.1f}, "
                              f"{opt['src_cy_m']:.1f}) m")
            if _cur_r - _opt_r > 0.5 * _opt_r:
                _parts.append(f"reduce size to **{_opt_r:.1f} m** "
                              f"(×{_area_gain:.0f} more signal hits per CPU hour)")
            st.info("💡 " + "  |  ".join(_parts) + ".", icon="💡")

        # Multi-detector breakdown (collapsed)
        if n_det > 1:
            with st.expander("Per-detector breakdown", expanded=False):
                st.dataframe(pd.DataFrame([{
                    "Det #":             i + 1,
                    "Axis (X,Y) [m]":    f"({g['cx']/100:.1f}, {g['cy']/100:.1f})",
                    "Depth top [m]":     f"{g['depth_top']/100:.1f}",
                    "R_eff [cm]":        f"{g['R_eff']:.1f}",
                    "E_eff [GeV]":       _fmt_E(opt['per_det_E_eff'][i]),
                    "MCS margin [cm]":   f"{opt['mcs_margins'][i]:.1f}",
                    "Min radius [m]":    f"{opt['per_det_r'][i]/100:.2f}",
                } for i, g in enumerate(geoms)]).set_index("Det #"),
                width='stretch')

    # ══════════════════════════════════════════════════════════════════════
    # B — Energy window
    # ══════════════════════════════════════════════════════════════════════
    st.divider()
    st.markdown("#### ⚡ Energy window")

    _e1, _e2, _e3, _e4 = st.columns(4)
    if opt["depth_top_m"] < 1.0:
        _e1.metric("E_min (enter top face)", "—",
                   help="Detector top face is at the surface — any muon enters it.")
    else:
        _e1.metric("E_min (enter top face)", f"{_fmt_E(opt['E_floor_enter'])} GeV",
                   help="1.1 × CSDA threshold to reach the shallowest detector face")
    _e2.metric("E_min (traverse)",        f"{_fmt_E(opt['E_floor_traverse'])} GeV",
               help="1.1 × CSDA threshold to reach the deepest detector face")
    if opt["muogr_window_ok"]:
        _e3.metric("Muography E window",
                   f"{_fmt_E(opt['muogr_window_E_lo'])} – {_fmt_E(opt['muogr_window_E_hi'])} GeV",
                   help="Range where transmission varies with rock density")
    else:
        _e3.metric("Muography E window", "—  (too narrow)")
    _mcs_str = f"{opt['E_mcs']:.0f} GeV" if not np.isinf(opt["E_mcs"]) else "flux counter only"
    _e4.metric("E_min (MCS pointing)",   _mcs_str,
               help="Min energy for σ_r < R_det (directional detector)")

    # Alerts
    if opt["depth_top_m"] >= 1.0 and e_min < opt["E_floor_enter"] * 0.85:
        st.error(f"🚨  E_min = {e_min:.0f} GeV is below the CSDA floor "
                 f"({_fmt_E(opt['E_floor_enter'])} GeV) — most generated muons stop "
                 f"before reaching the detector.  Set E_min ≥ **{_fmt_E(opt['E_floor_enter'])} GeV**.")
    if opt["muogr_window_ok"] and e_max > opt["muogr_window_E_hi"] * 1.2:
        st.warning(f"⚠️  E_max = {e_max:.0f} GeV exceeds the muography window "
                   f"({_fmt_E(opt['muogr_window_E_hi'])} GeV). "
                   f"High-energy muons add flux counts but no transmission contrast.")
    if (not np.isinf(opt["E_mcs"]) and opt["E_mcs"] > opt["E_floor_enter"] * 1.5
            and not (opt["muogr_window_ok"] and
                     opt["E_mcs"] < opt["muogr_window_E_hi"])):
        st.warning(f"⚠️  MCS pointing requires E > {opt['E_mcs']:.0f} GeV, "
                   f"but muography window ends at {_fmt_E(opt['muogr_window_E_hi'])} GeV. "
                   f"This detector works as a flux counter, not a directional tracker.")

    # Transmission table (collapsed)
    with st.expander("Transmission table T(E_min, depth)", expanded=False):
        _depths = sorted({round(g["depth_top"] / 100, 0) for g in geoms}
                       | {round(g["depth_bot"] / 100, 0) for g in geoms})[:5]
        _E_rows = sorted({round(float(e_min), 0),
                          round(float(opt["E_floor_enter"]), 0),
                          round(float(opt["muogr_window_E_lo"]), 0),
                          round(float(opt["E_thr_deepest_vert"]), 0),
                          round(float(opt["muogr_window_E_hi"]), 0)})
        _hdrs = ["E_min [GeV]"] + [f"D={d:.0f}m" for d in _depths] + ["Contrast"]
        _rows_T = []
        for _E in _E_rows:
            _row = [f"{_E:.0f}"]
            _Tvals = []
            for _D in _depths:
                _Ethr = _groom(rho * _D * 100)
                _T    = 1.0 if _E >= _Ethr else (_E / _Ethr) ** (_GAMMA - 1)
                _Tvals.append(_T)
                _row.append(f"{_T * 100:.1f}%")
            _row.append(f"{(max(_Tvals) - min(_Tvals)) * 100:.1f}%")
            _rows_T.append(_row)
        st.dataframe(pd.DataFrame(_rows_T, columns=_hdrs).set_index("E_min [GeV]"),
                     width='stretch')
        st.caption("T ≈ (E_min / E_thr)^(γ−1), γ=3.7. Contrast = max(T) − min(T) across depths.")

    # ══════════════════════════════════════════════════════════════════════
    # C — θ_max
    # ══════════════════════════════════════════════════════════════════════
    st.divider()
    st.markdown("#### 🧭 Recommended θ_max")

    if opt["theta_unconstrained"]:
        st.info(f"No MCS or geometric constraint on θ_max for this geometry — "
                f"choose θ_max from flux considerations. Current {theta_max:.0f}° is fine; "
                f"larger values only add near-horizontal flux.")
    else:
        _t1, _t2, _t3 = st.columns(3)
        _t1.metric("MCS limit",         f"{opt['theta_mcs_rec']:.1f}°",
                   help=f"Max θ where σ_r(E = {_fmt_E(opt['E_eff'])} GeV, θ, depth_top) ≤ R_det")
        _t2.metric("Geometric minimum", f"{opt['theta_geom_rec']:.1f}°",
                   help="Min θ to illuminate all detector edges from the optimised centre "
                        "(0° = entry face covered by the source area directly)")
        _t3.metric("Recommended θ_max", f"{opt['theta_rec']:.1f}°",
                   delta=f"{opt['theta_rec'] - theta_max:+.1f}° vs current",
                   delta_color="normal" if abs(opt["theta_rec"] - theta_max) < 5 else "inverse")

        if opt["theta_conflict"]:
            st.warning(f"⚠️  Geometric coverage needs θ_max ≥ {opt['theta_geom_rec']:.1f}° but "
                       f"the MCS limit is {opt['theta_mcs_rec']:.1f}° — edge coverage and "
                       f"directional pointing cannot both be satisfied. The recommendation "
                       f"favours coverage; treat the detector as a flux counter at large θ.")
        elif theta_max > opt["theta_mcs_rec"] + 5:
            st.warning(f"⚠️  θ_max = {theta_max:.0f}° exceeds MCS limit "
                       f"({opt['theta_mcs_rec']:.1f}°). Reduce to **{opt['theta_rec']:.1f}°**.")
        elif theta_max < opt["theta_geom_rec"] - 2:
            st.warning(f"⚠️  θ_max = {theta_max:.0f}° may not cover all detector edges "
                       f"from the optimised source centre (need ≥ {opt['theta_geom_rec']:.1f}°).")

        if n_det > 1:
            with st.expander("Per-detector θ breakdown", expanded=False):
                st.dataframe(pd.DataFrame([{
                    "Det #":        i + 1,
                    "Depth top [m]": f"{g['depth_top'] / 100:.1f}",
                    "R_eff [cm]":   f"{g['R_eff']:.1f}",
                    "θ_mcs [°]":    f"{opt['per_det_theta_mcs'][i]:.1f}",
                    "θ_geom [°]":   f"{opt['per_det_theta_geom'][i]:.1f}",
                } for i, g in enumerate(geoms)]).set_index("Det #"),
                width='stretch')

    # ══════════════════════════════════════════════════════════════════════
    # D — Summary
    # ══════════════════════════════════════════════════════════════════════
    st.divider()
    st.markdown("#### 📋 Settings at a glance")

    if source_mode == 1:
        _src_line = (f"Source disc    centre ({opt['src_cx_m']:.2f}, "
                     f"{opt['src_cy_m']:.2f}) m   radius {opt['src_r_m']:.2f} m")
    elif source_mode == 2:
        _src_line = (f"Source rect    centre ({opt['src_cx_m']:.2f}, "
                     f"{opt['src_cy_m']:.2f}) m   "
                     f"Lx {opt['src_lx_m']:.2f} m  Ly {opt['src_ly_m']:.2f} m")
    else:
        _src_line = "Source         hemisphere (centring not applicable)"

    _e_line = f"E_min          ≥ {_fmt_E(opt['E_floor_enter'])} GeV  (CSDA survival)"
    if opt["muogr_window_ok"]:
        _e_line += (f"\nE muography    {_fmt_E(opt['muogr_window_E_lo'])} – "
                    f"{_fmt_E(opt['muogr_window_E_hi'])} GeV")
    if not np.isinf(opt["E_mcs"]) and opt["E_mcs"] > opt["E_floor_enter"] * 1.5:
        _e_line += f"\nE directional  ≥ {opt['E_mcs']:.0f} GeV  (MCS < R_det)"

    if opt["theta_unconstrained"]:
        _th_line = f"θ_max          {theta_max:.1f}°  (unconstrained — user choice)"
    else:
        _th_line = f"θ_max          {opt['theta_rec']:.1f}°"

    st.code(f"{_src_line}\n{_e_line}\n{_th_line}", language="text")
