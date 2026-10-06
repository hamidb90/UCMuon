"""
basic_mode.py — the values the GUI's Basic mode derives instead of asking for.

Pure functions, no Streamlit, so they can be unit-tested
(tests/gui/test_basic_mode.py). What Basic derives is listed in CHANGELOG.md
(1.3.0, "Basic / Advanced mode").

Units follow the GUI: energies are total energies in GeV, detector geometry
is in cm (the per-detector dicts the generator receives), source sizes in m.
Depths are positive metres below the ground (z = 0).
"""
import math

import mcs_margin as _mm

# ── Recommended energy (and zenith) range per spectrum ──────────────────────
# From the README validity table ("Surface spectra"); Hamid's decision
# 2026-10-01: choosing a spectrum in Basic sets its range automatically.
# emin/emax None: the user's own range (power law, a sampling shape only).
RECOMMENDED = {
    1: dict(name="CosmoALEPH", emin=100.0, emax=2500.0, theta_max=85.0, normalised=True,
            note="fit to vertical data at 112-2239 GeV/c"),
    2: dict(name="Power law E⁻³·⁷", emin=None, emax=None, theta_max=85.0, normalised=False,
            note="no absolute normalisation: no rate or live time"),
    3: dict(name="PARMA/EXPACS", emin=1.0, emax=1000.0, theta_max=85.0, normalised=True,
            note="site-aware (altitude, latitude, solar activity of the date)"),
    4: dict(name="Guan 2015", emin=1.0, emax=2500.0, theta_max=85.0, normalised=True,
            note="fitted above 1 GeV, all zenith angles"),
    5: dict(name="Frosin 2025", emin=1.0, emax=2500.0, theta_max=85.0, normalised=True,
            note="fitted above 1 GeV, all zenith angles"),
    6: dict(name="Gaisser 1990", emin=100.0, emax=2500.0, theta_max=60.0, normalised=True,
            note="valid only above 100/cos θ GeV"),
    7: dict(name="Reyna 2006", emin=1.0, emax=2000.0, theta_max=85.0, normalised=True,
            note="valid for 1 < p < 2000/cos θ GeV/c"),
    8: dict(name="Cosmic e±", emin=0.01, emax=1.0, theta_max=85.0, normalised=False,
            note="no absolute normalisation: no rate or live time"),
}

BASIC_ANGULAR_MODE = 6        # joint J(p, θ) with the surface projection (the default)
BASIC_SOURCE_SIZE_M = 200.0   # source radius / half-width without a detector
DEFAULT_SPECTRUM = 4          # Guan 2015 (decision 2026-10-05, both modes)
DEFAULT_EMIN, DEFAULT_EMAX = 1.0, 2500.0


def recommended(spectrum):
    """The Basic range for a spectrum (a copy of its RECOMMENDED entry)."""
    return dict(RECOMMENDED.get(int(spectrum), RECOMMENDED[DEFAULT_SPECTRUM]))


# ── Detector geometry ────────────────────────────────────────────────────────
def top_depth_m(det):
    """Depth of the detector's top face below z = 0, in m (>= 0 above ground)."""
    z = max(det["az"], det["bz"]) if det["shape"] == 1 else max(det["zmin"], det["zmax"])
    return -float(z) / 100.0


def bottom_depth_m(det):
    """Depth of the detector's deepest point below z = 0, in m."""
    z = min(det["az"], det["bz"]) if det["shape"] == 1 else min(det["zmin"], det["zmax"])
    return -float(z) / 100.0


def detector_from_basic(shape, cx_m, cy_m, top_m, height_m, radius_m=1.0,
                        width_m=1.0, length_m=1.0, margin_cm=0.0):
    """A generator detector dict (cm) from the Basic inputs: a vertical
    cylinder (shape 1) or a box (shape 2) centred on (cx, cy), whose top face
    is top_m below the ground and which is height_m tall."""
    z_top, z_bot = -100.0 * float(top_m), -100.0 * (float(top_m) + float(height_m))
    if int(shape) == 1:
        return {"shape": 1, "margin": float(margin_cm),
                "ax": 100.0 * cx_m, "ay": 100.0 * cy_m, "az": z_bot,
                "bx": 100.0 * cx_m, "by": 100.0 * cy_m, "bz": z_top,
                "r": 100.0 * float(radius_m)}
    return {"shape": 2, "margin": float(margin_cm),
            "xmin": 100.0 * (cx_m - width_m / 2.0), "xmax": 100.0 * (cx_m + width_m / 2.0),
            "ymin": 100.0 * (cy_m - length_m / 2.0), "ymax": 100.0 * (cy_m + length_m / 2.0),
            "zmin": z_bot, "zmax": z_top}


def basic_from_detector(det):
    """The Basic inputs for a detector dict, or None when Basic cannot show it
    (a cylinder whose axis is not vertical)."""
    if det["shape"] == 1:
        if abs(det["ax"] - det["bx"]) > 1e-6 or abs(det["ay"] - det["by"]) > 1e-6:
            return None
        return dict(shape=1, cx_m=det["ax"] / 100.0, cy_m=det["ay"] / 100.0,
                    top_m=top_depth_m(det), height_m=abs(det["az"] - det["bz"]) / 100.0,
                    radius_m=det["r"] / 100.0)
    return dict(shape=2, cx_m=(det["xmin"] + det["xmax"]) / 200.0,
                cy_m=(det["ymin"] + det["ymax"]) / 200.0, top_m=top_depth_m(det),
                height_m=abs(det["zmax"] - det["zmin"]) / 100.0,
                width_m=abs(det["xmax"] - det["xmin"]) / 100.0,
                length_m=abs(det["ymax"] - det["ymin"]) / 100.0)


def detector_centre_m(det):
    """Horizontal centre (x, y) of a detector, in m."""
    if det["shape"] == 1:
        return (det["ax"] + det["bx"]) / 200.0, (det["ay"] + det["by"]) / 200.0
    return (det["xmin"] + det["xmax"]) / 200.0, (det["ymin"] + det["ymax"]) / 200.0


# ── Safety margin ────────────────────────────────────────────────────────────
def margin_cm(det, rho, e_min_total_gev):
    """2σ_r of multiple scattering at the detector's deepest point
    (gui/mcs_margin.py, the rule the Generator's margin warning checks).
    0 for a detector that does not reach below the ground."""
    depth = bottom_depth_m(det)
    if depth <= 0.0:
        return 0.0
    rec = _mm.suggested_margin(depth, float(rho), float(e_min_total_gev), 0.0)
    return float(rec["margin_cm"] or 0.0)


# ── Source size ──────────────────────────────────────────────────────────────
def source_reach(det, theta_max_deg, src_z_m, source_mode, source_plane, disk_tilt,
                 disk_cx_m, disk_cy_m, disk_r_m, u1_m, u2_m, v1_m, v2_m):
    """For a horizontal disk or rectangle source (not tilted), whether every
    straight path at zenith <= theta_max into the margin-inflated detector
    starts on the source. Returns None when it does (or the check does not
    apply), else (needed, have, unit_text): the radius, or the half-widths,
    the source needs. A muon reaching a point at depth d below the source at
    zenith θ starts d·tan θ away horizontally, so the deepest points decide."""
    if source_plane != 1 or source_mode not in (1, 2) or abs(disk_tilt) > 1e-9:
        return None
    t = math.tan(math.radians(min(float(theta_max_deg), 89.9)))
    z_src, m = float(src_z_m) * 100.0, float(det["margin"])
    if det["shape"] == 1:
        rr = float(det["r"]) + m
        pts = [(det["ax"], det["ay"], min(det["az"], det["bz"]) - m),
               (det["bx"], det["by"], min(det["az"], det["bz"]) - m),
               (det["ax"], det["ay"], det["az"]), (det["bx"], det["by"], det["bz"])]
        pts = [(x, y, z, rr) for x, y, z in pts]
    else:
        xs = (det["xmin"] - m, det["xmax"] + m)
        ys = (det["ymin"] - m, det["ymax"] + m)
        z0 = min(det["zmin"], det["zmax"]) - m
        pts = [(x, y, z0, 0.0) for x in xs for y in ys]
    pts = [(x, y, z, r) for x, y, z, r in pts if z < z_src]
    if not pts:
        return None
    if source_mode == 1:
        need = max(math.hypot(x - disk_cx_m * 100.0, y - disk_cy_m * 100.0) + r
                   + (z_src - z) * t for x, y, z, r in pts) / 100.0
        return None if need <= disk_r_m * 1.0001 else (need, disk_r_m, "radius")
    lo_u = min(x - r - (z_src - z) * t for x, y, z, r in pts) / 100.0
    hi_u = max(x + r + (z_src - z) * t for x, y, z, r in pts) / 100.0
    lo_v = min(y - r - (z_src - z) * t for x, y, z, r in pts) / 100.0
    hi_v = max(y + r + (z_src - z) * t for x, y, z, r in pts) / 100.0
    if lo_u >= u1_m - 1e-6 and hi_u <= u2_m + 1e-6 and lo_v >= v1_m - 1e-6 and hi_v <= v2_m + 1e-6:
        return None
    return ((lo_u, hi_u, lo_v, hi_v), (u1_m, u2_m, v1_m, v2_m), "extent")


def _round_up(x, step=10.0):
    return step * math.ceil(float(x) / step - 1e-9)


def required_source(dets, theta_max_deg, source_mode):
    """The Basic source for the (margin-inflated) detectors: horizontal, at
    z = 0, centred on the detectors, just large enough that every straight
    path at zenith <= theta_max into them starts on it (source_reach), rounded
    up to 10 m. Returns a dict with the generator's source parameters.

    Hemisphere (source_mode 3): the generator's hemisphere is centred on the
    z axis, so its radius is the disk radius measured from (0, 0) at z = 0: a
    path that crosses z = 0 inside it then leaves the sphere through its upper
    half."""
    if int(source_mode) == 3:
        cx = cy = 0.0
    else:
        cxs, cys = zip(*(detector_centre_m(d) for d in dets))
        cx, cy = sum(cxs) / len(cxs), sum(cys) / len(cys)
    out = dict(source_mode=int(source_mode), source_plane=1, src_w_m=0.0,
               disk_tilt=0.0, disk_tilt_az=0.0, disk_cx=cx, disk_cy=cy)
    if int(source_mode) in (1, 3):
        need = 0.0
        for d in dets:
            res = source_reach(d, theta_max_deg, 0.0, 1, 1, 0.0, cx, cy, 0.0,
                               0.0, 0.0, 0.0, 0.0)
            if res is not None:
                need = max(need, res[0])
        r = max(_round_up(need), 10.0)
        out.update(disk_r=r, src_u1_m=cx - r, src_u2_m=cx + r,
                   src_v1_m=cy - r, src_v2_m=cy + r)
        return out
    lo_u = lo_v = math.inf
    hi_u = hi_v = -math.inf
    for d in dets:
        res = source_reach(d, theta_max_deg, 0.0, 2, 1, 0.0, cx, cy, 0.0,
                           cx, cx, cy, cy)
        if res is None:
            continue
        (a, b, c, e), _, _ = res
        lo_u, hi_u, lo_v, hi_v = min(lo_u, a), max(hi_u, b), min(lo_v, c), max(hi_v, e)
    if not math.isfinite(lo_u):                    # nothing below ground
        lo_u, hi_u, lo_v, hi_v = cx - 10.0, cx + 10.0, cy - 10.0, cy + 10.0
    u1, u2 = cx - _round_up(cx - lo_u), cx + _round_up(hi_u - cx)
    v1, v2 = cy - _round_up(cy - lo_v), cy + _round_up(hi_v - cy)
    out.update(src_u1_m=u1, src_u2_m=u2, src_v1_m=v1, src_v2_m=v2,
               disk_cx=(u1 + u2) / 2.0, disk_cy=(v1 + v2) / 2.0,
               disk_r=min(u2 - u1, v2 - v1) / 2.0)
    return out


def free_source(source_mode, size_m):
    """The Basic source without a detector: horizontal, at z = 0, centred on
    the origin; size_m is the radius (disk, hemisphere) or the half-width of
    a square (rectangle). Same keys as required_source."""
    r = float(size_m)
    return dict(source_mode=int(source_mode), source_plane=1, src_w_m=0.0,
                disk_tilt=0.0, disk_tilt_az=0.0, disk_cx=0.0, disk_cy=0.0,
                disk_r=r, src_u1_m=-r, src_u2_m=r, src_v1_m=-r, src_v2_m=r)


# ── Controls Basic hides that are flagged for repair or removal ─────────────
# The 39 "Remove?" controls of the inventory (dead, duplicated or misleading).
# Decision 2026-10-01: hidden in Basic, not deleted; repaired or removed later
# one by one. (key or label, tab, reason)
REMOVE_CANDIDATES = [
    *[(k, "Generator", "manual PHITS export duplicating the automatic one")
      for k in ("phits_surf_enable", "phits_surf_src", "phits_surf_out",
                "btn_phits_surf", "dl_phits_surf")],
    *[(k, "Generator", "MCS-estimator material, ignores the Transport density")
      for k in ("mcs_material", "mcs_rho_custom", "mcs_X0_custom")],
    *[(k, "Generator", "source-optimizer material, ignores the Transport density")
      for k in ("srcopt_mat", "srcopt_rho", "srcopt_X0")],
    *[(k, "Transport", "density spread reaches no engine (B9c, now stated)")
      for k in ("density_mode", "music_rho (Gaussian mean)", "music_rho_sigma",
                "density_kde_file")],
    ("proposal_med_choice", "Transport", "duplicates the material preset"),
    ("phitsxs_mat_choice", "Transport", "duplicates the material preset"),
    ("music_omp_threads", "Transport", "only MUSIC uses it (B14a, now MUSIC only)"),
    ("terrain surface-file selector", "Terrain", "was ignored by the run (B10a, fixed)"),
    *[(k, "Terrain", "duplicate CSG single-ray check reading unset keys")
      for k in ("csg_test_az", "csg_test_ze", "csg_test_ray")],
    *[(k, "Results", "3D overlay that rewrites the run's detector metadata")
      for k in ("ov_sh", "ov_ax", "ov_ay", "ov_az", "ov_bx", "ov_by", "ov_bz", "ov_r",
                "ov_xn", "ov_yn", "ov_zn", "ov_xx", "ov_yx", "ov_zx", "ov_apply")],
    ("Download config JSON", "Config", "restored nothing up to v1.2.0 (B12, fixed)"),
    ("Reset autosave", "Config", "did not reset up to v1.2.0 (B8, fixed)"),
]
