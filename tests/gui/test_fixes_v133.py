"""Regression tests for the fixes after 1.3.2 (each fails on 1.3.2)."""
import re
from pathlib import Path

from conftest import exceptions, widget

ROOT = Path(__file__).resolve().parents[2]


def test_autosave_that_is_not_an_object_does_not_crash(new_app, repo_copy):
    """Valid JSON that is not an object (here a list) crashed the app on start
    with AttributeError on data.items(); it is now ignored."""
    (repo_copy / "ucmuon_autosave.json").write_text("[1, 2, 3]")
    at = new_app(keep_autosave=True)
    assert not exceptions(at)
    assert widget(at, "toggle", "ui_mode_toggle") is not None


def test_python_bethe_bloch_ice_radiation_length_matches_fortran():
    """The Python Bethe-Bloch driver used X0 = 33.10 g/cm2 for ice; the Fortran
    engine, UCMuon-MC and PDG give 36.08 g/cm2 (as water)."""
    import ucmuon_bb_driver as bb
    f90 = (ROOT / "src/transport/bethe_bloch/ucmuon_transport_bb_omp.f90").read_text()
    ice = re.search(r"case\(2\)\s*!\s*Ice.*?X0_gcm2=([0-9.]+)d0", f90, re.S)
    assert ice, "ice block not found in the Fortran Bethe-Bloch engine"
    assert bb._BB_MAT[2]["name"] == "Ice"
    assert bb._BB_MAT[2]["X0_gcm2"] == float(ice.group(1)) == 36.08


def test_phits_dump_with_fortran_d_exponents_loads(tmp_path):
    """load_file could not read a PHITS dump written with Fortran D exponents
    (the generator's own, the GUI's PHITS export, ucmuon_to_phits.f90): pandas
    left every column as text and the load failed. Electrons (kf 11) also got
    the muon charge test and the muon mass."""
    import ast
    import types

    import numpy as np
    import pandas as pd

    src = (ROOT / "gui" / "ucmuon_gui.py").read_text()
    fn = next(n for n in ast.parse(src).body
              if isinstance(n, ast.FunctionDef) and n.name == "load_file")
    fn.decorator_list = []                       # drop @st.cache_data
    msgs = []
    st = types.SimpleNamespace(info=msgs.append, warning=msgs.append, error=msgs.append)
    ns = {"pd": pd, "np": np, "st": st}
    exec(compile(ast.Module([fn], []), "load_file", "exec"), ns)

    rows = [(13, 2132.896254804029), (11, 605.55574), (-13, 4758.424802264717)]
    path = tmp_path / "muons_for_phits.dat"           # format 10(1pd24.15), as phits_module
    path.write_text("".join(
        "".join(f"{v:24.15E}".replace("E", "D")
                for v in (kf, 793.0, 162.7, 0.0, -0.30, -0.20, -0.93, ek, 1.0, 0.0)) + "\n"
        for kf, ek in rows))
    df = ns["load_file"](str(path))
    assert len(df) == 3, msgs
    assert df["charge"].tolist() == [-1, -1, 1]
    assert np.allclose(df["E"], [2.132896254804029 + 0.105658, 0.60555574 + 0.000511,
                                 4.758424802264717 + 0.105658])


def _beam(n, e_gev, z0=0.0):
    import numpy as np
    return {"Ekin_MeV": np.full(n, e_gev * 1e3), "cx": np.zeros(n), "cy": np.zeros(n),
            "cz": -np.ones(n), "x": np.linspace(-100, 100, n), "y": np.zeros(n),
            "z": np.full(n, z0), "depth_axis": 2}


def _pdg_csda_exit(key, e_gev, x_gcm2):
    """Kinetic energy [MeV] after x_gcm2 of material `key` (PDG 2024 table),
    integrated step by step, independently of the drivers' range tables."""
    import numpy as np
    import ucmuon_stochastic_driver as U
    t, T = U._PROC_TABLES[key], U._PROC_T
    tot = t["ION"] + t["BREMS"] + t["PAIR"] + t["NUCL"]
    tk, h = e_gev * 1e3, 0.5
    for _ in range(int(x_gcm2 / h)):
        tk -= np.interp(np.log(tk), np.log(T), tot) * h
    return tk


def test_python_bethe_bloch_water_follows_the_pdg_table_with_and_without_ms():
    """Water and ice rescaled the rock table by 1.046 (Z/A 0.523; water's is
    0.555) and the MS-off path divided where it should multiply: 100 GeV
    muons through 100 m of water lost 5 % (MS on) and 13 % (MS off) too
    little energy. A custom material gave different losses with MS on and off."""
    import numpy as np
    import ucmuon_bb_driver as bb
    ref = _pdg_csda_exit("water", 100.0, 10000.0)
    for ms in (True, False):
        r = bb.transport_bb(_beam(200, 100.0), 100.0, dict(bb._BB_MAT[3]), ms_enable=ms,
                            rng=np.random.default_rng(1))
        assert abs(r["E_kin_f_MeV"].mean() / ref - 1.0) < 0.003, (ms, r["E_kin_f_MeV"].mean(), ref)
    custom = {"name": "Custom", "Z": 11.0, "A": 20.0, "rho": 2.0, "I_eV": 136.4,
              "X0_gcm2": 26.0, "a_scale": (11.0 / 20.0) / 0.5,
              "b_scale": (121.0 / 20.0) / (121.0 / 22.0)}
    e_on, e_off = (bb.transport_bb(_beam(200, 10.0), 10.0, custom, ms_enable=ms,
                                   rng=np.random.default_rng(1))["E_kin_f_MeV"].mean()
                   for ms in (True, False))
    assert abs(e_on / e_off - 1.0) < 0.003, (e_on, e_off)


def test_exit_plane_is_measured_from_the_source():
    """UCMuon-MC and the Python Bethe-Bloch wrote survivors at z = -depth even
    when the source plane is at z0 != 0 (the Fortran engines, MUSIC and PUMAS
    write z0 - depth); the GUI's underground filter then shifted inclined
    muons laterally by z0 tan(theta)."""
    import numpy as np
    import ucmuon_bb_driver as bb
    import ucmuon_stochastic_driver as U
    r = U.transport(_beam(50, 100.0, z0=1000.0), 20.0, 2.65, dict(U._MAT_DB[1]),
                    ms_enable=False, rng=np.random.default_rng(1))
    assert np.allclose(r["z_f"][r["alive"] == 1], 1000.0 - 2000.0)
    r = bb.transport_bb(_beam(50, 100.0, z0=1000.0), 20.0, dict(bb._BB_MAT[1]),
                        ms_enable=False, rng=np.random.default_rng(1))
    assert np.allclose(r["z_f"][r["alive"] == 1], 1000.0 - 2000.0)


def test_point_source_depth_axis_from_the_momentum(tmp_path):
    """A pencil beam (every muon at one point) gave the PROPOSAL and PUMAS
    drivers depth along X (argmin of three zero spreads): no muon survived."""
    import numpy as np
    import proposal_driver as pd_
    import ucmuon_pumas_driver as pu
    muons = [{"xs": 0.0, "ys": 0.0, "zs": 0.0, "px": 0.0, "py": 0.0, "pz": -100.0}] * 20
    assert pd_.detect_source_plane(muons) == 2
    f = tmp_path / "beam.dat"
    f.write_text("# header\n" + "".join(
        f"{i+1} 0.0 0.0 0.0 100.0 0.0 0.0 -100.0 0.0 0.0 100.06 1 1\n" for i in range(20)))
    assert pu._detect_depth_axis(str(f)) == 2


def test_multiple_scattering_is_highland_for_the_layer_at_any_step():
    """UCMuon-MC and the Python Bethe-Bloch applied Highland per step and added
    the steps in quadrature: the angles depended on the step size (~10 % per
    decade) and one radiation length of rock gave 22-30 % less than Highland
    for that layer. The variance is now accumulated along the path, so a layer
    up to one X0 gives Highland's theta0 whatever the step."""
    import numpy as np
    import ucmuon_bb_driver as bb
    import ucmuon_stochastic_driver as U
    n, depth_m = 20000, 26.54 / 2.65 / 100.0            # one X0 of Standard Rock
    p = np.sqrt(10.0 ** 2 - 0.1056584 ** 2)
    beam = _beam(n, 10.0 - 0.1056584)                    # 10 GeV total energy
    theta0 = 13.6e-3 / (p * p / 10.0)                    # Highland, t = 1
    for steps in (30, 3000):
        r = U.transport(dict(beam), depth_m, 2.65, dict(U._MAT_DB[1]), n_steps=steps,
                        delta_rays=False, rng=np.random.default_rng(1))
        rms = np.sqrt(np.mean(np.arctan2(r["cx_f"], -r["cz_f"]) ** 2))
        assert abs(rms / theta0 - 1.0) < 0.03, ("UCMuon-MC", steps, rms, theta0)
    for steps in (50, 1000):
        r = bb.transport_bb(dict(beam), depth_m, dict(bb._BB_MAT[1]), n_steps=steps,
                            rng=np.random.default_rng(1))
        rms = np.sqrt(np.mean(np.arctan2(r["cx_f"], -r["cz_f"]) ** 2))
        assert abs(rms / theta0 - 1.0) < 0.03, ("Bethe-Bloch", steps, rms, theta0)


def test_terrain_ray_tracer_measures_the_rock_exactly():
    """The terrain ray tracer counted whole steps below the DEM: a detector
    100 m under a flat surface (underground deployment) got 50 m of rock
    vertically and 100 m instead of 141 m at 45 deg with the default 50 m
    step. Boundaries are now located by bisection."""
    import numpy as np
    import pytest
    pytest.importorskip("rasterio")
    import gui_geometry_engine as G
    import ucmuon_terrain_driver as td
    s = G.build_flat_slab(40.8, 14.4, -100.0, extent_km=3.0, res_m=50.0, z_surface=0.0)
    for ze in (0.0, 45.0, 70.0):
        _, L, sky = td.ray_overburden(s.elev, s.transform, 40.8, 14.4, -100.0, 37.0, ze,
                                      2.65, 50.0, underground=True)
        assert not sky and abs(L - 100.0 / np.cos(np.radians(ze))) < 0.01, (ze, L)


def test_dem_loaders_place_projected_and_asc_grids_correctly(tmp_path):
    """A DEM in a projected CRS was taken as degrees (GeoTIFF) or converted
    with a formula that ignored the UTM zone (XYZ, ASC): the detector fell
    outside it and every direction was open sky. An ASC in xllcenter
    convention was shifted half a cell, and one without the optional
    NODATA_value line lost its first row."""
    import numpy as np
    import pytest
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin
    from rasterio.warp import transform as wtransform
    import ucmuon_terrain_driver as td
    # A 600 m cone (apex 1000 m) on a UTM 33N grid of 20 m pixels
    (ex,), (ny,) = wtransform("EPSG:4326", "EPSG:32633", [14.4], [40.8])
    x = ex - 3000.0 + 20.0 * (np.arange(300) + 0.5)
    y = ny + 3000.0 - 20.0 * (np.arange(300) + 0.5)
    r = np.hypot(*np.meshgrid(x - ex, y - ny))
    z = np.maximum(1000.0 - r, 400.0).astype("float32")
    tif = tmp_path / "cone_utm.tif"
    with rasterio.open(tif, "w", driver="GTiff", height=300, width=300, count=1, dtype="float32",
                       crs="EPSG:32633", transform=from_origin(ex - 3000.0, ny + 3000.0, 20.0, 20.0)) as ds:
        ds.write(z, 1)
    elev, tr = td.load_dem(str(tif))
    assert abs(td.dem_elevation_at(elev, tr, 40.8, 14.4) - 1000.0) < 25.0
    xyz = tmp_path / "cone_utm.xyz"
    xx, yy = np.meshgrid(x, y)
    np.savetxt(xyz, np.column_stack([xx.ravel(), yy.ravel(), z.ravel()]), fmt="%.2f")
    with pytest.raises(ValueError, match="projected"):
        td.load_dem(str(xyz))
    # Geographic ASC: centre convention and a five-line header
    asc = tmp_path / "g.asc"
    asc.write_text("ncols 4\nnrows 3\nxllcenter 14.005\nyllcenter 40.005\ncellsize 0.01\n"
                   "1 2 3 4\n5 6 7 8\n9 10 11 12\n")
    elev, tr = td.load_dem(str(asc))
    assert elev.shape == (3, 4) and elev[0, 0] == 1.0
    assert abs(tr.c - 14.0) < 1e-9 and abs(tr.f - 40.03) < 1e-9


def test_backward_mc_materials_follow_their_pdg_tables():
    """The backward MC rescaled the rock range by a_scale where it should
    divide (and used one factor for all losses): water rates came out 29 %
    too high at 50 m, iron 25 % too low at 10 m."""
    import numpy as np
    import ucmuon_backward_mc as B
    import ucmuon_stochastic_driver as U

    def ref_rate(depth, rho, key):
        """CSDA rate (P_surv = 1) on the material's own PDG table, on the
        grids backward_mc_flux uses."""
        t, T = U._PROC_TABLES[key], U._PROC_T
        dedx = t["ION"] + t["BREMS"] + t["PAIR"] + t["NUCL"]
        tf = np.logspace(np.log10(T[0]), np.log10(2.0e6), 4000)
        df = np.exp(np.interp(np.log(tf), np.log(T), np.log(dedx)))
        rf = t["R10"] + np.concatenate(([0.0], np.cumsum(0.5 * (1 / df[1:] + 1 / df[:-1]) * np.diff(tf))))
        th_e = np.linspace(0.0, np.radians(60.0), 11)
        th, dth = 0.5 * (th_e[1:] + th_e[:-1]), np.diff(th_e)
        e_det = np.logspace(3.0, np.log10(2.0e6), 40)
        dlne = np.log(2000.0) / 39
        rate = 0.0
        for i, e in enumerate(e_det):
            de = e / 1000.0 * dlne * (0.5 if i in (0, 39) else 1.0)
            for a, da in zip(th, dth):
                r_s = np.interp(e, tf, rf) + depth * 100.0 * rho / max(np.cos(a), 0.02)
                if r_s > rf[-1]:
                    continue
                e_s = np.interp(r_s, rf, tf)
                if e_s / 1000.0 > 2000.0:
                    continue
                jac = np.interp(e_s, tf, df) / np.interp(e, tf, df)
                rate += (float(B._flux_surface(e_s / 1000.0, a, 3)) * jac
                         * 2.0 * np.pi * np.sin(a) * da * np.cos(a) * de)
        return rate

    for depth, rho, key, mat_id in ((50.0, 1.0, "water", 2), (10.0, 7.874, "iron", 4)):
        got = B.backward_mc_flux(depth, rho, mat_id, 3, E_min_GeV=1.0, E_max_GeV=2000.0,
                                 theta_max_deg=60.0, n_E=40, n_theta=10, mode=0)["rate_m2_s"]
        ref = ref_rate(depth, rho, key)
        assert abs(got / ref - 1.0) < 0.01, (key, got, ref)


def test_transmission_map_with_full_precision_angles_loads(tmp_path):
    """Grid keys were the parsed angles but lookups rounded them to 8
    decimals, so a map written at full precision (np.savetxt default) on the
    default 18 zenith bins loaded with 576 of 648 cells NaN, silently."""
    import numpy as np
    import ucmuon_density_analysis as D
    az = np.linspace(0, 360, 37)[:-1] + 5.0
    el = 90.0 - (np.arange(18) + 0.5) * 85.0 / 18.0
    f = tmp_path / "t.dat"
    with open(f, "w") as fh:
        fh.write("# Density: 2.65 g/cm3\n")
        np.savetxt(fh, [(a, e, 0.5) for a in az for e in el])
    _, _, T, _ = D.load_transmission_map(str(f))
    assert T.shape == (36, 18) and not np.isnan(T).any()


def test_stl_overburden_from_inside_the_mesh_and_errors_surface():
    """STL overburden: a detector inside the mesh (one exit hit per ray) got
    open sky, because the hit count was tested first, and any error (trimesh
    without rtree) was turned into open sky too."""
    import types
    import numpy as np
    import pytest
    import gui_csg_engine as C
    lo, hi = np.array([-100.0] * 3), np.array([100.0] * 3)

    def hits(ray_origins, ray_directions, **_):          # a 200 m cube
        o, d = np.asarray(ray_origins[0], float), np.asarray(ray_directions[0], float)
        t = [tt for tt in np.concatenate(((lo - o) / d, (hi - o) / d))
             if np.isfinite(tt) and tt >= 0 and np.all(np.abs(o + tt * d) <= 100.0 + 1e-9)]
        return np.array([o + tt * d for tt in t]).reshape(-1, 3), None, None

    g = C.CSGGeometry()
    g._stl_density = 2.65
    g._mesh = types.SimpleNamespace(
        ray=types.SimpleNamespace(intersects_location=hits),
        contains=lambda p: np.array([bool(np.all((np.asarray(q) > lo) & (np.asarray(q) < hi))) for q in p]))
    with np.errstate(divide="ignore", invalid="ignore"):
        assert abs(g.ray_overburden(0.0, 0.0, np.zeros(3))[1] - 100.0) < 1e-6
        assert abs(g.ray_overburden(0.0, 0.0, np.array([0.0, 0.0, -150.0]))[1] - 200.0) < 1e-6

    def broken(**_):
        raise ModuleNotFoundError("No module named 'rtree'")
    g._mesh.ray = types.SimpleNamespace(intersects_location=broken)
    with pytest.raises(ModuleNotFoundError):
        g.ray_overburden(0.0, 0.0, np.zeros(3))


def test_source_panel_never_advises_cutting_theta_max_for_scattering():
    """For the paper's underground example (r = 10 cm cylinder 90 m down,
    E_min 45 GeV) the source panel said "Reduce to 0.7 deg" (its MCS limit)
    and used a 102 cm MCS margin where the Generator's warning asks 231 cm."""
    import gui_source_optimizer as O
    import mcs_margin as M
    det = {"shape": 1, "ax": 0, "ay": 0, "az": -9000.0, "bx": 0, "by": 0, "bz": -9050.0,
           "r": 10.0, "margin": 0.0}
    opt = O._compute([det], 45.0, 70.0, 2.65, 26.7)
    assert opt["theta_rec"] >= 70.0
    want = M.suggested_margin(90.5, 2.65, 45.0, 0.0, 26.7)["margin_cm"]
    assert abs(opt["mcs_margins"][0] - want) < 1e-6


# ── Part 4 (GUI) ─────────────────────────────────────────────────────────────
from conftest import HAVE_GENERATOR, HAVE_RASTERIO, top_tab   # noqa: E402

_ROCK_WITH_CAVITY = """[ Surface ]
  1  pz  -5000.0
  2  pz      0.0
  3  cz  10000.0
  8  s   0.0 0.0 -2000.0 500.0
  9  rpp -200 200 -200 200 -3500 -3000
  7  so  100000.0
[ Cell ]
 100  1  -2.65   -3 1 -2 #102 #103
 102  0          -8
 103  2  -1.03   -9
 104  0          -7 #100 #102 #103
 105 -1           7
"""


def test_geant4_export_is_the_new_runs_muons(new_app, repo_copy):
    """The Geant4 export read the run's file through the cached loader with
    the path as the only key: a second run to the same file exported the
    first run's muons (200 rows after a 300-muon run)."""
    import time
    import pytest
    if not HAVE_GENERATOR:
        pytest.skip("needs bin/ucmuon_gen_omp (make ucmuon_gen_omp)")
    at = new_app(mode="Advanced")
    widget(at, "checkbox", "usedetector").uncheck().run()
    widget(at, "checkbox", "savegeant4").check().run()
    g4 = repo_copy / "output" / "muons_geant4.txt"
    for n in (200, 300):
        widget(at, "number_input", "nmuonsgen").set_value(n).run()
        [b for b in at.button if "Run UCMuon Surface Generator" in str(b.label)][0].click().run()
        for _ in range(300):
            time.sleep(1)
            at.run()
            if at.session_state["gen_ntry"]:
                break
        at.run()
        rows = [ln for ln in g4.read_text().splitlines() if ln and not ln.startswith("#")]
        assert len(rows) == n, (n, len(rows))
    assert not exceptions(at)


def test_margin_helper_uses_the_transport_material(new_app):
    """The MCS margin helper took the slant-path threshold at 2.65 g/cm3, and
    standard rock's X0 and range table for every material: through 100 m of
    water it suggested 96 cm where the Generator's warning asks 376 cm (281 cm
    with water's X0 and range)."""
    import mcs_margin as M
    at = new_app(mode="Advanced", state={"musicpreset": "Water", "music_preset_prev": "Water",
                                         "music_rho": 1.0, "music_rad": 36.08, "emin": 1.0,
                                         "hlp_mcs_depth": 100.0, "hlp_mcs_ze": 0.0})
    assert not exceptions(at)
    tab = [t for t in at.tabs if t.label == "🎯  MCS margin"][0]
    shown = float([m for m in tab.metric if m.label == "Suggested margin"][0].value.split()[0])
    want = M.suggested_margin(100.0, 1.0, 1.0, 0.0, 36.08, material="water")["margin_cm"]
    assert abs(shown - want) < 0.1 + 1e-3 * want, (shown, want)


def test_basic_margin_uses_the_materials_radiation_length(new_app):
    """Basic's safety margin used standard rock's X0 and range table for every
    material: for iron (X0 13.84 g/cm2) it came out 11-13 % below the
    scattering UCMuon-MC produces, not the envelope the 2 sigma rule needs."""
    import mcs_margin as M
    at = new_app(state={"musicpreset": "Iron", "music_preset_prev": "Iron", "music_rho": 7.874})
    assert not exceptions(at)
    card = {m.label: m.value for m in top_tab(at, "Generator").metric}
    want = M.suggested_margin(31.0, 7.874, 1.0, 0.0, 13.84, material="iron")["margin_cm"]  # default detector,
    rock = M.suggested_margin(31.0, 7.874, 1.0, 0.0)["margin_cm"]          # 31 m at its bottom
    assert abs(want - rock) > 10.0
    assert card["Safety margin"] == f"{round(want, 1):.0f} cm", (card["Safety margin"], want, rock)


def test_csg_cells_with_complements_and_parentheses(tmp_path):
    """The PHITS cell parser dropped '#N' complements and flattened
    parentheses: a void cavity cut out of the rock with #102 read as rock
    (39.5 m of rock above z = -40 m instead of 30 m), '-3 (1 : 2)' was
    read as '(-3 1) : 2', and a cell of an unsupported body (BOX) filled
    every point no earlier cell claimed."""
    import numpy as np
    import gui_csg_engine as C
    f = tmp_path / "cavity.inp"
    f.write_text(_ROCK_WITH_CAVITY)
    g = C.load_phits(str(f))
    assert g.density_at(np.array([0.0, 0.0, -20.0])) == 0.0          # in the cavity
    assert g.density_at(np.array([0.0, 0.0, -32.0])) == 1.03         # in the detector
    _, slant, sky = g.ray_overburden(0.0, 0.0, np.array([0.0, 0.0, -40.0]), 0.05, 1e4)
    assert not sky and abs(slant - 30.0) < 0.1, slant
    tree = C._parse_cell_expr("-3 (1 : 2)")
    assert tree == ("and", [("s", -1, 3), ("or", [("s", 1, 1), ("s", 1, 2)])])
    f2 = tmp_path / "box.inp"
    f2.write_text("[ Surface ]\n 1 box 0 0 0 100 0 0 0 100 0 0 0 100\n 2 so 1000\n"
                  "[ Cell ]\n 10 1 -2.65 -1\n 11 0 -2 #10\n 12 -1 2\n")
    g2 = C.load_phits(str(f2))
    assert any("Surface 1 (box) is not supported" in w for w in g2.warnings)
    # a cell bounded only by the unsupported box filled every point, the sky too
    assert g2.density_at(np.array([5.0, 5.0, 5.0])) == 0.0


def test_csg_preview_traces_from_the_detector_cell(tmp_path):
    """The CSG overburden preview, single-ray check and 3D marker read
    _csg_det_pos, which nothing set: they traced from the geometry's origin
    (open sky above a detector 30 m down). It is now a point inside the
    first selected detector cell."""
    import numpy as np
    from streamlit.testing.v1 import AppTest
    f = tmp_path / "cavity.inp"
    f.write_text(_ROCK_WITH_CAVITY)

    def page(deck, gui):
        import sys
        sys.path.insert(0, gui)
        import gui_csg_engine as C
        import gui_csg_transport as T
        T.render_detector_cell_selector(C.load_phits(deck))

    at = AppTest.from_function(page, args=(str(f), str(ROOT / "gui"))).run()
    ms = [m for m in at.multiselect if m.key == "csg_det_cell_ids"][0]
    ms.set_value([o for o in ms.options if o.startswith("Cell 103 ")]).run()
    assert "_csg_det_pos" in at.session_state
    assert np.allclose(at.session_state["_csg_det_pos"], [0.0, 0.0, -32.5])


def test_mcs_estimator_counts_detector_hits_after_transport():
    """The MCS estimator tested the scattered muons against the margin-
    inflated detector, at the surface momentum, counting muons that range
    out: for the Basic default detector it expected 5620 hits of 20000
    selected muons where transport gave 113."""
    import numpy as np
    import pandas as pd
    import gui_mcs_estimator as E
    det = {"shape": 1, "ax": 0.0, "ay": 0.0, "az": -3100.0, "bx": 0.0, "by": 0.0,
           "bz": -3000.0, "r": 50.0, "margin": 126.0}
    rng = np.random.default_rng(1)
    n = 4000
    rr = 176.0 * np.sqrt(rng.uniform(0, 1, n))
    a = rng.uniform(0, 2 * np.pi, n)

    def sample(e):
        p = np.sqrt(e * e - 0.10566 ** 2)
        return pd.DataFrame(dict(x=rr * np.cos(a), y=rr * np.sin(a), z=0.0, p=p, px=0.0, py=0.0,
                                 pz=-p, theta=0.0, phi=0.0))
    stop = E.estimate_mcs_acceptance(sample(5.0), [det], 2.65, 26.54, n_mc=200)
    assert stop["n_after_mcs"] == 0.0                      # 5 GeV stops in 30 m of rock
    hard = E.estimate_mcs_acceptance(sample(1000.0), [det], 2.65, 26.54, n_mc=200)
    want = n * (50.0 / 176.0) ** 2                        # the detector's share of the disc
    assert abs(hard["n_after_mcs"] / want - 1.0) < 0.1, (hard["n_after_mcs"], want)
    assert hard["acceptance_frac"] > 0.99


def test_terrain_warns_about_muons_beyond_the_zenith_bins(new_app, surface_run):
    """Muons steeper than the terrain grid's maximum zenith were dropped
    from the transport and its output without a word."""
    import numpy as np
    import pytest
    if not HAVE_RASTERIO:
        pytest.skip("needs rasterio")
    f, state = surface_run
    at = new_app(mode="Advanced", state={**state, "terrain_zemax": 40.0, "terrain_naz": 8,
                                         "terrain_nze": 4})
    [b for b in at.button if b.key == "terrain_run_btn"][0].click().run()
    assert not exceptions(at)
    d = np.loadtxt(f, comments="#")
    n_beyond = int(np.sum(np.degrees(np.arccos(np.clip(-d[:, 7] / d[:, 4], -1, 1))) >= 40.0))
    assert n_beyond > 0
    assert any(f"{n_beyond:,} of {len(d):,} surface muons" in w.value for w in at.warning), \
        [w.value[:120] for w in at.warning]


def _beam_df(theta_deg, energies, spread=True, seed=1):
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(seed)
    n = len(energies)
    th = np.radians(np.full(n, theta_deg) if not spread else rng.uniform(theta_deg - 2, theta_deg + 2, n))
    ph = rng.uniform(0, 2 * np.pi, n) if spread else np.zeros(n)
    E = np.asarray(energies, float)
    return pd.DataFrame(dict(EventID=np.arange(1, n + 1), x=rng.uniform(-500, 500, n),
                             y=rng.uniform(-500, 500, n), z=0.0, E=E, p=np.sqrt(E * E - 0.10566 ** 2),
                             theta=th, phi=ph, charge=1, cx=np.sin(th) * np.cos(ph),
                             cy=np.sin(th) * np.sin(ph), cz=-np.cos(th)))


def test_terrain_transport_puts_the_slant_path_in_front_of_each_muon():
    """The terrain forward transport set the engine depth to the slant
    opacity divided by cos θ instead of multiplied: 1/cos²θ too much rock.
    Through a flat 20 m slab seen at 60°, 1180 of 3000 muons survived where
    2391 do; at 80°, none."""
    import numpy as np
    import gui_terrain_engine as TE
    import ucmuon_stochastic_driver as U
    n, D = 3000, 20.0
    df = _beam_df(60.0, 10 ** np.random.default_rng(3).uniform(1, 2.7, n))
    az_c = np.arange(4) * 90 + 45.0
    ze_c = (np.arange(17) + 0.5) * 5.0
    ob = (2.65 * D * 100 / np.cos(np.radians(ze_c)))[None, :].repeat(4, 0)   # a flat slab
    out = TE.transport_stochastic_terrain(df, ob, az_c, ze_c, np.zeros((4, 17), bool), 2.65,
                                          4, 17, 85.0, 0.05, True, str(ROOT / "gui"))
    Ek = df.E.to_numpy() - U.M_MU_GEV
    mu = {k: df[k].to_numpy() for k in ("EventID", "theta", "phi", "charge", "cx", "cy", "cz")}
    mu.update(x=np.zeros(n), y=np.zeros(n), z=np.zeros(n), E_tot_GeV=df.E.to_numpy(),
              Ekin_GeV=Ek, Ekin_MeV=Ek * 1000.0)
    mat = dict(U._MAT_DB[1])
    mat["X0_cm"] = 26.54 / 2.65
    direct = int(U.transport(mu, D, 2.65, mat, rng=np.random.default_rng(2))["alive"].sum())
    terrain = int((out.alive == 1).sum())
    assert abs(terrain - direct) < 4 * np.sqrt(direct), (terrain, direct)


def test_terrain_bin_with_one_muon_keeps_its_depth_axis(repo_copy):
    """A direction bin holding one muon was written as a one-muon file, which
    the engines read as a point source: they took the depth axis from the
    momentum, so a 20 GeV muon at 81° (range 35 m) crossed 10 m of rock
    instead of 64 m and survived (Python Bethe-Bloch). The bin files were
    also created in the project root, and the engines' extra files stayed."""
    import numpy as np
    import gui_terrain_engine as TE

    def build_bb(cfg):
        return "\n".join([cfg["infile"], cfg["outfile"], "1", "13", str(cfg["depth_m"]),
                          f"1 {cfg['rho']}", "1"]) + "\n"
    ze_c = np.array([9.0, 27.0, 45.0, 63.0, 81.0])
    az_c = np.array([45.0, 135.0, 225.0, 315.0])
    ob = np.full((4, 5), 10.0 * 100 * 2.65 / np.cos(np.radians(81.0)))   # 64 m of rock at 81°
    before = set(p.name for p in repo_copy.iterdir())
    for E, alive in ((20.0, 0), (200.0, 1)):
        out = TE.transport_batched_terrain(
            _beam_df(81.0, [E], spread=False), ob, az_c, ze_c, np.zeros((4, 5), bool), 2.65,
            4, 5, 90.0, "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS",
            {"ncols": 13, "rho": 2.65}, ROOT / "gui", repo_copy, None, build_bb, None)
        assert out.alive.tolist() == [alive], (E, out.alive.tolist())
    # and nothing left in the project root (up to 1.3.2 the engines'
    # _stopped and _timing files of every bin stayed there)
    assert set(p.name for p in repo_copy.iterdir()) == before


def test_csg_run_transports_with_the_selected_engine(tmp_path):
    """The CSG run applied straight-line CSDA energy loss whatever engine was
    selected (no stochastic losses, no scattering). Through a 30 m rock layer
    with UCMuon-MC it must now survive and scatter as the engine does through
    a 30 m slab."""
    import numpy as np
    import gui_csg_engine as C
    import gui_csg_transport as T
    import ucmuon_stochastic_driver as U
    f = tmp_path / "slab.inp"
    f.write_text("[ Surface ]\n 1 pz -3000\n 2 pz 0\n 3 cz 100000\n 4 pz -3200\n 5 c/z 0 0 300\n"
                 " 6 pz -3100\n 7 so 1000000\n[ Cell ]\n 100 1 -2.65 -3 1 -2\n"
                 " 101 2 -1.03 -5 4 -6\n 104 0 -7 #100 #101\n 105 -1 7\n")
    g = C.load_phits(str(f))
    n = 3000
    rng = np.random.default_rng(5)
    df = _beam_df(20.0, 10 ** rng.uniform(np.log10(15), np.log10(300), n))
    t = 3100.0 / np.cos(df.theta)                                 # aim at the cell's top face
    rr, aa = 290 * np.sqrt(rng.uniform(0, 1, n)), rng.uniform(0, 2 * np.pi, n)
    df["x"], df["y"] = rr * np.cos(aa) - t * df.cx, rr * np.sin(aa) - t * df.cy
    r = T.transport_muons_csg_engine(df, g, [101], "UCMuon Stochastic", rho_ref=2.65,
                                     write_phits_dump=False, script_dir=str(ROOT / "gui"))[101]
    Ek = df.E.to_numpy() - U.M_MU_GEV
    mu = {k: df[k].to_numpy() for k in ("EventID", "theta", "phi", "charge", "cx", "cy", "cz")}
    mu.update(x=np.zeros(n), y=np.zeros(n), z=np.zeros(n), E_tot_GeV=df.E.to_numpy(),
              Ekin_GeV=Ek, Ekin_MeV=Ek * 1000.0)
    mat = dict(U._MAT_DB[1])
    mat["X0_cm"] = 26.54 / 2.65
    direct = int(U.transport(mu, 30.0, 2.65, mat, rng=np.random.default_rng(2))["alive"].sum())
    assert abs(int(r.alive.sum()) - direct) < 4 * np.sqrt(n - direct) + 5, (int(r.alive.sum()), direct)
    a = r.alive == 1
    cosang = (r[["cx", "cy", "cz"]].to_numpy() * df[["cx", "cy", "cz"]].to_numpy()).sum(1)
    assert np.degrees(np.arccos(np.clip(cosang[a], -1, 1))).mean() > 0.1   # scattered


def test_csg_engine_run_with_repeated_event_ids(tmp_path):
    """Files merged from MPI rank files repeated EventIDs (every rank counted
    from 1 up to 1.3.2), and the CSG run with a subprocess engine matched the
    engine's output back by EventID, which pandas refuses for repeated
    labels. The muons of a depth bin are numbered for the engine now."""
    import numpy as np
    import pandas as pd
    import gui_csg_engine as C
    import gui_csg_transport as T
    f = tmp_path / "slab.inp"
    f.write_text("[ Surface ]\n 1 pz -3000\n 2 pz 0\n 3 cz 100000\n 4 pz -3200\n 5 c/z 0 0 300\n"
                 " 6 pz -3100\n 7 so 1000000\n[ Cell ]\n 100 1 -2.65 -3 1 -2\n"
                 " 101 2 -1.03 -5 4 -6\n 104 0 -7 #100 #101\n 105 -1 7\n")
    g = C.load_phits(str(f))
    n = 400
    rng = np.random.default_rng(3)
    df = _beam_df(0.0, np.linspace(10.0, 200.0, n), spread=False)    # vertical
    rr, aa = 250 * np.sqrt(rng.uniform(0, 1, n)), rng.uniform(0, 2 * np.pi, n)
    df["x"], df["y"] = rr * np.cos(aa), rr * np.sin(aa)
    df["EventID"] = np.arange(n) % 100 + 1                            # each ID four times

    def run_bin(sub, depth_m):
        """A stand-in engine: loses 1 GeV per m, answers in reverse order."""
        E = sub.E.to_numpy() - depth_m
        out = pd.DataFrame(dict(
            EventID=sub.EventID, xs=sub.x, ys=sub.y, zs=0.0, Es=sub.E, theta_s=sub.theta,
            phi_s=sub.phi, charge=sub.charge, alive=(E > 0.2).astype(int), x=sub.x, y=sub.y,
            z=-depth_m * 100.0, E=np.where(E > 0.2, E, 0.0), cx=sub.cx, cy=sub.cy, cz=sub.cz,
            theta=sub.theta, phi=sub.phi))
        return out.iloc[::-1].reset_index(drop=True)

    r = T.transport_muons_csg_engine(df, g, [101], "MUSIC", rho_ref=2.65, run_bin=run_bin,
                                     write_phits_dump=False, script_dir=str(ROOT / "gui"))[101]
    assert len(r) == n and r.EventID.tolist() == df.EventID.tolist()
    loss = df.E.to_numpy() - r.E.to_numpy()                           # each muon's own
    a = r.alive.to_numpy() == 1
    assert a.sum() > 0.8 * n and np.ptp(loss[a]) < 1e-9, (a.sum(), np.ptp(loss[a]))
    assert np.all(df.E.to_numpy()[~a] < df.E.to_numpy()[a].min())    # the low-energy ones stop


def test_density_inversion_between_library_densities():
    """The library inversion interpolated T linearly between library
    densities; T falls almost exponentially with ρ, so ρ̂ came out high
    between the nodes (+2.6 % at 1.75 g/cm3 on the Vesuvius library). It
    interpolates ln T now."""
    import numpy as np
    import ucmuon_density_analysis as da
    rhos = (1.5, 2.0, 2.5, 2.65, 3.0)
    lam = np.linspace(0.4, 1.2, 12).reshape(3, 4)              # per-pixel attenuation
    lib = {r: np.exp(-r / lam) for r in rhos}
    for true in (1.75, 2.25, 2.8):
        rho_map, _, status = da.invert_density_map(np.exp(-true / lam), lib, min_sensitivity=0.0)
        assert np.all(status == 0)
        assert np.allclose(rho_map, true, rtol=1e-6), (true, rho_map.ravel()[:3])


def test_range_cut_muons_stop_where_the_stepped_ones_do():
    """A muon that cannot cross the slab is cut and placed at its range,
    which was rock's range times a_scale, the inverse of the material
    correction: 10 GeV muons stopped at 42 m in 45 m of water (stepped) but
    at 51 m in 200 m (cut); in iron at 7.0 m and 5.8 m."""
    import numpy as np
    import ucmuon_stochastic_driver as U
    n = 2000
    Ek = 10.0 - U.M_MU_GEV
    mu = dict(EventID=np.arange(n), x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
              theta=np.zeros(n), phi=np.zeros(n), E_tot_GeV=np.full(n, 10.0),
              Ekin_GeV=np.full(n, Ek), Ekin_MeV=np.full(n, Ek * 1000), charge=np.ones(n, int),
              cx=np.zeros(n), cy=np.zeros(n), cz=-np.ones(n))
    for key, rho, x0, thin, thick in ((2, 1.0, 36.08, 45.0, 200.0), (4, 7.874, 13.84, 8.0, 30.0)):
        mat = dict(U._MAT_DB[key])
        mat["X0_cm"] = x0 / rho
        depth = []
        for d in (thin, thick):
            r = U.transport(mu, d, rho, mat, rng=np.random.default_rng(4))
            depth.append(-r["z_stop"][r["alive"] == 0].mean() / 100.0)
        assert abs(depth[1] / depth[0] - 1.0) < 0.03, (mat["name"], depth)
