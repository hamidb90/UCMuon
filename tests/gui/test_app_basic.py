"""The GUI in Basic and Advanced mode (Streamlit AppTest, run in a copy of the
repository; see conftest.py)."""
import json

import pytest

from conftest import _TAB_LABELS, HAVE_GENERATOR, HAVE_RASTERIO, exceptions, top_tab, widget


# ── start-up ─────────────────────────────────────────────────────────────────
def test_fresh_install_starts_in_basic(new_app, repo_copy):
    at = new_app()
    assert not exceptions(at)
    assert at.session_state["ui_mode"] == "Basic"
    toggle = widget(at, "toggle", "ui_mode_toggle")                  # page header
    assert (toggle.label, toggle.value) == ("Advanced mode", False)
    assert widget(at, "selectbox", "gen_spectrum_mode").value == 4       # Guan 2015
    assert (widget(at, "number_input", "emin").value,
            widget(at, "number_input", "emax").value) == (1.0, 2500.0)
    assert widget(at, "selectbox", "angularmode").value == 6
    assert widget(at, "checkbox", "bd_use_detector").value is True
    card = {m.label: m.value for m in top_tab(at, "Generator").metric}
    assert {"Safety margin", "Source"} <= set(card)
    assert (repo_copy / "ucmuon_autosave.json").exists()
    # Basic: three tabs, all rendered (an exception or st.stop() in one would
    # blank the later ones)
    assert [t.label for t in at.tabs] == [_TAB_LABELS[n] for n in ("Generator", "Transport", "Results")]
    assert len(top_tab(at, "Results").markdown) + len(top_tab(at, "Results").info) > 0


def test_fresh_install_detector_is_underground(new_app):
    """Up to 1.3.1 the default detector was a 90 m tall, 5 cm column reaching
    up to the ground: Basic's first run took 100 000 hits (13 min) and then
    transported through no rock. Now its top face is 30 m deep, 10 000 hits."""
    at = new_app()
    assert not exceptions(at)
    assert (widget(at, "number_input", "bd_top").value,
            widget(at, "number_input", "bd_height").value,
            widget(at, "number_input", "bd_radius").value) == (30.0, 1.0, 0.5)
    assert widget(at, "number_input", "nmuonsgen").value == 10_000
    card = {m.label: m.value for m in top_tab(at, "Generator").metric}
    assert (card["Safety margin"], card["Source"]) == ("126 cm", "R = 380 m")


def test_fresh_install_advanced_renders(new_app):
    at = new_app(mode="Advanced")
    assert not exceptions(at)
    for key in ("gen_workflow", "density_mode"):
        assert widget(at, "radio", key) is not None
    for name in ("Terrain", "Density", "Config"):        # every tab rendered
        tab = top_tab(at, name)
        assert len(tab.button) + len(tab.radio) + len(tab.text_area) > 0, name


def test_guaranteed_hit_default_matches_the_standard_one(new_app):
    """Both workflows share the key nmuonsgen; 1.3.2 changed only the
    Standard default to 10 000 and left guaranteed-hit at 100 000."""
    at = new_app(mode="Advanced")
    widget(at, "radio", "gen_workflow").set_value("DAS-REM").run()
    assert not exceptions(at)
    assert widget(at, "number_input", "nmuonsgen").value == 10_000


def test_old_autosave_without_mode_opens_in_advanced(new_app, repo_copy):
    (repo_copy / "ucmuon_autosave.json").write_text(json.dumps({"emin": 3.0}))
    at = new_app(keep_autosave=True)
    assert at.session_state["ui_mode"] == "Advanced"
    assert widget(at, "toggle", "ui_mode_toggle").value is True


# ── Generator ────────────────────────────────────────────────────────────────
def test_basic_hides_generator_detail(new_app):
    at = new_app()
    for kind, key in (("radio", "gen_workflow"), ("checkbox", "gen_mono"),
                      ("number_input", "disk_r"), ("checkbox", "usedetector")):
        assert widget(at, kind, key) is None, key


def test_basic_output_files_and_exports(new_app):
    at = new_app()
    for kind, key in (("text_input", "outputall"), ("text_input", "outputsel"),
                      ("checkbox", "saveall"), ("checkbox", "savephits"),
                      ("checkbox", "savegeant4"), ("slider", "nthreads")):
        assert widget(at, kind, key) is not None, key
    widget(at, "checkbox", "savegeant4").check().run()
    assert widget(at, "radio", "g4_fmt_radio") is not None
    widget(at, "checkbox", "bd_use_detector").uncheck().run()
    assert widget(at, "text_input", "outputsel") is None      # no hits file without a detector


def test_basic_energy_range_and_angular_mode_are_inputs(new_app):
    at = new_app()
    widget(at, "selectbox", "gen_spectrum_mode").set_value(6).run()    # Gaisser
    assert (at.session_state["emin"], at.session_state["emax"],
            at.session_state["thetamax"]) == (100.0, 2500.0, 60.0)
    widget(at, "number_input", "emin").set_value(20.0).run()           # the user's own
    assert at.session_state["emin"] == 20.0
    [b for b in at.button if b.key == "btn_basic_rec_range"][0].click().run()
    assert at.session_state["emin"] == 100.0
    ang = widget(at, "selectbox", "angularmode")
    assert len(ang.options) == 6 and any("Vertical" in o for o in ang.options) \
        and any("cos²θ" in o for o in ang.options)
    ang.set_value(1).run()                                             # vertical
    assert not exceptions(at) and widget(at, "slider", "thetamax") is None
    card = {m.label: m.value for m in top_tab(at, "Generator").metric}
    assert card["Source"] == "R = 10 m"     # vertical: the detector's footprint only


def test_basic_without_a_detector(new_app, repo_copy):
    at = new_app()
    widget(at, "checkbox", "bd_use_detector").uncheck().run()
    assert not exceptions(at)
    assert widget(at, "number_input", "bd_top") is None
    assert widget(at, "number_input", "nmuonsgen").label == "Muons to generate"
    assert "Safety margin" not in {m.label for m in top_tab(at, "Generator").metric}
    widget(at, "number_input", "bd_src_size").set_value(50.0).run()
    saved = json.loads((repo_copy / "ucmuon_autosave.json").read_text())
    assert (saved["bd_use_detector"], saved["bd_src_size"]) == (False, 50.0)
    assert widget(at, "checkbox", "usedetector") is None     # Advanced's own filter untouched


def test_basic_detector_inputs_write_the_shared_keys(new_app):
    at = new_app()
    widget(at, "number_input", "bd_top").set_value(10.0).run()
    widget(at, "number_input", "bd_height").set_value(1.0).run()
    widget(at, "number_input", "bd_radius").set_value(0.5).run()
    ss = at.session_state
    assert (ss["az0"], ss["bz0"], ss["rr0"]) == (-1100.0, -1000.0, 50.0)
    card = {m.label: m.value for m in top_tab(at, "Generator").metric}
    assert card["Source"] == "R = 140 m"          # basic_mode.required_source, 2σ margin
    widget(at, "selectbox", "bd_shape").set_value(2).run()
    assert ss["xx0"] - ss["xn0"] == pytest.approx(100.0)      # box of 2 r


def test_power_law_keeps_the_users_range_in_basic(new_app):
    at = new_app()
    widget(at, "selectbox", "gen_spectrum_mode").set_value(2).run()
    assert widget(at, "number_input", "emin") is not None


def test_parma_shows_its_site_in_basic(new_app):
    at = new_app()
    widget(at, "selectbox", "gen_spectrum_mode").set_value(3).run()
    for key in ("parma_lat", "parma_lon", "parma_alt", "parma_year"):
        assert widget(at, "number_input", key) is not None, key
    assert widget(at, "number_input", "parma_sw") is None       # W from the date


def test_hidden_advanced_values_survive_basic_and_a_restart(new_app, repo_copy):
    at = new_app(mode="Advanced")
    widget(at, "checkbox", "usedetector").check().run()
    widget(at, "number_input", "emin").set_value(5.0).run()
    widget(at, "slider", "thetamax").set_value(70.0).run()
    widget(at, "number_input", "disk_r").set_value(321.0).run()
    widget(at, "number_input", "mg0").set_value(17.0).run()
    widget(at, "toggle", "ui_mode_toggle").set_value(False).run()
    at.run(); at.run(); at.run()                   # Streamlit drops hidden state here
    saved = json.loads((repo_copy / "ucmuon_autosave.json").read_text())
    assert (saved["emin"], saved["thetamax"], saved["disk_r"], saved["mg0"]) == (5.0, 70.0, 321.0, 17.0)
    widget(at, "toggle", "ui_mode_toggle").set_value(True).run()
    assert widget(at, "number_input", "emin").value == 5.0
    assert widget(at, "number_input", "mg0").value == 17.0
    restarted = new_app(keep_autosave=True)
    assert restarted.session_state["ui_mode"] == "Advanced"
    assert widget(restarted, "number_input", "disk_r").value == 321.0


@pytest.mark.skipif(not HAVE_GENERATOR, reason="needs bin/ucmuon_gen_omp (make ucmuon_gen_omp)")
def test_basic_generator_run_uses_derived_values_and_advanced_can_copy_them(new_app):
    import time
    at = new_app()
    widget(at, "number_input", "bd_top").set_value(5.0).run()
    widget(at, "number_input", "bd_height").set_value(1.0).run()
    widget(at, "number_input", "bd_radius").set_value(5.0).run()
    widget(at, "number_input", "nmuonsgen").set_value(500).run()
    [b for b in at.button if "Run UCMuon Surface Generator" in str(b.label)][0].click().run()
    ss = at.session_state
    assert (ss["gen_ui_mode"], ss["gen_emin"], ss["gen_emax"], ss["gen_angular_mode"]) == \
        ("Basic", 1.0, 2500.0, 6)
    assert ss["gen_use_detector"] and ss["gen_detectors"][0]["margin"] > 0
    for _ in range(300):                 # the generator runs in a thread
        time.sleep(1)
        at.run()
        if at.session_state["gen_ntry"]:     # set when the run has finished
            break
    assert at.session_state["gen_ntry"]
    widget(at, "toggle", "ui_mode_toggle").set_value(True).run()
    gen = top_tab(at, "Generator")
    assert any("used **Basic** settings" in i.value for i in gen.info)
    [b for b in gen.button if b.key == "btn_copy_basic"][0].click().run()
    assert widget(at, "selectbox", "angularmode").value == 6
    assert widget(at, "number_input", "emin").value == 1.0


@pytest.mark.skipif(not HAVE_GENERATOR, reason="needs bin/ucmuon_gen_omp (make ucmuon_gen_omp)")
def test_basic_generator_run_without_a_detector(new_app, repo_copy):
    import time
    at = new_app()
    widget(at, "checkbox", "bd_use_detector").uncheck().run()
    for f in ("muons_for_phits.dat", "muons_geant4.txt"):
        (repo_copy / "output" / f).unlink(missing_ok=True)
    widget(at, "checkbox", "savephits").check().run()
    widget(at, "checkbox", "savegeant4").check().run()
    widget(at, "number_input", "bd_src_size").set_value(10.0).run()
    widget(at, "selectbox", "angularmode").set_value(2).run()          # cos²θ
    widget(at, "number_input", "nmuonsgen").set_value(500).run()
    [b for b in at.button if "Run UCMuon Surface Generator" in str(b.label)][0].click().run()
    ss = at.session_state
    assert not ss["gen_use_detector"] and ss["gen_angular_mode"] == 2
    assert ss["gen_basic_used"]["usedetector"] is False
    for _ in range(300):
        time.sleep(1)
        at.run()
        if at.session_state["gen_ntry"]:
            break
    assert at.session_state["gen_ntry"] and not exceptions(at)
    for f in ("muons_for_phits.dat", "muons_geant4.txt"):            # the exports
        assert (repo_copy / "output" / f).stat().st_size > 0, f
    # no detector, so Basic Transport asks for the overburden depth
    assert widget(at, "number_input", "music_depth_m") is not None


# ── Transport ────────────────────────────────────────────────────────────────
def test_basic_transport_shows_the_essentials(new_app, surface_run):
    f, state = surface_run
    at = new_app(state=state)
    assert not exceptions(at)
    for kind, key in (("radio", "density_mode"), ("number_input", "music_rad")):
        assert widget(at, kind, key) is None, key
    for kind, key in (("selectbox", "transport_infile_select"),
                      ("text_input", "transport_infile_custom"),
                      ("text_input", "transport_outfile"), ("checkbox", "transport_all_chk")):
        assert widget(at, kind, key) is not None, key
    assert widget(at, "selectbox", "transport_infile_select").value == str(f)
    widget(at, "selectbox", "transport_engine").set_value("MUSIC").run()
    assert widget(at, "slider", "music_omp_threads") is not None


def test_basic_transport_of_your_own_file(new_app, surface_run, repo_copy):
    """A file that is not the last Generator run's output: the depth is asked
    for and the run's detector does not filter it."""
    import shutil
    f, state = surface_run
    own = repo_copy / "my_muons.dat"
    shutil.copy(f, own)
    det = {"shape": 1, "margin": 50.0, "ax": 0.0, "ay": 0.0, "az": -1100.0,
           "bx": 0.0, "by": 0.0, "bz": -1000.0, "r": 50.0}
    at = new_app(state=dict(state, gen_use_detector=True, gen_detectors=[det]))
    assert widget(at, "number_input", "music_depth_m") is None       # generator's own file
    widget(at, "text_input", "transport_infile_custom").input(str(own)).run()
    assert not exceptions(at)
    assert at.session_state["transport_infile"] == str(own)
    assert widget(at, "number_input", "music_depth_m") is not None
    tr = top_tab(at, "Transport")
    assert any("not from the last Generator run" in c.value for c in tr.caption)


def test_basic_transport_depth_is_the_detectors_top_face(new_app, surface_run):
    det = {"shape": 1, "margin": 50.0, "ax": 0.0, "ay": 0.0, "az": -1100.0,
           "bx": 0.0, "by": 0.0, "bz": -1000.0, "r": 50.0}
    at = new_app(state=dict(surface_run[1], gen_use_detector=True, gen_detectors=[det]))
    assert widget(at, "number_input", "music_depth_m") is None
    assert any("Overburden depth:** 10 m" in m.value for m in top_tab(at, "Transport").markdown)


@pytest.mark.parametrize("engine", ["MUSIC", "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS",
                                    "PROPOSAL", "PUMAS"])
def test_every_engine_renders_in_basic(new_app, engine):
    at = new_app(state={"transport_engine": engine})
    assert not exceptions(at)
    assert widget(at, "selectbox", "phitsxs_mat_choice") is None
    assert widget(at, "selectbox", "proposal_med_choice") is None


def test_basic_has_no_backward_mc_and_keeps_it_for_advanced(new_app, repo_copy):
    bmc, label = "Backward MC Flux Integrator", "⑤ Backward MC flux integrator"
    assert label not in widget(new_app(), "selectbox", "transport_engine").options
    assert label in widget(new_app(mode="Advanced"), "selectbox", "transport_engine").options
    at = new_app(mode="Advanced", state={"transport_engine": bmc})
    widget(at, "toggle", "ui_mode_toggle").set_value(False).run()
    at.run(); at.run(); at.run()
    assert not exceptions(at)
    assert widget(at, "selectbox", "transport_engine") is None
    basic_sel = widget(at, "selectbox", "transport_engine_basic")
    assert basic_sel.value == "UCMuon Stochastic (Python)" and label not in basic_sel.options
    assert json.loads((repo_copy / "ucmuon_autosave.json").read_text())["transport_engine"] == bmc
    widget(at, "toggle", "ui_mode_toggle").set_value(True).run()
    assert widget(at, "selectbox", "transport_engine").value == bmc


def test_transport_advanced_values_survive_basic(new_app):
    at = new_app(mode="Advanced")
    widget(at, "number_input", "music_rad").set_value(30.0).run()
    widget(at, "text_input", "transport_outfile").set_value("output/my_ug.dat").run()
    widget(at, "toggle", "ui_mode_toggle").set_value(False).run()
    at.run(); at.run()
    widget(at, "toggle", "ui_mode_toggle").set_value(True).run()
    assert widget(at, "number_input", "music_rad").value == 30.0
    assert widget(at, "text_input", "transport_outfile").value == "output/my_ug.dat"


# ── Results and Config ───────────────────────────────────────────────────────
def test_basic_results_after_transport(new_app, underground_run):
    ug, state = underground_run
    at = new_app(state=state)
    res = top_tab(at, "Results")
    assert not exceptions(at)
    assert res.selectbox[0].value == str(ug)
    assert any("Underground muons (transport)" in o for o in res.selectbox[0].options)
    labels = [m.label for m in res.metric]
    assert "Rate at depth [/s]" in labels and "Live time" in labels
    rate = next(m for m in res.metric if m.label == "Rate at depth [/s]")
    assert not rate.delta                         # 1.3.1: the error was a green "up" delta
    assert any("(statistical, 1σ)" in c.value for c in res.caption)
    assert [t for t in res.text_input if t.label == "Load any file"]
    assert any("Depth 5 m, ρ = 2.65" in c.value for c in res.caption)
    assert widget(at, "number_input", "res_depth") is None         # the last transport's
    # 3D trajectories, without the overlay editor (it rewrites the run's detector)
    assert [s for s in res.slider if s.label == "Trajectories per class"]
    assert widget(at, "selectbox", "ov_sh") is None
    assert [s for s in res.selectbox if s.label == "Variable"]


def test_basic_results_of_another_file(new_app, underground_run, repo_copy):
    import shutil
    ug, state = underground_run
    other = repo_copy / "other_ug.dat"
    shutil.copy(ug, other)
    at = new_app(state=state)
    res = top_tab(at, "Results")
    [t for t in res.text_input if t.label == "Load any file"][0].input(str(other)).run()
    assert not exceptions(at)
    assert widget(at, "number_input", "res_depth") is not None     # not the last transport


def test_advanced_results_keep_their_controls(new_app, underground_run):
    at = new_app(mode="Advanced", state=underground_run[1])
    res = top_tab(at, "Results")
    assert [t for t in res.text_input if t.label == "Load any file"]
    assert [n for n in res.number_input if n.key == "res_depth"]


def test_config_is_advanced_only(new_app):
    at = new_app()
    assert not at.json and not any("Reset autosave" in str(b.label) for b in at.button)
    adv = new_app(mode="Advanced")
    assert adv.json
    assert any("Reset autosave" in str(b.label) for b in top_tab(adv, "Config").button)


# ── Terrain and Density (Advanced only) ─────────────────────────────────────
def test_basic_has_no_terrain_or_density_controls(new_app):
    at = new_app()
    for kind, key in (("radio", "terrain_dem_mode"), ("number_input", "terrain_naz"),
                      ("radio", "da_method"), ("text_area", "da_lib_paths_text")):
        assert widget(at, kind, key) is None, key


def test_terrain_starts_on_the_bundled_dem(new_app):
    at = new_app(mode="Advanced")
    assert not exceptions(at)
    dem = widget(at, "radio", "terrain_dem_mode")
    assert dem.value == dem.options[0] and dem.options[0].startswith("🌋")
    assert (at.session_state["terrain_lat"], at.session_state["terrain_lon"]) == (40.8271, 14.4006)


def test_terrain_is_one_page(new_app):
    """No Setup / Overburden map / Run & Results sub-tabs (up to v1.2.0):
    sections 1-7 in order on one page."""
    at = new_app(mode="Advanced")
    assert not {"📋  Setup", "🗺️  Overburden map", "▶️  Run & Results"} & {t.label for t in at.tabs}
    heads = [m.value for m in top_tab(at, "Terrain").markdown if m.value.startswith("### ")]
    nums = [h[4] for h in heads if h[4].isdigit()]
    assert nums == sorted(nums) and {"1", "2", "3", "5", "6", "7"} <= set(nums), heads


@pytest.mark.skipif(not HAVE_RASTERIO, reason="needs rasterio")
def test_altitude_from_the_dem(new_app):
    at = new_app(mode="Advanced")
    widget(at, "number_input", "terrain_alt").set_value(100.0).run()
    [b for b in at.button if b.key == "terrain_alt_from_dem"][0].click().run()
    assert 400.0 < at.session_state["terrain_alt"] < 800.0


def test_terrain_values_survive_basic(new_app, repo_copy):
    at = new_app(mode="Advanced")
    widget(at, "number_input", "terrain_naz").set_value(72).run()
    widget(at, "number_input", "terrain_rho").set_value(2.4).run()
    widget(at, "toggle", "ui_mode_toggle").set_value(False).run()
    at.run(); at.run(); at.run()                   # Streamlit drops hidden state here
    saved = json.loads((repo_copy / "ucmuon_autosave.json").read_text())
    assert (saved["terrain_naz"], saved["terrain_rho"]) == (72, 2.4)
    widget(at, "toggle", "ui_mode_toggle").set_value(True).run()
    assert widget(at, "number_input", "terrain_naz").value == 72
    assert widget(at, "number_input", "terrain_rho").value == 2.4


def test_density_with_the_shipped_library(new_app):
    at = new_app(mode="Advanced")
    assert widget(at, "radio", "da_method").value == "library"
    assert "examples/vesuvius/tsim_library" in widget(at, "text_area", "da_lib_paths_text").value
    [b for b in at.button if b.key == "da_load_lib_btn"][0].click().run()
    assert len(at.session_state["da_tsim_lib"]) == 5
    widget(at, "radio", "da_tdata_mode").set_value("🔬 Generate synthetic (test)").run()
    widget(at, "slider", "da_synth_rho").set_value(2.5).run()
    [b for b in at.button if b.key == "da_gen_synth_btn"][0].click().run()
    [b for b in at.button if b.key == "da_run_inversion_btn"][0].click().run()
    assert not exceptions(at)
    med = [m.value for m in at.metric if m.label.startswith("Median ρ̂")]
    assert med and abs(float(med[0].split()[0]) - 2.5) < 0.05


def test_density_and_terrain_inputs_survive_basic(new_app, repo_copy):
    """Up to 1.3.0 the Terrain and Density inputs outside the autosave reset
    after Basic reruns (Basic does not render those tabs)."""
    def basic_cycle(at):
        widget(at, "toggle", "ui_mode_toggle").set_value(False).run()
        at.run(); at.run(); at.run()               # Streamlit drops hidden state here
        widget(at, "toggle", "ui_mode_toggle").set_value(True).run()

    at = new_app(mode="Advanced")
    synth = "🔬 Generate synthetic (test)"
    widget(at, "radio", "da_tdata_mode").set_value(synth).run()
    basic_cycle(at)
    assert widget(at, "radio", "da_tdata_mode").value == synth
    widget(at, "radio", "da_method").set_value("direct").run()
    basic_cycle(at)
    assert widget(at, "radio", "da_method").value == "direct"
    saved = json.loads((repo_copy / "ucmuon_autosave.json").read_text())
    assert (saved["da_method"], saved["da_tdata_mode"]) == ("direct", synth)
