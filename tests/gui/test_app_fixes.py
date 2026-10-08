"""Regression tests for the GUI bugs fixed before the Basic mode (B9-B15, N1
of the inventory of GUI controls; CHANGELOG.md, 1.3.0). Advanced mode, where
the controls concerned live."""
import json

from conftest import exceptions, top_tab, widget

BB = "Bethe-Bloch (PDG) + Groom radiative losses + Highland MS"


def test_b9a_bethe_bloch_uses_the_shared_density(new_app):
    at = new_app(mode="Advanced", state={"transport_engine": BB})
    widget(at, "number_input", "music_rho").set_value(2.71).run()
    assert all("ρ" not in o for o in widget(at, "selectbox", "phitsxs_mat_choice").options)
    assert any("ρ = 2.710 g/cm³ from Medium" in c.value for c in top_tab(at, "Transport").caption)


def test_b9a_driver_density_override(repo_copy, surface_run, tmp_path):
    """"1" and "1 2.65" give the same file; "1 2.71" a different one."""
    import subprocess, sys
    f, _ = surface_run
    outs = {}
    for line in ("1", "1 2.65", "1 2.71"):
        out = tmp_path / f"bb_{line.replace(' ', '_')}.dat"
        stdin = "\n".join([str(f), str(out), "1", "14", "5", line, "0"]) + "\n"
        subprocess.run([sys.executable, str(repo_copy / "gui" / "ucmuon_bb_driver.py")],
                       input=stdin, text=True, capture_output=True, check=True)
        outs[line] = out.read_text()
    assert outs["1"] == outs["1 2.65"]
    assert outs["1"] != outs["1 2.71"]


def test_b9b_old_proposal_label_is_mapped(new_app):
    at = new_app(mode="Advanced", state={"transport_engine": "PROPOSAL",
                                         "proposal_med_choice": "Standard Rock  (ρ=2.65 g/cm³)"})
    assert not exceptions(at)
    assert widget(at, "selectbox", "proposal_med_choice").value == "Standard Rock"


def test_b9c_sigma_is_said_not_to_be_propagated(new_app):
    at = new_app(mode="Advanced")
    widget(at, "radio", "density_mode").set_value("Gaussian prior").run()
    assert any("mean ρ only" in i.value for i in at.info)


def test_b14_threads_slider_and_worker_count(new_app):
    at = new_app(mode="Advanced")
    assert widget(at, "slider", "music_omp_threads") is None          # UCMuon-MC
    mode = [m.value for m in top_tab(at, "Transport").metric if m.label == "Mode"]
    assert mode and "Single-thread" not in mode[0]
    widget(at, "selectbox", "transport_engine").select("MUSIC").run()
    assert widget(at, "slider", "music_omp_threads") is not None


def test_b15_parma_w_from_the_date(new_app):
    at = new_app(mode="Advanced")
    widget(at, "selectbox", "gen_spectrum_mode").select(3).run()
    cb = widget(at, "checkbox", "parma_sw_from_date")
    assert cb is not None and cb.value is True
    assert widget(at, "number_input", "parma_sw") is None
    cb.uncheck().run()
    assert widget(at, "number_input", "parma_sw") is not None


def test_b10a_terrain_file_selector_is_keyed(new_app, surface_run):
    at = new_app(mode="Advanced", state=surface_run[1])
    assert widget(at, "selectbox", "terrain_infile_select") is not None


def test_b11_n1_results_after_transport(new_app, underground_run):
    ug, state = underground_run
    at = new_app(mode="Advanced", state=state)
    res = top_tab(at, "Results")
    assert res.selectbox[0].value == str(ug)                 # B11c
    assert "Rate at depth [/s]" in [m.label for m in res.metric]   # B11b


def test_b12_config_restore_applies_settings(new_app, repo_copy):
    """Restore applies a settings dict (the config JSON's "settings" part, or
    an autosave) before the widgets render."""
    at = new_app(mode="Advanced", state={"_cfg_restore_pending": {"emin": 7.0, "thetamax": 66.0}})
    assert widget(at, "number_input", "emin").value == 7.0
    assert widget(at, "slider", "thetamax").value == 66.0


def test_b13a_no_plain_gaisser_for_the_direct_inversion(new_app):
    at = new_app(mode="Advanced", state={"da_method": "direct"})
    assert "bugaev" not in widget(at, "selectbox", "di_model_sel").options


def test_b13b_measured_data_without_density_header(repo_copy, tmp_path):
    import sys
    sys.path.insert(0, str(repo_copy / "gui"))
    import ucmuon_density_analysis as da
    f = tmp_path / "t.dat"
    f.write_text("# az el T\n0 10 0.5\n0 20 0.6\n90 10 0.4\n90 20 0.7\n")
    *_, meta = da.load_transmission_map(f, require_density=False)
    assert meta["density"] is None


def test_underground_filter_uses_the_real_detector(new_app, repo_copy):
    """The detector-hit file after transport, and so the rate at depth Results
    gives it, is that of the real detector, not of the generator's
    margin-inflated one (v1.3.1). A 10 cm cylinder with a 2 m margin; one
    muon alive inside it, one alive 1 m off its axis, one stopped 1 m off it:
    only the first is a hit (the 1.3.0 filter counted all three)."""
    import numpy as np
    out = repo_copy / "output"
    ug, sel = out / "ug_margin_test.dat", out / "ug_margin_test_selected.dat"
    det = {"shape": 1, "margin": 200.0, "ax": 0.0, "ay": 0.0, "az": -600.0,
           "bx": 0.0, "by": 0.0, "bz": -550.0, "r": 10.0}
    rows = [  # EventID xs ys zs Es theta_s phi_s charge alive x y z E cx cy cz theta phi
        [1, 0, 0, 0, 20.0, 0, 0, 1, 1, 0.0, 0, -500, 18.0, 0, 0, -1, 0, 0],
        [2, 100, 0, 0, 20.0, 0, 0, 1, 1, 100.0, 0, -500, 18.0, 0, 0, -1, 0, 0],
        [3, 100, 0, 0, 2.0, 0, 0, -1, 0, 100.0, 0, -575, 0.0, 0, 0, 0, 0, 0],
    ]
    np.savetxt(ug, np.array(rows, float), fmt="%g",
               header="EventID xs ys zs Es theta_s phi_s charge alive x y z E cx cy cz theta phi")
    try:
        run_state = {f"{w}_{k}": v for w in ("gen", "music") for k, v in
                     (("lines", []), ("running", False), ("success", None),
                      ("stop_req", False), ("nmuons", 0), ("proc", None),
                      ("start_time", None), ("end_time", None))}
        run_state["music_success"] = True          # a finished transport
        at = new_app(mode="Advanced", state={
            "_state": run_state, "ug_use_filter": True, "ug_filter_done": False,
            "gen_use_detector": True, "gen_detectors": [det], "ug_file": str(ug),
            "ug_filter_file": str(sel), "ug_depth_m": 5.0})
        assert not exceptions(at)
        assert at.session_state["ug_filter_done"] is True
        hits = np.loadtxt(sel, comments="#", ndmin=2)
        assert [int(h[0]) for h in hits] == [1]
    finally:
        for f in (ug, sel):
            f.unlink(missing_ok=True)


def test_3d_title_names_the_file_kind(new_app, underground_run):
    """The 3D viewer said "UG Selected" for every file (v1.3.1)."""
    ug, state = underground_run
    at = new_app(mode="Advanced", state=state)
    specs = [c.proto.spec for c in at.get("plotly_chart")]
    assert any("3D Muon Trajectories: underground" in s for s in specs)
    assert not any("UG Selected" in s for s in specs)
