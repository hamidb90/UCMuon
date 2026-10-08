"""Regression tests for the bugs found by the pre-release audit of 1.3.1
(CHANGELOG.md, 1.3.1). Pure functions, no AppTest."""
import math
import subprocess
import sys

import numpy as np
import pytest


def test_fast_estimator_range_table_is_groom():
    """Groom (2001) Table IV-6, standard rock: the table was 1-11 % long from
    140 GeV to 1 TeV and ended at 1 TeV (deeper opacities gave no flux)."""
    import fast_flux_estimator as f
    for t_gev, r_gcm2 in ((100.0, 4.084e4), (200.0, 7.459e4), (1000.0, 2.453e5),
                          (2000.0, 3.616e5)):
        assert f._R_of_T(t_gev) == pytest.approx(r_gcm2, rel=1e-3)
    assert f.emin_from_opacity(2.453e5) == pytest.approx(1000.0, rel=1e-3)
    e = f.emin_from_opacity(5.0e5)                     # beyond the table
    assert e is not None and 2000.0 < e < 1.0e4
    assert f._R_of_T(f._T_of_R(5.0e5)) == pytest.approx(5.0e5, rel=1e-6)


def test_pumas_spectrum_integrates_to_the_rate(tmp_path):
    """Each backward event carries Φ/pdf, so dΦ/dE in a bin is Σfv/(NΔE); the
    bin mean gave Φ·E·ln(Emax/Emin). The binned spectrum must integrate to
    the event-level rate."""
    import ucmuon_pumas_driver as pd_
    rng = np.random.default_rng(3)
    n = 20000
    e = np.exp(rng.uniform(math.log(1.0), math.log(1000.0), n))       # log-uniform
    fv = 0.5 * e ** -1.7 * e * math.log(1000.0)                        # Φ(E)/pdf(E)
    rows = np.column_stack([np.arange(n), e, rng.uniform(0.5, 1.0, n),
                            rng.choice([-1, 1], n), e * 2.0, fv])
    f = tmp_path / "ev.dat"
    np.savetxt(f, rows, header="ev E_det cos charge E_surf flux_contribution")
    spec = pd_.compute_flux_spectrum(str(f), n_bins=40, cos_theta_min=1.0, cos_theta_max=0.5)
    w = 2.0 * math.pi * 0.5
    edges = np.exp(np.linspace(np.log(e.min()), np.log(e.max()), 41))
    integral = w * float(np.sum(spec["flux"] * np.diff(edges)))
    assert integral == pytest.approx(spec["rate_m2s"], rel=1e-9)
    assert int(spec["counts"].sum()) == n                             # last event kept
    # and dΦ/dE follows the input Φ ∝ E^-1.7 (not E^-0.7)
    ok = spec["counts"] > 200
    slope = np.polyfit(np.log(spec["E_det_GeV"][ok]), np.log(spec["flux"][ok]), 1)[0]
    assert slope == pytest.approx(-1.7, abs=0.1)


def test_bethe_bloch_driver_directions_are_finite(repo_copy, surface_run, tmp_path):
    """The Python Bethe-Bloch took the kinetic energy for the total in the
    Highland momentum: NaN directions below 105.7 MeV, too much scattering
    above. Survival and exit energies do not depend on it."""
    f, _ = surface_run
    ug = tmp_path / "ug.dat"
    stdin = "\n".join([str(f), str(ug), "1", "14", "8", "1 2.65", "1"]) + "\n"
    subprocess.run([sys.executable, str(repo_copy / "gui" / "ucmuon_bb_driver.py")],
                   input=stdin, text=True, capture_output=True, check=True)
    d = np.loadtxt(ug, comments="#")
    alive = d[:, 8] == 1
    assert alive.sum() > 100
    assert np.isfinite(d[alive, 13:18]).all()


def test_terrain_bins_by_arrival_direction():
    """The ray tracer's azimuth is where the muon comes from; the file holds
    the direction of travel. A muon travelling west came from the east."""
    import pandas as pd
    import gui_terrain_engine as te
    th = math.radians(40.0)
    df = pd.DataFrame({"cx": [-math.sin(th), 0.0], "cy": [0.0, -math.sin(th)],
                       "cz": [-math.cos(th)] * 2})          # travelling W, travelling S
    az, ze = te.assign_direction_bins(df, 36, 18, 90.0)
    assert list(az) == [9, 0]                               # from E (90°), from N (0°)


def test_density_library_nan_pixels_are_no_data():
    """Terrain-saved T_sim files leave bins without muons NaN: such pixels
    were "OK" with ρ̂ = NaN, and the synthetic generator crashed on them."""
    import ucmuon_density_analysis as da
    lib = {r: np.full((2, 2), t) for r, t in ((1.5, 0.30), (2.0, 0.20), (2.65, 0.10), (3.0, 0.05))}
    for v in lib.values():
        v[0, 0] = np.nan
    t, s = da.generate_synthetic_tdata(lib, true_rho=2.0, n_events=10000)
    rho, _, status = da.invert_density_map(t, lib, s)
    assert status[0, 0] == 1 and np.isnan(rho[0, 0])
    assert (status[1:, :] == 0).all() and np.isfinite(rho[1, 1])


@pytest.mark.parametrize("engine", ["ucmuon_bb_driver.py", "ucmuon_stochastic_driver.py"])
def test_stopped_muons_carry_their_stopping_point(repo_copy, surface_run, tmp_path, engine):
    """Every engine writes a stopped muon's (x, y, z) where it stopped (up to
    1.3.0 the Python engines wrote the surface x, y). Without scattering the
    stopping point lies on the surface direction: Δx = tanθ cosφ |Δz|."""
    f, _ = surface_run
    ug = tmp_path / "ug.dat"
    if engine == "ucmuon_bb_driver.py":
        stdin = "\n".join([str(f), str(ug), "1", "14", "8", "1 2.65", "0"]) + "\n"
    else:
        sys.path.insert(0, str(repo_copy / "gui"))
        import gui_stochastic_engine as ge
        stdin = ge.build_stochastic_stdin(dict(
            infile=str(f), outfile=str(ug), rho=2.65, rad=26.48, depth_m=8.0,
            stochastic_ms_enable=0, transport_all=True, ncols=14, stochastic_n_workers=1))
    subprocess.run([sys.executable, str(repo_copy / "gui" / engine)],
                   input=stdin, text=True, capture_output=True, check=True)
    d = np.loadtxt(ug, comments="#")
    st_ = (d[:, 8] == 0) & (np.abs(d[:, 11]) > 10.0) & (d[:, 5] > 0.3)
    assert st_.sum() > 20
    th, ph, dz = d[st_, 5], d[st_, 6], np.abs(d[st_, 11] - d[st_, 3])
    assert np.allclose(d[st_, 9] - d[st_, 1], np.tan(th) * np.cos(ph) * dz, rtol=0.02, atol=0.5)
    assert np.allclose(d[st_, 10] - d[st_, 2], np.tan(th) * np.sin(ph) * dz, rtol=0.02, atol=0.5)
