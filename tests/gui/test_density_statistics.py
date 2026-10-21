"""Regression tests for the density inversion's statistics (mock-review item 6;
each fails on 1.4.1)."""
import numpy as np
import pytest


def _exp_library(rhos=(1.5, 2.0, 2.65, 3.0)):
    """T = exp(-rho / lam), ln T exactly linear in rho: the inversion is exact
    without noise, so the pulls measure sigma_rho alone."""
    lam = np.linspace(0.25, 0.9, 40).reshape(5, 8)
    return {r: np.exp(-r / lam) for r in rhos}, lam


def test_pulls_have_unit_rms():
    """sigma_rho came from a binomial sigma_T and np.gradient of T at the
    library densities: pulls with an RMS of 1.3 on the Vesuvius example. With
    the Poisson sigma_T of the data and the slope of the ln-T interpolation
    that gives rho_hat, the RMS is 1."""
    import ucmuon_density_analysis as da
    lib, _ = _exp_library()
    pulls = []
    for seed in range(40):
        for true in (1.75, 2.0, 2.3):
            T, s = da.generate_synthetic_tdata(lib, true, n_events=100_000, seed=seed)
            rho, sr, st = da.invert_density_map(T, lib, s)
            ok = st == 0
            pulls.append((rho[ok] - true) / sr[ok])
    p = np.concatenate(pulls)
    assert p.size > 3000
    assert abs(np.sqrt(np.mean(p ** 2)) - 1.0) < 0.05, np.sqrt(np.mean(p ** 2))


def test_synthetic_sigma_is_poisson():
    """N_rock ~ Poisson(N T), so sigma_T = sqrt(T / N); it was sqrt(T(1-T)/N)."""
    import ucmuon_density_analysis as da
    lib, _ = _exp_library()
    T, s = da.generate_synthetic_tdata(lib, 2.0, n_events=10_000, seed=1)
    assert np.allclose(s, np.sqrt(lib[2.0] / 10_000))


def test_low_sensitivity_flag_uses_sigma_rho():
    """Status 4 was |dT/drho| < 0.005 whatever the statistics: under thick
    rock T and its slope are small, but so is sigma_T. With sigma_T it is
    sigma_rho > max_sigma_rho; without, the slope as before."""
    import ucmuon_density_analysis as da
    lib, lam = _exp_library()
    T = lib[2.0]
    slope = T / lam                                     # |dT/drho| at rho = 2
    deep = slope < 0.005
    assert deep.any() and (~deep).any()
    s = np.sqrt(T / 1e6)
    _, sr, st_sigma = da.invert_density_map(T, lib, s, max_sigma_rho=0.1)
    _, _, st_slope = da.invert_density_map(T, lib)
    assert np.all(st_slope[deep] == 4) and np.all(st_slope[~deep] == 0)
    assert np.all((st_sigma == 4) == (sr > 0.1))
    assert np.any((st_sigma == 0) & deep)               # measurable after all


def test_no_library_value_is_not_inverted_and_a_single_density_is_refused():
    import ucmuon_density_analysis as da
    lib, _ = _exp_library()
    lib[3.0][0, 0] = np.nan
    T, s = da.generate_synthetic_tdata(lib, 2.0, n_events=100_000, seed=2)
    _, _, st = da.invert_density_map(T, lib, s)
    assert st[0, 0] == 1 and np.isnan(T[0, 0])
    with pytest.raises(ValueError):
        da.invert_density_map(T, {2.0: lib[2.0]}, s)


def test_transmission_beyond_the_validated_opacity_is_nan():
    """The terrain driver wrote T for every direction; beyond 5000 m w.e. (the
    end of Engine 5's calibration) the flux is not validated, and under the
    thickest rock it was 0. Those directions are NaN now."""
    import ucmuon_backward_mc as bm
    import ucmuon_terrain_driver as td
    assert bm.opacity_limit(2, 3) == 5.0e5
    assert bm.opacity_limit(1, 3) == bm.RANGE_TABLE_R_MAX_GCM2
    assert bm.opacity_limit(2, 3, mat_id=2) == bm.RANGE_TABLE_R_MAX_GCM2
    flux = np.array([[1.0, 0.5, 1e-6, 0.0]])
    ob = np.array([[0.0, 1e5, 5.5e5, 7e5]])
    T = td.compute_transmission_map(flux, np.ones(4), ob, bm.opacity_limit(2, 3))
    assert T[0, 0] == 1.0 and T[0, 1] == 0.5 and np.isnan(T[0, 2]) and np.isnan(T[0, 3])
    assert td.compute_transmission_map(flux, np.ones(4))[0, 3] == 0.0   # unchanged without


def test_transmission_file_keeps_small_values(tmp_path):
    """Six decimals rounded T = 1e-6 by up to 40 %."""
    import ucmuon_density_analysis as da
    import ucmuon_terrain_driver as td
    T = np.array([[1.234567e-6, 0.5, np.nan]])
    f = tmp_path / "t.dat"
    td.write_transmission_map(np.array([0.5]), np.array([80.0, 70.0, 60.0]), T, f,
                              40.0, 14.0, 600.0, 2.0)
    _, _, back, _ = da.load_transmission_map(str(f))
    assert abs(back[0, 0] / 1.234567e-6 - 1) < 1e-6 and back[0, 1] == 0.5 and np.isnan(back[0, 2])


def test_direct_inversion_defaults_to_guan_and_warns_for_plain_gaisser():
    """Method 2 defaulted to Reyna, which reads the Vesuvius closure-test data
    as 1.70 g/cm3; the libraries use Guan. The plain Gaisser formula (0.15
    g/cm3 there) is not offered by the GUI and now warns in the API."""
    import inspect
    import warnings
    import ucmuon_density_analysis as da
    for f in (da.transmission_vs_opacity, da.invert_opacity_map, da.path_length_from_tsim):
        assert inspect.signature(f).parameters["model"].default == "guan_2015"
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        da.invert_opacity_map(np.array([[0.1]]), np.array([10.0]), model="bugaev")
        n_bugaev = len(rec)
        da.invert_opacity_map(np.array([[0.1]]), np.array([10.0]))
    assert n_bugaev == 1 and "Gaisser" in str(rec[0].message) and len(rec) == 1


def test_direct_inversion_matches_engine5_csda():
    """Method 2 integrated the open sky from 0.5 GeV and the rock from zero
    energy at the detector: -17 to +13 % from Engine 5's CSDA integral looking
    vertically. With one threshold (1 GeV) on both they are the same physics."""
    import ucmuon_backward_mc as bm
    import ucmuon_density_analysis as da
    xs = np.array([0.0, 1e3, 1e4, 5e4, 1e5, 2e5])
    for th in (0.0, 60.0, 80.0):
        _, T = da.transmission_vs_opacity(th, "guan_2015", 0.0, xs)
        kw = dict(E_min_GeV=1.0, E_max_GeV=2500.0, n_E=80, mode=0)
        sky = bm.directional_flux(0.0, np.radians(th), 3, **kw)
        e5 = np.array([bm.directional_flux(x, np.radians(th), 3, **kw) for x in xs]) / sky
        assert T[0] == 1.0
        assert np.allclose(T[1:] / e5[1:], 1.0, atol=0.006), (th, T[1:] / e5[1:])
