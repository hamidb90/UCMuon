"""Unit tests for gui/basic_mode.py: the values Basic mode derives."""
import math

import pytest

import basic_mode as BM


def test_recommended_ranges_follow_the_validity_table():
    # README "Surface spectra"; decision 2026-10-01.
    assert (BM.recommended(1)["emin"], BM.recommended(1)["emax"]) == (100.0, 2500.0)
    assert BM.recommended(2)["emin"] is None                 # power law: user's range
    assert (BM.recommended(3)["emin"], BM.recommended(3)["emax"]) == (1.0, 1000.0)
    for s in (4, 5):
        assert (BM.recommended(s)["emin"], BM.recommended(s)["emax"]) == (1.0, 2500.0)
    assert BM.recommended(6)["theta_max"] == 60.0            # Gaisser: > 100/cos θ only
    assert (BM.recommended(7)["emin"], BM.recommended(7)["emax"]) == (1.0, 2000.0)
    assert (BM.recommended(8)["emin"], BM.recommended(8)["emax"]) == (0.01, 1.0)
    assert not BM.recommended(2)["normalised"] and not BM.recommended(8)["normalised"]
    assert BM.DEFAULT_SPECTRUM == 4


@pytest.mark.parametrize("shape,kw", [(1, dict(radius_m=0.7)),
                                      (2, dict(width_m=2.0, length_m=3.0))])
def test_detector_round_trip(shape, kw):
    d = BM.detector_from_basic(shape, 12.0, -4.0, 25.0, 1.5, **kw)
    b = BM.basic_from_detector(d)
    assert b["shape"] == shape
    assert (b["cx_m"], b["cy_m"], b["top_m"], b["height_m"]) == pytest.approx((12.0, -4.0, 25.0, 1.5))
    for k, v in kw.items():
        assert b[k] == pytest.approx(v)
    assert BM.top_depth_m(d) == pytest.approx(25.0)
    assert BM.bottom_depth_m(d) == pytest.approx(26.5)


def test_a_face_at_the_ground_is_zero_not_minus_zero():
    """A top face at z = 0 printed as "-0 m" (1.3.1)."""
    d = {"shape": 1, "ax": 0.0, "ay": 0.0, "az": -9000.0,         # 1.3.1's default
         "bx": 0.0, "by": 0.0, "bz": 0.0, "r": 5.0}
    assert math.copysign(1.0, BM.top_depth_m(d)) == 1.0
    assert f"{BM.top_depth_m(d):g}" == "0"


def test_tilted_cylinder_is_not_basic():
    d = BM.detector_from_basic(1, 0.0, 0.0, 10.0, 2.0, radius_m=1.0)
    d["bx"] += 50.0
    assert BM.basic_from_detector(d) is None


def test_margin_is_two_sigma_at_the_deepest_point():
    import mcs_margin
    d = BM.detector_from_basic(1, 0.0, 0.0, 10.0, 1.0, radius_m=0.5)
    ref = mcs_margin.suggested_margin(11.0, 2.65, 1.0, 0.0)["margin_cm"]
    assert BM.margin_cm(d, 2.65, 1.0) == pytest.approx(ref)
    above = BM.detector_from_basic(1, 0.0, 0.0, -5.0, 1.0, radius_m=0.5)   # above ground
    assert BM.margin_cm(above, 2.65, 1.0) == 0.0


@pytest.mark.parametrize("theta", [60.0, 85.0])
@pytest.mark.parametrize("shape", [1, 2])
def test_required_source_reaches_every_path(theta, shape):
    d = BM.detector_from_basic(shape, 30.0, 10.0, 20.0, 2.0, radius_m=1.0, width_m=2.0, length_m=4.0)
    d["margin"] = BM.margin_cm(d, 2.65, 1.0)
    disk = BM.required_source([d], theta, 1)
    assert (disk["disk_cx"], disk["disk_cy"]) == pytest.approx((30.0, 10.0))
    assert disk["disk_r"] % 10.0 == pytest.approx(0.0)       # rounded up to 10 m
    # the GUI's own reach check finds nothing missing for the derived disk ...
    assert BM.source_reach(d, theta, 0.0, 1, 1, 0.0, disk["disk_cx"], disk["disk_cy"],
                           disk["disk_r"], 0, 0, 0, 0) is None
    # ... and 20 m less would be too small
    assert BM.source_reach(d, theta, 0.0, 1, 1, 0.0, disk["disk_cx"], disk["disk_cy"],
                           disk["disk_r"] - 20.0, 0, 0, 0, 0) is not None
    rect = BM.required_source([d], theta, 2)
    assert BM.source_reach(d, theta, 0.0, 2, 1, 0.0, 0, 0, 0, rect["src_u1_m"],
                           rect["src_u2_m"], rect["src_v1_m"], rect["src_v2_m"]) is None


def test_hemisphere_is_centred_on_the_z_axis():
    d = BM.detector_from_basic(1, 30.0, 0.0, 20.0, 2.0, radius_m=1.0)
    d["margin"] = 0.0
    hemi = BM.required_source([d], 85.0, 3)
    assert (hemi["disk_cx"], hemi["disk_cy"]) == (0.0, 0.0)
    # 22 m deep at 85°: 22 tan 85° + 30 m offset + 1 m radius ≈ 283 m -> 290 m
    assert hemi["disk_r"] == pytest.approx(10.0 * math.ceil((22.0 * math.tan(math.radians(85.0)) + 31.0) / 10.0))


def test_remove_candidates_are_the_39_of_the_inventory():
    assert len(BM.REMOVE_CANDIDATES) == 39
