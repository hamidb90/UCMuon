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
