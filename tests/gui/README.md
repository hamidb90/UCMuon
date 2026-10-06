# GUI tests

Headless tests of the Streamlit GUI with Streamlit's AppTest, run by
`.github/workflows/gui.yml`.

| file | what it checks |
|---|---|
| `test_basic_mode.py` | `gui/basic_mode.py`: recommended range per spectrum, detector re-parametrisation, the 2σ_r margin at the deepest point, the derived source reaching every path at θ ≤ θ_max |
| `test_app_basic.py` | Basic and Advanced mode, tab by tab: Basic's three tabs, what Basic shows and hides, values it derives, hidden Advanced settings (also of the Advanced-only tabs) surviving a mode switch and a restart, a Basic generator run (needs `bin/ucmuon_gen_omp`) and what Advanced then shows, Terrain on the bundled DEM, a density inversion with the shipped library |
| `test_app_fixes.py` | the GUI bugs B9-B15 / N1 found by the inventory of GUI controls (CHANGELOG.md, 1.3.0) |

```
python3 -m pip install pytest
python3 -m pytest tests/gui -v
```

The GUI changes to the project root of the copy it runs from and writes
`ucmuon_autosave.json` there, so `conftest.py` copies the repository to a
temporary directory once per session and every test runs the GUI from that
copy: running the tests never touches your own saved settings. Tests that need
rasterio or the Fortran generator are skipped without them.
