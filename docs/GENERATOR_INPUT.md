# Surface generator input

`bin/ucmuon_gen_omp` (OpenMP) and `bin/ucmuon_gen` (MPI + OpenMP) read their
parameters from standard input, one value per line (the detector's axis points
and box ranges take three or two values on one line). Which lines follow depends
on earlier values, as listed below. The GUI writes this stream for you
(`build_ucmuon_input` in `gui/ucmuon_gui.py`); `test_run/input_gen.dat` and
`hpc/input_params.dat` are commented examples for the two programs.

Every value is checked as it is read. A missing or unreadable line, or a value
out of range, stops the program with a message naming the parameter
(`INPUT ERROR: angular_mode must be 1 to 6, not 9`) and exit status 2; the
MPI program stops all its ranks. Up to UCMuon 1.4.1 most bad values were
replaced silently.

## Order

| # | Parameter | Unit | Valid values | Read when |
|---|---|---|---|---|
| 1 | `use_defaults` | | 0 = read the parameters below; 1 = built-in defaults (100–2500 GeV, CosmoALEPH, disk R = 800 m, cos²θ to 85°, 10⁵ muons, no detector), and nothing more is read | always |
| 2 | `e_min` | GeV | > 0 (muon spectra: above the muon mass) | `use_defaults` = 0 |
| 3 | `e_max` | GeV | > `e_min` | |
| 4 | `spectrum_mode` | | 1 CosmoALEPH, 2 power law E⁻³·⁷, 3 PARMA/EXPACS, 4 Guan 2015, 5 Frosin 2025, 6 Gaisser 1990, 7 Reyna 2006, 8 cosmic e± | |
| 4a | `parma_lat`, `parma_lon` | deg | −90…90, −360…360 | `spectrum_mode` = 3 |
| 4b | `parma_alt_km` | km | altitude of the site | |
| 4c | `parma_year`, `parma_month`, `parma_day` | | month 1…12, day 1…31 | |
| 4d | `parma_charge_mode` | | 0 both, 1 μ⁺ only, −1 μ⁻ only | |
| 4e | `parma_datapath` | | directory of the PARMA data (empty line = `.`) | |
| 4f | `parma_s_W` | | solar-modulation W index; ≤ −999 takes it from the date | |
| 5 | `source_mode` | | 1 disk, 2 rectangle, 3 hemisphere | |
| 6 | `source_plane` | | 1 XY (muons travel in −z), 2 XZ (−y), 3 YZ (−x) | `source_mode` = 1 or 2 |
| 7 | disk: `src_disk_cx_m`, `src_disk_cy_m`, `src_disk_r_m`, `src_w_m`, `src_tilt_deg`, `src_tilt_az_deg` | m, deg | radius > 0; tilt 0 ≤ t < 90 | `source_mode` = 1 |
| 7 | rectangle: `src_u1_m`, `src_u2_m`, `src_v1_m`, `src_v2_m`, `src_w_m`, `src_tilt_deg`, `src_tilt_az_deg` | m, deg | u2 > u1, v2 > v1; tilt 0 ≤ t < 90 | `source_mode` = 2 |
| 7 | hemisphere: `hemi_radius_m`, `hemi_cz_m` | m | radius > 0 | `source_mode` = 3 |
| 8 | `angular_mode` | | 1 vertical, 2 cos²θ, 3 uniform cone, 4 Guan/Frosin P(θ\|E), 5 cos³θ, 6 joint J(p, θ) through the surface (recommended: the only mode whose muons are flux-distributed, so that the live time needs no weights) | |
| 9 | `theta_max` | deg | 0 < θ ≤ 90 | `angular_mode` ≥ 2 |
| 10 | `nmuons` | | ≥ 1: muons to write, or with a detector the muons to select | |
| 11 | `use_detector` | | 0 no cut; 1 keep only muons whose straight line passes a detector (detector acceptance cut) | |
| 12 | `ndet` | | 1…10 | `use_detector` = 1 |
| 13 | per detector: `shape`, `safety margin` | —, cm | 1 cylinder, 2 box; margin ≥ 0 | |
| 13a | cylinder: `Ax Ay Az` (one line), `Bx By Bz` (one line), radius | cm | radius > 0; A ≠ B | `shape` = 1 |
| 13b | box: `xmin xmax`, `ymin ymax`, `zmin zmax` (one line each) | cm | min < max | `shape` = 2 |
| 14 | `save_all` | | 1 write every generated muon; 0 only the selected ones (needs a detector) | |
| 15 | `save_phits` | | 0 / 1 (OpenMP program only) | |
| 16 | surface output file | | empty line = `muons_surface.dat` | |
| 17 | selected-muon file | | empty line = `muons_selected.dat` | `use_detector` = 1 |
| 18 | PHITS file | | empty line = `muons_for_phits.dat` | `save_phits` = 1 (OpenMP) |
| 19 | (empty "Press Enter" line) | | may be omitted at the end of a piped input | OpenMP program |

## The MPI program

`bin/ucmuon_gen` reads the same stream on rank 0 and broadcasts it, except
that it has no `save_phits` line (15), no PHITS file (18) and no final empty
line, and that each rank writes its own files, named after the given ones
(default `ucmuon_surface.dat`, `ucmuon_selected.dat`) with `_RRRRR.dat` in place
of a final `.dat`; `hpc/run_ucmuon_gen.sh` concatenates them.

## Seeds and threads

The random seed comes from the environment variable `UCMUON_SEED` (otherwise
from the clock and the process ID); a run is reproduced event by event with
the same seed and `OMP_NUM_THREADS=1`.
