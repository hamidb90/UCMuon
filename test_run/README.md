# UCMuon — comprehensive test run

A small, fully seeded run of the complete two-stage pipeline, with reference
output committed so that a reviewer can confirm the package works after
installing it. It is the test run required by the Computer Physics
Communications program submission.

Everything here runs in a few seconds.

## Running it

From the repository root:

```bash
bash test_run/run_test.sh
```

Exit status 0 means the run reproduced the committed reference.

## What it exercises

| Stage | Program | What it covers |
|---|---|---|
| 1 | `bin/ucmuon_gen_omp` (Fortran 90, OpenMP) | CosmoALEPH spectrum sampling, momentum-dependent charge ratio, `cos^2(theta)` zenith sampling, disk source geometry, 13-column ASCII output |
| 2 | `gui/ucmuon_stochastic_driver.py` (UCMuon-MC, Engine 1) | per-process catastrophic energy loss, delta-ray straggling, Highland multiple scattering, muon decay, the range cut, 18-column ASCII output |

Stage 2 is pure Python and needs no compiler, so it can be checked on its own
even where `gfortran` is unavailable: if `bin/ucmuon_gen_omp` is missing,
`run_test.sh` falls back to the committed surface file and still runs and
grades the transport stage.

## Files

| File | Role |
|---|---|
| `input_gen.dat` | Stage 1 parameters, one documented value per line |
| `input_transport.dat` | Stage 2 parameters, one documented value per line |
| `run_test.sh` | runs both stages and grades the result |
| `expected/muons_surface.dat` | reference Stage 1 output |
| `expected/muons_underground.dat` | reference Stage 2 output |
| `expected/SUMMARY.txt` | reference physics numbers, SHA-256 sums, and the platform they were produced on |
| `output/` | created by the run; not tracked |

## Configuration

2000 muons, CosmoALEPH spectrum over 10–2500 GeV, disk source of radius 50 m,
`cos^2(theta)` up to a 70 degree zenith cut, transported through 25 m of
Standard Rock (2.65 g/cm^3) with UCMuon-MC.

Expected result: **574 of 2000 muons survive (28.70%)**, with a mean exit
kinetic energy of **13.562 GeV**. The depth is chosen so that a useful fraction
both survives and stops, which makes the survival fraction a meaningful check
rather than a formality.

## Reproducibility, and why the comparison has tolerances

`run_test.sh` sets `UCMUON_SEED=20260807` and `OMP_NUM_THREADS=1`, and Stage 2
uses seed 42 with a single worker. Both stages are then deterministic and
independent of the reviewer's core count, and on the reference platform the
output is byte-identical to `expected/`.

The test is nevertheless graded numerically rather than byte for byte: the muon
count must match exactly, the survival fraction to within a quarter of its own
binomial sigma (0.25 pp here), and the mean exit energy to within 0.5%. A
different `gfortran` or C library can change the last bits of a double without
anything being wrong, and a byte comparison would report that as a failure.
The script still reports byte identity when it holds.
