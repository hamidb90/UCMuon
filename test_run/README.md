# UCMuon — comprehensive test run

A small, seeded run of the complete two-stage pipeline, graded statistically
against a reference, so that a reviewer can confirm the package works after
installing it. It is the test run required by the Computer Physics
Communications program submission. It takes about 15 seconds.

## Running it

From the repository root:

```bash
bash test_run/run_test.sh
```

Exit status 0 means every statistic agrees with the reference (below).

## What it exercises

| Stage | Program | What it covers |
|---|---|---|
| 1 | `bin/ucmuon_gen_omp` (Fortran 90, OpenMP) | Guan et al. (2015) spectrum, joint sampling of momentum and direction through the source surface (angular mode 6), momentum-dependent charge ratio, disk source, the detector acceptance cut (only muons whose line passes within 1 m of a cylinder 25 m below are kept) and the live time those muons represent, 13-column ASCII output of the selected muons |
| 2 | `gui/ucmuon_stochastic_driver.py` (UCMuon-MC, Engine 1) | per-process catastrophic energy loss, delta-ray straggling, Highland multiple scattering, muon decay, the range cut, 18-column ASCII output |

Configuration: 10⁴ selected muons (about 3.8 million tried) of 5–2500 GeV, disk
source of radius 50 m, zenith angles up to 70°, transported through 25 m of
Standard Rock (2.65 g/cm³), the depth of the detector's top face.

Stage 2 is pure Python and needs no compiler: if `bin/ucmuon_gen_omp` is
missing, `run_test.sh` takes the reference platform's selected muons
(`expected/muons_selected.dat.gz`) and still runs and grades Stage 2.

## How it is graded

`grade.py` computes six statistics with their statistical errors from the
run's own files: the mean log₁₀ E, the μ⁺/μ⁻ ratio, the mean cos θ and the live
time per selected muon of Stage 1, and the survival fraction and mean exit
kinetic energy of Stage 2. Each must lie within 4 standard deviations of the
reference, the deviation combining the run's error with the reference's. The
reference (`expected/reference_stats.csv`) is a run 20 times larger with other
seeds, made by `make_reference.sh`.

| Statistic | Reference (2×10⁵ muons) | This run's error |
|---|---|---|
| mean log₁₀(E/GeV) | 1.0614 | 0.003 |
| μ⁺/μ⁻ | 1.253 | 0.025 |
| mean cos θ | 0.7602 | 0.0015 |
| live time per selected muon | 1.224 ms | 0.012 ms |
| survival fraction | 20.52 % | 0.40 % |
| mean exit kinetic energy | 20.39 GeV | 0.8 GeV |

A correct installation therefore fails a given check once in about 16 000
runs, and one of the six once in about 2600, whatever its compiler, C library or
random stream: five runs with other seeds (generator 1001–1005, transport
43–47) all pass, the largest deviation being 2.0σ. A real change is caught:
rock 10 % denser fails the survival check by 8.5σ.

`run_test.sh` sets `UCMUON_SEED=20260807` and `OMP_NUM_THREADS=1`, and Stage 2
uses seed 42 with one worker, so on the reference platform the output is also
byte-identical to the reference run's (`expected/sha256.txt`). This is
reported as a regression check but not required. Up to UCMuon 1.4.1 the test
used 2000 muons, the legacy cos²θ mode without a detector, and tolerances a
quarter of the statistical error, so it passed only on an identical random
stream.

## Files

| File | Role |
|---|---|
| `input_gen.dat` | Stage 1 parameters, one documented value per line (all parameters: `docs/GENERATOR_INPUT.md`) |
| `input_transport.dat` | Stage 2 parameters, one documented value per line |
| `run_test.sh` | runs both stages, reports byte identity and grades the result |
| `grade.py` | the six statistics and the 4σ comparison |
| `make_reference.sh` | regenerates `expected/` (a few minutes) |
| `expected/reference_stats.csv` | the reference statistics |
| `expected/muons_selected.dat.gz` | the reference platform's Stage 1 output, for Stage 2 without a compiler |
| `expected/sha256.txt` | SHA-256 of the reference platform's output |
| `expected/SUMMARY.txt` | the reference platform and its results |
| `output/` | created by the run; not tracked |
