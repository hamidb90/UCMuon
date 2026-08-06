# UCMuGen validation

The Fortran generator in `src/generator` is the reference implementation. This
suite exists to prove that `UCMuGen.h` reproduces it, and to state precisely
how large an error it would have caught if it did not.

## Why this is not a normal test suite

The obvious approach, "run a KS test and require p > 0.05", does not work for a
random generator:

- Two samples from *identical* distributions give p-values uniform on [0, 1].
  A `p > 0.05` gate therefore fails 5% of the time by construction, and across
  a 23-configuration matrix something fails on essentially every run.
- At the sample sizes needed here (1e5 to 1e6), the test also becomes sensitive
  to differences far below anything physically meaningful, so a correct port
  fails on a 1e-4 relative shift.

Hand-picked tolerances are no better. The right bound depends on the statistic,
the sample size, and the shape of the distribution. The standard deviation of a
steeply falling momentum spectrum is dominated by its rare high-momentum tail
and fluctuates by more than 1% between runs of the *same* model at N = 5e4, so
a 0.5% tolerance on it fails 100% of the time on correct code. (For
dN/dp ~ p^-2.7 the second moment diverges outright, so the sample mean of p has
no 1/sqrt(N) scaling at all. Median and IQR carry the statistical power here;
`mean:p` is retained but is deliberately weak.)

## What it does instead

Every threshold is **measured, not chosen**:

1. `stats.calibrate` runs the same model `n_runs` times with different seeds and
   compares every distinct pair, giving `n_runs*(n_runs-1)/2` null observations
   from `n_runs` generator runs.
2. Each probe's threshold is `3 x median` of what it observed under the null.
3. A configuration passes when every probe falls inside its measured band.

The pass criterion then has a concrete meaning: *the port differs from the
Fortran by no more than two Fortran runs differ from each other.*

The safety factor of 3.0 is itself calibrated, against the two-sample
Kolmogorov distribution: over 4000 replicates at n = m = 2e4,
`P(D > 3 x median) < 1/4000`, while `P(D > 2 x median) = 7e-3` (too tight for a
23-configuration matrix). The median is used rather than a high quantile
because a 99th percentile estimated from a handful of trials is essentially the
maximum of those trials: measured at N = 5e4 and 2e5, such thresholds failed to
scale as 1/sqrt(N), varying between 1.2x and 3.3x instead of the expected 2.0x.

## Validating the validator

Null calibration can only ever make a suite too loose, never too tight, so it is
worthless on its own. `stats.resolving_power` is the counterpart: it bisects the
detection boundary for injected distortions standing in for the porting mistakes
that actually happen (a units or mass slip, a mis-transcribed spectral index, a
wrong angular exponent, a broken charge-ratio parametrisation).

Reproduce with `python3 measure_power.py`:

| distortion       |  N = 50k |  N = 200k |  N = 800k |
|------------------|---------:|----------:|----------:|
| momentum scale   |   0.67%  |    0.53%  |    0.28%  |
| spectral index   |   0.71%  |    0.70%  |    0.31%  |
| angular exponent |   1.00%  |    0.60%  |    0.29%  |
| charge fraction  |   0.94%  |    0.45%  |    0.23%  |

Read as: the smallest error of that kind the suite would have caught. The
1/sqrt(N) scaling is visible and is the reason the reference matrix runs at
N = 1e6, where the guarantee is roughly 0.25% across the board.

Two details matter for this number to be honest, and both were wrong in the
first version:

- The distortion must be applied to a sample **independent** of the reference.
  Perturbing a sample against itself shares its fluctuations, which cancel and
  make the suite look more sensitive than it is.
- Every distortion must be exactly the identity at zero magnitude. The spectral
  index shift therefore uses rejection sampling; resampling with replacement is
  a bootstrap and perturbs the sample even at shift = 0, which fakes an
  arbitrarily good resolving power under bisection.

## Numbers in the documentation

A number lives in three places: the program that computes it, the README that
quotes it, and the paper that quotes it. Only the first executes. The other two
are transcriptions, and transcriptions rot as the code changes underneath them.

Both known instances were found by hand rather than by a test, which is the
reason this check exists:

- the Geant4 example's README advertised `99928` muons and `325.283 Hz`, from
  before that example became multithreaded;
- `feature_tour.cc` printed a ratio of *acceptances* under a column headed
  `speed-up`, with a code comment asserting the two were interchangeable. They
  differ by a factor of two to three, because a directed proposal costs more
  than a blind one.

Neither changed any physics. Both would have been quoted by a reader as fact.

```
feature_tour --numbers | check_numbers.py --exe <feature_tour>
```

`--numbers` emits every documented value as `key value`; `check_numbers.py`
diffs that against `reference/tour_numbers.txt` and, when a value has moved,
names the documents that quote it so the fix is a list to follow rather than a
memory test. Refresh with `--update` after an intended change, and read the
resulting diff as the checklist of documents to update.

The configurations behind these values are declared once in `feature_tour.cc`
and shared between the printed tables and `--numbers`, so the check cannot
verify one set of configurations while the tour prints another.

Timings are deliberately excluded: they are properties of the machine rather
than of the physics, so they are quoted with the machine named beside them and
cannot be diffed. Tolerances are set by what the platform can do to a value,
and are documented next to them in `check_numbers.py`.

## Reproducibility

The Fortran generator seeds its RNG from `UCMUON_SEED` when set. Runs are made
with `OMP_NUM_THREADS=1`, which is required: the RNG is seeded per thread, so a
multi-threaded run is only reproducible at fixed thread count, and output rows
are written under an OpenMP CRITICAL whose ordering is not deterministic.

Reference samples are **not** committed. At 1e6 events x 14 columns x 8 bytes
that would be ~2.5 GB across the matrix. Since the Fortran is reproducible, the
suite commits only the configuration matrix, the calibrated thresholds, and
small golden summaries, and regenerates full samples on demand. This also means
the reference can never go stale relative to the Fortran.

## Running it

```
./run_tests.sh          # header hygiene + test_projection + test_detector + test_geant4
python3 compare_legacy.py --n 20000     # exact event equality vs the Fortran
```

`run_tests.sh` skips the Geant4 test rather than failing when `geant4-config` is
absent, because the header's whole claim is that it needs no dependencies.
`CXX` sweeps the compiler; the Geant4 test deliberately ignores it and uses the
system default, since `geant4-config --cflags` returns the flags Geant4 was
itself configured with (a clang build emits `-Qunused-arguments`, which gcc
rejects) and the prebuilt libraries must be ABI-compatible with whatever links
them. Override with `G4CXX` if needed.

The first check is header hygiene, and it runs before any physics because a
header-only library fails in two ways no single-file test can see: a function
missing `inline` links fine until a *second* translation unit includes the
header, and a broken include guard only shows on the second include within one
file. Both are one compile to rule out. This is not hypothetical here: the
Geant4 helpers had to be moved outside the main include guard after a real bug
in which defining `UCMUGEN_WITH_GEANT4` after a first plain include silently
left `FireG4` undefined.

## Continuous integration

`.github/workflows/ucmugen.yml` runs the suite on four compiler/platform pairs
(ubuntu and macos, gcc and clang) with warnings fatal, plus a second job that
builds the Fortran generator with gfortran and runs the exact event comparison.
A single differing event fails the build. The Geant4 test is not run in CI (it
needs a multi-GB Geant4 installation) and is a local pre-release step.

## Layout

```
run_tests.sh         build and run the C++ suites
test_projection.cc   closed-form checks of the surface projection
test_detector.cc     closed-form checks of detector-directed sampling
test_geant4.cc       the Geant4 hand-off, needs a Geant4 install
check_numbers.py     documented numbers vs the code that produces them
ucmuref/fortran.py   drive the Fortran generator, parse its output
ucmuref/stats.py     probes, null calibration, resolving power
ucmuref/cases.py     the 23-configuration spanning set
measure_power.py     produce the resolving-power table above
```

`cases.py` uses a spanning set rather than a full cross product: spectrum,
angular model and surface are independent code paths, so 7 x 5 x 3 = 105
configurations is redundant. 23 exercise every path, plus the edge cases most
likely to break a port (mono-energetic input, narrow and wide energy windows,
near-horizon acceptance, and each model's natural angular pairing).

## Known deliberate divergence

The Fortran applies no projection weighting on tilted or vertical generation
surfaces: `ucmuon_gen_omp.f90:830-838` transforms positions only, and direction
cosines are bit-identical between `tilt=0` and `tilt=30`. Those configurations
therefore sample the sky *intensity* at points on the surface, not the flux
*through* it. UCMuGen computes `J x max(0, -n.d)` per surface and will not
reproduce the Fortran there by design. Those cases are excluded from the
must-match set and validated analytically instead.
