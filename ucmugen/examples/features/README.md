# Feature tour

Every UCMuGen option in one runnable file, with the numbers each one produces.

The Geant4 example next door is deliberately minimal, because its job is to be
copied into your application. This one is the opposite: it is meant to be read
and run, not copied, and it exists so you can see what a choice does before you
make it.

```bash
c++ -std=c++17 -O2 -I../../include feature_tour.cc -o feature_tour
./feature_tour
```

Runs in a few seconds. No Geant4, no external data. To include the
site-aware PARMA spectrum, generate `UCMuGen_PARMA.h` first (see
`ucmugen/tools/make_parma_header.py`) and add `-DUCMUGEN_TOUR_PARMA`.

`./feature_tour --numbers` prints the same values as machine-readable
`key value` pairs. That is what CI diffs against
`ucmugen/validation/reference/tour_numbers.txt`, so the numbers quoted below
cannot go stale without the build failing. See
`ucmugen/validation/README.md`, "Numbers in the documentation".

## What it covers

| Section | Feature |
|---|---|
| 1 | All eight spectra, and which of them carry an absolute normalisation |
| 2 | All four generation surfaces: `Plane` (flat and tilted), `Disk`, `HSphere`, `Cylinder` |
| 3 | All three detector shapes, directed vs blind sampling, and the speed-up |
| 4 | `setThetaRange` and `setPhiRange` |
| 5 | `setSeed`, `setRng`, reproducibility, `provenance()` |
| 6 | `rate()`, `rateAndError()`, and converting a muon count into live time |
| 7 | The legacy Fortran-compatible driver and its five angular models |

## Three things in the output worth reading twice

**Rates add up.** Zenith 0-30 deg gives 173.5 Hz and 30-70 deg gives 149.6 Hz,
and the full 0-70 deg range gives 323.1 Hz. Halving the azimuth range halves the
rate; taking a quadrant quarters it. Nothing enforces this: the restricted
generators each run their own integral over their own region. That they sum
correctly is a check that the angular ranges are honoured by the normalisation
and not only by `generate()`.

**Directed and blind sampling agree.** For each detector shape the two modes are
compared against their combined Monte Carlo error. Directed sampling is a
variance-reduction technique, not an approximation, so agreement here is a
correctness requirement rather than a nice result.

**The speed-up is a property of your geometry, not of the generator.** It grows
as the detector covers less of the sky:

```
  plate (full size)     acc dir  acc blind  dir us/mu blind us/mu  speed-up
  200 x 200 cm         8.03e-02   1.43e-02       1.69        4.8        3x
  60 x 60 cm           9.60e-02   1.52e-03       1.53       34.6       23x
  20 x 20 cm           9.44e-02   2.42e-04       1.58      208.3      132x
  6 x 6 cm             4.26e-02   4.63e-05       2.67     1085.4      407x
  2 x 2 cm             1.33e-02   1.18e-05       6.07     4189.7      691x
```

(Timings: Apple clang 21, `-O2`, Apple Silicon, one of three runs on
2026-10-16 that agree within a few per cent; they depend on the machine.)

A plate that fills the view from the source plane gains almost nothing, because
blind sampling was already hitting it most of the time. Real muography detectors
sit at the bottom of that table.

The acceptance columns and the speed-up column are deliberately separate
numbers. The ratio of the acceptances is the saving in *proposals*; the
wall-clock speed-up is 1.6 to 3 times smaller, because a directed proposal
costs more than a blind one (a cone has to be built and a ray-geometry test
run). Quoting the acceptance ratio as a speed-up would overstate the win.

## The surfaces are the part with no equivalent elsewhere

The projection cosine lives in the `Surface`, not in a per-geometry formula, so
a tilted plane, a dome and a cylinder are all first-class and none of them needs
a hand-derived integrator. Section 2 prints a flat and a tilted plane of
identical area side by side; the tilted one has a lower rate because it presents
less of itself to a sky that is brightest overhead.

This is also where UCMuGen and the Fortran generator genuinely differ: the
Fortran applies no projection weighting on a tilted or vertical surface, so its
sky intensity there is not a through-surface flux.
