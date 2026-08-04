# UCMuGen

A header-only C++17 cosmic muon generator for Geant4 and standalone use.
Drop `include/UCMuGen.h` in your project and include it. No dependencies.

```cpp
#include "UCMuGen.h"
using namespace ucmugen;

Generator gen;
gen.setSpectrum(Spectrum::Guan)
   .setSurface(std::make_shared<Plane>(100.0, 100.0, Vec3{0, 0, 500}))
   .setEnergyRange(1.0, 1000.0)
   .setThetaRange(0.0, 70.0 * kPi / 180.0);

Muon m = gen.generate();     // position, direction, momentum, kinetic, pdg
double hz = gen.rate();      // absolute rate through the surface, s^-1
double t  = gen.liveTime(1000000);   // seconds of real time 1e6 muons represent
```

Attach a detector and only muons that can reach it are generated, with `rate()`
reporting the rate into it:

```cpp
gen.setDetector(BoxDetector::centred(Vec3{0, 0, -20000}, 100, 100, 50));
```

## What it offers over EcoMug

**Seven spectrum models, not one.** CosmoALEPH, truncated power law, Guan 2015,
Frosin 2025, Gaisser/Bugaev, Reyna-Bugaev, and cosmic electrons. Swapping
between them is one call, which makes a spectrum-model systematic something you
can actually quote.

**Correct through-surface weighting on every surface.** EcoMug hardcodes a
separate `J'` expression per geometry (`EcoMug.h:975`, `:980`, `:1026`), each
needing its own hand-derived MC integrator plus a magic normalisation constant.
All three are the same thing: `J(p, theta) * max(0, -n.d) * sin(theta)`.
Factoring the projection into the surface means one sampling loop, one rate
integrator, every spectrum working with every surface, and a new surface costing
about twenty lines and no new physics.

**Absolute rates for any combination, including custom fluxes.** EcoMug's
`GetEstimatedTime` returns `0.` when a user supplies their own `J`
(`EcoMug.h:664`). Here the same estimator covers every case.

**Momentum-dependent charge ratio.** EcoMug uses a single constant.

**Cylinder caps.** EcoMug's cylinder is lateral-surface-only, so a vertical
muon has `J' = 0` and can never be generated. Measured with
`validation/test_projection.cc`, lateral-only misses **33%** of the rate for a
tall cylinder (h/r = 4) and **89%** for a squat one (h/r = 0.25).

**Detector-directed sampling.** Attach a box, sphere or arbitrarily oriented
cylinder and directions are drawn only inside the cone it subtends, instead of
over the whole sky with almost every proposal thrown away. For a 2 m detector
under 200 m of rock, generated over a 40 m disk, that is a measured **5800x**
speed-up (`validation/test_detector.cc`), and `rate()` becomes the rate *into
the detector*. EcoMug has no equivalent: you generate over the sky and discard
the misses yourself, and the rate estimate does not know you did.

The correctness argument is that the sampled cone is the one subtended by the
detector's bounding sphere, so it provably contains every direction that could
reach the volume, and the per-trial cone solid angle enters the same estimator
that produces the rate. Directed sampling and a plain acceptance cut share no
sampling code, and the suite checks that they agree.

**A direction you can use directly.** `Muon::direction` is a Cartesian unit
vector and `Muon::kinetic` is kinetic energy, so there is no trigonometry to get
wrong. EcoMug returns a theta already flipped by pi, and every caller re-derives
the same three lines.

## Validation

The Fortran generator of [UCMuon](https://github.com/hamidb90/UCMuon) is the
reference implementation. UCMuGen reproduces its RNG bit-for-bit, so agreement
is checked as **exact event equality**, not a statistical tolerance:

    20 of 20 applicable configurations produce bit-identical events
    (3 excluded by design, see below), at 20k events each and
    spot-checked at 500k.

The surface projection has no Fortran counterpart, so it is validated against
closed forms instead (`validation/test_projection.cc`), including the tilted
plane law `R(alpha)/R(0) = (1 + cos alpha)/2`. That is deliberately *not* the
naive projected-area law `cos(alpha)`, which holds only for a parallel beam or a
full-sphere isotropic field; at 45 degrees the two differ by 21%.

Detector-directed sampling is validated the same way, in
`validation/test_detector.cc`, because the way variance reduction fails is by
staying plausible while changing the normalisation. In the far field the rate
into a detector is fixed by its projected area alone, which gives a closed form
for four shapes at once and exercises all three intersection routines: a sphere
(`pi R^2`), a box (`4 hx hy`), an upright cylinder (`pi r^2`, so its length must
not matter) and the same cylinder on its side (`2 r L`). All four agree to
better than 1e-3, at MC precision. Two further tests pin the machinery rather
than the geometry: directed sampling must reproduce the plain acceptance cut,
with which it shares no sampling code, and a detector large enough to intercept
every ray must reproduce the open-sky rate.

The Geant4 hand-off is validated against a real Geant4 build
(`validation/test_geant4.cc`, checked on 11.4.1). It is the one boundary the
dependency-free tests cannot reach, and it is where the conversions live: cm to
mm, GeV to MeV, total to kinetic energy, and a PDG code whose sign is the
opposite of the naive guess. Two generators run on the same seed, one through
`FireG4` and one calling `generate()` directly, so the primary vertex Geant4
ends up holding can be compared against the muon that produced it. Over 5000
events every conversion is exact, the primaries read back out of `G4Event` are
still aimed at the detector, and the same Geant4 seed reproduces the sequence.

Run everything with `validation/run_tests.sh`, which skips the Geant4 test
rather than failing when `geant4-config` is absent.

One trap worth recording, since it cost a debugging cycle: a detector has depth,
so the face that intercepts the rays sits nearer than its centre by its own
half-thickness. At 200 m depth a 6 m cylinder is off from the naive
centre-distance formula by `(20000/19700)^2 = 3%`. The test places detectors at
20 km so the effect drops to 3e-4. A sphere is the one shape immune to it.

See `validation/README.md` for the methodology, including why the statistical
suite calibrates every threshold by measurement rather than choosing tolerances,
and the measured resolving power (~0.3% at N = 800k, scaling as 1/sqrt(N)).

## Deliberate divergence from the Fortran

The Fortran applies no projection weighting on tilted or vertical generation
surfaces: its tilt transform moves positions only, and direction cosines are
bit-identical between `tilt=0` and `tilt=30`. Those modes therefore sample the
sky *intensity* at points on the surface rather than the flux *through* it.
UCMuGen is correct-only here and does not reproduce them.

## Site-aware flux (PARMA)

This is the capability no other drop-in muon generator has, and the reason to
prefer UCMuGen over EcoMug for muography: flux that knows where and when you
are.

```cpp
#include "UCMuGen.h"
#include "UCMuGen_PARMA.h"

ucmugen::parma::Site site;
site.cutoff_GV  = 3.25;                                    // geomagnetic cutoff
site.depth_gcm2 = ucmugen::parma::depth_from_altitude(3.0, 45.0);   // 3 km
ucmugen::parma::install(site);

gen.setSpectrum(ucmugen::Spectrum::Parma);
```

Altitude, geomagnetic cutoff rigidity and solar epoch all move the answer:
measured with `validation/test_parma.cc`, vertical 10 GeV/c intensity is
**1.30x** sea level at 3 km and **1.44x** at 5 km, and the rate through a
horizontal detector rises **1.64x** at 3 km once the whole spectrum is
integrated. A sea-level parametrisation cannot express any of that, so for a
volcano observatory it is not a refinement, it is the difference between a
right and a wrong normalisation.

`UCMuGen_PARMA.h` is generated by `tools/make_parma_header.py` from JAEA's
official `parma_cpp` release, taking the routines verbatim and applying three
mechanical edits (file reads become embedded-table reads; two species loops are
restricted to muons so only the muon tables need embedding). The generated
header is 269 kB with 232 kB of embedded tables and needs no data files at
runtime. Cutoff rigidity from latitude and longitude, and the W index from a
date, need the two tables too large to embed (1.8 MB and 233 kB), so those are
loaded from a data directory you point at; the second is also revised by JAEA
over time, which is a reason not to freeze a copy.

The restriction to muons is the one edit that could silently change the physics,
so it is checked rather than assumed: at generation time the header is compared
against the stock JAEA build over **129,600 points** covering every W index,
cutoff, depth, geometry mode, energy and angle, and agrees **bit-for-bit**.

> **Licence:** `UCMuGen_PARMA.h` contains JAEA code and data and is **not MIT**.
> The EXPACS conditions grant redistribution and modification for
> **non-commercial** use and forbid commercial use without prior agreement, so
> including this header makes your build non-commercial-only. `UCMuGen.h` alone
> is unaffected. Cite Sato 2015 and Sato 2016. Delete the file and everything
> else still works; the suite skips its test.

## Layout

```
include/UCMuGen.h              the library
examples/geant4/               a complete, multithreaded Geant4 application
examples/features/             every option, with the numbers each produces
comparison/                    measured UCMuGen vs EcoMug
validation/                    the test suite (see its README)
```

Start with `examples/geant4/` to integrate, and `examples/features/` to decide
what to integrate: the first is the minimum you need and is meant to be copied,
the second exercises all eight spectra, all four surfaces, all three detector
shapes, the angular ranges, seeding and normalisation, and prints what each
choice does to the rate. It needs no Geant4 and runs in two seconds.

## Citing

Please cite the UCMuon suite, Zenodo concept DOI
[10.5281/zenodo.20826984](https://doi.org/10.5281/zenodo.20826984), plus the
paper for whichever spectrum model you use.

## Licence

MIT, except PARMA as noted above.
