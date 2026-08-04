# UCMuGen in Geant4: a complete, runnable example

Builds and runs on its own. Tested against Geant4 11.4.1.

```bash
cmake -S . -B build && cmake --build build -j
cd build
./ucmugen_example run.mac       # batch: prints a rate and a live time
./ucmugen_example tracks.mac    # writes trajectories for a figure
./ucmugen_example               # interactive, with visualisation
```

## What you actually need to copy

`PrimaryGeneratorAction.hh` and `PrimaryGeneratorAction.cc`. That is the whole
integration. `ucmugen_example.cc` is scaffolding so that this directory runs by
itself: a plate of scintillator, `FTFP_BERT`, a step counter and a `main()`. You
already have those.

Inside `GeneratePrimaries` the integration is one line:

```cpp
ucmugen::FireG4(fGen, fGun, event);
```

which converts cm to mm, GeV to MeV, total energy to kinetic energy, and the
charge to a PDG code whose sign is the opposite of the naive guess.

Define `UCMUGEN_WITH_GEANT4` before including `UCMuGen.h` to get `FireG4`. The
header pulls in the Geant4 headers it needs, so the include order in your own
file does not matter.

## What the run prints

```
[UCMuGen] rate into the detector : 325.517 Hz
[UCMuGen] 1e6 muons correspond to: 3072.03 s of live time
---------------- UCMuGen example ----------------
  generated muons        : 100000
  muons entering the plate: 99928  (99.928 %)
  rate into the detector : 325.517 Hz
  live time of this run  : 307.203 s
  measured hit rate      : 325.283 Hz
-------------------------------------------------
```

Two things in that output are worth dwelling on, because they are the reason to
use this generator rather than a spectrum sampler.

**The predicted rate and the measured rate agree to 0.07%.** `rate()` is a Monte
Carlo integral computed before the run; the measured value is Geant4 counting
muons that actually crossed into the plate. Nothing forces them to agree, and
that they do is an end-to-end check of the flux normalisation, the surface
projection and the unit conversions at once. It also means the run has an
*exposure*: 10⁵ muons here are 307 seconds of real sky, so simulated counts
convert into a rate.

**99.93% of generated muons reach the plate, not 100%.** Detector-directed
sampling guarantees the *ray* intersects the detector, not that the muon
survives to it. The missing 0.07% is decay in flight across 3 m of air, which is
the right order for these momenta. If you disable the detector (comment out
`setDetector`) the run still works and gives the same rate, but generates
hundreds of times more muons to get there.

## Options

Location- and date-aware PARMA flux, off by default because it carries JAEA's
non-commercial licence:

```bash
cmake -S . -B build -DUCMUGEN_PARMA=ON
```

Change the site in `PrimaryGeneratorAction.cc` with
`parma::depth_from_altitude(altitude_km, latitude_deg)`. At 3 km the vertical
10 GeV/c intensity is 1.30x its sea-level value, so for a mountain or volcano
site this is not a refinement.

## Figures

`tracks.mac` renders muon trajectories and their secondaries to a PNG through
Geant4's offscreen renderer, so it needs no display and works on a batch node.
It is there to let you look at your own setup; it is not used for any figure in
the paper.

For the published figure of blind versus directed sampling, see
`ucmugen/validation/dump_rays.cc`, which writes the generated rays out to be
plotted directly and needs neither Geant4 nor a screenshot.
