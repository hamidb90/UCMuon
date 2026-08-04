# UCMuGen vs EcoMug

Two measured comparisons behind Section 6 of the paper. Both give the two
generators the **identical differential flux**, so nothing here depends on their
default parametrisations differing: what is compared is the machinery, not the
physics input.

EcoMug is GPL-3.0 and is **not redistributed** with UCMuon. Download `EcoMug.h`
separately (see `the EcoMug project page`) and point
the include path at it. These programs live here rather than in
`ucmugen/validation/` for exactly that reason: the validation suite must stay
runnable with no Geant4 and no EcoMug.

## rate_crosscheck.cc

Two independently written rate integrators, one differential flux. UCMuGen's
`rate()` and EcoMug's `GetAverageGenRateAndError` with a custom `J` are
integrating the same thing, so they must agree; agreement validates both.

```
c++ -std=c++17 -O2 -I../include -I<EcoMug dir> \
    rate_crosscheck.cc -o rate_crosscheck && ./rate_crosscheck
```

EcoMug's custom-J estimator samples uniformly in momentum over a steeply
falling spectrum, so it is noisy: the program averages over eight seeds and
uses the spread of the means as the uncertainty rather than trusting the quoted
per-run error.

It also shows the live-time gap: with a user-supplied `J`, EcoMug's
`GetEstimatedTime` returns `0` by construction (`EcoMug.h`, `if (mCustomJ)
return 0.;`), so the exposure has to be derived externally.

## cylinder_g4.cc

Runs inside Geant4 (tested with 11.4.1). A 1 m x 1 m x 10 cm plate at the
origin, primaries generated through four different surfaces, tracked with
transportation only through a vacuum world, counting muons that enter the plate.
Physics processes are deliberately absent: what is under test is geometric
acceptance, and adding scattering would only add noise to it.

```
c++ -std=c++17 -O2 $(geant4-config --cflags) \
    -I../include -I<EcoMug dir> \
    cylinder_g4.cc $(geant4-config --libs) -o cylinder_g4
for i in 0 1 2 3; do ./cylinder_g4 $i 2000000; done
```

The four cases are chosen so the experiment identifies a cause rather than only
a discrepancy:

| # | generation surface | role |
|---|---|---|
| 0 | UCMuGen, flat plane above the detector | reference: a horizontal plane has no cap ambiguity |
| 1 | UCMuGen, cylinder **with** caps | must reproduce 0, since the rate into a detector cannot depend on the surface used to generate through |
| 2 | UCMuGen, cylinder, caps **disabled** | isolates the missing-caps effect within one codebase |
| 3 | EcoMug, cylinder (lateral surface only) | the tool under comparison |

Case 2 is what makes this conclusive. If case 3 matched case 1 the caps would be
irrelevant; if case 3 matches case 2, the missing caps are the entire
explanation and no other difference between the two tools contributes.

### Two traps worth knowing

**The reference plane must not truncate the acceptance.** A 70 degree muon
starting at height *z* travels *z*·tan(70) laterally, so the plane's half-width
must exceed that plus the detector half-width or the reference is silently too
low. At the first attempt it was not, and the "reference" came out below the
cylinder.

**The surface-rate estimator has its own Monte Carlo error.** Reporting only the
binomial error on the hit fraction understates the total and turns a 1%
agreement into an apparent 2 sigma discrepancy. Both terms are propagated.

### Fairness

For a flat horizontal plate, EcoMug's `Sky` mode is the appropriate generation
surface and would not show this deficit. The point is narrower and still worth
making: the cylindrical surface is incomplete, it is the natural choice when
enclosing a detector, and selecting it produces a badly wrong normalisation with
no diagnostic.
