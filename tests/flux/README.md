# Flux tests

Checks on the surface muon spectra, the generator's sampling and the live
time, added with the September 2026 flux audit
(`docs/FLUX_NORMALISATION_AUDIT.md`).

| script | what it checks | needs |
|---|---|---|
| `test_flux_consistency.py` | the same spectrum gives the same I(p, θ), vertical flux and horizontal-surface rate in the Fortran generator, UCMuGen and `gui/fast_flux_estimator.py`; UCMuGen's `rateAndError` matches the quadrature; UCMuGen and Python warn alike | `c++`, `gfortran`, numpy |
| `test_flux_reference.py` | each spectrum against PDG, CosmoALEPH Table 1 and the digitised Reyna 2006 Fig. 3 data, with the tolerances in its header (PARMA when `UCMuGen_PARMA.h` is present) | numpy; `c++` for PARMA |
| `test_joint_sampling.py` | angular mode 6 of `bin/ucmuon_gen_omp` against `ucmugen::Generator` in bins of energy and zenith angle, six configurations; legacy mode 4 as a control that must fail | `make ucmuon_gen_omp`, `c++` |
| `test_backward_mc_range.py` | the backward MC warns, with an estimate of what is left out, when its range table (10 TeV, 6889 m w.e.; 2 TeV up to 1.3.2) cuts a flux, instead of returning a silent zero; the terrain engine prints one summary line per map and warns below a spectrum's fitted range | numpy |
| `test_live_time.py` | the GUI live time (`gui/live_time.py`, T = N_tried / R) against UCMuGen's `liveTime` for plane, disk and hemisphere sources, with and without a detector filter; weighted legacy modes | `make ucmuon_gen_omp`, `c++` |
| `test_generator_fixes.py` | the generator fixes after 1.3.2: PARMA "W from the date" after the last year of its table falls back to the latest W with a warning (it was W = 0); the MPI generator keeps E_min below 1 GeV, generates electrons for spectrum 8, and writes exactly N rows without a detector | `make ucmuon_gen_omp`; `make ucmuon_gen` and `mpirun` for the MPI part (skipped without them) |

```
python3 tests/flux/test_flux_consistency.py
python3 tests/flux/test_flux_reference.py
python3 tests/flux/test_backward_mc_range.py
python3 tests/flux/test_joint_sampling.py     # --quick: CI size
python3 tests/flux/test_live_time.py          # --quick: CI size
python3 tests/flux/test_generator_fixes.py
```

All exit 0 on success and 1 on failure, and run in the `flux` job of
`.github/workflows/ucmugen.yml`.

Helpers: `intensity_dump.cc` / `intensity_dump.f90` print the same grid and
integrals from UCMuGen and from `src/generator/ucmuon_source_module.f90`
(`spectrum_intensity`); `parma_dump.cc` does it for PARMA; `ucmugen_sample.cc`
generates events and rates with `ucmugen::Generator` on a given surface and
detector; `_harness.py` drives the Fortran generator (OpenMP and MPI, detector
filter included) and bins the events.

The detector-filter intersection test is in `tests/geometry/`.
