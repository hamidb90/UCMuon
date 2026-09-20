# Bethe–Bloch + Highland MSC — UCMuon Engine 2

Deterministic CSDA energy loss (PDG/Groom) with Highland multiple scattering.
**Part of the main UCMuon repo.** Binary: `bin/ucmuon_transport_bb_omp`
(source: `src/transport/bethe_bloch/`); a pure-Python path also exists
(`gui/ucmuon_bb_driver.py`). No external files needed.

## Produce the benchmark output
```bash
python3 compare_engines.py --engine bb --source benchmark_surface.dat \
        --depths 1,10,25,50,100,200
```
Writes `BB_bench_<depth>m.dat` to the scratch tree (`geant4_muon_rock_v5/BB/`).
Timing: `../../results/BB_timing.txt`.

> Best **energy** agreement with Geant4 (±2 % at all depths). As a deterministic
> CSDA engine it carries no range straggling, so it returns 100 % survival up to
> the sharp CSDA cutoff (it brackets the stochastic engines from above) and
> underestimates MCS angles by ~40 % — expected, do not tune. See
> `../../reports/BENCHMARK_FEEDBACK.md`.
