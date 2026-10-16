# Bethe–Bloch + Highland MSC — UCMuon Engine 3

Deterministic CSDA energy loss (PDG/Groom) with Highland multiple scattering.
**Part of the main UCMuon repo.** Binary: `bin/ucmuon_transport_bb_omp`
(source: `src/transport/bethe_bloch/`); a pure-Python path also exists
(`gui/ucmuon_bb_driver.py`), and it is the one the benchmark runs. No external
files needed.

## Produce the benchmark output
```bash
bash benchmark/analysis/run_mono_bench.sh "1 10 25 50 100 200"
```
runs the four UCMuon engines (this one included) at every depth and writes
`BB_bench_<depth>m.dat` to the scratch tree (`geant4_muon_rock_v5/BB/`).
Timing: `../../results/BB_timing.txt`.

> As a deterministic CSDA engine it carries no range straggling, so it returns
> 100 % survival up to the sharp CSDA cutoff and brackets the stochastic engines
> from above (survival 0.3 to 0.6 pp above Geant4 from 25 to 200 m). Its
> survivors are not the stochastic engines' survivors, so their mean exit
> energy is 1 to 2 % below Geant4 down to 50 m and 3.9 and 6.6 % below at 100
> and 200 m, and their mean scattering angle 3 to 14 % below (1 to 200 m;
> `../../results/benchmark_summary_v2.csv`, October 2026 run with the
> scattering variance accumulated along the path).
