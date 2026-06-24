# UCMuon-MC — UCMuon Engine 1 (flagship)

Native pure-Python stochastic engine: Groom dE/dx + Poisson-sampled hard
radiative events + Highland MSC + muon decay. **Part of the main UCMuon repo.**
Driver: `gui/ucmuon_stochastic_driver.py` (numpy/scipy only — no compiler, no
external tables, cross-platform).

## Produce the benchmark output
```bash
python3 compare_engines.py --engine ucmuon --source benchmark_surface.dat \
        --depths 1,10,25,50,100,200
```
Writes `UCMuon_bench_<depth>m.dat` to the scratch tree
(`geant4_muon_rock_v5/UCMuon/`). Timing: `../../results/UCMuon_*_timing.txt`.

> Tracks the full-MC codes within a few percent (aggregate survival bias −1.2 %,
> exit KE ±3.5 %) and is systematically the lowest survivor in high-loss cells.
> **Use the v2 (post-fix) outputs** — the flat `output/*_bench_*.dat` in the main
> repo are an earlier v1 run with inflated exit KE (126.0 vs 123.74 GeV at 100 m).
> See `../../reports/BENCHMARK_FEEDBACK.md` §4 and `manuscript/PUBLICATION_TODO.md`.
