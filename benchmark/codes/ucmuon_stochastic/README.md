# UCMuon-MC — UCMuon Engine 1 (flagship)

Native pure-Python stochastic engine: PDG 2024 per-process energy-loss tables,
sampled bremsstrahlung and photonuclear losses (pair production continuous),
δ-ray straggling, Highland MCS accumulated along the path, muon decay.
**Part of the main UCMuon repo.**
Driver: `gui/ucmuon_stochastic_driver.py` (numpy/scipy only — no compiler, no
external tables, cross-platform).

## Produce the benchmark output
```bash
bash benchmark/analysis/run_mono_bench.sh "1 10 25 50 100 200"
```
runs the four UCMuon engines (this one included) at every depth and writes
`UCMuon_bench_<depth>m.dat` to the scratch tree
(`geant4_muon_rock_v5/UCMuon/`). Timing: `../../results/UCMuon_*_timing.txt`.

> Against Geant4 (`../../results/benchmark_summary_v2.csv`, October 2026 run):
> survival within 0.25 pp at every depth, between MUSIC and PROPOSAL in the
> survival matrix; mean exit energy 0.9 to 4.6 % low from 25 to 200 m, as
> MUSIC (0.9 to 4.1 %) and PROPOSAL (0.8 to 4.3 %); mean scattering angle
> 0.8 to 5.0 % low from 1 to 200 m. Use the `run_mono_bench.sh` outputs: the
> flat `output/*_bench_*.dat` in the main repo are an earlier v1 run with
> inflated exit KE (126.0 vs 123.74 GeV at 100 m).
