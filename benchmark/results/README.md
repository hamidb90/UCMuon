# Benchmark results

What each file is, and from which run. The monoenergetic population is six
vertical beams (5 to 300 GeV, 10⁵ muons each) through Standard Rock, scored at
1, 10, 25, 50, 100 and 200 m.

| File | Contents | Run |
|---|---|---|
| `benchmark_summary_v2.csv` | six-code summary per depth: transmission, energies, energy loss, angles, lateral displacement | current: Geant4 `run_stdrock_20261019_mono` (Standard Rock; up to 1.4.1 `run_file_20260522_210114`, sodium rock) and PHITS (May 2026), MUSIC and PROPOSAL (July 2026), Bethe-Bloch and UCMuon-MC (re-run 2026-10-16); written by `analysis/run_benchmark.sh`; all six depths (10 m was excluded up to 1.4.1). Equals `../figures/v2_six_code/benchmark_summary.csv` |
| `comparison_metrics.csv` | each code's quantities and their ratio to Geant4 | current, written by `analysis/run_benchmark.sh` (copy of `../figures/v2_six_code/comparison_metrics.csv`) |
| `engine_timing_20261018.csv` | engine times per depth and per 10⁵ muons, one process each, one session, one machine, with Geant4 and PHITS as reference rows; the paper's engine times and Fig. 10 | current (after 1.4.1: new random-number routines); read its header |
| `engine_timing_20261016.csv` | the same, 2026-10-16, v1.4.0 | superseded: its UCMuon-MC times were for 10 worker processes, not one, and MUSIC predates the new routines |
| `engine_timing_20260730.csv` | the same, July 2026 | superseded |
| `*_timing.txt` | wall times of the individual June and July runs | superseded by `engine_timing_20261018.csv` |
| `backward_validation_20261019.csv` | Engine 5 (backward CSDA flux integrator, mode 1) against UCMuon-MC, MUSIC and PROPOSAL: Guan, zenith 0/60/80/85°, 200 to 3000 m w.e.; panel (a) of the paper's `fig_bmc_validation` | `analysis/validate_backward.py`, 2026-10-19 |
| `backward_calibration_run_20261020.csv`, `backward_calibration_20261020.csv` | the run behind Engine 5's mode-2 calibration (four spectra, zenith 0/60/80/85°, 100 to 5000 m w.e., UCMuon-MC and MUSIC) and the table built from it | `analysis/validate_backward.py`, `analysis/build_backward_calibration.py`, 2026-10-20 |
| `backward_validation_offnode_20261020.csv` | mode 2 between the calibration's nodes (zenith 40/70/75/82°, 350 to 4500 m w.e., three engines); panel (b) of `fig_bmc_validation` | `analysis/validate_backward.py`, 2026-10-20 |
| `benchmark_summary_4code.csv`, `muon_depth_stats.csv` | early tests (May 2026) with earlier sources (10⁵ and 10⁴ muons, not the population above) | kept for the record; not comparable with the files above |
