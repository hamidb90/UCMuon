# Benchmark results

What each file is, and from which run. The monoenergetic population is six
vertical beams (5 to 300 GeV, 10⁵ muons each) through Standard Rock, scored at
1, 10, 25, 50, 100 and 200 m.

| File | Contents | Run |
|---|---|---|
| `benchmark_summary_v2.csv` | six-code summary per depth: transmission, energies, energy loss, angles, lateral displacement | current: Geant4 `run_stdrock_xy1km_20261022_mono` (Standard Rock, slab ±1 km wide; up to 1.4.1 `run_file_20260522_210114`, sodium rock, and in between `run_stdrock_20261019_mono`, ±25 m, which lost muons through its sides) and PHITS (3.37, re-run 2026-10-25 with its rock as near Standard Rock as PHITS allows and a slab ±1 km wide; up to 1.4.1 PHITS 3.36, May 2026, Mg/O rock, ±25 m), MUSIC and PROPOSAL (July 2026), Bethe-Bloch (re-run 2026-10-22 with the thick-target scattering) and UCMuon-MC (re-run 2026-10-23 with its model after 1.4.1: tabulated radiative spectra, spin-1/2 knock-ons, v_cut = 0.01); Geant4's lateral displacement from the source positions (empty up to 1.4.1); written by `analysis/run_benchmark.sh`; all six depths (10 m was excluded up to 1.4.1). Equals `../figures/v2_six_code/benchmark_summary.csv` |
| `comparison_metrics.csv` | each code's quantities and their ratio to Geant4 | current, written by `analysis/run_benchmark.sh` (copy of `../figures/v2_six_code/comparison_metrics.csv`) |
| `engine_timing_20261027.csv` | engine times per depth and per 10⁵ muons, one process each, one session, one machine, with Geant4 and PHITS as reference rows (one core each, timed on every 32nd muon and scaled); the paper's engine times and Fig. 10 | current (UCMuon-MC's model after 1.4.1; PHITS 3.37 on one core instead of the CPU time of 10 OpenMP threads); read its header |
| `engine_timing_20261018.csv` | the same, 2026-10-18 | superseded: UCMuon-MC's model of 1.4.1, and PHITS as the CPU time of the 10 OpenMP threads of a PHITS 3.36 run |
| `engine_timing_20261016.csv` | the same, 2026-10-16, v1.4.0 | superseded: its UCMuon-MC times were for 10 worker processes, not one, and MUSIC predates the new routines |
| `engine_timing_20260730.csv` | the same, July 2026 | superseded |
| `*_timing.txt` | wall times of the individual June and July runs | superseded by `engine_timing_20261027.csv` |
| `backward_validation_20261019.csv` | Engine 5 (backward CSDA flux integrator, mode 1) against UCMuon-MC, MUSIC and PROPOSAL: Guan, zenith 0/60/80/85°, 200 to 3000 m w.e.; panel (a) of the paper's `fig_bmc_validation` | `analysis/validate_backward.py`, 2026-10-19; UCMuon-MC re-run 2026-10-23 |
| `backward_calibration_run_20261020.csv`, `backward_calibration_20261020.csv` | the run behind Engine 5's mode-2 calibration (four spectra, zenith 0/60/80/85°, 100 to 5000 m w.e., UCMuon-MC and MUSIC) and the table built from it | `analysis/validate_backward.py`, `analysis/build_backward_calibration.py`, 2026-10-20; UCMuon-MC re-run 2026-10-23 |
| `backward_validation_offnode_20261020.csv` | mode 2 between the calibration's nodes (zenith 40/70/75/82°, 350 to 4500 m w.e., three engines); panel (b) of `fig_bmc_validation` | `analysis/validate_backward.py`, 2026-10-20; UCMuon-MC re-run 2026-10-23 |
| `mcs_widths_20261022.csv`, `mcs_histograms_20261022.csv` | multiple scattering of every engine per beam and depth: core width as E_s, tail percentiles, lateral core and its model; binned angular and lateral distributions for the paper's `fig_mcs` | `analysis/mcs_distributions.py`, 2026-10-22, after the thick-target limit (UCMuon-MC and Bethe-Bloch re-run); 2026-10-23 with the Geant4 ±1 km slab and UCMuon-MC after its energy-loss changes |
| `geant4_survival_by_beam.csv` | Geant4 survival per beam and depth of the monoenergetic suite (the per-crossing Geant4 output, 400 MB, is not shipped); read by the paper's survival table | `analysis/geant4_survival_by_beam.py`, from `run_stdrock_xy1km_20261022_mono` |
| `tev_beams_20261023.csv` | 1, 3 and 10 TeV beams through 1, 2 and 3 km w.e. of Standard Rock, 10⁴ muons per cell: survival and exit energy of UCMuon-MC, MUSIC and PROPOSAL | `analysis/run_tev_bench.py`, 2026-10-23 |
| `benchmark_summary_4code.csv`, `muon_depth_stats.csv` | early tests (May 2026) with earlier sources (10⁵ and 10⁴ muons, not the population above) | kept for the record; not comparable with the files above |
