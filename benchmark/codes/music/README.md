# MUSIC — UCMuon Engine 1

Kudryavtsev stochastic Monte Carlo. **Part of the main UCMuon repo**, not
duplicated here. Binary: `bin/ucmuon_transport_music_omp`
(source: `src/transport/music/`).

> ⚠️ MUSIC source (`music.f`, `music-crosssections.f`) is **not redistributed** —
> obtain it from Prof. V. Kudryavtsev (`v.kudryavtsev@sheffield.ac.uk`) and place
> it in `src/transport/music/`, then `bash setup.sh`. See `docs/MUSIC_FILES.md`.

## Produce the benchmark output
From the repo root, `compare_engines.py` feeds the shared `benchmark_surface.dat`
through this engine at every depth and writes `MUSIC_bench_<depth>m.dat`:
```bash
python3 compare_engines.py --engine music --source benchmark_surface.dat \
        --depths 1,10,25,50,100,200
```
Raw `.dat` outputs land in the git-ignored scratch tree (`geant4_muon_rock_v5/MUSIC/`).
Timing: `../../results/MUSIC_timing.txt`. Result: agrees with PROPOSAL to <0.7 pp
on survival and <0.4 % on exit KE (see `../../reports/BENCHMARK_FEEDBACK.md`).
