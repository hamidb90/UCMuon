# PROPOSAL — UCMuon Engine 3

Full stochastic Monte Carlo (Koehne 2013; Alameddine 2024). **Part of the main
UCMuon repo.** Driver: `gui/proposal_driver.py`. Install: `pip install proposal`
(a *system* Python venv — not miniforge/anaconda, due to a pybind11 ABI clash;
`setup.sh` sets up `~/venvs/ucmuon` automatically).

## Produce the benchmark output
```bash
python3 compare_engines.py --engine proposal --source benchmark_surface.dat \
        --depths 1,10,25,50,100,200
```
Writes `PROPOSAL_bench_<depth>m.dat` to the scratch tree
(`geant4_muon_rock_v5/PROPOSAL/`). Timing: `../../results/PROPOSAL_*_timing.txt`.

> ⚠️ Open item: with MCS enabled, PROPOSAL overestimates the scattering angle by
> +15–19 % at ≥66 MWE — likely the medium radiation length X₀ (should be
> ≈26.7 g cm⁻² for Z=11/A=22/ρ=2.65). Survival and exit KE are unaffected and
> match Geant4/MUSIC. See `../../reports/BENCHMARK_FEEDBACK.md` §4.3.
