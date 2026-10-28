# PHITS — external reference code

PHITS input decks and helper scripts for the second external code of the
transport benchmark (PHITS 3.37).
**Physics:** ATIMA dE/dx, Lynch–Molière MSC, with muon brem/pair/nuclear
processes enabled (see header table in `../../reports/BENCHMARK_FEEDBACK.md` §3).

**Rock:** PHITS sets the mean excitation energy only for water, so the rock is
the single nuclide ²²Na: Standard Rock's Z = 11 and A = 22 (Z/A = 0.5, and
Standard Rock's Z²/A for the radiative losses), but sodium's I = 149 eV
instead of 136.4 eV, about 0.7 % less ionisation loss. Up to UCMuon 1.4.1 it
was Mg/O 75:25 by mass (Z/A = 0.495).

**Geometry:** a 200 m slab, ±1 km wide as in the Geant4 reference, so that no
muon leaves through its sides (the beams start within ±25 m; with a slab only
±25 m wide, up to 1.1 % of a beam was lost at 100 to 200 m).

## Files here
| File | Purpose |
|---|---|
| `muon_rock.inp` | monoenergetic benchmark: six beams, 6×10⁵ muons, 6 depth tallies |
| `muon_rock_realistic.inp` | realistic-source benchmark (10⁵ muons, up to 2.5 TeV); tallies run to 10 TeV |
| `convert_source_to_phits.py` | converts the common source → PHITS s-type=17 dump |
| `read_phits_output.py` | parses the raw tally `.out` files → `phits_summary.csv` |
| `run_phits_mpi.sh` | runs a deck under MPI on a Slurm cluster (one master, N−1 workers) |
| `run_phits_chunks.sh`, `merge_chunks.py` | run a deck on one machine as K independent single-thread chunks (every K-th muon, its own seed), P at a time, and merge their summaries |
| `phits_summary.csv`, `phits_timing.txt` | parsed monoenergetic results + timing |
| `phits_summary_realistic.csv` | parsed realistic-source results, read by `../../analysis/analyze_sixcode.py` |
| `batch.out`, `muon_source_phits_info.txt` | run metadata |

Up to 1.4.1 the realistic-source summary was not shipped
(`analyze_sixcode.py` read it from the scratch tree); its numbers are also in
`benchmark/run_realistic_20260614/benchmark_sixcode_summary.csv`.

Both summaries are from runs of 2026-10-25 on Lemaitre4 with `run_phits_mpi.sh`:
the monoenergetic deck with `maxcas = 1000`, `maxbch = 600` on 60 workers, the
realistic one with `maxcas = 200`, `maxbch = 500` on 50 workers.

## How to run the benchmark
```bash
# 1. Convert the shared source to PHITS format (once)
python3 convert_source_to_phits.py <benchmark_surface>.dat muons_for_phits.dat

# 2a. Run PHITS on one machine (set file(1) in [Parameters] to the PHITS directory;
#     maxcas x maxbch must equal the number of source records)
phits muon_rock.inp        # or: /path/to/phits/bin/phits.sh muon_rock.inp

# 2b. or on a Slurm cluster (see the script's header for the sizes used)
sbatch --ntasks=61 --export=ALL,PHITSPATH=/path/to/phits \
       run_phits_mpi.sh muon_rock.inp muons_for_phits.dat RUN_DIR

# 2c. or on one machine: 32 single-thread chunks, 8 at a time, parsed and merged
#     into RUN_DIR/phits_summary.csv (no step 3 needed)
PHITSPATH=/path/to/phits bash run_phits_chunks.sh muon_rock.inp muons_for_phits.dat RUN_DIR 32 8

# 3. Parse the tallies into a summary CSV (in the run directory)
python3 read_phits_output.py --n 600000
```

The inclined-beam decks (`muon_rock_th60.inp`, `_th80`, `_th85`: the planes at
L cos(theta)) and their dumps are written by
`../../analysis/make_inclined_inputs.py`; their summaries of 2026-10-27 were made
with `run_phits_chunks.sh` on an Apple M4 (32 chunks, 8 at a time). Runs of
the same decks with `run_phits_mpi.sh` on Lemaitre4 (other seeds) agree with
them within 0.04 percentage points in survival and 0.2 % in mean exit energy.
PHITS's T-Cross counts every crossing, so near the horizon a muon that
crosses a plane twice counts twice; `../../analysis/range_probe_50gev.py`
dumps every crossing and finds that at most 0.3 % of the muons do.

## Timing
The paper's PHITS time is one core of the Apple M4, PHITS's own single-thread
binary, on every 32nd muon of the source (18 750 of the monoenergetic, 3 125
of the realistic), scaled by 32, as for Geant4. Cluster times are not PHITS's
speed: PHITS reads its CPU clock about 10⁵ times per muon for the
time-by-process table in `phits.out`, and on Lemaitre4 each read costs several
system calls (`getrusage`, `uname`, and with `TZ` unset a `stat` of
`/etc/localtime`), so the same muons took 6 to 24 times longer per core there
than on the M4. Up to 2026-10-24 the paper quoted the CPU time summed over the
10 threads of an OpenMP run of PHITS 3.36 (371 146 s).

> PHITS produces aggregate tallies only — there is no per-muon output, so MCS
> angles and lateral displacement are unavailable for PHITS (shown as `—` in the
> comparison).
