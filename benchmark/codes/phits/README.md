# PHITS — external reference code

PHITS input deck and helper scripts for the second external reference.
**Physics:** ATIMA dE/dx, Lynch–Molière MSC, with muon brem/pair/nuclear
processes enabled (see header table in `../../reports/BENCHMARK_FEEDBACK.md` §3).

## Files here
| File | Purpose |
|---|---|
| `muon_rock.inp` | PHITS input (slab geometry, 6 depth tallies, `maxcas`×`maxbch`=6e5) |
| `convert_source_to_phits.py` | converts the common source → PHITS s-type=17 dump |
| `read_phits_output.py` | parses the raw tally `.out` files → `phits_summary.csv` |
| `phits_summary.csv`, `phits_timing.txt` | parsed results + timing |
| `batch.out`, `muon_source_phits_info.txt` | run metadata |

## How to run the benchmark
```bash
# 1. Convert the shared source to PHITS format (once)
python3 convert_source_to_phits.py <benchmark_surface>.dat muons_for_phits.dat

# 2. Run PHITS
phits muon_rock.inp        # or: /path/to/phits/bin/phits.sh muon_rock.inp

# 3. Parse the tallies into a summary CSV
python3 read_phits_output.py
```

> PHITS produces aggregate tallies only — there is no per-muon output, so MCS
> angles and lateral displacement are unavailable for PHITS (shown as `—` in the
> comparison). Known open item: exit-KE diverges to −12 % at 200 m (see reports).
