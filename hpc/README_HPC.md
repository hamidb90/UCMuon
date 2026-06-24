# UCMuon — HPC Cluster Guide
### Lemaitre4 / CECI (Louvain-la-Neuve, Belgium)

**Hamid Basiri** - UCLouvain Muography Group - hamid.basiri@uclouvain.be

---

## Directory structure

```
UCMuon/                              <- project root; always run sbatch from here
|
|-- Makefile
|-- README.md
|-- requirements.txt                 <- Python deps (GUI only)
|
|-- ucmuon_gen                       <- binary: MPI+OMP generator
|-- ucmuon_transport_music           <- binary: MPI+OMP MUSIC transport
|-- ucmuon_transport_bb              <- binary: MPI+OMP Bethe-Bloch transport
|-- ucmuon_to_phits                  <- binary: PHITS converter (no MPI needed)
|-- music-eloss-rock.dat             <- generated on first MUSIC run (must be in root)
|-- music-cross-sections-rock.dat    <- generated on first MUSIC run (must be in root)
|
|-- src/
|   |-- generator/                   <- generator Fortran source
|   |-- parma/                       <- PARMA spectrum model
|   |-- music/                       <- MUSIC transport source
|   |-- bethe_bloch/                 <- Bethe-Bloch transport source
|   |-- common/                      <- shared Fortran support files
|   `-- phits/
|       |-- ucmuon_to_phits.f90      <- PHITS converter (Fortran, fast)
|       `-- ucmuon_to_phits.py       <- PHITS converter (Python, local use)
|
|-- hpc/
|   |-- README_HPC.md                <- this file
|   |-- run_ucmuon_gen.sh            <- SLURM: generator
|   |-- run_ucmuon_transport.sh      <- SLURM: transport (MUSIC + BB)
|   |-- input_params.dat             <- generator input (annotated)
|   |-- input_transport_music.dat    <- MUSIC transport input (annotated)
|   `-- input_transport_bb.dat       <- BB transport input (annotated)
|
|-- data/                            <- static data (never modified by runs)
|   |-- music-double-diff-rock.dat   <- required for MUSIC init_tables=0
|   |-- music-double-diff-water.dat
|   `-- EXPACS/parma/input/          <- standalone PARMA v4.10 data files
|
|-- logs/                            <- SLURM .out/.err files (auto-created)
|
|-- output_<JOBID>/                  <- one folder per job (auto-created)
|   |-- ucmuon_selected.dat          <- generator: detector-aimed muons
|   |-- ucmuon_selected_phits.dat    <- PHITS dump of above (auto-generated)
|   |-- ucmuon_underground.dat       <- transport: all muons with alive flag
|   |-- ucmuon_underground_phits.dat <- PHITS dump of above (alive==1 only)
|   `-- input_used.dat               <- copy of input for reproducibility
|
`-- build/                           <- object files (auto-created by make)
```

**Note on MUSIC table files:** `music.f` opens energy-loss tables by plain
filename with no path. They must stay in the project root (`UCMuon/`), not
in a subdirectory. The code copies them from `data/` if needed automatically.

---

## One-time setup

### 1. Load MPI toolchain

Add to `~/.bashrc` so it loads automatically on every login:

```bash
echo "module load releases/2023b" >> ~/.bashrc
echo "module load foss/2023b"     >> ~/.bashrc
source ~/.bashrc
```

### 2. Build binaries

```bash
cd ~/UCMuon
make hpc                    # builds all MPI+OMP binaries + PHITS converter
```

Individual targets:
```bash
make ucmuon_gen
make ucmuon_transport_music
make ucmuon_transport_bb
make ucmuon_to_phits        # no MPI needed; fast build
```

### 3. MUSIC data files

Place `music-double-diff-rock.dat` in `data/` (from the MUSIC distribution).
The energy-loss tables (`music-eloss-rock.dat`, `music-cross-sections-rock.dat`)
are generated automatically on the first run with `init_tables=0`.

---

## How to choose the number of MPI ranks

### Hard rule

```
nranks <= N_muons
```

If you request more ranks than muons, some ranks get zero muons. They still
initialise (wasting a CPU slot) but produce no output. This is not a crash
and the total output is still correct, just inefficient.

### Practical rule

Each rank should have at least ~1000 muons to amortise MPI overhead:

```
nranks <= N_muons / 1000
```

### Main rule: match free CPUs

```
nranks = free CPUs on best node   (check with `status`)
```

### Quick reference

| N_muons | Max useful nranks | nranks to set |
|---|---|---|
| 10,000 | 10 | 10 (1000/rank) |
| 100,000 | 100 | 26 (free CPUs) |
| 1,000,000 | 1000 | 26-56 |
| 5,000,000 | 5000 | 26-56 |
| 10,000,000 | 10000 | 56-128 (multi-node) |

### Does it always produce exactly N muons? YES.

The work split is exact in both codes:

```
base    = floor(N / nranks)
rank 0  = base + (N mod nranks)    <- absorbs remainder
ranks 1..nranks-1 = base each
sum     = N exactly
```

For the **generator**: each rank runs rejection sampling until it saves its
exact quota. The acceptance rate (geometry, detector filter) affects wall time
but never the final muon count.

For **transport**: rank 0 reads the full input file and scatters exactly
`base` muons to each rank via `MPI_Scatterv`. Every input muon is transported.

### OMP threads

For Lemaitre4, pure MPI (`--cpus-per-task=1`) is recommended for both
generation and transport. Use hybrid only when the number of free CPUs on
one node exceeds ~100:

```bash
# Pure MPI (recommended):
#SBATCH --ntasks=26
#SBATCH --cpus-per-task=1

# Hybrid (use when many cores, few ranks makes sense):
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=4    # 7 x 4 = 28 CPUs total
```

---

## Running jobs

Always check free CPUs and set `--ntasks` in the script before submitting:

```bash
status       # shows free CPUs per node and best node recommendation
squeue --me  # your current jobs
```

### Step 1 - Generate surface muons

Edit `hpc/input_params.dat` (all parameters annotated with # comments).

```bash
cd ~/UCMuon
# Edit --ntasks in hpc/run_ucmuon_gen.sh
sbatch hpc/run_ucmuon_gen.sh hpc/input_params.dat
tail -f logs/ucmuon_gen_<JOBID>.out
```

Output in `output_<JOBID>/`:
- `ucmuon_selected.dat` -- detector-aimed muons (14 col)
- `ucmuon_selected_phits.dat` -- PHITS dump (auto-converted)
- `ucmuon_surface.dat` -- all muons (only if `save_all=1`)

### Step 2A - Transport with MUSIC (reference engine)

Edit `hpc/input_transport_music.dat`:
- **Line 1**: path to generator output, e.g. `output_6307xxx/ucmuon_selected.dat`
- **Line 10** (`init_tables`): `0` on first ever run, `1` for all subsequent runs

```bash
sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_music.dat
```

Output in `output_<JOBID>/`:
- `ucmuon_underground.dat` -- all muons with alive flag (18 col)
- `ucmuon_underground_phits.dat` -- PHITS dump, alive muons only (auto-converted)
- `input_used.dat` -- input copy for reproducibility

### Step 2B - Transport with Bethe-Bloch (cross-validation)

Use the same depth and density as Step 2A.

```bash
sbatch hpc/run_ucmuon_transport.sh hpc/input_transport_bb.dat
```

Both engines produce identical 18-column format. Compare survival rates directly.

---

## PHITS conversion

### Automatic (runs after every job)

Both SLURM scripts automatically call `ucmuon_to_phits` after merging.
No manual step needed.

### Why post-processing is better than in-loop PHITS writing

Writing to a PHITS file inside the generator hot loop forces all OMP threads
to serialise for every muon (`!$OMP CRITICAL` on each write). With 56 threads
and 5M muons this causes ~5M forced serialisations and reduces parallel
efficiency dramatically. The Fortran post-processor reads the merged file
once sequentially after the job -- no synchronisation overhead.

### Performance

| File size | Fortran converter | Python converter |
|---|---|---|
| 160 MB (1M muons) | ~3 s | ~20 s |
| 800 MB (5M muons) | ~15 s | ~100 s |
| 1.6 GB (10M muons) | ~30 s | ~200 s |

Always use the compiled Fortran converter for production files.

### Manual conversion

```bash
# Generator output -> PHITS (all muons in file):
./ucmuon_to_phits gen < output_6307xxx/ucmuon_selected.dat \
                      > output_6307xxx/ucmuon_selected_phits.dat

# Transport output -> PHITS (alive muons only):
./ucmuon_to_phits transport < output_6307xxx/ucmuon_underground.dat \
                            > output_6307xxx/ucmuon_underground_phits.dat
```

Progress is printed to stderr:
```
  ucmuon_to_phits: written        3436
  ucmuon_to_phits: skipped       96564
  Use in PHITS: s-type=17, dump=-10, 1 2 3 4 5 6 7 8 9 10
```

### PHITS output format (10 columns, D-exponent)

```
kf   x[cm]  y[cm]  z[cm]  u  v  w  Ekin[MeV]  weight  time[ns]
```

- `kf`: PDG code -- mu+ = -13, mu- = +13
- `u v w`: unit direction cosines
- `Ekin`: kinetic energy in MeV = (E_GeV - 0.10566) * 1000
- weight = 1.0, time = 0.0

### Using in a PHITS simulation

```
[ Source ]
  s-type = 17
  file   = output_6307xxx/ucmuon_underground_phits.dat
  dump   = -10
  1 2 3 4 5 6 7 8 9 10
```

For surface muons as source (before transport):
```
[ Source ]
  s-type = 17
  file   = output_6307xxx/ucmuon_selected_phits.dat
  dump   = -10
  1 2 3 4 5 6 7 8 9 10
```

---

## Output column formats

### Generator: ucmuon_selected.dat (14 columns)

```
EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV
         theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask
```

Coordinate convention: z=0 at surface, z<0 underground.
charge: +1=mu+, -1=mu-. hit_flag=1 means muon intersects a detector volume.

### Transport: ucmuon_underground.dat (18 columns)

```
EventID  x_srf_cm  y_srf_cm  z_srf_cm  E_srf_GeV  theta_srf  phi_srf  charge
         alive  x_ug_cm  y_ug_cm  z_ug_cm  E_ug_GeV  cx  cy  cz
         theta_ug  phi_ug
```

alive=1: muon reached the specified depth. alive=0: muon stopped in rock.
Only alive=1 rows are included in the PHITS dump.

---

## Expected wall times

Standard Rock, 90 m depth, 26 MPI ranks (pure MPI, one node):

| Task | N_muons | Wall time |
|---|---|---|
| Generator | 100k | ~1 min |
| Generator | 1M | ~7 min |
| Generator | 5M | ~35 min |
| MUSIC transport | 100k | ~1 min |
| MUSIC transport | 1M | ~8 min |
| BB transport | 100k | ~1 min |
| BB transport | 1M | ~12 min |
| PHITS conversion | 1M | ~3 s |
| PHITS conversion | 5M | ~15 s |

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `mpif90: not found` | MPI not loaded | `module load releases/2023b && module load foss/2023b` |
| `Cannot open file 'music-eloss-rock.dat'` | Tables not in root | Run with `init_tables=0` first |
| `No output files found` | Binary crashed | Check `logs/*.err` |
| Exit code 137 | SIGKILL from crash | Check `logs/*.err` for the real error |
| `Circular dependency dropped` | Make warning | Harmless, ignore |
| Some ranks idle | N_muons < nranks | Reduce `--ntasks` or increase `nmuons` |
| PHITS skipped 0, written 0 | File not found or all stopped | Check output_JOBID/ exists and has data |
