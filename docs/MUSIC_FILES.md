# MUSIC Source Files

The MUSIC muon transport code is **not included** in this repository.
It is distributed separately by its author and must be obtained directly.

## Files needed

Copy these two files into this directory (`src/music/`):

| File | Description |
|------|-------------|
| `music.f` | MUSIC core transport engine |
| `music-crosssections.f` | Cross-section calculation routines |

## How to obtain

Contact the author:

**Prof. Vitaly Kudryavtsev**  
Department of Physics & Astronomy, University of Sheffield  
✉️  v.kudryavtsev@sheffield.ac.uk

Request the MUSIC code for muography research.

**Reference (please cite):**  
Kudryavtsev, V.A. (2009). *Muon simulation codes MUSIC and MUSUN for underground physics.*  
Comput. Phys. Commun. 180, 339–346. https://doi.org/10.1016/j.cpc.2008.10.013

## After receiving the files

```bash
cp music.f               src/music/
cp music-crosssections.f src/music/
bash setup.sh            # or: make music
```

## About music-eloss.dat

`music-eloss.dat` is **not a static file** — it is generated automatically
during the first run of `cosmoaleph_music_driver_omp` for your specific
rock composition. It will appear in the project root after the first
successful MUSIC transport run.

## No MUSIC? Use these alternatives

| Engine | How to enable |
|--------|--------------|
| **Bethe-Bloch + MS** | Place `cosmoaleph_phitsxs_omp.f90` in `src/bethe_bloch/`, run `make bb` |
| **PROPOSAL** | `pip install proposal` (system Python venv required — not miniforge) |
