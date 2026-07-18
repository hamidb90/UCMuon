# MUSIC Source Files

The MUSIC muon transport code is **not included** in this repository.
It is distributed separately by its author and must be obtained directly.

## Files needed

Copy these two files into `src/transport/music/`:

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
cp music.f               src/transport/music/
cp music-crosssections.f src/transport/music/
bash setup.sh            # or: make music
```

## About music-eloss.dat

`music-eloss.dat` is **not a static file** — it is generated automatically
during the first run of the MUSIC transport driver for your specific
rock composition. It will appear in the project root after the first
successful MUSIC transport run.

## No MUSIC? Use these alternatives

| Engine | How to enable |
|--------|--------------|
| **Bethe-Bloch CSDA + MS** | Built automatically by `setup.sh` / `install.ps1` (needs gfortran only) |
| **PROPOSAL** | `pip install proposal` (system Python venv required — not miniforge) |
