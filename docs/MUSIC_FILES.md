# Enabling the MUSIC transport engine (Engine 2)

The MUSIC muon transport code is **not included** in this repository and is
**not** on Zenodo. It is the property of its author and is distributed only
on request. UCMuon works fully without it: Engines 1, 3, 4, 5, 6, and 7
build and run with no MUSIC files present. Only Engine 2 (MUSIC stochastic
MC) needs the files below.

## Step 1: Request the code from the author

Contact:

**Prof. Vitaly Kudryavtsev**
Department of Physics and Astronomy, University of Sheffield
v.kudryavtsev@sheffield.ac.uk

Ask for the MUSIC code for muography / muon transport research.

**Please cite:**
Kudryavtsev, V.A. (2009). *Muon simulation codes MUSIC and MUSUN for
underground physics.* Comput. Phys. Commun. 180, 339-346.
https://doi.org/10.1016/j.cpc.2008.10.013

## Step 2: What you will receive, and where to put it

The author sends a small set of files. Place each one as shown below
(paths are relative to the UCMuon project root):

| File you receive          | Copy it to                     | Purpose                                        |
|---------------------------|--------------------------------|------------------------------------------------|
| `music.f`                 | `src/transport/music/`         | MUSIC core transport engine                    |
| `music-crosssections.f`   | `src/transport/music/`         | Cross-section / energy-loss initialisation     |
| `music-double-diff-rock.dat`  | `data/`                    | Double-differential table, **rock** (required) |
| `music-double-diff-water.dat` | `data/`                    | Double-differential table, water (only for water / sea overburden) |

Notes:

- The two `.f` files go in `src/transport/music/` (next to the UCMuon
  wrapper files that are already there).
- The `music-double-diff-*.dat` tables go in `data/`. The build links
  them into `bin/` automatically.
- The author may also send standard CERN random-number routines
  (`ranlux`, `ranmar`, `corset`, `corgen`, `rnorml`). UCMuon already
  bundles equivalent CERN routines in `src/common/`, so you normally do
  **not** need to add these. Only add them if a build error reports one
  as missing.

## Step 3: Build

```bash
bash setup.sh          # Linux / macOS   (or: install.ps1 on Windows)
```

`setup.sh` detects `src/transport/music/music.f` and enables Engine 2
automatically. If the file is absent it simply skips MUSIC and builds the
other engines.

## Files that are generated for you (do not request these)

You do **not** receive these from the author. They are computed on the
first MUSIC run for your rock composition and written next to the binary:

- `music-eloss-rock.dat` (continuous energy-loss table)
- `music-cross-sections-rock.dat` (integral cross-sections)

This matches the author's instructions: cross-sections and energy losses
are *initialised (calculated)* for your medium from the double-differential
table, they are not shipped.

## No MUSIC? Use another transport engine

| Engine                    | How to enable                                              |
|---------------------------|-----------------------------------------------------------|
| **Bethe-Bloch CSDA + MS** | Built automatically by `setup.sh` / `install.ps1` (needs gfortran only) |
| **PROPOSAL**              | `pip install proposal` (system Python venv, not miniforge) |
| **PUMAS**                 | Auto-downloaded by `setup.sh`, or unzip the pumas repo to `external/pumas-master` |
