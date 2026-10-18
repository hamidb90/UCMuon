# Third-party components and licences

Everything in this repository is UCMuon's own code and data, under the MIT
License (see `LICENSE`), apart from the public-domain elevation model listed at
the end. The components below are **not** in this repository: the installers
download them, or you obtain them from their authors, and they keep their own
terms.

---

## Downloaded by the installers

### PARMA 4.10 / EXPACS (JAEA): spectrum mode 3, optional

`tools/install_parma.py`, which `setup.sh` and `install.ps1` run (the default
answer is yes), downloads JAEA's `parma.zip` and `parma_cpp.zip` from
https://phits.jaea.go.jp/expacs/, checks the two code files against the
release UCMuon was tested with (`tools/parma_v4.10.sha256`), and writes three
git-ignored items:

- `data/EXPACS/parma/`: JAEA's package, unpacked (the data tables);
- `src/parma/parma_subroutines.f90`: JAEA's Fortran with UCMuon's edits
  (configurable data directory; date lookup limited to each table's years),
  listed in its header and in the script;
- `ucmugen/include/UCMuGen_PARMA.h`: generated from JAEA's C++ by
  `ucmugen/tools/make_parma_header.py`.

Terms (JAEA's conditions for use):

- Copyright (c) 2006 Japan Atomic Energy Agency (JAEA); PARMA was developed by
  Dr. Tatsuhiko Sato.
- Free for **non-commercial** use, redistribution and modification. Commercial
  use requires a prior agreement with JAEA.
- Publications using it must cite T. Sato, *PLoS ONE* **10**(12): e0144679
  (2015), doi:10.1371/journal.pone.0144679, and T. Sato, *PLoS ONE* **11**(8):
  e0160390 (2016), doi:10.1371/journal.pone.0160390, acknowledge
  https://phits.jaea.go.jp/expacs, and be sent to nsed-expacs@jaea.go.jp.
- JAEA's package also contains a Mersenne Twister (M. Matsumoto and
  T. Nishimura; Fortran translation by H. Takano) under the GNU LGPL, version 2
  or later, in its example program `main-generator.f90`, which UCMuon does not
  use.

A build that includes PARMA (a generator linked with `parma_subroutines.f90`,
or C++ code that includes `UCMuGen_PARMA.h`) is for non-commercial use only.
Without PARMA every other spectrum and engine works, the generator's mode 3
stops with a message saying how to install it, and `UCMuGen.h` stays MIT.

### PUMAS (V. Niess): Engine 7, optional

Downloaded into `external/pumas-master/` (a pinned commit) by the installers.
LGPL-3.0: see `LICENSE` and `COPYING.LESSER` there. Reference: V. Niess,
*Comput. Phys. Commun.* **279** (2022) 108438.

### Python packages

Installed with pip from `requirements.txt` (Streamlit, pandas, NumPy, Plotly,
SciPy, Matplotlib), and optionally rasterio (Engine 6) and PROPOSAL (Engine 4,
LGPL-3.0), each under its own licence.

---

## Obtained from their authors

- **MUSIC** (V. A. Kudryavtsev): Engine 2. Not redistributable; UCMuon ships
  only the wrapper. See `docs/MUSIC_FILES.md`.
- **EcoMug** (D. Pagano, GPL-3.0): used only by the optional Geant4 benchmark
  application; download it separately
  (`benchmark/codes/geant4/include/EcoMug_NOTICE.txt`).
- **Geant4** and **PHITS**: used only to run the benchmarks; not included.

---

## Random number routines

`src/common/ucmuon_rng.f90` is UCMuon's own code (MIT): RANLUX (Lüscher 1994;
James 1994), RANMAR (Marsaglia, Zaman and Tsang 1990; James 1990) and Leva's
normal generator (Leva 1992), written from the published descriptions, with
one state per OpenMP thread. `tests/rng/` checks that it reproduces, number for
number, the CERN Program Library routines that UCMuon bundled up to 1.4.1.

---

## Licence notice for releases 0.9.0 to 1.4.1

Besides UCMuon's MIT code, those releases contained:

- PARMA/EXPACS (`src/parma/parma_subroutines.f90`, `data/EXPACS/`,
  `ucmugen/include/UCMuGen_PARMA.h`), under JAEA's terms above, as this file
  then stated;
- CERN Program Library routines, distributed by CERN under the GNU General
  Public License, which this file did not list: `src/common/ranlux.f` and
  `ranlux_omp.f` (RANLUX, Mathlib V115, F. James), `ranmar_omp.f` (RANMAR,
  Kernlib V113), `rnorml.f` (RNORML, Mathlib) and `corset.f` (Mathlib, unused);
  `ranlux_omp.f` and `ranmar_omp.f` were modified by UCMuon to keep one
  generator state per OpenMP thread;
- the LGPL Mersenne Twister above, inside `data/EXPACS/parma/`, also not
  listed.

In those releases the MIT licence covers UCMuon's own code only; the files
above stay under their own licences.

---

## Bundled data

- `examples/vesuvius/vesuvius_dem.tif`: SRTM GL1 (NASA/USGS), public domain;
  see `examples/vesuvius/DEM_SOURCE.md`.
