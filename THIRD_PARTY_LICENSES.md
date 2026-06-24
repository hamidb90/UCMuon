# Third-party components and licenses

UCMuon's own source code is released under the MIT License (see `LICENSE`).
However, this repository bundles third-party components that are **NOT** covered
by the MIT License and carry their own terms. By using those components you
agree to their respective conditions.

---

## PARMA / EXPACS  —  `src/parma/`, `data/EXPACS/`

Spectrum mode 3 ("PARMA / EXPACS") uses the PARMA v4.10 model and its data
tables developed by **Dr. Tatsuhiko Sato, Japan Atomic Energy Agency (JAEA)**.

- Copyright (c) Japan Atomic Energy Agency (JAEA).
- **Non-commercial use only.** Commercial use is NOT allowed without a prior
  agreement with JAEA.
- Any published use **must cite**:
  - T. Sato, *PLoS ONE* **10**(12): e0144679 (2015). doi:10.1371/journal.pone.0144679
  - T. Sato, *PLoS ONE* **11**(8): e0160390 (2016). doi:10.1371/journal.pone.0160390
  - and acknowledge the URL https://phits.jaea.go.jp/expacs
- Contact: nsed-expacs@jaea.go.jp
- PARMA bundles a Mersenne Twister RNG by Prof. M. Matsumoto (Hiroshima Univ.).

Full conditions: `data/EXPACS/EXPACS_CONDITIONS_FOR_USE.txt`.
The UCMuon copy is unmodified physics; only the data-directory path was made
configurable (see header of `src/parma/parma_subroutines.f90`).

> **Implication:** because PARMA/EXPACS is non-commercial-only, the repository
> as a whole cannot be used commercially without removing `src/parma/` +
> `data/EXPACS/` (and obtaining JAEA agreement). All other engines are MIT.

---

## PUMAS  —  `external/pumas-master/`

The PUMAS muon/tau transport library by **V. Niess** is licensed **LGPL-3.0**
(see `external/pumas-master/LICENSE` and `COPYING.LESSER`). Reference:
V. Niess, *Comput. Phys. Commun.* **279** (2022) 108438.

---

## NOT included in this repository

- **MUSIC** (V. A. Kudryavtsev, University of Sheffield) — not redistributable.
  UCMuon ships only the wrapper. To obtain MUSIC, see `docs/MUSIC_FILES.md`.
- **EcoMug** (D. Pagano, GPL-3.0) — used only by the optional Geant4 benchmark
  app; download it separately
  (`benchmark/codes/geant4/include/EcoMug_NOTICE.txt`).
