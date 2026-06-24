# UCMuon — directory map & roadmap

A single orientation page: what every folder is, and the ordered steps to a
publication-ready version. Updated 2026-06-11.

---

## 1. Directory map

```
UCMuon/
├── README.md              project overview + install instructions
├── ROADMAP.md             ← this file
├── CITATION.cff           how to cite
├── requirements.txt       Python dependencies
├── setup.sh / install.ps1 installers (auto-detect what to build)
├── run_gui.sh             launch the Streamlit GUI
├── tools/                 standalone CLI tools
│   └── compare_engines.py   run all engines on one source — physics audit + plots
│
├── src/                   Fortran/C engine sources
│   ├── generator/           surface muon generator (8 spectra)
│   ├── transport/           music/ · bethe_bloch/ · pumas/   (transport engines)
│   ├── parma/ · common/ · converters/
│
├── gui/                   Streamlit GUI (entry: ucmuon_gui.py) + Python drivers
│   └── ucmuon_stochastic_driver.py   ← UCMuon-MC (Engine ①, flagship)
├── data/                  physics tables (PARMA, MUSIC tables, …)
├── bin/                   compiled binaries (built by setup.sh; git-ignored)
├── docs/                  user docs (incl. MUSIC_FILES.md = how to obtain MUSIC)
├── examples/              worked examples (e.g. Vesuvius)
├── hpc/                   SLURM scripts + HPC input decks
│
├── manuscript/           the paper  ← see manuscript/PUBLICATION_TODO.md
│   ├── ucmuon_cpc_paper.tex / .pdf   + cas-sc template, .bib
│   ├── tab_survival_matrix.tex / tab_exitke.tex   generated validation tables
│   ├── scripts/             paper-asset generators (make_all_figs.sh,
│   │                        make_survival_table.py, make_fig06.py, …)
│   └── PUBLICATION_TODO.md       the granular pre-submission checklist
│
├── benchmark/            ← CURATED 6-code validation (codes + analysis + results)
│   ├── see benchmark/README.md   (per-code subfolders under benchmark/codes/)
│   └── geant4_muon_rock_v5/   ← git-ignored scratch tree (~8 GB raw .dat)
│       └── sources/benchmark_surface.dat   ← the authoritative v2 source file
│
├── external/             vendored / cloned third-party libraries
│   └── pumas-master/        PUMAS C library (used by the PUMAS engine / Makefile)
│
├── misc/                 loose regenerable scratch (see misc/README.md)
│
├── references/           literature PDFs + MURAVES slides (kept private for now)
│
└── (scratch — regenerable, kept on disk, NOT part of the clean repo)
    ├── output/                4.6 GB engine outputs
    ├── build/                 Fortran build artifacts
    └── *.dat, *.tif, ucmuon_autosave.json   loose regenerable files
```

**Rule of thumb:** everything above the "scratch" block is the real project;
everything in the scratch block is regenerable and can be deleted/ignored.

---

## 2. Where things stand (2026-06-11)

- **Code:** compiles and runs; GUI verified. **UCMuon-MC (Engine ①) now has the
  v2 physics**: per-process hard-event spectra (brems (1−v)/v, pair 1/v³,
  photonuclear 1/v), explicit δ-ray straggling, and the *deterministic-bound*
  pre-filter that removes the old mean-loss CSDA survival bias.
- **Manuscript:** compiles clean (42 pp, no undefined refs). Structure, novelty
  framing, density-inversion narrative, and the §8.4 "disabled scattering"
  rewording are done. **BUT the validation numbers (Table 1, survival matrix,
  exit-KE table, engine timings, quoted −1.2 % bias) are still from the
  2026-05-25 v2 benchmark run of the OLD engine** (v1 physics + mean-loss
  pre-filter). Two dated TODO comments in the .tex (lines ~952, ~1106) mark
  this; §8.4 currently *explains* a bias that §4.1 says no longer exists.
- **Benchmark:** curated `benchmark/` folder in place; raw 8 GB scratch in
  git-ignored `geant4_muon_rock_v5/`. Geant4/PHITS reference rows are verified
  and do NOT need re-running — only the four UCMuon engines do.
- **Repo:** private, not pushed. MUSIC + reference PDFs intentionally kept.

---

## 3. WHAT TO RUN NEXT (the critical path, in order)

The single decisive task is the **benchmark re-run with the new UCMuon-MC
physics**, then refreshing everything derived from it.

### Step 1 — re-run the four UCMuon engines on the v2 source
```bash
# all six depths; same source population the Geant4/PHITS rows used
for d in 1 25 50 100 200; do
  python3 tools/compare_engines.py \
      benchmark/geant4_muon_rock_v5/sources/benchmark_surface.dat  $d  2.65
done
# place/rename outputs as benchmark/geant4_muon_rock_v5/<ENGINE>/<ENGINE>_bench_<d>m.dat
# (ENGINE ∈ BB, UCMuon, MUSIC, PROPOSAL — see make_survival_table.py header)
```

### Step 2 — regenerate the paper tables
```bash
cd manuscript/scripts
python3 make_survival_table.py > ../tab_survival_matrix.tex
# tab_exitke.tex is hand-transcribed: recompute mean exit-KE per engine/depth
# from the new .dat files and update it + BENCHMARK_FEEDBACK.md §4.2
```

### Step 3 — regenerate figures and recompile
```bash
bash manuscript/scripts/make_all_figs.sh        # incl. fig06 survival curve
cd manuscript && latexmk -pdf ucmuon_cpc_paper
```

### Step 4 — update the prose that depends on the numbers
- §8.4: replace the "residual bias stems from the CSDA pre-filter" paragraph
  (now contradicts §4.1's deterministic-bound pre-filter) with the new measured
  bias; delete both `% TODO (2026-06-10)` comments.
- §9.4 (multi-engine cross-validation): still calls UCMuon-MC a
  "mean-energy-loss engine" and quotes the old 36–47 % spread — re-measure and
  reword (UCMuon-MC is now a full stochastic engine).
- Table `tab:engines` + Table 1 timing column: refresh wall-times.

### Step 5 — HPC scaling rows (can run in parallel with 1–4)
```bash
# on a CECI Lemaitre4 node:
bash manuscript/scripts/run_scaling_hpc.sh      # → scaling_hpc.csv
python3 manuscript/scripts/make_fig07_scaling.py
# then fill the k = 8/16/32 rows in tab:scaling (TODO at tex line ~1816)
```

---

## 4. Remaining roadmap (after the re-run)

Work top-to-bottom; details live in `manuscript/PUBLICATION_TODO.md`.

### Stage A — finish the science (decisive for acceptance)
1. **Benchmark re-run + table/prose refresh** — Section 3 above. ← **DO FIRST**
2. **Fix the PUMAS engine inconsistency.** The paper says "six engines" with
   Terrain as Engine 6 and never describes PUMAS, yet the intro (§1.4) claims
   "PROPOSAL and PUMAS as interchangeable cross-check engines" and the GUI
   ships "⑥ PUMAS" in the Transport selector (Terrain lives in Tab 3, not the
   engine list; `docs/ENGINE6_USAGE_GUIDE.md` calls Terrain Engine 6). Decide:
   either add a short PUMAS subsection + table row (7 engines), or drop the
   §1.4 PUMAS claim and renumber/relabel the GUI selector. A code-running CPC
   referee WILL notice.
3. **Add a "validity domain" paragraph** (§4.1 or Conclusions limitations):
   hard-event shapes are mean-conserving asymptotic spectra (not full
   Kelner/Kokoulin–Petrukhin — tail approximate beyond ~2 km.w.e., use the
   MUSIC engine there); no LPM suppression; δ-ray spectrum without the spin
   term; Highland Gaussian MS without Molière tails. State each as a bounded
   limitation with the in-package fallback. (Decided 2026-06-11: state, don't
   implement; the δ-ray spin factor is the only optional cheap fix.)
4. **Close the two open benchmark items**: PROPOSAL MCS X₀ fix (+15–19 % angle
   overestimate; correct X₀ ≈ 26.7 g cm⁻²) and the PHITS −12 % exit-KE
   single-layer test (monoenergetic 300 GeV, 200 m).
5. **Re-screenshot Fig. 4** (GUI) showing the new tab order and the
   "① ★ UCMuon-MC" engine selector
   (`manuscript/scripts/make_fig04_gui_screenshot.py`).
6. ~~Reword the §8.4 "disabled scattering" sentence~~ — **done** (now says
   survival is energy-loss dominated; Molière MCS doesn't alter it). Keep the
   guard: no sentence may claim MCS-*angle* agreement.
7. Rename the stale `\label{subsec:pumas}` on the UCMuon-MC subsection to
   `subsec:ucmuonmc` (cosmetic, but confusing given item 2).

### Stage B — journal & frame
8. **Target Computers & Geosciences** (decision made; CPC as hedge). The .tex
   still says `\journal{Computer Physics Communications}` — switch it (or
   consciously keep CPC) and re-check the abstract leads with the
   terrain→density-inversion workflow.

### Stage C — documentation & examples
9. Fill the Zenodo DOI placeholders in `README.md` + `CITATION.cff` once minted.
10. Add a second worked example (or a synthetic-inversion tutorial) beyond Vesuvius.
11. Verify `setup.sh` / `install.ps1` on a clean machine (the macOS PROPOSAL
    venv step is the fragile one).

### Stage D — public release (only when ready to submit)
12. Derive a clean public copy and **strip the restricted files**: MUSIC
    sources (`src/transport/music/music.f`, `music-crosssections.f`,
    `data/music-*.dat`) and copyrighted PDFs (`references/`). Keep the UCMuon
    MUSIC wrappers. Quick scripted step — see PUBLICATION_TODO.md.
13. Archive the tagged release on Zenodo; for CPC also prepare the Program
    Summary + Mendeley Data deposit, for C&G the software-availability
    statement.

---

## 5. Quick links
- Pre-submission checklist → `manuscript/PUBLICATION_TODO.md`
- Benchmark reproduction → `benchmark/README.md`
- Obtaining MUSIC → `docs/MUSIC_FILES.md`
