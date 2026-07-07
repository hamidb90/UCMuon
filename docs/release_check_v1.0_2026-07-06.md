# UCMuon v1.0 pre-release verification — 2026-07-06

Run against branch `fix/ucmuon-mc-decay-length` @ `ddf4cd4` (+ uncommitted
benchmark CSV refreshes). Expected values pinned to the 2026-07-05 baseline
(commit 021faac). Full sweep: 31 checks, findings collected without fixing.

## Results

| # | Check | Expected | Measured | Verdict |
|---|-------|----------|----------|---------|
| 1 | git clean | clean (manuscript artifacts excepted) | benchmark summary CSVs + run.log modified uncommitted; untracked: misc/, tools/, references/, CPC-guide HTML | **FAIL** (release hygiene) |
| 2 | clean build + MPI source fixes | 3 OMP binaries; b_rad_shape & Rayleigh in MPI src | all 3 built; MPI src has b_rad_shape (13×), Rayleigh sampler, no old cos-form; MPI link fails (known local gfortran/mpi.mod mismatch) | **PASS** (MPI build SKIP, cluster-only) |
| 3 | py_compile gui/*.py | pass | pass | **PASS** |
| 4 | mode-1 index | −3.195 ± 0.05 | MLE γ = 3.200; KS D=0.0017 < 0.0043; flux printed in cm⁻² | **PASS** |
| 5 | mode-2 endpoint | no spike; χ²/dof < 2 | top-1% z = −0.62σ; χ²/dof = 0.88 | **PASS** |
| 6 | PARMA | runs, no NaNs | N=100k, mean 7.44 GeV, 0 NaNs | **PASS** |
| 7 | modes 4–7 | run, falling spectra | all four: finite means, monotonic falls >10 GeV | **PASS** |
| 8 | angular modes | cos² KS; cone ⟨cosθ⟩=0.750±0.002; cos³ ±0.5% | KS D=0.0025 ✓; cone 0.7498 ✓; cos³ dev 0.04% ✓ | **PASS** |
| 9 | HEPEvt format | 8 fields, correct order | NHEP + ISTHEP IDHEP JDA1 JDA2 px py pz m; μ⁻→13; momenta exact | **PASS** |
| 10 | 100 m survival | 5 engines ± 0.4 pp | MUSIC 85.80 / UCMuon 85.82 / PROPOSAL 85.14 / py-BB 87.71 / f-BB 88.19; ǀUC−MUSICǀ=0.03 pp | **PASS** |
| 11 | 200 m exit E | 4 engines ± 0.6 GeV | 97.00 / 96.82 / 94.73 / 95.07; ǀpy−fǀ=0.34 GeV | **PASS** |
| 12 | MCS RMS 10 m/100 GeV | 0.110° ± 10 %; corr < 0.05 | f-BB 0.1101°, UCMuon-MC 0.1094°; corr −0.009 / −0.001 | **PASS** |
| 13 | 1.2 GeV stop depth | 228 ± 3 cm; InitKE 1.094 | both engines 228.2 cm; InitKE 1.0943 | **PASS** |
| 14 | E-column conventions | total for alive, 0 for stopped | all 5 engines OK; PROPOSAL: 46/85k survivors at E = m_μ exactly (KE=0 boundary crawlers, path-distance alive criterion) — cosmetic | **PASS** (note) |
| 15 | Gen tab + flux tracking | rate follows current run | run 156 ms; I_vert 0.00566→0.00046, band 58.6%→4.8% after E_min change | **PASS** |
| 16 | Transport tab 50 m | ≈ 98.1 %; Ice X₀ 36.08 | 98.03 % (10k); Ice → ρ 0.917, X₀ 36.08 g/cm², opacity consistent | **PASS** |
| 17 | Terrain tab | DEM loads, T=1 open sky, n_E-independent | tab renders, **graceful "rasterio not installed"** in the streamlit conda env; full pipeline verified headless in venv (open-sky T=1.000000 exact, ρ-monotone, n_E <1%) | **PARTIAL** (env, not code) |
| 18 | Results tab | φ-wrap clean; ρ changes curve | θ–φ map clean at ±180° (benign φ=0 bin for exactly-vertical muons); ρ 2.65→1.0 moves survival knee 155→420 m (×2.7 ≈ ρ ratio) | **PASS** |
| 19 | Density tab | renders | 5-map T_sim library loads: ρ∈[1.5,3.0], 360×85 | **PASS** |
| 20 | Config round-trip | restore w/o exception | tab renders, autosave live & matches session; upload untestable in automation harness (tool limit); restore fix verified at commit time (26f8ecb) | **PARTIAL** |
| 21 | Backward-MC panel | ≈1.14e-4 m⁻²s⁻¹, n_E-indep. | GUI 1.132e-4 @ n_E=80; CLI n_E 20→160 within 0.8% | **PASS** |
| 22 | CSG spot-check | 1.2 GeV stops 228 cm; E labeled total | verified numerically (228.2 cm exact vs PDG-2024; labels updated bce58f0); GUI drive skipped (needs example geometry setup) | **PARTIAL** |
| 23 | PROPOSAL custom medium | Sternheimer line, no fallback | API-level: custom Fe vs stock Iron dE/dx ≤1.1 %, custom rock vs StandardRock 0.2 % (2026-07-05); full driver run skipped (hours of PROPOSAL table building) | **PASS** (API-verified) |
| 24 | mono benchmark vs G4 | engines within Phase-2 tol. | v3 outputs (2026-07-05, this code): UCMuon 32.97 / MUSIC 32.92 / PROPOSAL 32.82 / BB 33.33 vs G4 32.72 % — all within tolerance | **PASS** |
| 25 | figures 9/10 visual | G4 ~100 % at 1 m; panels on-scale | fig06 loads refs from distilled summary ([SUMMARY] lines confirm); both PDFs visually verified post-021faac | **PASS** |
| 26 | paper builds | no errors/undefined | latexmk clean, 0 errors, 0 undefined refs; table numbers match regenerated CSVs (checked during 029484f) | **PASS** |
| 27 | CITATION.cff | version 1.0.0 | **version "0.9.0"**; DOI ✓ 10.5281/zenodo.20826984; title ✓ | **FAIL** (needs bump — release action) |
| 28 | README/quickstart | commands work; version strings | make ✓, streamlit run ✓, rasterio documented (§Install/§DEM); **no global app version string in GUI** | **PARTIAL** |
| 29 | PARMA carve-out | present | LICENSE §52 + THIRD_PARTY_LICENSES.md ✓ | **PASS** |
| 30 | public-repo delta | enumerate | ../UCMuon-public @ 21a40be (pre-July): needs full sync of all July fixes (src 2, gui 1+, benchmark analysis, examples/tsim, docs); private-only items to exclude: manuscript/, references/, misc/, tools/, benchmark run outputs, .claude/, autosaves | **FAIL** (sync pending — release action) |
| 31 | no hidden artifact deps | regenerable or shipped | T_sim library SHIPS (committed ddf4cd4) ✓; MUSIC tables self-generate (init=1) ✓; **misc/dem_site.tif is gitignored but make_tsim_library.py references it** — regeneration (not use) of the example library needs the DEM; README documents DEM download | **PARTIAL** |

## Known open items (status)

- ~~T_sim libraries~~ regenerated with fixed pipeline, ship in examples/vesuvius/ ✓
- ~~G4 angle comparison~~ closed 2026-07-05: UCMuon-MC/BB 0.90–1.03× G4 (BENCHMARK_FEEDBACK §8) ✓
- MPI binaries: cluster rebuild pending; fixed source verified in check 2.
- Raw mono Geant4 CSV pruned; distilled `figures_benchmark/benchmark_summary.csv` is authoritative → release-notes item.

## Verdict: **NO-GO — but only release mechanics remain**

Zero physics or code failures. All 26 substantive checks pass (3 partial for
environment/tooling reasons with the underlying physics verified numerically).

### Blockers (must do before tagging v1.0)
1. **Commit or discard the working tree** (benchmark CSVs, run.log; decide fate
   of untracked misc/, tools/, references/ — none of these ship).
2. **Version bump**: CITATION.cff → 1.0.0; add a visible app version string to
   the GUI header/footer.
3. **Public-repo sync**: UCMuon-public is at the pre-July state — every physics
   fix in this branch must be pushed; exclude the private dirs listed in #30.
4. **DEM availability**: either ship a DEM download script/pointer next to
   `examples/vesuvius/make_tsim_library.py` or document that regeneration
   requires the SRTM tile (README §"Get a DEM file" covers the how).

### Non-blockers (release notes)
- MPI binaries require cluster rebuild (source is fixed).
- PROPOSAL boundary muons written alive with KE=0 (0.05 %, cosmetic).
- Terrain tab requires `pip install rasterio` (documented; error is graceful).
- Raw mono Geant4 event CSV pruned — distilled summary is the reference.
- Streamlit env for full GUI needs rasterio if Terrain tab is used.
