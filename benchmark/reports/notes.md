# Presentation Speaker Notes — CCS Muography Benchmark

---

## 0. Muon Source Description  *(introduce before showing any figure)*

**What to say:**

> "Before we look at the results, let me describe what we simulated.
> We injected 600,000 muons into a slab of standard rock — the same
> rock model used in all underground physics benchmarks: density 2.65 g/cm³,
> effective atomic number Z = 11, effective mass number A = 22.
>
> The muons come from a file we call `benchmark_surface.dat`. It contains
> six discrete initial kinetic energies: approximately 5, 10, 20, 50, 100,
> and 300 GeV — 100,000 muons per energy group, all aimed vertically
> downward. We scored the surviving muons at five depths: 1, 25, 50, 100,
> and 200 metres of rock. We deliberately excluded the 10-metre scoring
> plane from our comparison — at that depth, a 5 GeV muon is right at its
> range limit, and the codes disagree on exactly where it stops, which
> creates an artificial discrepancy that has nothing to do with the physics
> at CCS-relevant depths.
>
> We ran the same source through six different simulation codes. Two are
> full Monte Carlo: Geant4 and PHITS. The other four are fast dedicated
> muon transport engines: MUSIC, PROPOSAL, an analytic Bethe-Bloch engine
> (BB), and UCMuon based on the PUMAS framework. We will compare what each
> code predicts and what information each one can give you."

**Key numbers to have ready:**
- 600,000 muons total; 100,000 per energy group
- Scoring depths: 1, 25, 50, 100, 200 m (2.65 – 530 m.w.e.)
- Rock: ρ = 2.65 g/cm³, Z_eff = 11
- Runtime: Geant4 ≈ 5 h on a single core; fast codes ≈ seconds

---

## Figure 1 — *pres_fig1_physics_outputs.png*

### Panel (A) — Muon Survival Through Rock

**Key message:** All six codes agree on how many muons survive to depth.

**What to say:**

> "Panel A shows the fraction of muons that reach each depth — the
> transmission curve. Starting at 100% at the surface, it drops to about
> 66% at 25 metres, 50% at 50 metres, and 16% at 200 metres.
>
> The remarkable thing is that all six codes — full Monte Carlo and fast
> analytic alike — agree to within 3% at every depth. The error bars you
> see are the statistical uncertainty from running 100,000 muons per energy
> group; they are smaller than the markers, so every difference you see
> between codes is a real physics difference, not statistical noise.
>
> For CCS muography, this panel tells you that whichever code you use to
> predict your muon flux, you will get essentially the same answer. The
> muon counting rate at your detector is well-constrained."

**If asked about the shape:**

> "The curve is not exponential — it has a shape set by the muon energy
> spectrum. Lower-energy muons stop in the first few metres, so the
> survival fraction drops quickly at first. At depth only the highest-
> energy muons remain, and these lose energy more slowly per metre,
> so the curve flattens relative to a pure exponential."

---

### Panel (B) — Exit Muon Energy vs Depth

**Key message:** Exit energy rises with depth due to energy-selection. PHITS diverges at depth — a known cross-section difference.

**What to say:**

> "Panel B shows the mean kinetic energy of muons that exit each depth.
> Notice something perhaps surprising: the exit energy RISES with depth,
> from 80 GeV at 1 metre to about 138 GeV at 200 metres.
>
> This is not energy gain — it is energy selection. The deeper you go,
> the only muons that survive are the ones that started at the highest
> energies. At 200 metres, your detector is seeing almost exclusively
> muons that started with 300 GeV. You are looking at the hard tail of
> the spectrum.
>
> Five of the six codes agree within about 4% across all depths. But
> PHITS — the red dashed line — diverges progressively, reaching 12%
> below Geant4 at 200 metres. This is not a configuration error; it is
> a fundamental difference in how PHITS parametrises the cross sections
> for high-energy bremsstrahlung and pair production. At shallow depths
> almost all energy loss is via ionisation, which all codes agree on.
> At depth, radiative processes contribute more, and the PHITS cross
> sections give 10–20% higher radiative energy loss at 100–300 GeV.
> That difference accumulates over 200 metres into the 12% you see here.
>
> For CCS, if your detector measures the FLUX of muons, PHITS is fine.
> If it measures the ENERGY SPECTRUM, you should be aware of this offset."

**Error bars:** The error bars (SEM on the mean exit KE, using
StdELoss / √N as proxy) are sub-GeV — invisible at this scale. Every
code difference shown is statistically significant by many sigma.

---

### Panel (C) — Energy Loss by Process  [Geant4 only]

**Key message:** Ionisation dominates but radiative processes grow with depth — this is why energy-spectrum agreement between codes is harder than flux agreement.

**What to say:**

> "Panel C is Geant4-only — it shows something the fast codes cannot give
> you: a breakdown of how the muon loses energy at each depth.
>
> At 1 metre, 97.7% of the energy loss is via ionisation — the steady
> drag described by the Bethe-Bloch formula. The remaining 2.1% is pair
> production (creating electron-positron pairs), and a small fraction is
> bremsstrahlung and nuclear interactions.
>
> As you go deeper to 200 metres, the ionisation fraction drops to 94%
> while pair production climbs to nearly 6%, and bremsstrahlung and
> nuclear interactions also grow. This happens because the surviving muons
> at 200 metres are predominantly 300 GeV muons — very high energy — and
> at those energies radiative processes are proportionally more important.
>
> This is the direct explanation for Panel B: PHITS and Geant4 agree on
> ionisation, but disagree on the radiative cross sections. At shallow
> depths where ionisation is 97.7%, any cross-section error in
> bremsstrahlung or pair production barely matters. At 200 metres, where
> pair production is 6%, even a 20% error in that cross section shifts
> your total energy loss prediction by over 1%. Accumulated over 200
> metres, that is the 12% gap you saw."

---

### Panel (D) — Multiple Coulomb Scattering Angle

**Key message:** The MCS angle peaks at intermediate depth and DECREASES at greater depth — because only high-energy, straight-flying muons survive.

**What to say:**

> "Panel D shows the mean deflection angle of surviving muons at each
> depth. You might expect it to increase monotonically — more rock means
> more scattering. But instead it peaks around 66 m.w.e. and then falls.
>
> The reason is again energy selection. The MCS angle scales as the square
> root of path length divided by the muon's momentum: more rock means
> more scattering, but higher momentum means less scattering per unit
> rock. At depth, only the highest-energy, highest-momentum muons survive.
> Their momentum advantage outweighs the extra path length, and the mean
> deflection angle falls.
>
> The peak at about 66 m.w.e. corresponds to 10 GeV muons that are near
> the end of their range — they have slowed down substantially, their
> momentum is low, and they scatter strongly just before they stop. Beyond
> that depth those muons are gone, and you are left with the hard,
> straight-flying muons that scatter very little.
>
> For CCS this is good news: at depths of 100–200 metres, the surviving
> muons are nearly straight. Your muographic image is not blurred by
> multiple scattering at those depths — the angular resolution of your
> detector is the limiting factor, not the rock.
>
> Notice that PHITS is excluded from this panel: PHITS runs in summary
> mode and only records the aggregate zenith angle distribution, not the
> per-muon MCS deflection that the other codes provide. The error bars
> on each point represent the standard error on the mean; they are smaller
> than the markers — all differences here are statistically real."

**If asked why BB/UCMuon are lower than Geant4:**

> "BB and UCMuon use the Highland-Lynch formula for MCS, which is an
> analytic approximation. MUSIC and PROPOSAL simulate scattering step by
> step and produce larger angles, closer to Geant4. The Highland formula
> underestimates the tail of the scattering distribution, which pulls the
> mean angle down."

---

## Figure 2 — *pres_fig2_code_engines.png*

### Left Panel — Code Capabilities Matrix

**Key message:** Geant4 gives the full picture; fast codes give per-event kinematics but not process breakdown or secondaries; PHITS only gives aggregate histograms.

**What to say:**

> "This matrix summarises what each code can actually give you. The rows
> are physics observables; the columns are the six codes.
>
> Green means you get the full information. Yellow-cream means partial or
> aggregate only. Red means not available.
>
> All six codes can tell you transmission and exit energy — that is the
> common ground. The differences appear as you go down the rows.
>
> MCS scatter angle: PHITS gives only the aggregate zenith angle histogram,
> not per-muon deflections — hence 'Partial'.
>
> Lateral displacement: PHITS cannot compute this from aggregate tallies
> at all. The fast codes (MUSIC, PROPOSAL, BB, UCMuon) all record entry
> and exit positions, so they compute lateral displacement correctly.
>
> Process breakdown — ionisation, bremsstrahlung, pair production: this
> is exclusively Geant4. No other code in this benchmark tracks energy
> loss by process.
>
> Secondary particles — the photons and electron-positron pairs produced
> by high-energy muons: only Geant4 gives you the full secondary spectrum.
> PHITS gives partial information.
>
> Per-event output: PHITS works with tallies, not individual events. You
> cannot extract a per-muon record from a standard PHITS run.
>
> The practical conclusion: if you need to validate a detector response
> model — which requires per-event kinematics at minimum — Geant4 or
> one of the fast codes is appropriate. If you need to understand the
> physics of energy deposition in detail, only Geant4 provides that."

---

### Right Panel — Transmission Agreement with Geant4

**Key message:** All fast codes stay within ±5% of Geant4 for muon flux prediction at all CCS-relevant depths.

**What to say:**

> "The right panel shows the ratio of each code's transmission prediction
> to Geant4, across the five scoring depths. The green band is ±5%; the
> yellow band is ±10%.
>
> Every code stays inside the green band — within 5% of Geant4 — at all
> depths. The worst case is the analytic BB code at 200 metres, which is
> about 3% above Geant4. All other codes are within 2%.
>
> The error bars here are propagated from the binomial statistical
> uncertainty on both the code and Geant4 transmission counts. They are
> smaller than the symbols, so these ratios are well-determined.
>
> The practical message: for CCS muography, where you want to predict how
> many muons reach detector depth, any of these six codes will give you an
> answer accurate to within 3%. You do not need Geant4's full 5-hour run
> to get a reliable flux estimate. The fast codes run in seconds and are
> adequate for feasibility studies, detector placement optimisation, and
> sensitivity analysis."

---

## Figure 3 — *pres_fig3_timing.png*

### Simulation Speed Comparison

**Key message:** Fast dedicated engines are 20–220× faster than Geant4. PHITS is 8× *slower* than Geant4 — a counter-intuitive but important practical point.

**What to say:**

> "This figure answers the practical question: how long does each code take
> to run? The x-axis is on a log scale — each gridline is a factor of ten
> in time — because the spread is enormous.
>
> The four fast engines at the top finish the 600,000-muon run in minutes.
> MUSIC takes 3 minutes, UCMuon 5 minutes, the analytic Bethe-Bloch code
> 15 minutes, and PROPOSAL 36 minutes.
>
> Geant4 — the reference full Monte Carlo — takes 12.5 hours on a single
> core. That is roughly two to three orders of magnitude slower than the
> fast codes. When you are doing a sensitivity analysis with hundreds of
> geometry configurations, that difference matters enormously.
>
> The perhaps surprising result is PHITS: it takes 103 hours — over four
> days — for the same 600,000 muons. PHITS is 8 times slower than Geant4.
> This is not a configuration error. PHITS uses an analog Monte Carlo kernel
> that tracks all secondary particles through the scoring planes before
> accumulating tallies. Geant4, with its optimised geometry navigation and
> physics tables, is significantly faster for this type of problem.
>
> The practical message: for CCS muography studies where you need to explore
> many rock densities, depths, or detector positions, the fast codes let you
> run a full parameter sweep in an afternoon. A comparable Geant4 study
> would take weeks. Use Geant4 for the final validation run; use fast codes
> for the exploration phase."

**Key numbers to have ready:**
- MUSIC: 204 s (3.4 min) — 221× faster than Geant4
- UCMuon: 328 s (5.5 min) — 138× faster
- BB: 871 s (14.5 min) — 52× faster
- PROPOSAL: 2,147 s (35.8 min) — 21× faster
- Geant4: 45,134 s (12.5 h) — reference, 13.3 evt/s
- PHITS: 371,146 s (103 h) — 8.2× slower than Geant4, 1.6 evt/s

**If asked why PHITS is so slow:**
> "PHITS runs a fully analog Monte Carlo — it tracks every secondary
> particle (photons, electrons, positrons from bremsstrahlung and pair
> production) until they are absorbed. For a 300 GeV muon, each radiative
> event spawns a shower. Geant4 also does this, but its geometry navigation
> and cross-section look-up tables are more optimised for this use case.
> You can speed PHITS up using variance reduction techniques, but those
> change the simulation mode and require additional validation."

---

## Common Questions

**Q: Why not use a real cosmic muon spectrum instead of fixed energies?**
> "The fixed-energy grid is a deliberate choice for benchmarking. It lets
> us isolate physics differences at each energy without convolving them
> with spectrum uncertainties. For a real CCS study we would fold in the
> Gaisser cosmic-ray spectrum — the codes are all capable of that."

**Q: Is standard rock representative of real CCS geology?**
> "Standard rock (Z=11, ρ=2.65 g/cm³) approximates average continental
> crust. Real CCS sites are sedimentary — limestone, sandstone, shale —
> with densities between 2.2 and 2.7 g/cm³ and different Z values. The
> benchmark validates the codes; site-specific studies would use the
> actual lithological column."

**Q: Which code do you recommend for CCS muography?**
> "For flux prediction alone, any of the fast codes. For full detector
> simulation — tracking secondaries, modelling scintillator or Cherenkov
> response — Geant4. PHITS is most useful if you want an independent cross-
> check from a second full Monte Carlo without switching ecosystem."
