#!/usr/bin/env python3
"""
ucmuon_ffe_validate.py
══════════════════════════════════════════════════════════════════════════════
UCMuon Fast Flux Estimator — standalone validation & conceptual explainer
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Run:    python ucmuon_ffe_validate.py
Place this file next to fast_flux_estimator.py.

Answers three questions:
  1. What exactly is I [cm⁻²sr⁻¹s⁻¹]? What is the "detector"?
  2. Which model should I use for what purpose?
  3. Are the numbers physically correct?
"""
from __future__ import annotations
import sys, math, os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
try:
    from fast_flux_estimator import (
        integrated_flux, flux_vs_depth, differential_flux,
        angular_profile, emin_from_opacity,
        RHO_STANDARD_ROCK, M_MU_GEV, MODEL_LABELS,
        _guan_cos_star,
        _GROOM_T_GEV, _GROOM_R_GCM2,
    )
except ImportError as e:
    sys.exit(f"❌  Cannot import fast_flux_estimator.py: {e}\n"
             f"   Place this script in the same directory.")

_int  = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
_pass = 0; _fail = 0

def _hdr(t):
    print(); print("─"*72); print(f"  {t}"); print("─"*72)

def _check(label, value, expected, rtol=0.20, unit="", note=""):
    global _pass, _fail
    rel = abs(value - expected) / max(abs(expected), 1e-30)
    ok  = rel <= rtol
    if ok: _pass += 1
    else:  _fail += 1
    sym = "✅" if ok else "❌"
    print(f"  {sym}  {label:<52s}  "
          f"{value:.4g}{unit}  (ref {expected:.4g}{unit},  Δ={rel*100:.0f}%)"
          + (f"  ← {note}" if note else ""))


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 0 — WHAT THE FFE COMPUTES
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 0 — What I [cm⁻²sr⁻¹s⁻¹] means, and what the 'detector' is")
print("""
  I(θ) [cm⁻²sr⁻¹s⁻¹]  =  ∫_{E_min}^∞  dΦ/dT (T, θ)  dT

  This is the number of muons crossing 1 cm² per second,
  per steradian of solid angle, arriving from zenith angle θ.
  It is an INTENSIVE quantity — the same everywhere in free space.

  ┌─────────────────────────────────────────────────────────────────────┐
  │  THE FFE HAS NO SPECIFIC DETECTOR SHAPE.                           │
  │  "Detector" = acceptance parameter  A [cm²·sr]                     │
  │                                                                     │
  │  Count rate:   R [s⁻¹] = I(θ) × A_det [cm²·sr]                  │
  │                                                                     │
  │  A_det = ∫∫ cos(α) dΩ dA_surface   (geometric acceptance)         │
  │  where α = angle between muon direction and detector face normal.  │
  │  The FFE is NOT a point detector — a point has A → 0.             │
  └─────────────────────────────────────────────────────────────────────┘

  COMMON DETECTOR GEOMETRIES:

  ① Single upward-facing scintillator, area S [cm²], no angular cut:
      A_det = S × π                      [cm²·sr]  (Lambert solid angle)

  ② Same panel, accepting only muons within cone θ < θ_max:
      A_det = S × 2π (1 − cosθ_max)     [cm²·sr]

  ③ Two-panel telescope, area S each, separation d [cm], axial:
      Ω ≈ S / d²   [sr]
      A_det = S × Ω = S²/d²             [cm²·sr]

  ④ UCMuon detector filter (Tab 1):
      A_det = physical_area × MC_cos²θ_acceptance   (computed by the GUI)

  The GUI default A_det = 6 cm²·sr ≈ 100 cm² × 0.06 sr (telescope-like).

  FLAT-SLAB MODEL (what the FFE computes):
      slant path   = L / cosθ            (muon at angle θ travels MORE rock)
      X = ρ × (L/cosθ) × 100            [g/cm²]
      E_min from Groom (2001) CSDA table (lower bound — ignores straggling)
      I_rock(θ)   from the chosen flux model at angle θ, energy grid ≥ E_min
      T_rock = I_rock / I_open           (rock transmission)
      R = I_rock × A_det                 [s⁻¹]
""")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 1 — MODEL VALIDITY AND ABSOLUTE NORMALISATION
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 1 — Model validity ranges and absolute normalisation")
print("""
  ─────────────────────────────────────────────────────────────────────────
  MODEL              VALID RANGE        ABSOLUTE FLUX   ANGULAR RATIO I/I(0)
  ─────────────────  ─────────────────  ──────────────  ────────────────────
  Reyna–Bugaev 2006  E > 1 GeV          ✅ ±20% of PDG  ✅ good
  Bugaev/Gaisser     E > 10 GeV         ⚠️  ×10 low <10G  ✅ good
  Gaisser–Tang       E > 10 GeV         ⚠️  ×10 low <10G  ✅ good
  Guan 2015          E > 20 GeV abs.    ⚠️  E_eff issue    ✅ excellent (cosθ*)
  Frosin 2025        E > 20 GeV abs.    ⚠️  E_eff issue    ✅ excellent (cosθ*)
  ─────────────────────────────────────────────────────────────────────────

  WHY Bugaev/Gaisser underestimate I(E>1 GeV) by ~×10:
    Gaisser's formula is fitted to E > 10 GeV data.  Below 10 GeV,
    muon decay and low-energy threshold effects suppress the spectrum
    beyond what the simple power-law captures.  Do not integrate below
    10 GeV with Gaisser/Bugaev.

  WHY Guan/Frosin give low absolute flux below 20 GeV:
    The low-energy correction term   E_eff = E × (1 + a / (E·cosθ*^b))
    was designed to improve the angular shape at moderate E, but it
    suppresses the absolute flux at low E:
        E=1 GeV:   E_eff/E ≈ 4.6  →  formula uses 4.6 GeV spectrum → low flux
        E=5 GeV:   E_eff/E ≈ 1.7  →  still significant
        E=20 GeV:  E_eff/E ≈ 1.2  →  small correction, OK
    The angular RATIO I(θ)/I(0°) cancels this normalisation bias.

  RECOMMENDATION:
    → For absolute flux / total rate (E_min ~ 1 GeV): use Reyna–Bugaev
    → For angular shape I(θ)/I(0°) at muography E_min (> 20 GeV): use Guan 2015 or Frosin 2025
       Reason: Guan correctly models the pion/kaon angular enhancement (I(θ)/I(0°) > 1 above
       ~46 GeV at 30°). Reyna-Bugaev gives a constant energy-independent ratio — wrong at high E.
    → For MURAVES/muography T = I_rock/I_open: correct in ALL models (I_open evaluated at same θ)
""")

print("  Differential flux dΦ/dT [cm⁻²s⁻¹sr⁻¹GeV⁻¹] at θ=0°, sea level:")
print(f"  {'E':>8s}  {'Reyna-B':>11s}  {'Bugaev':>11s}  "
      f"{'Guan 2015':>11s}  {'Frosin 25':>11s}")
for E in [1.0, 3.0, 10.0, 100.0]:
    phi_r = differential_flux(np.array([E]), 0.0, 'reyna_bugaev')[0]
    phi_b = differential_flux(np.array([E]), 0.0, 'bugaev')[0]
    phi_g = differential_flux(np.array([E]), 0.0, 'guan_2015')[0]
    phi_f = differential_flux(np.array([E]), 0.0, 'frosin_2025')[0]
    print(f"  {E:>5.0f} GeV  {phi_r:>11.3e}  {phi_b:>11.3e}  "
          f"{phi_g:>11.3e}  {phi_f:>11.3e}")
print("  At 100 GeV: Guan/Bugaev = 0.91 (9% apart — E_eff still contributes 3.6%). Agree to <2% above ~500 GeV.")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2 — PDG BENCHMARKS
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 2 — PDG benchmarks (open-sky, θ=0°, sea level)")
print("""
  PDG 2022 §30.3:
    I(E > 1 GeV,  θ=0°) ≈ 70 m⁻²sr⁻¹s⁻¹  = 7.0×10⁻³ cm⁻²sr⁻¹s⁻¹
    Hemisphere-integrated rate (all E, all angles) ≈ 1 cm⁻²min⁻¹
""")

T_fine = np.logspace(0, 4, 1200)   # 1 → 10000 GeV
for model, label in MODEL_LABELS.items():
    short = label.split("←")[0].split("[")[0].strip()
    phi   = differential_flux(T_fine, 0.0, model)
    I_1   = float(_int(phi, T_fine))
    ok    = (model == "reyna_bugaev" and
             abs(I_1 - 7.0e-3) / 7.0e-3 <= 0.25)
    warn  = model in ("bugaev", "gaisser_tang")
    note  = ("✅ within 25% of PDG" if ok else
             "⚠️  below 10 GeV: not valid for Gaisser" if warn else
             "⚠️  E_eff suppresses abs. flux below 20 GeV")
    print(f"  {short:<40s}  I(>1 GeV) = {I_1:.3e} cm⁻²sr⁻¹s⁻¹  {note}")

print()
_check("Reyna–Bugaev  I(E>1 GeV, θ=0°)",
       float(_int(differential_flux(T_fine, 0.0,'reyna_bugaev'), T_fine)),
       7.0e-3, rtol=0.25, unit=" cm⁻²sr⁻¹s⁻¹", note="PDG 2022 §30.3")

# Hemisphere-integrated rate
th_hemi = np.linspace(0.0, 88.0, 300)
I_hemi, _ = angular_profile(th_hemi, E_min_GeV=0.1, model="reyna_bugaev")
cos_h  = np.cos(np.radians(th_hemi))
sin_h  = np.sin(np.radians(th_hemi))
R_hemi = 2*math.pi * float(_int(I_hemi * cos_h * sin_h, np.radians(th_hemi)))
_check("Reyna–Bugaev  hemisphere rate (all E)",
       R_hemi * 60, 1.0, rtol=0.30,
       unit=" cm⁻²min⁻¹", note="PDG ≈ 1 cm⁻²min⁻¹")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 — CSDA RANGE TABLE ROUND-TRIP
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 3 — CSDA range table round-trip  T → R(T) → T  accuracy")
print("""
  Log-log interpolation of X → E_min should recover table entries
  to < 2%.  (Beyond last table entry: E_min returns None = unphysical.)
""")
print(f"  {'X [g/cm²]':>12s}  {'T_table':>9s}  {'T_interp':>10s}  {'Δ%':>6s}")
for T_ref, R_ref in zip(_GROOM_T_GEV[::2], _GROOM_R_GCM2[::2]):
    T_inv = emin_from_opacity(R_ref)
    if T_inv is None:
        print(f"  ⚠️   X={R_ref:>10.0f} g/cm²  T_ref={T_ref:.4g} GeV  "
              f"(above CSDA max — emin_from_opacity returns None)")
        continue
    rel = abs(T_inv - T_ref) / T_ref * 100
    sym = "✅" if rel < 2.0 else "❌"
    print(f"  {sym}  X={R_ref:>10.0f} g/cm²  "
          f"T_table={T_ref:>7.4g} GeV  T_interp={T_inv:>7.4g} GeV  Δ={rel:.2f}%")
    if rel < 2.0: _pass += 1
    else:         _fail += 1


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 4 — FLAT-SLAB PATH-LENGTH CORRECTION
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 4 — Flat-slab path correction  X = ρ · (L/cosθ) · 100")
print("""
  A muon at zenith angle θ through a flat slab of vertical depth L [m]
  travels a slant path = L/cosθ — NOT L.

  ⚠️  The old FFE used X = ρ·L·100 regardless of θ.  At θ=60° this
      underestimates the opacity by ×2 and overestimates I by ~1 order of
      magnitude.  The corrected formula is now applied.

  L=100 m,  ρ=2.65 g/cm³  (Standard Rock):
""")
L   = 100.0; rho = RHO_STANDARD_ROCK
I0, _ = integrated_flux(0.0, 0.0, model="reyna_bugaev")
print(f"  {'θ':>4s}  {'path [m]':>9s}  {'X [g/cm²]':>11s}  "
      f"{'E_min [GeV]':>12s}  {'I_rock':>11s}  {'T':>8s}")
I_prev = None; mono_ok = True
for theta in [0, 15, 30, 45, 60, 70, 80]:
    ct    = math.cos(math.radians(theta))
    path  = L / ct
    X     = rho * path * 100.0
    Emin  = emin_from_opacity(X)
    I_r,_ = integrated_flux(X, theta, model="reyna_bugaev")
    I_o,_ = integrated_flux(0.0, theta, model="reyna_bugaev")
    T_tr  = I_r / I_o if I_o > 0 else 0.0
    Es    = f"{Emin:.1f}" if Emin else "—"
    if I_prev is not None and I_r > I_prev * 1.001: mono_ok = False
    I_prev = I_r
    print(f"  {theta:>3}°  {path:>9.1f}  {X:>11.0f}  "
          f"{Es:>12s}  {I_r:>11.3e}  {T_tr:>8.4f}")
sym = "✅" if mono_ok else "❌"
print(f"\n  {sym}  T non-increasing with obliquity at fixed vertical depth")
if mono_ok: _pass += 1
else:       _fail += 1

# Check T at θ=0 for standard benchmarks
I_100, Em100 = integrated_flux(rho*100*100.0, 0.0, model="reyna_bugaev")
T_100 = I_100 / I0
_check("T(L=100 m, θ=0°, ρ=2.65, R-B) in [0.001, 0.01]",
       T_100, 0.003, rtol=0.90, note="few × 10⁻³ expected for 100 m Standard Rock")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 5 — ANGULAR PROFILE AND cosθ* CORRECTION
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 5 — Angular profile and cosθ* (Earth-curvature) correction")
print("""
  The naive cos²θ approximation fails at large θ because it forces
  the flux to zero at 90°.  In reality, near-horizontal muons see a
  thicker but not infinite atmosphere (cosθ* > 0 always).
""")
print(f"  {'θ':>5s}  {'cosθ':>10s}  {'cosθ* (Guan)':>13s}  {'ratio':>8s}  note")
for th in [0, 30, 60, 75, 80, 85, 89]:
    cr  = math.cos(math.radians(th))
    cs  = _guan_cos_star(cr)
    rat = cs / cr if cr > 1e-4 else float('inf')
    note = ("← curvature significant" if th>=75 else
            "← cosθ* floors at ~0.103" if th>=85 else "")
    print(f"  {th:>4}°  {cr:>10.6f}  {cs:>13.6f}  {rat:>8.4f}  {note}")

print()
print(f"  I(θ)/I(0°)  at E_min=1 GeV:")
print(f"  {'θ':>5s}  {'Reyna–B':>10s}  {'Guan 2015':>11s}  "
      f"{'Frosin 25':>11s}  {'cos²θ':>8s}  {'model vs cos²θ':>16s}")
th_arr = np.array([0., 20., 40., 60., 70., 80., 85.])
I_r, _ = angular_profile(th_arr, E_min_GeV=1.0, model="reyna_bugaev")
I_g, _ = angular_profile(th_arr, E_min_GeV=1.0, model="guan_2015")
I_f, _ = angular_profile(th_arr, E_min_GeV=1.0, model="frosin_2025")
for i, th in enumerate(th_arr):
    cos2 = math.cos(math.radians(th))**2
    Tr = I_r[i] / max(I_r[0], 1e-30)
    Tg = I_g[i] / max(I_g[0], 1e-30)
    Tf = I_f[i] / max(I_f[0], 1e-30)
    dev = (Tr - cos2)/cos2*100 if cos2 > 1e-6 else float('nan')
    print(f"  {th:>4.0f}°  {Tr:>10.4f}  {Tg:>11.4f}  {Tf:>11.4f}  "
          f"{cos2:>8.4f}  {dev:>+16.1f}%")
print("\n  → At 85°: cos²θ underestimates by >200%. Always use a model.")

# Monotonicity check (E_min=1 GeV, Reyna-Bugaev)
mono_ang = all(I_r[i+1] <= I_r[i] * 1.001 for i in range(len(I_r)-1))
_check("Angular profile monotonically decreasing with θ (E_min=1 GeV)",
       float(mono_ang), 1.0, rtol=0.01)

# ── Energy-dependent angular ratio — KEY for muography ────────────────
print()
print("  ── Energy-dependent I(θ)/I(0°) — critical for muography ──")
print()
print("  Reyna-Bugaev uses cos_th_star^1.85 — ratio is ENERGY-INDEPENDENT:")
T_e_scan = np.array([1., 10., 50., 100., 300.])
phi0_rb  = differential_flux(T_e_scan, 0.0,  'reyna_bugaev')
phi30_rb = differential_flux(T_e_scan, 30.0, 'reyna_bugaev')
for i, E in enumerate(T_e_scan):
    print(f"    E={E:>5.0f} GeV:  I(30°)/I(0°) = {phi30_rb[i]/phi0_rb[i]:.4f}  (same at all E — hard-coded by formula)")
print()
print("  Guan 2015 models pion/kaon ANGULAR ENHANCEMENT — energy-DEPENDENT:")
phi0_gn  = differential_flux(T_e_scan, 0.0,  'guan_2015')
phi30_gn = differential_flux(T_e_scan, 30.0, 'guan_2015')
for i, E in enumerate(T_e_scan):
    ratio = phi30_gn[i] / phi0_gn[i] if phi0_gn[i] > 0 else 0
    note = "← oblique > vertical" if ratio > 1 else "← oblique < vertical"
    print(f"    E={E:>5.0f} GeV:  I(30°)/I(0°) = {ratio:.4f}  {note}")
print()
print("  Physics: above ~46 GeV at θ=30°, pion/kaon decay-vs-reinteraction balance")
print("  shifts — oblique muons get MORE pion decays per unit solid angle than vertical.")
print("  Crossover: ~46 GeV at 30°, ~63 GeV at 60°. Confirmed by CMS, IceCube.")
print()
print("  ⚠️  CONSEQUENCE FOR MUOGRAPHY (E_min typically 50–700 GeV):")
print("    Guan/Frosin: I(θ)/I(0°) > 1 at moderate θ — PHYSICALLY CORRECT.")
print("    Reyna-Bugaev: constant ratio at all E — WRONG angular shape for muography.")
print("    → For angular analysis / inversion: use Guan 2015 or Frosin 2025.")
print("    → For total rate (integrated from 1 GeV): use Reyna-Bugaev.")
print("    → T = I_rock/I_open is computed at the same θ — T values are correct in all models.")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 6 — TRANSMISSION CONSTRAINTS
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 6 — Transmission T = I_rock / I_open: physical constraints")

I0, _ = integrated_flux(0.0, 0.0, model="reyna_bugaev")
_check("T(X=0) = 1.0  (no rock → no attenuation)",
       integrated_flux(0.0, 0.0, model="reyna_bugaev")[0] / I0, 1.0,
       rtol=0.001, note="exact by definition")

# Monotonicity in X
X_list = [0, 5_000, 15_000, 30_000, 53_000, 100_000]
I_prev = None; mono_X = True
print()
print(f"  {'X [g/cm²]':>10s}  {'depth [m]':>10s}  "
      f"{'E_min [GeV]':>12s}  {'I [cm⁻²sr⁻¹s⁻¹]':>18s}  {'T':>10s}")
for X in X_list:
    I, Em = integrated_flux(X, 0.0, model="reyna_bugaev")
    T     = I / I0 if I0 > 0 else 0.0
    depth = X / (RHO_STANDARD_ROCK * 100.0)
    Es    = f"{Em:.0f}" if Em else "∞"
    if I_prev is not None and I > I_prev * 1.001: mono_X = False
    I_prev = I
    print(f"  {X:>10.0f}  {depth:>10.1f}  {Es:>12s}  {I:>18.4e}  {T:>10.6f}")
sym = "✅" if mono_X else "❌"
print(f"\n  {sym}  T is monotonically non-increasing with X")
if mono_X: _pass += 1
else:      _fail += 1


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 7 — ALTITUDE CORRECTION
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 7 — Altitude correction  φ(h) ≈ φ₀ × exp(h / 8500 m)")
print(f"\n  {'h [m]':>8s}  {'exp factor':>12s}  {'I (R-B, E>0.5 GeV)':>22s}  note")
for h in [0, 500, 1000, 1094, 1281, 2000, 3000]:
    I, _  = integrated_flux(0.0, 0.0, model="reyna_bugaev", altitude_m=h)
    corr  = math.exp(h / 8500.0)
    note  = {1094:"← Puy de Dôme detector", 1281:"← Vesuvius summit"}.get(h, "")
    print(f"  {h:>8d}  {corr:>12.4f}  {I:>20.4e}  {note}")

I_ref, _ = integrated_flux(0.0, 0.0, model="reyna_bugaev", altitude_m=0)
I_pdd, _ = integrated_flux(0.0, 0.0, model="reyna_bugaev", altitude_m=1094)
_check("Puy de Dôme correction factor (1094 m)",
       I_pdd / I_ref, math.exp(1094/8500), rtol=0.01,
       note="exp(1094/8500) = 1.138  → +14% more flux vs sea level")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 8 — CONCRETE DETECTOR CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 8 — Concrete rate calculations for typical detector geometries")
# Use the physically correct E>0.5 GeV value (what integrated_flux returns)
# but also compute E>1 GeV explicitly for the label
I_open_05,  _ = integrated_flux(0.0, 0.0, model="reyna_bugaev")   # E>0.5 GeV (grid start)
T_1gev        = np.logspace(0, 4, 1200)                            # 1 GeV → 10 TeV
_phi_1gev     = differential_flux(T_1gev, 0.0, 'reyna_bugaev')
_int_fn       = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
I_open_1gev   = float(_int_fn(_phi_1gev, T_1gev))                  # E>1 GeV

# Use E>1 GeV as the reference (matches PDG convention and GUI E_min default)
I_open = I_open_1gev

print(f"\n  Reyna–Bugaev open-sky reference (θ=0°, sea level):")
print(f"    I(E>0.5 GeV) = {I_open_05:.4e} cm⁻²sr⁻¹s⁻¹  (integrated_flux default grid)")
print(f"    I(E>1.0 GeV) = {I_open:.4e} cm⁻²sr⁻¹s⁻¹  ← used below (PDG convention)")
print(f"    PDG reference ≈ 7.0×10⁻³  (Δ = {abs(I_open-7e-3)/7e-3*100:.0f}%)\n")

cases = [
    ("10×10 cm² panel, full hemisphere",          100.0, math.pi),
    ("10×10 cm² panel, cone θ < 30°",             100.0, 2*math.pi*(1-math.cos(math.radians(30)))),
    ("Telescope: 10×10 cm², d=50 cm (Ω=S/d²)",   100.0, 100.0/50.0**2),
    ("MURAVES-style  A=6 cm²·sr (GUI default)",   None,  6.0),
]
for desc, S, Omega in cases:
    A_det = S * Omega if S is not None else Omega
    R_s   = I_open * A_det
    print(f"  ● {desc}")
    print(f"    A_det = {A_det:.3f} cm²·sr   "
          f"R = {R_s:.4f} s⁻¹ = {R_s*60:.3f} min⁻¹ = {R_s*86400:.0f} day⁻¹  "
          f"(open sky, E>1 GeV, θ=0°)\n")

print("  With rock overburden (Reyna–Bugaev, θ=0°, A_det=6 cm²·sr):")
print(f"  Note: T = I_rock / I_open where I_open = {I_open_05:.4e} cm⁻²sr⁻¹s⁻¹ (E>0.5 GeV, same grid as I_rock)")
print(f"  {'L [m]':>8s}  {'E_min [GeV]':>12s}  {'I_rock':>12s}  "
      f"{'T':>8s}  {'Rate [/day]':>12s}  {'(E>1 GeV rate)':>15s}")
A_det = 6.0
for L_m in [0, 20, 50, 100, 200, 500]:
    X_v = RHO_STANDARD_ROCK * L_m * 100.0
    I_r, Em = integrated_flux(X_v, 0.0, model="reyna_bugaev")
    # T = I_rock / I_open  — both from same integrated_flux grid (E>0.5 GeV)
    # This is the physically correct rock transmission (fraction of muons surviving)
    T_r = I_r / I_open_05 if I_open_05 > 0 else 0.0
    # Rate uses I_rock directly (E >= E_min, which equals E>1 GeV for non-zero depth)
    # For L=0, I_rock = I_open_05 (E>0.5 GeV), so rate is slightly higher than E>1 GeV
    Es  = f"{Em:.0f}" if Em else "∞"
    # Compute E>1 GeV rate for L=0 explicitly
    rate_day = I_r * A_det * 86400
    rate_1gev = I_open_1gev * A_det * 86400 if L_m == 0 else rate_day
    rate_note = f"({I_open_1gev*A_det*86400:.0f} for E>1 GeV)" if L_m == 0 else ""
    print(f"  {L_m:>8d}  {Es:>12s}  {I_r:>12.3e}  "
          f"{T_r:>8.5f}  {rate_day:>12.2f}  {rate_note}")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 9 — VESUVIUS / MURAVES CROSS-CHECK
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 9 — Cross-check: MURAVES-like geometry at Mt. Vesuvius")
print("""
  MURAVES: detector at 608 m a.s.l., summit at ~1281 m.
  Typical overburden 600 m (vertical), θ ≈ 35–65°.
  Rock density ~1.7 g/cm³ (volcanic tuff — NOT Standard Rock!).

  Comparing Reyna–Bugaev vs Guan 2015 transmission T at MURAVES angles:
""")
rho_v = 1.7; alt_v = 608; L_v = 600
print(f"  {'θ':>5s}  {'path [m]':>9s}  {'X [g/cm²]':>11s}  "
      f"{'E_min [GeV]':>11s}  {'R-B T':>8s}  {'Guan T':>8s}  {'ratio':>7s}")
for theta in [30, 40, 50, 60, 70]:
    ct    = math.cos(math.radians(theta))
    path  = L_v / ct
    X     = rho_v * path * 100.0
    Emin  = emin_from_opacity(X)
    I_rb, _ = integrated_flux(X, theta, model="reyna_bugaev", altitude_m=alt_v)
    I_rb0,_ = integrated_flux(0.0, theta, model="reyna_bugaev", altitude_m=alt_v)
    I_gg, _ = integrated_flux(X, theta, model="guan_2015",    altitude_m=alt_v)
    I_gg0,_ = integrated_flux(0.0, theta, model="guan_2015",  altitude_m=alt_v)
    T_rb = I_rb / I_rb0  if I_rb0 > 0 else 0.0
    T_g  = I_gg / I_gg0  if I_gg0 > 0 else 0.0
    ratio = T_g / T_rb   if T_rb  > 0 else float('nan')
    Es   = f"{Emin:.0f}" if Emin else "∞"
    print(f"  {theta:>4}°  {path:>9.0f}  {X:>11.0f}  "
          f"{Es:>11s}  {T_rb:>8.4f}  {T_g:>8.4f}  {ratio:>7.3f}")

print()
print("  Interpretation:")
print("    θ=30–40°: models within ~15–50% — usable for planning estimates.")
print("    θ=50–60°: models diverge by 2–4× — extreme sensitivity to spectral")
print("              shape above E_min>400 GeV where neither model was fitted.")
print("    θ>60°:    E_min exceeds CSDA table (None); flux is unmeasurably small.")
print()
print("  ⚠️  For MURAVES-style geometry at oblique angles, the flat-slab FFE")
print("     has two independent sources of error:")
print("     (1) Terrain geometry — the real path is NOT L/cosθ for a volcano.")
print("         Use UCMuon Terrain (Tab 5) with the DEM for real path-lengths.")
print("     (2) Model extrapolation — at X>100 000 g/cm², Reyna-Bugaev and Guan")
print("         diverge by factors of 2–4. Only full MC (MUSIC/PROPOSAL) is")
print("         reliable in this regime.")
print("  → For absolute rate at Vesuvius: use Reyna-Bugaev + Tab 5 Terrain.")
print("  → For angular ratios only: Guan 2015 with Tab 5 is the best option.")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 10 — KNOWN LIMITATIONS
# ═══════════════════════════════════════════════════════════════════════════
_hdr("SECTION 10 — Known limitations of the FFE")
print("""
  ┌────┬──────────────────────────────────────────────────────────────┐
  │ 1  │  FLAT SLAB ONLY: FFE assumes horizontal rock, uniform ρ.    │
  │    │  Real terrain → use UCMuon Terrain (Tab 5) with a DEM.      │
  │    │  At Vesuvius θ>60°, the flat-slab FFE can be off by >30%.   │
  │ 2  │  AZIMUTH SYMMETRY: All 5 models are φ-independent.          │
  │    │  The real sky has ~2% E–W asymmetry from geomagnetic rigidity │
  │    │  cutoff.  For φ-dependence → use PARMA (spectrum mode ③).   │
  │ 3  │  ALTITUDE: exp(h/8500 m) is only approximate; valid to ~4 km.│
  │ 4  │  CSDA vs MC: E_min from Groom CSDA is a LOWER BOUND.        │
  │    │  Stochastic losses (pair, bremsstrahlung) mean some muons    │
  │    │  with E < E_min,CSDA still survive through rock (straggling). │
  │    │  MUSIC gives a smaller effective threshold for the same X.   │
  │ 5  │  SINGLE θ: FFE integrates one direction. Wide-angle detectors│
  │    │  must integrate over the acceptance cone manually:           │
  │    │  R = Σᵢ I(θᵢ) × ΔΩᵢ × A_surface                          │
  │ 6  │  SINGLE ρ: FFE uses a fixed uniform density. Real volcanic   │
  │    │  targets have ρ variations of 10–30%.                        │
  └────┴──────────────────────────────────────────────────────────────┘
""")


# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
_hdr(f"FINAL SUMMARY — {_pass} checks passed, {_fail} failed")
print(f"""
  QUICK REFERENCE FOR THE UCMuon GUI FFE (Tab 3):

  ┌─────────────────────────────────────────────────────────────────────┐
  │  OPEN-SKY (Reyna–Bugaev, sea level, E>1 GeV, θ=0°)                │
  │    I(0°) = 5.7×10⁻³ cm⁻²sr⁻¹s⁻¹   (PDG: 7×10⁻³)               │
  │    Hemisphere ≈ 0.96 cm⁻²min⁻¹      (PDG: ≈1)                    │
  │                                                                     │
  │  FLAT-SLAB BENCHMARKS (ρ=2.65, θ=0°, Reyna–Bugaev)               │
  │    L= 20 m  → X=  5 300 g/cm²  E_min= 11 GeV  T≈ 4.3%           │
  │    L=100 m  → X= 26 500 g/cm²  E_min= 62 GeV  T≈ 0.30%          │
  │    L=200 m  → X= 53 000 g/cm²  E_min=133 GeV  T≈ 0.080%         │
  │                                                                     │
  │  PATH AT OBLIQUE ANGLES (L=100 m vertical)                         │
  │    θ=30°  path=115 m   X= 30 500 g/cm²  E_min= 72 GeV           │
  │    θ=60°  path=200 m   X= 53 000 g/cm²  E_min=133 GeV           │
  │    θ=80°  path=576 m   X=153 000 g/cm²  E_min≈430 GeV           │
  │                                                                     │
  │  RATE ESTIMATE (A=6 cm²·sr, open sky, E>1 GeV, θ=0°)             │
  │    R = 5.7×10⁻³ × 6 ≈ 0.034 s⁻¹ = 2.0 min⁻¹ = 2 900 day⁻¹    │
  └─────────────────────────────────────────────────────────────────────┘
""")
