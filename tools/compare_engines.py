#!/usr/bin/env python3
"""
UCMuon engine comparison diagnostic
====================================
Runs all available engines (BB+MS, MUSIC, UCMuon-MC, PROPOSAL) on the
same input with the same depth and material, then compares survival rates,
energy loss, and lateral displacement.

Usage:
    python3 tools/compare_engines.py [input_file] [depth_m] [rho_gcm3]

Defaults:
    input_file = output/muons_surface.dat
    depth_m    = 90.0
    rho        = 2.65  (Standard Rock)
"""

import sys, os, subprocess, time, shutil
from pathlib import Path
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE      = Path(__file__).resolve().parent.parent   # tools/ → repo root
BIN        = _HERE / "bin"
GUI_DIR    = _HERE / "gui"
OUTPUT_DIR = _HERE / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# ── Parameters from CLI ────────────────────────────────────────────────────────
INFILE  = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else _HERE / "output/muons_surface.dat"
DEPTH_M = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
RHO     = float(sys.argv[3]) if len(sys.argv) > 3 else 2.65

# Standard Rock (Groom 2001 / PDG)
X0_GCM2  = 26.54   # radiation length [g/cm²]
X0_CM    = X0_GCM2 / RHO  # correct X0 in cm = 10.02 for Standard Rock
MAT_ID   = 1        # Standard Rock for all engines

OPACITY  = DEPTH_M * 100.0 * RHO   # [g/cm²]
MWE      = OPACITY / 100.0          # m.w.e.

print(f"\n{'='*64}")
print(f"  UCMuon Engine Comparison Diagnostic  (v2 — source-plane aware)")
print(f"{'='*64}")
print(f"  Input file : {INFILE}")
print(f"  Depth      : {DEPTH_M} m")
print(f"  Density    : {RHO} g/cm³")
print(f"  Opacity    : {OPACITY:.1f} g/cm²  ({MWE:.1f} m.w.e.  =  {MWE/1000:.4f} km.w.e.)")
print(f"{'='*64}\n")

# ── Source-plane auto-detection (mirrors driver logic) ─────────────────────────
def detect_source_plane(path, ncols):
    xs, ys, zs = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 4:
                continue
            xs.append(float(p[1])); ys.append(float(p[2])); zs.append(float(p[3]))
    xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)
    stds = [xs.std(), ys.std(), zs.std()]
    axis = int(np.argmin(stds))
    names = {0: "YZ (depth=X)", 1: "XZ (depth=Y)", 2: "XY (depth=Z)"}
    return axis, names[axis], stds

# ── Sanity-check input file ────────────────────────────────────────────────────
if not INFILE.exists():
    print(f"ERROR: input file not found: {INFILE}")
    print("       Run the generator first (Tab 1) or pass a path as argv[1].")
    sys.exit(1)

# ── Probe input file ──────────────────────────────────────────────────────────
def detect_ncols(path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            return len(line.split())
    return 0

def count_muons(path, ncols):
    n_total = n_hit = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < ncols:
                continue
            n_total += 1
            if ncols == 14:
                if int(float(parts[12])) == 1:
                    n_hit += 1
            else:
                n_hit += 1
    return n_total, n_hit

ncols = detect_ncols(INFILE)
if ncols == 0:
    print("ERROR: input file is empty or all comments.")
    sys.exit(1)

n_total, n_hit = count_muons(INFILE, ncols)
print(f"  Input  : {n_total} muons total, {ncols} columns")
if ncols == 14:
    print(f"           {n_hit} with hit_flag=1  (will be transported)")

src_axis, src_name, src_stds = detect_source_plane(INFILE, ncols)
print(f"  Source plane: {src_name}")
print(f"  Pos. σ(x,y,z): {src_stds[0]:.1f}  {src_stds[1]:.1f}  {src_stds[2]:.1f}  cm")
if src_axis != 2:
    print(f"\n  *** NOTE: Non-XY source plane detected. ***")
    print(f"  Both BB+MS (fixed) and UCMuon-MC (fixed) now use the")
    print(f"  correct depth-direction component for slant computation.")
    print(f"  Older engine runs (before this fix) gave incorrect results.")
print()

# ── Run a transport process ────────────────────────────────────────────────────
def run_binary(label, cmd, stdin_str, timeout=600):
    cmd = [str(c) for c in cmd]
    if not Path(cmd[0]).exists():
        return False, f"not found: {cmd[0]}"
    print(f"  [{label}] running ...")
    t0 = time.time()
    try:
        res = subprocess.run(
            cmd,
            input=stdin_str,
            capture_output=True, text=True,
            cwd=str(BIN),
            timeout=timeout,
            env={**os.environ, "OMP_NUM_THREADS": "4",
                 "PROPOSAL_LOG_LEVEL": "error", "SPDLOG_LEVEL": "error"},
        )
        elapsed = time.time() - t0
        if res.returncode != 0:
            snippet = (res.stderr or res.stdout)[:600]
            return False, f"exit {res.returncode}\n{snippet}"
        # Print last summary lines
        lines = [l for l in res.stdout.splitlines() if l.strip()]
        for l in lines[-8:]:
            print(f"      {l}")
        print(f"  [{label}] done  ({elapsed:.1f} s)\n")
        return True, res.stdout
    except subprocess.TimeoutExpired:
        return False, f"timeout (>{timeout} s)"
    except Exception as e:
        return False, str(e)

# ── Load 18-column output ─────────────────────────────────────────────────────
# 18-col layout:
#  0=eventid  1=xs  2=ys  3=zs  4=E_srf  5=theta_srf  6=phi_srf  7=charge
#  8=alive    9=x_ug 10=y_ug 11=z_ug 12=E_ug
#  13=cx_ug  14=cy_ug 15=cz_ug  16=theta_ug  17=phi_ug
def load_output(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 13:
                rows.append([float(x) for x in parts[:18]])
    return np.array(rows) if rows else None

# ── Build stdin strings ────────────────────────────────────────────────────────

# 1) BB+MS
BB_OUT = OUTPUT_DIR / "compare_bb.dat"
bb_stdin = "\n".join([
    str(INFILE),     # infile (absolute)
    str(BB_OUT),     # outfile (absolute)
    "0",             # transport_all: hit_flag=1 only
    str(ncols),      # ncols hint
    str(DEPTH_M),    # depth [m]
    str(MAT_ID),     # material (1=Standard Rock)
    "1",             # multiple scattering ON
]) + "\n"

# 2) MUSIC
MUSIC_OUT  = OUTPUT_DIR / "compare_music.dat"
_eloss_ok  = (BIN / "music-eloss-rock.dat").exists()
_xsec_ok   = (BIN / "music-cross-sections-rock.dat").exists()
music_init = 1 if (_eloss_ok and _xsec_ok) else 0
music_stdin = "\n".join([
    str(INFILE),     # infile
    str(MUSIC_OUT),  # outfile
    str(RHO),        # density [g/cm³]
    str(X0_GCM2),    # radiation length [g/cm²]
    str(DEPTH_M),    # depth [m]
    "1",             # idim=1 (3-D lateral scattering ON)
    "1",             # idim1=1
    "-30",           # minv
    str(music_init), # init: 1=load tables, 0=recalculate
    str(MAT_ID),     # mat_id (1=Standard Rock)
    "0",             # transport_all
    "",              # press Enter
]) + "\n"

# 3) UCMuon-MC (ucmuon_stochastic_driver.py)
STOCHASTIC_OUT = OUTPUT_DIR / "compare_stochastic.dat"
stochastic_stdin = "\n".join([
    str(INFILE),     # infile
    str(STOCHASTIC_OUT),  # outfile
    str(DEPTH_M),    # depth [m]
    str(RHO),        # density [g/cm³]
    str(X0_CM),      # X0 [cm] — correct value: X0_gcm2 / rho = 26.54/2.65 = 10.02
    str(MAT_ID),     # mat_id (1=Standard Rock)
    "0",             # transport_all
    str(ncols),      # ncols hint
    "0",             # n_steps (0=auto)
    "0.05",          # v_cut
    "1",             # ms_enable
]) + "\n"

# 4) PROPOSAL (proposal_driver.py)
PROP_OUT  = OUTPUT_DIR / "compare_proposal.dat"
prop_stdin = "\n".join([
    str(INFILE),     # infile
    str(PROP_OUT),   # outfile
    str(DEPTH_M),    # depth [m]
    str(MAT_ID),     # medium_type (1=StandardRock)
    "0",             # transport_all
    "500",           # e_cut [MeV]
    "0.001",         # v_cut
    "2",             # scattering: HighlandIntegral
    "",              # tables_dir (default)
]) + "\n"

def _find_proposal_python():
    """
    Find a Python executable that can actually import proposal without crashing.
    Tests candidates via subprocess so a segfaulting install doesn't kill us.
    """
    home = Path.home()
    candidates = [sys.executable]
    # Common venv / conda locations
    for pattern in [
        home / "venvs/*/bin/python",
        home / "venvs/*/bin/python3",
        home / "miniforge3/bin/python3",
        home / "miniforge3/bin/python",
        home / "miniforge3/envs/*/bin/python3",
        home / "miniforge3/envs/*/bin/python",
        home / "anaconda3/bin/python3",
        home / "anaconda3/envs/*/bin/python3",
        home / "opt/anaconda3/bin/python3",
    ]:
        candidates.extend(sorted(home.glob(str(pattern.relative_to(home)))))
    for name in ["python3", "python"]:
        found = shutil.which(name)
        if found:
            candidates.append(found)

    seen = set()
    for py in candidates:
        py = str(py)
        if py in seen or not Path(py).exists():
            continue
        seen.add(py)
        try:
            r = subprocess.run(
                [py, "-c", "import proposal; print('ok')"],
                capture_output=True, text=True, timeout=20,
            )
            if r.returncode == 0 and "ok" in r.stdout:
                return py
        except Exception:
            pass
    return None

# ── Run all engines ───────────────────────────────────────────────────────────
print("=" * 64)
print("  RUNNING ENGINES")
print("=" * 64 + "\n")

raw_results = {}

# 1) BB+MS
BB_BIN = BIN / "ucmuon_transport_bb_omp"
ok, msg = run_binary("BB+MS", [BB_BIN], bb_stdin)
if ok and BB_OUT.exists():
    raw_results["BB+MS"] = BB_OUT
else:
    print(f"  [BB+MS] FAILED: {msg}\n")

# 2) MUSIC
MUSIC_BIN = BIN / "ucmuon_transport_music_omp"
if not MUSIC_BIN.exists():
    print("  [MUSIC] binary not found — skipping.\n")
elif music_init == 0 and not (BIN / "music-double-diff-rock.dat").exists():
    print("  [MUSIC] music-double-diff-rock.dat missing (needed for init=0) — skipping.\n")
else:
    ok, msg = run_binary("MUSIC", [MUSIC_BIN], music_stdin)
    if ok and MUSIC_OUT.exists():
        raw_results["MUSIC"] = MUSIC_OUT
    else:
        print(f"  [MUSIC] FAILED: {msg}\n")

# 3) UCMuon-MC
STOCHASTIC_DRV = GUI_DIR / "ucmuon_stochastic_driver.py"
if not STOCHASTIC_DRV.exists():
    print("  [UCMuon-MC] ucmuon_stochastic_driver.py not found — skipping.\n")
else:
    ok, msg = run_binary("UCMuon-MC", [sys.executable, str(STOCHASTIC_DRV)], stochastic_stdin)
    if ok and STOCHASTIC_OUT.exists():
        raw_results["UCMuon-MC"] = STOCHASTIC_OUT
    else:
        print(f"  [UCMuon-MC] FAILED: {msg}\n")

# 4) PROPOSAL
PROP_DRV = GUI_DIR / "proposal_driver.py"
if not PROP_DRV.exists():
    print("  [PROPOSAL] proposal_driver.py not found — skipping.\n")
else:
    print("  [PROPOSAL] locating a working Python+proposal install ...", flush=True)
    _prop_python = _find_proposal_python()
    if _prop_python is None:
        print("  [PROPOSAL] no Python found that can import proposal — skipping.\n"
              "             Install with: pip install proposal\n")
    else:
        print(f"  [PROPOSAL] using {_prop_python}\n")
        ok, msg = run_binary("PROPOSAL", [_prop_python, str(PROP_DRV)], prop_stdin, timeout=900)
        if ok and PROP_OUT.exists():
            raw_results["PROPOSAL"] = PROP_OUT
        else:
            print(f"  [PROPOSAL] FAILED: {msg}\n")

if not raw_results:
    print("No engines ran successfully. Check binary paths and input file.")
    sys.exit(1)

# ── Load outputs ───────────────────────────────────────────────────────────────
data = {}
for name, path in raw_results.items():
    d = load_output(path)
    if d is not None and len(d) > 0:
        data[name] = d
    else:
        print(f"  WARNING: {name} output empty or unreadable at {path}")

if not data:
    print("No output data loaded.")
    sys.exit(1)

# ── Per-engine statistics ──────────────────────────────────────────────────────
print("=" * 64)
print("  RESULTS SUMMARY")
print(f"  depth={DEPTH_M} m | ρ={RHO} g/cm³ | X={OPACITY:.0f} g/cm² | {MWE:.1f} m.w.e.")
print("=" * 64)

stats = {}
for name, d in data.items():
    alive     = d[:, 8].astype(int)
    E_srf     = d[:, 4]
    E_ug      = d[:, 12]
    n_tot     = len(d)
    n_surv    = int((alive == 1).sum())
    surv_frac = n_surv / n_tot if n_tot > 0 else 0.0
    mask_s    = alive == 1
    E_srf_s   = E_srf[mask_s]
    E_ug_s    = E_ug[mask_s]
    dE_s      = E_srf_s - E_ug_s
    dx        = (d[:, 9]  - d[:, 1]) / 100.0
    dy        = (d[:, 10] - d[:, 2]) / 100.0
    dr        = np.sqrt(dx**2 + dy**2)
    dr_s      = dr[mask_s]
    stats[name] = dict(
        n_tot=n_tot, n_surv=n_surv, surv_frac=surv_frac,
        alive=alive, E_srf=E_srf, E_ug=E_ug,
        E_srf_s=E_srf_s, E_ug_s=E_ug_s, dE_s=dE_s,
        dr=dr, dr_s=dr_s,
    )
    print(f"\n  ── {name} {'─'*(44-len(name))}")
    print(f"  Transported          : {n_tot:,}")
    print(f"  Survived             : {n_surv:,}  ({surv_frac*100:.2f} %)")
    if n_surv > 0:
        print(f"  E_surface (survived) : mean={np.mean(E_srf_s):.2f}  "
              f"median={np.median(E_srf_s):.2f}  std={np.std(E_srf_s):.2f}  GeV")
        print(f"  E_underground        : mean={np.mean(E_ug_s):.2f}  "
              f"median={np.median(E_ug_s):.2f}  std={np.std(E_ug_s):.2f}  GeV")
        print(f"  ΔE (mean loss)       : {np.mean(dE_s):.2f} GeV")
        print(f"  Lateral disp (surv)  : mean={np.mean(dr_s):.2f} m  p95={np.percentile(dr_s,95):.2f} m")

# ── Survival rate table ────────────────────────────────────────────────────────
print(f"\n{'='*64}")
print(f"  SURVIVAL RATE TABLE")
print(f"{'='*64}")
ref_name = next(iter(stats))
ref_surv = stats[ref_name]["surv_frac"]
print(f"  {'Engine':<26}  {'N_survived':>12}  {'Rate':>8}  {'Delta vs '+ref_name:>20}")
print(f"  {'─'*26}  {'─'*12}  {'─'*8}  {'─'*20}")
for name, s in stats.items():
    delta = (s["surv_frac"] - ref_surv) * 100.0
    delta_str = f"{delta:+.2f} pp" if name != ref_name else "  (reference)"
    print(f"  {name:<26}  {s['n_surv']:>12,}  {s['surv_frac']*100:>7.2f}%  {delta_str:>20}")

# ── Muon-by-muon pairwise comparison ──────────────────────────────────────────
names = list(stats.keys())
if len(names) >= 2:
    print(f"\n{'='*64}")
    print(f"  MUON-BY-MUON COMPARISON (requires same N muons in each output)")
    print(f"{'='*64}")
    ref = stats[names[0]]
    for other_name in names[1:]:
        s = stats[other_name]
        if ref["n_tot"] != s["n_tot"]:
            print(f"\n  {names[0]} vs {other_name}: different counts "
                  f"({ref['n_tot']} vs {s['n_tot']}) — skipping per-muon compare.")
            print(f"  Likely cause: one engine uses transport_all differently.")
            continue
        a0, a1 = ref["alive"], s["alive"]
        both_live = (a0 == 1) & (a1 == 1)
        both_dead = (a0 == 0) & (a1 == 0)
        only_ref  = (a0 == 1) & (a1 == 0)
        only_oth  = (a0 == 0) & (a1 == 1)
        n = len(a0)
        print(f"\n  {names[0]}  vs  {other_name}:")
        print(f"    Both survive          : {both_live.sum():>7,}  ({100*both_live.sum()/n:.1f}%)")
        print(f"    Both stopped          : {both_dead.sum():>7,}  ({100*both_dead.sum()/n:.1f}%)")
        print(f"    Only {names[0]:20s}: {only_ref.sum():>7,}  ({100*only_ref.sum()/n:.1f}%)")
        print(f"    Only {other_name:20s}: {only_oth.sum():>7,}  ({100*only_oth.sum()/n:.1f}%)")
        if both_live.sum() > 0:
            dE = ref["E_ug"][both_live] - s["E_ug"][both_live]
            print(f"    E_ug diff (co-surv)  : mean={np.mean(dE):+.3f}  "
                  f"median={np.median(dE):+.3f}  std={np.std(dE):.3f}  GeV")

# ── Physics parameter audit ────────────────────────────────────────────────────
print(f"\n{'='*64}")
print(f"  PHYSICS PARAMETER AUDIT")
print(f"  (explains systematic differences between engines)")
print(f"{'='*64}")
print(f"""
  Parameter        BB+MS        MUSIC        UCMuon Stoc.  PROPOSAL
  ─────────────────────────────────────────────────────────────────
  b_rad [cm²/g]    3.08e-6      internal     3.475e-6      internal
  X0 used [g/cm²]  26.54        internal     {X0_CM:.2f}×{RHO}={X0_CM*RHO:.1f}     N/A (Molière)
  Energy loss      CSDA mean    stochastic   stochastic    stochastic
  Landau fluct.    NO           YES          YES (Poisson) YES
  MS model         Highland     Molière      Highland      HighlandInt
  Decay            NO           NO           YES           YES
  Source plane     auto-detect  (see MUSIC)  auto-detect   (see driver)
  ─────────────────────────────────────────────────────────────────

  KEY DIFFERENCES (post source-plane fix):
  1. Source-plane convention (FIXED in BB+MS and UCMuon-MC):
     Input was generated with '{src_name}' source plane.
     All engines now use the correct depth-direction component for slant.
     Before this fix, ~50% of muons were killed incorrectly by BB+MS, causing
     BB+MS ≈ 18.9% vs UCMuon-MC ≈ 32.8% (14 pp spurious discrepancy).

  2. b_rad: BB+MS (3.08e-6) vs UCMuon-MC (3.475e-6) — +12.8%.
     This causes UCMuon-MC to have higher radiative losses, leading to
     LOWER survival at depth relative to BB+MS after the source-plane fix.

  3. CSDA vs stochastic: BB+MS (no Landau fluctuations) OVERestimates survival
     at depth because hard events are not sampled. UCMuon-MC correctly
     applies Poisson-sampled catastrophic losses.
     Combined with b_rad difference, the net effect depends on depth and energy.

  4. X0 unit conversion (FIXED in GUI and standalone driver):
     gui_stochastic_engine.py now converts rad_gcm2 / rho before passing to driver.
     ucmuon_stochastic_driver.py _MAT_DB default corrected to X0_cm=10.015 (was 26.48).
     This script uses X0={X0_CM:.4f} cm → X0_gcm2={X0_CM*RHO:.2f} g/cm² ✓
""")

# ── Expected literature context ────────────────────────────────────────────────
print(f"{'='*64}")
print(f"  EXPECTED BB+MS vs MUSIC  (literature benchmarks)")
print(f"{'='*64}")
if MWE < 2000:
    exp = f"~2–5 pp  (shallow — engines should nearly agree)"
elif MWE < 10000:
    exp = f"~5–15 pp  (medium — BB+MS overestimates survival)"
else:
    exp = f"~15–30 pp  (deep — Landau fluctuations dominate)"
print(f"  Depth {DEPTH_M:.0f} m → {MWE:.0f} m.w.e.")
print(f"  BB+MS vs MUSIC survival delta: {exp}")

# ── Optional matplotlib plot ───────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    COLORS = {
        "BB+MS":             "#2196F3",
        "MUSIC":             "#FF5722",
        "UCMuon-MC": "#4CAF50",
        "PROPOSAL":          "#9C27B0",
    }

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for ax_i, (title, xcol_key, xlabel) in enumerate([
        ("Surface energy of survivors", "E_srf_s", "E_surface [GeV]"),
        ("Underground energy of survivors", "E_ug_s",  "E_underground [GeV]"),
        ("Lateral displacement of survivors", "dr_s",   "Lateral disp. [m]"),
    ]):
        ax = axes[ax_i]
        for name, s in stats.items():
            if s["n_surv"] > 0:
                vals = s[xcol_key]
                ax.hist(vals, bins=60, density=True, histtype="step", linewidth=2,
                        color=COLORS.get(name, "gray"),
                        label=f"{name}  ({s['surv_frac']*100:.1f}%)")
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.set_title(title, fontsize=10)
        if ax_i < 2:
            ax.set_yscale("log")
        ax.legend(fontsize=8)

    surv_str = "  |  ".join(f"{n}: {s['surv_frac']*100:.1f}%" for n, s in stats.items())
    fig.suptitle(
        f"Engine comparison  —  depth={DEPTH_M} m, ρ={RHO} g/cm³, X={OPACITY:.0f} g/cm²  "
        f"({MWE:.0f} m.w.e.)\n"
        f"input: {INFILE.name}   |   survival: {surv_str}",
        fontsize=9
    )
    fig.tight_layout()

    plot_path = OUTPUT_DIR / "engine_comparison.png"
    fig.savefig(str(plot_path), dpi=150, bbox_inches="tight")
    print(f"\n  Plot saved → {plot_path}")

except ImportError:
    print("\n  (install matplotlib for plots: pip install matplotlib)")
except Exception as e:
    print(f"\n  Plot error: {e}")

print()
