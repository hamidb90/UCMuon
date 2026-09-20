#!/usr/bin/env python3
"""
read_phits_output.py
Parse PHITS tally output files from muon_rock.inp → phits_summary.csv + phits_timing.txt

Usage:
    python3 read_phits_output.py [--n N] [--outdir DIR]

Options:
    --n N      Override number of source particles (default: read from muon_rock.out)
    --outdir   Directory to write output CSVs (default: current directory)

Reads (in current directory):
    out_muon_transmission.out  → transmission at 6 depths (mesh=reg, axis=reg)
    out_muon_KE_spectrum.out   → exit KE spectra (mesh=reg, axis=eng, samepage=reg)
    out_muon_angle.out         → angle distribution (mesh=r-z, axis=cos, samepage=z)
    out_muon_edep_total.out    → T-Deposit layer traversal counts (mesh=reg, axis=reg)
    muon_rock.out              → for n_source and timing

Writes:
    phits_summary.csv   per-depth statistics (DepthCm, MWE, N_transmitted,
                        Transmission_%, MeanExitKE_GeV, MeanAngle_deg, ...)
    phits_timing.txt    elapsed time in Geant4 summary format
"""

import os, re, argparse
import numpy as np

# ── Geometry constants (must match muon_rock.inp) ────────────────────────────
SCORING_DEPTHS_CM = [100, 1000, 2500, 5000, 10000, 20000]
ROCK_DENSITY      = 2.65
SLAB_AREA_CM2     = 5000.0 * 5000.0   # 25 000 000 cm²

# Region pair index (1-6 in T-Cross) → scoring depth [cm]
REG_TO_DEPTH = {i + 1: d for i, d in enumerate(SCORING_DEPTHS_CM)}

# Angle tally z-column index (0-5) → depth [cm]
# samepage=z → 7 z-boundary columns in order: -20000,-10000,-5000,-2500,-1000,-100, 0
# col 6 (z=0) is the rock entry surface → skip
ANGLE_COL_TO_DEPTH = {0: 20000, 1: 10000, 2: 5000, 3: 2500, 4: 1000, 5: 100}

# KE spectrum: region column (0-5, samepage=reg) → depth [cm]
KE_COL_TO_DEPTH = {i: d for i, d in enumerate(SCORING_DEPTHS_CM)}

# T-Deposit region number (1-6) → layer-bottom depth [cm]
DEP_REG_TO_DEPTH = {1: 100, 2: 1000, 3: 2500, 4: 5000, 5: 10000, 6: 20000}


# ═══════════════════════════════════════════════════════════════════════════════
#  Low-level helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _is_data_line(line):
    """True if line is a numeric data row (not a header/comment)."""
    s = line.strip()
    if not s:
        return False
    if s[0] in ('#', '$'):
        return False
    if s[0].isalpha():
        return False
    return True


def _read_reg_cross_1d(path):
    """
    Parse T-Cross with mesh=reg, axis=reg.
    Actual column layout:  num  area  muon+  err  muon-  err
    Returns {pair_no: (muon_plus, muon_minus)} in units [1/cm²/src].
    """
    result = {}
    if not os.path.exists(path):
        print(f"  [WARN] not found: {path}")
        return result
    with open(path) as fh:
        for raw in fh:
            if not _is_data_line(raw):
                continue
            parts = raw.split()
            if len(parts) < 6:
                continue
            try:
                no      = int(parts[0])
                muplus  = float(parts[2])
                muminus = float(parts[4])
                result[no] = (muplus, muminus)
            except (ValueError, IndexError):
                continue
    return result


def _read_multicolumn_2d(path, n_cols):
    """
    Parse T-Cross output where samepage puts all regions/z-levels as columns.

    Format (per data row):
        x_lower  x_upper  val1  err1  val2  err2  ...  valN  errN

    Particle-type blocks are separated by comment lines:
      - KE spectrum:  "# no. = 1  muon+"  (particle type on same line as no.)
      - Angle tally:  "# muon+"            (particle type on its own line)

    Returns {particle_type: (x_lower[M], x_upper[M], data[M, n_cols])}
    """
    blocks = {}
    current_part = None
    x_lowers, x_uppers, data_rows = [], [], []

    if not os.path.exists(path):
        print(f"  [WARN] not found: {path}")
        return blocks

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith('#'):
                has_plus  = 'muon+' in line
                has_minus = 'muon-' in line
                if has_plus and not has_minus:
                    new_part = 'muon+'
                elif has_minus and not has_plus:
                    new_part = 'muon-'
                else:
                    new_part = None

                if new_part and new_part != current_part:
                    if current_part is not None and x_lowers:
                        blocks[current_part] = (
                            np.array(x_lowers),
                            np.array(x_uppers),
                            np.array(data_rows))
                    current_part = new_part
                    x_lowers, x_uppers, data_rows = [], [], []
                continue

            if not _is_data_line(raw) or current_part is None:
                continue

            parts = line.split()
            if len(parts) < 2 + 2 * n_cols:
                continue
            try:
                x_lowers.append(float(parts[0]))
                x_uppers.append(float(parts[1]))
                data_rows.append([float(parts[2 + 2 * i]) for i in range(n_cols)])
            except (ValueError, IndexError):
                continue

    if current_part is not None and x_lowers:
        blocks[current_part] = (
            np.array(x_lowers),
            np.array(x_uppers),
            np.array(data_rows))

    return blocks


def _read_reg_deposit(path):
    """
    Read T-Deposit with mesh=reg, output=deposit, unit=3.
    Actual column layout:  num  reg  volume  muon+  err  muon-  err  all  err

    NOTE: With axis=reg and ne=1, PHITS counts deposit EVENTS per source particle,
    NOT energy in MeV. The 'all' value equals the fraction of source particles
    that enter (traverse or stop in) each layer.  Volume is always 1.0.

    Returns {region (1-6): traverse_fraction_per_source}
    """
    result = {}
    if not os.path.exists(path):
        print(f"  [WARN] not found: {path}")
        return result
    with open(path) as fh:
        for raw in fh:
            if not _is_data_line(raw):
                continue
            parts = raw.split()
            if len(parts) < 8:
                continue
            try:
                reg = int(parts[1])    # region number (col 1)
                val = float(parts[7])  # 'all' particles (col 7)
                result[reg] = val
            except (ValueError, IndexError):
                continue
    return result


# ═══════════════════════════════════════════════════════════════════════════════
#  Run-log parser
# ═══════════════════════════════════════════════════════════════════════════════

def parse_run_log(path="muon_rock.out"):
    n_source = None
    timing_s = None

    # Try muon_rock.out first, then phits.out
    for try_path in [path, "phits.out", "batch.out"]:
        if not os.path.exists(try_path):
            continue
        with open(try_path) as fh:
            content = fh.read()

        if n_source is None:
            m = re.search(r'maxcas\s*=\s*(\d+)', content)
            if m:
                maxcas = int(m.group(1))
                m2 = re.search(r'maxbch\s*=\s*(\d+)', content)
                maxbch = int(m2.group(1)) if m2 else 1
                n_source = maxcas * maxbch

        if timing_s is None:
            # "total cpu time = 674.29" (seconds)
            m = re.search(r'total cpu time\s*=\s*([\d.]+)', content, re.IGNORECASE)
            if m:
                timing_s = float(m.group(1))
            else:
                # "cpu time = 11 m. 14.14 s."
                m = re.search(r'cpu time\s*=\s*(?:(\d+)\s*m\.\s*)?([\d.]+)\s*s', content, re.IGNORECASE)
                if m:
                    minutes = int(m.group(1)) if m.group(1) else 0
                    timing_s = minutes * 60 + float(m.group(2))
                else:
                    m = re.search(r'elapsed time\s*=\s*([\d.]+)', content, re.IGNORECASE)
                    if m:
                        timing_s = float(m.group(1))

    return n_source, timing_s


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="Parse PHITS output → benchmark CSVs")
    ap.add_argument("--n",      type=int, default=None, help="Override n_source")
    ap.add_argument("--outdir", default=".",            help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=" * 60)
    print("  PHITS output reader → benchmark CSVs")
    print("=" * 60)

    # ── n_source and timing ────────────────────────────────────────────────────
    n_source, timing_s = parse_run_log("muon_rock.out")
    if args.n:
        n_source = args.n
    if n_source is None and os.path.exists("muon_rock.inp"):
        with open("muon_rock.inp") as fh:
            inp = fh.read()
        m  = re.search(r'^\s*maxcas\s*=\s*(\d+)', inp, re.MULTILINE)
        m2 = re.search(r'^\s*maxbch\s*=\s*(\d+)',  inp, re.MULTILINE)
        if m:
            n_source = int(m.group(1)) * (int(m2.group(1)) if m2 else 1)
    if n_source is None:
        n_source = 10003
        print(f"  [WARN] n_source unknown — defaulting to {n_source}")
    else:
        print(f"  n_source  = {n_source:,}")
    if timing_s:
        print(f"  Elapsed   = {timing_s:.1f} s")

    # ── 1. Transmission (mesh=reg, axis=reg) ──────────────────────────────────
    print("\n[1] Transmission: out_muon_transmission.out")
    trans_raw = _read_reg_cross_1d("out_muon_transmission.out")
    trans_at_depth = {}
    for pair_no, (muplus, muminus) in trans_raw.items():
        depth = REG_TO_DEPTH.get(pair_no)
        if depth is None:
            continue
        total_curr   = muplus + muminus       # [1/cm²/src]
        transmission = total_curr * SLAB_AREA_CM2   # dimensionless fraction
        trans_at_depth[depth] = max(0.0, transmission)
        print(f"    pair {pair_no}  depth={depth:6d} cm "
              f"mu+={muplus:.4e}  mu-={muminus:.4e}  T={transmission*100:.3f}%")

    # ── 2. Exit KE spectrum (mesh=reg, axis=eng, samepage=reg → 6 region cols) ─
    print("\n[2] KE spectrum: out_muon_KE_spectrum.out")
    ke_blocks = _read_multicolumn_2d("out_muon_KE_spectrum.out", n_cols=6)
    mean_ke_at_depth = {}
    std_ke_at_depth  = {}

    if ke_blocks:
        ref_xl, ref_xu, _ = next(iter(ke_blocks.values()))
        combined_ke = None
        for pt in ('muon+', 'muon-'):
            if pt in ke_blocks:
                _, _, data = ke_blocks[pt]
                combined_ke = data if combined_ke is None else combined_ke + data

        if combined_ke is not None:
            de  = ref_xu - ref_xl
            e_c = 0.5 * (ref_xl + ref_xu)
            for col_idx, depth in KE_COL_TO_DEPTH.items():
                if col_idx >= combined_ke.shape[1]:
                    continue
                spec  = combined_ke[:, col_idx]   # [1/cm²/(MeV/n)/src]
                w     = spec * de                  # weight [1/cm²/src] per bin
                w_sum = w.sum()
                if w_sum > 0:
                    mean_MeV = (e_c * w).sum() / w_sum
                    var_MeV  = ((e_c - mean_MeV)**2 * w).sum() / w_sum
                    mean_ke_at_depth[depth] = mean_MeV / 1000.0
                    std_ke_at_depth[depth]  = np.sqrt(max(var_MeV, 0)) / 1000.0
                    print(f"    depth={depth:6d} cm  "
                          f"mean KE = {mean_ke_at_depth[depth]:.3f} ± "
                          f"{std_ke_at_depth[depth]:.3f} GeV")
                else:
                    print(f"    depth={depth:6d} cm  zero spectrum")
    else:
        print("  [WARN] no blocks parsed")

    # ── 3. Angle distribution (mesh=r-z, axis=cos, samepage=z → 7 z-cols) ─────
    print("\n[3] Angle: out_muon_angle.out")
    ang_blocks = _read_multicolumn_2d("out_muon_angle.out", n_cols=7)
    mean_angle_at_depth = {}
    std_angle_at_depth  = {}

    if ang_blocks:
        ref_al, ref_au, _ = next(iter(ang_blocks.values()))
        combined_ang = None
        for pt in ('muon+', 'muon-'):
            if pt in ang_blocks:
                _, _, data = ang_blocks[pt]
                combined_ang = data if combined_ang is None else combined_ang + data

        if combined_ang is not None:
            da    = ref_au - ref_al
            cos_c = 0.5 * (ref_al + ref_au)
            # iangform=3: cosine is from the +z axis.
            # Downward muons have cos_z ≈ -1. Scattering angle from -z:
            #   theta_scatter = arccos(-cos_z)
            scatter_cos_c = np.clip(-cos_c, -1.0, 1.0)
            scatter_deg_c = np.degrees(np.arccos(scatter_cos_c))

            for col_idx, depth in ANGLE_COL_TO_DEPTH.items():
                if col_idx >= combined_ang.shape[1]:
                    continue
                spec  = combined_ang[:, col_idx]   # [1/cm²/sr/src]
                w     = spec * da
                w_sum = w.sum()
                if w_sum > 0:
                    mean_a = (scatter_deg_c * w).sum() / w_sum
                    var_a  = ((scatter_deg_c - mean_a)**2 * w).sum() / w_sum
                    mean_angle_at_depth[depth] = mean_a
                    std_angle_at_depth[depth]  = np.sqrt(max(var_a, 0))
                    mean_cos = (cos_c * w).sum() / w_sum
                    print(f"    depth={depth:6d} cm  "
                          f"mean cos_z={mean_cos:.4f}  "
                          f"mean scatter={mean_a:.3f}° ± "
                          f"{std_angle_at_depth[depth]:.3f}°")
                else:
                    print(f"    depth={depth:6d} cm  zero flux")
    else:
        print("  [WARN] no blocks parsed")

    # ── 4. T-Deposit (traversal counts, NOT energy in MeV) ────────────────────
    print("\n[4] T-Deposit traversal fractions: out_muon_edep_total.out")
    print("    (NOTE: unit=3 with output=deposit counts traverse events/source,")
    print("     not energy. Values ≈ fraction of muons entering each layer.)")
    dep_all = _read_reg_deposit("out_muon_edep_total.out")
    for reg in sorted(dep_all):
        depth = DEP_REG_TO_DEPTH.get(reg, reg)
        print(f"    layer {reg}  (layer bottom at {depth:6d} cm)"
              f"  entry fraction = {dep_all[reg]:.4f}")

    # ── Write phits_summary.csv ────────────────────────────────────────────────
    # Note on MeanZenithAngle_deg: this is the mean angle of surviving muons
    # from the downward (-z) direction, including the initial beam divergence.
    # It is NOT the MCS scattering angle from the initial muon direction
    # (Geant4's AngleScat_Deg). The two quantities are not directly comparable:
    # PHITS shows zenith-angle hardening with depth; Geant4 shows MCS growth.
    summary_path = os.path.join(args.outdir, "phits_summary.csv")
    with open(summary_path, "w") as fout:
        fout.write("DepthCm,MWE,N_transmitted,Transmission_%,"
                   "MeanExitKE_GeV,StdExitKE_GeV,"
                   "MeanZenithAngle_deg,StdZenithAngle_deg\n")
        for d in SCORING_DEPTHS_CM:
            mwe     = d * ROCK_DENSITY / 100.0
            tr      = trans_at_depth.get(d, float("nan"))
            n_trans = tr * n_source if tr == tr else float("nan")
            tr_pct  = tr * 100.0   if tr == tr else float("nan")
            ke      = mean_ke_at_depth.get(d, float("nan"))
            ke_std  = std_ke_at_depth.get(d, float("nan"))
            ang     = mean_angle_at_depth.get(d, float("nan"))
            ang_std = std_angle_at_depth.get(d, float("nan"))
            fout.write(f"{d},{mwe:.3f},{n_trans:.1f},{tr_pct:.4f},"
                       f"{ke:.4f},{ke_std:.4f},"
                       f"{ang:.4f},{ang_std:.4f}\n")
    print(f"\nSaved {summary_path}")

    # ── Write phits_timing.txt ─────────────────────────────────────────────────
    timing_path = os.path.join(args.outdir, "phits_timing.txt")
    with open(timing_path, "w") as fout:
        fout.write(f"Events  : {n_source}\n")
        if timing_s:
            fout.write(f"Elapsed : {timing_s:.2f} s\n")
            fout.write(f"Rate    : {n_source / timing_s:.1f} evt/s\n")
    print(f"Saved {timing_path}")

    print("\n" + "=" * 60)
    print("  Done.  Run benchmark comparison from the parent directory:")
    print("    python3 benchmark_analysis.py --geant4 outputs/ --phits phits/")
    print("=" * 60)


if __name__ == "__main__":
    main()
