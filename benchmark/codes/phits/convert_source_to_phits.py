#!/usr/bin/env python3
"""
convert_source_to_phits.py
Convert the Geant4 muon source file (muons_geant4.txt) into a PHITS
ASCII dump file (s-type=17) that muon_rock.inp can read directly.

Usage:
    python3 convert_source_to_phits.py [input_file] [output_file]

Defaults:
    input_file  = ../muons_geant4.txt
    output_file = muon_source_phits.dat

Geant4 source format (space-separated, comment lines start with #):
    PDG  x_mm  y_mm  z_mm  px_MeV  py_MeV  pz_MeV  Ekin_MeV

PHITS dump format (8 columns, ASCII):
    kf   x_cm   y_cm   z_cm   u   v   w   ek_MeV
    kf   = PDG particle code  (13 = µ−,  -13 = µ+)
    x,y  = lateral position in cm (converted from mm)
    z    = 0.0 cm  (all muons placed at the rock entry surface)
    u,v,w= direction cosines  (NOTE: z-axis is FLIPPED relative to Geant4)
           Geant4: pz < 0 means downward (−z is down)
           PHITS:  w  > 0 means downward (+z is down in muon_rock.inp)
           Conversion: w = −pz / |p|
    ek   = kinetic energy in MeV  (same units, no conversion needed)

Coordinate check:
    Geant4 muon going straight down: pz = −P, px = py = 0
    → u = 0,  v = 0,  w = −(−P)/P = +1   (downward in PHITS)  ✓

Rock slab clamping:
    The Geant4 source can place muons outside ±25 m (the slab XY boundary).
    Geant4's PrimaryGeneratorAction clamps XY to ±24 990 mm automatically.
    This script applies the same clamp: x,y are capped at ±2499.0 cm.
"""

import sys
import os
import math

# ── Configurable limits (must match muon_rock.inp geometry) ─────────────────
SLAB_HALF_XY_CM = 2499.0    # cap at 0.1 cm inside the wall (matching Geant4)
Z_ENTRY_CM      = 0.0       # all muons start at the entry surface

def convert(inpath, outpath):
    n_read = n_written = n_clamped = 0

    # PHITS s-type=17 ASCII dump: purely numeric, NO header/comment lines.
    # Write a separate info file alongside the dump so the format is documented.
    info_path = outpath.replace(".dat", "_info.txt")
    with open(info_path, "w") as fi:
        fi.write(f"PHITS ASCII dump file: {os.path.basename(outpath)}\n")
        fi.write(f"Source: {os.path.basename(inpath)}\n")
        fi.write("Columns (10 items, matches: dump = -10  1 2 3 4 5 6 7 8 9 10):\n")
        fi.write("  1: kf   PDG particle code (13=mu-  -13=mu+)\n")
        fi.write("  2: x    cm  (from mm, clamped to slab)\n")
        fi.write("  3: y    cm\n")
        fi.write("  4: z    cm  (always 0 = rock entry surface)\n")
        fi.write("  5: u    direction cosine x\n")
        fi.write("  6: v    direction cosine y\n")
        fi.write("  7: w    direction cosine z  (w<0 = downward = into rock, matches z<0 geometry)\n")
        fi.write("  8: ek   kinetic energy (MeV)\n")
        fi.write("  9: wgt  statistical weight (always 1.0)\n")
        fi.write(" 10: time arrival time (always 0.0 ns)\n")

    with open(inpath, "r") as fin, open(outpath, "w") as fout:
        for raw in fin:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("$"):
                continue

            parts = line.split()
            if len(parts) < 8:
                continue
            n_read += 1

            try:
                pdg    = int(parts[0])
                x_mm   = float(parts[1])
                y_mm   = float(parts[2])
                # parts[3] = z_mm  (ignored; all muons placed at entry surface)
                px     = float(parts[4])   # MeV/c
                py     = float(parts[5])
                pz     = float(parts[6])
                ek     = float(parts[7])   # MeV
            except ValueError:
                continue

            # Position: mm → cm, clamp to slab
            x_cm = x_mm / 10.0
            y_cm = y_mm / 10.0
            x_cm = max(-SLAB_HALF_XY_CM, min(SLAB_HALF_XY_CM, x_cm))
            y_cm = max(-SLAB_HALF_XY_CM, min(SLAB_HALF_XY_CM, y_cm))
            if abs(x_mm / 10.0) > SLAB_HALF_XY_CM or abs(y_mm / 10.0) > SLAB_HALF_XY_CM:
                n_clamped += 1

            # Direction cosines: keep Geant4 sign convention.
            # Geant4: pz < 0 for downward muons (+z points up).
            # PHITS geometry (muon_rock.inp): z < 0 is deeper into rock,
            # so w < 0 also means downward — same convention, NO sign flip.
            p_mag = math.sqrt(px*px + py*py + pz*pz)
            if p_mag < 1e-9:
                continue  # skip zero-momentum particles
            u = px / p_mag
            v = py / p_mag
            w = pz / p_mag    # pz<0 for downward → w<0 downward in z<0 geometry

            fout.write(
                f"  {pdg:4d}  {x_cm:12.6f}  {y_cm:12.6f}  {Z_ENTRY_CM:8.4f}"
                f"  {u:12.8f}  {v:12.8f}  {w:12.8f}  {ek:14.6f}"
                f"  {1.0:8.4f}  {0.0:8.4f}\n"
            )
            n_written += 1

    return n_read, n_written, n_clamped


def main():
    inpath  = sys.argv[1] if len(sys.argv) > 1 else "../muons_geant4.txt"
    outpath = sys.argv[2] if len(sys.argv) > 2 else "muon_source_phits.dat"

    if not os.path.exists(inpath):
        print(f"[ERROR] Input file not found: {inpath}")
        print(f"  Usage: python3 {os.path.basename(__file__)} [input] [output]")
        sys.exit(1)

    print(f"Converting: {inpath}  →  {outpath}")
    n_read, n_written, n_clamped = convert(inpath, outpath)

    print(f"  Lines read    : {n_read:,}")
    print(f"  Particles out : {n_written:,}")
    if n_clamped:
        print(f"  XY clamped    : {n_clamped:,} (outside ±{SLAB_HALF_XY_CM} cm, clipped to slab boundary)")
    print(f"  Done → {outpath}")
    print()
    print("  In muon_rock.inp, set:  maxcas = " + str(n_written))
    print("  (maxcas must match the number of lines in the dump file)")


if __name__ == "__main__":
    main()
