#!/usr/bin/env python3
"""
ucmuon_to_phits.py  —  UCMuon → PHITS s-type=17 dump converter
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Converts UCMuon output files to PHITS ASCII dump format (s-type=17).
Works on both generator output (ucmuon_selected.dat) and transport
output (ucmuon_underground.dat).

Usage
-----
  python3 ucmuon_to_phits.py <input.dat> [output.dat] [--mode gen|transport]

  mode=gen        (default) reads ucmuon_selected.dat / ucmuon_surface.dat
                  13-col or 14-col generator format
                  → uses surface position and direction

  mode=transport  reads ucmuon_underground.dat (18-col transport output)
                  → uses underground position and direction, alive==1 only

Output format
-------------
  PHITS s-type=17 ASCII dump, 10 columns per muon:
    kf   x[cm]  y[cm]  z[cm]  u  v  w  Ekin[MeV]  wt  time[ns]

  kf codes (PDG): mu- = 13,  mu+ = -13

  Use in PHITS input:
    [Source]
    s-type = 17
    file   = ucmuon_phits.dat
    dump   = -10
    1 2 3 4 5 6 7 8 9 10

Column layout reference
-----------------------
  Generator 13-col: id  x  y  z  p  px  py  pz  theta  phi  E  charge  det_mask
  Generator 14-col: id  x  y  z  p  px  py  pz  theta  phi  E  charge  hit_flag  det_mask
  Transport 18-col: id  x_srf  y_srf  z_srf  E_srf  theta_srf  phi_srf  charge
                    alive  x_ug  y_ug  z_ug  E_ug  cx_ug  cy_ug  cz_ug
                    theta_ug  phi_ug
"""

import sys
import os
import math
import argparse

MUON_MASS_GEV = 0.105658370   # GeV/c²
MUON_MASS_MEV = 105.658370    # MeV/c²


def detect_format(line):
    """Return (ncols, mode) from the first data line."""
    parts = line.split()
    n = len(parts)
    if n == 13:
        return 13, 'gen'
    elif n == 14:
        return 14, 'gen'
    elif n == 18:
        return 18, 'transport'
    else:
        raise ValueError(f"Unrecognised column count: {n} (expected 13, 14, or 18)")


def direction_from_momentum(px, py, pz, p):
    """Return unit direction cosines (u, v, w) from momentum components."""
    if p > 0.0:
        return px / p, py / p, pz / p
    return 0.0, 0.0, -1.0   # fallback: straight down


def direction_from_angles(theta, phi):
    """Return unit direction cosines from zenith/azimuth angles [rad]."""
    st = math.sin(theta)
    return st * math.cos(phi), st * math.sin(phi), -math.cos(theta)


def energy_to_kf(charge, e_gev):
    """Return (kf, ekin_MeV). charge: +1=mu+, -1=mu-."""
    kf = -13 if charge == 1 else 13          # PDG: mu+=−13, mu−=+13
    ekin_mev = max(0.0, (e_gev - MUON_MASS_GEV) * 1000.0)
    return kf, ekin_mev


def fmt(v):
    """PHITS dump-a format: 1pd24.15 (D-style exponent, 24 wide, 15 decimals)."""
    s = f"{v:24.15E}"
    # Python uses 'E', PHITS expects 'D' — both are accepted by most PHITS versions
    return s


def convert(infile, outfile, mode):
    n_written = 0
    n_skipped = 0

    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for raw in fin:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            ncols = len(parts)

            if mode == 'auto':
                mode = detect_format(line)[1]

            # ----------------------------------------------------------------
            # Generator format (13 or 14 columns)
            # ----------------------------------------------------------------
            if mode == 'gen':
                if ncols not in (13, 14):
                    continue
                try:
                    x     = float(parts[1])
                    y     = float(parts[2])
                    z     = float(parts[3])
                    p     = float(parts[4])
                    px    = float(parts[5])
                    py    = float(parts[6])
                    pz    = float(parts[7])
                    theta = float(parts[8])
                    phi   = float(parts[9])
                    e_gev = float(parts[10])
                    charge= int(parts[11])
                except (ValueError, IndexError):
                    n_skipped += 1
                    continue

                u, v, w = direction_from_momentum(px, py, pz, p)
                kf, ekin = energy_to_kf(charge, e_gev)

            # ----------------------------------------------------------------
            # Transport format (18 columns) — only alive muons
            # ----------------------------------------------------------------
            elif mode == 'transport':
                if ncols != 18:
                    continue
                try:
                    charge = int(parts[7])
                    alive  = int(parts[8])
                    x      = float(parts[9])
                    y      = float(parts[10])
                    z      = float(parts[11])
                    e_gev  = float(parts[12])
                    cx     = float(parts[13])
                    cy     = float(parts[14])
                    cz     = float(parts[15])
                except (ValueError, IndexError):
                    n_skipped += 1
                    continue

                if alive != 1:
                    n_skipped += 1
                    continue

                u, v, w = cx, cy, cz
                # Normalise direction (should already be unit but just in case)
                norm = math.sqrt(u*u + v*v + w*w)
                if norm > 1e-14:
                    u, v, w = u/norm, v/norm, w/norm

                kf, ekin = energy_to_kf(charge, e_gev)

            else:
                raise ValueError(f"Unknown mode: {mode}")

            # Write PHITS record: kf x y z u v w ekin wt time
            fout.write(
                fmt(float(kf)) +
                fmt(x) + fmt(y) + fmt(z) +
                fmt(u) + fmt(v) + fmt(w) +
                fmt(ekin) +
                fmt(1.0) +   # weight = 1
                fmt(0.0) +   # time = 0 ns
                '\n'
            )
            n_written += 1

    return n_written, n_skipped


def main():
    parser = argparse.ArgumentParser(
        description='Convert UCMuon output to PHITS s-type=17 dump')
    parser.add_argument('input',  help='Input file (ucmuon_selected.dat or ucmuon_underground.dat)')
    parser.add_argument('output', nargs='?', help='Output PHITS file (default: <input_stem>_phits.dat)')
    parser.add_argument('--mode', choices=['gen', 'transport', 'auto'],
                        default='auto',
                        help='Input format: gen=generator, transport=18-col, auto=detect')
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file not found: {args.input}")
        sys.exit(1)

    if args.output is None:
        stem = os.path.splitext(args.input)[0]
        args.output = stem + '_phits.dat'

    print(f"UCMuon → PHITS converter")
    print(f"  Input  : {args.input}  (mode={args.mode})")
    print(f"  Output : {args.output}")

    n_written, n_skipped = convert(args.input, args.output, args.mode)

    print(f"  Written: {n_written:,} muons")
    if n_skipped:
        print(f"  Skipped: {n_skipped:,} (stopped muons or bad lines)")
    print()
    print("Use in PHITS input section:")
    print("  [Source]")
    print("  s-type = 17")
    print(f"  file   = {args.output}")
    print("  dump   = -10")
    print("  1 2 3 4 5 6 7 8 9 10")


if __name__ == '__main__':
    main()
