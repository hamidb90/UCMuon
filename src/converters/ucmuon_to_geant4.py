#!/usr/bin/env python3
"""
ucmuon_to_geant4.py  —  UCMuon → Geant4 primary-vertex file converter
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Converts UCMuon output files to the UCMuon ASCII Geant4 source format.
Works on both generator output and transport output.

Usage
-----
  python3 ucmuon_to_geant4.py <input.dat> [output.txt] [--mode gen|transport|auto]
                              [--electrons]

  mode=auto       (default) picks gen or transport from the column count
  mode=gen        reads 13- or 14-col generator format
                  → surface positions and directions
  mode=transport  reads 18-col transport output (alive muons only)
                  → underground positions and directions
  --electrons     the input is a spectrum-8 (cosmic e±) generator file:
                  PDG ±11 and the electron mass. A generator file has no
                  particle column, so this has to be said; a row whose own E
                  and p contradict the particle being written stops the run.

Output format (UCMuon ASCII)
----------------------------
  Header comments (lines starting with #), then one muon per line:
    PDG  x[mm]  y[mm]  z[mm]  px[MeV/c]  py[MeV/c]  pz[MeV/c]  Ekin[MeV]

  PDG codes: mu- = 13,  mu+ = -13   (--electrons: e- = 11,  e+ = -11)

Unit conversions applied
------------------------
  Positions  : cm  → mm  (×10)
  Momenta    : GeV/c → MeV/c  (×1000)
  Energy     : GeV  → MeV kinetic  (×1000 − 105.658, or − 0.511 for e±)

Read in Geant4 (PrimaryGeneratorAction.cc example)
---------------------------------------------------
  while (std::getline(fFile, line)) {
      if (line.empty() || line[0] == '#') continue;
      G4int pdg; G4double x,y,z,px,py,pz,Ekin;
      std::istringstream(line) >> pdg >> x >> y >> z >> px >> py >> pz >> Ekin;
      fGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(pdg));
      fGun->SetParticlePosition(G4ThreeVector(x*mm, y*mm, z*mm));
      fGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz).unit());
      fGun->SetParticleEnergy(Ekin * MeV);
      fGun->GeneratePrimaryVertex(event); break;
  }

Column layout reference
-----------------------
  Generator 13-col: id  x  y  z  p  px  py  pz  theta  phi  E  charge  det_mask
  Generator 14-col: id  x  y  z  p  px  py  pz  theta  phi  E  charge  hit_flag  det_mask
  Transport 18-col: id  xs ys zs Es theta_s phi_s charge alive x y z E cx cy cz theta phi
"""

import sys
import os
import math
import argparse

MUON_MASS_MEV = 105.658370   # MeV/c²
MUON_MASS_GEV = 0.105658370  # GeV/c²
ELECTRON_MASS_MEV = 0.51099895   # MeV/c²


class SpeciesError(ValueError):
    """A row's own E and p say it is not the particle being written."""


def check_species(p_gev, e_gev, electrons, event_id):
    """Below 20 GeV the file's six decimals give sqrt(E^2 - p^2) to better than
    1 %, and a muon (0.106 GeV) and an electron (0.0005 GeV) are far apart; above
    it the check is skipped."""
    if e_gev >= 20.0:
        return
    m = math.sqrt(max(e_gev * e_gev - p_gev * p_gev, 0.0))
    if electrons and m > 0.05:
        raise SpeciesError(f"event {event_id}: E = {e_gev} GeV, p = {p_gev} GeV/c give a mass "
                           f"of {m:.4f} GeV, a muon, not an electron. Drop --electrons.")
    if not electrons and m < 0.05:
        raise SpeciesError(f"event {event_id}: E = {e_gev} GeV, p = {p_gev} GeV/c give a mass "
                           f"of {m:.4f} GeV, not a muon: a spectrum-8 (cosmic e±) file? "
                           f"Use --electrons.")


def direction_from_momentum(px, py, pz, p):
    if p > 0.0:
        return px / p, py / p, pz / p
    return 0.0, 0.0, -1.0


def convert(infile, outfile, mode, electrons=False):
    n_written = 0
    n_skipped = 0
    base = 11 if electrons else 13
    mass_mev = ELECTRON_MASS_MEV if electrons else MUON_MASS_MEV

    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        fout.write("# UCMuon Geant4 source file — UCLouvain Muography Group\n")
        fout.write(f"# Source: {os.path.basename(infile)}\n")
        fout.write("# Format: PDG  x[mm]  y[mm]  z[mm]  px[MeV/c]  py[MeV/c]  pz[MeV/c]  Ekin[MeV]\n")
        fout.write("# PDG: 11 = e-   -11 = e+\n" if electrons else "# PDG: 13 = mu-   -13 = mu+\n")

        for raw in fin:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            ncols = len(parts)

            # ── Generator format (13, 14, or 15 columns) ─────────────────────
            if mode in ('gen', 'auto') and ncols in (13, 14, 15):
                try:
                    x_cm  = float(parts[1])
                    y_cm  = float(parts[2])
                    z_cm  = float(parts[3])
                    p_gev = float(parts[4])
                    px_g  = float(parts[5])
                    py_g  = float(parts[6])
                    pz_g  = float(parts[7])
                    e_gev = float(parts[10])
                    charge = int(parts[11])
                except (ValueError, IndexError):
                    n_skipped += 1
                    continue

                check_species(p_gev, e_gev, electrons, parts[0])
                pdg = base if charge < 0 else -base
                # Use stored momentum components — more precise than recomputing from angles
                px_mev = px_g * 1000.0
                py_mev = py_g * 1000.0
                pz_mev = pz_g * 1000.0
                ekin   = max(0.0, e_gev * 1000.0 - mass_mev)

            # ── Transport format (18 columns, alive==1 only) ──────────────────
            elif mode in ('transport', 'auto') and ncols == 18:
                if electrons:
                    raise SpeciesError("an 18-column transport file holds muons only; "
                                       "--electrons applies to generator files")
                try:
                    charge = int(parts[7])
                    alive  = int(parts[8])
                    x_cm   = float(parts[9])
                    y_cm   = float(parts[10])
                    z_cm   = float(parts[11])
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

                # Normalise direction (safety)
                norm = math.sqrt(cx*cx + cy*cy + cz*cz)
                if norm > 1e-14:
                    cx /= norm; cy /= norm; cz /= norm

                pdg = 13 if charge < 0 else -13
                p_mev = math.sqrt(max(0.0, e_gev**2 - MUON_MASS_GEV**2)) * 1000.0
                px_mev = cx * p_mev
                py_mev = cy * p_mev
                pz_mev = cz * p_mev
                ekin   = max(0.0, e_gev * 1000.0 - MUON_MASS_MEV)

            else:
                n_skipped += 1
                continue

            # Positions: cm → mm
            x_mm = x_cm * 10.0
            y_mm = y_cm * 10.0
            z_mm = z_cm * 10.0

            fout.write(
                f"{pdg:4d}  {x_mm:12.4f}  {y_mm:12.4f}  {z_mm:12.4f}"
                f"  {px_mev:14.6f}  {py_mev:14.6f}  {pz_mev:14.6f}  {ekin:14.6f}\n"
            )
            n_written += 1

    return n_written, n_skipped


def main():
    parser = argparse.ArgumentParser(
        description='Convert UCMuon output to Geant4 ASCII primary-vertex format')
    parser.add_argument('input',  help='Input file (UCMuon generator or transport output)')
    parser.add_argument('output', nargs='?',
                        help='Output Geant4 file (default: <input_stem>_geant4.txt)')
    parser.add_argument('--mode', choices=['gen', 'transport', 'auto'],
                        default='auto',
                        help='Input format: gen=generator (13/14-col), '
                             'transport=18-col, auto=detect')
    parser.add_argument('--electrons', action='store_true',
                        help='spectrum-8 (cosmic e±) generator file: PDG ±11, electron mass')
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.output is None:
        stem = os.path.splitext(args.input)[0]
        args.output = stem + '_geant4.txt'

    print(f"UCMuon → Geant4 converter")
    print(f"  Input  : {args.input}  (mode={args.mode})")
    print(f"  Output : {args.output}")
    print(f"  Units  : positions cm→mm (×10), momenta GeV/c→MeV/c (×1000), Ekin MeV")

    try:
        n_written, n_skipped = convert(args.input, args.output, args.mode, args.electrons)
    except SpeciesError as exc:
        os.remove(args.output)
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    print(f"  Written: {n_written:,} {'electrons' if args.electrons else 'muons'}")
    if n_skipped:
        print(f"  Skipped: {n_skipped:,} (stopped muons or unrecognised lines)")
    print()
    print("Geant4 PrimaryGeneratorAction snippet:")
    print("  fGun->SetParticlePosition(G4ThreeVector(x*mm, y*mm, z*mm));")
    print("  fGun->SetParticleEnergy(Ekin * MeV);")


if __name__ == '__main__':
    main()
