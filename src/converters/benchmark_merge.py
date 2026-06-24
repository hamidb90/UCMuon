#!/usr/bin/env python3
"""
benchmark_merge.py  —  Merge per-energy UCMuon benchmark source files
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Scans an output directory for  muons_surface_<E>.dat,  muons_for_phits_<E>.dat,
and  muons_geant4_<E>.txt  files, sorts them by energy, and writes three merged
benchmark source files ready for transport runs:

  benchmark_surface.dat   —  UCMuon 14-col format  (input to UCMuon transport)
  benchmark_phits.dat     —  PHITS s-type=17 dump  (input to PHITS)
  benchmark_geant4.txt    —  UCMuon ASCII Geant4   (input to Geant4)

Usage
-----
  python3 benchmark_merge.py [output_dir] [--outdir DIR]

  output_dir   directory containing the per-energy files  (default: ./output)
  --outdir     where to write the three merged files       (default: same as output_dir)

Example
-------
  python3 src/converters/benchmark_merge.py output/
"""

import os
import re
import sys
import argparse


# ── helpers ──────────────────────────────────────────────────────────────────

def energy_from_name(path):
    """Extract numeric energy (GeV) from filenames like muons_surface_300.dat."""
    m = re.search(r'_(\d+)\.', os.path.basename(path))
    return int(m.group(1)) if m else 0


def energy_label(e_gev):
    if e_gev >= 1000:
        return f"{e_gev // 1000} TeV"
    return f"{e_gev} GeV"


def find_files(directory, pattern):
    """Return files matching pattern, sorted by embedded energy number."""
    regex = re.compile(pattern)
    files = [os.path.join(directory, f)
             for f in os.listdir(directory)
             if regex.match(f)]
    return sorted(files, key=energy_from_name)


# ── surface merge (UCMuon 14-col, renumber EventIDs) ─────────────────────────

def merge_surface(files, outpath):
    print(f"\n  Merging {len(files)} surface files → {outpath}")
    event_id = 0
    batch_log = []

    with open(outpath, 'w') as fout:
        fout.write(
            "# UCMuon benchmark — merged surface source\n"
            "# Format: EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV"
            "  theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask\n"
            "#\n"
        )
        # placeholder lines for batch table — will be filled after writing
        table_pos = fout.tell()
        # write a fixed-width block we can seek back to
        PLACEHOLDER = "# {:<72s}\n"
        MAX_BATCHES = 10
        for _ in range(MAX_BATCHES):
            fout.write(PLACEHOLDER.format(""))
        fout.write("#\n")

        for fpath in files:
            e_gev  = energy_from_name(fpath)
            first  = event_id + 1
            n_rows = 0
            with open(fpath) as fin:
                for raw in fin:
                    line = raw.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    event_id += 1
                    parts[0] = str(event_id)
                    fout.write('  '.join(parts) + '\n')
                    n_rows += 1
            batch_log.append((e_gev, first, event_id, n_rows))
            print(f"    {energy_label(e_gev):>8s}  {n_rows:>8,} muons  "
                  f"EventID {first}–{event_id}")

    # Rewrite the placeholder block in-place via full file rewrite (simplest)
    with open(outpath, 'r') as f:
        content = f.read()

    table_lines = []
    for e_gev, first, last, n in batch_log:
        table_lines.append(
            f"# {energy_label(e_gev):>8s}  {n:>8,} muons  EventID {first}–{last}"
        )
    # Pad to MAX_BATCHES lines so total file structure is identical
    while len(table_lines) < MAX_BATCHES:
        table_lines.append("#")

    old_block = "\n".join([PLACEHOLDER.format("").rstrip('\n')] * MAX_BATCHES)
    new_block  = "\n".join(table_lines)
    content = content.replace(old_block, new_block, 1)

    with open(outpath, 'w') as f:
        f.write(content)

    total = sum(n for *_, n in batch_log)
    print(f"    Total: {total:,} muons")
    return total


# ── PHITS merge (no header — pure data concatenation) ────────────────────────

def merge_phits(files, outpath):
    print(f"\n  Merging {len(files)} PHITS files → {outpath}")
    total = 0
    with open(outpath, 'w') as fout:
        for fpath in files:
            e_gev = energy_from_name(fpath)
            n_rows = 0
            with open(fpath) as fin:
                for raw in fin:
                    line = raw.strip()
                    if not line or line.startswith('#'):
                        continue
                    fout.write(raw)   # keep original Fortran D-format spacing
                    n_rows += 1
            total += n_rows
            print(f"    {energy_label(e_gev):>8s}  {n_rows:>8,} muons")
    print(f"    Total: {total:,} muons")
    return total


# ── Geant4 merge (one header block, then data from all files) ─────────────────

def merge_geant4(files, outpath):
    print(f"\n  Merging {len(files)} Geant4 files → {outpath}")
    total = 0
    with open(outpath, 'w') as fout:
        fout.write(
            "# UCMuon benchmark — merged Geant4 source\n"
            "# Format: PDG  x[mm]  y[mm]  z[mm]  px[MeV/c]  py[MeV/c]  pz[MeV/c]  Ekin[MeV]\n"
            "# PDG: 13 = mu-   -13 = mu+\n"
        )
        for fpath in files:
            e_gev = energy_from_name(fpath)
            n_rows = 0
            with open(fpath) as fin:
                for raw in fin:
                    line = raw.strip()
                    if not line or line.startswith('#'):
                        continue
                    fout.write(raw)
                    n_rows += 1
            total += n_rows
            print(f"    {energy_label(e_gev):>8s}  {n_rows:>8,} muons")
    print(f"    Total: {total:,} muons")
    return total


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Merge per-energy UCMuon benchmark files into three merged sources')
    parser.add_argument('directory', nargs='?', default='output',
                        help='Directory containing per-energy source files (default: output)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory for merged files (default: same as directory)')
    parser.add_argument('--exclude', nargs='+', type=int, default=[1000, 10000],
                        metavar='E_GEV',
                        help='Energies (GeV) to exclude from the merge '
                             '(default: 1000 10000)')
    args = parser.parse_args()

    src_dir  = args.directory
    out_dir  = args.outdir if args.outdir else src_dir
    excluded = set(args.exclude)

    if not os.path.isdir(src_dir):
        print(f"ERROR: directory not found: {src_dir}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(out_dir, exist_ok=True)

    def _keep(path):
        return energy_from_name(path) not in excluded

    surface_files = [f for f in find_files(src_dir, r'muons_surface_\d+\.dat$')    if _keep(f)]
    phits_files   = [f for f in find_files(src_dir, r'muons_for_phits_\d+\.dat$')  if _keep(f)]
    geant4_files  = [f for f in find_files(src_dir, r'muons_geant4_\d+\.txt$')     if _keep(f)]

    if not surface_files:
        print(f"ERROR: no muons_surface_<E>.dat files found in {src_dir}", file=sys.stderr)
        sys.exit(1)

    energies = sorted(set(energy_from_name(f) for f in surface_files))
    print("=" * 60)
    print("  UCMuon benchmark merge")
    print("=" * 60)
    print(f"  Source directory : {src_dir}")
    print(f"  Output directory : {out_dir}")
    if excluded:
        print(f"  Excluded         : {', '.join(energy_label(e) for e in sorted(excluded))}")
    print(f"  Energies merged  : {', '.join(energy_label(e) for e in energies)}")
    print(f"  Surface files    : {len(surface_files)}")
    print(f"  PHITS files      : {len(phits_files)}")
    print(f"  Geant4 files     : {len(geant4_files)}")
    print("=" * 60)

    n_surf = merge_surface(surface_files,
                           os.path.join(out_dir, 'benchmark_surface.dat'))

    if phits_files:
        merge_phits(phits_files,
                    os.path.join(out_dir, 'benchmark_phits.dat'))
    else:
        print("\n  WARNING: no PHITS files found — skipping benchmark_phits.dat")

    if geant4_files:
        merge_geant4(geant4_files,
                     os.path.join(out_dir, 'benchmark_geant4.txt'))
    else:
        print("\n  WARNING: no Geant4 files found — skipping benchmark_geant4.txt")

    print()
    print("=" * 60)
    print("  Done.  Files written:")
    print(f"    {os.path.join(out_dir, 'benchmark_surface.dat')}")
    print(f"    {os.path.join(out_dir, 'benchmark_phits.dat')}")
    print(f"    {os.path.join(out_dir, 'benchmark_geant4.txt')}")
    print()
    print("  Use benchmark_surface.dat as input for UCMuon transport runs.")
    print("  Use benchmark_phits.dat   with:  s-type=17 / dump=-10 / 1 2 3 4 5 6 7 8 9 10")
    print("  Use benchmark_geant4.txt  with your Geant4 PrimaryGeneratorAction reader.")
    print("=" * 60)


if __name__ == '__main__':
    main()
