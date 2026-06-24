#!/usr/bin/env python3
"""
merge_sources.py  —  Merge UCMuon monoenergetic source files into one
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Usage
-----
  python3 merge_sources.py source_10GeV.dat source_20GeV.dat ... -o benchmark_source.dat

  Or use the --energies shorthand to auto-name input files:
  python3 merge_sources.py --energies 10 20 50 100 300 1000 10000 -o benchmark_source.dat

Output
------
  Standard UCMuon 14-column generator format with:
    - EventIDs renumbered 1..N (sequential across all batches)
    - Header comments listing each energy batch and its EventID range
    - A companion index file  <output>.index  mapping EventID → initial_E_GeV
      (used during analysis to split results by initial energy)
"""

import argparse
import os
import sys


HEADER = (
    "# UCMuon benchmark source — merged monoenergetic batches\n"
    "# UCLouvain Muography Group\n"
    "# Format: EventID  x_cm  y_cm  z_cm  p_GeV  px_GeV  py_GeV  pz_GeV"
    "  theta_rad  phi_rad  E_GeV  charge  hit_flag  det_mask\n"
)


def parse_energy_label(filename):
    """Extract a human-readable energy label from the filename."""
    base = os.path.splitext(os.path.basename(filename))[0]
    # e.g.  source_10GeV  →  10 GeV
    for suffix in ("tev", "TeV", "gev", "GeV"):
        if suffix.lower() in base.lower():
            idx = base.lower().index(suffix.lower())
            num = base[:idx].split("_")[-1]
            unit = "TeV" if "t" in suffix.lower() else "GeV"
            return f"{num} {unit}"
    return base


def numeric_energy_GeV(filename):
    """Return energy in GeV inferred from filename for the index file."""
    base = os.path.splitext(os.path.basename(filename))[0].lower()
    for suffix, scale in (("tev", 1000.0), ("gev", 1.0)):
        if suffix in base:
            idx = base.index(suffix)
            try:
                return float(base[:idx].split("_")[-1]) * scale
            except ValueError:
                pass
    return 0.0


def merge(input_files, output_path):
    index_path = output_path + ".index"

    event_id   = 0
    batch_info = []   # (label, first_id, last_id)

    with open(output_path, "w") as fout, open(index_path, "w") as fidx:
        fidx.write("# EventID  initial_E_GeV\n")

        # Write header placeholder — will prepend batch table after we know ranges
        fout.write(HEADER)

        for fpath in input_files:
            label   = parse_energy_label(fpath)
            e_GeV   = numeric_energy_GeV(fpath)
            first   = event_id + 1
            n_written = 0

            if not os.path.isfile(fpath):
                print(f"  WARNING: file not found, skipping: {fpath}", file=sys.stderr)
                continue

            with open(fpath) as fin:
                for raw in fin:
                    line = raw.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 13:
                        continue

                    event_id += 1
                    # Replace original EventID (col 0) with new sequential ID
                    parts[0] = str(event_id)
                    fout.write("  ".join(parts) + "\n")
                    fidx.write(f"{event_id}  {e_GeV:.1f}\n")
                    n_written += 1

            last = event_id
            batch_info.append((label, first, last, n_written))
            print(f"  {label:>10s}  {n_written:>8,} muons  EventID {first}–{last}")

    # Prepend batch table into a fresh copy of the file
    with open(output_path, "r") as f:
        original = f.read()

    table_lines = ["#\n", "# Batch summary:\n",
                   "#   Energy       N_muons    EventID range\n"]
    for label, first, last, n in batch_info:
        table_lines.append(f"#   {label:<12s}  {n:>8,}    {first}–{last}\n")
    table_lines.append("#\n")

    # Insert table after the first 3 header lines
    lines = original.split("\n", 3)
    with open(output_path, "w") as f:
        f.write("\n".join(lines[:3]) + "\n")
        f.writelines(table_lines)
        f.write(lines[3])

    total = sum(n for _, _, _, n in batch_info)
    print(f"\n  Total: {total:,} muons → {output_path}")
    print(f"  Index: {index_path}  (maps EventID → initial energy for analysis)")


def main():
    parser = argparse.ArgumentParser(
        description="Merge UCMuon monoenergetic source files into one benchmark source")
    parser.add_argument("inputs", nargs="*",
                        help="Input .dat files (in order, lowest to highest energy)")
    parser.add_argument("--energies", nargs="+", type=float, metavar="E_GeV",
                        help="Energy list in GeV; auto-names files as source_<E>GeV.dat "
                             "(use e.g. 1000 for 1 TeV)")
    parser.add_argument("-o", "--output", default="benchmark_source.dat",
                        help="Output merged file (default: benchmark_source.dat)")
    args = parser.parse_args()

    if args.energies:
        files = []
        for e in args.energies:
            if e >= 1000:
                name = f"source_{int(e/1000)}TeV.dat"
            else:
                name = f"source_{int(e)}GeV.dat"
            files.append(name)
    elif args.inputs:
        files = args.inputs
    else:
        parser.print_help()
        sys.exit(1)

    print(f"Merging {len(files)} source files → {args.output}\n")
    merge(files, args.output)


if __name__ == "__main__":
    main()
