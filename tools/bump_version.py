#!/usr/bin/env python3
"""Move the version, everywhere it is written, in one step.

The version is stated in five places that must agree, and in the paper when the
paper is in this tree. Changing four of five is the easiest mistake to make and
the least visible, so this does all of them, prints the diff it made, and leaves
a dated stub at the top of the CHANGELOG for the entry to be written into.

    python3 tools/bump_version.py 1.1.3
    python3 tools/bump_version.py 1.1.3 --also ../UCMuon-public
    python3 tools/bump_version.py 1.1.3 --dry-run

It does not commit, tag, push or touch Zenodo DOIs: the version DOI for a
release does not exist until the release is archived, and is added afterwards.
Run tools/check_consistency.py when it is done.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def edits_for(version: str) -> list[tuple[str, str, str]]:
    """(file, regex, replacement). The regex must match exactly once."""
    return [
        ("CITATION.cff", r'^version:\s*"[^"]+"', f'version: "{version}"'),
        ("CITATION.cff", r'^date-released:\s*"[^"]+"',
         f'date-released: "{dt.date.today().isoformat()}"'),
        (".zenodo.json", r'"version"\s*:\s*"[^"]+"', f'"version": "{version}"'),
        ("gui/ucmuon_gui.py", r'^__version__\s*=\s*"[^"]+"', f'__version__ = "{version}"'),
        ("README.md", r'^## Status & scope \(v[0-9][^)]*\)', f"## Status & scope (v{version})"),
    ]


def apply_to(root: Path, version: str, dry: bool) -> list[str]:
    changed: list[str] = []
    for rel, pattern, repl in edits_for(version):
        p = root / rel
        if not p.is_file():
            print(f"  -- {rel}: not in this tree, skipped")
            continue
        text = p.read_text()
        new, n = re.subn(pattern, repl, text, count=1, flags=re.M)
        if n != 1:
            print(f"  !! {rel}: pattern not found; fix it by hand", file=sys.stderr)
            continue
        if new == text:
            print(f"  == {rel}: already there")
            continue
        if not dry:
            p.write_text(new)
        changed.append(rel)
        old = re.search(pattern, text, re.M).group(0).strip()
        print(f"  -> {rel}: {old}  ==>  {repl.strip()}")

    # The paper names the version it describes; in the public tree there is none.
    paper = root / "manuscript" / "ucmuon_cpc_paper.tex"
    if paper.is_file():
        text = paper.read_text()
        others = sorted(set(re.findall(r"v(\d+\.\d+\.\d+)", text)) - {version})
        if others:
            new = text
            for old in others:
                new = new.replace(f"v{old}", f"v{version}")
            if not dry:
                paper.write_text(new)
            changed.append("manuscript/ucmuon_cpc_paper.tex")
            print(f"  -> manuscript/ucmuon_cpc_paper.tex: v{', v'.join(others)}  ==>  v{version}")
            print("     (rebuild the PDF: bash manuscript/make_submission.sh)")

    # A stub, not an entry: what changed is for a person to write.
    ch = root / "CHANGELOG.md"
    if ch.is_file():
        text = ch.read_text()
        if f"## [{version}]" in text:
            print("  == CHANGELOG.md: already has an entry for this version")
        else:
            stub = (f"## [{version}] — {dt.date.today().isoformat()}\n\n"
                    "### Added\n- \n\n### Fixed\n- \n\n")
            if not dry:
                ch.write_text(text.replace("# Changelog\n\n", "# Changelog\n\n" + stub, 1))
            changed.append("CHANGELOG.md")
            print(f"  -> CHANGELOG.md: stub entry added for {version}, fill it in")
    return changed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("version", help="the new version, e.g. 1.1.3")
    ap.add_argument("--also", metavar="PATH", action="append", default=[],
                    help="another tree to apply the same bump to (repeatable)")
    ap.add_argument("--dry-run", action="store_true", help="print, change nothing")
    args = ap.parse_args()

    if not SEMVER.match(args.version):
        print(f"error: '{args.version}' is not a MAJOR.MINOR.PATCH version", file=sys.stderr)
        return 2

    root = Path(__file__).resolve().parent.parent
    trees = [root] + [Path(p).resolve() for p in args.also]
    for tree in trees:
        print(f"\n{tree}")
        if not (tree / "CITATION.cff").is_file():
            print("  !! no CITATION.cff here; is this a UCMuon tree?", file=sys.stderr)
            continue
        apply_to(tree, args.version, args.dry_run)

    print("\nNext: write the CHANGELOG entry, then")
    print("  python3 tools/check_consistency.py")
    if args.dry_run:
        print("\n(dry run: nothing was written)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
