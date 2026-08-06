#!/usr/bin/env python3
"""Check the numbers the documentation quotes against the code that produces them.

The problem this solves is drift. A number lives in three places: the program
that computes it, the README that quotes it, and the paper that quotes it. Only
the first one executes, so the other two rot silently as the code changes
underneath them. That is not hypothetical here: the Geant4 example's README
advertised a muon count from before the example went multithreaded, and the
feature tour reported a ratio of acceptances under a column headed "speed-up",
overstating the real figure by a factor of three. Neither was caught by a test,
because no test looked at the numbers.

So: `feature_tour --numbers` emits every quoted value as `key value`, this
script diffs that against reference/tour_numbers.txt, and a value that has moved
fails the build *and names the documents that quote it*, so the fix is a matter
of following the list rather than remembering where a number was copied to.

Usage:
    python3 check_numbers.py --exe /path/to/feature_tour
    python3 check_numbers.py --exe /path/to/feature_tour --update

`--update` rewrites the reference from the current run. It is the right move
when a value changes on purpose; the diff it produces is then the checklist of
documents to update, and reviewing that diff is the point of the exercise.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "tour_numbers.txt"

# Tolerances by key prefix, relative. These are set by what the platform can do
# to a value, not by what looks tidy:
#
#   areas       pure arithmetic, identical everywhere.
#   rates       a mean over a fixed, deterministically sampled set of points, so
#               only libm ULPs and the occasional flipped ray-geometry test move
#               it. Observed spread is ~1e-12; 1e-4 is far looser than that and
#               still catches any change visible at the four significant figures
#               the documents quote.
#   acceptance  a count, and a single flipped accept/reject shifts the whole
#               subsequent random stream, so cross-platform it behaves like an
#               independent sample: ~0.7% relative at these statistics. 5% is
#               the non-flaky bound, and the failure this guards against (a
#               mislabelled quantity) is a factor of two or more.
TOLERANCES = (
    ("surface_area.", 1e-9),
    ("acceptance.", 5e-2),
    ("", 1e-4),
)


def tolerance_for(key: str) -> float:
    for prefix, tol in TOLERANCES:
        if key.startswith(prefix):
            return tol
    return 1e-4


def run(exe: Path) -> dict[str, float]:
    out = subprocess.run([str(exe), "--numbers"], capture_output=True, text=True,
                         check=True).stdout
    values: dict[str, float] = {}
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        key, _, value = line.partition(" ")
        values[key] = float(value)
    if not values:
        sys.exit("check_numbers: --numbers produced nothing")
    return values


def load_reference() -> dict[str, tuple[float, str]]:
    ref: dict[str, tuple[float, str]] = {}
    for line in REFERENCE.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("|")]
        key, value = parts[0], float(parts[1])
        quoted_in = parts[2] if len(parts) > 2 else ""
        ref[key] = (value, quoted_in)
    return ref


def write_reference(values: dict[str, float],
                    previous: dict[str, tuple[float, str]]) -> None:
    lines = [
        "# Numbers the documentation quotes, as produced by",
        "#     feature_tour --numbers",
        "#",
        "# Format:  key | value | documents that quote it",
        "#",
        "# Regenerate with check_numbers.py --update. Do not hand-edit a value:",
        "# the point of the file is that every value in it came out of a run.",
        "# The third column is the reason the file exists. When a value moves,",
        "# it is the list of places that have to move with it.",
        "",
    ]
    width = max(len(k) for k in values)
    for key, value in values.items():
        # Carry the "quoted in" annotation across an update; it is curated by
        # hand and has nothing to do with the measured value.
        quoted_in = previous.get(key, (0.0, ""))[1]
        row = f"{key:<{width}} | {value:<16.6f} |"
        lines.append(f"{row} {quoted_in}" if quoted_in else row)
    REFERENCE.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--update", action="store_true")
    args = ap.parse_args()

    values = run(args.exe)

    if args.update:
        previous = load_reference() if REFERENCE.exists() else {}
        write_reference(values, previous)
        print(f"check_numbers: wrote {len(values)} values to "
              f"{REFERENCE.relative_to(HERE.parent.parent)}")
        return 0

    if not REFERENCE.exists():
        sys.exit(f"check_numbers: no reference at {REFERENCE}; "
                 "create it with --update")

    reference = load_reference()
    failures: list[str] = []

    for key, (expected, quoted_in) in reference.items():
        if key not in values:
            failures.append(f"  {key}: gone from the program output")
            continue
        got = values[key]
        tol = tolerance_for(key)
        scale = max(abs(expected), 1e-12)
        rel = abs(got - expected) / scale
        if rel > tol:
            msg = (f"  {key}\n"
                   f"      reference {expected:.6g}, now {got:.6g} "
                   f"(rel {rel:.2g}, tolerance {tol:g})")
            if quoted_in:
                msg += f"\n      quoted in: {quoted_in}"
            failures.append(msg)

    for key in values:
        if key not in reference:
            failures.append(f"  {key}: new, not in the reference "
                            "(run --update and review the diff)")

    if failures:
        print(f"\n  {len(failures)} number(s) no longer match the reference:\n")
        print("\n".join(failures))
        print("\n  Either the change is a regression, or it is intended and the"
              "\n  documents listed above need updating. Once they are, refresh"
              "\n  the reference with:"
              "\n      python3 ucmugen/validation/check_numbers.py "
              "--exe <feature_tour> --update\n")
        return 1

    print(f"  ok: {len(reference)} documented numbers match the code")
    return 0


if __name__ == "__main__":
    sys.exit(main())
