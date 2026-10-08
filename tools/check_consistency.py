#!/usr/bin/env python3
"""Check that the things which have to agree, agree.

A release, a paper and a repository drift apart in small ways: a version bumped
in one file and not another, a path the paper cites that nobody published, a
figure with no script left to draw it. Each is cheap to fix and expensive for a
referee to find. This script looks for all of them in one pass, so anything can
be changed and the damage is caught before it ships.

    python3 tools/check_consistency.py              # everything, ~10 s
    python3 tools/check_consistency.py --fast       # skip the test run
    python3 tools/check_consistency.py --public ../UCMuon-public

It runs in either tree. Checks whose inputs are absent (no manuscript in the
public tree, no sibling tree to compare against) are reported as skipped, not
as failures, so the public copy of this file is as useful as the private one.

Exit status is 0 when nothing failed, 1 otherwise. Warnings never fail the run:
they are the things that are normal in one phase and wrong in another, such as
a version with no Zenodo DOI yet, which is expected right up until the release.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

OK, WARN, FAIL, SKIP = "ok", "warn", "FAIL", "skip"
results: list[tuple[str, str, str]] = []


def record(status: str, check: str, detail: str = "") -> None:
    results.append((status, check, detail))


def read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


# ---------------------------------------------------------------- versions ---
def declared_versions(root: Path) -> dict[str, str | None]:
    """The version as each file states it. None means the file is unreadable or
    the pattern is gone, which is itself worth reporting."""
    def grab(rel: str, pattern: str) -> str | None:
        text = read(root / rel)
        if text is None:
            return None
        m = re.search(pattern, text, re.M)
        return m.group(1) if m else None

    return {
        "CITATION.cff": grab("CITATION.cff", r'^version:\s*"([^"]+)"'),
        ".zenodo.json": grab(".zenodo.json", r'"version"\s*:\s*"([^"]+)"'),
        "gui/ucmuon_gui.py": grab("gui/ucmuon_gui.py", r'^__version__\s*=\s*"([^"]+)"'),
        "README.md": grab("README.md", r'^## Status & scope \(v([0-9][^)]*)\)'),
        "CHANGELOG.md": grab("CHANGELOG.md", r'^## \[([0-9][^\]]*)\]'),
    }


def check_versions(root: Path) -> str | None:
    declared = declared_versions(root)
    missing = [f for f, v in declared.items() if v is None]
    if missing:
        record(FAIL, "version strings readable", "no version found in: " + ", ".join(missing))
        return None
    distinct = sorted(set(declared.values()))
    if len(distinct) != 1:
        detail = "; ".join(f"{f}={v}" for f, v in declared.items())
        record(FAIL, "version strings agree", detail)
        return None
    version = distinct[0]
    record(OK, "version strings agree", f"all five say {version}")
    return version


def check_doi(root: Path, version: str | None) -> None:
    text = read(root / "CITATION.cff")
    if text is None or version is None:
        record(SKIP, "Zenodo DOI for this version")
        return
    if re.search(r'description:\s*"Version DOI for v%s"' % re.escape(version), text):
        record(OK, "Zenodo DOI for this version", f"v{version} identifier present")
    else:
        record(WARN, "Zenodo DOI for this version",
               f"no identifier for v{version} yet; expected until the release is archived")


# ------------------------------------------------------------------ paper ---
PAPER = "manuscript/ucmuon_cpc_paper.tex"


def check_paper_version(root: Path, version: str | None) -> None:
    text = read(root / PAPER)
    if text is None:
        record(SKIP, "paper version matches the code")
        return
    quoted = sorted(set(re.findall(r"v(\d+\.\d+\.\d+)", text)))
    if not quoted:
        record(WARN, "paper version matches the code", "the paper names no version")
    elif version is not None and quoted == [version]:
        record(OK, "paper version matches the code", f"v{version} in {text.count('v' + version)} places")
    else:
        record(FAIL, "paper version matches the code",
               f"paper says {', '.join(quoted)}, the code says {version}")


def cited_paths(text: str) -> list[str]:
    raw = re.findall(r"\\texttt\{([^}]*)\}", text)
    out = set()
    for item in raw:
        item = item.replace("\\_", "_").replace("\\&", "&").strip()
        if "/" not in item or item.startswith("/") or item.startswith("foss/"):
            continue
        # A repository path, not a shell command, a URL or a module name.
        if " " in item or "://" in item or any(c in item for c in "\\${}"):
            continue
        out.add(item)
    return sorted(out)


def check_cited_paths(root: Path, public: Path | None) -> None:
    text = read(root / PAPER)
    if text is None:
        record(SKIP, "paths the paper cites exist")
        return
    target = public if public is not None else root
    where = "the public tree" if public is not None else "this tree"
    # bin/ is created by make, so it is absent from a fresh checkout by design.
    missing = [p for p in cited_paths(text)
               if p != "bin/" and not (target / p).exists()]
    if missing:
        record(FAIL, "paths the paper cites exist",
               f"absent from {where}: " + ", ".join(missing))
    else:
        record(OK, "paths the paper cites exist", f"all resolve in {where}")


def check_figures_have_scripts(root: Path) -> None:
    text = read(root / PAPER)
    scripts = root / "manuscript" / "scripts"
    if text is None or not scripts.is_dir():
        record(SKIP, "every figure has a script")
        return
    figs = sorted({m.rsplit("/", 1)[-1]
                   for m in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]*)\}", text)})
    body = "\n".join(filter(None, (read(p) for p in scripts.glob("*"))))
    orphans = [f for f in figs if f not in body]
    if orphans:
        record(WARN, "every figure has a script",
               "no script mentions: " + ", ".join(orphans))
    else:
        record(OK, "every figure has a script", f"{len(figs)} figures accounted for")


def check_highlights(root: Path) -> None:
    hl = read(root / "manuscript" / "highlights.txt")
    if hl is None:
        record(SKIP, "highlights within the Elsevier limit")
        return
    bullets = [l.strip() for l in hl.splitlines() if l.strip().startswith("•")]
    if not bullets:
        record(FAIL, "highlights within the Elsevier limit", "no bullets found")
        return
    over = [b for b in bullets if len(b.lstrip("• ")) > 85]
    if over:
        record(FAIL, "highlights within the Elsevier limit",
               "; ".join(f"{len(b.lstrip('• '))} chars: {b[:50]}…" for b in over))
    elif not 3 <= len(bullets) <= 5:
        record(FAIL, "highlights within the Elsevier limit",
               f"{len(bullets)} bullets, Elsevier wants 3 to 5")
    else:
        longest = max(len(b.lstrip("• ")) for b in bullets)
        record(OK, "highlights within the Elsevier limit",
               f"{len(bullets)} bullets, longest {longest} characters")


def check_cover_letter(root: Path, version: str | None) -> None:
    text = read(root / "manuscript" / "cover_letter.md")
    if text is None:
        record(SKIP, "cover letter ready to send")
        return
    problems = []
    if "to be filled" in text.lower() or "TODO" in text:
        problems.append("still has a placeholder")
    if version and f"v{version}" not in text:
        problems.append(f"does not name v{version}")
    if problems:
        record(FAIL, "cover letter ready to send", "; ".join(problems))
    else:
        record(OK, "cover letter ready to send")


# ----------------------------------------------------------------- public ---
PRIVATE_MARKERS = (
    ("the manuscript itself", re.compile(r"^manuscript/(?!scripts/)")),
    ("submission material", re.compile(r"submission_cpc|cover_letter|highlights\.txt|declaration_of_interest|credit_author_statement")),
    ("release checks", re.compile(r"release_check|release_notes")),
    ("MUSIC sources or tables", re.compile(r"src/transport/music/.*\.(f|for)$|data/music-.*\.dat$")),
    ("superseded benchmark material", re.compile(r"benchmark/figures/four_code/|benchmark/reports/(notes|benchmark_plan|COMPARISON_ALL_CODES)")),
    # The pre-publication working roadmap (June 2026); published up to 1.3.1.
    ("the internal roadmap", re.compile(r"^ROADMAP\.md$")),
)


def tracked_files(tree: Path) -> list[str] | None:
    try:
        out = subprocess.run(["git", "-C", str(tree), "ls-files"],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.splitlines()


def check_public_clean(public: Path | None) -> None:
    if public is None:
        record(SKIP, "nothing private is published")
        return
    files = tracked_files(public)
    if files is None:
        record(SKIP, "nothing private is published", "not a git tree")
        return
    hits = []
    for label, pattern in PRIVATE_MARKERS:
        found = [f for f in files if pattern.search(f)]
        if found:
            hits.append(f"{label}: {found[0]}" + (f" (+{len(found)-1} more)" if len(found) > 1 else ""))
    if hits:
        record(FAIL, "nothing private is published", "; ".join(hits))
    else:
        record(OK, "nothing private is published", f"{len(files)} tracked files checked")


def check_tree_parity(root: Path, public: Path | None) -> None:
    """Shared files must be identical. README.md and .gitignore differ by design."""
    if public is None:
        record(SKIP, "shared files identical in both trees")
        return
    files = tracked_files(public)
    if files is None:
        record(SKIP, "shared files identical in both trees", "not a git tree")
        return
    allowed = {"README.md", ".gitignore"}
    differing = []
    for rel in files:
        a, b = root / rel, public / rel
        if not a.is_file() or not b.is_file():
            continue
        if a.read_bytes() != b.read_bytes() and rel not in allowed:
            differing.append(rel)
    if differing:
        record(FAIL, "shared files identical in both trees",
               ", ".join(differing[:5]) + (f" (+{len(differing)-5} more)" if len(differing) > 5 else ""))
    else:
        record(OK, "shared files identical in both trees",
               "apart from README.md and .gitignore, which differ on purpose")


# ------------------------------------------------------------------ tests ---
def check_test_run(root: Path) -> None:
    script = root / "test_run" / "run_test.sh"
    if not script.is_file():
        record(SKIP, "test run passes")
        return
    try:
        proc = subprocess.run(["bash", str(script)], cwd=root,
                              capture_output=True, text=True, timeout=900)
    except (OSError, subprocess.TimeoutExpired) as exc:
        record(FAIL, "test run passes", str(exc))
        return
    if proc.returncode == 0 and "RESULT: PASS" in proc.stdout:
        tail = [l for l in proc.stdout.splitlines() if "[PASS]" in l]
        record(OK, "test run passes", f"{len(tail)} graded checks")
    else:
        record(FAIL, "test run passes", (proc.stdout or proc.stderr).strip().splitlines()[-1:] or ["no output"])


def check_ucmugen_numbers(root: Path) -> None:
    checker = root / "ucmugen" / "validation" / "check_numbers.py"
    tour = root / "ucmugen" / "examples" / "features" / "feature_tour.cc"
    if not checker.is_file() or not tour.is_file():
        record(SKIP, "UCMuGen documented numbers")
        return
    record(SKIP, "UCMuGen documented numbers",
           "run ucmugen/validation/run_tests.sh; it needs the tour compiled")


# ------------------------------------------------------------------- main ---
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--public", metavar="PATH", default="../UCMuon-public",
                    help="the published tree to compare against (default: ../UCMuon-public)")
    ap.add_argument("--fast", action="store_true", help="skip the test run")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    public = Path(args.public)
    if not public.is_absolute():
        public = (root / public).resolve()
    public_tree = public if (public.is_dir() and public != root) else None

    version = check_versions(root)
    check_doi(root, version)
    check_paper_version(root, version)
    check_cited_paths(root, public_tree)
    check_figures_have_scripts(root)
    check_highlights(root)
    check_cover_letter(root, version)
    check_public_clean(public_tree)
    check_tree_parity(root, public_tree)
    check_ucmugen_numbers(root)
    if not args.fast:
        check_test_run(root)

    width = max(len(c) for _, c, _ in results)
    print()
    for status, check, detail in results:
        mark = {OK: "  ok  ", WARN: " warn ", FAIL: " FAIL ", SKIP: " skip "}[status]
        line = f"[{mark}] {check.ljust(width)}"
        if detail:
            line += f"  {detail if isinstance(detail, str) else ' '.join(detail)}"
        print(line)

    failed = [c for s, c, _ in results if s == FAIL]
    warned = [c for s, c, _ in results if s == WARN]
    print()
    if failed:
        print(f"FAILED: {len(failed)} of {len(results)} checks ({', '.join(failed)})")
        return 1
    print(f"All {len(results)} checks passed" + (f", {len(warned)} warning(s)" if warned else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
