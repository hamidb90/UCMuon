# Changing things safely

Most of what breaks a release is not a bug in the physics. It is a version
bumped in four files out of five, a path quoted in one place and moved in
another, or a number that a document repeats and no longer matches the code.
Each is cheap to fix and expensive for someone else to find. The tools here
exist so that anything can be changed and the damage is caught in the same
sitting.

## Before you push anything

```bash
python3 tools/check_consistency.py          # about ten seconds
python3 tools/check_consistency.py --fast   # without the test run
```

It checks that the five files stating the version agree, that the release has a
Zenodo DOI recorded (a warning until it is archived, which is normal), that the
test run reproduces its reference, that nothing private is tracked in the
published tree, and that the two trees hold identical copies of every shared
file. Where this repository also holds the paper, it checks that the paper names
the same version, that every repository path the paper cites exists in the
published tree, that every figure still has a script that draws it, and that the
highlights are within the Elsevier limit.

Checks whose inputs are missing are reported as skipped, not failed, so the
published copy of the script is as useful as the private one. Warnings never
fail the run; failures exit non-zero.

## Changing the version

```bash
python3 tools/bump_version.py 1.1.3 --also ../UCMuon-public
python3 tools/bump_version.py 1.1.3 --dry-run     # look first
```

It moves the version in `CITATION.cff` (with today's release date),
`.zenodo.json`, `gui/ucmuon_gui.py` and the README status heading, opens a dated
stub in `CHANGELOG.md` to write the entry into, and, where the paper is present,
updates the version the paper says it describes. It does not commit, tag, push
or invent a Zenodo DOI: the version DOI does not exist until the release is
archived, and is added afterwards.

## Changing a figure of the paper

The figures themselves are not published; the scripts that draw them are, in
`manuscript/scripts/`, one per figure. Regenerate with
`python3 manuscript/scripts/make_figNN_*.py`, or `bash scripts/make_all_figs.sh`
for all of them.

If the change is only to the drawing, nothing else follows. If it changes what
the figure shows, or the numbers behind it, then the published scripts and
benchmark products no longer match the paper, and that wants a release.

## Making a release

1. `python3 tools/bump_version.py <version> --also ../UCMuon-public`, then write
   the CHANGELOG entry.
2. `python3 tools/check_consistency.py` until it passes.
3. Commit in both trees. The published tree is `../UCMuon-public`; this one has
   no remote and must never be pushed.
4. Push, tag and release. These steps are run by hand, deliberately:

   ```bash
   git -C ../UCMuon-public push origin main
   git -C ../UCMuon-public tag -a v<version> -m "UCMuon <version>"
   git -C ../UCMuon-public push origin v<version>
   git tag -a v<version> -m "UCMuon <version>"
   gh release create v<version> --title "UCMuon v<version>" --notes-file <notes>
   ```
5. When Zenodo has minted the version DOI, add it to `CITATION.cff` in both
   trees, after the previous entry, and push. The concept DOI
   10.5281/zenodo.20826984 stays the one to cite.
6. Download the published archive and check it: that it matches the tag, that
   `bash test_run/run_test.sh` passes from it, and that it carries nothing that
   was meant to stay private.

## What must never be published

The manuscript, its figures, the cover letter, the highlights file, the
declaration of interest, the submission package, the release-check notes, and
the MUSIC sources and data tables, which are not ours to redistribute. The
published tree's `.gitignore` enforces the manuscript part, and
`tools/check_consistency.py` fails if any of the rest appears in it.
