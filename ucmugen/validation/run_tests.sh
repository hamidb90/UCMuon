#!/usr/bin/env bash
# Build and run the C++ validation programs.
#
# The Geant4 test is skipped, not failed, when geant4-config is not on the
# path: the other two must stay runnable on a machine with no Geant4, which is
# the whole point of the header being dependency-free.
#
# The Python suite (compare_legacy.py) is separate because it needs the Fortran
# generator built; see README.md.

set -u
cd "$(dirname "$0")"

CXX=${CXX:-c++}
FLAGS="-std=c++17 -O2 -Wall -Wextra -Wpedantic -I../include"

# The Geant4 test is built with the system default compiler, not $CXX. Setting
# CXX is how the portability matrix sweeps compilers, but `geant4-config
# --cflags` hands back the flags Geant4 itself was configured with (on a clang
# build that includes -Qunused-arguments, which gcc rejects outright), and the
# prebuilt libraries have to be ABI-compatible with whatever links them. Export
# G4CXX to override.
G4CXX=${G4CXX:-c++}
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

status=0

# Header hygiene, before the physics: a header-only library breaks in two ways
# that no single-file test can see. A function missing `inline` links fine until
# a second translation unit includes the header, and a broken include guard only
# shows up on the second include in one file. Both are cheap to rule out.
echo "=== header hygiene (double include, two-TU link)"
# Written into $OUT rather than via `mktemp -t`: BSD mktemp treats -t's
# argument as a prefix, GNU treats it as a template and fails unless it ends in
# XXXXXX. On GNU the substitution came back empty, both names collapsed to the
# literal ".cc", and the one file was compiled twice, which fails as a
# duplicate main() rather than as anything to do with the header.
hdr_a="$OUT/ucmugen_tu_a.cc"
hdr_b="$OUT/ucmugen_tu_b.cc"
printf '#include "UCMuGen.h"\n#include "UCMuGen.h"\nint tu_a(){ return int(ucmugen::flux::intensity(ucmugen::Spectrum::Guan,10.0,1.0)>0); }\n' > "$hdr_a"
printf '#include "UCMuGen.h"\nint tu_a();\nint main(){ ucmugen::Generator g; (void)g; return tu_a()?0:1; }\n' > "$hdr_b"
if $CXX $FLAGS "$hdr_a" "$hdr_b" -o "$OUT/hdr_hygiene" && "$OUT/hdr_hygiene"; then
  echo "  ok"
else
  echo "  FAILED"; status=1
fi
rm -f "$hdr_a" "$hdr_b"
echo

for t in test_projection test_detector; do
  echo "=== $t"
  if ! $CXX $FLAGS "$t.cc" -o "$OUT/$t"; then
    echo "  BUILD FAILED"; status=1; continue
  fi
  "$OUT/$t" || status=1
  echo
done

# PARMA is optional: the generated header is JAEA's, non-commercial-only, and a
# user who has removed it should still get a passing suite for the MIT core.
echo "=== test_parma"
if [ -f ../include/UCMuGen_PARMA.h ]; then
  if $CXX $FLAGS test_parma.cc -o "$OUT/test_parma"; then
    "$OUT/test_parma" || status=1
  else
    echo "  BUILD FAILED"; status=1
  fi
else
  echo "  skipped: UCMuGen_PARMA.h not present"
fi
echo

# The numbers the documentation quotes, against the code that produces them.
# Kept here rather than only in CI because the failure it catches is one a
# developer creates locally: change a default, and every README and paper line
# quoting a rate is silently wrong until someone re-measures. Needs python3
# only for the comparison; the values come from the example itself.
echo "=== documented numbers"
if $CXX $FLAGS ../examples/features/feature_tour.cc -o "$OUT/feature_tour"; then
  if command -v python3 >/dev/null 2>&1; then
    python3 check_numbers.py --exe "$OUT/feature_tour" || status=1
  else
    echo "  skipped: python3 not found"
  fi
else
  echo "  BUILD FAILED"; status=1
fi
echo

echo "=== test_geant4"
if command -v geant4-config >/dev/null 2>&1; then
  # geant4-config --cflags brings its own -std and warning set.
  if $G4CXX -std=c++17 -O2 $(geant4-config --cflags) -I../include \
          test_geant4.cc $(geant4-config --libs) -o "$OUT/test_geant4"; then
    DYLD_LIBRARY_PATH="$(geant4-config --prefix)/lib:${DYLD_LIBRARY_PATH:-}" \
    LD_LIBRARY_PATH="$(geant4-config --prefix)/lib:${LD_LIBRARY_PATH:-}" \
      "$OUT/test_geant4" || status=1
  else
    echo "  BUILD FAILED"; status=1
  fi
else
  echo "  skipped: geant4-config not found"
fi

echo
[ $status -eq 0 ] && echo "ALL SUITES PASSED" || echo "SOME SUITES FAILED"
exit $status
