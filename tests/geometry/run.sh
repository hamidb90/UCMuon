#!/usr/bin/env bash
# Build and run the detector-filter intersection test (see the header of
# test_ray_intersections.f90). Exit status is the test's.
set -eu
cd "$(dirname "$0")/../.."
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
gfortran -O2 -J"$tmp" src/generator/geom_module.f90 \
  tests/geometry/test_ray_intersections.f90 -o "$tmp/test_ray_intersections"
"$tmp/test_ray_intersections"
