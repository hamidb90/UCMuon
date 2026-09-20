#!/bin/bash
# =============================================================================
# make_all_figs.sh — regenerate every figure of the CPC paper.
# Run from anywhere:  bash manuscript/scripts/make_all_figs.sh
#
# fig01  pipeline schematic            (pure matplotlib)
# fig02  flux spectra                  (analytic + PARMA helper)
# fig03  zenith distributions          (analytic)
# fig04  GUI screenshot                (needs manual capture, see script)
# fig05  surface/underground spectra   (runs gen + MUSIC, ~30 s)
# fig06  survival curves               (benchmark data, make_fig06.py)
# fig07  OpenMP scaling                (runs scaling benchmark, ~5 min
#                                       on first run; cached afterwards)
# fig08+09 source-spectrum validation  (reads benchmark/sources CSVs)
# fig10  six-code timing comparison    (reads benchmark timing.txt files)
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"

fail=0
for script in make_fig01_pipeline.py make_fig02_flux_spectra.py \
              make_fig03_zenith_dist.py make_fig05_energy_spectra.py \
              make_fig06.py make_fig07_scaling.py \
              make_fig08_source_validation.py make_fig10_timing.py; do
  echo "=== $script ==="
  python3 "$script" "$@" || { echo "[FAIL] $script"; fail=1; }
  echo
done

echo "=== make_fig04_gui_screenshot.py ==="
if [ -f ../figs/gui_screenshot.png ]; then
  python3 make_fig04_gui_screenshot.py || fail=1
else
  echo "[TODO] capture the GUI Results tab to manuscript/figs/gui_screenshot.png"
  echo "       then run: python3 manuscript/scripts/make_fig04_gui_screenshot.py"
fi

exit $fail
