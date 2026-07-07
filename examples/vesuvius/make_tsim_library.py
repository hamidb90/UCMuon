#!/usr/bin/env python3
# make_tsim_library.py — regenerate the Vesuvius/MURAVES T_sim transmission library.
#
# Run from the repo root:   python3 examples/vesuvius/make_tsim_library.py
#
# ── DEM dependency ────────────────────────────────────────────────────────────
# This script needs a Vesuvius DEM at  misc/dem_site.tif  (a GeoTIFF covering the
# summit region around the MURAVES detector at 40.8271 N, 14.4006 E). GeoTIFFs are
# gitignored (see .gitignore: *.tif), so the file is NOT shipped in the repo — the
# tsim_library/ output it produces IS shipped, so you only need the DEM to
# *regenerate* the library, not to use it.
#
# To obtain the DEM (see README §"Get a DEM file" for full options):
#   • GUI:  Tab 2 → UCMuon Terrain → Section 1 → "Auto-download" (no account), OR
#   • CLI:  pip install elevation && \
#           eio clip -o misc/dem_site.tif --bounds 14.30 40.75 14.50 40.90
#           (bounds = LON_W LAT_S LON_E LAT_N around Vesuvius), OR
#   • OpenTopography SRTM GL1 30 m → export GeoTIFF → save as misc/dem_site.tif
# ──────────────────────────────────────────────────────────────────────────────
import sys, time, numpy as np
sys.path.insert(0, 'gui')
import ucmuon_terrain_driver as td

DET = (40.8271, 14.4006, 608.0)          # MURAVES SW-flank detector
RHO_REF = 2.65
RHOS = [1.5, 2.0, 2.5, 2.65, 3.0]
OUT = 'examples/vesuvius/tsim_library'

DEM_PATH = 'misc/dem_site.tif'           # gitignored; see header note to obtain it
elev, transform = td.load_dem(DEM_PATH)
t0 = time.time()
az_c, ze_c, ob_ref, sky = td.compute_overburden_map(
    elev, transform, *DET, RHO_REF, 360, 85, 85.0, 25.0,
    progress_cb=lambda d, t: print(f"overburden {d}/{t}", flush=True) if d % 3400 == 0 else None)
print(f"OVERBURDEN DONE in {(time.time()-t0)/60:.1f} min; blocked {np.sum(~sky)}/{sky.size}", flush=True)
td.write_overburden_map(az_c, ze_c, ob_ref, sky, f'{OUT}/terrain_overburden_ref2.65.dat', *DET, RHO_REF)

prev_T = None
for rho in RHOS:
    t1 = time.time()
    ob = ob_ref * (rho / RHO_REF)        # uniform-density scaling, no re-trace
    fmap, osky = td.compute_flux_map(az_c, ze_c, ob, sky, rho,
                                     spectrum_mode=1, mode=1, n_E=40,
                                     script_dir='gui')
    T = td.compute_transmission_map(fmap, osky)
    # validation
    open_ok = np.allclose(T[sky], 1.0)
    blocked = ~sky & (ob > 1.0)
    mono_ok = True
    if prev_T is not None:
        mono_ok = np.all(T[blocked] <= prev_T[blocked] + 1e-12)
    fname = f'{OUT}/terrain_transmission_{rho:.2f}.dat'
    td.write_transmission_map(az_c, ze_c, T, fname, *DET, rho)
    print(f"RHO {rho:.2f}: open-sky T=1 {'OK' if open_ok else 'FAIL'}; "
          f"monotone-vs-prev-rho {'OK' if mono_ok else 'FAIL'}; "
          f"min T {np.nanmin(T[blocked]):.3e}; blocked-bin mean T {np.nanmean(T[blocked]):.4f}; "
          f"{time.time()-t1:.0f}s -> {fname}", flush=True)
    prev_T = T
print("LIBRARY COMPLETE", flush=True)
