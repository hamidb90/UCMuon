#!/usr/bin/env python3
"""
ucmuon_terrain_example.py
─────────────────────────────────────────────────────────────────────────────
Standalone example: generate the correct surface muon file for the
UCMuon Terrain tab (Puy de Dôme, south flank detector).

Run this from the UCMuon project root (not from gui/):
    python ucmuon_terrain_example.py

It reads the existing muons_surface.dat (if already generated in Tab 1)
and tells you whether it is suitable.  If not, it prints the exact Tab 1
settings you need and writes a minimal test file for quick validation.

Author: UCLouvain Muography Group
"""

import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Puy de Dôme example parameters
# ─────────────────────────────────────────────────────────────────────────────
DET_LAT   = 45.760000   # south flank detector
DET_LON   =  2.955000
DET_ALT_M = 1094.0      # from DEM query

SUMMIT_LAT  = 45.7716
SUMMIT_ALT  = 1462.0    # DEM max
H_HORIZ     = (SUMMIT_LAT - DET_LAT) * 111320   # ~1291 m
H_DIFF      = SUMMIT_ALT - DET_ALT_M             # 368 m
ZE_CRIT_DEG = np.degrees(np.arctan2(H_HORIZ, H_DIFF))

print("═" * 65)
print("  UCMuon Terrain — Puy de Dôme example")
print("═" * 65)
print()
print(f"Detector position:")
print(f"  lat = {DET_LAT}°N   lon = {DET_LON}°E   alt = {DET_ALT_M} m a.s.l.")
print()
print(f"Summit geometry:")
print(f"  horizontal distance to summit  = {H_HORIZ:.0f} m")
print(f"  height difference              = {H_DIFF:.0f} m")
print(f"  critical blocking angle        = {ZE_CRIT_DEG:.1f}°")
print()
print(f"  → You MUST set Max zenith > {ZE_CRIT_DEG:.0f}° in Tab 5")
print(f"  → Recommended: ze_max = 85°  (last bin centre at 82.6°)")
print()

# ─────────────────────────────────────────────────────────────────────────────
# Check existing surface file
# ─────────────────────────────────────────────────────────────────────────────
surf_file = Path("muons_surface.dat")
if surf_file.exists():
    with open(surf_file) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    ncols = len(lines[0].split()) if lines else 0
    n_muons = len(lines)
    print(f"Found {surf_file}: {n_muons:,} muons, {ncols} columns")

    # Check angular coverage: do we have muons at ze > 74°?
    # In 14-col format: col 8 = theta [rad]
    if ncols in (13, 14):
        thetas = []
        for line in lines[:min(10000, len(lines))]:
            parts = line.split()
            thetas.append(float(parts[8]))
        thetas = np.degrees(np.array(thetas))
        n_near_horiz = np.sum(thetas > ZE_CRIT_DEG)
        frac = n_near_horiz / len(thetas) * 100
        print(f"  Muons with ze > {ZE_CRIT_DEG:.0f}° (sampled): {n_near_horiz} / {len(thetas)} = {frac:.2f}%")
        if n_muons >= 100_000 and thetas.max() > ZE_CRIT_DEG:
            print(f"  ✅ File is suitable for Tab 5 (has {n_muons:,} muons and covers ze>{ZE_CRIT_DEG:.0f}°)")
        elif thetas.max() <= ZE_CRIT_DEG:
            print(f"  ⚠️  Max zenith in file = {thetas.max():.1f}° — does not reach blocking angle!")
            print(f"     Regenerate with theta_max = 85° in Tab 1.")
        else:
            print(f"  ⚠️  Only {n_muons:,} muons — consider generating 500,000+ for better statistics.")
else:
    print(f"No {surf_file} found.")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Exact Tab 1 settings needed
# ─────────────────────────────────────────────────────────────────────────────
print()
print("═" * 65)
print("  EXACT TAB 1 SETTINGS for this example")
print("═" * 65)
print()
settings = [
    ("E min [GeV]",          "100",       "CSDA threshold for ~200 m basalt"),
    ("E max [GeV]",          "2500",      "CosmoALEPH upper limit"),
    ("Spectrum model",       "① CosmoALEPH", "recommended"),
    ("Source shape",         "💿 Circular disk", ""),
    ("Radius [m]",           "800",       "covers ze_max=85° from 1300m slant"),
    ("z plane [m]",          "0",         "surface"),
    ("N muons to generate",  "500,000",   "~750 muons per direction bin"),
    ("Angular mode",         "② cos²θ",  "realistic distribution"),
    ("θ_max [°]",            "85",        "MUST exceed critical angle 74.3°"),
    ("Detector filter",      "OFF",       "terrain tab handles all directions"),
    ("Output file",          "muons_surface.dat", ""),
]
for name, val, note in settings:
    comment = f"  ← {note}" if note else ""
    print(f"  {name:<30} {val:<20}{comment}")

print()
print("═" * 65)
print("  EXACT TAB 5 SETTINGS")
print("═" * 65)
print()
tab5 = [
    ("Surface muon file",    "muons_surface.dat", "generated above"),
    ("Output file",          "muons_terrain_ug.dat", ""),
    ("DEM file",             "puydedome_dem.tif", "SRTM GL1 30m, 0.6 MB"),
    ("Detector latitude",    f"{DET_LAT}",     "south flank"),
    ("Detector longitude",   f"{DET_LON}",     ""),
    ("Detector altitude",    f"{DET_ALT_M} m", "from DEM query"),
    ("Transport engine",     "UCMuon-MC", "recommended"),
    ("Rock density ρ",       "2.65 g/cm³",    "standard rock (or 2.85 for basalt)"),
    ("Max zenith",           "85°",           f"MUST be > {ZE_CRIT_DEG:.0f}°"),
    ("Azimuth bins",         "36",            "10° steps"),
    ("Zenith bins",          "18",            "~4.7° steps → last centre 82.6°"),
    ("Ray-trace step",       "50 m",          "2× SRTM 30m pixel"),
]
for name, val, note in tab5:
    comment = f"  ← {note}" if note else ""
    print(f"  {name:<30} {val:<25}{comment}")

print()
print("═" * 65)
print("  EXPECTED RESULTS")
print("═" * 65)
print()
print("  Overburden map (log scale):")
print("    North sector ze≈77–83°: 50,000–342,695 g/cm²  (red on Jet scale)")
print("    All other directions:    open sky (gray/masked)")
print("    38 of 648 bins blocked  (5.9% of directions)")
print()
print("  Overall survival rate:    ~99.9%  (physically correct)")
print("  Why: cos²θ puts <1.5% of muons in blocked directions ze>74°")
print()
print("  Survival rate map:")
print("    North blocked bins: 0.001–1%  survival  (almost all stopped)")
print("    Open-sky bins:      ~100%     survival  (no rock)")
print()
print("  This directional contrast IS the muographic signal.")
print("  Low survival in north = lava dome blocking muons.")
print()
print("═" * 65)
print("  GET THE DEM FILE")
print("═" * 65)
print()
print("  Method A — eio command (requires: pip install elevation):")
print("    eio clip -o puydedome_dem.tif --bounds 2.75 45.65 3.15 45.90")
print()
print("  Method B — OpenTopography website (no install):")
print("    https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1")
print("    Bounds: West=2.75  South=45.65  East=3.15  North=45.90")
print("    Format: GeoTIFF  →  Save as puydedome_dem.tif")
print()
print("  Method C — Tab 5 auto-download (in the GUI):")
print("    Section 2 → Auto-download tab")
print("    S=45.65  N=45.90  W=2.75  E=3.15")
print("    Product: SRTM GL1 (30m)")
print()

# ─────────────────────────────────────────────────────────────────────────────
# Quick validation: reproduce what the DEM ray-tracer would find
# ─────────────────────────────────────────────────────────────────────────────
print("═" * 65)
print("  QUICK GEOMETRY VALIDATION  (no DEM file needed)")
print("═" * 65)
print()
print("  Simulating the critical north-sector ray analytically:")
print()
for ze_deg in [72.9, 74.3, 77.9, 82.6]:
    ze_rad = np.radians(ze_deg)
    cos_ze = np.cos(ze_rad)
    sin_ze = np.sin(ze_rad)
    # Ray altitude at horizontal distance H_HORIZ
    # alt(dist) = DET_ALT_M + cos_ze * dist
    # dist when horizontal projection = H_HORIZ: dist = H_HORIZ / sin_ze
    dist = H_HORIZ / sin_ze
    ray_alt_at_summit = DET_ALT_M + cos_ze * dist
    terrain_alt = SUMMIT_ALT
    blocked = ray_alt_at_summit < terrain_alt
    if blocked:
        slant_cm  = dist * 100.0
        overburden = slant_cm * 2.65
        print(f"  ze={ze_deg:5.1f}°: ray alt at summit = {ray_alt_at_summit:6.0f} m  "
              f"< terrain {terrain_alt:.0f} m  → BLOCKED  "
              f"overburden={overburden/1000:.0f}k g/cm²  ({dist:.0f} m slant)")
    else:
        print(f"  ze={ze_deg:5.1f}°: ray alt at summit = {ray_alt_at_summit:6.0f} m  "
              f"> terrain {terrain_alt:.0f} m  → open sky")

print()
print("  Last bin centre with ze_max=75° (OLD default): 72.9° → open sky ✗")
print("  Last bin centre with ze_max=85° (FIXED):       82.6° → BLOCKED ✓")
print()
