#!/usr/bin/env python3
"""make_fig_ucmugen_directed.py — blind vs detector-directed sampling.

Draws the rays the generator actually produced, not a re-enactment of them:
ucmugen/validation/dump_rays.cc writes one row per generated muon and this
script plots those rows. The two panels use the same generator, the same
spectrum, the same surface and the same seed; the only difference is whether a
detector is attached.

Two panels rather than one plot with two colours: the quantity being compared is
where a whole population of rays goes, and overlaying the populations would hide
exactly the thing worth seeing.

Usage:
    c++ -std=c++17 -O2 -I ucmugen/include ucmugen/validation/dump_rays.cc \
        -o /tmp/dump_rays && /tmp/dump_rays > /tmp/rays.csv
    python3 manuscript/scripts/make_fig_ucmugen_directed.py /tmp/rays.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_PDF = SCRIPT_DIR.parent / "figs" / "fig_ucmugen_directed.pdf"
OUT_PNG = SCRIPT_DIR.parent / "figs" / "fig_ucmugen_directed.png"

# Geometry, in metres, mirroring dump_rays.cc.
DET_HX, DET_HY, DET_HZ = 0.5, 0.5, 0.1
SKY_HALF, SKY_Z = 5.0, 5.0
Z_FLOOR = -1.0

# Identity encoding. The hit colour is the teal already used for UCMuGen in the
# other figures; misses are neutral ink, not a second hue, because "missed" is
# an absence rather than a competing category. Width and opacity carry the same
# distinction, so the panels survive greyscale printing and colour-blind
# readers without relying on hue.
C_HIT = "#00838F"
C_MISS = "#9E9E9E"
C_DET = "#C62828"
C_SKY = "#607D8B"


def load(path: Path):
    rows = {"blind": [], "directed": []}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows[r["mode"]].append(
                (float(r["x"]) / 100.0, float(r["y"]) / 100.0,
                 float(r["z"]) / 100.0, float(r["dx"]), float(r["dy"]),
                 float(r["dz"]), int(r["hit"])))
    return rows


def entry_t(x, y, z, dx, dy, dz):
    """Ray parameter at which the ray first enters the detector box.

    Slab method, the same one UCMuGen's BoxDetector uses. Returns None if the
    ray misses, which lets the caller fall back to drawing it to the floor.
    """
    lo = (-DET_HX, -DET_HY, -DET_HZ)
    hi = (DET_HX, DET_HY, DET_HZ)
    t_near, t_far = -1e300, 1e300
    for o, d, a, b in ((x, dx, lo[0], hi[0]), (y, dy, lo[1], hi[1]),
                       (z, dz, lo[2], hi[2])):
        if abs(d) < 1e-15:
            if o < a or o > b:
                return None
            continue
        ta, tb = (a - o) / d, (b - o) / d
        if ta > tb:
            ta, tb = tb, ta
        t_near = max(t_near, ta)
        t_far = min(t_far, tb)
        if t_near > t_far:
            return None
    return t_near if t_far >= max(t_near, 0.0) and t_near > 0 else None


def box_faces(hx, hy, hz):
    x0, x1, y0, y1, z0, z1 = -hx, hx, -hy, hy, -hz, hz
    c = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    idx = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
           (2, 3, 7, 6), (1, 2, 6, 5), (0, 3, 7, 4)]
    return [[c[i] for i in f] for f in idx]


def draw(ax, rays, title, subtitle):
    # Sky plane: a thin outline, deliberately recessive.
    ax.add_collection3d(Poly3DCollection(
        [[(-SKY_HALF, -SKY_HALF, SKY_Z), (SKY_HALF, -SKY_HALF, SKY_Z),
          (SKY_HALF, SKY_HALF, SKY_Z), (-SKY_HALF, SKY_HALF, SKY_Z)]],
        facecolors=C_SKY, alpha=0.06, edgecolors=C_SKY, linewidths=0.8,
        zorder=1))

    # Rays. A ray that reaches the detector is truncated at the face it enters,
    # not continued to the floor: drawn through, the population appears to
    # converge on a point below the box and the box itself disappears behind
    # the bundle. Misses run to a floor below the plate so the sky coverage of
    # panel (a) is visible.
    misses, hits = [], []
    for x, y, z, dx, dy, dz, hit in rays:
        if dz >= 0:
            continue
        t = entry_t(x, y, z, dx, dy, dz) if hit else None
        if t is None:
            t = (Z_FLOOR - z) / dz
        seg = [(x, y, z), (x + dx * t, y + dy * t, z + dz * t)]
        (hits if hit else misses).append(seg)

    ax.add_collection3d(Line3DCollection(
        misses, colors=C_MISS, linewidths=0.4, alpha=0.35, zorder=2))
    ax.add_collection3d(Line3DCollection(
        hits, colors=C_HIT, linewidths=0.9, alpha=0.95, zorder=3))

    ax.add_collection3d(Poly3DCollection(
        box_faces(DET_HX, DET_HY, DET_HZ), facecolors=C_DET, alpha=0.95,
        edgecolors="#7F1D1D", linewidths=0.6, zorder=4))

    ax.set_xlim(-SKY_HALF, SKY_HALF)
    ax.set_ylim(-SKY_HALF, SKY_HALF)
    ax.set_zlim(Z_FLOOR, SKY_Z)
    ax.set_box_aspect((1, 1, 0.72))
    ax.view_init(elev=16, azim=-58)

    # Recessive frame: the geometry is the message, the axes are scaffolding.
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.0)
        axis._axinfo["grid"]["color"] = (0.85, 0.85, 0.85, 0.5)
        axis._axinfo["grid"]["linewidth"] = 0.4
    ax.set_xlabel("x [m]", labelpad=-6, fontsize=8)
    ax.set_ylabel("y [m]", labelpad=-6, fontsize=8)
    ax.set_zlabel("z [m]", labelpad=-6, fontsize=8)
    ax.tick_params(labelsize=7, pad=-2)

    n_hit = len(hits)
    n_tot = len(hits) + len(misses)
    ax.set_title(f"{title}\n{subtitle}: {n_hit} of {n_tot} reach the detector",
                 fontsize=9.5, pad=-2)
    return n_hit, n_tot


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/rays.csv")
    if not src.exists():
        sys.exit(f"no ray file: {src}\nbuild and run ucmugen/validation/"
                 "dump_rays.cc first (see the docstring)")
    rows = load(src)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })
    fig = plt.figure(figsize=(7.2, 3.5))
    ax1 = fig.add_subplot(121, projection="3d", computed_zorder=False)
    ax2 = fig.add_subplot(122, projection="3d", computed_zorder=False)

    b_hit, b_tot = draw(ax1, rows["blind"], "(a) sampling over the whole sky",
                        "blind")
    d_hit, d_tot = draw(ax2, rows["directed"], "(b) detector-directed sampling",
                        "directed")

    handles = [
        plt.Line2D([], [], color=C_HIT, lw=1.6, label="reaches the detector"),
        plt.Line2D([], [], color=C_MISS, lw=1.0, alpha=0.6, label="misses"),
        plt.Line2D([], [], color=C_DET, lw=3.0, alpha=0.6,
                   label="detector, $1\\times1\\times0.2$ m"),
        plt.Line2D([], [], color=C_SKY, lw=1.0,
                   label="generation surface, $10\\times10$ m"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.02))

    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.14,
                        wspace=0.02)
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, bbox_inches="tight")
    print(f"[OK] wrote {OUT_PDF}")
    print(f"[OK] wrote {OUT_PNG}")
    print(f"     blind    {b_hit}/{b_tot}")
    print(f"     directed {d_hit}/{d_tot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
