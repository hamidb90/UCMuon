#!/usr/bin/env python3
"""
make_fig04_gui_screenshot.py — Convert a GUI screenshot PNG into
fig04_gui_screenshot.pdf for the CPC paper.

The screenshot itself must be captured manually (the figure shows the
live Results tab with loaded data):

  1.  bash run_gui.sh                       # launches http://localhost:8501
  2.  Run a generator + transport pass and open the Results tab.
  3.  Capture the browser window, e.g. on macOS:
        cmd-shift-4, space, click the window
      and save it as manuscript/figs/gui_screenshot.png
  4.  python3 manuscript/scripts/make_fig04_gui_screenshot.py
        [path/to/screenshot.png]

The script trims uniform border pixels and embeds the image in a
correctly sized vector PDF canvas so LaTeX scaling stays sharp.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

FIGS_DIR = Path(__file__).resolve().parent.parent / "figs"
DEFAULT_PNG = FIGS_DIR / "gui_screenshot.png"
OUT_PDF = FIGS_DIR / "fig04_gui_screenshot.pdf"


def autocrop(img: np.ndarray, tol: float = 0.02) -> np.ndarray:
    """Trim uniform-colour borders (window margins) from the screenshot."""
    rgb = img[..., :3]
    ref = rgb[0, 0]
    mask = (np.abs(rgb - ref).sum(axis=2) > tol * 3)
    if not mask.any():
        return img
    rows = np.flatnonzero(mask.any(axis=1))
    cols = np.flatnonzero(mask.any(axis=0))
    return img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]


def main():
    png = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PNG
    if not png.exists():
        sys.exit(f"[ERROR] screenshot not found: {png}\n"
                 f"Capture the GUI Results tab first (see module docstring).")

    img = mpimg.imread(str(png))
    img = autocrop(img)
    h, w = img.shape[:2]

    width_in = 6.5                       # \linewidth of the cas-sc class
    fig = plt.figure(figsize=(width_in, width_in * h / w))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis("off")
    fig.savefig(OUT_PDF, dpi=300)
    print(f"[OK] wrote {OUT_PDF}  ({w}x{h} px source)")


if __name__ == "__main__":
    main()
