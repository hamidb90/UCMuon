#!/usr/bin/env python3
"""Digitize Fig. 3 of Reyna (2006), arXiv:hep-ph/0604145, from the vector PDF.

Fig. 3 shows the surface muon intensity of 11 experimental data sets scaled
to the vertical with 1/cos^3(theta), plotted against zeta = p*cos(theta),
together with Reyna's best-fit curve (the parametrization used as UCMuon
spectrum Mode 7).  The figure in references/source/Reyna-B.pdf is pure
vector graphics, so marker centres are recovered exactly from the PDF
drawing commands — no raster digitization error.

Axis calibration is a least-squares fit of the decade tick-label centres
(5 labels in x, 8 in y), giving sub-percent accuracy in both axes.

Outputs (relative to this script):
    data/reyna2006_fig3_data.csv     zeta [GeV], I_vertical [cm^-2 s^-1 sr^-1 GeV^-1]
                                     per experiment, with vertical errors where
                                     an error bar could be associated
    data/reyna2006_fig3_bestfit.csv  the published best-fit curve, digitized
    figures/qa_reyna_fig3_overlay.png  QA overlay: extracted points on the
                                       original figure

Requires: pymupdf  (pip install pymupdf)
"""

import csv
import os

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "..", "references", "source", "Reyna-B.pdf")
PAGE = 4          # 0-based page index holding Fig. 3
TOL = 0.30        # pt tolerance for coincidence tests

# ---------------------------------------------------------------------------
# Marker signatures: (stroke colour, fill colour, item-orientation string)
# Orientation string: one letter per path item, H/V/D for horizontal /
# vertical / diagonal line segments, 'c' for a Bezier arc.
# ---------------------------------------------------------------------------
SIGNATURES = {
    # name                      stroke           fill             pattern
    "nandi_sinha_0deg":   ((0.0, 0.0, 0.0), None,            "HDDDD"),
    "mars_0deg":          ((0.0, 0.0, 0.0), None,            "HVHVHVHVHVHV"),
    "okayama_0deg":       ((0.0, 0.0, 0.0), None,            "DD"),
    "okayama_30deg":      (None,            (1.0, 0.0, 0.0), "HVHV"),
    "kellogg_30deg":      ((1.0, 0.0, 0.0), None,            "DDDD"),
    "okayama_60deg":      (None,            (0.0, 1.0, 0.0), "DH"),
    "okayama_75deg":      (None,            (0.0, 0.0, 1.0), "DH"),
    "kiel_desy_75deg":    ((0.0, 0.0, 1.0), None,            "HVHVH"),
    "okayama_80deg":      (None,            (0.0, 1.0, 1.0), "HDDD"),
    "mutron_89deg":       ((1.0, 0.0, 1.0), None,            "cccc"),
}
# okayama_60/75 share the same signature; disambiguated by fill colour
# (green vs blue), already part of the key.  kellogg_75deg ('+' marker) is
# a single blue two-segment path (PAW renders this marker in blue); its
# centre is the intersection of the two strokes, handled separately below.


def orientation(item):
    if item[0] != "l":
        return item[0]
    a, b = item[1], item[2]
    dx, dy = abs(b.x - a.x), abs(b.y - a.y)
    return "H" if dy < 0.05 else "V" if dx < 0.05 else "D"


def path_signature(path):
    return "".join(orientation(it) for it in path["items"])


def calibrate(page):
    """Least-squares page-coordinate -> log10(data) mapping from the major
    tick marks (more precise than the tick-label text bounding boxes, whose
    visual centres are offset by ~0.5 pt)."""
    xticks, yticks = [], []
    for p in page.get_drawings():
        if p.get("color") != (0.0, 0.0, 0.0) or p.get("fill") is not None \
                or len(p["items"]) != 1 or p["items"][0][0] != "l":
            continue
        a, b = p["items"][0][1], p["items"][0][2]
        dx, dy = abs(b.x - a.x), abs(b.y - a.y)
        # x-axis major ticks: long vertical strokes rising from the bottom frame
        if dx < 0.05 and dy > 4.0 and max(a.y, b.y) > 349.0 and min(a.y, b.y) > 330:
            xticks.append(a.x)
        # y-axis major ticks: long horizontal strokes at the left frame edge
        if dy < 0.05 and 4.5 < dx < 10.0 and 160 < min(a.x, b.x) < 165:
            yticks.append(a.y)
    xticks = sorted(set(round(t, 2) for t in xticks))
    yticks = sorted(set(round(t, 2) for t in yticks))
    assert len(xticks) == 5, f"expected 5 x major ticks, got {xticks}"
    assert len(yticks) == 8, f"expected 8 y major ticks, got {yticks}"
    ax, bx = np.polyfit(xticks, [-1, 0, 1, 2, 3], 1)          # 10^-1 .. 10^3
    ay, by = np.polyfit(yticks, [-3, -4, -5, -6, -7, -8, -9, -10], 1)
    return (lambda x: ax * x + bx), (lambda y: ay * y + by)


def main():
    doc = fitz.open(PDF)
    page = doc[PAGE]
    paths = page.get_drawings()
    to_logx, to_logy = calibrate(page)

    # Plot frame: the largest stroked rectangle-ish extent of thin black lines
    # (known from inspection; clip slightly inside to drop axis/tick artwork).
    frame = fitz.Rect(163.0, 147.0, 447.0, 352.0)
    # Legend box (top-right): exclude its sample markers.
    legend = fitz.Rect(345.0, 165.0, 412.0, 250.0)

    def keep(rect):
        c = ((rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2)
        return frame.contains(fitz.Point(*c)) and not legend.contains(fitz.Point(*c))

    datasets = {name: [] for name in SIGNATURES}
    datasets["kellogg_75deg"] = []
    segs_v = []                      # bare black segments, for error bars
    bestfit = []

    for p in paths:
        r = p["rect"]
        sig = path_signature(p)
        col = p.get("color")
        fil = p.get("fill")

        if max(r.width, r.height) <= 8:
            for name, (sc, fc, pat) in SIGNATURES.items():
                if sig == pat and col == sc and fil == fc and keep(r):
                    if name == "kellogg_30deg" and not (2.0 < r.width < 2.8):
                        continue   # diamond is 2.4 x 3.6; reject strays
                    datasets[name].append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
                    break

        # Kellogg 75deg '+': one blue path of two crossing strokes
        if col == (0.0, 0.0, 1.0) and fil is None and sig in ("VH", "HV"):
            cx = cy = None
            for it in p["items"]:
                a, b = it[1], it[2]
                if abs(b.y - a.y) < 0.05:
                    cy = a.y                       # horizontal stroke -> y
                else:
                    cx = a.x                       # vertical stroke -> x
            if cx is not None and cy is not None and keep(fitz.Rect(cx, cy, cx, cy)):
                datasets["kellogg_75deg"].append((cx, cy))

        if col == (0.0, 0.0, 0.0) and fil is None and len(p["items"]) == 1 \
                and sig == "V":
            a, b = p["items"][0][1], p["items"][0][2]
            segs_v.append((min(a.x, b.x), min(a.y, b.y),
                           max(a.x, b.x), max(a.y, b.y)))

        # Best-fit curve: two long black polylines (~50 segments each)
        if col == (0.0, 0.0, 0.0) and fil is None and len(p["items"]) > 20 \
                and r.width > 50:
            for it in p["items"]:
                bestfit.append((it[1].x, it[1].y))
            bestfit.append((p["items"][-1][2].x, p["items"][-1][2].y))
    bestfit = sorted(set(bestfit))

    # ---- vertical error bars ------------------------------------------------
    # PAW draws them as two half-segments that stop at the marker edge (or as
    # one full segment through open markers).  Collect, per marker, the upper
    # and lower bar ends among V segments sharing the marker's x.
    def vbar(center):
        x, y = center
        hi_end = lo_end = None          # page-y of the bar extremities
        for s in segs_v:
            if abs((s[0] + s[2]) / 2 - x) > TOL:
                continue
            top, bot = s[1], s[3]       # top = smaller page-y = larger value
            if bot <= y + 0.6 and bot >= y - 4.0:          # upper half-bar
                hi_end = top if hi_end is None else min(hi_end, top)
            elif top >= y - 0.6 and top <= y + 4.0:        # lower half-bar
                lo_end = bot if lo_end is None else max(lo_end, bot)
            elif top < y - 0.6 < y + 0.6 < bot:            # full bar through
                hi_end = top if hi_end is None else min(hi_end, top)
                lo_end = bot if lo_end is None else max(lo_end, bot)
        return hi_end, lo_end

    # ---- write outputs -----------------------------------------------------
    theta = {"nandi_sinha_0deg": 0, "mars_0deg": 0, "okayama_0deg": 0,
             "okayama_30deg": 30, "kellogg_30deg": 30, "okayama_60deg": 60,
             "okayama_75deg": 75, "kellogg_75deg": 75, "kiel_desy_75deg": 75,
             "okayama_80deg": 80, "mutron_89deg": 89}

    out = os.path.join(HERE, "data", "reyna2006_fig3_data.csv")
    n_total = 0
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["experiment", "zenith_deg", "zeta_gev",
                    "intensity_cm2_s_sr_gev", "err_lo", "err_hi"])
        for name in theta:
            pts = sorted(datasets[name])
            n_total += len(pts)
            for (x, y) in pts:
                zeta = 10 ** to_logx(x)
                inten = 10 ** to_logy(y)
                hi_end, lo_end = vbar((x, y))
                err_hi = max(0.0, 10 ** to_logy(hi_end) - inten) if hi_end else ""
                err_lo = max(0.0, inten - 10 ** to_logy(lo_end)) if lo_end else ""
                w.writerow([name, theta[name], f"{zeta:.6g}", f"{inten:.6g}",
                            err_lo and f"{err_lo:.3g}", err_hi and f"{err_hi:.3g}"])
            print(f"{name:20s} {len(pts):3d} points")
    print(f"total {n_total} data points -> {out}")

    out_fit = os.path.join(HERE, "data", "reyna2006_fig3_bestfit.csv")
    with open(out_fit, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["zeta_gev", "intensity_cm2_s_sr_gev"])
        for (x, y) in bestfit:
            w.writerow([f"{10 ** to_logx(x):.6g}", f"{10 ** to_logy(y):.6g}"])
    print(f"best-fit curve: {len(bestfit)} samples -> {out_fit}")

    # ---- QA overlay --------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dpi = 200
    pix = page.get_pixmap(dpi=dpi)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    fig, ax = plt.subplots(figsize=(pix.width / dpi, pix.height / dpi), dpi=dpi)
    ax.imshow(img, extent=(0, pix.width * 72 / dpi, pix.height * 72 / dpi, 0))
    for name, pts in datasets.items():
        if pts:
            arr = np.array(pts)
            ax.plot(arr[:, 0], arr[:, 1], "o", ms=1.2, mew=0, color="orange")
    if bestfit:
        arr = np.array(bestfit)
        ax.plot(arr[:, 0], arr[:, 1], ",", color="lime")
    ax.set_axis_off()
    qa = os.path.join(HERE, "figures", "qa_reyna_fig3_overlay.png")
    fig.savefig(qa, bbox_inches="tight", pad_inches=0)
    print(f"QA overlay -> {qa}")


if __name__ == "__main__":
    main()
