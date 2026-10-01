#!/usr/bin/env python3
"""Angular mode 6 of the Fortran generator samples J(p, θ) · max(0, -n·d).

Each case runs the Fortran generator (bin/ucmuon_gen_omp) with angular mode 6
and ucmugen::Generator, the reference sampler, on the same spectrum, surface
and detector, and compares the event fractions in bins of total energy and
zenith angle (the bins of the production check: E 15-30, 30-100, 100-1500 GeV;
θ 0-20, 20-40, 40-60, 60-85°), plus the direction components that a rotated
or tilted source would get wrong.

Pass criterion per case: every bin within 4σ (pull), and a χ² p-value above
1e-3 over the bins. The threshold is set for CI: χ²/ndf < 2, the first choice,
fails a correct sampler 2.4 % of the time at 11 degrees of freedom, which is
about one run in eight over the five cases; with p > 1e-3, and the 4σ pull
checks counted too, a correct sampler fails about one run in a hundred.
The legacy mode 4 is run on the first case as a control, and must FAIL the
same criterion: it is the sampler this mode replaces (χ² ≈ 900 for 11).

    python3 tests/flux/test_joint_sampling.py [--quick]
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _harness as H  # noqa: E402

QUICK = "--quick" in sys.argv
N = 150_000 if QUICK else 600_000
failures: list[str] = []


def chi2_pvalue(chi2: float, ndf: int) -> float:
    """Upper-tail χ² probability, Wilson-Hilferty (good to ~1e-4 here; no scipy)."""
    k = max(ndf, 1)
    z = ((chi2 / k) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k))) / math.sqrt(2.0 / (9.0 * k))
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def check_bins(name, fr, uc, e_bins=H.E_BINS, theta_bins=H.THETA_BINS, expect_fail=False):
    Hf, Hu = H.binned(fr, e_bins, theta_bins), H.binned(uc, e_bins, theta_bins)
    ratio, err, pull = H.compare_bins(Hf, Hu)
    ok_bins = Hu > 25
    chi2 = float(np.sum(pull[ok_bins] ** 2))
    ndf = int(ok_bins.sum()) - 1
    worst = float(np.max(np.abs(pull[ok_bins])))
    pval = chi2_pvalue(chi2, ndf)
    ok = worst < 4.0 and pval > 1e-3
    print(f"\n{name}")
    print("   E \\ θ      " + "".join(f"{a:>4.0f}-{b:<4.0f}  " for a, b in zip(theta_bins[:-1], theta_bins[1:])))
    for i in range(len(e_bins) - 1):
        cells = "".join(f"{ratio[i,j]:6.3f}±{err[i,j]:.3f} " if ok_bins[i, j] else "    --      "
                        for j in range(len(theta_bins) - 1))
        print(f"   {e_bins[i]:5.0f}-{e_bins[i+1]:<5.0f} {cells}")
    tag = "ok  " if ok != expect_fail else "FAIL"
    what = "fails, as it must" if expect_fail and not ok else ("agrees" if ok else "DISAGREES")
    print(f"  [{tag}] Fortran/UCMuGen: max |pull| {worst:.1f}, χ²/ndf {chi2:.1f}/{ndf} (p = {pval:.2g})  -> {what}")
    if ok == expect_fail:
        failures.append(name)
    return ratio


def check_components(name, fr, uc, comp, bins):
    a = np.histogram(fr.d[:, comp], bins=bins)[0]
    b = np.histogram(uc.d[:, comp], bins=bins)[0]
    fa, fb = a / a.sum(), b / b.sum()
    m = b > 25
    pull = (fa - fb)[m] / np.sqrt(fa[m] * (1 - fa[m]) / a.sum() + fb[m] * (1 - fb[m]) / b.sum())
    ok = np.max(np.abs(pull)) < 4.0
    print(f"  [{'ok  ' if ok else 'FAIL'}] {name}: direction component {'xyz'[comp]} "
          f"in {len(bins)-1} bins, max |pull| {np.max(np.abs(pull)):.1f}")
    if not ok:
        failures.append(f"{name} component {'xyz'[comp]}")


def main() -> int:
    # 1. The production configuration: Guan, 15-1500 GeV, θ <= 85°, R = 500 m disk.
    prod = H.Config(spectrum=4, emin=15.0, emax=1500.0, theta_max_deg=85.0,
                    source_mode=1, radius_m=500.0, n=N)
    uc = H.run_ucmugen(prod)
    check_bins("1. Guan, disk R = 500 m, mode 6 (production configuration)",
               H.run_fortran(prod), uc)
    legacy = H.run_fortran(H.Config(**{**prod.__dict__, "angular_mode": 4}))
    check_bins("1b. same, legacy mode 4 (control: must fail)", legacy, uc, expect_fail=True)

    # 2. With a cylinder detector filter: a borehole-like vertical cylinder
    #    under a smaller disk, so that the hit statistics stay reasonable.
    det = H.Detector("cyl", (0.0, 0.0, -2400.0, 0.0, 0.0, -2000.0, 150.0), margin=50.0)
    cyl = H.Config(spectrum=4, emin=15.0, emax=1500.0, theta_max_deg=85.0,
                   source_mode=1, radius_m=40.0, detector=det, n=N // 3)
    check_bins("2. Guan, disk R = 40 m, cylinder filter at 20-24 m depth, mode 6",
               H.run_fortran(cyl), H.run_ucmugen(cyl))

    # 3. Reyna on a hemisphere: position-dependent normal.
    hemi = H.Config(spectrum=7, emin=2.0, emax=1000.0, theta_max_deg=80.0,
                    source_mode=3, radius_m=20.0, n=N)
    eb = (2.0, 10.0, 50.0, 1000.0)
    fr, uc = H.run_fortran(hemi), H.run_ucmugen(hemi)
    check_bins("3. Reyna, hemisphere R = 20 m, mode 6", fr, uc, e_bins=eb,
               theta_bins=(0, 20, 40, 60, 80))

    # 4. Frosin on a vertical plane (XZ): the sky must stay fixed, so muons
    #    travel downward and only those with d_y < 0 cross the +y face.
    vert = H.Config(spectrum=5, emin=2.0, emax=1000.0, theta_max_deg=85.0,
                    source_mode=2, source_plane=2, half_m=10.0, n=N)
    fr, uc = H.run_fortran(vert), H.run_ucmugen(vert)
    check_bins("4. Frosin, vertical XZ rectangle, mode 6", fr, uc, e_bins=eb)
    check_components("4", fr, uc, 1, np.linspace(-1, 0, 11))
    check_components("4", fr, uc, 2, np.linspace(-1, 0, 11))
    if np.any(fr.d[:, 1] > 1e-9) or np.any(fr.d[:, 2] > 1e-9):
        print("  [FAIL] 4: muons going up or leaving through the face")
        failures.append("4 directions")

    # 5. CosmoALEPH on a disk tilted by 30°: projection onto the tilted normal.
    tilt = H.Config(spectrum=1, emin=100.0, emax=2500.0, theta_max_deg=85.0,
                    source_mode=1, radius_m=10.0, tilt_deg=30.0, tilt_az_deg=40.0, n=N)
    fr, uc = H.run_fortran(tilt), H.run_ucmugen(tilt)
    check_bins("5. CosmoALEPH, disk tilted 30° (az 40°), mode 6", fr, uc,
               e_bins=(100.0, 200.0, 500.0, 2500.0))
    check_components("5", fr, uc, 0, np.linspace(-1, 1, 11))
    check_components("5", fr, uc, 1, np.linspace(-1, 1, 11))

    # 6. PARMA (when the JAEA data and UCMuGen_PARMA.h are present).
    if H.PARMA_DIR.exists() and (H.ROOT / "ucmugen/include/UCMuGen_PARMA.h").exists():
        pc = H.Config(spectrum=3, emin=2.0, emax=1000.0, theta_max_deg=80.0,
                      source_mode=1, radius_m=10.0, n=N // 2)
        fr = H.run_fortran(pc)
        import re
        rc = float(re.search(r"PARMA cutoff rigid:\s*([0-9.]+)", fr.stdout).group(1))
        dd = float(re.search(r"PARMA atm. depth:\s*([0-9.]+)", fr.stdout).group(1))
        w = float(re.search(r"PARMA W index:\s*([-0-9.]+)", fr.stdout).group(1))
        uc = H.run_ucmugen(pc, parma_site=(w, rc, dd, 0.0))
        check_bins("6. PARMA (sea level, Louvain), disk, mode 6", fr, uc, e_bins=eb,
                   theta_bins=(0, 20, 40, 60, 80))
    else:
        print("\n6. PARMA: skipped (data/EXPACS/parma or UCMuGen_PARMA.h absent)")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {failures}")
        return 1
    print("ALL AGREE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
