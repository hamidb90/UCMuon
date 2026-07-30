"""The matrix of configurations the C++ port must reproduce.

A full cross product of 7 spectra x 5 angular models x 3 surfaces is 105
configurations, which is both slow and redundant: the spectrum, the angular
model and the surface are independent code paths in the generator, so a
spanning set exercises every path with far fewer runs than a cross product.

The set below is built from four groups:

* **spectra**    every spectrum against a fixed angular model and surface,
* **angular**    every angular model against a fixed spectrum and surface,
* **surfaces**   every surface (and plane orientation) against a fixed spectrum,
* **edges**      the cases most likely to break a port: mono-energetic input,
                 a narrow energy window, a tilted plane, and the natural
                 spectrum/angular pairings the models were fitted with.

PARMA (spectrum 3) is excluded here. It needs the EXPACS data tables, is
licensed separately, and is validated by its own module.
"""

from __future__ import annotations

from .fortran import GenConfig

# Sample size for the reference runs. At this N the validation suite detects
# momentum-scale, spectral-index, angular-exponent and charge-fraction errors
# at or above roughly 0.25%; see validation/reference/resolving_power.json and
# the measurement in validation/measure_power.py.
DEFAULT_N = 1_000_000

# Fast setting for pre-commit and smoke runs. Resolving power degrades to
# roughly 1%, which still catches gross porting mistakes.
SMOKE_N = 50_000

NON_PARMA_SPECTRA = (1, 2, 4, 5, 6, 7, 8)
ANGULAR_MODES = (1, 2, 3, 4, 5)


def build(n: int = DEFAULT_N) -> list[GenConfig]:
    """Return the spanning set of configurations, de-duplicated."""
    base = GenConfig(nmuons=n, emin=1.0, emax=1000.0, spectrum=4,
                     angular_mode=2, source_mode=1, disk_r_m=2.0,
                     theta_max_deg=70.0)

    cases: list[GenConfig] = []

    # Every spectrum, holding the angular model and surface fixed.
    for s in NON_PARMA_SPECTRA:
        cfg = base.replace(spectrum=s)
        if s == 8:
            # Cosmic electrons: 10 MeV - 1 GeV, well below the muon range.
            cfg = cfg.replace(emin=0.01, emax=1.0)
        cases.append(cfg)

    # Every angular model, holding the spectrum and surface fixed.
    for a in ANGULAR_MODES:
        cases.append(base.replace(angular_mode=a))

    # Every surface, plus the two vertical plane orientations. The projection
    # convention differs per surface, so these are the highest-risk paths.
    cases += [
        base.replace(source_mode=1),                        # horizontal disk
        base.replace(source_mode=2),                        # horizontal rect
        base.replace(source_mode=3, hemi_radius_m=2.0),     # hemisphere
        base.replace(source_mode=2, source_plane=2),        # XZ (vertical)
        base.replace(source_mode=2, source_plane=3),        # YZ (vertical)
    ]

    # Edge cases and the pairings each model was actually fitted with.
    cases += [
        # Mono-energetic: emin == emax skips the CDF entirely (see
        # build_cosmoaleph_cdf's early return), a separate code path.
        base.replace(emin=10.0, emax=10.0),
        # Narrow window: the CDF grid is 300 points over [p_min, p_max], so a
        # narrow range stresses interpolation near the grid edges.
        base.replace(emin=100.0, emax=110.0),
        # Wide window: four decades, the opposite interpolation stress.
        base.replace(emin=0.5, emax=5000.0),
        # Tilted generation plane.
        base.replace(source_mode=1, tilt_deg=30.0, tilt_az_deg=45.0),
        # Near-horizon acceptance, where cos(theta) -> 0 and the angular
        # samplers are least well behaved.
        base.replace(angular_mode=2, theta_max_deg=89.0),
        # Reyna-Bugaev with its own cos^3 angular model.
        base.replace(spectrum=7, angular_mode=5),
        # Frosin with the self-consistent P(theta|E) it was fitted with.
        base.replace(spectrum=5, angular_mode=4),
        # Guan with the self-consistent angular model, the flagship pairing.
        base.replace(spectrum=4, angular_mode=4),
        # Electrons with cos^3, their documented pairing.
        base.replace(spectrum=8, angular_mode=5, emin=0.01, emax=1.0),
    ]

    # De-duplicate while preserving order; several groups overlap at the base.
    seen: set[str] = set()
    unique: list[GenConfig] = []
    for c in cases:
        key = c.digest
        if key not in seen:
            seen.add(key)
            unique.append(c)
    return unique


def label(cfg: GenConfig) -> str:
    """Human-readable label including the details `name` omits."""
    bits = [cfg.name]
    if cfg.emin == cfg.emax:
        bits.append(f"mono E={cfg.emin:g}")
    else:
        bits.append(f"E=[{cfg.emin:g},{cfg.emax:g}]")
    if cfg.angular_mode != 1:
        bits.append(f"thmax={cfg.theta_max_deg:g}deg")
    return " ".join(bits)


if __name__ == "__main__":
    cs = build()
    print(f"{len(cs)} configurations at N={DEFAULT_N:,}\n")
    for i, c in enumerate(cs, 1):
        print(f"{i:3d}. {label(c)}")
