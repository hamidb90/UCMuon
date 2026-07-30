"""Two-sample agreement tests for validating the C++ port against Fortran.

Why not p-values
----------------
At the sample sizes needed to validate a generator (1e5 to 1e7 events), a
Kolmogorov-Smirnov p-value is not a useful pass criterion.  Two samples drawn
from *identical* distributions produce p-values uniform on [0, 1], so a
"p > 0.05" gate fails 5% of the time by construction, and across a matrix of 35
configurations something fails on essentially every run.  At large N the test
also becomes sensitive to differences far below anything physically meaningful,
so a correct port can fail on a 1e-4 relative shift.

Why not hand-picked tolerances either
-------------------------------------
The obvious alternative, "require D < 0.005", is just as bad: the right bound
depends on the statistic, the sample size, and the shape of the distribution.
The standard deviation of a steeply falling power-law momentum spectrum is
dominated by its rare high-momentum tail, and fluctuates by more than 1% between
independent runs of the *same* model at N = 5e4.  A 0.5% tolerance on it fails
100% of the time on correct code.

What this module does instead
-----------------------------
Every statistic is calibrated against the null hypothesis by measurement:

1. :func:`calibrate` runs the *same* model many times with different seeds and
   records each probe's statistic over independent pairs.
2. The threshold for that probe is a high quantile of what it observed, with
   headroom.
3. A comparison passes when every probe falls inside its measured band.

The pass criterion is then a statement with a concrete meaning: *the port
differs from the Fortran by no more than two Fortran runs differ from each
other.*  No tolerance in this file is chosen by hand, and
:func:`sensitivity_check` verifies the resulting suite still detects injected
errors, so the thresholds cannot be silently loose.

p-values are reported alongside, because they are what a reviewer looks for.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy import stats as sps

Sample = dict[str, np.ndarray]

# Scalar quantities compared for every configuration.
SCALAR_QUANTITIES = ("p", "cos_theta", "phi", "x", "y")

# Calibration settings.
#
# The threshold for each probe is SAFETY times the *median* null statistic.
# The median is used rather than a high quantile because a 99th percentile
# estimated from a handful of trials is essentially the maximum of those
# trials, which is far too noisy: measured at N = 5e4 and N = 2e5, such
# thresholds failed to scale as 1/sqrt(N), varying between 1.2x and 3.3x
# instead of the expected 2.0x.
#
# SAFETY = 3.0 is calibrated against the two-sample Kolmogorov distribution:
# over 4000 replicates at n = m = 2e4, P(D > 3 x median) was below 1/4000,
# while P(D > 2 x median) was 7e-3, which would be far too tight across a
# 35-configuration matrix.
DEFAULT_RUNS = 12
SAFETY = 3.0
MAX_PAIRS = 60


# ----------------------------------------------------------------------
# Probes
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Probe:
    """One named statistic computed from a pair of samples.

    ``fn`` returns ``(statistic, pvalue, note)``. The statistic must be
    non-negative and increase with disagreement, so that a single upper
    threshold is meaningful.
    """

    name: str
    fn: Callable[[Sample, Sample], tuple[float, float, str]]

    def __call__(self, a: Sample, b: Sample) -> tuple[float, float, str]:
        return self.fn(a, b)


def _ks(quantity: str) -> Probe:
    def fn(a, b):
        res = sps.ks_2samp(a[quantity], b[quantity], method="asymp")
        return float(res.statistic), float(res.pvalue), ""
    return Probe(f"ks:{quantity}", fn)


def _relative(quantity: str, stat_name: str,
              estimator: Callable[[np.ndarray], float]) -> Probe:
    """Relative difference of a location or scale estimator."""
    def fn(a, b):
        va, vb = estimator(a[quantity]), estimator(b[quantity])
        denom = max(abs(va), abs(vb), 1e-300)
        return (float(abs(va - vb) / denom), float("nan"),
                f"{va:.6g} vs {vb:.6g}")
    return Probe(f"{stat_name}:{quantity}", fn)


def _iqr(v: np.ndarray) -> float:
    q1, q3 = np.percentile(v, [25.0, 75.0])
    return float(q3 - q1)


def _charge() -> Probe:
    """Two-proportion z-test on the mu+ fraction.

    A KS test on a two-valued variable has essentially no power, so the charge
    ratio needs its own probe. This matters more here than for EcoMug, whose
    charge ratio is a constant; UCMuon's is momentum dependent, so an error in
    it shows up only as a subtle distortion of the charge split with momentum.
    """
    def fn(a, b):
        pa_arr, pb_arr = np.asarray(a["charge"]) > 0, np.asarray(b["charge"]) > 0
        na, nb = pa_arr.size, pb_arr.size
        pa, pb = pa_arr.mean(), pb_arr.mean()
        pooled = (pa_arr.sum() + pb_arr.sum()) / (na + nb)
        se = np.sqrt(pooled * (1 - pooled) * (1 / na + 1 / nb))
        z = 0.0 if se == 0 else (pa - pb) / se
        return float(abs(z)), float(2 * sps.norm.sf(abs(z))), f"{pa:.5f} vs {pb:.5f}"
    return Probe("charge:mu+ fraction", fn)


def _chi2_2d(qx: str = "log10_p", qy: str = "cos_theta",
             bins: int = 24, min_count: int = 10) -> Probe:
    """Chi-squared on the joint 2D density of (log10 p, cos theta).

    A generator can match every 1D marginal and still get the *correlation*
    between momentum and zenith angle wrong. That correlation is the entire
    content of the Guan and Frosin self-consistent angular models, so it needs
    a probe that a pair of 1D tests cannot substitute for.
    """
    def fn(a, b):
        ax, ay, bx, by = a[qx], a[qy], b[qx], b[qy]
        ex = np.histogram_bin_edges(np.concatenate([ax, bx]), bins=bins)
        ey = np.histogram_bin_edges(np.concatenate([ay, by]), bins=bins)
        ha, _, _ = np.histogram2d(ax, ay, bins=[ex, ey])
        hb, _, _ = np.histogram2d(bx, by, bins=[ex, ey])
        na, nb = ha.sum(), hb.sum()
        keep = (ha + hb) >= min_count
        ha, hb = ha[keep], hb[keep]
        if ha.size < 2:
            return 0.0, 1.0, "too few populated bins"
        # Two-sample chi2 for unequal sample sizes (PDG / Bityukov form).
        ka, kb = np.sqrt(nb / na), np.sqrt(na / nb)
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = float(np.nansum((ka * ha - kb * hb) ** 2 / (ha + hb)))
        dof = int(ha.size - 1)
        return chi2 / dof, float(sps.chi2.sf(chi2, dof)), f"dof={dof}"
    return Probe(f"chi2:({qx},{qy})", fn)


def default_probes(angular_mode: int = 2) -> list[Probe]:
    """The probe set applied to every configuration.

    Angular mode 1 is vertical-only: cos(theta) is identically 1 and phi is
    arbitrary, so the probes that depend on them are dropped rather than
    left to compare two constants.
    """
    quantities = list(SCALAR_QUANTITIES)
    if angular_mode == 1:
        quantities = [q for q in quantities if q not in ("cos_theta", "phi")]

    probes = [_ks(q) for q in quantities]
    probes += [
        _relative("p", "median", lambda v: float(np.median(v))),
        _relative("p", "iqr", _iqr),
        _relative("p", "mean", lambda v: float(np.mean(v))),
        _charge(),
    ]
    if angular_mode != 1:
        probes.append(_chi2_2d())
        probes.append(_relative("cos_theta", "mean", lambda v: float(np.mean(v))))
    return probes


# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
@dataclass
class Result:
    probe: str
    statistic: float
    pvalue: float
    threshold: float | None
    note: str = ""

    @property
    def passed(self) -> bool:
        if self.threshold is None:
            return True                       # informational only
        return bool(self.statistic <= self.threshold)

    def __str__(self) -> str:
        mark = "ok  " if self.passed else "FAIL"
        thr = "  n/a  " if self.threshold is None else f"{self.threshold:.5f}"
        p = "   -   " if np.isnan(self.pvalue) else f"{self.pvalue:.3g}"
        return (f"  [{mark}] {self.probe:<22s} stat={self.statistic:<10.5f} "
                f"thr={thr}  p={p}" + (f"   ({self.note})" if self.note else ""))


@dataclass
class Report:
    label: str
    results: list[Result] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(r.passed for r in self.results)

    @property
    def failures(self) -> list[Result]:
        return [r for r in self.results if not r.passed]

    def __str__(self) -> str:
        head = f"{'PASS' if self.passed else 'FAIL'}  {self.label}"
        return "\n".join([head, *(str(r) for r in self.results)])


# ----------------------------------------------------------------------
# Calibration and comparison
# ----------------------------------------------------------------------
def calibrate(sample_fn: Callable[[int], Sample],
              probes: Sequence[Probe],
              n_runs: int = DEFAULT_RUNS,
              safety: float = SAFETY,
              max_pairs: int = MAX_PAIRS,
              seed0: int = 1_000_000,
              return_observations: bool = False,
              return_samples: bool = False):
    """Measure the null band for every probe.

    ``sample_fn(seed)`` returns one sample of the model under test, as produced
    by :func:`ucmuref.fortran.derived`. ``n_runs`` independent samples are drawn
    and every distinct pair of them is compared, giving ``n_runs*(n_runs-1)/2``
    null observations from only ``n_runs`` generator runs. The pairs are not
    mutually independent, but the median of the resulting statistics is a far
    more stable estimator than one built from ``n_runs/2`` disjoint pairs, and
    it costs half as many runs.

    Each probe's threshold is ``safety`` times the median observed statistic.
    See the notes on SAFETY at the top of this module.

    Cost is ``n_runs`` generator runs, so callers should cache the result per
    (model, sample size) rather than recalibrating per comparison.
    """
    if n_runs < 4:
        raise ValueError("n_runs must be >= 4 to estimate a median")

    samples = [sample_fn(seed0 + i) for i in range(n_runs)]
    pairs = [(i, j) for i in range(n_runs) for j in range(i + 1, n_runs)]
    if len(pairs) > max_pairs:                 # keep calibration cost bounded
        step = len(pairs) / max_pairs
        pairs = [pairs[int(k * step)] for k in range(max_pairs)]

    observed: dict[str, list[float]] = {p.name: [] for p in probes}
    for i, j in pairs:
        for probe in probes:
            stat, _, _ = probe(samples[i], samples[j])
            observed[probe.name].append(stat)

    thresholds = {name: float(np.median(v) * safety)
                  for name, v in observed.items() if v}
    extras: list = []
    if return_observations:
        extras.append(observed)
    if return_samples:
        extras.append(samples)
    return (thresholds, *extras) if extras else thresholds


def compare(a: Sample, b: Sample, label: str,
            thresholds: dict[str, float],
            probes: Sequence[Probe] | None = None,
            angular_mode: int = 2) -> Report:
    """Run every probe on a pair of samples and judge against the null band."""
    probes = probes if probes is not None else default_probes(angular_mode)
    rep = Report(label)
    for probe in probes:
        stat, pval, note = probe(a, b)
        rep.results.append(Result(probe.name, stat, pval,
                                  thresholds.get(probe.name), note))
    return rep


# ----------------------------------------------------------------------
# Validating the validator
# ----------------------------------------------------------------------
def perturb(sample: Sample, p_scale: float = 1.0, p_index_shift: float = 0.0,
            cos_theta_power: float = 1.0, charge_shift: float = 0.0,
            seed: int = 0) -> Sample:
    """Apply a known distortion to a sample, for sensitivity testing.

    ``p_scale``          multiplies all momenta (a calibration error).
    ``p_index_shift``    reweights by p**shift (a spectral index error, the most
                         physically likely mistake in a port).
    ``cos_theta_power``  raises cos(theta) to a power (an angular model error).
    ``charge_shift``     moves the mu+ fraction by an absolute amount.

    Every distortion is exactly the identity at zero magnitude. That matters for
    :func:`resolving_power`, which bisects towards zero: a distortion that
    injected noise of its own at zero magnitude would appear detectable at
    arbitrarily small settings and report a resolving power far better than the
    suite really has. In particular the spectral-index shift uses rejection
    sampling rather than resampling with replacement, since the latter is a
    bootstrap and perturbs the sample even at shift = 0.
    """
    rng = np.random.default_rng(seed)
    out = {k: np.array(v, copy=True) for k, v in sample.items()}

    if p_index_shift != 0.0:
        w = out["p"] ** p_index_shift
        w /= w.max()                       # accept probability, max 1
        keep = rng.random(w.size) < w
        out = {k: v[keep] for k, v in out.items()}

    if p_scale != 1.0:
        for k in ("p", "E"):
            if k in out:
                out[k] = out[k] * p_scale
        out["log10_p"] = np.log10(out["p"])

    if cos_theta_power != 1.0 and "cos_theta" in out:
        out["cos_theta"] = np.clip(out["cos_theta"], 0.0, 1.0) ** cos_theta_power

    if charge_shift != 0.0 and "charge" in out:
        flip = rng.random(out["charge"].size) < abs(charge_shift)
        target = 1.0 if charge_shift > 0 else -1.0
        out["charge"] = np.where(flip, target, out["charge"])

    if "p" in out:
        out["log10_p"] = np.log10(out["p"])
    return out


def sensitivity_check(reference: Sample, other: Sample,
                      thresholds: dict[str, float],
                      angular_mode: int = 2,
                      distortions: dict[str, dict] | None = None) -> dict[str, bool]:
    """Confirm the calibrated thresholds still catch injected errors.

    Null calibration can only ever make a suite too loose, never too tight, so a
    calibrated suite is worthless without this counterpart: it verifies that
    realistic porting mistakes are still detected.

    ``other`` must be *statistically independent* of ``reference``. Perturbing a
    sample against itself would share its fluctuations, which cancel in the
    comparison and make the suite look more sensitive than it is.

    Returns a mapping from distortion name to whether the suite caught it.
    """
    distortions = distortions or DEFAULT_DISTORTIONS
    probes = default_probes(angular_mode)
    caught = {}
    for name, kw in distortions.items():
        bad = perturb(other, seed=abs(hash(name)) % (2**32), **kw)
        rep = compare(reference, bad, name, thresholds, probes, angular_mode)
        caught[name] = not rep.passed
    return caught


# Distortions standing in for the porting mistakes that actually happen:
# a units or mass slip (scale), a mis-transcribed spectral index, a wrong
# angular exponent, and a broken charge-ratio parametrisation.
DEFAULT_DISTORTIONS = {
    "momentum scale": dict(p_scale=1.005),
    "spectral index": dict(p_index_shift=0.05),
    "angular exponent": dict(cos_theta_power=1.05),
    "charge fraction": dict(charge_shift=0.01),
}

# Which perturb() keyword each distortion varies, and how a magnitude maps onto
# it, for the bisection in resolving_power().
_MAGNITUDE_MAP: dict[str, Callable[[float], dict]] = {
    "momentum scale": lambda m: dict(p_scale=1.0 + m),
    "spectral index": lambda m: dict(p_index_shift=m),
    "angular exponent": lambda m: dict(cos_theta_power=1.0 + m),
    "charge fraction": lambda m: dict(charge_shift=m),
}


def resolving_power(pairs: Sequence[tuple[Sample, Sample]],
                    thresholds: dict[str, float],
                    angular_mode: int = 2,
                    lo: float = 1e-4, hi: float = 0.5,
                    iterations: int = 12,
                    return_spread: bool = False):
    """Smallest injected distortion of each kind that the suite still catches.

    This is the number that makes the validation claim concrete. Instead of
    "the port agrees with the Fortran", the suite can state "at this sample
    size, any momentum-scale error above 0.5%, or angular-exponent error above
    0.7%, would have been detected". Anything smaller is below the resolving
    power, and is reported as such rather than claimed as agreement.

    ``pairs`` is a list of ``(reference, other)`` sample pairs, each
    statistically independent of the other member of its pair. The detection
    boundary is bisected separately for each pair and the **median** across
    pairs is returned: a single pair carries its own fluctuation, which makes
    the measurement vary by tens of percent and destroys the expected
    1/sqrt(N) scaling. Passing the samples already drawn by :func:`calibrate`
    costs no extra generator runs.

    Returned magnitudes are as defined by ``_MAGNITUDE_MAP``. ``inf`` means the
    distortion escaped detection even at ``hi``.
    """
    if not pairs:
        raise ValueError("need at least one (reference, other) pair")
    probes = default_probes(angular_mode)
    per_pair: dict[str, list[float]] = {}

    for k, (reference, other) in enumerate(pairs):
        for name, make_kw in _MAGNITUDE_MAP.items():
            if name == "angular exponent" and angular_mode == 1:
                continue                   # no angular spread to distort

            def caught(mag: float, _n=name, _m=make_kw, _o=other, _r=reference,
                       _k=k) -> bool:
                bad = perturb(_o, seed=(abs(hash(_n)) + _k) % (2**32), **_m(mag))
                return not compare(_r, bad, _n, thresholds, probes,
                                   angular_mode).passed

            if not caught(hi):
                mag = float("inf")
            elif caught(lo):
                mag = lo                   # already resolved at the floor
            else:
                a, b = lo, hi
                for _ in range(iterations):
                    mid = float(np.sqrt(a * b))   # geometric: spans decades
                    if caught(mid):
                        b = mid
                    else:
                        a = mid
                mag = b
            per_pair.setdefault(name, []).append(mag)

    median = {n: float(np.median(v)) for n, v in per_pair.items()}
    if return_spread:
        spread = {n: (float(np.min(v)), float(np.max(v)))
                  for n, v in per_pair.items()}
        return median, spread
    return median
