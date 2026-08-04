// =============================================================================
//  UCMuGen.h -- header-only cosmic muon generator for Geant4 and standalone C++
//
//  A C++17 port of the UCMuon Fortran generator (src/generator/*.f90).
//  No dependencies. Drop this file in and include it.
//
//  Physics: 7 sea-level spectrum models and 5 angular models, a
//  momentum-dependent charge ratio, and generation on disk, rectangle,
//  hemisphere or cylinder surfaces with correct through-surface projection.
//
//  Optionally aim at a detector (box, sphere, or oriented cylinder): sampling
//  is then restricted to the cone it subtends, which for a small detector is
//  orders of magnitude faster, and rate() becomes the rate into it.
//
//  Reference implementation: UCMuon (https://github.com/hamidb90/UCMuon),
//  Zenodo concept DOI 10.5281/zenodo.20826984.
//
//  Layout of this file
//    1. units
//    2. Rng            random number sources, incl. a Fortran-exact LCG
//    3. flux           the spectrum models, as pure functions
//    4. MomentumCdf    momentum sampling (analytic or tabulated)
//    5. angular        the zenith-angle samplers
//    6. charge         the momentum-dependent mu+/mu- ratio
//    7. legacy         a driver reproducing the Fortran bit-for-bit
//
//  MIT licensed. PARMA (spectrum 3) is not in this file: it is JAEA's and is
//  non-commercial-only, so it lives in UCMuGen_PARMA.h, which installs itself
//  through flux::parma_provider(). Including that header is what takes on the
//  non-commercial restriction; this file alone never does.
// =============================================================================
#ifndef UCMUGEN_H
#define UCMUGEN_H

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <vector>

#define UCMUGEN_VERSION_MAJOR 0
#define UCMUGEN_VERSION_MINOR 1
#define UCMUGEN_VERSION_PATCH 0

namespace ucmugen {

// =============================================================================
// 1. Units and constants
// =============================================================================
namespace units {
inline constexpr double GeV = 1.0;
inline constexpr double MeV = 1.0e-3;
inline constexpr double cm = 1.0;
inline constexpr double m = 100.0;
inline constexpr double mm = 0.1;
inline constexpr double s = 1.0;
inline constexpr double hertz = 1.0;
}  // namespace units

// Values are those of the Fortran module, kept to the digit so that the port
// reproduces it exactly. MUON_MASS in particular is the Fortran's 0.10566,
// not the PDG 0.1056583755; changing it would shift every kinematic quantity.
inline constexpr double kPi = 3.141592654;
inline constexpr double kMuonMass = 0.10566;      // GeV/c^2
inline constexpr double kElectronMass = 0.000511; // GeV/c^2

enum class Spectrum {
  CosmoALEPH = 1,
  PowerLaw = 2,
  Parma = 3,          // needs UCMuGen_PARMA.h, licensed separately by JAEA
  Guan = 4,
  Frosin = 5,
  GaisserBugaev = 6,
  ReynaBugaev = 7,
  Electron = 8,
};

enum class Angular {
  Vertical = 1,
  Cos2 = 2,
  UniformCone = 3,
  GuanSelfConsistent = 4,
  Cos3 = 5,
};

inline bool is_electron(Spectrum s) { return s == Spectrum::Electron; }

inline double particle_mass(Spectrum s) {
  return is_electron(s) ? kElectronMass : kMuonMass;
}

// =============================================================================
// 2. Random number sources
// =============================================================================
class Rng {
 public:
  virtual ~Rng() = default;
  virtual double flat() = 0;                       // in (0, 1)
  double flat(double a, double b) { return a + (b - a) * flat(); }
};

/// The generator used by the UCMuon Fortran code, reproduced exactly.
///
/// Despite its Fortran name (`par_ranlux`) this is not RANLUX: it is a 64-bit
/// LCG with the MMIX constants, from which bits 40..63 are taken and scaled by
/// 2^-24 in *single* precision. All three details matter for bit-exact
/// agreement, so the arithmetic below is deliberately not "cleaned up":
///
///   - the state is unsigned, giving defined wraparound with the same bit
///     pattern as the Fortran signed integer(8) overflow;
///   - the 24-bit draw is assembled as a float, not a double, because the
///     Fortran buffer is real(4) and the narrowing is observable;
///   - a zero draw is forced to 1, which is what keeps the result inside the
///     open interval (0, 1).
///
/// Bit-exactness is what lets the validation suite compare event streams
/// directly rather than only their distributions.
class UCMuonLcg : public Rng {
 public:
  explicit UCMuonLcg(std::int32_t base_seed = 20260729, int thread_id = 0) {
    seed(base_seed, thread_id);
  }

  /// Reproduces par_init_rng() from src/generator/rng_parallel.f90.
  void seed(std::int32_t base_seed, int thread_id = 0) {
    const std::uint64_t tid = static_cast<std::uint64_t>(thread_id);
    std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::int64_t>(base_seed)) * kMul
                    + (tid + 1u) * 2654435761u
                    + kInc;
    s = s * kMul + kInc;
    s = s * kMul + kInc;
    s = s * kMul + kInc + tid;
    s = s * kMul + kInc;
    state_ = s;
  }

  /// One draw, matching par_ranlux() exactly.
  float flat32() {
    state_ = state_ * kMul + kInc;
    std::uint64_t bits = (state_ >> 40) & 0xFFFFFFull;   // ibits(seed, 40, 24)
    if (bits == 0ull) bits = 1ull;
    return static_cast<float>(bits) * kInv24;
  }

  double flat() override { return static_cast<double>(flat32()); }

  std::uint64_t state() const { return state_; }

 private:
  static constexpr std::uint64_t kMul = 6364136223846793005ull;
  static constexpr std::uint64_t kInc = 1442695040888963407ull;
  // The Fortran writes 5.96046448e-8 as a real(4) literal; 2^-24 exactly is
  // 5.9604644775390625e-8, and the two round to the same float.
  static constexpr float kInv24 = 5.96046448e-8f;

  std::uint64_t state_ = 0;
};

/// Wraps any callable returning a uniform double, so Geant4 users can hand the
/// generator G4UniformRand and keep one reproducible stream per run:
///     gen.setRng(std::make_unique<ucmugen::CallbackRng>(
///         [] { return G4UniformRand(); }));
class CallbackRng : public Rng {
 public:
  explicit CallbackRng(std::function<double()> f) : f_(std::move(f)) {}
  double flat() override { return f_(); }

 private:
  std::function<double()> f_;
};

// =============================================================================
// 3. Flux models
//
// Every function returns a differential intensity in
// cm^-2 s^-1 sr^-1 (GeV/c)^-1, per unit solid angle. No surface projection is
// folded in: that is the generation surface's responsibility, which is what
// lets any spectrum combine with any surface.
// =============================================================================
namespace flux {

// CosmoALEPH (Schmelling 2013): dN/dp = 10^3.8467 * p^-3.1952, originally in
// m^-2; the factor 1e-4 puts it in cm^-2 to match the other models.
inline constexpr double kCosmoA = 3.8467;
inline constexpr double kCosmoB = -3.1952;
inline constexpr double kCosmoM2ToCm2 = 1.0e-4;

// Guan et al. 2015 (arXiv:1509.06176), and the cos* parametrisation shared
// with Reyna-Bugaev.
inline constexpr double kP1 = 0.102573, kP2 = -0.068287, kP3 = 0.958633;
inline constexpr double kP4 = 0.0407253, kP5 = 0.817285;
inline constexpr double kDenom = 0.99144315;
inline constexpr double kEpi = 115.0, kEk = 850.0;
inline constexpr double kKaonFrac = 0.054, kPrefactor = 0.14, kIndex = -2.7;
inline constexpr double kGuanA = 3.64, kGuanB = 1.29;      // Guan
inline constexpr double kFrosinA = 3.512, kFrosinB = 1.388;  // Frosin 2025
inline constexpr double kBugaevA = 0.0, kBugaevB = 1.0;    // no atm. correction

// Reyna-Bugaev 2006 (arXiv:hep-ph/0604145), Eq. 6-7.
inline constexpr double kReynaC0 = 0.00253, kReynaC1 = 0.2455;
inline constexpr double kReynaC2 = 1.288, kReynaC3 = -0.2555, kReynaC4 = 0.0209;

inline constexpr double kElectronIndex = -3.0;

/// Effective zenith cosine accounting for Earth curvature.
inline double cos_star(double cos_theta) {
  const double numer = cos_theta * cos_theta + kP1 * kP1
                     + kP2 * std::pow(cos_theta, kP3)
                     + kP4 * std::pow(cos_theta, kP5);
  double cs = std::sqrt(numer > 0.0 ? numer : 0.0) / kDenom;
  if (cs < 0.0) cs = 0.0;
  if (cs > 1.0) cs = 1.0;
  return cs;
}

/// Gaisser pion+kaon form with the low-energy correction of Guan et al.
/// `a` and `b` select the variant: Guan, Frosin, or uncorrected Gaisser.
inline double guan(double E_GeV, double cos_theta, double a, double b) {
  const double cs = cos_star(cos_theta);
  const double E_eff = E_GeV * (1.0 + a / (E_GeV * std::pow(cs, b)));
  const double pion = 1.0 / (1.0 + 1.1 * E_GeV * cs / kEpi);
  const double kaon = kKaonFrac / (1.0 + 1.1 * E_GeV * cs / kEk);
  return kPrefactor * std::pow(E_eff, kIndex) * (pion + kaon);
}

/// Reyna-Bugaev log-polynomial vertical flux.
inline double reyna(double p_GeV, double cos_theta) {
  const double p_eff = p_GeV * cos_star(cos_theta);
  if (p_eff <= 0.0) return 0.0;
  const double lp = std::log10(p_eff);
  const double n = kReynaC1 + kReynaC2 * lp + kReynaC3 * lp * lp
                 + kReynaC4 * lp * lp * lp;
  const double phi = kReynaC0 * std::pow(p_eff, -n);
  return phi > 0.0 ? phi : 0.0;
}

inline double cosmoaleph(double p_GeV) {
  return std::pow(10.0, kCosmoA) * kCosmoM2ToCm2 * std::pow(p_GeV, kCosmoB);
}

/// The (a, b) pair a given spectrum feeds to `guan`.
inline bool guan_family(Spectrum s, double* a, double* b) {
  switch (s) {
    case Spectrum::Guan:          *a = kGuanA;   *b = kGuanB;   return true;
    case Spectrum::Frosin:        *a = kFrosinA; *b = kFrosinB; return true;
    case Spectrum::GaisserBugaev: *a = kBugaevA; *b = kBugaevB; return true;
    default: return false;
  }
}

/// Hook for a spectrum whose implementation cannot live in this file.
///
/// PARMA is the only one, and the reason is licensing rather than size: it is
/// JAEA's, non-commercial-only, so putting it here would make this header
/// non-commercial too. Including UCMuGen_PARMA.h installs a provider; without
/// it, Spectrum::Parma throws rather than silently returning a wrong number.
///
/// Threading: this is process-wide state, not per-generator. Install the
/// provider once, before any worker thread starts, and treat it as read-only
/// afterwards. Installing it from a per-thread constructor is a data race:
/// every thread would write the same std::function while the others read it
/// through intensity(). Everything else in this header is per-instance, so one
/// Generator per thread needs no further care.
using IntensityFn = std::function<double(double p_GeV, double cos_theta)>;
inline IntensityFn& parma_provider() {
  static IntensityFn f;
  return f;
}

/// Differential intensity for any spectrum, as a function of momentum.
inline double intensity(Spectrum s, double p_GeV, double cos_theta) {
  double a = 0.0, b = 0.0;
  if (guan_family(s, &a, &b)) {
    const double E = std::sqrt(p_GeV * p_GeV + kMuonMass * kMuonMass);
    // dN/dE -> dN/dp via the Jacobian dE/dp = p/E.
    return guan(E, cos_theta, a, b) * p_GeV / E;
  }
  switch (s) {
    case Spectrum::CosmoALEPH:  return cosmoaleph(p_GeV);
    case Spectrum::ReynaBugaev: return reyna(p_GeV, cos_theta);
    case Spectrum::PowerLaw:    return std::pow(p_GeV, -3.7);
    case Spectrum::Electron:    return std::pow(p_GeV, kElectronIndex);
    case Spectrum::Parma:
      if (parma_provider()) return parma_provider()(p_GeV, cos_theta);
      throw std::invalid_argument(
          "ucmugen: Spectrum::Parma needs a provider. Include "
          "UCMuGen_PARMA.h and call ucmugen::parma::install(site).");
    default:
      throw std::invalid_argument("ucmugen: spectrum has no built-in "
                                  "intensity");
  }
}

}  // namespace flux

// =============================================================================
// 4. Momentum sampling
// =============================================================================

/// Momentum sampler: analytic inverse CDF where one exists, otherwise a
/// tabulated CDF on a 300-point logarithmic grid.
///
/// The grid size, the trapezoidal integration and the *linear* interpolation
/// in the CDF all reproduce build_cosmoaleph_cdf / sample_momentum from
/// ucmuon_source_module.f90. Linear interpolation on a log grid is not the
/// most accurate choice, but changing it would change the sampled distribution
/// at the 1e-3 level, which the validation suite resolves; it is kept
/// deliberately.
class MomentumCdf {
 public:
  static constexpr int kGrid = 300;

  MomentumCdf() = default;
  MomentumCdf(double p_min, double p_max, Spectrum s) { build(p_min, p_max, s); }

  void build(double p_min, double p_max, Spectrum s) {
    p_min_ = p_min;
    p_max_ = p_max;
    spectrum_ = s;
    analytic_ = false;
    integrated_flux_ = 0.0;
    p_.clear();
    cdf_.clear();

    if (p_min >= p_max) return;            // mono-energetic: no CDF needed

    // Spectra with an exactly invertible CDF skip the table entirely.
    if (s == Spectrum::PowerLaw || s == Spectrum::Electron) {
      analytic_ = true;
      return;
    }
    // CosmoALEPH is a pure power law, so it is analytic too. The Fortran only
    // takes that branch when p_min < 10 GeV/c and otherwise builds a table;
    // both give the same distribution, and the threshold is reproduced here so
    // that the sampled stream matches draw for draw.
    if (s == Spectrum::CosmoALEPH && p_min < 10.0) {
      analytic_ = true;
      integrated_flux_ = std::pow(10.0, flux::kCosmoA) * flux::kCosmoM2ToCm2
          * (std::pow(p_max, flux::kCosmoB + 1.0)
             - std::pow(p_min, flux::kCosmoB + 1.0))
          / (flux::kCosmoB + 1.0);
      return;
    }

    p_.resize(kGrid);
    cdf_.resize(kGrid);
    std::vector<double> f(kGrid), partial(kGrid);

    const double lr = std::log(p_max / p_min);
    for (int j = 0; j < kGrid; ++j)
      p_[j] = p_min * std::exp(double(j) / double(kGrid - 1) * lr);
    for (int j = 0; j < kGrid; ++j)
      f[j] = flux::intensity(s, p_[j], 1.0);   // vertical, as the Fortran does

    partial[0] = 0.0;
    for (int j = 1; j < kGrid; ++j)
      partial[j] = partial[j - 1] + 0.5 * (f[j] + f[j - 1]) * (p_[j] - p_[j - 1]);

    integrated_flux_ = partial[kGrid - 1];
    if (integrated_flux_ <= 0.0)
      throw std::runtime_error("ucmugen: spectrum integrates to zero");
    for (int j = 0; j < kGrid; ++j) cdf_[j] = partial[j] / integrated_flux_;
  }

  /// Map a uniform deviate to a momentum.
  double sample(double u) const {
    if (p_min_ >= p_max_) return p_min_;      // mono-energetic

    if (analytic_) {
      double alpha;
      switch (spectrum_) {
        // dN/dp ~ p^-3.7 truncated. Sampling the unbounded law and clamping
        // piles the overflow into a delta spike at p_max (about 15% of events
        // for a factor-2 window), so the truncated inverse CDF is used.
        case Spectrum::PowerLaw:   alpha = -2.7; break;
        case Spectrum::Electron:   alpha = flux::kElectronIndex + 1.0; break;
        default:                   alpha = flux::kCosmoB + 1.0; break;
      }
      const double lo = std::pow(p_min_, alpha);
      const double hi = std::pow(p_max_, alpha);
      return std::pow(lo + u * (hi - lo), 1.0 / alpha);
    }

    if (u <= cdf_[0]) return p_[0];
    for (int j = 1; j < kGrid; ++j) {
      if (u <= cdf_[j]) {
        const double d = cdf_[j] - cdf_[j - 1];
        const double frac = d > 0.0 ? (u - cdf_[j - 1]) / d : 0.0;
        return p_[j - 1] + frac * (p_[j] - p_[j - 1]);
      }
    }
    return p_[kGrid - 1];
  }

  double p_min() const { return p_min_; }
  double p_max() const { return p_max_; }
  bool analytic() const { return analytic_; }
  /// Integral of the vertical intensity over [p_min, p_max]; 0 when the
  /// spectrum carries no absolute normalisation (PowerLaw, Electron).
  double integrated_flux() const { return integrated_flux_; }

 private:
  double p_min_ = 0.0, p_max_ = 0.0;
  Spectrum spectrum_ = Spectrum::Guan;
  bool analytic_ = false;
  double integrated_flux_ = 0.0;
  std::vector<double> p_, cdf_;
};

// =============================================================================
// 5. Angular samplers
//
// Each returns a zenith angle in [0, theta_max]. The number of Rng draws each
// consumes matters: it is part of what makes the port stream-compatible with
// the Fortran.
// =============================================================================
namespace angular {

/// theta ~ cos^2(theta) per unit solid angle, by rejection. Two draws per
/// trial, matching sample_cos2.
inline double cos2(Rng& rng, double theta_max) {
  constexpr double kMaxPdf = 0.385;
  for (;;) {
    const double theta = rng.flat() * theta_max;
    const double pdf = std::cos(theta) * std::cos(theta) * std::sin(theta);
    if (rng.flat() * kMaxPdf < pdf) return theta;
  }
}

/// theta ~ cos^3(theta) per unit solid angle, by exact inverse CDF.
inline double cos3(Rng& rng, double theta_max) {
  const double u = rng.flat();
  const double c4 = std::pow(std::cos(theta_max), 4.0);
  return std::acos(std::pow(1.0 - u * (1.0 - c4), 0.25));
}

/// Uniform in solid angle within the cone.
inline double uniform_cone(Rng& rng, double theta_max) {
  const double u = rng.flat();
  return std::acos(1.0 - u * (1.0 - std::cos(theta_max)));
}

/// P(theta | E) built from the Guan/Frosin flux itself, so that the angular
/// distribution stays consistent with the spectrum instead of being imposed
/// independently. This is the pairing the models were fitted with.
///
/// The extra factor of cos(theta) in the pdf is the projection onto a
/// horizontal plane, which is why this sampler is specific to that geometry.
inline double guan_self_consistent(Rng& rng, double E_GeV, double theta_max,
                                   Spectrum spectrum) {
  constexpr int kN = 100;
  double a = flux::kGuanA, b = flux::kGuanB;
  if (spectrum == Spectrum::Frosin) { a = flux::kFrosinA; b = flux::kFrosinB; }

  const double cos_min = std::cos(theta_max);
  std::array<double, kN> c{}, pdf{}, cdf{};
  for (int k = 0; k < kN; ++k)
    c[k] = cos_min + double(k) / double(kN - 1) * (1.0 - cos_min);
  for (int k = 0; k < kN; ++k)
    pdf[k] = flux::guan(E_GeV, c[k], a, b) * c[k];

  cdf[0] = 0.0;
  for (int k = 1; k < kN; ++k)
    cdf[k] = cdf[k - 1] + 0.5 * (pdf[k] + pdf[k - 1]) * (c[k] - c[k - 1]);

  const double total = cdf[kN - 1];
  if (!(total > 0.0)) return 0.0;
  for (int k = 0; k < kN; ++k) cdf[k] /= total;

  const double u = rng.flat();
  double cc = c[kN - 1];
  for (int k = 1; k < kN; ++k) {
    if (u <= cdf[k]) {
      const double d = cdf[k] - cdf[k - 1];
      const double frac = d > 0.0 ? (u - cdf[k - 1]) / d : 0.0;
      cc = c[k - 1] + frac * (c[k] - c[k - 1]);
      break;
    }
  }
  if (cc < cos_min) cc = cos_min;
  if (cc > 1.0) cc = 1.0;
  return std::acos(cc);
}

/// Dispatch on the angular model. `E_GeV` is used only by
/// Angular::GuanSelfConsistent.
inline double sample(Rng& rng, Angular mode, double theta_max, double E_GeV,
                     Spectrum spectrum) {
  switch (mode) {
    case Angular::Vertical:           return 0.0;
    case Angular::Cos2:               return cos2(rng, theta_max);
    case Angular::UniformCone:        return uniform_cone(rng, theta_max);
    case Angular::GuanSelfConsistent:
      return guan_self_consistent(rng, E_GeV, theta_max, spectrum);
    case Angular::Cos3:               return cos3(rng, theta_max);
  }
  return 0.0;
}

}  // namespace angular

// =============================================================================
// 6. Charge ratio
// =============================================================================

/// Momentum-dependent mu+/mu- ratio, as a step function of p in GeV/c.
///
/// These are measured values, not a smooth fit, which is why the ratio is not
/// monotonic (the 0.648 at 1413-1778 GeV/c is a real downward fluctuation in
/// the input data, not a transcription error). EcoMug uses a single constant
/// ratio, so this is one of the places the two generators genuinely differ.
/// Hook for a spectrum that models the two charges separately and therefore
/// knows its own ratio.
///
/// PARMA is the case: it computes mu+ and mu- fluxes independently, so its
/// ratio varies with atmospheric depth and geomagnetic cutoff as well as with
/// momentum, where the table below is a sea-level fit and flat below 112 GeV/c.
/// Using a site-aware flux and then overriding its charge split with a
/// sea-level constant discards half of what the model provides. Installed by
/// UCMuGen_PARMA.h; empty otherwise.
inline std::function<double(double p_GeV)>& charge_ratio_provider() {
  static std::function<double(double p_GeV)> f;
  return f;
}

inline double charge_ratio(double p_GeV) {
  if (p_GeV <=  112.0) return 1.252;
  if (p_GeV <=  141.0) return 1.293;
  if (p_GeV <=  178.0) return 1.259;
  if (p_GeV <=  224.0) return 1.271;
  if (p_GeV <=  282.0) return 1.239;
  if (p_GeV <=  355.0) return 1.348;
  if (p_GeV <=  447.0) return 1.541;
  if (p_GeV <=  562.0) return 1.373;
  if (p_GeV <=  708.0) return 1.243;
  if (p_GeV <=  891.0) return 1.547;
  if (p_GeV <= 1122.0) return 1.785;
  if (p_GeV <= 1413.0) return 1.361;
  if (p_GeV <= 1778.0) return 0.648;
  return 1.495;
}

/// The ratio to use for a given spectrum: the spectrum's own where it has one,
/// the sea-level table otherwise. This is what the generators call, so a
/// spectrum that knows its charge composition is never overridden by the fit.
inline double charge_ratio(Spectrum s, double p_GeV) {
  if (s == Spectrum::Parma && charge_ratio_provider())
    return charge_ratio_provider()(p_GeV);
  return charge_ratio(p_GeV);
}

// =============================================================================
// 7. Legacy driver
//
// Reproduces src/generator/ucmuon_gen_omp.f90 exactly, including the order in
// which random numbers are consumed. Its purpose is validation and backward
// compatibility with existing UCMuon runs.
//
// NOTE: like the Fortran, this applies no projection weighting to the
// generation surface. For a horizontal surface that is the conventional
// choice; for tilted or vertical surfaces it samples the sky *intensity*
// rather than the flux *through* the surface. The projection-correct API
// (ucmugen::Generator) is the one to use for new work.
// =============================================================================
namespace legacy {

enum class Source { Disk = 1, Rectangle = 2, Hemisphere = 3 };

struct Config {
  double e_min = 1.0;              // total energy, GeV
  double e_max = 1000.0;
  Spectrum spectrum = Spectrum::Guan;
  Angular angular = Angular::Cos2;
  double theta_max_deg = 70.0;
  Source source = Source::Disk;
  double radius_cm = 200.0;        // disk and hemisphere
  double half_lx_cm = 200.0;       // rectangle
  double half_ly_cm = 200.0;
  double source_z_cm = 0.0;
  std::int32_t seed = 20260729;
};

struct Muon {
  double x = 0, y = 0, z = 0;             // cm
  double px = 0, py = 0, pz = 0;          // GeV/c
  double p = 0, energy = 0;               // GeV/c, GeV (total)
  double theta = 0, phi = 0;              // rad
  int charge = 0;                         // +1 / -1
  int pdg = 0;                            // mu-=13, mu+=-13, e-=11, e+=-11
};

class Generator {
 public:
  explicit Generator(const Config& cfg) : cfg_(cfg), rng_(cfg.seed) {
    const double m = particle_mass(cfg.spectrum);
    // The Fortran takes energy limits and converts to momentum once.
    p_min_ = std::sqrt(std::max(cfg.e_min * cfg.e_min - m * m, 0.0));
    p_max_ = std::sqrt(std::max(cfg.e_max * cfg.e_max - m * m, 0.0));
    cdf_.build(p_min_, p_max_, cfg.spectrum);
    theta_max_ = cfg.theta_max_deg * kPi / 180.0;
  }

  const MomentumCdf& cdf() const { return cdf_; }
  UCMuonLcg& rng() { return rng_; }

  /// One muon. The order of draws below is load-bearing: momentum, charge,
  /// position, direction, azimuth. Reordering it would still be correct
  /// physics but would break stream compatibility with the Fortran.
  Muon generate() {
    Muon mu;
    const double m = particle_mass(cfg_.spectrum);

    // --- momentum ---
    const double p = cdf_.sample(rng_.flat());

    // --- charge and total energy ---
    if (is_electron(cfg_.spectrum)) {
      mu.charge = (rng_.flat() < 0.5) ? 1 : -1;
      mu.energy = std::sqrt(p * p + kElectronMass * kElectronMass);
      mu.pdg = (mu.charge == 1) ? -11 : 11;
    } else {
      const double r = charge_ratio(p);
      const double pos_frac = r / (1.0 + r);
      mu.charge = (rng_.flat() < pos_frac) ? 1 : -1;
      mu.energy = std::sqrt(p * p + kMuonMass * kMuonMass);
      mu.pdg = (mu.charge == 1) ? -13 : 13;
    }

    // --- position ---
    switch (cfg_.source) {
      case Source::Rectangle: {
        mu.x = cfg_.half_lx_cm * (2.0 * rng_.flat() - 1.0);
        mu.y = cfg_.half_ly_cm * (2.0 * rng_.flat() - 1.0);
        mu.z = cfg_.source_z_cm;
        break;
      }
      case Source::Hemisphere: {
        const double ph = 2.0 * kPi * rng_.flat();
        const double th = std::acos(1.0 - rng_.flat());
        const double st = std::sin(th);
        mu.x = cfg_.radius_cm * st * std::cos(ph);
        mu.y = cfg_.radius_cm * st * std::sin(ph);
        mu.z = cfg_.radius_cm * std::cos(th) + cfg_.source_z_cm;
        break;
      }
      case Source::Disk:
      default: {
        const double tc = 2.0 * kPi * rng_.flat();
        const double r = cfg_.radius_cm * std::sqrt(rng_.flat());
        mu.x = r * std::cos(tc);
        mu.y = r * std::sin(tc);
        mu.z = cfg_.source_z_cm;
        break;
      }
    }

    // --- direction ---
    const double theta = angular::sample(rng_, cfg_.angular, theta_max_,
                                         mu.energy, cfg_.spectrum);
    const double phi = 2.0 * kPi * rng_.flat();

    const double cx = std::sin(theta) * std::cos(phi);
    const double cy = std::sin(theta) * std::sin(phi);
    const double cz = -std::cos(theta);            // downward

    // --- derived kinematics ---
    mu.p = std::sqrt(std::max(mu.energy * mu.energy - m * m, 0.0));
    mu.px = mu.p * cx;
    mu.py = mu.p * cy;
    mu.pz = mu.p * cz;

    // theta and phi are recomputed from the direction, as the Fortran does
    // after its plane rotation, rather than reused from the sampler.
    double c = -cz;
    if (c < -1.0) c = -1.0;
    if (c > 1.0) c = 1.0;
    mu.theta = std::acos(c);
    if (std::fabs(std::sin(mu.theta)) < 1.0e-9) {
      mu.phi = 0.0;
    } else {
      mu.phi = std::atan2(cy, cx);
      if (mu.phi < 0.0) mu.phi += 2.0 * kPi;
    }
    return mu;
  }

 private:
  Config cfg_;
  UCMuonLcg rng_;
  MomentumCdf cdf_;
  double p_min_ = 0.0, p_max_ = 0.0, theta_max_ = 0.0;
};

}  // namespace legacy

// =============================================================================
// 8. Vectors
// =============================================================================
struct Vec3 {
  double x = 0.0, y = 0.0, z = 0.0;

  Vec3() = default;
  Vec3(double a, double b, double c) : x(a), y(b), z(c) {}

  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  double norm() const { return std::sqrt(dot(*this)); }
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  Vec3 unit() const {
    const double n = norm();
    return n > 0.0 ? Vec3{x / n, y / n, z / n} : Vec3{0.0, 0.0, 1.0};
  }
  Vec3 cross(const Vec3& o) const {
    return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
  }
};

/// Any unit vector orthogonal to n, chosen without a degenerate axis.
inline Vec3 any_perpendicular(const Vec3& n) {
  const Vec3 a = (std::fabs(n.z) < 0.9) ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
  return n.cross(a).unit();
}

// =============================================================================
// 9. Generation surfaces
//
// A surface supplies a point, an outward unit normal, and its area. The
// projection factor max(0, -n.d) is the whole of the convention: it converts a
// differential *intensity* (per unit solid angle, per unit area perpendicular
// to the muon) into a *rate through this surface*, and it rejects muons that
// would be leaving rather than entering.
//
// This is the single design decision that separates UCMuGen from EcoMug, which
// hardcodes one J' expression per geometry (EcoMug.h:975, 980, 1026) plus a
// bespoke MC integrator for each. Here the projection lives in one place, so
// every spectrum works with every surface and adding a surface needs no new
// physics.
// =============================================================================
struct SurfacePoint {
  Vec3 position;   // cm
  Vec3 normal;     // outward unit normal
};

class Surface {
 public:
  virtual ~Surface() = default;
  virtual SurfacePoint sample(Rng& rng) const = 0;
  virtual double area() const = 0;              // cm^2
  virtual const char* name() const = 0;

  /// Rate per unit intensity for a muon travelling along `d` through a face
  /// whose outward normal is `n`. Zero for outgoing muons.
  static double projection(const Vec3& n, const Vec3& d) {
    const double c = -n.dot(d);
    return c > 0.0 ? c : 0.0;
  }
};

/// Flat rectangle. Equivalent to EcoMug's "Sky" when the normal is +Z, but the
/// normal is free, so a tilted detector face is a first-class case rather than
/// a position-only transform.
class Plane : public Surface {
 public:
  Plane(double half_x_cm, double half_y_cm, Vec3 centre,
        Vec3 normal = Vec3{0.0, 0.0, 1.0})
      : hx_(half_x_cm), hy_(half_y_cm), centre_(centre), n_(normal.unit()) {
    u_ = any_perpendicular(n_);
    v_ = n_.cross(u_).unit();
  }

  SurfacePoint sample(Rng& rng) const override {
    const double a = hx_ * (2.0 * rng.flat() - 1.0);
    const double b = hy_ * (2.0 * rng.flat() - 1.0);
    return {centre_ + u_ * a + v_ * b, n_};
  }
  double area() const override { return 4.0 * hx_ * hy_; }
  const char* name() const override { return "Plane"; }

 private:
  double hx_, hy_;
  Vec3 centre_, n_, u_, v_;
};

/// Flat disk. UCMuon's default generation surface; EcoMug has no equivalent.
class Disk : public Surface {
 public:
  Disk(double radius_cm, Vec3 centre, Vec3 normal = Vec3{0.0, 0.0, 1.0})
      : r_(radius_cm), centre_(centre), n_(normal.unit()) {
    u_ = any_perpendicular(n_);
    v_ = n_.cross(u_).unit();
  }

  SurfacePoint sample(Rng& rng) const override {
    const double phi = 2.0 * kPi * rng.flat();
    const double rr = r_ * std::sqrt(rng.flat());   // sqrt for uniform area
    return {centre_ + u_ * (rr * std::cos(phi)) + v_ * (rr * std::sin(phi)), n_};
  }
  double area() const override { return kPi * r_ * r_; }
  const char* name() const override { return "Disk"; }

 private:
  double r_;
  Vec3 centre_, n_, u_, v_;
};

/// Upper hemisphere, outward radial normal.
class HSphere : public Surface {
 public:
  HSphere(double radius_cm, Vec3 centre) : r_(radius_cm), centre_(centre) {}

  SurfacePoint sample(Rng& rng) const override {
    const double phi = 2.0 * kPi * rng.flat();
    const double ct = rng.flat();                  // uniform in cos: equal area
    const double st = std::sqrt(1.0 - ct * ct);
    const Vec3 n{st * std::cos(phi), st * std::sin(phi), ct};
    return {centre_ + n * r_, n};
  }
  double area() const override { return 2.0 * kPi * r_ * r_; }
  const char* name() const override { return "HSphere"; }

 private:
  double r_;
  Vec3 centre_;
};

/// Closed vertical cylinder: lateral surface plus both caps.
///
/// EcoMug's cylinder is lateral-only (GetGenSurfaceArea returns dphi*r*h, and
/// its J' carries sin^2(theta)), so a vertical muon has J' = 0 and can never be
/// generated. That is fine for a tall thin detector and wrong for a squat one.
/// Including the caps costs nothing here because the projection factor already
/// rejects muons leaving through a face.
class Cylinder : public Surface {
 public:
  Cylinder(double radius_cm, double height_cm, Vec3 centre, bool caps = true)
      : r_(radius_cm), h_(height_cm), centre_(centre), caps_(caps) {
    const double lateral = 2.0 * kPi * r_ * h_;
    const double cap = kPi * r_ * r_;
    area_ = caps_ ? lateral + 2.0 * cap : lateral;
    p_lateral_ = lateral / area_;
  }

  SurfacePoint sample(Rng& rng) const override {
    const double phi = 2.0 * kPi * rng.flat();
    if (!caps_ || rng.flat() < p_lateral_) {
      const double z = h_ * (rng.flat() - 0.5);
      const Vec3 n{std::cos(phi), std::sin(phi), 0.0};
      return {centre_ + Vec3{n.x * r_, n.y * r_, z}, n};
    }
    const double rr = r_ * std::sqrt(rng.flat());
    const bool top = rng.flat() < 0.5;
    const double z = top ? 0.5 * h_ : -0.5 * h_;
    return {centre_ + Vec3{rr * std::cos(phi), rr * std::sin(phi), z},
            Vec3{0.0, 0.0, top ? 1.0 : -1.0}};
  }
  double area() const override { return area_; }
  const char* name() const override { return "Cylinder"; }

 private:
  double r_, h_;
  Vec3 centre_;
  bool caps_;
  double area_ = 0.0, p_lateral_ = 0.0;
};

// =============================================================================
// 10. Detectors
//
// A detector is a closed volume the user actually cares about. Attaching one
// switches the generator into *detector-directed sampling*: directions are
// drawn only inside the cone the detector subtends from the current surface
// point, rather than over the whole allowed sky. For a small detector under a
// large generation surface that is where the orders-of-magnitude speed-up
// lives, because the fraction of blind proposals that would have hit is tiny.
//
// The two shapes of the UCMuon paper are here (axis-aligned box and
// arbitrarily oriented cylinder), plus a sphere, which is the natural bounding
// volume and costs three lines.
//
// Every detector must supply a bounding sphere, and that is a correctness
// requirement rather than a convenience. Directed sampling is unbiased only if
// the cone that is sampled provably contains *every* direction that could
// reach the volume; the cone subtended by a bounding sphere is exactly such a
// cone, and no tighter construction is safe for a general shape.
// =============================================================================
class Detector {
 public:
  virtual ~Detector() = default;

  /// True if the forward ray `origin + t*dir` with t > 0 reaches the volume.
  /// An origin already inside counts as a hit.
  virtual bool intersects(const Vec3& origin, const Vec3& dir) const = 0;

  /// Centre and radius of a sphere that fully contains the volume. Must be a
  /// true bound: the directed sampler relies on it.
  virtual Vec3 boundCentre() const = 0;
  virtual double boundRadius() const = 0;

  virtual const char* name() const = 0;
};

/// Axis-aligned box, given by two opposite corners. Intersection is the slab
/// method: the ray parameter intervals of the three slabs must have a common
/// overlap that extends to positive t.
class BoxDetector : public Detector {
 public:
  BoxDetector(Vec3 corner_a, Vec3 corner_b) {
    lo_ = {std::fmin(corner_a.x, corner_b.x), std::fmin(corner_a.y, corner_b.y),
           std::fmin(corner_a.z, corner_b.z)};
    hi_ = {std::fmax(corner_a.x, corner_b.x), std::fmax(corner_a.y, corner_b.y),
           std::fmax(corner_a.z, corner_b.z)};
  }

  /// Convenience: centre plus half-sizes, which is how detectors are usually
  /// quoted.
  static std::shared_ptr<BoxDetector> centred(Vec3 centre, double hx,
                                              double hy, double hz) {
    const Vec3 h{std::fabs(hx), std::fabs(hy), std::fabs(hz)};
    return std::make_shared<BoxDetector>(centre - h, centre + h);
  }

  bool intersects(const Vec3& o, const Vec3& d) const override {
    double t_near = -kInf, t_far = kInf;
    if (!slab(o.x, d.x, lo_.x, hi_.x, t_near, t_far)) return false;
    if (!slab(o.y, d.y, lo_.y, hi_.y, t_near, t_far)) return false;
    if (!slab(o.z, d.z, lo_.z, hi_.z, t_near, t_far)) return false;
    return t_far >= std::fmax(t_near, 0.0);
  }

  Vec3 boundCentre() const override { return (lo_ + hi_) * 0.5; }
  double boundRadius() const override { return (hi_ - lo_).norm() * 0.5; }
  const char* name() const override { return "BoxDetector"; }

 private:
  /// One slab of the AABB. Returns false if the ray misses it outright.
  /// The parallel case is handled explicitly rather than by relying on
  /// infinities, which would produce 0*inf NaN when the origin sits exactly on
  /// a face plane.
  static bool slab(double o, double d, double lo, double hi, double& t_near,
                   double& t_far) {
    if (std::fabs(d) < 1e-15) return o >= lo && o <= hi;
    const double inv = 1.0 / d;
    double ta = (lo - o) * inv, tb = (hi - o) * inv;
    if (ta > tb) { const double s = ta; ta = tb; tb = s; }
    if (ta > t_near) t_near = ta;
    if (tb < t_far) t_far = tb;
    return t_near <= t_far;
  }

  static constexpr double kInf = 1e300;
  Vec3 lo_, hi_;
};

/// Sphere. Also the cheapest way to wrap a detector whose exact shape does not
/// matter for the acceptance study.
class SphereDetector : public Detector {
 public:
  SphereDetector(Vec3 centre, double radius_cm)
      : c_(centre), r_(std::fabs(radius_cm)) {}

  bool intersects(const Vec3& o, const Vec3& d) const override {
    const Vec3 m = o - c_;
    const double c0 = m.dot(m) - r_ * r_;
    if (c0 <= 0.0) return true;                 // origin inside
    const double b = m.dot(d);
    if (b > 0.0) return false;                  // pointing away
    const double disc = b * b - c0;
    return disc >= 0.0;
  }

  Vec3 boundCentre() const override { return c_; }
  double boundRadius() const override { return r_; }
  const char* name() const override { return "SphereDetector"; }

 private:
  Vec3 c_;
  double r_;
};

/// Finite cylinder of arbitrary orientation, given by its two axis endpoints
/// and a radius. Lateral surface by the quadratic, then both cap planes, which
/// is what makes a squat cylinder work as well as a tall one.
class CylinderDetector : public Detector {
 public:
  CylinderDetector(Vec3 end_a, Vec3 end_b, double radius_cm)
      : a_(end_a), r_(std::fabs(radius_cm)) {
    const Vec3 ab = end_b - end_a;
    len_ = ab.norm();
    u_ = (len_ > 0.0) ? ab * (1.0 / len_) : Vec3{0.0, 0.0, 1.0};
  }

  bool intersects(const Vec3& o, const Vec3& d) const override {
    const Vec3 w = o - a_;
    const double wu = w.dot(u_), du = d.dot(u_);
    const Vec3 wp = w - u_ * wu;               // components across the axis
    const Vec3 dp = d - u_ * du;

    // Lateral surface: |wp + t*dp|^2 = r^2.
    const double A = dp.dot(dp);
    if (A > 1e-15) {
      const double B = 2.0 * wp.dot(dp);
      const double C = wp.dot(wp) - r_ * r_;
      const double disc = B * B - 4.0 * A * C;
      if (disc >= 0.0) {
        const double sq = std::sqrt(disc);
        for (const double t : {(-B - sq) / (2.0 * A), (-B + sq) / (2.0 * A)}) {
          if (t <= 0.0) continue;
          const double s = wu + t * du;         // axial coordinate of the hit
          if (s >= 0.0 && s <= len_) return true;
        }
      }
    }

    // Cap planes at s = 0 and s = len.
    if (std::fabs(du) > 1e-15) {
      for (const double s_plane : {0.0, len_}) {
        const double t = (s_plane - wu) / du;
        if (t <= 0.0) continue;
        const Vec3 rel = w + d * t - u_ * s_plane;
        if (rel.dot(rel) <= r_ * r_) return true;
      }
    }

    // Origin already inside.
    return wu >= 0.0 && wu <= len_ && wp.dot(wp) <= r_ * r_;
  }

  Vec3 boundCentre() const override { return a_ + u_ * (0.5 * len_); }
  double boundRadius() const override {
    return std::sqrt(0.25 * len_ * len_ + r_ * r_);
  }
  const char* name() const override { return "CylinderDetector"; }

 private:
  Vec3 a_, u_;
  double len_ = 0.0, r_ = 0.0;
};

// =============================================================================
// 11. The projection-correct generator
// =============================================================================
struct Muon {
  Vec3 position;          // cm
  Vec3 direction;         // unit
  double momentum = 0.0;  // GeV/c
  double energy = 0.0;    // total, GeV
  double kinetic = 0.0;   // GeV
  int charge = 0;         // +1 / -1
  int pdg = 0;            // mu-=13, mu+=-13, e-=11, e+=-11
};

/// Samples muons crossing a generation surface, with the correct
/// through-surface weighting, for any spectrum and any surface.
///
/// Sampling strategy: momentum from the tabulated CDF of the *vertical*
/// intensity (which absorbs the steep momentum dependence and keeps the
/// acceptance high), direction uniform in solid angle over the allowed cone,
/// position uniform on the surface, then accept with probability proportional
/// to `[J(p, cos t) / J(p, 1)] * max(0, -n.d)`.
///
/// The ratio is genuinely greater than one in places: well above the pion
/// critical energy the flux rises towards the horizon roughly as 1/cos(theta),
/// so the acceptance envelope must be measured rather than assumed to be 1.
///
/// Detector-directed sampling
/// --------------------------
/// With a detector attached, directions are drawn from the cone C(x) that its
/// bounding sphere subtends from the sampled surface point x, instead of from
/// the whole allowed sky. Writing w = J(p,c)/J(p,1) and P = max(0, -n.d), the
/// rate is in both cases a single expectation over (p, x, d):
///
///   blind :  R = A * dOmega * I_p * E[ w * P * 1{hit} ]
///   direct:  R = A *          I_p * E[ Omega_c(x) * w * P * 1{hit} ]
///
/// The second line is the first with the direction integral restricted to
/// C(x): the constant dOmega is replaced by the per-trial Omega_c(x), which is
/// legitimate precisely because C(x) contains every direction that can reach
/// the volume. Both the event distribution and the rate therefore come from
/// one per-trial quantity f, computed in one place (`sampleTrial`) and used by
/// both `generate()` and `rateAndError()`. Keeping a single f is the whole
/// safety argument: an estimator that a fast path and a normalisation path
/// compute separately is one refactor away from disagreeing silently.
///
/// Directed sampling can be turned off with `setDirectedSampling(false)`,
/// which keeps the detector as a plain geometric acceptance cut. The two modes
/// must agree on `rate()` to within Monte Carlo error, and the validation
/// suite checks exactly that.
class Generator {
 public:
  Generator() : rng_(new UCMuonLcg(20260729)) {}

  // --- configuration -------------------------------------------------------
  Generator& setSpectrum(Spectrum s) { spectrum_ = s; dirty_ = true; return *this; }
  Generator& setSurface(std::shared_ptr<Surface> s) {
    surface_ = std::move(s); return *this;
  }
  Generator& setRng(std::unique_ptr<Rng> r) { rng_ = std::move(r); return *this; }
  Generator& setSeed(std::int32_t seed) {
    rng_.reset(new UCMuonLcg(seed)); seed_ = seed; return *this;
  }
  /// Total-energy limits, GeV, matching the Fortran's convention.
  Generator& setEnergyRange(double e_min, double e_max) {
    e_min_ = e_min; e_max_ = e_max; dirty_ = true; return *this;
  }
  Generator& setThetaRange(double t_min_rad, double t_max_rad) {
    theta_min_ = t_min_rad; theta_max_ = t_max_rad; dirty_ = true; return *this;
  }
  Generator& setPhiRange(double p_min_rad, double p_max_rad) {
    phi_min_ = p_min_rad; phi_max_ = p_max_rad; dirty_ = true; return *this;
  }
  /// Attach a detector. Muons whose ray misses it are never returned, and
  /// `rate()` reports the rate *into the detector* rather than through the
  /// whole surface. Pass nullptr to remove it.
  Generator& setDetector(std::shared_ptr<Detector> d) {
    detector_ = std::move(d); dirty_ = true; return *this;
  }
  /// Directed sampling is on by default whenever a detector is set. Turning it
  /// off falls back to a plain acceptance cut: same answer, much slower for a
  /// small detector, and useful as a cross-check.
  Generator& setDirectedSampling(bool on) {
    directed_ = on; dirty_ = true; return *this;
  }

  // --- generation ----------------------------------------------------------
  Muon generate() {
    prepare();
    const double m = particle_mass(spectrum_);
    const double env = envelopeTotal();
    for (;;) {
      ++tried_;
      const Trial t = sampleTrial(*rng_);
      if (t.f > env) ++envelope_violations_;
      // The draw happens unconditionally, including for f = 0, so that the
      // random stream does not depend on how many proposals miss.
      if (rng_->flat() * env >= t.f) continue;        // rejected

      const double p = t.momentum;
      Muon mu;
      mu.position = t.position;
      mu.direction = t.direction;
      mu.momentum = p;
      mu.energy = std::sqrt(p * p + m * m);
      mu.kinetic = mu.energy - m;
      if (is_electron(spectrum_)) {
        mu.charge = (rng_->flat() < 0.5) ? 1 : -1;
        mu.pdg = (mu.charge == 1) ? -11 : 11;
      } else {
        // Spectrum-aware: PARMA supplies its own site-dependent ratio, every
        // other spectrum falls back to the sea-level table. The legacy driver
        // deliberately does not do this, because its contract is bit-exact
        // reproduction of the Fortran, which uses the table unconditionally.
        const double r = charge_ratio(spectrum_, p);
        mu.charge = (rng_->flat() < r / (1.0 + r)) ? 1 : -1;
        mu.pdg = (mu.charge == 1) ? -13 : 13;
      }
      ++accepted_;
      return mu;
    }
  }

  // --- normalisation -------------------------------------------------------
  /// Absolute rate through the surface, in s^-1.
  ///
  /// Works for every (spectrum, surface) pair from a single estimator, because
  /// the projection lives in the surface rather than in a per-geometry
  /// formula. EcoMug needs one hand-derived integrator per geometry and gives
  /// up entirely for a user-supplied flux (GetEstimatedTime returns 0).
  ///
  /// Returns a negative value when the spectrum carries no absolute
  /// normalisation (PowerLaw and Electron are shapes only), rather than
  /// silently reporting a meaningless number.
  double rate(int npoints = 200000) const {
    double r, err;
    rateAndError(r, err, npoints);
    return r;
  }

  void rateAndError(double& out_rate, double& out_error,
                    int npoints = 200000) const {
    prepare();
    const double I_p = cdf_.integrated_flux();
    if (!(I_p > 0.0)) { out_rate = -1.0; out_error = -1.0; return; }

    // The same per-trial f that drives generate(); see the class comment for
    // the two forms of the estimator and the derivation note in
    // validation/README.md.
    UCMuonLcg local(seed_ ^ 0x5eed1234);
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < npoints; ++i) {
      const double f = sampleTrial(local).f;
      sum += f;
      sum2 += f * f;
    }
    const double mean = sum / npoints;
    const double var = std::max(sum2 / npoints - mean * mean, 0.0);
    // Omega_c is already inside f in directed mode, so only the blind form
    // carries the constant solid angle.
    const double norm = surface_->area() * I_p *
                        (directedActive() ? 1.0 : solidAngle());
    out_rate = norm * mean;
    out_error = norm * std::sqrt(var / npoints);
  }

  /// Real time, in seconds, that `n` generated muons correspond to.
  double liveTime(long long n) const {
    const double r = rate();
    return r > 0.0 ? double(n) / r : -1.0;
  }

  double surfaceArea() const { return surface_ ? surface_->area() : 0.0; }
  double solidAngle() const {
    return (phi_max_ - phi_min_) * (std::cos(theta_min_) - std::cos(theta_max_));
  }
  /// Fraction of proposals accepted; a diagnostic on sampling efficiency.
  double acceptance() const {
    return tried_ > 0 ? double(accepted_) / double(tried_) : 0.0;
  }
  double envelope() const { prepare(); return envelope_; }
  const MomentumCdf& cdf() const { prepare(); return cdf_; }
  const Detector* detector() const { return detector_.get(); }
  bool directedSampling() const { return directedActive(); }
  /// Largest cone solid angle over the surface, sr. Meaningful only in
  /// directed mode, where it is the envelope's angular factor.
  double coneSolidAngleMax() const { prepare(); return omega_max_; }
  /// Proposals whose weight exceeded the envelope. Must be 0 for a clean run;
  /// a nonzero count means the sampled maximum missed the true one and those
  /// events are slightly under-weighted.
  long long envelopeViolations() const { return envelope_violations_; }

  /// Everything needed to reproduce the run, for a logfile or a methods
  /// section.
  std::string provenance() const {
    prepare();
    char buf[640];
    std::snprintf(buf, sizeof(buf),
                  "UCMuGen %d.%d.%d | spectrum=%d | surface=%s(A=%.4g cm^2) | "
                  "detector=%s(%s) | "
                  "E=[%g,%g] GeV | theta=[%g,%g] rad | seed=%d | "
                  "envelope=%.4g",
                  UCMUGEN_VERSION_MAJOR, UCMUGEN_VERSION_MINOR,
                  UCMUGEN_VERSION_PATCH, int(spectrum_),
                  surface_ ? surface_->name() : "none", surfaceArea(),
                  detector_ ? detector_->name() : "none",
                  detector_ ? (directed_ ? "directed" : "acceptance cut")
                            : "-",
                  e_min_, e_max_, theta_min_, theta_max_, seed_, envelope_);
    return std::string(buf);
  }

 private:
  /// J(p, cos t) / J(p, 1): the angular shape, with the momentum dependence
  /// already carried by the CDF.
  double weight(double p, double cos_theta) const {
    if (cos_theta <= 0.0) return 0.0;
    const double num = flux::intensity(spectrum_, p, cos_theta);
    const double den = flux::intensity(spectrum_, p, 1.0);
    return den > 0.0 ? num / den : 0.0;
  }

  bool directedActive() const { return detector_ && directed_; }

  /// One proposal: momentum, position, direction, and the figure of merit f
  /// that both the acceptance test and the rate integral are built from.
  ///
  /// RNG consumption per call is fixed in blind mode (1 momentum + surface +
  /// 2 direction), so attaching nothing leaves the existing stream untouched.
  struct Trial {
    Vec3 position, direction;
    double momentum = 0.0;
    double f = 0.0;
  };

  Trial sampleTrial(Rng& rng) const {
    Trial t;
    t.momentum = cdf_.sample(rng.flat());
    const SurfacePoint sp = surface_->sample(rng);
    t.position = sp.position;

    double omega = 1.0;
    if (directedActive()) {
      const Cone c = coneTo(sp.position);
      t.direction = sampleInCone(rng, c.axis, c.cos_alpha);
      omega = c.solid_angle;
      // The cone is not aware of the user's angular window, so directions
      // outside it are dropped here rather than by shrinking the cone, which
      // would break the containment guarantee.
      if (!inAngularWindow(t.direction)) return t;
    } else {
      t.direction = sampleDirectionWith(rng);
    }

    if (detector_ && !detector_->intersects(t.position, t.direction)) return t;

    t.f = omega * weight(t.momentum, -t.direction.z) *
          Surface::projection(sp.normal, t.direction);
    return t;
  }

  /// The cone subtended by the detector's bounding sphere, seen from `from`.
  /// Degenerates to the full sphere when the point is inside the bound, which
  /// is correct but slow: it means the detector reaches the generation surface.
  struct Cone {
    Vec3 axis{0.0, 0.0, -1.0};
    double cos_alpha = -1.0;                 // -1 => whole sphere
    double solid_angle = 4.0 * kPi;
  };

  Cone coneTo(const Vec3& from) const {
    Cone c;
    const Vec3 to = detector_->boundCentre() - from;
    const double dist = to.norm();
    const double rad = detector_->boundRadius();
    if (!(dist > rad)) return c;
    c.axis = to * (1.0 / dist);
    const double s = rad / dist;
    c.cos_alpha = std::sqrt(std::max(1.0 - s * s, 0.0));
    c.solid_angle = 2.0 * kPi * (1.0 - c.cos_alpha);
    return c;
  }

  static Vec3 sampleInCone(Rng& rng, const Vec3& axis, double cos_alpha) {
    const double c = cos_alpha + rng.flat() * (1.0 - cos_alpha);
    const double s = std::sqrt(std::max(1.0 - c * c, 0.0));
    const double phi = 2.0 * kPi * rng.flat();
    const Vec3 u = any_perpendicular(axis);
    const Vec3 v = axis.cross(u).unit();
    return (axis * c + u * (s * std::cos(phi)) + v * (s * std::sin(phi))).unit();
  }

  /// Whether a direction falls inside the configured theta and phi window.
  /// Only needed in directed mode; blind sampling cannot leave the window.
  bool inAngularWindow(const Vec3& d) const {
    constexpr double kTol = 1e-12;
    const double c = -d.z;                                  // sky zenith cosine
    if (c < std::cos(theta_max_) - kTol) return false;
    if (c > std::cos(theta_min_) + kTol) return false;
    if (phi_max_ - phi_min_ >= 2.0 * kPi - kTol) return true;
    double phi = std::atan2(d.y, d.x);
    while (phi < phi_min_) phi += 2.0 * kPi;
    while (phi >= phi_min_ + 2.0 * kPi) phi -= 2.0 * kPi;
    return phi <= phi_max_;
  }

  Vec3 sampleDirection() const { return sampleDirectionWith(*rng_); }

  Vec3 sampleDirectionWith(Rng& rng) const {
    const double c = std::cos(theta_max_)
        + rng.flat() * (std::cos(theta_min_) - std::cos(theta_max_));
    const double st = std::sqrt(std::max(1.0 - c * c, 0.0));
    const double phi = phi_min_ + rng.flat() * (phi_max_ - phi_min_);
    return {st * std::cos(phi), st * std::sin(phi), -c};   // downward
  }

  void prepare() const {
    if (!dirty_) return;
    if (!surface_)
      throw std::runtime_error("ucmugen::Generator: no surface set");
    const double m = particle_mass(spectrum_);
    const double p_min = std::sqrt(std::max(e_min_ * e_min_ - m * m, 0.0));
    const double p_max = std::sqrt(std::max(e_max_ * e_max_ - m * m, 0.0));
    cdf_.build(p_min, p_max, spectrum_);
    computeEnvelope(p_min, p_max);
    computeOmegaMax();
    dirty_ = false;
  }

  /// Envelope for the full figure of merit. In directed mode f carries the
  /// per-trial cone solid angle, so the bound has to carry its maximum.
  double envelopeTotal() const { prepare(); return envelope_ * omega_max_; }

  /// Largest cone solid angle over the generation surface, i.e. the closest
  /// approach to the detector.
  ///
  /// Measured by sampling rather than derived, because Surface is an open
  /// interface and a user-supplied shape has no analytic closest point. The
  /// result is clamped at 4*pi, which is a true bound, so the safety margin can
  /// never push it past what is physically possible. If a proposal still
  /// exceeds the envelope it is always accepted, which mildly under-weights it
  /// relative to its peers rather than truncating it; `envelopeViolations()`
  /// reports how often that happened so the effect is visible instead of
  /// silent.
  void computeOmegaMax() const {
    omega_max_ = 1.0;
    if (!directedActive()) return;
    constexpr int kN = 20000;
    constexpr double kMargin = 1.15;
    UCMuonLcg probe(seed_ ^ 0x0c0ffee);
    double best = 0.0;
    for (int i = 0; i < kN; ++i) {
      const double w = coneTo(surface_->sample(probe).position).solid_angle;
      if (w > best) best = w;
    }
    omega_max_ = std::min(best * kMargin, 4.0 * kPi);
  }

  /// Maximum of the acceptance weight over the sampled domain.
  ///
  /// Scanned on a grid rather than assumed, because J(p,c)/J(p,1) exceeds 1
  /// near the horizon at high momentum. The margin covers the grid being
  /// finite; the acceptance test would silently truncate the distribution if
  /// the envelope were ever too small.
  void computeEnvelope(double p_min, double p_max) const {
    constexpr int kNp = 128, kNc = 128;
    constexpr double kMargin = 1.15;
    const double c_lo = std::cos(theta_max_), c_hi = std::cos(theta_min_);
    double best = 0.0;
    const double lr = (p_max > p_min) ? std::log(p_max / p_min) : 0.0;
    for (int i = 0; i < kNp; ++i) {
      const double p = (lr > 0.0)
          ? p_min * std::exp(double(i) / double(kNp - 1) * lr)
          : p_min;
      for (int j = 0; j < kNc; ++j) {
        const double c = c_lo + double(j) / double(kNc - 1) * (c_hi - c_lo);
        const double w = weight(p, c);          // projection <= 1 separately
        if (w > best) best = w;
      }
      if (lr == 0.0) break;
    }
    envelope_ = (best > 0.0) ? best * kMargin : 1.0;
  }

  Spectrum spectrum_ = Spectrum::Guan;
  std::shared_ptr<Surface> surface_;
  std::shared_ptr<Detector> detector_;
  std::unique_ptr<Rng> rng_;
  std::int32_t seed_ = 20260729;
  bool directed_ = true;

  double e_min_ = 1.0, e_max_ = 1000.0;
  double theta_min_ = 0.0, theta_max_ = kPi / 2.0;
  double phi_min_ = 0.0, phi_max_ = 2.0 * kPi;

  mutable bool dirty_ = true;
  mutable MomentumCdf cdf_;
  mutable double envelope_ = 1.0;
  mutable double omega_max_ = 1.0;
  long long tried_ = 0, accepted_ = 0, envelope_violations_ = 0;
};

// =============================================================================
// 12. Geant4 integration
//
// See the block below, which sits outside the main include guard on purpose.
// =============================================================================

}  // namespace ucmugen

#endif  // UCMUGEN_H

// =============================================================================
// Geant4 helpers
//
// Header-only and Geant4-free unless UCMUGEN_WITH_GEANT4 is defined, so the
// core still compiles without a Geant4 installation.
//
// This block deliberately lives *outside* UCMUGEN_H and carries its own guard.
// Otherwise the definition order matters: a header that includes UCMuGen.h
// before the translation unit defines UCMUGEN_WITH_GEANT4 would set UCMUGEN_H,
// the later include would expand to nothing, and FireG4 would silently not
// exist. The error that produces ("use of undeclared identifier 'FireG4'")
// points at the call site and says nothing about include order, so it is worth
// making structurally impossible rather than documenting.
//
// It exists at all because the unit and convention conversions are where
// hand-rolled integrations go wrong: cm to mm, GeV to MeV, total energy to
// kinetic energy, and the PDG sign convention (mu+ is -13, not +13). Returning
// a Cartesian direction rather than a (theta, phi) pair removes a further class
// of error: EcoMug hands back a theta already flipped by pi, so every user
// re-derives the same three lines of trigonometry and some get the sign wrong.
// =============================================================================
#if defined(UCMUGEN_WITH_GEANT4) && !defined(UCMUGEN_GEANT4_H)
#define UCMUGEN_GEANT4_H

// Pulled in here rather than assumed, so this block does not depend on the
// order in which the including file happens to have arranged its own includes.
// Defining UCMUGEN_WITH_GEANT4 *before* including this header would otherwise
// compile FireG4 with none of these types declared, which is a confusing error
// a long way from its cause. Requiring only that the macro be defined is a
// simpler contract than requiring a particular include order.
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

namespace ucmugen {

/// Fill a G4ParticleGun from a generated muon and fire it.
///
///     void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
///       ucmugen::FireG4(fGen, fGun, evt);
///     }
inline void FireG4(Generator& gen, G4ParticleGun* gun, G4Event* event) {
  const Muon m = gen.generate();
  gun->SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle(m.pdg));
  gun->SetParticlePosition(G4ThreeVector(m.position.x * CLHEP::cm,
                                         m.position.y * CLHEP::cm,
                                         m.position.z * CLHEP::cm));
  gun->SetParticleMomentumDirection(
      G4ThreeVector(m.direction.x, m.direction.y, m.direction.z));
  gun->SetParticleEnergy(m.kinetic * CLHEP::GeV);
  gun->GeneratePrimaryVertex(event);
}

}  // namespace ucmugen

#endif  // UCMUGEN_WITH_GEANT4 && !UCMUGEN_GEANT4_H
