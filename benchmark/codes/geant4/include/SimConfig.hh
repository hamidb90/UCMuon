#ifndef SimConfig_hh
#define SimConfig_hh

#include <string>
#include <vector>

namespace SimConfig {

// ============================================================================
// SOURCE MODE — override at runtime:
//   ./MuonRock -f /path/to/file.txt 10000
//   ./MuonRock -m powerlaw 50000
//   ./MuonRock -m ecomug   50000
// ============================================================================
enum class SourceMode { FILE, POWERLAW, ECOMUG };

inline SourceMode  SOURCE_MODE = SourceMode::FILE;
inline std::string SOURCE_FILE = "";
inline std::string OUTPUT_DIR  = "outputs";   // set via -o flag

// ============================================================================
// FILE MODE column layout
// ============================================================================
struct CSVColumnConfig {
    bool   has_header    = false;
    bool   pdg_is_string = false;
    int    col_pdg  = 0,  col_x  = 1,  col_y  = 2,  col_z  = 3;
    int    col_px   = 4,  col_py = 5,  col_pz = 6,  col_ke = 7;
    double pos_to_mm = 1.0;
    double ke_to_MeV = 1.0;
};
inline const CSVColumnConfig CSV_COLS;

// ============================================================================
// POWERLAW MODE
// ============================================================================
inline double SPEC_EMIN_GEV              = 100.0;
inline double SPEC_EMAX_GEV              = 10000.0;
constexpr double SPEC_GAMMA                 = 2.7;
constexpr double POWERLAW_MU_MINUS_FRACTION = 0.43;
constexpr double SOURCE_Z_POWERLAW_CM       = 1.0;
inline   double POWERLAW_MAX_THETA_RAD     = 0.5236;   // 30 deg

/// ============================================================================
// ECOMUG MODE
// ============================================================================
// ── ECOMUG ────────────────────────────────────────────────────────────────────
inline double ECOMUG_PMIN_MEV          = 0.0;    // set from SPEC_EMIN_GEV at startup
inline double ECOMUG_PMAX_MEV          = 0.0;    // set from SPEC_EMAX_GEV at startup
inline   double ECOMUG_MAX_THETA_RAD     = 0.5236;      // 30 deg

// Sky = flat rectangular source plane above rock entry face
// SetSkySize receives {full_width_mm, full_height_mm}
constexpr double ECOMUG_PLANE_HALFX_MM    = 24990.0; // just inside ±25 m slab
constexpr double ECOMUG_PLANE_HALFY_MM    = 24990.0;
constexpr double ECOMUG_PLANE_Z_MM        = 1.0;     // 1 mm above rock entry z=0
// ============================================================================
// DETECTOR GEOMETRY
//   Slab: 50 m × 50 m lateral  (halfXY = 25 m = 2500 cm)
//         200 m deep            (scoring planes up to 20000 cm)
//   World: z = +505 cm (5 m headroom) → z = -20010 cm
//   Rock entry face: z = 0   (muons enter here, travel in -z)
// ============================================================================
constexpr double SLAB_HALF_XY_CM = 2500.0;           // 25 m
inline const std::vector<double> SCORING_DEPTHS_CM = {
    100.0, 1000.0, 2500.0, 5000.0, 10000.0, 20000.0  // 1m 10m 25m 50m 100m 200m
};

// ============================================================================
// ROCK MATERIAL
// ============================================================================
constexpr double ROCK_DENSITY = 2.65;
constexpr double ROCK_Z_EFF   = 11.0;
constexpr double ROCK_A_EFF   = 22.0;

// ============================================================================
// PHYSICS
// ============================================================================
constexpr double PRODUCTION_CUT_MM = 10.0;

}  // namespace SimConfig
#endif
