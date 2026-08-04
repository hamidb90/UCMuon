#!/usr/bin/env python3
"""Generate ucmugen/include/UCMuGen_PARMA.h from an official PARMA C++ release.

Why generate rather than hand-port
----------------------------------
PARMA's muon path is ~800 lines of interpolation over fitted tables. Retyping it
would be a transcription exercise with no upside and a large downside, so this
script lifts the routines *verbatim* from JAEA's own `subroutines.cpp` and applies
exactly three mechanical edits, each of which is easy to audit and is listed in
the generated header:

  1. every `ifstream X(dname, ios::in)` becomes an `istringstream` over a table
     embedded in the header, falling back to the real file when a data directory
     has been configured;
  2. the angular initialiser, which loops over all six particle species, is
     restricted to muons, so only the muon tables need embedding;
  3. the black-hole-factor initialiser is restricted the same way.

Edits 2 and 3 are not taken on trust: `validation/test_parma.cc` compares the
generated header against the stock build over a dense grid and requires exact
agreement.

Licence
-------
The generated header contains JAEA code and data. The EXPACS "Conditions for
Use" grant redistribution and modification for any non-commercial purpose
(clause 2), which is what makes this legitimate, and forbid commercial use
without prior agreement. The generated header is therefore NOT MIT, which is
precisely why it is a separate file from UCMuGen.h.

Usage
-----
    python3 make_parma_header.py /path/to/parma_cpp [-o ../include/UCMuGen_PARMA.h]

where the argument is an unpacked parma_cpp.zip from
https://phits.jaea.go.jp/expacs/download-eng.htm
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# The muon path, and nothing else. Order is the order of definition in the
# source, which is also a valid definition order for the header.
NEEDED = [
    "getFFPfromWCpp",
    "getPowCpp",
    "getAmuon",
    "getMuonSpecCpp",
    "funcAngCpp",
    "getintCpp",
    "adjustParaAdepCpp",
    "getParaAdepCpp",
    "BHfactorCpp",
    "getGneutCpp",
    "getGmuonCpp",
    "getSpecAngCpp",
    "getSpecAngFinalCpp",
    "getHPcpp",
    "getrcpp",
    "getdcpp",
]

# Tables small enough to embed. The muon physics needs the first nine; the tenth
# converts altitude to atmospheric depth, which is the headline feature and only
# 48 kB. Deliberately NOT embedded: input/CORdata.inp (1.8 MB, latitude and
# longitude to cutoff rigidity) and input/FFPtable.day (233 kB, date to solar
# modulation, and revised by JAEA over time so a frozen copy would go stale).
# Those load from a user-supplied data directory instead.
EMBED = [
    "input/muon--/final135.plus",
    "input/muon--/final135.mins",
    "input/muon--/final2467.plus",
    "input/muon--/final2467.mins",
    "input/muon--/solar-dep.plus",
    "input/muon--/solar-dep.mins",
    "input/angle/muon--.out",
    "input/angle/muon---Eint.out",
    "input/angle/BkH-muon--.inp",
    "input/angle/NeutronGround.out",
    "input/AtomDepth.inp",
]


def split_functions(src: str) -> dict[str, str]:
    """Map function name -> verbatim source text, including its leading comment."""
    lines = src.split("\n")
    heads = []
    for i, line in enumerate(lines):
        m = re.match(r"^(?:double|void|int)\s+(\w+)\s*\(", line)
        if m:
            heads.append((i, m.group(1)))
    out = {}
    for k, (i, name) in enumerate(heads):
        end = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        # Walk back over the banner comment that precedes the definition.
        start = i
        while start > 0 and (lines[start - 1].startswith("//")
                             or lines[start - 1].strip() == ""):
            if lines[start - 1].strip() == "" and start - 1 < i - 1:
                break
            start -= 1
        body = "\n".join(lines[start:end]).rstrip()
        out[name] = body
    return out


def patch(text: str, name: str) -> tuple[str, list[str]]:
    """Apply the three mechanical edits. Returns the text and what was done."""
    notes = []

    # 1. File reads become reads from the embedded tables.
    def repl(m):
        return f"std::istringstream {m.group(1)}(table({m.group(2)}));"

    text, n = re.subn(r"ifstream\s+(\w+)\s*\(\s*(\w+)\s*,\s*ios::in\s*\);",
                      repl, text)
    if n:
        notes.append(f"{n} file open(s) -> embedded table lookup")

    # 2 and 3. Restrict the species initialisers to muons. Both loops fill
    # per-species arrays, but getSpecAngCpp's also fills the depth and rigidity
    # grids that the interpolation reads, and those grids are shared across
    # species rather than per-species. Leaving the loop at 1..6 with only the
    # muon tables embedded would let the five empty reads corrupt them.
    if name == "getSpecAngCpp":
        text, n = re.subn(r"for\(ip1=1;ip1<=npart;ip1\+\+\)",
                          "for(ip1=4;ip1<=4;ip1++)  /* muons only */", text)
        if n != 1:
            sys.exit("getSpecAngCpp: species loop not found; source changed?")
        notes.append("species loop restricted to muons (ip1=4)")
    if name == "BHfactorCpp":
        text, n = re.subn(r"for\(ip2=1;ip2<=npart;ip2\+\+\)",
                          "for(ip2=4;ip2<=4;ip2++)  /* muons only */", text)
        if n != 1:
            sys.exit("BHfactorCpp: species loop not found; source changed?")
        notes.append("species loop restricted to muons (ip2=4)")

    return text, notes


def embed_tables(root: Path) -> tuple[str, int]:
    """Emit the tables as raw string literals.

    Raw literals rather than parsed arrays on purpose: the bytes reach
    istringstream exactly as they would have reached ifstream, so the parse is
    the same one PARMA has always done and no float round-trip is introduced.
    """
    chunks, total = [], 0
    chunks.append("inline const std::map<std::string, std::string>& tables() {\n"
                  "  static const std::map<std::string, std::string> t = {\n")
    for rel in EMBED:
        path = root / rel
        if not path.exists():
            sys.exit(f"missing table: {path}")
        data = path.read_text(encoding="utf-8", errors="strict")
        if ')PARMA"' in data:
            sys.exit(f"{rel}: contains the raw-literal delimiter")
        total += len(data)
        chunks.append(f'    {{"{rel}", R"PARMA({data})PARMA"}},\n')
    chunks.append("  };\n  return t;\n}\n")
    return "".join(chunks), total


HEADER = '''// =============================================================================
//  UCMuGen_PARMA.h -- PARMA/EXPACS muon flux for UCMuGen. GENERATED FILE.
//
//  Regenerate with:
//      python3 ucmugen/tools/make_parma_header.py /path/to/parma_cpp
//
//  ---------------------------------------------------------------------------
//  LICENCE: THIS FILE IS NOT MIT.
//
//  It contains code and data from PARMA/EXPACS v4.10, Copyright (c) 2006 Japan
//  Atomic Energy Agency. The published "Conditions for Use"
//  (https://phits.jaea.go.jp/expacs/download-eng.htm) state:
//
//    "These program are free softwares; you can redistribute them and/or modify
//     them for any purposes except for commercial use."
//
//  so redistribution and modification are granted for NON-COMMERCIAL use, and
//  COMMERCIAL USE IS NOT PERMITTED without prior agreement with JAEA. Including
//  this header therefore makes the resulting work non-commercial-only, which is
//  why it is a separate file: UCMuGen.h on its own stays MIT.
//
//  Any published use must cite:
//    T. Sato, PLoS ONE 11(8): e0160390 (2016).
//    T. Sato, PLoS ONE 10(12): e0144679 (2015).
//  and acknowledge http://phits.jaea.go.jp/expacs/
//  Send a copy of publications to nsed-expacs@jaea.go.jp
//  ---------------------------------------------------------------------------
//
//  Usage:
//      #include "UCMuGen.h"
//      #include "UCMuGen_PARMA.h"
//
//      ucmugen::parma::Site site;              // Louvain-la-Neuve, solar min
//      site.cutoff_GV  = 3.25;                 // vertical cutoff rigidity
//      site.depth_gcm2 = ucmugen::parma::depth_from_altitude(0.0, 50.67);
//      ucmugen::parma::install(site);          // wires Spectrum::Parma
//
//      gen.setSpectrum(ucmugen::Spectrum::Parma);
//
//  Cutoff rigidity from latitude and longitude, and the solar W index from a
//  date, need the two large tables that are NOT embedded here. Point at an
//  unpacked parma_cpp (or the parma/ directory of an EXPACS install) to use
//  them:
//      ucmugen::parma::set_data_dir("/path/to/parma_cpp");
//      site.cutoff_GV = ucmugen::parma::cutoff_from_location(50.67, 4.62);
//      site.w_index   = ucmugen::parma::w_index_from_date(2026, 1, 1);
//
//  What was changed from the stock JAEA source, and nothing else:
{EDITS}
//
//  Embedded tables: {NTABLES} files, {KB} kB.
// =============================================================================
#ifndef UCMUGEN_PARMA_H
#define UCMUGEN_PARMA_H

#ifndef UCMUGEN_H
#error "include UCMuGen.h before UCMuGen_PARMA.h"
#endif

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace ucmugen {
namespace parma {
namespace detail {

// The stock subroutines.cpp opens with `using namespace std;`. Reproducing
// that wholesale in a header would be antisocial, so the names it actually
// uses are pulled in one by one.
//
// std::abs is not optional here. The source calls unqualified abs() on a
// double; without this line it resolves to C's int abs() and silently
// truncates, which changes getFFPfromW(-135.4) from 0.1988 to 1.7818. libc++
// hoists the double overload into the global namespace so clang gets it right
// by accident, while libstdc++ does not. It was caught by comparing against
// the stock build under both compilers, which is why that check is worth
// keeping.
using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::ceil;
using std::cos;
using std::cout;
using std::exp;
using std::fabs;
using std::floor;
using std::getline;
using std::istringstream;
using std::log;
using std::log10;
using std::max;
using std::min;
using std::pow;
using std::sin;
using std::sqrt;
using std::string;
using std::tan;

'''


FOOTER = '''
}  // namespace detail

// =============================================================================
// Public interface
// =============================================================================

/// Atmospheric depth in g/cm^2 from altitude, using the embedded table.
/// Pass lat = 100 to use the US Standard Atmosphere 1976 instead of the
/// latitude-dependent profile.
inline double depth_from_altitude(double alt_km, double lat_deg) {
  return detail::getdcpp(alt_km, lat_deg);
}

/// Where to find the two tables too large to embed. Only needed for
/// cutoff_from_location and w_index_from_date.
inline std::string& data_dir() { return detail::data_dir_ref(); }
inline void set_data_dir(const std::string& d) { data_dir() = d; }

/// Vertical cutoff rigidity in GV. Needs set_data_dir (input/CORdata.inp).
inline double cutoff_from_location(double lat_deg, double lon_deg) {
  if (data_dir().empty())
    throw std::runtime_error("ucmugen::parma: cutoff_from_location needs "
                             "set_data_dir(); CORdata.inp is not embedded");
  return detail::getrcpp(lat_deg, lon_deg);
}

/// Solar modulation W index for a date. Needs set_data_dir (input/FFPtable.*).
inline double w_index_from_date(int year, int month, int day) {
  if (data_dir().empty())
    throw std::runtime_error("ucmugen::parma: w_index_from_date needs "
                             "set_data_dir(); FFPtable is not embedded");
  return detail::getHPcpp(year, month, day);
}

/// A place and a moment in the solar cycle.
struct Site {
  double w_index = 0.0;        ///< solar modulation W index (0 = mean activity)
  double cutoff_GV = 3.0;      ///< vertical cutoff rigidity, GV
  double depth_gcm2 = 1033.0;  ///< atmospheric depth, g/cm^2 (sea level)
  /// Local geometry: 0..1 is the ground water fraction (0.15 is the JAEA
  /// recommendation when unknown), 10 is the ideal atmosphere with no ground,
  /// 100 is "black hole" mode with no albedo.
  double local_g = 0.15;
};

/// Differential intensity of mu+ plus mu- at the site,
/// in cm^-2 s^-1 sr^-1 (GeV/c)^-1, to match ucmugen::flux::intensity.
///
/// PARMA works in total energy per MeV and splits the two charges, so this
/// converts on both counts: it sums the charges and applies dE/dp = p/E along
/// with the MeV to GeV factor.
inline double intensity(const Site& s, double p_GeV, double cos_theta) {
  if (cos_theta <= 0.0 || p_GeV <= 0.0) return 0.0;
  const double E_GeV = std::sqrt(p_GeV * p_GeV + kMuonMass * kMuonMass);
  const double Ek_MeV = (E_GeV - kMuonMass) * 1.0e3;   // PARMA wants kinetic
  if (Ek_MeV <= 0.0) return 0.0;
  const double plus = detail::getMuonSpecCpp(1, s.w_index, s.cutoff_GV,
                                             s.depth_gcm2, Ek_MeV);
  const double mins = detail::getMuonSpecCpp(2, s.w_index, s.cutoff_GV,
                                             s.depth_gcm2, Ek_MeV);
  const double ang = detail::getSpecAngFinalCpp(4, s.w_index, s.cutoff_GV,
                                                s.depth_gcm2, Ek_MeV,
                                                s.local_g, cos_theta);
  const double per_MeV_per_sr = std::max(0.0, plus + mins) * std::max(0.0, ang);
  return per_MeV_per_sr * 1.0e3 * (p_GeV / E_GeV);
}

/// PARMA's own mu+/mu- ratio, which is a function of the site as well as the
/// momentum. UCMuGen's built-in ratio is a fit to sea-level data; this one
/// carries the altitude and geomagnetic dependence.
inline double charge_ratio(const Site& s, double p_GeV) {
  const double E_GeV = std::sqrt(p_GeV * p_GeV + kMuonMass * kMuonMass);
  const double Ek_MeV = (E_GeV - kMuonMass) * 1.0e3;
  if (Ek_MeV <= 0.0) return 1.0;
  const double plus = detail::getMuonSpecCpp(1, s.w_index, s.cutoff_GV,
                                             s.depth_gcm2, Ek_MeV);
  const double mins = detail::getMuonSpecCpp(2, s.w_index, s.cutoff_GV,
                                             s.depth_gcm2, Ek_MeV);
  return (mins > 0.0) ? plus / mins : 1.0;
}

/// Wire this site into Spectrum::Parma, so the generator, the momentum CDF and
/// rate() all use it. Call again to change site; the caller must mark the
/// generator dirty (any setter does) for the CDF to rebuild.
inline void install(const Site& s) {
  flux::parma_provider() = [s](double p_GeV, double cos_theta) {
    return intensity(s, p_GeV, cos_theta);
  };
  // PARMA models the two charges separately, so it also supplies the charge
  // ratio. Without this the generator would use a site-independent sea-level
  // fit alongside a site-aware spectrum, which is inconsistent: at 5 km the
  // fit is high by about 12% at 1 GeV/c.
  charge_ratio_provider() = [s](double p_GeV) { return charge_ratio(s, p_GeV); };
}

/// Remove both providers, so Spectrum::Parma throws again and the charge ratio
/// reverts to the built-in table.
inline void uninstall() {
  flux::parma_provider() = nullptr;
  charge_ratio_provider() = nullptr;
}

}  // namespace parma
}  // namespace ucmugen

#endif  // UCMUGEN_PARMA_H
'''


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("parma_cpp", type=Path,
                    help="unpacked parma_cpp.zip directory")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path(__file__).resolve().parents[1]
                    / "include" / "UCMuGen_PARMA.h")
    args = ap.parse_args()

    src_path = args.parma_cpp / "subroutines.cpp"
    if not src_path.exists():
        sys.exit(f"not a parma_cpp directory: {src_path} missing")
    funcs = split_functions(src_path.read_text(encoding="utf-8",
                                               errors="replace"))

    missing = [n for n in NEEDED if n not in funcs]
    if missing:
        sys.exit(f"routines not found in subroutines.cpp: {missing}")

    bodies, all_notes = [], []
    for name in NEEDED:
        text, notes = patch(funcs[name], name)
        bodies.append(text)
        for n in notes:
            all_notes.append(f"{name}: {n}")

    tables, nbytes = embed_tables(args.parma_cpp)

    edits = "\n".join(f"//    - {n}" for n in all_notes)
    header = (HEADER.replace("{EDITS}", edits)
                    .replace("{NTABLES}", str(len(EMBED)))
                    .replace("{KB}", f"{nbytes / 1024:.0f}"))

    # Forward declarations: the routines call each other in both directions.
    fwd = ["// Forward declarations: these routines are mutually recursive.",
           "inline double getFFPfromWCpp(double);",
           "inline double getPowCpp(int,double,double);",
           "inline void getAmuon(int,int,double,double,double*);",
           "inline double getMuonSpecCpp(int,double,double,double,double);",
           "inline double funcAngCpp(double,double*);",
           "inline double getintCpp(double,double,double,double,double);",
           "inline void adjustParaAdepCpp(double*);",
           "inline double getParaAdepCpp(int,double,double*);",
           "inline double BHfactorCpp(int,double,double);",
           "inline double getGneutCpp(double,int);",
           "inline double getGmuonCpp(double,double);",
           "inline double getSpecAngCpp(int,double,double,double,double,double,double);",
           "inline double getSpecAngFinalCpp(int,double,double,double,double,double,double);",
           "inline double getHPcpp(int,int,int);",
           "inline double getrcpp(double,double);",
           "inline double getdcpp(double,double);",
           ""]

    lookup = ""

    lookup_impl = '''
// Where the two tables too large to embed are found. Held here rather than in
// the public interface so `table()` below can reach it.
inline std::string& data_dir_ref() { static std::string d; return d; }

// Table lookup. Embedded tables win; anything else falls back to a real file
// under the configured data directory, so the two large tables still work when
// the user points at an EXPACS install. A miss returns empty, which parses
// exactly as a missing file did, so an unconfigured lookup degrades the way
// PARMA already degrades rather than crashing.
inline const std::string& table(const std::string& name) {
  static const std::string empty;
  const auto& t = tables();
  const auto it = t.find(name);
  if (it != t.end()) return it->second;
  static std::map<std::string, std::string> cache;
  const auto c = cache.find(name);
  if (c != cache.end()) return c->second;
  const std::string& dir = data_dir_ref();
  if (dir.empty()) return empty;
  std::ifstream f(dir + "/" + name);
  if (!f) return empty;
  std::ostringstream ss;
  ss << f.rdbuf();
  return cache.emplace(name, ss.str()).first->second;
}
'''

    # `inline` every definition: this is a header.
    text_bodies = []
    for b in bodies:
        b = re.sub(r"^(double|void|int)\s+(\w+)\s*\(", r"inline \1 \2(",
                   b, flags=re.M)
        text_bodies.append(b)

    out = (header
           + tables + "\n"
           + lookup + lookup_impl + "\n"
           + "\n".join(fwd) + "\n\n"
           + "\n\n".join(text_bodies) + "\n"
           + FOOTER)

    args.out.write_text(out, encoding="utf-8")
    print(f"wrote {args.out}  ({len(out) / 1024:.0f} kB source, "
          f"{nbytes / 1024:.0f} kB tables, {len(NEEDED)} routines)")
    for n in all_notes:
        print(f"  edit: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
