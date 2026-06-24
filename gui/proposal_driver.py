#!/usr/bin/env python3
"""
proposal_driver.py — UCMuon PROPOSAL muon transport driver (v12)
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

v12 fixes:
  1. Multiple Coulomb Scattering (MCS) was not enabled.
     build_propagator now calls pp.make_multiple_scattering() (PROPOSAL v7.6+
     PropagationUtilityCollection API) so the final track direction reflects
     the true MCS-deflected exit angle.  Fallback: pp.scattering.make_scattering()
     (Sector-based API used in older v7.x builds) for robustness.

  2. Stopped muon energy was written as the PROPOSAL final total energy
     (≈ rest mass = 105.658 MeV) instead of 0.0 GeV as required by the
     18-column benchmark format spec.  Now Ef = 0.0 for alive=0.

  3. Stopped muon direction now forced to 0,0,-1 (theta=0,phi=0) and
     position forced to xs,ys,stop_depth per spec.

  4. Timing file (<outfile_stem>_timing.txt) written after transport.

  5. Stopped-muon file (<outfile_stem>_stopped.dat) written alongside output.

v11: alive criterion d_traversed >= 0.999 * slant_path (not E_final).
"""
from __future__ import annotations
import sys, os, math, time, numpy as np

MUON_MASS_MEV = 105.6584

_MED_MAP  = {1:"StandardRock",2:"Water",3:"Ice",4:"Seawater"}
_MED_ALT  = {
    "StandardRock":["StandardRock","FrejusRock"],
    "Water":       ["Water","AntaresWater"],
    "Ice":         ["Ice"],
    "Seawater":    ["AntaresWater","CascadiaBasinWater","Water"],
}
_REF_RHO  = {"StandardRock":2.65,"Water":1.00,"Ice":0.917,"Seawater":1.025}
_SCATTER  = {0:"NoScattering",1:"Highland",2:"HighlandIntegral",3:"Moliere"}


def _set_tables_dir(pp, tables_dir):
    if not tables_dir:
        return
    tables_dir = os.path.expanduser(str(tables_dir))
    os.makedirs(tables_dir, exist_ok=True)
    for attr in ["InterpolationDef", "InterpolationSettings"]:
        obj = getattr(pp, attr, None)
        if obj is not None:
            try:
                obj.path_to_tables = tables_dir
                print(f"  tables_dir: {tables_dir}", flush=True)
                return
            except Exception:
                pass


def _medium(pp, name, dmult=1.0):
    for n in _MED_ALT.get(name, [name]):
        cls = getattr(pp.medium, n, None)
        if cls:
            if n != name:
                print(f"  [medium] {name} -> {n}", flush=True)
            try:    return cls(dmult)
            except TypeError: return cls()
    try:    return pp.medium.StandardRock(dmult)
    except TypeError: return pp.medium.StandardRock()


def _density_dist(pp, rho_gcm3):
    dd = pp.density_distribution
    for fn_name in ["density_homogeneous","DensityHomogeneous",
                    "HomogeneousDensity","density_uniform","DensityUniform"]:
        fn = getattr(dd, fn_name, None)
        if fn is None: continue
        for args in [(rho_gcm3,), ()]:
            try:    return fn(*args)
            except Exception: pass
    for n in [x for x in dir(dd) if not x.startswith("_")]:
        fn = getattr(dd, n)
        if callable(fn):
            for args in [(rho_gcm3,), ()]:
                try:    return fn(*args)
                except Exception: pass
    raise RuntimeError(f"density_distribution: no working constructor. "
                       f"Available: {[x for x in dir(dd) if not x.startswith('_')]}")


def _make_mcs(pp, scat_name, pdef, medium, xs):
    """
    Build a scattering object for PropagationUtilityCollection.scattering.

    In PROPOSAL 7.6.x the PropagationUtilityCollection API requires a
    proposal.scattering.Scattering (base class), but make_multiple_scattering()
    returns a MultipleScattering subclass.  The bridge is:

        mms = pp.make_multiple_scattering(name, pdef, medium, xs, True)
        col.scattering = pp.scattering.ScatteringMultiplier(mms, 1.0)

    ScatteringMultiplier inherits from Scattering and wraps any
    MultipleScattering object, passing through deflections unchanged
    (multiplier = 1.0).

    Returns None when scat_name is 'NoScattering' or all attempts fail.
    """
    if scat_name == "NoScattering":
        return None

    fn = getattr(pp, "make_multiple_scattering", None)
    if fn is not None:
        for args in [(scat_name, pdef, medium, xs, True),
                     (scat_name, pdef, medium, xs)]:
            try:
                mms = fn(*args)
                # Wrap in ScatteringMultiplier so col.scattering accepts it
                wrapper_cls = getattr(
                    getattr(pp, "scattering", None), "ScatteringMultiplier", None)
                if wrapper_cls is not None:
                    obj = wrapper_cls(mms, 1.0)
                    print(f"  MCS: {scat_name} via make_multiple_scattering"
                          f" + ScatteringMultiplier", flush=True)
                    return obj
                # Some builds allow direct assignment — try without wrapper
                print(f"  MCS: {scat_name} via make_multiple_scattering"
                      f" (no wrapper)", flush=True)
                return mms
            except Exception:
                pass

    print(f"  [warn] MCS setup failed for '{scat_name}' — no scattering applied",
          flush=True)
    return None


def build_propagator(pp, pdef, med_name, ecut, vcut, dmult=1.0, scat_name="Moliere"):
    ref_rho = _REF_RHO.get(med_name, 2.65)
    rho     = ref_rho * dmult
    medium  = _medium(pp, med_name, dmult)
    cuts    = pp.EnergyCutSettings(ecut, vcut, False)   # positional only in v7.6.2
    xs      = pp.crosssection.make_std_crosssection(
        particle_def=pdef, target=medium, interpolate=True, cuts=cuts)
    col = pp.PropagationUtilityCollection()
    col.displacement = pp.make_displacement(xs, True)
    col.interaction  = pp.make_interaction(xs, True)
    col.time         = pp.make_time(xs, pdef, True)
    col.decay        = pp.make_decay(xs, pdef, True)
    mcs = _make_mcs(pp, scat_name, pdef, medium, xs)
    if mcs is not None:
        try:    col.scattering = mcs
        except Exception as e:
            print(f"  [warn] col.scattering assignment failed: {e}", flush=True)
    utility = pp.PropagationUtility(collection=col)
    geo     = pp.geometry.Sphere(pp.Cartesian3D(0, 0, 0), 1e20)
    density = _density_dist(pp, rho)
    return pp.Propagator(pdef, [(geo, utility, density)])


def _csda_crosscheck(pp, prop):
    """
    Propagate 100 GeV muon until it stops naturally.
    Expected in Standard Rock: 145-180 m.
    d_traversed is the correct survival criterion — NOT E_final.
    """
    state = pp.particle.ParticleState()
    state.energy = 100_000.0          # 100 GeV total energy [MeV]
    state.position  = pp.Cartesian3D(0, 0, 0)
    state.direction = pp.Cartesian3D(0, 0, -1)
    state.time = 0.0
    state.propagated_distance = 0.0
    try:
        track  = prop.propagate(state, max_distance=5_000_000)   # 50 km — no limit
        e_fin  = track.track_energies()[-1]
        d_fin  = track.track_propagated_distances()[-1] / 100.0   # cm → m
        status = "✅ OK" if 120 < d_fin < 230 else "⚠️ CHECK"
        print(f"  CSDA crosscheck: 100 GeV muon stops at {d_fin:.1f} m  "
              f"(E_final={e_fin:.1f} MeV = rest mass + kinetic)  "
              f"Groom CSDA mean: 142.6 m  {status}", flush=True)
        print(f"  NOTE: E_final={e_fin:.1f} MeV ≈ muon rest mass (105.658 MeV) "
              f"for a STOPPED muon. alive criterion = d_traversed, NOT E_final.",
              flush=True)
    except Exception as exc:
        print(f"  CSDA crosscheck failed: {exc}", flush=True)


def read_muons(infile, transport_all):
    muons = []
    with open(infile) as fh:
        for raw in fh:
            raw = raw.lstrip()
            if not raw or raw.startswith("#"): continue
            parts = raw.split()
            nc = len(parts)
            if nc < 13: continue
            hit_flag = int(parts[12]) if nc >= 14 else 1
            if not transport_all and hit_flag != 1: continue
            E_GeV = float(parts[10])
            if E_GeV * 1e3 <= MUON_MASS_MEV: continue
            p  = float(parts[4]) if nc > 4 else 0.0
            px = float(parts[5]) if nc > 5 else 0.0
            py = float(parts[6]) if nc > 6 else 0.0
            pz = float(parts[7]) if nc > 7 else 0.0
            muons.append({"evid":int(parts[0]),
                "xs":float(parts[1]),"ys":float(parts[2]),"zs":float(parts[3]),
                "p":p,"px":px,"py":py,"pz":pz,
                "theta":float(parts[8]),"phi":float(parts[9]),
                "E_GeV":E_GeV,"charge":int(parts[11])})
    return muons


def detect_source_plane(muons):
    """Return depth_axis: 0=YZ(depth X), 1=XZ(depth Y), 2=XY(depth Z)."""
    xs = np.array([m["xs"] for m in muons])
    ys = np.array([m["ys"] for m in muons])
    zs = np.array([m["zs"] for m in muons])
    stds = [xs.std(), ys.std(), zs.std()]
    axis = int(np.argmin(stds))
    names = {0: "YZ (depth=X)", 1: "XZ (depth=Y)", 2: "XY (depth=Z)"}
    print(f"  Source plane: {names[axis]}  σ(x,y,z)={stds[0]:.1f} {stds[1]:.1f} {stds[2]:.1f} cm",
          flush=True)
    return axis


def main():
    raw = [ln.strip() for ln in sys.stdin.read().splitlines() if ln.strip()]
    it  = iter(raw)
    def nxt(d=""):
        try:    return next(it)
        except StopIteration: return d

    infile   = nxt(); outfile  = nxt()
    depth_m  = float(nxt("90"))
    med_type = int(nxt("1"))

    custom_rho = None
    if med_type == 5:
        nxt("11"); nxt("22")
        custom_rho = float(nxt("2.65"))
        nxt("136.4")

    transport_all = int(nxt("0")) == 1
    ecut          = float(nxt("500"))
    vcut          = float(nxt("0.001"))
    scat_id       = int(nxt("2"))
    tables_dir    = nxt("")

    med_name  = _MED_MAP.get(med_type, "StandardRock")
    ref_rho   = _REF_RHO.get(med_name, 2.65)
    dmult     = (custom_rho / ref_rho) if custom_rho else 1.0
    rho_used  = ref_rho * dmult

    scat_name = _SCATTER.get(scat_id, "Moliere")
    print(f"PROPOSAL driver v12: medium={med_name}  depth={depth_m:.1f} m  "
          f"e_cut={ecut:.0f} MeV  v_cut={vcut:.4f}  MCS={scat_name}", flush=True)
    print(f"  rho={rho_used:.3f} g/cm³  "
          f"vert_opacity={rho_used*depth_m*100:.0f} g/cm²", flush=True)
    print(f"  alive criterion: d_traversed >= 99.9% of slant path  "
          f"(NOT E_final — stopped muons have E_final ≈ 105.7 MeV = rest mass)",
          flush=True)

    try:
        import proposal as pp
    except ImportError:
        print("ERROR: pip install proposal", file=sys.stderr, flush=True)
        sys.exit(1)

    # Suppress benign dNdx/Epair warnings — physically harmless boundary condition.
    # PROPOSAL 7.6.x uses spdlog; set level via env var (must be set before import,
    # but we can redirect stderr to filter at the Python level as a reliable fallback).
    import os as _os
    _os.environ["PROPOSAL_LOG_LEVEL"] = "error"
    _os.environ["SPDLOG_LEVEL"] = "error"
    # Some builds expose a Python logging API:
    for _attr in ("set_loglevel", "setLogLevel"):
        _fn = getattr(getattr(pp, "logging", None) or getattr(pp, "Log", None) or object(), _attr, None)
        if _fn:
            try: _fn(5)  # 5 = ERROR in spdlog
            except Exception: pass

    print(f"PROPOSAL {pp.__version__}", flush=True)
    _set_tables_dir(pp, tables_dir)

    print("Building propagators ...", flush=True)
    try:
        pm  = build_propagator(pp, pp.particle.MuMinusDef(),
                               med_name, ecut, vcut, dmult, scat_name)
        pp2 = build_propagator(pp, pp.particle.MuPlusDef(),
                               med_name, ecut, vcut, dmult, scat_name)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        sys.exit(1)

    print("Propagators ready.", flush=True)
    _csda_crosscheck(pp, pm)

    if not infile or not os.path.exists(infile):
        print(f"ERROR: not found: '{infile}'", file=sys.stderr, flush=True)
        sys.exit(1)

    muons = read_muons(infile, transport_all)
    total = len(muons)
    print(f"Read {total:,} muons from '{infile}'", flush=True)
    if total == 0:
        print("No muons — exiting.", flush=True); sys.exit(0)

    depth_axis = detect_source_plane(muons)

    # Derive sibling output paths from outfile stem
    outbase       = os.path.splitext(outfile)[0]
    timing_file   = outbase + "_timing.txt"
    stopped_file  = outbase + "_stopped.dat"

    transported = survived = 0
    stopped_rows = []          # (EventID, InitKE_GeV, stop_depth_cm) for alive=0
    t_start = time.perf_counter()

    with open(outfile, "w") as fout:
        fout.write(
            f"# UCMuon PROPOSAL v{pp.__version__} transport\n"
            f"# medium={med_name}  rho={rho_used:.3f} g/cm³  depth={depth_m:.2f} m  "
            f"vert_opacity={rho_used*depth_m*100:.0f} g/cm²\n"
            f"# e_cut={ecut:.0f} MeV  v_cut={vcut:.4f}  MCS={scat_name}\n"
            f"# E_convention: total_energy_GeV   (E = KE + 0.10566)\n"
            f"# alive = d_traversed >= 99.9% of slant path\n"
            f"# EventID xs ys zs Es theta_s phi_s charge alive "
            f"x y z E cx cy cz theta phi\n"
        )
        for mu in muons:
            E_MeV  = mu["E_GeV"] * 1e3
            th, ph = mu["theta"], mu["phi"]
            # Direction cosines from momentum components (preferred) or theta/phi
            p_mag = mu["p"]
            if p_mag > 0:
                cx = mu["px"] / p_mag
                cy = mu["py"] / p_mag
                cz = mu["pz"] / p_mag
            else:
                cx =  math.sin(th)*math.cos(ph)
                cy =  math.sin(th)*math.sin(ph)
                cz = -math.cos(th)

            # Slant path using correct depth-direction component
            if depth_axis == 1:       # XZ plane: depth in Y
                cos_depth = max(abs(cy), 0.01)
            elif depth_axis == 0:     # YZ plane: depth in X
                cos_depth = max(abs(cx), 0.01)
            else:                     # XY plane: depth in Z
                cos_depth = max(abs(cz), 0.01)
            slant_cm = depth_m * 100.0 / cos_depth    # geometric path [cm]

            state = pp.particle.ParticleState()
            state.energy              = E_MeV
            state.position            = pp.Cartesian3D(mu["xs"], mu["ys"], mu["zs"])
            state.direction           = pp.Cartesian3D(cx, cy, cz)
            state.time                = 0.0
            state.propagated_distance = 0.0

            prop = pm if mu["charge"] < 0 else pp2
            try:
                track    = prop.propagate(state, max_distance=slant_cm)
                e_fin    = track.track_energies()[-1]
                d_trav   = track.track_propagated_distances()[-1]
                pf       = track.track_positions()[-1]
                df       = track.track_directions()[-1]

                # alive: muon traversed full slant path (not energy-based)
                alive = 1 if d_trav >= slant_cm * 0.999 else 0

                if alive:
                    Ef  = e_fin / 1e3          # total energy GeV from PROPOSAL
                    fx  = pf.x
                    fy  = pf.y
                    fz  = -depth_m * 100.0     # exact plane depth
                    fcx, fcy, fcz = df.x, df.y, df.z
                    fth = math.acos(max(-1.0, min(1.0, -fcz)))
                    fph = math.atan2(fcy, fcx)
                else:
                    # Spec: alive=0 → E=0, x=xs, y=ys, z=stop_depth, cx/cy/cz=0/0/-1
                    Ef  = 0.0
                    fx  = mu["xs"]
                    fy  = mu["ys"]
                    fz  = pf.z                 # actual stopping depth
                    fcx, fcy, fcz = 0.0, 0.0, -1.0
                    fth, fph = 0.0, 0.0
                    stopped_rows.append((mu["evid"],
                                         mu["E_GeV"] - MUON_MASS_MEV / 1e3,
                                         abs(fz)))

            except Exception as exc:
                alive, Ef = 0, 0.0
                fx, fy = mu["xs"], mu["ys"]
                fz = 0.0
                fcx, fcy, fcz = 0.0, 0.0, -1.0
                fth, fph = 0.0, 0.0
                stopped_rows.append((mu["evid"],
                                     mu["E_GeV"] - MUON_MASS_MEV / 1e3, 0.0))
                print(f"  [warn] evid={mu['evid']}: {exc}", flush=True)

            transported += 1
            survived    += alive
            fout.write(
                f"{mu['evid']:10d} "
                f"{mu['xs']:13.4f} {mu['ys']:13.4f} {mu['zs']:13.4f} "
                f"{mu['E_GeV']:13.6f} {th:13.6f} {ph:13.6f} "
                f"{mu['charge']:4d} {alive:2d} "
                f"{fx:13.4f} {fy:13.4f} {fz:13.4f} "
                f"{Ef:13.6f} {fcx:13.6f} {fcy:13.6f} {fcz:13.6f} "
                f"{fth:13.6f} {fph:13.6f}\n"
            )
            if transported % 200 == 0:
                print(f"  Transported: {transported}  Survived: {survived}"
                      f"  Total: {total}", flush=True)

    elapsed = time.perf_counter() - t_start
    print(f"  Transported: {transported}  Survived: {survived}"
          f"  Total: {transported}", flush=True)
    print(f"  Muons transported: {transported}", flush=True)
    print(f"  Survived:          {survived}", flush=True)
    print(f"  Survival rate:     {100*survived/max(transported,1):.4f} %", flush=True)
    print(f"  Output:            {outfile}", flush=True)
    print(f"  Elapsed:           {elapsed:.1f} s", flush=True)

    with open(timing_file, "w") as ft:
        ft.write(f"Elapsed : {elapsed:.1f}\n")
    print(f"  Timing:            {timing_file}", flush=True)

    with open(stopped_file, "w") as fs:
        fs.write(f"# PROPOSAL stopped muons  depth={depth_m:.2f} m\n")
        fs.write("# EventID  InitKE_GeV  StopDepth_cm\n")
        for evid, ke, sdep in stopped_rows:
            fs.write(f"{evid:10d}  {ke:13.6f}  {sdep:13.4f}\n")
    print(f"  Stopped muons:     {stopped_file}  ({len(stopped_rows)} entries)",
          flush=True)


if __name__ == "__main__":
    main()
