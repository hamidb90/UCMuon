#!/usr/bin/env python3
"""
ucmuon_pumas_driver.py — UCMuon PUMAS transport driver
UCLouvain Muography Group | Hamid Basiri <hamid.basiri@uclouvain.be>

Wraps the ucmuon_transport_pumas C binary. Supports two modes:

  forward  — reads muons_surface.dat, transports each muon through a flat rock
             slab via PUMAS, writes muons_underground.dat (18-col, same format
             as BB/MUSIC/PROPOSAL engines).

  backward — backward MC flux integration starting from the detector. No input
             muon file needed. Samples (E_det, theta, charge) at the detector
             and transports backward via PUMAS to the generation surface.
             Writes a per-event file and returns a binned flux spectrum.

Stdin protocol (called by ucmuon_gui.py via start_run):
  mode           forward | backward
  infile         [forward only] path to surface muon file
  outfile        path for output file
  depth_m
  mat_id         1=StandardRock  2=Water  3=Seawater  4=Ice  5=Custom
  rho_gcm3       density override (0 = material default)
  energy_loss    0=CSDA  1=MIXED  2=STRAGGLED
  scattering     0=disabled  1=mixed
  transport_all  [forward only] 0 or 1
  E_min_GeV      [backward only]
  E_max_GeV      [backward only]
  theta_max_deg  [backward only]
  n_events       [backward only]
  spectrum_id    [backward only] 0=GCCLY  1=Gaisser
  seed           [backward only] 0=time-based
"""
from __future__ import annotations
import math
import os
import subprocess
import sys
import time
import numpy as np

# ── Material database ────────────────────────────────────────────────────────
# mat_id → (pumas_material_name, default_rho_gcm3)
# Water is used for Ice/Seawater with a density override (same composition).
_MAT_DB = {
    1: ("StandardRock", 2.65),
    2: ("Water",        1.00),
    3: ("Water",        1.025),   # Seawater
    4: ("Water",        0.917),   # Ice
    5: ("StandardRock", 2.65),    # Custom: name stays StandardRock, rho from user
}

# ── Path helpers ─────────────────────────────────────────────────────────────
_HERE        = os.path.dirname(os.path.abspath(__file__))
_PROJECT     = os.path.dirname(_HERE)
_BIN         = os.path.join(_PROJECT, "bin")
_PUMAS_SRC   = os.path.join(_PROJECT, "external", "pumas-master")
_MDF_PATH    = os.path.join(_PUMAS_SRC, "examples", "data", "materials.xml")
_BINARY      = os.path.join(_BIN, "ucmuon_transport_pumas")


def _dump_path(mat_name: str) -> str:
    """Physics dump cache per material — saved once in bin/."""
    return os.path.join(_BIN, f"pumas_{mat_name}.pumas")


def _detect_depth_axis(infile: str) -> int:
    """Return 0/1/2 for which of X/Y/Z has the smallest std (= depth axis)."""
    xs, ys, zs = [], [], []
    try:
        with open(infile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                xs.append(float(parts[1]))
                ys.append(float(parts[2]))
                zs.append(float(parts[3]))
                if len(xs) >= 2000:
                    break
    except Exception:
        return 2
    if not xs:
        return 2
    stds = [np.std(xs), np.std(ys), np.std(zs)]
    axis = int(np.argmin(stds))
    names = {0: "YZ (depth=X)", 1: "XZ (depth=Y)", 2: "XY (depth=Z)"}
    print(f"  Source plane: {names[axis]}  σ(x,y,z)="
          f"{stds[0]:.1f} {stds[1]:.1f} {stds[2]:.1f} cm", flush=True)
    return axis


def _check_binary() -> bool:
    if os.path.exists(_BINARY):
        return True
    print(f"  [warn] Binary not found: {_BINARY}", flush=True)
    print("  Run: make pumas   (or make all) to build it.", flush=True)
    return False


# ── Binned flux spectrum from per-event file ─────────────────────────────────

def compute_flux_spectrum(event_file: str, n_bins: int = 50,
                          cos_theta_min: float = 1.0,
                          cos_theta_max: float = 0.0) -> dict:
    """
    Read the per-event backward MC output file and compute a binned flux
    spectrum  dΦ/dE_det  [m-2 s-1 GeV-1 sr-1]  as a function of E_det.

    Returns a dict with keys:
      E_det_GeV, flux, flux_err, E_surf_mean, cos_theta_mean,
      rate_m2s, rate_err, n_events, n_reached
    """
    evs = []
    try:
        with open(event_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                evs.append([float(p) for p in parts])
    except Exception as exc:
        print(f"  [warn] Cannot read event file: {exc}", flush=True)
        return {}

    if not evs:
        return {}

    arr       = np.array(evs)
    # columns: ev  E_det  cos_theta  charge  E_surf  flux_contribution
    E_det     = arr[:, 1]
    cos_theta = arr[:, 2]
    # charge  = arr[:, 3]
    E_surf    = arr[:, 4]
    flux_val  = arr[:, 5]

    n_reached = len(evs)

    # Solid angle weight for rate integration
    w_angle = 2. * math.pi * abs(cos_theta_min - cos_theta_max)

    # We need to know total N (including events that didn't reach surface)
    # Unfortunately we only see the survived events here.
    # The driver prints "Survived: N" — the caller should pass n_total.
    # For now use n_reached as a lower bound; rate will be a lower estimate
    # if many muons stop before surface. In practice for backward MC every
    # muon should reach surface (that's the point of backward MC).
    n_total = n_reached  # correct if all events reach surface

    # Binned spectrum: log-uniform bins in E_det
    E_min = E_det.min()
    E_max = E_det.max()
    if E_max <= E_min:
        return {}
    log_edges = np.linspace(np.log(E_min), np.log(E_max), n_bins + 1)
    edges = np.exp(log_edges)
    E_centers = np.exp(0.5 * (log_edges[:-1] + log_edges[1:]))
    dE        = np.diff(edges)

    flux_binned    = np.zeros(n_bins)
    flux_err_binned = np.zeros(n_bins)
    E_surf_mean    = np.zeros(n_bins)
    cos_mean       = np.zeros(n_bins)
    counts         = np.zeros(n_bins, dtype=int)

    for b in range(n_bins):
        mask = (E_det >= edges[b]) & (E_det < edges[b + 1])
        n_b  = mask.sum()
        if n_b == 0:
            continue
        counts[b]      = n_b
        fv             = flux_val[mask]
        flux_binned[b] = fv.mean()          # dΦ/dE per sr [m-2 s-1 GeV-1 sr-1]
        flux_err_binned[b] = fv.std() / math.sqrt(n_b) if n_b > 1 else 0.
        E_surf_mean[b]     = E_surf[mask].mean()
        cos_mean[b]        = cos_theta[mask].mean()

    # Total integrated rate [m-2 s-1]:
    #   Rate = w_angle × (1/N) × Σ flux_val_i × dE_i (element-wise by bin)
    # Use event-level sum for accuracy:
    rate     = w_angle * flux_val.sum() / n_total
    var_rate = w_angle**2 * ((flux_val**2).mean() - flux_val.mean()**2) / n_total
    rate_err = math.sqrt(max(0., var_rate))

    return {
        "E_det_GeV":    E_centers,
        "flux":         flux_binned,           # dΦ/dE_det/dΩ [m-2 s-1 GeV-1 sr-1]
        "flux_err":     flux_err_binned,
        "E_surf_mean":  E_surf_mean,
        "cos_theta_mean": cos_mean,
        "counts":       counts,
        "rate_m2s":     rate,
        "rate_err":     rate_err,
        "n_events":     n_total,
        "n_reached":    n_reached,
    }


# ── Main entry point ─────────────────────────────────────────────────────────

def main():
    raw   = [ln.strip() for ln in sys.stdin.read().splitlines() if ln.strip()]
    it    = iter(raw)

    def nxt(default=""):
        try:    return next(it)
        except StopIteration: return default

    mode   = nxt("backward")          # "forward" or "backward"
    infile = nxt("")                   # surface file (forward only)
    outfile = nxt("pumas_output.dat")
    depth_m = float(nxt("100.0"))
    mat_id  = int(nxt("1"))
    rho_gcm3_user = float(nxt("0.0"))
    energy_loss   = int(nxt("0"))
    scattering    = int(nxt("0"))
    transport_all = int(nxt("0"))

    # Backward-only params
    E_min_GeV    = float(nxt("1.0"))
    E_max_GeV    = float(nxt("1000.0"))
    theta_max_deg = float(nxt("85.0"))
    n_events     = int(nxt("50000"))
    spectrum_id  = int(nxt("0"))
    seed         = int(nxt("0"))

    # Resolve material
    mat_name, rho_default = _MAT_DB.get(mat_id, ("StandardRock", 2.65))
    rho_gcm3 = rho_gcm3_user if rho_gcm3_user > 0. else rho_default

    print(f"PUMAS driver: mode={mode}  mat={mat_name}  rho={rho_gcm3:.3f} g/cm3"
          f"  depth={depth_m:.1f} m", flush=True)

    if not _check_binary():
        sys.exit(1)
    if not os.path.exists(_MDF_PATH):
        print(f"ERROR: materials XML not found: {_MDF_PATH}", file=sys.stderr)
        print("  Expected at: external/pumas-master/examples/data/materials.xml",
              file=sys.stderr, flush=True)
        sys.exit(1)

    # Derived geometry params for backward mode
    cos_theta_max_val = math.cos(math.radians(theta_max_deg))  # e.g. cos(85°)≈0.087
    cos_theta_min_val = 1.0  # vertical (theta=0)

    # Build stdin for the C binary
    if mode == "forward":
        if not infile or not os.path.exists(infile):
            print(f"ERROR: surface file not found: '{infile}'",
                  file=sys.stderr, flush=True)
            sys.exit(1)
        depth_axis = _detect_depth_axis(infile)

        stdin_lines = [
            "0",                         # mode=forward
            _MDF_PATH,
            _dump_path(mat_name),
            mat_name,
            f"{rho_gcm3:.6f}",
            f"{depth_m:.6f}",
            str(energy_loss),
            str(scattering),
            infile,
            outfile,
            str(transport_all),
            str(depth_axis),
        ]
        n_total = _count_muons(infile, transport_all)
        print(f"  Forward: {n_total} muons  axis={depth_axis}", flush=True)

    else:
        # Determine output file for per-event data
        bwd_event_file = os.path.splitext(outfile)[0] + "_bwd_events.dat"

        stdin_lines = [
            "1",                             # mode=backward
            _MDF_PATH,
            _dump_path(mat_name),
            mat_name,
            f"{rho_gcm3:.6f}",
            f"{depth_m:.6f}",
            str(energy_loss),
            str(scattering),
            bwd_event_file,
            f"{E_min_GeV:.6f}",
            f"{E_max_GeV:.6f}",
            f"{cos_theta_min_val:.6f}",
            f"{cos_theta_max_val:.6f}",
            str(n_events),
            str(spectrum_id),
            str(seed),
        ]
        n_total = n_events
        print(f"  Backward: N={n_events}  E=[{E_min_GeV}, {E_max_GeV}] GeV"
              f"  theta_max={theta_max_deg}°", flush=True)

    stdin_str = "\n".join(stdin_lines) + "\n"
    t_start   = time.perf_counter()

    # Spawn C binary
    try:
        proc = subprocess.Popen(
            [_BINARY],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except OSError as exc:
        print(f"ERROR: cannot launch {_BINARY}: {exc}", file=sys.stderr, flush=True)
        sys.exit(1)

    stdout_data, _ = proc.communicate(input=stdin_str)
    elapsed = time.perf_counter() - t_start

    # Echo C binary output to our stdout (GUI reads this for progress)
    for line in stdout_data.splitlines():
        print(line, flush=True)

    if proc.returncode != 0:
        print(f"ERROR: binary exited with code {proc.returncode}",
              file=sys.stderr, flush=True)
        sys.exit(proc.returncode)

    print(f"  Elapsed: {elapsed:.1f} s", flush=True)

    # Post-process backward MC: build flux spectrum file alongside outfile
    if mode == "backward":
        bwd_event_file = os.path.splitext(outfile)[0] + "_bwd_events.dat"
        if os.path.exists(bwd_event_file):
            spec = compute_flux_spectrum(
                bwd_event_file,
                n_bins=50,
                cos_theta_min=cos_theta_min_val,
                cos_theta_max=cos_theta_max_val,
            )
            if spec:
                _write_flux_spectrum(outfile, spec, depth_m, mat_name,
                                     rho_gcm3, theta_max_deg, n_events,
                                     spectrum_id)
                print(f"  Rate:   {spec['rate_m2s']:.4e} +/- "
                      f"{spec['rate_err']:.4e} m-2 s-1", flush=True)
                print(f"  Output: {outfile}", flush=True)


def _count_muons(infile: str, transport_all: int) -> int:
    n = 0
    try:
        with open(infile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 12:
                    continue
                if transport_all:
                    n += 1
                else:
                    hit = int(parts[12]) if len(parts) >= 14 else 1
                    if hit == 1:
                        n += 1
    except Exception:
        pass
    return n


def _write_flux_spectrum(outfile: str, spec: dict, depth_m: float,
                         mat_name: str, rho: float, theta_max: float,
                         n_events: int, spectrum_id: int):
    """Write the binned flux spectrum to outfile (5-col ASCII)."""
    with open(outfile, "w") as fh:
        fh.write(
            f"# PUMAS backward MC flux spectrum\n"
            f"# depth={depth_m:.2f} m  mat={mat_name}  rho={rho:.3f} g/cm3\n"
            f"# theta_max={theta_max:.1f} deg  N={n_events}"
            f"  spectrum={'GCCLY' if spectrum_id == 0 else 'Gaisser'}\n"
            f"# Rate: {spec['rate_m2s']:.4e} +/- {spec['rate_err']:.4e} m-2 s-1\n"
            f"# E_det_GeV  flux[m-2 s-1 GeV-1 sr-1]  flux_err"
            f"  E_surf_mean_GeV  n_events_in_bin\n"
        )
        E   = spec["E_det_GeV"]
        fl  = spec["flux"]
        err = spec["flux_err"]
        Es  = spec["E_surf_mean"]
        cnt = spec["counts"]
        for i in range(len(E)):
            if cnt[i] > 0:
                fh.write(f"{E[i]:13.6f} {fl[i]:15.6e} {err[i]:15.6e}"
                         f" {Es[i]:13.6f} {cnt[i]:8d}\n")


if __name__ == "__main__":
    main()
