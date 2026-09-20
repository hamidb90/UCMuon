#!/usr/bin/env python3
"""
analyze_output.py — Post-process muon_output.csv from MuonRock v3.

Usage:
    python3 analyze_output.py muon_output.csv [total_events]

    total_events = number you passed to ./MuonRock (for transmission fraction).
                   If omitted, the count at the shallowest depth is used.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

csv_file     = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("muon_output.csv")
total_events = int(sys.argv[2]) if len(sys.argv) > 2 else None

df = pd.read_csv(csv_file, comment="#")
print(f"Loaded {len(df)} muon-depth records from {csv_file}")

depths_cm = sorted(df["Depth_cm"].unique())
if total_events is None:
    total_events = df[df["Depth_cm"] == depths_cm[0]].shape[0]
    print(f"  (No total_events given — using n_transmitted at shallowest depth: {total_events})")

# ── Per-depth statistics ───────────────────────────────────────────────────
rows = []
for d in depths_cm:
    sub = df[df["Depth_cm"] == d]
    n = len(sub)
    rows.append({
        "depth_cm"      : d,
        "mwe"           : round(d * 2.65 / 100.0, 4),
        "n_transmitted" : n,
        "transmission"  : round(n / total_events, 6),
        "mean_initKE"   : round(sub["InitKE_GeV"].mean(), 3),
        "mean_exitKE"   : round(sub["ExitKE_GeV"].mean(), 3),
        "mean_eLoss"    : round(sub["EnergyLoss_GeV"].mean(), 3),
        "std_eLoss"     : round(sub["EnergyLoss_GeV"].std(), 3),
        "mean_angle"    : round(sub["AngleDeg"].mean(), 4),
        "std_angle"     : round(sub["AngleDeg"].std(), 4),
        "p95_angle"     : round(sub["AngleDeg"].quantile(0.95), 4),
    })

stats = pd.DataFrame(rows)
print("\n── Per-depth summary ─────────────────────────────────────────────")
print(stats.to_string(index=False, float_format="{:.4f}".format))
stats.to_csv("muon_depth_stats.csv", index=False)
print("\nSaved: muon_depth_stats.csv")

# ── 4-panel figure ─────────────────────────────────────────────────────────
colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(depths_cm)))
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("Geant4 v3 — Muon Transport through Standard Rock\n"
             "FTFP_BERT + G4EmStandardPhysics_option4", fontsize=13, fontweight="bold")

# Panel 1: Survival curve
ax = axes[0, 0]
ax.semilogy(stats["depth_cm"], stats["transmission"],
            "o-", color="#1a6fa8", lw=2, ms=7)
ax.set_xlabel("Depth from entry face (cm)")
ax.set_ylabel("Transmission probability")
ax.set_title("Muon Survival Fraction")
ax.grid(True, alpha=0.3)
ax2 = ax.twiny()
ax2.set_xlim(np.array(ax.get_xlim()) * 2.65 / 100.0)
ax2.set_xlabel("m.w.e.", color="gray", fontsize=10)

# Panel 2: Mean energy loss
ax = axes[0, 1]
ax.errorbar(stats["depth_cm"], stats["mean_eLoss"], yerr=stats["std_eLoss"],
            fmt="s-", color="#c0392b", lw=2, ms=7, capsize=4, label="Mean ± std")
ax.set_xlabel("Depth (cm)")
ax.set_ylabel("Energy loss (GeV)")
ax.set_title("Cumulative Energy Loss")
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 3: Angular deflection
ax = axes[1, 0]
ax.errorbar(stats["depth_cm"], stats["mean_angle"], yerr=stats["std_angle"],
            fmt="^-", color="#27ae60", lw=2, ms=7, capsize=4, label="Mean ± std")
ax.plot(stats["depth_cm"], stats["p95_angle"],
        "^--", color="#27ae60", alpha=0.5, lw=1.5, label="95th percentile")
ax.set_xlabel("Depth (cm)")
ax.set_ylabel("Deflection angle (deg)")
ax.set_title("Multiple Scattering Angle")
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 4: Angle distributions per depth
ax = axes[1, 1]
for d, c in zip(depths_cm, colors):
    sub = df[df["Depth_cm"] == d]
    ax.hist(sub["AngleDeg"], bins=40, alpha=0.55, color=c,
            label=f"{int(d)} cm", density=True)
ax.set_xlabel("Deflection angle (deg)")
ax.set_ylabel("Probability density")
ax.set_title("Angle Distributions by Depth")
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("muon_depth_analysis.png", dpi=150, bbox_inches="tight")
print("Saved: muon_depth_analysis.png")

# ── Energy loss histograms per depth ──────────────────────────────────────
n = len(depths_cm)
fig2, axes2 = plt.subplots(1, n, figsize=(3.8*n, 4), sharey=False)
if n == 1: axes2 = [axes2]
for ax, d, c in zip(axes2, depths_cm, colors):
    sub = df[df["Depth_cm"] == d]
    ax.hist(sub["EnergyLoss_GeV"], bins=40, color=c, edgecolor="white", lw=0.3)
    ax.axvline(sub["EnergyLoss_GeV"].mean(), color="black", ls="--", lw=1.5,
               label=f"Mean: {sub['EnergyLoss_GeV'].mean():.1f} GeV")
    ax.set_title(f"{int(d)} cm ({d*2.65/100:.2f} m.w.e.)", fontsize=10)
    ax.set_xlabel("Energy loss (GeV)")
    ax.set_ylabel("Counts")
    ax.legend(fontsize=8)
    ax.set_yscale("log")
fig2.suptitle("Energy Loss Distributions at Each Scoring Depth",
              fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("muon_eloss_by_depth.png", dpi=150, bbox_inches="tight")
print("Saved: muon_eloss_by_depth.png")
