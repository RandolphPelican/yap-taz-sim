"""
Simulation 1: Feedback Knockout — Is the Hill loop actually necessary?
======================================================================
Compares the wild-type bistable model (Hill cooperativity n=4) against
a feedback-deficient variant (sigma(R) = 1, i.e., n→0 or n=1, theta→0).

Prediction: removing positive feedback collapses bistability into a graded
monotonic response — confirming the Hill loop is sine qua non.

Author: John D. Stabler | MIT License
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import PARAMS, steady_states, f_stiffness

S_range = np.linspace(0, 10, 300)

# ── Wild-type: full feedback ─────────────────────────────────────────────────
S_wt, R_lo_wt, R_hi_wt = steady_states(S_range, p=PARAMS)

# ── Knockout: n=1, theta→0 (linear, no ultrasensitivity) ────────────────────
p_ko = {**PARAMS, "n": 1, "theta": 0.01}
S_ko, R_lo_ko, R_hi_ko = steady_states(S_range, p=p_ko)

# ── Partial knockout: n=2 (intermediate cooperativity) ──────────────────────
p_pk = {**PARAMS, "n": 2}
S_pk, R_lo_pk, R_hi_pk = steady_states(S_range, p=p_pk)

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
fig.suptitle(
    "Simulation 1: Feedback Knockout — Necessity of Ultrasensitive Hill Loop",
    fontsize=13, fontweight="bold", y=1.02
)

configs = [
    (axes[0], S_wt,  R_lo_wt,  R_hi_wt,  "Wild-Type (n=4)\nFull Bistability",   "#2196F3"),
    (axes[1], S_pk,  R_lo_pk,  R_hi_pk,  "Partial KO (n=2)\nWeakened Bistability", "#FF9800"),
    (axes[2], S_ko,  R_lo_ko,  R_hi_ko,  "Full KO (n=1)\nMonotonic Response",   "#F44336"),
]

for ax, S, R_lo, R_hi, title, color in configs:
    ax.plot(S, R_hi, "-",  color=color, lw=2.5, label="High IC (R₀=0.99)")
    ax.plot(S, R_lo, "--", color=color, lw=2.5, label="Low IC (R₀=0.01)", alpha=0.7)

    gap = R_hi - R_lo
    bistable_mask = gap > 0.1
    if bistable_mask.any():
        ax.fill_between(S, R_lo, R_hi, where=bistable_mask,
                        alpha=0.15, color=color, label="Bistable zone")
        lo_edge = S[bistable_mask][0]
        hi_edge = S[bistable_mask][-1]
        ax.axvline(lo_edge, color=color, lw=1, ls=":", alpha=0.8)
        ax.axvline(hi_edge, color=color, lw=1, ls=":", alpha=0.8)
        ax.text(
            (lo_edge + hi_edge) / 2, 0.05,
            f"{lo_edge:.2f}–{hi_edge:.2f} kPa",
            ha="center", fontsize=8, color=color
        )

    ax.axvspan(2, 8, alpha=0.06, color="green", label="Physiological gingiva")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("ECM Stiffness (kPa)", fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0, 10)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

axes[0].set_ylabel("Steady-State Regenerative Fraction R*", fontsize=10)

plt.tight_layout()
plt.savefig("figures/sim1_feedback_knockout.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim1_feedback_knockout.png")

# ── Summary stats ─────────────────────────────────────────────────────────────
gap_wt = R_hi_wt - R_lo_wt
gap_ko = R_hi_ko - R_lo_ko
print(f"\nWild-type   max hysteresis gap : {gap_wt.max():.3f}")
print(f"Full KO     max hysteresis gap : {gap_ko.max():.3f}")
print(f"WT bistable window             : "
      f"{S_range[gap_wt>0.1][0]:.2f}–{S_range[gap_wt>0.1][-1]:.2f} kPa")
bistable_ko = S_range[gap_ko > 0.1]
print(f"KO bistable window             : "
      f"{'NONE — bistability abolished' if len(bistable_ko)==0 else str(bistable_ko[[0,-1]])}")
#..
