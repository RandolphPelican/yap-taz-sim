"""
Simulation 2: Imperfect Washout — Does hysteresis survive sloppy drug clearance?
==================================================================================
Models exponential drug decay post-pulse with variable half-lives.
Identifies the maximum tolerable clearance time for successful state lock-in.

Clinical relevance: real pharmacokinetics aren't step functions.

Author: John D. Stabler | MIT License
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from model import PARAMS, f_stiffness, sigma_feedback

K_L       = 0.02    # TRULI IC₅₀ (μM)
L0        = 2.0     # initial drug dose (μM)
T_PULSE   = 48.0    # pulse duration (h)
S_sim     = 4.0     # scaffold stiffness (kPa) — within therapeutic window
T_MAX     = 300     # total sim time (h)
N_STEPS   = 6000

# Half-lives to compare (hours)
HALF_LIVES = [1, 6, 12, 24, 48, 96]
COLORS     = ["#1a237e", "#1565C0", "#42A5F5", "#81D4FA", "#FFB74D", "#E53935"]


def drug_concentration(t, t_pulse=T_PULSE, L0=L0, half_life=1.0):
    """Exponential decay after pulse ends."""
    if t <= t_pulse:
        return L0
    k_clear = np.log(2) / half_life
    return L0 * np.exp(-k_clear * (t - t_pulse))


def full_ode(state, t, S, p, half_life):
    R = np.clip(state[0], 0, 1)
    L = drug_concentration(t, half_life=half_life)
    fS = f_stiffness(S, p["S0"], p["sigma"])
    fb = sigma_feedback(R, p["theta"], p["n"])
    g  = K_L / (L + K_L)
    dR = p["alpha"] * fS * fb * (1 - R) - p["beta"] * g * (1 - fS) * R
    return [dR]


t_span = np.linspace(0, T_MAX, N_STEPS)

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(
    "Simulation 2: Imperfect Washout — Hysteresis Robustness to Drug Clearance Kinetics",
    fontsize=13, fontweight="bold", y=1.01
)

results_summary = []

for idx, (hl, color, ax) in enumerate(zip(HALF_LIVES, COLORS, axes.flatten())):
    sol = odeint(full_ode, [0.05], t_span, args=(S_sim, PARAMS, hl))
    R_traj = sol[:, 0]

    # Drug concentration trajectory for overlay
    L_traj = np.array([drug_concentration(t, half_life=hl) for t in t_span])

    ax2 = ax.twinx()
    ax2.fill_between(t_span, 0, L_traj, alpha=0.15, color="#FF7043", label="[Drug]")
    ax2.set_ylabel("[TRULI] μM", fontsize=8, color="#FF7043")
    ax2.tick_params(axis="y", labelcolor="#FF7043", labelsize=7)
    ax2.set_ylim(0, L0 * 1.3)

    ax.plot(t_span, R_traj, color=color, lw=2.5, zorder=5)
    ax.axhline(0.5, color="gray", ls="--", lw=1, alpha=0.7, label="Threshold R=0.5")
    ax.axhline(0.7, color="green", ls=":", lw=1, alpha=0.7, label="Target R=0.7")
    ax.axvline(T_PULSE, color="black", ls="--", lw=1, alpha=0.5)

    R_final = R_traj[-1]
    success = R_final >= 0.6
    status  = "✓ SUCCESS" if success else "✗ FAILURE"
    ax.set_facecolor("#E8F5E9" if success else "#FFEBEE")
    ax.set_title(
        f"Drug t½ = {hl}h  |  {status}\nR(final) = {R_final:.3f}",
        fontsize=10, fontweight="bold",
        color="#1B5E20" if success else "#B71C1C"
    )
    ax.set_xlabel("Time (h)", fontsize=9)
    ax.set_ylabel("Regenerative Fraction R(t)", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0, T_MAX)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")

    results_summary.append((hl, R_final, success))

plt.tight_layout()
plt.savefig("figures/sim2_imperfect_washout.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim2_imperfect_washout.png")

# ── Critical clearance window summary plot ────────────────────────────────────
fig2, ax = plt.subplots(figsize=(8, 4))
hls   = [r[0] for r in results_summary]
R_fin = [r[1] for r in results_summary]
cols  = ["#2E7D32" if r[2] else "#C62828" for r in results_summary]
bars  = ax.bar(range(len(hls)), R_fin, color=cols, edgecolor="white", lw=1.5)
ax.axhline(0.6, color="black", ls="--", lw=1.5, label="Success threshold (R=0.6)")
ax.set_xticks(range(len(hls)))
ax.set_xticklabels([f"t½={h}h" for h in hls])
ax.set_ylabel("Final Regenerative Fraction R(300h)")
ax.set_title("Critical Washout Window: Maximum Allowable Drug Half-Life", fontweight="bold")
ax.set_ylim(0, 1.05)
ax.legend()
ax.grid(axis="y", alpha=0.3)
for bar, val in zip(bars, R_fin):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.02,
            f"{val:.2f}", ha="center", fontsize=9, fontweight="bold")
plt.tight_layout()
plt.savefig("figures/sim2_clearance_summary.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim2_clearance_summary.png")

print("\n── Washout Robustness Summary ──────────────────────────────")
for hl, Rf, ok in results_summary:
    print(f"  t½ = {hl:>3}h  →  R_final = {Rf:.3f}  {'✓' if ok else '✗'}")
#..
