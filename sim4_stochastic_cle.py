"""
Simulation 4: Stochastic (CLE) — Is the switch noise-robust?
=============================================================
Chemical Langevin Equation adds multiplicative Gaussian noise to the ODE.
Simulates N=1000 independent cells to get population distributions.

Key question: does intrinsic noise cause unwanted random state-flipping?

Prediction: bimodal distribution at bistable stiffness; near-zero spontaneous
switching rate confirms hysteresis is not noise-driven drift.

Author: John D. Stabler | MIT License
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import PARAMS, f_stiffness, sigma_feedback

N_CELLS   = 1000          # independent cells
DT        = 0.1           # time step (h)
T_MAX     = 500           # total sim time (h)
ETA       = 0.05          # noise amplitude (tune to ~5% of typical dR)
SEED      = 2025
np.random.seed(SEED)

K_L = 0.02
L0  = 2.0

def deterministic_dRdt(R, S, p, L=0.0):
    R = np.clip(R, 1e-4, 1 - 1e-4)
    fS = f_stiffness(S, p["S0"], p["sigma"])
    fb = sigma_feedback(R, p["theta"], p["n"])
    g  = K_L / (L + K_L)
    return p["alpha"] * fS * fb * (1 - R) - p["beta"] * g * (1 - fS) * R


def run_cle(S, R0_arr, p, eta=ETA, drug_pulse_end=None):
    """
    Stochastic Euler-Maruyama integration (CLE).
    R0_arr : array of initial conditions, shape (N,)
    Returns: R trajectory, shape (n_steps, N)
    """
    n_steps = int(T_MAX / DT)
    R = R0_arr.copy().astype(float)
    trajectories = np.zeros((n_steps, len(R)))

    for step in range(n_steps):
        t = step * DT
        L = (L0 if (drug_pulse_end and t < drug_pulse_end) else 0.0)
        drift     = deterministic_dRdt(R, S, p, L)
        diffusion = eta * np.sqrt(np.abs(drift) + 1e-6)
        noise     = diffusion * np.random.randn(len(R))
        R         = np.clip(R + DT * drift + np.sqrt(DT) * noise, 0, 1)
        trajectories[step] = R

    return trajectories


# ── 1. Bimodal distribution test at S=4 kPa (bistable zone) ─────────────────
print("Running CLE @ S=4 kPa with mixed ICs (500 fibrotic / 500 regenerative)...")
R0_mixed = np.concatenate([
    np.random.uniform(0.01, 0.15, 500),   # fibrotic
    np.random.uniform(0.80, 0.99, 500),   # regenerative
])
traj_bistable = run_cle(S=4.0, R0_arr=R0_mixed, p=PARAMS)
R_final_bistable = traj_bistable[-1]

# ── 2. Monostable test at S=0.5 kPa (below bistable window) ─────────────────
print("Running CLE @ S=0.5 kPa (below window, should be monostable fibrotic)...")
traj_low = run_cle(S=0.5, R0_arr=R0_mixed.copy(), p=PARAMS)
R_final_low = traj_low[-1]

# ── 3. High stiffness test at S=9.5 kPa (above window, monostable fibrotic) ─
print("Running CLE @ S=9.5 kPa (above window, should be monostable fibrotic)...")
traj_high = run_cle(S=9.5, R0_arr=R0_mixed.copy(), p=PARAMS)
R_final_high = traj_high[-1]

# ── 4. Spontaneous switching rate ─────────────────────────────────────────────
print("Computing spontaneous switching rate @ S=4 kPa (all fibrotic start)...")
R0_all_fibrotic = np.random.uniform(0.01, 0.1, N_CELLS)
traj_switch = run_cle(S=4.0, R0_arr=R0_all_fibrotic, p=PARAMS)
# Count cells that crossed R=0.5 spontaneously
crossings = (traj_switch > 0.5).any(axis=0).sum()
switch_rate_per_day = crossings / N_CELLS / (T_MAX / 24) * 100

# ── Figures ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle(
    "Simulation 4: Stochastic CLE — Noise Robustness of Bistable Switch",
    fontsize=13, fontweight="bold"
)

# Panel A: Bimodal histogram at bistable stiffness
ax = axes[0, 0]
ax.hist(R_final_bistable, bins=50, color="#1565C0", edgecolor="white", lw=0.5, alpha=0.85)
ax.axvline(0.5, color="red", ls="--", lw=1.5, label="State threshold")
ax.set_title("A. S=4 kPa (Bistable Zone)\nExpected: Bimodal Distribution", fontweight="bold")
ax.set_xlabel("Final R"); ax.set_ylabel("Cell Count")
ax.legend(); ax.grid(alpha=0.3)
bimodal_pct = (R_final_bistable > 0.5).mean() * 100
ax.text(0.7, ax.get_ylim()[1]*0.85, f"Regen: {bimodal_pct:.0f}%\nFib: {100-bimodal_pct:.0f}%",
        fontsize=10, bbox=dict(boxstyle="round", fc="white", alpha=0.8))

# Panel B: Monostable at S=0.5 kPa
ax = axes[0, 1]
ax.hist(R_final_low, bins=50, color="#C62828", edgecolor="white", lw=0.5, alpha=0.85)
ax.axvline(0.5, color="black", ls="--", lw=1.5, label="State threshold")
ax.set_title("B. S=0.5 kPa (Below Window)\nExpected: Fibrotic Monostable", fontweight="bold")
ax.set_xlabel("Final R"); ax.set_ylabel("Cell Count")
ax.legend(); ax.grid(alpha=0.3)

# Panel C: Monostable at S=9.5 kPa
ax = axes[1, 0]
ax.hist(R_final_high, bins=50, color="#E65100", edgecolor="white", lw=0.5, alpha=0.85)
ax.axvline(0.5, color="black", ls="--", lw=1.5, label="State threshold")
ax.set_title("C. S=9.5 kPa (Above Window)\nExpected: Fibrotic Monostable (high stiffness)", fontweight="bold")
ax.set_xlabel("Final R"); ax.set_ylabel("Cell Count")
ax.legend(); ax.grid(alpha=0.3)

# Panel D: Spontaneous switching trajectories (sample of 50 cells)
ax = axes[1, 1]
t_axis = np.linspace(0, T_MAX, traj_switch.shape[0])
for i in range(50):
    ax.plot(t_axis, traj_switch[:, i], alpha=0.15, lw=0.7, color="#1B5E20")
ax.axhline(0.5, color="red", ls="--", lw=1.5, label="State threshold")
ax.set_title(
    f"D. Spontaneous Switching (all fibrotic start, S=4 kPa)\n"
    f"Switch rate: {switch_rate_per_day:.2f}% per day — {'STABLE ✓' if switch_rate_per_day < 0.5 else 'UNSTABLE ✗'}",
    fontweight="bold",
    color="#1B5E20" if switch_rate_per_day < 0.5 else "#B71C1C"
)
ax.set_xlabel("Time (h)"); ax.set_ylabel("R(t)")
ax.set_ylim(-0.05, 1.05)
ax.legend(); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figures/sim4_stochastic_cle.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim4_stochastic_cle.png")

# ── Mean first passage time estimate ──────────────────────────────────────────
# (analytical approximation via Kramers' theory, illustrative)
# For a double-well with barrier height ΔU and noise η, MFPT ~ exp(2ΔU/η²)/ω
# We approximate ΔU from the bistability gap
S_test   = 4.0
gap_est  = 0.75 - 0.20   # rough bistability gap
delta_U  = 0.5 * gap_est ** 2 / (ETA ** 2 + 1e-6)
mfpt_h   = np.exp(min(delta_U, 50)) / (2 * np.pi)   # cap to avoid overflow
mfpt_yr  = mfpt_h / (24 * 365)

print(f"\n── Stochastic Summary ─────────────────────────────────────")
print(f"  N cells simulated              : {N_CELLS}")
print(f"  Noise amplitude (η)            : {ETA}")
print(f"  S=4 kPa bimodal distribution?  : {'YES ✓' if (R_final_bistable > 0.5).mean() > 0.2 else 'NO ✗'}")
print(f"  % cells regenerative (S=4 kPa): {bimodal_pct:.1f}%")
print(f"  Spontaneous switch rate        : {switch_rate_per_day:.3f}% per day")
print(f"  Estimated MFPT (Kramers)       : {mfpt_yr:.2e} years")
print(f"  Noise-robustness verdict       : {'ROBUST — switch stable against thermal noise ✓' if switch_rate_per_day < 1 else 'FRAGILE ✗'}")
#..
