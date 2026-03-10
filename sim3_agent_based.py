"""
Simulation 3: Agent-Based Model — Does bistability survive spatial coupling?
=============================================================================
2D numpy grid of cells. Each cell runs the ODE internally.
Cells share local stiffness via ECM secretion from neighbors.
No external frameworks needed — pure numpy + matplotlib.

Prediction: local cooperative coupling stabilizes regenerative clusters,
showing bistability is a tissue-level property, not just single-cell.

Author: John D. Stabler | MIT License
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from model import PARAMS, f_stiffness, sigma_feedback

# ── Grid configuration ────────────────────────────────────────────────────────
GRID_SIZE   = 30          # 30×30 cells
DT          = 0.5         # time step (h)
T_MAX       = 200         # total simulation time (h)
T_PULSE_END = 48          # LATS inhibitor pulse end (h)
L0          = 2.0         # drug concentration during pulse (μM)
K_L         = 0.02        # IC₅₀
SEED        = 42
np.random.seed(SEED)

# Stiffness: uniform scaffold at 4 kPa + small noise
S_BASE       = 4.0
S_NOISE      = 0.3
S_GRID       = np.clip(np.random.normal(S_BASE, S_NOISE, (GRID_SIZE, GRID_SIZE)), 1, 9)

# ECM coupling: cells stiffen local environment proportional to R
ECM_COUPLING = 0.8        # how much a regenerative cell stiffens neighbors (kPa)
NEIGHBOR_K   = 0.15       # spatial coupling weight (Moore neighborhood average)


def moore_average(grid):
    """Average of 8 Moore neighbors (with zero-padding at edges)."""
    pad = np.pad(grid, 1, mode="edge")
    avg = np.zeros_like(grid)
    for di in range(3):
        for dj in range(3):
            if di == 1 and dj == 1:
                continue
            avg += pad[di:di+GRID_SIZE, dj:dj+GRID_SIZE]
    return avg / 8.0


def cell_dRdt(R, S_eff, p, L):
    """Single-cell ODE derivative."""
    R = np.clip(R, 0.001, 0.999)
    fS = f_stiffness(S_eff, p["S0"], p["sigma"])
    fb = sigma_feedback(R, p["theta"], p["n"])
    g  = K_L / (L + K_L)
    return p["alpha"] * fS * fb * (1 - R) - p["beta"] * g * (1 - fS) * R


# ── Initialise: small "seed" patch of pre-conditioned cells in centre ─────────
R = np.full((GRID_SIZE, GRID_SIZE), 0.05)           # almost all fibrotic
cx, cy = GRID_SIZE // 2, GRID_SIZE // 2
R[cx-2:cx+2, cy-2:cy+2] = 0.85                     # 4×4 regenerative seed

n_steps    = int(T_MAX / DT)
snapshots  = {}
snap_times = [0, 24, 48, 72, 100, 150, 200]

# ── Time evolution ────────────────────────────────────────────────────────────
for step in range(n_steps):
    t = step * DT

    # Drug concentration
    L = L0 if t < T_PULSE_END else 0.0

    # Local stiffness = base + ECM secretion from neighbors
    neighbor_R = moore_average(R)
    S_eff = S_GRID + ECM_COUPLING * neighbor_R * NEIGHBOR_K

    # Euler step for all cells simultaneously
    dR = cell_dRdt(R, S_eff, PARAMS, L)
    R  = np.clip(R + DT * dR, 0, 1)

    t_round = round(t)
    if t_round in snap_times and t_round not in snapshots:
        snapshots[t_round] = R.copy()

# Ensure final snapshot
snapshots[200] = R.copy()

# ── Figure: spatial maps ──────────────────────────────────────────────────────
snap_keys = sorted(snapshots.keys())
n_panels  = len(snap_keys)
fig, axes = plt.subplots(2, 4, figsize=(16, 9))
fig.suptitle(
    "Simulation 3: Agent-Based Model — Spatial Propagation of Regenerative State",
    fontsize=13, fontweight="bold"
)
axes = axes.flatten()

for idx, tk in enumerate(snap_keys):
    ax = axes[idx]
    im = ax.imshow(snapshots[tk], vmin=0, vmax=1, cmap="RdYlGn", origin="lower")
    mean_R = snapshots[tk].mean()
    frac_regen = (snapshots[tk] > 0.5).mean() * 100
    drug_label = " [DRUG ON]" if tk < T_PULSE_END else " [POST-WASHOUT]"
    ax.set_title(f"t = {tk}h{drug_label}\nMean R={mean_R:.2f} | {frac_regen:.0f}% regen",
                 fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="R")

# Hide unused panel
axes[-1].axis("off")
plt.tight_layout()
plt.savefig("figures/sim3_abm_spatial.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim3_abm_spatial.png")

# ── Figure: mean R over time with spatial variance ────────────────────────────
# Re-run to collect time series
R2 = np.full((GRID_SIZE, GRID_SIZE), 0.05)
R2[cx-2:cx+2, cy-2:cy+2] = 0.85
ts, means, stds, frac_regens = [], [], [], []

for step in range(n_steps):
    t = step * DT
    L = L0 if t < T_PULSE_END else 0.0
    neighbor_R2 = moore_average(R2)
    S_eff2 = S_GRID + ECM_COUPLING * neighbor_R2 * NEIGHBOR_K
    dR2 = cell_dRdt(R2, S_eff2, PARAMS, L)
    R2 = np.clip(R2 + DT * dR2, 0, 1)
    if step % 4 == 0:
        ts.append(t)
        means.append(R2.mean())
        stds.append(R2.std())
        frac_regens.append((R2 > 0.5).mean())

ts = np.array(ts); means = np.array(means); stds = np.array(stds)
frac_regens = np.array(frac_regens)

fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
fig2.suptitle("ABM: Population Dynamics — Spatial Cooperative Bistability",
              fontsize=12, fontweight="bold")

ax1.fill_between(ts, means - stds, means + stds, alpha=0.2, color="#43A047")
ax1.plot(ts, means, color="#2E7D32", lw=2.5, label="Mean R across tissue")
ax1.axvline(T_PULSE_END, color="red", ls="--", lw=1.5, label="Drug washout")
ax1.axhline(0.5, color="gray", ls=":", lw=1, label="Threshold")
ax1.set_ylabel("Mean Regenerative Fraction ⟨R⟩")
ax1.legend(fontsize=9); ax1.grid(alpha=0.3)

ax2.fill_between(ts, 0, frac_regens * 100, alpha=0.3, color="#1E88E5")
ax2.plot(ts, frac_regens * 100, color="#1565C0", lw=2.5)
ax2.axvline(T_PULSE_END, color="red", ls="--", lw=1.5, label="Drug washout")
ax2.set_xlabel("Time (h)"); ax2.set_ylabel("% Cells in Regenerative State (R>0.5)")
ax2.legend(fontsize=9); ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("figures/sim3_abm_dynamics.png", dpi=150, bbox_inches="tight")
print("✓ Saved: figures/sim3_abm_dynamics.png")

print(f"\n── ABM Summary ────────────────────────────────────────────")
print(f"  Grid size         : {GRID_SIZE}×{GRID_SIZE} = {GRID_SIZE**2} cells")
print(f"  Seed patch        : 4×4 (regenerative) in fibrotic background")
print(f"  Final mean R      : {means[-1]:.3f}")
print(f"  Final % regen     : {frac_regens[-1]*100:.1f}%")
print(f"  Spatial spread    : {'YES — bistability propagated' if frac_regens[-1]>0.5 else 'NO — patch contained'}")
#..
