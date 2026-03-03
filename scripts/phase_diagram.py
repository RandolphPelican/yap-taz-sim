import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import odeint

# Parameters (same as phase_sim.py)
k_act, k_inact = 0.5, 1.0
K_E, K_L = 5.0, 0.02
k_remodel, k_relax = 0.05, 0.01
E_min, E_max, E_opt, sigma_E = 1, 15, 5, 3
k_form, k_decay = 0.1, 0.02
k_clear = 0.03
W_base, alpha_Y, beta_B = 0.2, 3, 0.5

def f_E(E): return E / (E + K_E)
def g_L(L): return K_L / (L + K_L)
def h_E(E): return np.exp(-((E - E_opt)**2) / (2 * sigma_E**2))
def E_target(B): return E_min + (E_max - E_min) * B
def I_drug(t, dose, T_pulse=48): return dose if t < T_pulse else 0.0

def system(state, t, dose, E0):
    Y, E, B, L = state
    W = W_base * (1 + alpha_Y * Y)
    dY_dt = k_act * f_E(E) * (1 - Y) - k_inact * (g_L(L) + beta_B * B) * Y
    dE_dt = k_remodel * Y * (E_target(B) - E) + k_relax * (E0 - E)
    dB_dt = k_form * Y * h_E(E) * W * (1 - B) - k_decay * (1 - Y) * B
    dL_dt = -k_clear * L + I_drug(t, dose)
    return [dY_dt, dE_dt, dB_dt, dL_dt]

print("\n" + "="*60)
print("  Phase Diagram Generator")
print("="*60 + "\n")

# Create phase space grid
Y_range = np.linspace(0, 1, 50)
E_range = np.linspace(0.5, 16, 50)
Y_grid, E_grid = np.meshgrid(Y_range, E_range)

# Calculate boundary formation rate at each point (dB/dt with B=0)
B_rate = np.zeros_like(Y_grid)
for i in range(len(E_range)):
    for j in range(len(Y_range)):
        Y_val, E_val = Y_grid[i,j], E_grid[i,j]
        W = W_base * (1 + alpha_Y * Y_val)
        B_rate[i,j] = k_form * Y_val * h_E(E_val) * W

# Define regions
def classify_region(Y, E):
    if Y < 0.2 and E > 8: return 1  # Fibrotic
    elif Y > 0.6 and E > 9: return 2  # Hyperplastic
    elif Y < 0.2 and E < 3: return 3  # Quiescent
    elif 0.2 <= Y <= 0.7 and 3 <= E <= 9: return 4  # Regenerative
    else: return 0  # Transition

# Simulate trajectories
t = np.linspace(0, 200, 2000)
trajectories = [
    ("Success", 2.0, 4, 'green'),
    ("Failed", 0.0, 12, 'red'),
    ("Drug-only", 2.0, 12, 'orange'),
    ("Scaffold-only", 0.0, 4, 'blue')
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Phase diagram with trajectories
contour = ax1.contourf(Y_grid, E_grid, B_rate, levels=20, cmap='viridis', alpha=0.3)
ax1.contour(Y_grid, E_grid, B_rate, levels=[0.01], colors='white', linewidths=2)

# Region boundaries (approximate)
ax1.add_patch(Rectangle((0, 8), 0.2, 8, fc='lightcoral', alpha=0.2, label='Region I (Fibrotic)'))
ax1.add_patch(Rectangle((0.6, 9), 0.4, 7, fc='orange', alpha=0.2, label='Region II (Hyperplastic)'))
ax1.add_patch(Rectangle((0, 0.5), 0.2, 2.5, fc='lightblue', alpha=0.2, label='Region III (Quiescent)'))
ax1.add_patch(Rectangle((0.2, 3), 0.5, 6, fc='lightgreen', alpha=0.3, label='Region IV (Regenerative)'))

# Plot trajectories
for label, dose, E0_val, color in trajectories:
    sol = odeint(system, [0.05, 12 if E0_val==12 else 4, 0, 0], t, args=(dose, E0_val))
    Y_traj, E_traj, B_traj = sol[:, 0], sol[:, 1], sol[:, 2]
    
    ax1.plot(Y_traj, E_traj, color=color, linewidth=2.5, label=label, alpha=0.8)
    ax1.scatter(Y_traj[0], E_traj[0], color=color, s=100, marker='o', edgecolors='black', zorder=5)
    ax1.scatter(Y_traj[-1], E_traj[-1], color=color, s=150, marker='*', edgecolors='black', zorder=5)
    
    # Annotate final B value
    ax1.annotate(f'B={B_traj[-1]:.2f}', 
                xy=(Y_traj[-1], E_traj[-1]), 
                xytext=(10, -10), 
                textcoords='offset points',
                fontsize=9, 
                bbox=dict(boxstyle='round,pad=0.3', fc=color, alpha=0.3))

ax1.set_xlabel('YAP Activity (nuclear fraction)', fontsize=13, fontweight='bold')
ax1.set_ylabel('ECM Stiffness (kPa)', fontsize=13, fontweight='bold')
ax1.set_title('Phase Diagram: Tissue State Space', fontsize=14, fontweight='bold')
ax1.legend(loc='upper left', fontsize=9, framealpha=0.9)
ax1.grid(alpha=0.3)
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 16)

# Right: Time evolution of key metrics
for label, dose, E0_val, color in trajectories[:3]:  # Only show first 3 for clarity
    sol = odeint(system, [0.05, 12 if E0_val==12 else 4, 0, 0], t, args=(dose, E0_val))
    B_traj = sol[:, 2]
    ax2.plot(t/24, B_traj, color=color, linewidth=2.5, label=label, alpha=0.8)

ax2.axhline(y=0.7, color='green', linestyle='--', linewidth=2, alpha=0.5, label='Repair Threshold (B>0.7)')
ax2.axvline(x=48/24, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Drug Washout (48h)')
ax2.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
ax2.set_ylabel('Boundary Integrity (B)', fontsize=13, fontweight='bold')
ax2.set_title('Temporal Dynamics: Does Repair Persist?', fontsize=14, fontweight='bold')
ax2.legend(loc='lower right', fontsize=10, framealpha=0.9)
ax2.grid(alpha=0.3)
ax2.set_xlim(0, 8)
ax2.set_ylim(0, 1)

plt.tight_layout()
plt.savefig('figures/phase_diagram.png', dpi=300, bbox_inches='tight')
print("✓ Phase diagram saved: figures/phase_diagram.png")

# Calculate key metrics
print("\n" + "="*60)
print("  Quantitative Analysis")
print("="*60)

for label, dose, E0_val, _ in trajectories:
    sol = odeint(system, [0.05, 12 if E0_val==12 else 4, 0, 0], t, args=(dose, E0_val))
    B_traj = sol[:, 2]
    
    # Time to reach B > 0.7
    idx_threshold = np.where(B_traj > 0.7)[0]
    time_to_repair = t[idx_threshold[0]]/24 if len(idx_threshold) > 0 else np.inf
    
    # Final B at day 8
    B_final = B_traj[-1]
    
    # Peak YAP (cancer risk proxy)
    Y_peak = np.max(sol[:, 0])
    
    print(f"\n{label:20s}: B_final={B_final:.3f}  |  Time to B>0.7: {time_to_repair:5.2f}d  |  Y_peak={Y_peak:.3f}")

print("\n" + "="*60)
print("Done.")
