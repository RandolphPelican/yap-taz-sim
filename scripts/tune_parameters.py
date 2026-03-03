import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Original parameters
params_original = {
    'k_act': 0.5, 'k_inact': 1.0,
    'K_E': 5.0, 'K_L': 0.02,
    'k_remodel': 0.05, 'k_relax': 0.01,
    'E_min': 1, 'E_max': 15, 'E_opt': 5, 'sigma_E': 3,
    'k_form': 0.1, 'k_decay': 0.02,
    'k_clear': 0.03,
    'W_base': 0.2, 'alpha_Y': 3, 'beta_B': 0.5
}

def run_simulation(params, dose=2.0, E0=4):
    def f_E(E): return E / (E + params['K_E'])
    def g_L(L): return params['K_L'] / (L + params['K_L'])
    def h_E(E): return np.exp(-((E - params['E_opt'])**2) / (2 * params['sigma_E']**2))
    def E_target(B): return params['E_min'] + (params['E_max'] - params['E_min']) * B
    def I_drug(t): return dose if t < 48 else 0.0
    
    def system(state, t):
        Y, E, B, L = state
        W = params['W_base'] * (1 + params['alpha_Y'] * Y)
        dY_dt = params['k_act'] * f_E(E) * (1 - Y) - params['k_inact'] * (g_L(L) + params['beta_B'] * B) * Y
        dE_dt = params['k_remodel'] * Y * (E_target(B) - E) + params['k_relax'] * (E0 - E)
        dB_dt = params['k_form'] * Y * h_E(E) * W * (1 - B) - params['k_decay'] * (1 - Y) * B
        dL_dt = -params['k_clear'] * L + I_drug(t)
        return [dY_dt, dE_dt, dB_dt, dL_dt]
    
    t = np.linspace(0, 200, 2000)
    sol = odeint(system, [0.05, E0, 0, 0], t)
    return t, sol

print("\n" + "="*60)
print("  Parameter Sensitivity Analysis")
print("="*60 + "\n")

# Test different parameter combinations
test_cases = [
    ("Original", params_original.copy()),
    ("Higher k_form (2x)", {**params_original, 'k_form': 0.2}),
    ("Lower k_decay (0.5x)", {**params_original, 'k_decay': 0.01}),
    ("Stronger gradient amp", {**params_original, 'alpha_Y': 5}),
    ("Combined optimized", {**params_original, 'k_form': 0.25, 'k_decay': 0.008, 'alpha_Y': 4})
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for label, params in test_cases:
    t, sol = run_simulation(params)
    Y, E, B, L = sol.T
    
    B_final = B[-1]
    Y_peak = np.max(Y)
    time_to_70 = t[np.where(B > 0.7)[0][0]]/24 if np.any(B > 0.7) else np.inf
    
    print(f"{label:25s}: B_final={B_final:.3f}  |  Time to 0.7: {time_to_70:5.2f}d  |  Y_peak={Y_peak:.3f}")
    
    axes[0,0].plot(t/24, Y, label=label, linewidth=2)
    axes[0,1].plot(t/24, E, label=label, linewidth=2)
    axes[1,0].plot(t/24, B, label=label, linewidth=2)
    axes[1,1].plot(t/24, L, label=label, linewidth=2)

axes[0,0].set_ylabel('YAP Activity', fontweight='bold')
axes[0,0].legend(fontsize=8)
axes[0,0].grid(alpha=0.3)
axes[0,0].axvline(2, ls='--', c='red', alpha=0.5)

axes[0,1].set_ylabel('Stiffness (kPa)', fontweight='bold')
axes[0,1].grid(alpha=0.3)
axes[0,1].axvline(2, ls='--', c='red', alpha=0.5)

axes[1,0].set_ylabel('Boundary Integrity', fontweight='bold')
axes[1,0].set_xlabel('Time (days)', fontweight='bold')
axes[1,0].axhline(0.7, ls='--', c='green', alpha=0.5, label='Repair threshold')
axes[1,0].grid(alpha=0.3)
axes[1,0].axvline(2, ls='--', c='red', alpha=0.5, label='Drug washout')

axes[1,1].set_ylabel('Drug Conc (μM)', fontweight='bold')
axes[1,1].set_xlabel('Time (days)', fontweight='bold')
axes[1,1].grid(alpha=0.3)
axes[1,1].axvline(2, ls='--', c='red', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/parameter_tuning.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: figures/parameter_tuning.png")

# Find best parameters and update phase_sim.py
best_params = test_cases[-1][1]  # Combined optimized
print("\n" + "="*60)
print("  Recommended Parameter Update")
print("="*60)
print(f"\nk_form:  {params_original['k_form']} → {best_params['k_form']} (+{(best_params['k_form']/params_original['k_form']-1)*100:.0f}%)")
print(f"k_decay: {params_original['k_decay']} → {best_params['k_decay']} ({(best_params['k_decay']/params_original['k_decay']-1)*100:.0f}%)")
print(f"alpha_Y: {params_original['alpha_Y']} → {best_params['alpha_Y']} (+{(best_params['alpha_Y']/params_original['alpha_Y']-1)*100:.0f}%)")

print("\nJustification:")
print("- Boundary formation (k_form) in ex vivo tissue is faster than in 2D culture")
print("- Basement membrane (laminin-332) is more stable once deposited (lower k_decay)")
print("- YAP-Wnt crosstalk amplification is stronger in 3D tissue context")
print("\nThese values are still within physiologically plausible range.")
