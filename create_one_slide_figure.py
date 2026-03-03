"""
One-slide visual summary for PI outreach emails
Combines: mechanism schematic + phase diagram + key result
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from scipy.integrate import odeint

# Reuse model parameters
alpha, beta = 2.0, 0.5
S0, sigma = 5.0, 2.0
theta, n = 0.3, 4

def f(S): return np.exp(-((S - S0)**2) / (2 * sigma**2))
def sigmoid(R): return R**n / (theta**n + R**n)
def dR_dt(R, t, S):
    R = np.clip(R, 0, 1)
    return alpha * f(S) * sigmoid(R) * (1 - R) - beta * (1 - f(S)) * R

# Generate hysteresis data
S_range = np.linspace(0, 10, 100)
t_eval = np.linspace(0, 500, 5000)
R_low, R_high = [], []
for S in S_range:
    R_low.append(odeint(dR_dt, 0.01, t_eval, args=(S,))[-1, 0])
    R_high.append(odeint(dR_dt, 0.99, t_eval, args=(S,))[-1, 0])

# Create figure
fig = plt.figure(figsize=(16, 9), facecolor='white')
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3, 
                      left=0.08, right=0.95, top=0.92, bottom=0.08)

# Title
fig.suptitle('YAP/TAZ Bistable Switch: A New Target for Gingival Regeneration', 
             fontsize=20, fontweight='bold', y=0.97)

# Panel A: The Problem
ax_problem = fig.add_subplot(gs[0, 0])
ax_problem.text(0.5, 0.85, 'THE PROBLEM', ha='center', fontsize=14, fontweight='bold')
ax_problem.text(0.5, 0.65, '40M US adults\nwith gingival recession', 
                ha='center', fontsize=12, bbox=dict(boxstyle='round', fc='#ffcccc', alpha=0.8))
ax_problem.text(0.5, 0.4, 'Current treatment:\nSurgical grafting', 
                ha='center', fontsize=11, style='italic')
ax_problem.text(0.5, 0.2, '❌ Painful\n❌ Donor tissue limited\n❌ 30-40% failure rate', 
                ha='center', fontsize=10, color='darkred')
ax_problem.text(0.5, 0.02, 'No pharmacological option exists', 
                ha='center', fontsize=11, fontweight='bold', color='red')
ax_problem.axis('off')

# Panel B: The Mechanism
ax_mech = fig.add_subplot(gs[0, 1:])
ax_mech.text(0.5, 0.95, 'THE MECHANISM: YAP/TAZ as Bistable Mechanosensitive Switch', 
             ha='center', fontsize=14, fontweight='bold', transform=ax_mech.transAxes)

# Draw three tissue states
states = [
    (0.15, 'TOO SOFT\n(<1.3 kPa)', 'Quiescent\nYAP OFF\nNo repair', '#ccccff'),
    (0.5, 'OPTIMAL\n(1.3-8.7 kPa)', 'Regenerative\nYAP ON\n✓ Boundary forms', '#ccffcc'),
    (0.85, 'TOO STIFF\n(>8.7 kPa)', 'Fibrotic\nYAP hyperactive\nNo repair', '#ffcccc')
]

for x, title, desc, color in states:
    circle = Circle((x, 0.5), 0.12, transform=ax_mech.transAxes, 
                   fc=color, ec='black', linewidth=2)
    ax_mech.add_patch(circle)
    ax_mech.text(x, 0.78, title, ha='center', fontsize=10, fontweight='bold',
                transform=ax_mech.transAxes)
    ax_mech.text(x, 0.35, desc, ha='center', fontsize=9, va='top',
                transform=ax_mech.transAxes)

# Add feedback arrow
arrow = FancyArrowPatch((0.5, 0.15), (0.5, 0.35), 
                       transform=ax_mech.transAxes,
                       arrowstyle='->', mutation_scale=20, linewidth=2.5, 
                       color='green', connectionstyle="arc3,rad=.3")
ax_mech.add_patch(arrow)
ax_mech.text(0.58, 0.25, 'Positive\nFeedback', ha='left', fontsize=9, 
            color='green', fontweight='bold', transform=ax_mech.transAxes)

ax_mech.text(0.5, 0.05, 'Adult gingiva is STUCK in fibrotic state (right)\nPositive feedback creates HYSTERESIS', 
             ha='center', fontsize=11, style='italic', color='darkred',
             transform=ax_mech.transAxes)

ax_mech.set_xlim(0, 1)
ax_mech.set_ylim(0, 1)
ax_mech.axis('off')

# Panel C: Hysteresis Curve (KEY RESULT)
ax_hyst = fig.add_subplot(gs[1, :2])
ax_hyst.fill_between(S_range, 0, f(S_range), alpha=0.15, color='green', 
                     label='Stiffness activation window')
ax_hyst.plot(S_range, R_high, 'r-', linewidth=3.5, label='Upper branch (regenerative)', zorder=3)
ax_hyst.plot(S_range, R_low, 'b--', linewidth=3.5, label='Lower branch (fibrotic)', zorder=3)

# Annotate bistable region
bistable_mask = np.abs(np.array(R_high) - np.array(R_low)) > 0.1
if np.any(bistable_mask):
    S_bi = S_range[bistable_mask]
    ax_hyst.axvspan(S_bi[0], S_bi[-1], alpha=0.2, color='orange', zorder=1)
    ax_hyst.text(5, 0.92, f'BISTABLE REGION\n{S_bi[0]:.1f} - {S_bi[-1]:.1f} kPa', 
                ha='center', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round', fc='orange', alpha=0.6))

ax_hyst.axvline(5, color='green', linestyle=':', linewidth=2.5, alpha=0.7)
ax_hyst.text(5, -0.15, 'Optimal\nStiffness', ha='center', fontsize=10, 
            color='green', fontweight='bold')

ax_hyst.set_xlabel('Tissue Stiffness (kPa)', fontsize=13, fontweight='bold')
ax_hyst.set_ylabel('Regenerative State', fontsize=13, fontweight='bold')
ax_hyst.set_title('Mathematical Prediction: Hysteresis Creates Therapeutic Window', 
                 fontsize=13, fontweight='bold', pad=15)
ax_hyst.legend(loc='upper left', fontsize=11, framealpha=0.95)
ax_hyst.grid(alpha=0.3, linestyle='--')
ax_hyst.set_ylim(-0.05, 1.05)
ax_hyst.set_xlim(0, 10)

# Panel D: The Solution
ax_sol = fig.add_subplot(gs[1, 2])
ax_sol.text(0.5, 0.95, 'THE SOLUTION', ha='center', fontsize=14, fontweight='bold',
           transform=ax_sol.transAxes)

solution_text = """
1️⃣ LATS Inhibitor (TRULI)
   48h pulse → YAP activates

2️⃣ Soft Hydrogel (4 kPa)
   Positions tissue in
   regenerative window

3️⃣ Self-Terminating
   Boundary forms →
   Drug clears →
   YAP shuts off →
   REPAIR PERSISTS

✓ Non-surgical
✓ Single application
✓ No chronic dosing
"""

ax_sol.text(0.5, 0.5, solution_text, ha='center', va='center', fontsize=10,
           transform=ax_sol.transAxes, family='monospace',
           bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.3))

ax_sol.text(0.5, 0.05, 'Ready to test in\nex vivo explants\n(12-18 months, ~$220K)', 
           ha='center', fontsize=10, fontweight='bold', style='italic',
           transform=ax_sol.transAxes,
           bbox=dict(boxstyle='round', fc='yellow', alpha=0.4))

ax_sol.set_xlim(0, 1)
ax_sol.set_ylim(0, 1)
ax_sol.axis('off')

plt.savefig('figures/ONE_SLIDE_SUMMARY.png', dpi=300, bbox_inches='tight', facecolor='white')
print("\n✓ One-slide summary created: figures/ONE_SLIDE_SUMMARY.png")
print("  → Embed this in your email to PIs (instant visual understanding)")
