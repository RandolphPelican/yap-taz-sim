"""
run_all_sims.py — Execute all 4 validation simulations in sequence
===================================================================
Usage: python run_all_sims.py

Outputs figures to ./figures/
Runtime: ~2-5 minutes on standard hardware.

Author: John D. Stabler | MIT License
"""

import os, time, subprocess, sys

os.makedirs("figures", exist_ok=True)

SIMS = [
    ("sim1_feedback_knockout.py",  "Sim 1: Feedback Knockout"),
    ("sim2_imperfect_washout.py",  "Sim 2: Imperfect Washout"),
    ("sim3_agent_based.py",        "Sim 3: Agent-Based (Spatial)"),
    ("sim4_stochastic_cle.py",     "Sim 4: Stochastic CLE"),
]

print("=" * 60)
print("YAP/TAZ Bistable Switch — Validation Simulation Suite")
print("=" * 60)

total_start = time.time()
for script, label in SIMS:
    print(f"\n▶  Running {label}...")
    t0 = time.time()
    result = subprocess.run([sys.executable, script], capture_output=True, text=True)
    elapsed = time.time() - t0
    if result.returncode == 0:
        print(result.stdout)
        print(f"   ✓ Completed in {elapsed:.1f}s")
    else:
        print(f"   ✗ ERROR:\n{result.stderr}")

total_elapsed = time.time() - total_start
print("\n" + "=" * 60)
print(f"All simulations complete in {total_elapsed:.1f}s")
print("Figures saved to ./figures/")
print("=" * 60)
#..
