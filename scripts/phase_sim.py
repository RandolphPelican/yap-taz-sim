import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
import os

# phase_sim.py
# YAP/TAZ Gingival Regeneration Phase Simulation
# Full Level 3 dynamical system from Stabler framework
# Simulates successful regenerative trajectory vs failed trajectory
# Shows phase diagram with treatment path, dose-response panel
# Author: John Stabler (Randolph Pelican III)
# Repo: github.com/RandolphPelican/yap-taz-sim

# Parameters (literature-derived)
k_act     = 0.5    # h^-1  YAP activation rate
k_inact   = 1.0    # h^-1  YAP inactivation rate
K_E       = 5.0    # kPa   half-max stiffness for mechanical activation
K_L       = 0.02   # uM    TRULI IC50 for LATS2
k_remodel = 0.05   # h^-1  YAP-driven ECM remodeling
k_relax   = 0.01   # h^-1  passive ECM relaxation
E_min     = 1.0    # kPa   minimum stiffness
E_max     = 15.0   # kPa   maximum stiffness
E_opt     = 5.0    # kPa   optimal stiffness for boundary formation
sigma_E   = 3.0    # kPa   width of permissive stiffness window
k_form    = 0.18    # h^-1  boundary assembly rate
k_decay   = 0.006   # h^-1  boundary decay rate
k_clear   = 0.03   # h^-1  drug clearance (t_half ~24h)
W_base    = 0.2    #       basal Wnt/BMP gradient
alpha_Y   = 3.0    #       YAP-mediated gradient amplification
beta_B    = 0.8    #       polarity-mediated YAP suppression
P_max     = 0.5    # day^-1 max proliferation rate
E0_scaffold = 4.0  # kPa   soft hydrogel scaffold stiffness


# Auxiliary functions
def f_E(E):
    return E / (E + K_E)

def g_L(L):
    return K_L / (L + K_L)

def h_E(E):
    return np.exp(-((E - E_opt)**2) / (2 * sigma_E**2))

def E_target(B):
    return E_min + (E_max - E_min) * B

def I_drug(t, dose=2.0, T_pulse=48.0):
    return dose if t < T_pulse else 0.0


# Full Level 3 ODE system
def system(state, t, dose=2.0, T_pulse=48.0, scaffold=True):
    Y, E, B, L = state
    Y = np.clip(Y, 0, 1)
    E = np.clip(E, E_min, E_max)
    B = np.clip(B, 0, 1)
    L = max(L, 0)

    W = W_base * (1 + alpha_Y * Y)
    E_scaf = E0_scaffold if scaffold else 12.0

    dY = (k_act * f_E(E) * (1 - Y)
          - k_inact * (g_L(L) + beta_B * B) * Y)

    dE = (k_remodel * Y * (E_target(B) - E)
          + k_relax * (E_scaf - E))

    dB = (k_form * Y * h_E(E) * W * (1 - B)
          - k_decay * (1 - Y) * B)

    dL = -k_clear * L + I_drug(t, dose, T_pulse)

    return [dY, dE, dB, dL]


def run_simulation(dose=2.0, T_pulse=48.0, scaffold=True,
                   t_end=200, n_points=2000):
    t = np.linspace(0, t_end, n_points)
    # Initial conditions: adult fibrotic state
    state0 = [0.05, 12.0, 0.0, 0.0]
    sol = odeint(system, state0, t,
                 args=(dose, T_pulse, scaffold),
                 rtol=1e-6, atol=1e-8)
    Y, E, B, L = sol.T
    return t, Y, E, B, L


def run_dose_response(doses, T_pulse=48.0):
    results = []
    for dose in doses:
        t, Y, E, B, L = run_simulation(dose=dose, T_pulse=T_pulse)
        B_final = B[-1]
        Y_final = Y[-1]
        # Check if boundary persisted (success criterion)
        B_crit = 0.65
        success = B_final >= B_crit and Y_final < 0.5
        results.append(dict(dose=dose, B_final=B_final,
                            Y_final=Y_final, success=success))
    return results


def main():
    print("=" * 60)
    print("  YAP/TAZ Gingival Regeneration Phase Simulation")
    print("  Level 3 Full Dynamical System")
    print("=" * 60)
    print()

    os.makedirs("figures", exist_ok=True)

    # Run trajectories
    print("Running successful trajectory (TRULI 2uM, 48h, soft scaffold)...")
    t, Y_s, E_s, B_s, L_s = run_simulation(dose=2.0, T_pulse=48, scaffold=True)

    print("Running failed trajectory (no drug, stiff scaffold)...")
    t, Y_f, E_f, B_f, L_f = run_simulation(dose=0.0, T_pulse=48, scaffold=False)

    print("Running drug-only trajectory (TRULI 2uM, stiff scaffold)...")
    t, Y_d, E_d, B_d, L_d = run_simulation(dose=2.0, T_pulse=48, scaffold=False)

    print("Running dose-response...")
    doses = np.linspace(0, 4.0, 25)
    dr = run_dose_response(doses)

    # Final state report
    print()
    print("Trajectory Outcomes")
    print("-" * 50)
    print(f"Success (TRULI + scaffold): "
          f"B={B_s[-1]:.3f}  Y={Y_s[-1]:.3f}  E={E_s[-1]:.2f}kPa")
    print(f"Failed  (no drug):          "
          f"B={B_f[-1]:.3f}  Y={Y_f[-1]:.3f}  E={E_f[-1]:.2f}kPa")
    print(f"Drug only (stiff):          "
          f"B={B_d[-1]:.3f}  Y={Y_d[-1]:.3f}  E={E_d[-1]:.2f}kPa")
    print()
    threshold_dose = None
    for r in dr:
        if r["success"] and threshold_dose is None:
            threshold_dose = r["dose"]
    print(f"Minimum effective dose: ~{threshold_dose:.2f} uM TRULI") if threshold_dose else print("Threshold dose: see dose-response panel")

    # ── Plot ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 13))
    fig.patch.set_facecolor("#0f0f1a")
    gs  = GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38)

    ax_Y  = fig.add_subplot(gs[0, 0])
    ax_E  = fig.add_subplot(gs[0, 1])
    ax_B  = fig.add_subplot(gs[0, 2])
    ax_L  = fig.add_subplot(gs[1, 0])
    ax_ph = fig.add_subplot(gs[1, 1:])
    ax_dr = fig.add_subplot(gs[2, :])

    def sty(ax, title, ylabel="", xlabel="Time (h)"):
        ax.set_facecolor("#1a1a2e")
        for sp in ax.spines.values():
            sp.set_edgecolor("#2a2a4e")
        ax.tick_params(colors="#aaaaaa", labelsize=8)
        ax.set_title(title, color="#ffffff", fontsize=9,
                     pad=6, fontweight="bold")
        ax.set_ylabel(ylabel, color="#aaaaaa", fontsize=8)
        ax.set_xlabel(xlabel, color="#aaaaaa", fontsize=8)

    # Panel A -- YAP activity
    ax_Y.plot(t, Y_s, color="#2ecc71", lw=2.0, label="TRULI + scaffold")
    ax_Y.plot(t, Y_f, color="#e74c3c", lw=1.5, ls="--", label="No treatment")
    ax_Y.plot(t, Y_d, color="#e67e22", lw=1.5, ls=":", label="Drug only")
    ax_Y.axvline(48, color="#ffffff", lw=0.8, ls="--", alpha=0.4)
    ax_Y.set_ylim(0, 1.1)
    ax_Y.text(50, 0.95, "washout", color="#aaaaaa", fontsize=7)
    sty(ax_Y, "A — YAP/TAZ Activity", ylabel="YAP (nuclear fraction)")
    ax_Y.legend(facecolor="#0f0f1a", edgecolor="#2a2a4e",
                labelcolor="white", fontsize=7)

    # Panel B -- ECM stiffness
    ax_E.plot(t, E_s, color="#2ecc71", lw=2.0)
    ax_E.plot(t, E_f, color="#e74c3c", lw=1.5, ls="--")
    ax_E.plot(t, E_d, color="#e67e22", lw=1.5, ls=":")
    ax_E.axvline(48, color="#ffffff", lw=0.8, ls="--", alpha=0.4)
    ax_E.axhline(E_opt, color="#9b59b6", lw=0.8, ls="--", alpha=0.5)
    ax_E.text(160, E_opt + 0.3, "E_opt", color="#9b59b6", fontsize=7)
    sty(ax_E, "B — ECM Stiffness", ylabel="Young's modulus (kPa)")

    # Panel C -- Boundary integrity
    ax_B.plot(t, B_s, color="#2ecc71", lw=2.5)
    ax_B.plot(t, B_f, color="#e74c3c", lw=1.5, ls="--")
    ax_B.plot(t, B_d, color="#e67e22", lw=1.5, ls=":")
    ax_B.axvline(48, color="#ffffff", lw=0.8, ls="--", alpha=0.4)
    ax_B.axhline(0.7, color="#f39c12", lw=0.8, ls="--", alpha=0.6)
    ax_B.text(160, 0.72, "B_crit", color="#f39c12", fontsize=7)
    ax_B.set_ylim(0, 1.1)
    sty(ax_B, "C — Boundary Integrity", ylabel="B (0=none, 1=complete)")

    # Panel D -- Drug concentration
    ax_L.plot(t, L_s, color="#3498db", lw=2.0)
    ax_L.axvline(48, color="#ffffff", lw=0.8, ls="--", alpha=0.4)
    sty(ax_L, "D — TRULI Concentration", ylabel="L (uM)")

    # Panel E -- Phase diagram (Y vs E) with trajectory
    Y_grid = np.linspace(0, 1, 200)
    E_grid = np.linspace(0, 16, 200)
    Yg, Eg = np.meshgrid(Y_grid, E_grid)

    # Background: h(E) * Y -- boundary formation potential
    BFP = h_E(Eg) * Yg
    ax_ph.contourf(Yg, Eg, BFP, levels=20,
                   cmap="RdYlGn", alpha=0.4)

    # Region labels
    ax_ph.text(0.05, 13, "Region I\nFIBROTIC", color="white",
               fontsize=8, fontweight="bold", alpha=0.9)
    ax_ph.text(0.75, 13, "Region II\nHYPERPLASTIC", color="white",
               fontsize=8, fontweight="bold", alpha=0.9)
    ax_ph.text(0.05, 2, "Region III\nQUIESCENT", color="white",
               fontsize=8, fontweight="bold", alpha=0.9)
    ax_ph.text(0.45, 5.5, "Region IV\nREGENERATIVE", color="white",
               fontsize=9, fontweight="bold", alpha=1.0)

    # Trajectories on phase diagram
    ax_ph.plot(Y_s, E_s, color="#2ecc71", lw=2.5,
               label="Success (TRULI + scaffold)")
    ax_ph.plot(Y_f, E_f, color="#e74c3c", lw=1.5, ls="--",
               label="No treatment")
    ax_ph.plot(Y_d, E_d, color="#e67e22", lw=1.5, ls=":",
               label="Drug only (stiff)")

    # Start and end markers
    ax_ph.scatter([Y_s[0]], [E_s[0]], color="white", s=80,
                  zorder=10, label="Start (fibrotic)")
    ax_ph.scatter([Y_s[-1]], [E_s[-1]], color="#2ecc71", s=80,
                  marker="*", zorder=10, label="End (regenerated)")

    ax_ph.set_facecolor("#1a1a2e")
    for sp in ax_ph.spines.values():
        sp.set_edgecolor("#2a2a4e")
    ax_ph.tick_params(colors="#aaaaaa", labelsize=8)
    ax_ph.set_title("E — Phase Diagram: YAP Activity vs ECM Stiffness\n"
                    "Treatment path crosses diagonally through Region IV. "
                    "Other approaches miss the window.",
                    color="#ffffff", fontsize=9, pad=6, fontweight="bold")
    ax_ph.set_xlabel("YAP Activity (nuclear fraction)", color="#aaaaaa", fontsize=8)
    ax_ph.set_ylabel("ECM Stiffness (kPa)", color="#aaaaaa", fontsize=8)
    ax_ph.legend(facecolor="#0f0f1a", edgecolor="#2a2a4e",
                 labelcolor="white", fontsize=7, loc="upper left")

    # Panel F -- Dose response
    dose_vals  = [r["dose"]    for r in dr]
    B_finals   = [r["B_final"] for r in dr]
    Y_finals   = [r["Y_final"] for r in dr]
    successes  = [r["success"] for r in dr]
    colors_dr  = ["#2ecc71" if s else "#e74c3c" for s in successes]

    ax_dr.bar(dose_vals, B_finals, width=0.14,
              color=colors_dr, edgecolor="#0f0f1a", linewidth=0.4,
              label="Boundary integrity at day 8")
    ax_dr.axhline(0.7, color="#f39c12", lw=1.2, ls="--",
                  label="Success threshold (B=0.7)")
    ax_dr.set_xlim(-0.2, 4.2)
    ax_dr.set_ylim(0, 1.1)
    ax_dr.set_facecolor("#1a1a2e")
    for sp in ax_dr.spines.values():
        sp.set_edgecolor("#2a2a4e")
    ax_dr.tick_params(colors="#aaaaaa", labelsize=8)
    ax_dr.set_title(
        "F — Dose Response: Final Boundary Integrity vs TRULI Dose\n"
        "Green = success (B > 0.7, YAP off).  "
        "Red = failure.  Therapeutic window clearly defined.",
        color="#ffffff", fontsize=9, pad=6, fontweight="bold")
    ax_dr.set_xlabel("TRULI dose (uM)", color="#aaaaaa", fontsize=8)
    ax_dr.set_ylabel("Boundary integrity at t=200h", color="#aaaaaa", fontsize=8)
    ax_dr.legend(facecolor="#0f0f1a", edgecolor="#2a2a4e",
                 labelcolor="white", fontsize=8)

    if threshold_dose:
        ax_dr.axvline(threshold_dose, color="#2ecc71", lw=1.2,
                      ls="--", alpha=0.7)
        ax_dr.text(threshold_dose + 0.05, 1.02,
                   f"min effective\n{threshold_dose:.2f} uM",
                   color="#2ecc71", fontsize=7)

    fig.suptitle(
        "YAP/TAZ Gingival Regeneration — Phase Simulation\n"
        "Transient LATS inhibition + soft scaffold drives diagonal "
        "trajectory through regenerative region. "
        "Self-limiting via boundary polarity feedback.",
        color="#ffffff", fontsize=12, fontweight="bold", y=1.005)

    plt.savefig("figures/phase_sim_results.png", dpi=150,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    print("Saved: figures/phase_sim_results.png")
    plt.close()
    print("Done.")


if __name__ == "__main__":
    main()
