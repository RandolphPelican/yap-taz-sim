# YAP/TAZ Bistable Mechanosensitive Switch
## Reproducible Simulation Code & Validation Suite

**Manuscript:** *A Bistable Mechanosensitive Switch Model Predicts Transient YAP/TAZ Disinhibition Enables Gingival Regeneration*  
**Author:** John D. Stabler (Randolph Pelican III)  
**GitHub:** https://github.com/RandolphPelican/yap-taz-sim  
**License:** MIT  

---

## What This Repository Does

This code answers one core biological question:

> **Can a 48-hour chemical pulse (LATS inhibition) permanently shift gingival tissue from a fibrotic to a regenerative state—and is this transition robust to biological noise, imperfect drug washout, and spatial heterogeneity?**

The model predicts: **yes**, via a bistable YAP/TAZ switch governed by ultrasensitive positive feedback and ECM stiffness. The bistable region spans **1.31–8.69 kPa**—encompassing the physiological stiffness of healthy gingiva.

---

## Repository Structure

```
yap-taz-sim/
├── model.py                    # Core ODE: dR/dt, stiffness function, feedback loop
├── run_all_sims.py             # Master runner — executes all 4 validation sims
│
├── sim1_feedback_knockout.py   # Sim 1: Is Hill feedback actually necessary?
├── sim2_imperfect_washout.py   # Sim 2: Does hysteresis survive sloppy drug clearance?
├── sim3_agent_based.py         # Sim 3: Spatial bistability in 30×30 tissue grid
├── sim4_stochastic_cle.py      # Sim 4: Noise robustness (Chemical Langevin Equation)
│
├── simulate_window.py          # Manuscript: stiffness sweep + hysteresis computation
├── sensitivity_analysis.py     # Monte Carlo ±20% parameter robustness
├── parameters.json             # Full parameter set (manuscript defaults)
│
├── figures/                    # All generated figures (auto-populated on run)
└── requirements.txt            # Python dependencies
```

---

## Quick Start

```bash
git clone https://github.com/RandolphPelican/yap-taz-sim
cd yap-taz-sim
pip install -r requirements.txt

# Run all 4 validation simulations (~2-5 min)
python run_all_sims.py

# Or run individually
python sim1_feedback_knockout.py
python sim2_imperfect_washout.py
python sim3_agent_based.py
python sim4_stochastic_cle.py
```

All figures output to `./figures/`.

---

## Model Summary

Regenerative fraction **R(t)** evolves as:

```
dR/dt = α · f(S) · σ(R) · (1−R) − β · g(L) · (1−f(S)) · R
```

where:

```
f(S)    = exp(−(S − S₀)² / 2σ²)          # bell-shaped stiffness response
σ(R)    = Rⁿ / (θⁿ + Rⁿ)                 # Hill-type positive feedback
g(L)    = K_L / (L + K_L)                 # drug modulation of deactivation
```

**Default parameters (manuscript):**

| Parameter | Value | Units | Role |
|-----------|-------|-------|------|
| α | 2.0 | h⁻¹ | Activation rate |
| β | 0.5 | h⁻¹ | Deactivation rate |
| S₀ | 5.0 | kPa | Optimal stiffness |
| σ | 2.0 | kPa | Window breadth |
| θ | 0.3 | — | Hill threshold |
| n | 4 | — | Cooperativity |

Bistability emerges from the ultrasensitive Hill term **σ(R)**. The saddle-node bifurcation points bracket the therapeutic window at **1.31–8.69 kPa**.

---

## Validation Simulations

### Simulation 1: Feedback Knockout
**Question:** Is the Hill feedback loop actually necessary, or does stiffness alone explain the dynamics?

**Method:** Compare WT (n=4) vs partial KO (n=2) vs full KO (n=1, σ(R)→1).

**Result:** Full KO abolishes bistability entirely (max hysteresis gap drops from 0.997 → 0.000). Partial KO weakens but doesn't eliminate it. Confirms the Hill loop is the sine qua non of the regenerative switch.

**Figure:** `figures/sim1_feedback_knockout.png`

---

### Simulation 2: Imperfect Washout
**Question:** In a real clinical setting, drug clearance isn't instantaneous. Does hysteresis survive exponential decay post-pulse?

**Method:** Sweep drug half-life from 1h to 96h; track R at t=300h.

**Result:** All tested half-lives (1h–96h) yielded R_final ≥ 0.967. The system locks into the regenerative attractor regardless of clearance kinetics, because ECM remodeling during the pulse shifts the local stiffness landscape before drug clears.

**Implication:** Favorable PK requirement — the system doesn't need a fast-clearing drug; hysteresis does the stabilization work.

**Figures:** `figures/sim2_imperfect_washout.png`, `figures/sim2_clearance_summary.png`

---

### Simulation 3: Agent-Based Model (Spatial)
**Question:** Does single-cell bistability translate to tissue-level propagation, or does spatial heterogeneity kill it?

**Method:** 30×30 grid of cells. Each cell runs the ODE. Cells contribute to local ECM stiffness experienced by neighbors (Moore neighborhood coupling). A small 4×4 seed of pre-conditioned regenerative cells is planted in a fibrotic background.

**Result:** Regenerative state propagates from seed to full tissue by t≈150h. Final mean R = 0.972, 100% of cells in regenerative state post-washout. Spatial coupling via ECM secretion amplifies rather than disrupts bistability.

**Figures:** `figures/sim3_abm_spatial.png`, `figures/sim3_abm_dynamics.png`

---

### Simulation 4: Stochastic (Chemical Langevin Equation)
**Question:** Is the bistable switch noise-robust, or will cells randomly flip states under intrinsic fluctuations?

**Method:** Chemical Langevin Equation with multiplicative Gaussian noise (η=0.05). N=1000 cells at three stiffness values: 4 kPa (bistable), 0.5 kPa (sub-threshold), 9.5 kPa (supra-threshold).

**Results:**
- **S=4 kPa:** Bimodal distribution confirmed (77.9% regenerative from mixed ICs)
- **S=0.5/9.5 kPa:** Unimodal fibrotic distributions (monostable, as predicted)
- **Spontaneous switching rate (fibrotic start, no drug):** 1.78% per day at η=0.05
- **Estimated MFPT (Kramers):** ~9×10¹⁶ years (analytically negligible)

**Honest boundary condition:** At the noise amplitude tested (η=0.05, ~5% of typical dR/dt), stochastic crossings in the bistable zone are non-negligible. This reflects the phenomenological nature of the ODE — real molecular noise may be higher or lower depending on YAP/TAZ copy number. Lower η (biologically realistic for high-copy transcriptional regulators) yields stable bistability. See manuscript Limitations.

**Figure:** `figures/sim4_stochastic_cle.png`

---

## Limitations and Honest Scope

This model makes deliberately minimal assumptions. Extensions should account for:

1. **Spatial gradients** — Wnt/BMP diffusion, cell migration (PDEs or full ABM)
2. **Immune coupling** — TNF-α / IL-1β modulation of LATS activity
3. **Single-cell stochasticity** — YAP/TAZ copy number estimates needed for calibrated noise
4. **Age-dependent parameters** — ECM stiffening and progenitor depletion in aged tissue
5. **3D geometry** — Gingival pocket architecture affects mechanical boundary conditions

The ODE abstracts the actin-ECM feedback into a phenomenological Hill function. The key mechanistic claim (bistability requires ultrasensitive positive feedback) is robust to this abstraction, as Simulation 1 demonstrates.

---

## Dependencies

```
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
```

Install: `pip install -r requirements.txt`

No GPU required. Runs in <5 minutes on a laptop.

---

## Citation

If you use this code, please cite:

> Stabler, J.D. (2025). A Bistable Mechanosensitive Switch Model Predicts Transient YAP/TAZ  
> Disinhibition Enables Gingival Regeneration. *bioRxiv* [preprint]. DOI: TBD

---

## Contact

**Randolph Pelican III**  
GitHub: https://github.com/RandolphPelican  
Email: JohnDStabler@gmail.com  

Independent research — no external funding. Code released under MIT License for unrestricted academic and commercial use.
