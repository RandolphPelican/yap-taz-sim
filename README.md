# YAP/TAZ Bistable Mechanosensitive Switch: Gingival Regeneration Model

**A mathematical framework for understanding and targeting tissue state transitions in periodontal repair**

## Overview

This repository contains a **bistable ODE model** demonstrating that adult gingival tissue is trapped in a fibrotic state via dysregulated YAP/TAZ mechanotransduction.

**Key Finding:** YAP/TAZ acts as a stiffness-gated bistable switch with hysteresis, creating a therapeutic window at **1.31–8.69 kPa**.

## Quick Start

```bash
pip install numpy scipy matplotlib
cd scripts
python3 bistable_yap_model.py
Key Results
Bistable region: 1.31–8.69 kPa (matches gingival tissue)
Hysteresis confirmed (upper/lower branches separated)
Sharp state transitions at boundary stiffness values
Repository Structure
yap-taz-sim/
├── scripts/           # Simulation code
├── figures/           # Generated plots
├── docs/              # Documentation & concept notes
└── README.md
Mathematical Model
dR/dt = α·f(S)·[R^n/(θ^n + R^n)]·(1-R) - β·(1-f(S))·R
Parameters grounded in:
Dupont et al., Nature 2011 (YAP mechanotransduction)
Elosegui-Artola et al., Cell 2017 (YAP-actin feedback)
License
MIT License - see LICENSE file
Contact
John D. Stabler
GitHub: @RandolphPelican
