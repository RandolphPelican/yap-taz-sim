"""
YAP/TAZ Bistable Mechanosensitive Switch — Core ODE Model
Author: John D. Stabler (Randolph Pelican III)
GitHub: https://github.com/RandolphPelican/yap-taz-sim
License: MIT
"""

import numpy as np
from scipy.integrate import odeint

# ── Default parameters ──────────────────────────────────────────────────────
PARAMS = {
    "alpha": 2.0,   # activation rate (h⁻¹)
    "beta":  0.5,   # deactivation rate (h⁻¹)
    "S0":    5.0,   # optimal stiffness (kPa)
    "sigma": 2.0,   # stiffness window width (kPa)
    "theta": 0.3,   # Hill threshold (unitless)
    "n":     4,     # Hill cooperativity
}


def f_stiffness(S, S0, sigma):
    """Bell-shaped stiffness activation function."""
    return np.exp(-((S - S0) ** 2) / (2 * sigma ** 2))


def sigma_feedback(R, theta, n):
    """Hill-type positive feedback (YAP-actin-ECM loop)."""
    Rn = R ** n
    return Rn / (theta ** n + Rn)


def dRdt(R, t, S, p, feedback=True, L=0.0, K_L=0.02):
    """
    ODE for regenerative fraction R.

    Parameters
    ----------
    R       : float  regenerative fraction [0,1]
    t       : float  time (h)
    S       : float  ECM stiffness (kPa)
    p       : dict   parameter dict
    feedback: bool   if False, disables Hill feedback (Sim 1 knockout)
    L       : float  LATS inhibitor concentration (μM)
    K_L     : float  drug IC₅₀ (μM)
    """
    fS = f_stiffness(S, p["S0"], p["sigma"])
    fb = sigma_feedback(R, p["theta"], p["n"]) if feedback else 1.0
    g  = K_L / (L + K_L)           # drug modulation of deactivation
    activation   = p["alpha"] * fS * fb * (1 - R)
    deactivation = p["beta"]  * g  * (1 - fS) * R
    return activation - deactivation


def simulate(S, R0=0.01, t_max=500, n_steps=5000, p=None, **kw):
    """Integrate ODE; returns (t, R)."""
    p = p or PARAMS
    t = np.linspace(0, t_max, n_steps)
    R = odeint(dRdt, R0, t, args=(S, p), **kw)
    return t, R.flatten()


def steady_states(S_range=None, p=None):
    """
    Sweep stiffness; return (S, R_low, R_high) for bistability diagram.
    """
    p = p or PARAMS
    S_range = S_range if S_range is not None else np.linspace(0, 10, 200)
    R_low, R_high = [], []
    for S in S_range:
        _, r_lo = simulate(S, R0=0.01, p=p)
        _, r_hi = simulate(S, R0=0.99, p=p)
        R_low.append(r_lo[-1])
        R_high.append(r_hi[-1])
    return S_range, np.array(R_low), np.array(R_high)
#..
