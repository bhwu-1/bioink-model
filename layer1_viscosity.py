"""
layer1_viscosity.py
===================
Layer 1 — Apparent viscosity of the fibrin–alginate bioink.

Uses the Power Law model:  eta(gamma_dot) = K * gamma_dot^(n-1)
K and n are bivariate functions of both alginate and fibrinogen concentration,
evaluated via K_func() and n_func() from parameters.py.

Wall shear rate (Rabinowitsch correction for power-law fluid in a tube):
    gamma_dot_w = (3n+1)/(4n) * 4Q / (pi * R^3)
"""

import numpy as np
from parameters import K_func, n_func, R_channel, Q_print, c_alg_nom, c_fib_nom


def apparent_viscosity(gamma_dot: float, K: float, n: float) -> float:
    """Return apparent viscosity [Pa·s] at a given shear rate [1/s]."""
    return K * gamma_dot ** (n - 1)


def wall_shear_rate(Q: float = Q_print, R: float = R_channel, n: float = 1.0) -> float:
    """
    Corrected (Rabinowitsch) wall shear rate [1/s] for a power-law fluid
    in a cylindrical tube.
    """
    gamma_apparent = 4 * Q / (np.pi * R**3)
    correction = (3 * n + 1) / (4 * n)
    return correction * gamma_apparent


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom,
        Q: float = Q_print, R: float = R_channel) -> dict:
    """
    Compute wall shear rate and apparent viscosity for the fibrin–alginate bioink.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    Q     : volumetric flow rate [m³/s]
    R     : nozzle radius [m]
    """
    K = K_func(c_alg, c_fib)
    n = n_func(c_alg, c_fib)
    gamma_w = wall_shear_rate(Q, R, n)
    eta_w = apparent_viscosity(gamma_w, K, n)
    return {"gamma_wall": gamma_w, "eta_wall": eta_w, "K": K, "n": n}


def plot(save_path: str = "figures/layer1_viscosity.png"):
    """
    Viscosity vs shear rate curves for four representative formulations.
    A vertical dashed line marks the operating shear rate inside the RX1 nozzle.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25% alg / 10 mg/mL fib"),
        (0.50, 20,  "0.50% alg / 20 mg/mL fib  (nominal)"),
        (1.00, 20,  "1.00% alg / 20 mg/mL fib"),
        (1.50, 30,  "1.50% alg / 30 mg/mL fib"),
    ]
    colors = ["#a8d1e7", "#5ba3c9", "#2171b5", "#08306b"]
    gamma_range = np.logspace(0, 3, 300)   # 1 – 1000 1/s

    fig, ax = plt.subplots(figsize=(8, 5))

    for (c_alg, c_fib, label), color in zip(test_points, colors):
        K = K_func(c_alg, c_fib)
        n = n_func(c_alg, c_fib)
        eta = apparent_viscosity(gamma_range, K, n)
        ax.loglog(gamma_range, eta, color=color, linewidth=2, label=label)

    # Mark operating shear rate at nominal formulation
    r_nom = run(c_alg_nom, c_fib_nom)
    gamma_op = r_nom["gamma_wall"]
    eta_op   = r_nom["eta_wall"]
    ax.axvline(gamma_op, color="k", linestyle="--", linewidth=1.2, alpha=0.7,
               label=f"RX1 nozzle  ({gamma_op:.0f} 1/s)")
    ax.plot(gamma_op, eta_op, "k*", markersize=13, zorder=5)

    ax.set_xlabel("Shear Rate [1/s]", fontsize=12, fontweight="bold")
    ax.set_ylabel("Apparent Viscosity [Pa·s]", fontsize=12, fontweight="bold")
    ax.set_title("Layer 1: Bioink Viscosity\nPower Law model, fibrin–alginate", fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, labelsize=11)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    test_points = [(0.25, 10), (0.5, 20), (1.0, 20), (1.5, 30)]
    for c_alg, c_fib in test_points:
        r = run(c_alg, c_fib)
        print(f"c_alg={c_alg} wt%, c_fib={c_fib:4.0f} mg/mL | "
              f"K={r['K']:.4f} Pa·s^n  n={r['n']:.3f}  "
              f"gamma_w={r['gamma_wall']:.1f} 1/s  eta_w={r['eta_wall']:.4f} Pa·s")
    plot()
