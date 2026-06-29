"""
layer7_ecm.py
=============
Layer 7 — Neural ECM deposition and dynamic stiffness evolution.

Models post-printing accumulation of cell-secreted ECM (fibronectin, laminin,
tenascin-C) over a 30-day culture period and its effect on effective gel
stiffness experienced by embedded GBM / neural cells.

Design intent:
  G'_initial ≈ 30 Pa  (freshly printed nominal formulation)
  G'_eff(30d) ≈ 278 Pa → approaching G'_opt = 300 Pa as ECM accumulates
"""

import numpy as np
from scipy.integrate import solve_ivp
from parameters import (k_dep, alpha_ecm, G_prime_func, G_prime_opt,
                        t_culture, c_alg_nom, c_fib_nom)
from layer6_stiffness import stiffness_viability


def ecm_ode(t, y, k):
    """
    First-order ECM deposition ODE.

    d(rho_ecm)/dt = k_dep * (1 - rho_ecm)

    rho_ecm ∈ [0, 1]: normalised ECM density (0 = none, 1 = fully deposited).
    Half-life ≈ 8 days at k_dep = 1e-6 1/s.
    """
    return [k * (1.0 - y[0])]


def effective_modulus(G_gel: float, rho_ecm: float,
                      alpha: float = alpha_ecm) -> float:
    """
    Effective substrate modulus including ECM contribution [Pa].

    G'_eff = G_gel * (1 + alpha_ecm * rho_ecm)
    """
    return G_gel * (1.0 + alpha * rho_ecm)


def ecm_evolution(
    c_alg: float = c_alg_nom,
    c_fib: float = c_fib_nom,
    gel_alpha: float = 1.0,
    t_span: tuple = None,
    n_points: int = 300,
) -> tuple:
    """
    Simulate ECM deposition and effective stiffness over the 30-day culture.

    Parameters
    ----------
    c_alg     : alginate concentration [wt%]
    c_fib     : fibrinogen concentration [mg/mL]
    gel_alpha : initial gel conversion fraction (0–1)
    t_span    : (t_start, t_end) [s]; defaults to (0, t_culture = 30 days)
    n_points  : number of output time points

    Returns
    -------
    t         : time array [s]
    rho_ecm   : normalised ECM density
    G_eff     : effective modulus [Pa]
    viability : stiffness-based viability score (0–1)
    """
    if t_span is None:
        t_span = (0, t_culture)

    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(ecm_ode, t_span, [0.0], args=(k_dep,), t_eval=t_eval,
                    method="RK45", rtol=1e-6)

    rho = sol.y[0]
    G_gel = G_prime_func(c_alg, c_fib, gel_alpha)
    G_eff = effective_modulus(G_gel, rho)
    viability = np.vectorize(stiffness_viability)(G_eff)

    return sol.t, rho, G_eff, viability


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom) -> dict:
    """Return ECM evolution results dict over the full 30-day culture period."""
    t, rho, G_eff, viability = ecm_evolution(c_alg=c_alg, c_fib=c_fib)
    return {
        "t": t,
        "rho_ecm": rho,
        "G_eff": G_eff,
        "stiffness_viability": viability,
        "G_eff_final": float(G_eff[-1]),
        "viability_final": float(viability[-1]),
    }


def plot(save_path: str = "figures/layer7_ecm.png"):
    """
    Left panel  — G'_eff(t) over 30-day culture for all test formulations.
                  Dashed red line marks G'_opt = 300 Pa.
    Right panel — Stiffness-based viability score over the same period.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25% alg / 10 mg/mL fib"),
        (0.50, 20,  "0.50% alg / 20 mg/mL fib  (nominal)"),
        (1.00, 20,  "1.00% alg / 20 mg/mL fib"),
        (1.50, 30,  "1.50% alg / 30 mg/mL fib"),
    ]
    colors = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for (c_alg, c_fib, label), color in zip(test_points, colors):
        r = run(c_alg, c_fib)
        t_days = r["t"] / 86400
        ax1.plot(t_days, r["G_eff"],              color=color, linewidth=2, label=label)
        ax2.plot(t_days, r["stiffness_viability"], color=color, linewidth=2, label=label)

    ax1.axhline(G_prime_opt, color="red", linestyle="--", linewidth=1.5,
                label=f"G_opt = {G_prime_opt:.0f} Pa  (brain tissue)")
    ax1.set_xlabel("Culture Time [days]", fontsize=11)
    ax1.set_ylabel("Effective Modulus G'_eff [Pa]", fontsize=11)
    ax1.set_title("ECM-Driven Stiffness Evolution\n(30-day culture)", fontsize=11)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.axhline(1.0, color="gray", linestyle=":", linewidth=1, alpha=0.5)
    ax2.set_xlabel("Culture Time [days]", fontsize=11)
    ax2.set_ylabel("Stiffness Viability Score", fontsize=11)
    ax2.set_title("Mechanosensing Score Over 30 Days", fontsize=11)
    ax2.legend(fontsize=8)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Layer 7 — Neural ECM Deposition & Stiffness Evolution", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    test_points = [(0.25, 10), (0.5, 20), (1.0, 20), (1.5, 30)]
    for c_alg, c_fib in test_points:
        r = run(c_alg, c_fib)
        print(f"c_alg={c_alg} wt%, c_fib={c_fib:4.0f} mg/mL | "
              f"G'_eff(30d)={r['G_eff_final']:.1f} Pa  "
              f"viability(30d)={r['viability_final']:.3f}")
    plot()
