"""
layer3_kinetics.py
==================
Layer 3 — Ionic crosslinking kinetics (Ca-alginate gelation).

Uses an Arrhenius rate model to compute the crosslinking rate constant,
then integrates a first-order conversion equation to track gel fraction
over time during the printing / bath immersion process.
"""

import numpy as np
from scipy.integrate import solve_ivp
from parameters import Ea_cross, A_cross, R_gas, T_process, T_bath


def crosslink_rate_constant(T: float = T_process) -> float:
    """
    Arrhenius crosslinking rate constant [1/s].

    k(T) = A_cross * exp(-Ea_cross / (R_gas * T))
    """
    return A_cross * np.exp(-Ea_cross / (R_gas * T))


def gel_fraction_ode(t, y, k):
    """ODE: d(alpha)/dt = k * (1 - alpha),  alpha = gel conversion fraction."""
    alpha = y[0]
    return [k * (1.0 - alpha)]


def crosslinking_profile(t_span: tuple = (0, 300), T: float = T_bath, n_points: int = 300):
    """
    Integrate gel conversion from t=0 to t_span[1] seconds at temperature T.

    Returns
    -------
    t      : time array [s]
    alpha  : gel conversion fraction (0 = sol, 1 = fully gelled)
    k      : rate constant used [1/s]
    """
    k = crosslink_rate_constant(T)
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(gel_fraction_ode, t_span, [0.0], args=(k,), t_eval=t_eval, dense_output=True)
    return sol.t, sol.y[0], k


def run(T: float = T_bath, t_end: float = 300.0):
    """Return crosslinking profile dict for downstream layers."""
    t, alpha, k = crosslinking_profile(t_span=(0, t_end), T=T)
    return {"t": t, "alpha": alpha, "k_cross": k, "alpha_final": alpha[-1]}


def plot(save_path: str = "figures/layer3_kinetics.png"):
    """
    Gel conversion alpha(t) at bath temperature (37 C).
    t_gel (where alpha = 0.9) is marked with a red dashed line.
    """
    import matplotlib.pyplot as plt

    result = run(t_end=120.0)
    t     = result["t"]
    alpha = result["alpha"]
    k     = result["k_cross"]
    t_gel = -np.log(1.0 - 0.9) / k

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t, alpha, "b-", linewidth=2.5, label="Gel conversion alpha(t)")
    ax.axhline(0.9, color="gray",  linestyle=":",  linewidth=1.5,
               label="alpha = 0.9  (gel threshold)")
    ax.axvline(t_gel, color="red", linestyle="--", linewidth=1.5,
               label=f"t_gel = {t_gel:.2f} s")
    ax.plot(t_gel, 0.9, "r*", markersize=14, zorder=5)

    ax.set_xlabel("Time [s]", fontsize=12)
    ax.set_ylabel("Gel Conversion Fraction alpha", fontsize=12)
    ax.set_title(
        f"Layer 3 — Ca-Alginate Crosslinking Kinetics\n"
        f"Arrhenius model  |  T = {T_bath - 273.15:.0f} C  |  "
        f"k = {k:.3e} 1/s",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, t[-1])
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    result = run()
    print(f"Rate constant k = {result['k_cross']:.4e} 1/s")
    print(f"Final gel conversion after {result['t'][-1]:.0f} s: alpha = {result['alpha_final']:.4f}")
    plot()
