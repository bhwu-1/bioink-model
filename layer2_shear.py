"""
layer2_shear.py
===============
Layer 2 — Cell damage assessment from wall shear stress.

Computes wall shear stress and evaluates whether it exceeds tau_crit (100 Pa
for fragile neural / GBM cells). A sigmoidal damage function smoothly models
viability loss above the threshold.
"""

import numpy as np
from parameters import tau_crit, K_func, n_func, R_channel, Q_print, c_alg_nom, c_fib_nom
from layer1_viscosity import wall_shear_rate, apparent_viscosity


def wall_shear_stress(K: float, n: float,
                      Q: float = Q_print, R: float = R_channel) -> float:
    """Return wall shear stress [Pa] for a power-law fluid in a cylindrical tube."""
    gamma_w = wall_shear_rate(Q, R, n)
    eta_w = apparent_viscosity(gamma_w, K, n)
    return eta_w * gamma_w


def cell_viability_shear(tau: float, tau_c: float = tau_crit,
                         k: float = 0.05) -> float:
    """
    Sigmoidal cell viability factor (0–1) as a function of wall shear stress.

    viability = 1 / (1 + exp(k * (tau - tau_c)))

    Parameters
    ----------
    tau   : wall shear stress [Pa]
    tau_c : critical threshold [Pa]  (100 Pa for neural cells)
    k     : steepness [1/Pa]  — tune to experimental viability data
    """
    return 1.0 / (1.0 + np.exp(k * (tau - tau_c)))


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom,
        Q: float = Q_print, R: float = R_channel) -> dict:
    """
    Compute shear stress and viability factor for the fibrin–alginate bioink.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    Q     : volumetric flow rate [m³/s]
    R     : nozzle radius [m]
    """
    K = K_func(c_alg, c_fib)
    n = n_func(c_alg, c_fib)
    tau_w = wall_shear_stress(K, n, Q, R)
    viability = cell_viability_shear(tau_w)
    return {"tau_wall": tau_w, "viability_shear": viability}


def plot(save_path: str = "figures/layer2_shear.png"):
    """
    Left panel  — bar chart of wall shear stress per formulation vs tau_crit.
    Right panel — sigmoidal viability function with formulations marked.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25%\n10 mg/mL"),
        (0.50, 20,  "0.50%\n20 mg/mL\n(nominal)"),
        (1.00, 20,  "1.00%\n20 mg/mL"),
        (1.50, 30,  "1.50%\n30 mg/mL"),
    ]
    colors = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c"]

    results = [run(c_alg, c_fib) for c_alg, c_fib, _ in test_points]
    taus      = [r["tau_wall"]        for r in results]
    viabs     = [r["viability_shear"] for r in results]
    bar_colors = ["#4daf4a" if t < tau_crit else "#e41a1c" for t in taus]
    labels    = [p[2] for p in test_points]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: bar chart ---
    bars = ax1.bar(labels, taus, color=bar_colors, edgecolor="k", linewidth=0.8)
    ax1.axhline(tau_crit, color="red", linestyle="--", linewidth=1.5,
                label=f"tau_crit = {tau_crit:.0f} Pa  (neural cells)")
    for bar, tau in zip(bars, taus):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 f"{tau:.1f}", ha="center", va="bottom", fontsize=9)
    ax1.set_ylabel("Wall Shear Stress [Pa]", fontsize=11)
    ax1.set_title("Wall Shear Stress by Formulation", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.set_ylim(0, max(taus) * 1.35)

    # --- Right: sigmoidal viability curve ---
    tau_range = np.linspace(0, 200, 400)
    viab_curve = np.array([cell_viability_shear(t) for t in tau_range])
    ax2.plot(tau_range, viab_curve, "b-", linewidth=2.5)
    ax2.axvline(tau_crit, color="red", linestyle="--", linewidth=1.5,
                label=f"tau_crit = {tau_crit:.0f} Pa")
    for (c_alg, c_fib, _), color, tau, viab in zip(test_points, colors, taus, viabs):
        ax2.plot(tau, viab, "o", color=color, markersize=9, zorder=5)
    ax2.set_xlabel("Wall Shear Stress [Pa]", fontsize=11)
    ax2.set_ylabel("Cell Viability Score", fontsize=11)
    ax2.set_title("Sigmoidal Viability Function\n(neural / GBM cells)", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Layer 2 — Shear-Induced Cell Damage", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    test_points = [(0.25, 10), (0.5, 20), (1.0, 20), (1.5, 30)]
    for c_alg, c_fib in test_points:
        r = run(c_alg, c_fib)
        print(f"c_alg={c_alg} wt%, c_fib={c_fib:4.0f} mg/mL | "
              f"tau_wall={r['tau_wall']:.2f} Pa  viability={r['viability_shear']:.3f}")
    plot()
