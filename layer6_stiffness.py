"""
layer6_stiffness.py
===================
Layer 6 — Gel stiffness and cell mechanosensing response.

Predicts the storage modulus G' of the crosslinked fibrin–alginate gel using
the bivariate model G_prime_func(c_alg, c_fib) from parameters.py.

Computes a Gaussian viability/functionality score based on how closely G'
matches G_prime_opt = 300 Pa (brain / GBM tissue target).

Design intent: the freshly printed gel starts at ~30 Pa; ECM deposition
(layer7) drives G'_eff toward 300 Pa over the 30-day culture period.
"""

import numpy as np
from parameters import (G_prime_func, G_prime_opt, sigma_stiff,
                        c_alg_nom, c_fib_nom)


def stiffness_viability(G_prime: float, G_opt: float = G_prime_opt,
                        sigma: float = sigma_stiff) -> float:
    """
    Gaussian cell viability/functionality score (0–1) based on substrate stiffness.

    score = exp( -(G' - G'_opt)^2 / (2 * sigma^2) )

    Score = 1 at G' = G'_opt = 300 Pa; falls off with deviation.
    """
    return float(np.exp(-0.5 * ((G_prime - G_opt) / sigma) ** 2))


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom,
        alpha_gel: float = 1.0) -> dict:
    """
    Compute gel modulus and stiffness-based viability.

    Parameters
    ----------
    c_alg     : alginate concentration [wt%]
    c_fib     : fibrinogen concentration [mg/mL]
    alpha_gel : gel conversion fraction (0–1); 1 = fully crosslinked
    """
    G_prime = G_prime_func(c_alg, c_fib, alpha_gel)
    score = stiffness_viability(G_prime)
    return {"G_prime": G_prime, "stiffness_viability": score}


def plot(save_path: str = "figures/layer6_stiffness.png"):
    """
    Left panel  — G' bar chart with G'_opt = 300 Pa reference line.
    Right panel — Gaussian mechanosensing score curve with formulations marked.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25%\n10 mg/mL"),
        (0.50, 20,  "0.50%\n20 mg/mL\n(nominal)"),
        (1.00, 20,  "1.00%\n20 mg/mL"),
        (1.50, 30,  "1.50%\n30 mg/mL"),
    ]
    colors  = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c"]
    results = [run(c_alg, c_fib) for c_alg, c_fib, _ in test_points]
    G_vals  = [r["G_prime"]             for r in results]
    scores  = [r["stiffness_viability"] for r in results]
    labels  = [p[2]                     for p in test_points]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: G' bars ---
    bars = ax1.bar(labels, G_vals, color=colors, edgecolor="k", linewidth=0.8)
    ax1.axhline(G_prime_opt, color="red", linestyle="--", linewidth=1.5,
                label=f"G_opt = {G_prime_opt:.0f} Pa  (brain tissue)")
    for bar, G in zip(bars, G_vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 f"{G:.1f}", ha="center", va="bottom", fontsize=9)
    ax1.set_ylabel("Storage Modulus G' [Pa]", fontsize=11)
    ax1.set_title("Initial Gel Stiffness\nby Formulation", fontsize=11)
    ax1.legend(fontsize=9)

    # --- Right: Gaussian mechanosensing curve ---
    G_range = np.linspace(0, 900, 500)
    M_curve = np.array([stiffness_viability(g) for g in G_range])
    ax2.plot(G_range, M_curve, "b-", linewidth=2.5, label="Mechanosensing score M(G')")
    ax2.axvline(G_prime_opt, color="red", linestyle="--", linewidth=1.5,
                label=f"G_opt = {G_prime_opt:.0f} Pa")
    for (c_alg, c_fib, lbl), color, G, score in zip(test_points, colors, G_vals, scores):
        ax2.plot(G, score, "o", color=color, markersize=9, zorder=5,
                 label=lbl.replace("\n", " "))
    ax2.set_xlabel("Storage Modulus G' [Pa]", fontsize=11)
    ax2.set_ylabel("Mechanosensing Score M", fontsize=11)
    ax2.set_title("Gaussian Mechanosensing Score\nvs Substrate Stiffness", fontsize=11)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Layer 6 — Gel Stiffness & Cell Mechanosensing", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    test_points = [(0.25, 10), (0.5, 20), (1.0, 20), (1.5, 30)]
    for c_alg, c_fib in test_points:
        r = run(c_alg, c_fib)
        print(f"c_alg={c_alg} wt%, c_fib={c_fib:4.0f} mg/mL | "
              f"G'={r['G_prime']:.2f} Pa  stiffness_viability={r['stiffness_viability']:.3f}")
    plot()
