"""
process_window.py
=================
Process window analysis for the GBM fibrin–alginate bioink.

Sweeps over:
  X-axis — fibrinogen concentration  (10–30 mg/mL)
  Y-axis — alginate concentration    (0.25–1.5 wt%)

Evaluates layers 2, 5, 6, 7 at each grid point, computes a weighted composite
printability / biocompatibility score, and saves two figures:
  figures/process_window.png          — composite score heatmap
  figures/process_window_subscores.png — four-panel sub-score breakdown
"""

import numpy as np
import matplotlib.pyplot as plt
from parameters import (
    R_channel, T_bath,
    c_fib_min, c_fib_max, c_fib_points,
    c_alg_min, c_alg_max, c_alg_points,
    w_shear, w_shape, w_stiff, w_ecm,
    c_fib_nom, c_alg_nom,
)

import layer2_shear as L2
import layer5_shape as L5
import layer6_stiffness as L6
import layer7_ecm as L7


def composite_score(viability_shear: float, shape_fidelity: float,
                    stiffness_viability: float, ecm_viability: float) -> float:
    """
    Weighted composite printability / biocompatibility score (0–1).

    Weights (w_shear, w_shape, w_stiff, w_ecm) are loaded from parameters.py.
    Shape fidelity index (low = good) is converted to a 0–1 score via exp(-SF).
    """
    sf_score = np.exp(-shape_fidelity)
    return (w_shear * viability_shear
            + w_shape * sf_score
            + w_stiff * stiffness_viability
            + w_ecm   * ecm_viability)


def sweep(c_fib_arr=None, c_alg_arr=None,
          R: float = R_channel, T: float = T_bath):
    """
    Evaluate all layers across the fibrinogen × alginate design space.

    Parameters
    ----------
    c_fib_arr : fibrinogen concentrations [mg/mL]; defaults to linspace from parameters
    c_alg_arr : alginate concentrations   [wt%];   defaults to linspace from parameters
    R         : nozzle radius [m]
    T         : crosslinking bath temperature [K]

    Returns
    -------
    FIB_grid : 2-D fibrinogen array  (shape: n_alg × n_fib)
    ALG_grid : 2-D alginate array    (shape: n_alg × n_fib)
    scores   : 2-D composite score   (shape: n_alg × n_fib)
    details  : dict of 2-D sub-score arrays
    """
    if c_fib_arr is None:
        c_fib_arr = np.linspace(c_fib_min, c_fib_max, c_fib_points)
    if c_alg_arr is None:
        c_alg_arr = np.linspace(c_alg_min, c_alg_max, c_alg_points)

    # meshgrid: rows index alginate, cols index fibrinogen
    FIB_grid, ALG_grid = np.meshgrid(c_fib_arr, c_alg_arr)
    scores  = np.zeros_like(FIB_grid)
    v_shear = np.zeros_like(FIB_grid)
    v_shape = np.zeros_like(FIB_grid)
    v_stiff = np.zeros_like(FIB_grid)
    v_ecm   = np.zeros_like(FIB_grid)

    n_total = len(c_alg_arr) * len(c_fib_arr)
    count = 0
    for i, c_alg in enumerate(c_alg_arr):
        for j, c_fib in enumerate(c_fib_arr):
            shear_res = L2.run(c_alg, c_fib, R=R)
            shape_res = L5.run(c_alg, c_fib, R=R, T=T)
            stiff_res = L6.run(c_alg, c_fib)
            ecm_res   = L7.run(c_alg, c_fib)

            vs = shear_res["viability_shear"]
            sf = shape_res["shape_fidelity"]
            sv = stiff_res["stiffness_viability"]
            ev = ecm_res["viability_final"]

            v_shear[i, j] = vs
            v_shape[i, j] = sf
            v_stiff[i, j] = sv
            v_ecm[i, j]   = ev
            scores[i, j]  = composite_score(vs, sf, sv, ev)

            count += 1
            if count % 50 == 0 or count == n_total:
                print(f"  sweep progress: {count}/{n_total}")

    details = {
        "viability_shear":     v_shear,
        "shape_fidelity":      v_shape,
        "stiffness_viability": v_stiff,
        "ecm_viability":       v_ecm,
    }
    return FIB_grid, ALG_grid, scores, details


def plot_process_window(FIB_grid, ALG_grid, scores,
                        save_path: str = "figures/process_window.png"):
    """Composite score heatmap with nominal formulation marker."""
    fig, ax = plt.subplots(figsize=(8, 6))
    cf = ax.contourf(FIB_grid, ALG_grid, scores, levels=20,
                     cmap="RdYlGn", vmin=0, vmax=1)
    cs = ax.contour(FIB_grid, ALG_grid, scores, levels=[0.4, 0.6, 0.8],
                    colors="k", linewidths=0.8)
    ax.clabel(cs, fmt="%.1f", fontsize=8)
    plt.colorbar(cf, ax=ax, label="Composite Score (0–1)")

    ax.plot(c_fib_nom, c_alg_nom, "w*", markersize=14,
            label=f"Nominal ({c_fib_nom:.0f} mg/mL, {c_alg_nom} wt%)")
    ax.legend(fontsize=8, loc="upper left")

    ax.set_xlabel("Fibrinogen Concentration [mg/mL]")
    ax.set_ylabel("Alginate Concentration [wt%]")
    ax.set_title("GBM Bioink Process Window\n"
                 "Aspect Biosystems RX1 — fibrin–alginate composite")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Process window saved to {save_path}")
    plt.close(fig)


def plot_subscores(FIB_grid, ALG_grid, details,
                   save_path: str = "figures/process_window_subscores.png"):
    """Four-panel breakdown of individual score components."""
    panels = [
        ("viability_shear",     "Shear Viability",              "RdYlGn", 0, 1),
        ("shape_fidelity",      "Shape Fidelity (SF index)\nlow = good", "RdYlGn_r", None, None),
        ("stiffness_viability", "Initial Stiffness Score\n(vs G\'_opt = 300 Pa)", "RdYlGn", 0, 1),
        ("ecm_viability",       "ECM Stiffness Score (30 d)\n(vs G\'_opt = 300 Pa)", "RdYlGn", 0, 1),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (key, title, cmap, vmin, vmax) in zip(axes.flat, panels):
        cf = ax.contourf(FIB_grid, ALG_grid, details[key], levels=20,
                         cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(cf, ax=ax)
        ax.plot(c_fib_nom, c_alg_nom, "w*", markersize=10)
        ax.set_xlabel("Fibrinogen [mg/mL]")
        ax.set_ylabel("Alginate [wt%]")
        ax.set_title(title, fontsize=10)

    fig.suptitle("GBM Bioink — Process Window Sub-scores", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Sub-score plot saved to {save_path}")
    plt.close(fig)


def run():
    print("Running process window sweep...")
    FIB_grid, ALG_grid, scores, details = sweep()
    plot_process_window(FIB_grid, ALG_grid, scores)
    plot_subscores(FIB_grid, ALG_grid, details)

    best_idx = np.unravel_index(np.argmax(scores), scores.shape)
    print(f"\nBest score : {scores[best_idx]:.3f}")
    print(f"  at c_fib = {FIB_grid[best_idx]:.1f} mg/mL, "
          f"c_alg = {ALG_grid[best_idx]:.3f} wt%")
    return FIB_grid, ALG_grid, scores, details


if __name__ == "__main__":
    run()
