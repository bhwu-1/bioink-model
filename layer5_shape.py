"""
layer5_shape.py
===============
Layer 5 — Printed filament shape fidelity.

Estimates whether the extruded filament retains its cylindrical shape before
gelation arrests viscous spreading. Uses a dimensionless shape-fidelity index:

    SF = t_spread / t_gel

SF < 1 → good shape retention; SF >> 1 → filament collapses before gelling.

Note: gelation here refers to Ca-alginate crosslinking (fast, ionic).
Fibrin polymerisation is modelled separately in layer3_kinetics.py.
"""

import numpy as np
from parameters import (K_func, n_func, R_channel, Q_print, T_bath,
                        c_alg_nom, c_fib_nom)
from layer1_viscosity import wall_shear_rate, apparent_viscosity
from layer3_kinetics import crosslink_rate_constant


def gelation_time(k: float, alpha_gel: float = 0.9) -> float:
    """
    Time [s] to reach gel conversion alpha_gel (first-order kinetics).

    t_gel = -ln(1 - alpha_gel) / k
    """
    return -np.log(1.0 - alpha_gel) / k


def viscous_spreading_time(eta: float, R: float = R_channel) -> float:
    """
    Characteristic viscous spreading time [s] after filament deposition.

    t_spread ~ eta * R / sigma_s

    sigma_s is the surface tension of the bioink.
    Value below is a placeholder — REPLACE with measured value for fibrin-alginate.
    """
    sigma_s = 0.04  # [N/m]  Surface tension — REPLACE; protein lowers vs pure water
    return eta * R / sigma_s


def shape_fidelity_index(t_gel: float, t_spread: float) -> float:
    """Dimensionless SF = t_spread / t_gel."""
    return t_spread / t_gel


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom,
        Q: float = Q_print, R: float = R_channel, T: float = T_bath) -> dict:
    """
    Compute shape fidelity metrics for the fibrin–alginate bioink.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    Q     : volumetric flow rate [m³/s]
    R     : nozzle radius [m]
    T     : crosslinking bath temperature [K]
    """
    K = K_func(c_alg, c_fib)
    n = n_func(c_alg, c_fib)
    gamma_w = wall_shear_rate(Q, R, n)
    eta_w = apparent_viscosity(gamma_w, K, n)
    k = crosslink_rate_constant(T)
    t_gel = gelation_time(k)
    t_spread = viscous_spreading_time(eta_w, R)
    SF = shape_fidelity_index(t_gel, t_spread)
    return {"eta_wall": eta_w, "t_gel": t_gel, "t_spread": t_spread, "shape_fidelity": SF}


def plot(save_path: str = "figures/layer5_shape.png"):
    """
    Left panel  — SF index bar chart (SF < 1 = good shape retention).
    Right panel — grouped bars comparing t_gel vs t_spread per formulation.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25%\n10 mg/mL"),
        (0.50, 20,  "0.50%\n20 mg/mL\n(nominal)"),
        (1.00, 20,  "1.00%\n20 mg/mL"),
        (1.50, 30,  "1.50%\n30 mg/mL"),
    ]
    colors = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c"]

    results  = [run(c_alg, c_fib) for c_alg, c_fib, _ in test_points]
    sfs      = [r["shape_fidelity"] for r in results]
    t_gels   = [r["t_gel"]          for r in results]
    t_spreads= [r["t_spread"]       for r in results]
    labels   = [p[2]                for p in test_points]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: SF index bars ---
    bar_colors = ["#4daf4a" if sf < 1 else "#e41a1c" for sf in sfs]
    bars = ax1.bar(labels, sfs, color=bar_colors, edgecolor="k", linewidth=0.8)
    ax1.axhline(1.0, color="red", linestyle="--", linewidth=1.5,
                label="SF = 1.0  (spreading = gelation)")
    for bar, sf in zip(bars, sfs):
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + max(sfs) * 0.02,
                 f"{sf:.4f}", ha="center", va="bottom", fontsize=9)
    ax1.set_ylabel("Shape Fidelity Index  SF = t_spread / t_gel", fontsize=10)
    ax1.set_title("Shape Fidelity Index\n(lower is better; green = safe)", fontsize=11)
    ax1.legend(fontsize=9)

    # --- Right: t_gel vs t_spread grouped bars (log scale) ---
    x = np.arange(len(labels))
    w = 0.35
    ax2.bar(x - w / 2, t_gels,    w, label="t_gel [s]",    color="#377eb8", edgecolor="k")
    ax2.bar(x + w / 2, t_spreads, w, label="t_spread [s]", color="#ff7f00", edgecolor="k")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_yscale("log")
    ax2.set_ylabel("Time [s]  (log scale)", fontsize=11)
    ax2.set_title("Gelation Time vs Spreading Time", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, which="both", axis="y", alpha=0.3)

    fig.suptitle("Layer 5 — Filament Shape Fidelity", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    test_points = [(0.25, 10), (0.5, 20), (1.0, 20), (1.5, 30)]
    for c_alg, c_fib in test_points:
        r = run(c_alg, c_fib)
        print(f"c_alg={c_alg} wt%, c_fib={c_fib:4.0f} mg/mL | "
              f"t_gel={r['t_gel']:.1f} s  t_spread={r['t_spread']:.4f} s  "
              f"SF={r['shape_fidelity']:.4f}")
    plot()
