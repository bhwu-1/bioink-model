"""
layer5_shape.py
===============
Layer 5 — Printed filament shape fidelity.

Estimates whether the extruded filament retains its cylindrical shape before
gelation arrests viscous spreading.  Uses a dimensionless shape-fidelity index:

    SF = t_spread / t_gel

where:
  t_spread = eta * R / sigma_s   (viscous resistance to surface-tension spreading)
  t_gel    = R² / (2 * D_Ca_eff) (diffusion-limited gelation time from Layer 3)

Convention:
  SF > 1  → t_spread > t_gel → filament GELS before spreading significantly → GOOD
  SF < 1  → t_spread < t_gel → filament SPREADS before gelling              → POOR

The tension that creates the process window boundary:
  Low alginate  → low eta → short t_spread → small SF → poor shape (collapses)
  High alginate → longer t_gel (slower D_Ca_eff) → also lowers SF, but slower
  Nominal       → eta high enough AND D_Ca fast enough → SF maximised near (0.5%, 20 mg/mL)
"""

import numpy as np
from parameters import (K_func, n_func, R_channel, Q_print,
                        c_alg_nom, c_fib_nom)
from layer1_viscosity import wall_shear_rate, apparent_viscosity
from layer3_kinetics import t_gel_diffusion


def viscous_spreading_time(eta: float, R: float = R_channel) -> float:
    """
    Characteristic time [s] for surface-tension-driven viscous spreading.

    t_spread = eta * R / sigma_s

    Physical basis: surface tension (sigma_s) drives spreading; viscosity (eta)
    resists it.  High eta → large t_spread → better shape retention.

    sigma_s = 0.050 N/m — dilute aqueous hydrogel, reduced from water (72 mN/m)
    by dissolved polymer.  Ref: alginate surface tension literature ~40–70 mN/m.
    """
    sigma_s = 0.050  # [N/m]  Surface tension of fibrin-alginate bioink  [EST]
                     #  — REPLACE with pendant-drop tensiometry measurement
    return eta * R / sigma_s


def shape_fidelity_index(t_gel: float, t_spread: float) -> float:
    """
    Dimensionless SF = t_spread / t_gel.

    SF > 1: gels before spreading → GOOD shape retention
    SF < 1: spreads before gelling → POOR shape retention
    """
    return t_spread / t_gel


def run(c_alg: float = c_alg_nom, c_fib: float = c_fib_nom,
        Q: float = Q_print, R: float = R_channel, **kwargs) -> dict:
    """
    Compute shape fidelity metrics for the fibrin–alginate bioink.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    Q     : volumetric flow rate [m³/s]
    R     : nozzle radius [m]

    Note: t_gel is diffusion-limited (from layer3_kinetics.t_gel_diffusion),
    NOT from Arrhenius kinetics.  Da >> 1 everywhere in this system.
    """
    K = K_func(c_alg, c_fib)
    n = n_func(c_alg, c_fib)
    gamma_w  = wall_shear_rate(Q, R, n)
    eta_w    = apparent_viscosity(gamma_w, K, n)
    t_gel    = t_gel_diffusion(c_alg, R)      # diffusion-limited gelation time
    t_spread = viscous_spreading_time(eta_w, R)
    SF       = shape_fidelity_index(t_gel, t_spread)
    return {"eta_wall": eta_w, "t_gel": t_gel, "t_spread": t_spread, "shape_fidelity": SF}


def plot(save_path: str = "figures/layer5_shape.png"):
    """
    Left panel  — SF index bars (SF > 1 = good; SF < 1 = poor shape retention).
    Right panel — t_gel vs t_spread on log scale showing the 3-order-of-magnitude gap.
    """
    import matplotlib.pyplot as plt

    test_points = [
        (0.25, 10,  "0.25%\n10 mg/mL"),
        (0.50, 20,  "0.50%\n20 mg/mL\n(nominal)"),
        (1.00, 20,  "1.00%\n20 mg/mL"),
        (1.50, 30,  "1.50%\n30 mg/mL"),
    ]
    colors = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c"]

    results   = [run(c_alg, c_fib) for c_alg, c_fib, _ in test_points]
    sfs       = [r["shape_fidelity"] for r in results]
    t_gels    = [r["t_gel"]          for r in results]
    t_spreads = [r["t_spread"]       for r in results]
    labels    = [p[2]                for p in test_points]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: SF index bars (SF > 1 = good = green) ---
    bar_colors = ["#4daf4a" if sf > 1 else "#e41a1c" for sf in sfs]
    bars = ax1.bar(labels, sfs, color=bar_colors, edgecolor="k", linewidth=0.8)
    ax1.axhline(1.0, color="red", linestyle="--", linewidth=1.5,
                label="SF = 1.0  (gelation = spreading)")
    for bar, sf in zip(bars, sfs):
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + max(sfs) * 0.02,
                 f"{sf:.4f}", ha="center", va="bottom", fontsize=9)
    ax1.set_ylabel("Shape Fidelity Index  SF = t_spread / t_gel", fontsize=10)
    ax1.set_title("Shape Fidelity Index\n(SF > 1 = good; green = gels before spreading)", fontsize=11)
    ax1.legend(fontsize=9)

    # --- Right: t_gel vs t_spread (log scale) ---
    x = np.arange(len(labels))
    w = 0.35
    ax2.bar(x - w / 2, t_gels,    w, label="t_gel [s]  (diffusion-limited)",
            color="#377eb8", edgecolor="k")
    ax2.bar(x + w / 2, t_spreads, w, label="t_spread [s]  (viscous spreading)",
            color="#ff7f00", edgecolor="k")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_yscale("log")
    ax2.set_ylabel("Time [s]  (log scale)", fontsize=11)
    ax2.set_title("Diffusion t_gel vs Viscous t_spread\n(gap narrows at high alginate)", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, which="both", axis="y", alpha=0.3)

    fig.suptitle("Layer 5 — Filament Shape Fidelity  (diffusion-limited gelation)", fontsize=13)
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
