"""
layer3_kinetics.py
==================
Layer 3 — Crosslinking regime analysis via the Damköhler number.

Ca-alginate ionic crosslinking is chemistry-fast but diffusion-limited:
Ca²⁺ must diffuse through the alginate network before binding occurs.
The Damköhler number quantifies this:

    Da = t_diffusion / t_reaction = (R² / D_Ca_eff) * k_cross

When Da >> 1, gelation is diffusion-controlled (always true here, Da ~ 20–50).
The physically meaningful gelation time is therefore the diffusion time:

    t_gel = R² / (2 * D_Ca_eff(c_alg))

This t_gel depends on alginate concentration through D_Ca_eff and is exported
to layer5_shape.py as the governing timescale for shape fidelity.
"""

import numpy as np
from parameters import (Ea_cross, A_cross, R_gas, T_process, T_bath,
                        D_Ca_eff, R_channel, c_alg_nom, c_alg_min, c_alg_max)


def crosslink_rate_constant(T: float = T_process) -> float:
    """
    Arrhenius Ca-alginate crosslinking rate constant [1/s].

    k(T) = A_cross * exp(-Ea_cross / (R_gas * T))

    Used for Da calculation only — not the governing timescale for gelation.
    """
    return A_cross * np.exp(-Ea_cross / (R_gas * T))


def t_gel_diffusion(c_alg: float, R: float = R_channel) -> float:
    """
    Diffusion-limited gelation time [s].

    t_gel = R² / (2 * D_Ca_eff(c_alg))

    Increases with alginate concentration as D_Ca_eff decreases (denser network).
    This is the t_gel used by layer5_shape.py.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    R     : filament radius [m]
    """
    return R**2 / (2.0 * D_Ca_eff(c_alg))


def damkohler_number(c_alg: float, T: float = T_bath, R: float = R_channel) -> float:
    """
    Damköhler number Da = t_diffusion / t_reaction.

    Da = t_gel_diffusion * k_cross

    Da >> 1 confirms diffusion-limited regime (chemistry is not rate-limiting).
    """
    k = crosslink_rate_constant(T)
    return t_gel_diffusion(c_alg, R) * k


def run(c_alg: float = c_alg_nom, T: float = T_bath) -> dict:
    """Return crosslinking regime diagnostics for a given alginate concentration."""
    k    = crosslink_rate_constant(T)
    t_gel = t_gel_diffusion(c_alg)
    Da   = damkohler_number(c_alg, T)
    D    = D_Ca_eff(c_alg)
    return {"k_cross": k, "D_Ca_eff": D, "t_gel_diffusion": t_gel, "Da": Da}


def plot(save_path: str = "figures/layer3_kinetics.png"):
    """
    Left panel  — Damköhler number vs alginate concentration.
                  Da >> 1 everywhere confirms diffusion-limited regime.
    Right panel — Diffusion-limited t_gel vs alginate concentration.
                  Shows t_gel increasing with c_alg (slower diffusion).
    """
    import matplotlib.pyplot as plt

    c_alg_range = np.linspace(c_alg_min, c_alg_max, 100)
    Da_vals    = [damkohler_number(c) for c in c_alg_range]
    t_gel_vals = [t_gel_diffusion(c)  for c in c_alg_range]

    test_points = [(0.25, "0.25%"), (0.50, "0.50% (nom.)"), (1.00, "1.00%"), (1.50, "1.50%")]
    colors = ["#a8d1e7", "#5ba3c9", "#2171b5", "#08306b"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: Da vs c_alg ---
    ax1.plot(c_alg_range, Da_vals, color="#2171b5", linewidth=2.5)
    ax1.axhline(1.0, color="red", linestyle="--", linewidth=1.2,
                label="Da = 1  (reaction-limited threshold)")
    for (c, lbl), col in zip(test_points, colors):
        ax1.plot(c, damkohler_number(c), "o", color=col, markersize=9,
                 zorder=5, label=f"{lbl}  Da={damkohler_number(c):.0f}")
    ax1.set_xlabel("Alginate Concentration [wt%]", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Damkohler Number  Da", fontsize=12, fontweight="bold")
    ax1.set_title("Crosslinking Regime\nDa >> 1 = diffusion-limited everywhere", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.5)
    ax1.tick_params(width=1.5, labelsize=11)

    # --- Right: t_gel_diffusion vs c_alg ---
    ax2.plot(c_alg_range, t_gel_vals, color="#2171b5", linewidth=2.5)
    for (c, lbl), col in zip(test_points, colors):
        ax2.plot(c, t_gel_diffusion(c), "o", color=col, markersize=9, zorder=5,
                 label=f"{lbl}  t_gel={t_gel_diffusion(c):.1f} s")
    ax2.set_xlabel("Alginate Concentration [wt%]", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Diffusion-Limited t_gel [s]", fontsize=12, fontweight="bold")
    ax2.set_title("Gelation Time (diffusion-limited)\nIncreases with alginate concentration", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.5)
    ax2.tick_params(width=1.5, labelsize=11)

    fig.suptitle("Layer 3: Ca-Alginate Crosslinking: Damkohler Analysis", fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    for c in [0.25, 0.5, 1.0, 1.5]:
        r = run(c)
        print(f"c_alg={c} wt% | D_Ca_eff={r['D_Ca_eff']:.2e} m2/s  "
              f"t_gel={r['t_gel_diffusion']:.1f} s  Da={r['Da']:.1f}")
    plot()
