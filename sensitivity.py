"""
sensitivity.py
==============
One-at-a-time (OAT) parameter sensitivity analysis for the GBM bioink model.

For each uncertain parameter, runs two full 2-D process window sweeps
(+20 % / -20 % from nominal) and records:
  (1) Composite score at the nominal formulation (20 mg/mL fib, 0.5 wt% alg)
  (2) Predicted-optimum location (c_fib*, c_alg*) and its score

The swing |score(+20%) - score(-20%)| at the nominal point is used as the
sensitivity metric for ranking in the tornado plot.

Parameters analysed
-------------------
tau_crit   [LIT]          Neural cell shear damage threshold — low-sensitivity control
D_Ca       [LIT/central]  Ca2+ diffusion coefficient at nominal alginate
SF_ref     [EST]          Shape fidelity score reference value
Ea_cross   [EST]          Ca-alginate crosslinking activation energy
A_cross    [PLACEHOLDER]  Arrhenius pre-exponential factor
w_shear    [EST]          Composite score weight — shear viability
w_shape    [EST]          Composite score weight — shape fidelity
w_stiff    [EST]          Composite score weight — initial stiffness
w_ecm      [EST]          Composite score weight — ECM trajectory

Note on Ea_cross / A_cross
--------------------------
These parameters enter only through the Damkohler number Da = t_gel * k_cross.
Because Da >> 1 everywhere (diffusion-limited regime), t_gel = R^2 / (2 * D_Ca_eff)
is the governing timescale — independent of Arrhenius kinetics.  Consequently,
Ea_cross and A_cross are expected to show near-zero sensitivity in the composite
score.  This is itself an important finding: it confirms that these PLACEHOLDER
values do not need urgent experimental validation.

Weight perturbations
--------------------
Weights are constrained to sum to 1.0.  Each weight is perturbed by +/-W_DELTA
(default 0.05) and the complement is redistributed proportionally among the
remaining three weights.

Outputs
-------
figures/sensitivity_tornado.png       — parameters ranked by score swing
figures/sensitivity_optimum_shift.png — scatter of optimum location shifts
Console summary table
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import parameters as params
import layer2_shear as L2
import layer3_kinetics as L3
import layer5_shape as L5
import layer6_stiffness as L6
import layer7_ecm as L7

from parameters import (
    c_fib_nom, c_alg_nom,
    c_fib_min, c_fib_max,
    c_alg_min, c_alg_max,
)

# ── Tunable constants ──────────────────────────────────────────────────────────
_N_FIB  = 15     # grid points along fibrinogen axis (increase for higher fidelity)
_N_ALG  = 15     # grid points along alginate axis
_DELTA  = 0.20   # fractional perturbation for scalar params (±20 %)
_W_DELTA = 0.05  # absolute weight shift (reallocated among remaining weights)

# ── Tag colour palette (matching existing blue-palette style) ─────────────────
_TAG_COLORS = {
    "[LIT]":         "#2171b5",   # blue  — well-supported
    "[EST]":         "#41ab5d",   # green — estimated, needs validation
    "[PLACEHOLDER]": "#d94801",   # orange — rough order of magnitude
}


# =============================================================================
# Core sweep engine with explicit parameter overrides
# =============================================================================

def _sweep(c_fib_arr, c_alg_arr,
           tau_crit_val, D_Ca_val, SF_ref_val,
           Ea_cross_val, A_cross_val,
           ws, wsh, wst, we):
    """
    Full 2-D grid sweep with explicit parameter overrides.

    D_Ca is temporarily patched in the parameters module because D_Ca_eff()
    reads it at call-time (not at import time).  Ea_cross and A_cross are
    patched in both parameters and layer3_kinetics (which imported them as
    values at load time), though they do not affect the composite score.

    tau_crit is passed explicitly to cell_viability_shear() to bypass its
    default-argument binding.  SF_ref and weights are used inline below.

    Returns
    -------
    FIB_grid, ALG_grid : 2-D meshgrid arrays  (n_alg x n_fib)
    scores             : 2-D composite score array
    """
    # Save and patch module-level parameters used at call-time
    saved = {
        "params.D_Ca":     params.D_Ca,
        "params.Ea_cross": params.Ea_cross,
        "params.A_cross":  params.A_cross,
        "L3.Ea_cross":     L3.Ea_cross,
        "L3.A_cross":      L3.A_cross,
    }
    params.D_Ca     = D_Ca_val
    params.Ea_cross = Ea_cross_val
    params.A_cross  = A_cross_val
    L3.Ea_cross     = Ea_cross_val
    L3.A_cross      = A_cross_val

    try:
        FIB_grid, ALG_grid = np.meshgrid(c_fib_arr, c_alg_arr)
        scores = np.zeros_like(FIB_grid)

        for i, c_alg in enumerate(c_alg_arr):
            for j, c_fib in enumerate(c_fib_arr):
                # Layer 2: get tau_wall, apply overridden tau_crit explicitly
                tau_w = L2.run(c_alg, c_fib)["tau_wall"]
                vs = L2.cell_viability_shear(tau_w, tau_c=tau_crit_val)

                # Layer 5: shape fidelity (uses patched D_Ca via D_Ca_eff)
                sf = L5.run(c_alg, c_fib)["shape_fidelity"]

                # Layers 6 & 7: stiffness scores (independent of kinetics params)
                sv = L6.run(c_alg, c_fib)["stiffness_viability"]
                ev = L7.run(c_alg, c_fib)["viability_final"]

                # Composite score with overridden SF_ref and weights
                sf_score = 1.0 - np.exp(-sf / SF_ref_val)
                scores[i, j] = ws * vs + wsh * sf_score + wst * sv + we * ev

    finally:
        # Always restore — even if a sweep raises an exception
        params.D_Ca     = saved["params.D_Ca"]
        params.Ea_cross = saved["params.Ea_cross"]
        params.A_cross  = saved["params.A_cross"]
        L3.Ea_cross     = saved["L3.Ea_cross"]
        L3.A_cross      = saved["L3.A_cross"]

    return FIB_grid, ALG_grid, scores


def _nominal_score(scores, FIB_grid, ALG_grid):
    """Score at the grid point closest to the nominal formulation."""
    j = int(np.argmin(np.abs(FIB_grid[0, :] - c_fib_nom)))
    i = int(np.argmin(np.abs(ALG_grid[:, 0] - c_alg_nom)))
    return float(scores[i, j])


def _optimum(scores, FIB_grid, ALG_grid):
    """(c_fib, c_alg, score) at the global grid maximum."""
    idx = np.unravel_index(np.argmax(scores), scores.shape)
    return float(FIB_grid[idx]), float(ALG_grid[idx]), float(scores[idx])


# =============================================================================
# Weight perturbation helper
# =============================================================================

def _perturb_weights(base_w, idx, delta):
    """
    Shift weight[idx] by delta; redistribute -delta proportionally among
    the other three weights.  Sum is preserved at 1.0.

    Parameters
    ----------
    base_w : sequence of 4 nominal weights (ws, wsh, wst, we)
    idx    : index of weight to perturb (0=shear, 1=shape, 2=stiff, 3=ecm)
    delta  : signed shift (+/- _W_DELTA)
    """
    w = list(base_w)
    others = [i for i in range(4) if i != idx]
    others_sum = sum(base_w[i] for i in others)
    w[idx] += delta
    for i in others:
        w[i] -= delta * (base_w[i] / others_sum)
    return tuple(w)


# =============================================================================
# Main analysis loop
# =============================================================================

def run(n_fib=_N_FIB, n_alg=_N_ALG):
    """
    Execute the full OAT sensitivity analysis.

    Parameters
    ----------
    n_fib : grid points along fibrinogen axis
    n_alg : grid points along alginate axis

    Returns
    -------
    results : list of dicts, one per parameter, sorted by swing (descending)
    base    : dict with base-run arrays and scalar summary values
    """
    c_fib_arr = np.linspace(c_fib_min, c_fib_max, n_fib)
    c_alg_arr = np.linspace(c_alg_min, c_alg_max, n_alg)

    # Nominal parameter values
    nom = dict(
        tau_crit = params.tau_crit,
        D_Ca     = params.D_Ca,
        SF_ref   = params.SF_ref,
        Ea_cross = params.Ea_cross,
        A_cross  = params.A_cross,
        w_shear  = params.w_shear,
        w_shape  = params.w_shape,
        w_stiff  = params.w_stiff,
        w_ecm    = params.w_ecm,
    )
    base_weights = (nom["w_shear"], nom["w_shape"], nom["w_stiff"], nom["w_ecm"])

    def _run_with(overrides):
        p = {**nom, **overrides}
        return _sweep(
            c_fib_arr, c_alg_arr,
            tau_crit_val = p["tau_crit"],
            D_Ca_val     = p["D_Ca"],
            SF_ref_val   = p["SF_ref"],
            Ea_cross_val = p["Ea_cross"],
            A_cross_val  = p["A_cross"],
            ws  = p["w_shear"],
            wsh = p["w_shape"],
            wst = p["w_stiff"],
            we  = p["w_ecm"],
        )

    # --- Base sweep ---
    print("Running base sweep...")
    F0, A0, S0 = _run_with({})
    base = {
        "FIB_grid":  F0,
        "ALG_grid":  A0,
        "scores":    S0,
        "nom_score": _nominal_score(S0, F0, A0),
        "opt_fib":   _optimum(S0, F0, A0)[0],
        "opt_alg":   _optimum(S0, F0, A0)[1],
        "opt_score": _optimum(S0, F0, A0)[2],
    }

    # --- Parameter definitions ---
    # (display_label, param_key, nominal_value, tag, is_weight, weight_index)
    param_defs = [
        ("tau_crit",  "tau_crit",  nom["tau_crit"],  "[LIT]",         False, None),
        ("D_Ca",      "D_Ca",      nom["D_Ca"],       "[LIT]",         False, None),
        ("SF_ref",    "SF_ref",    nom["SF_ref"],     "[EST]",         False, None),
        ("Ea_cross",  "Ea_cross",  nom["Ea_cross"],   "[EST]",         False, None),
        ("A_cross",   "A_cross",   nom["A_cross"],    "[PLACEHOLDER]", False, None),
        ("w_shear",   "w_shear",   nom["w_shear"],    "[EST]",         True,  0),
        ("w_shape",   "w_shape",   nom["w_shape"],    "[EST]",         True,  1),
        ("w_stiff",   "w_stiff",   nom["w_stiff"],    "[EST]",         True,  2),
        ("w_ecm",     "w_ecm",     nom["w_ecm"],      "[EST]",         True,  3),
    ]

    results = []
    total = len(param_defs)
    for k, (label, key, val, tag, is_weight, widx) in enumerate(param_defs):
        print(f"  [{k + 1}/{total}] Perturbing {label} ...")

        if is_weight:
            w_hi = _perturb_weights(base_weights, widx, +_W_DELTA)
            w_lo = _perturb_weights(base_weights, widx, -_W_DELTA)
            ov_hi = dict(zip(["w_shear", "w_shape", "w_stiff", "w_ecm"], w_hi))
            ov_lo = dict(zip(["w_shear", "w_shape", "w_stiff", "w_ecm"], w_lo))
        else:
            ov_hi = {key: val * (1.0 + _DELTA)}
            ov_lo = {key: val * (1.0 - _DELTA)}

        F_hi, A_hi, S_hi = _run_with(ov_hi)
        F_lo, A_lo, S_lo = _run_with(ov_lo)

        s_hi = _nominal_score(S_hi, F_hi, A_hi)
        s_lo = _nominal_score(S_lo, F_lo, A_lo)
        swing = abs(s_hi - s_lo)

        opt_fib_hi, opt_alg_hi, opt_score_hi = _optimum(S_hi, F_hi, A_hi)
        opt_fib_lo, opt_alg_lo, opt_score_lo = _optimum(S_lo, F_lo, A_lo)

        results.append({
            "label":        label,
            "key":          key,
            "nominal_val":  val,
            "tag":          tag,
            "is_weight":    is_weight,
            "score_hi":     s_hi,
            "score_lo":     s_lo,
            "swing":        swing,
            "opt_fib_hi":   opt_fib_hi,
            "opt_alg_hi":   opt_alg_hi,
            "opt_score_hi": opt_score_hi,
            "opt_fib_lo":   opt_fib_lo,
            "opt_alg_lo":   opt_alg_lo,
            "opt_score_lo": opt_score_lo,
        })

    results.sort(key=lambda r: r["swing"], reverse=True)
    return results, base


# =============================================================================
# Figures
# =============================================================================

def plot_tornado(results, base,
                 save_path="figures/sensitivity_tornado.png"):
    """
    Horizontal tornado chart: parameters sorted by score swing, largest at top.
    Bars are coloured by parameter tag ([LIT] / [EST] / [PLACEHOLDER]).
    """
    # Sort ascending so largest ends up at top after barh
    sorted_r = sorted(results, key=lambda r: r["swing"])

    labels  = [r["label"] for r in sorted_r]
    swings  = [r["swing"] for r in sorted_r]
    tags    = [r["tag"]   for r in sorted_r]
    colors  = [_TAG_COLORS.get(t, "#888888") for t in tags]

    fig, ax = plt.subplots(figsize=(9, 6))

    bars = ax.barh(labels, swings, color=colors, edgecolor="k",
                   linewidth=0.8, height=0.6)

    # Value labels on bars
    x_max = max(swings) if max(swings) > 0 else 1e-6
    for bar, swing in zip(bars, swings):
        ax.text(bar.get_width() + x_max * 0.02,
                bar.get_y() + bar.get_height() / 2,
                f"{swing:.4f}", va="center", ha="left", fontsize=9)

    # Tag legend
    legend_patches = [
        mpatches.Patch(color=_TAG_COLORS["[LIT]"],
                       label="[LIT] — literature (included as sensitivity control)"),
        mpatches.Patch(color=_TAG_COLORS["[EST]"],
                       label="[EST] — estimated, needs experimental validation"),
        mpatches.Patch(color=_TAG_COLORS["[PLACEHOLDER]"],
                       label="[PLACEHOLDER] — rough order of magnitude"),
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc="lower right", framealpha=0.85)

    ax.set_xlabel(
        "Score Swing  |score(+20%) - score(-20%)|  at nominal formulation\n"
        "(weights: |score(+0.05) - score(-0.05)| with proportional redistribution)",
        fontsize=10, fontweight="bold",
    )
    ax.set_title(
        "GBM Bioink Model — Parameter Sensitivity (OAT)\n"
        f"Base score at nominal formulation: {base['nom_score']:.4f}",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlim(0, x_max * 1.30)
    ax.grid(True, axis="x", alpha=0.3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, labelsize=11)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Tornado plot saved to {save_path}")
    plt.close(fig)


def plot_optimum_shift(results, base,
                       save_path="figures/sensitivity_optimum_shift.png"):
    """
    Scatter showing how the predicted-optimum location shifts under each
    perturbation, overlaid on the base process window contour.

    Upward triangle (^) = +20 % perturbation
    Downward triangle (v) = -20 % perturbation
    Connecting line shows the direction and magnitude of the shift.
    """
    F0, A0, S0 = base["FIB_grid"], base["ALG_grid"], base["scores"]

    fig, ax = plt.subplots(figsize=(9, 7))

    # Background: base process window
    cf = ax.contourf(F0, A0, S0, levels=100, cmap="RdYlGn",
                     vmin=S0.min(), vmax=S0.max(), alpha=0.50)
    plt.colorbar(cf, ax=ax, label="Base Composite Score")
    ax.contour(F0, A0, S0, levels=6, colors="k", linewidths=0.5, alpha=0.35)

    cmap_p = plt.cm.tab10
    n = len(results)
    for k, r in enumerate(results):
        color = cmap_p(k / n)
        ax.plot(r["opt_fib_hi"], r["opt_alg_hi"], "^",
                color=color, markersize=8,
                markeredgecolor="k", markeredgewidth=0.6, zorder=4)
        ax.plot(r["opt_fib_lo"], r["opt_alg_lo"], "v",
                color=color, markersize=8,
                markeredgecolor="k", markeredgewidth=0.6, zorder=4)
        ax.plot([r["opt_fib_hi"], r["opt_fib_lo"]],
                [r["opt_alg_hi"], r["opt_alg_lo"]],
                "-", color=color, linewidth=1.5, alpha=0.75,
                label=r["label"], zorder=3)

    # Base optimum and nominal formulation markers
    ax.plot(base["opt_fib"], base["opt_alg"], "*",
            markersize=18, color="gold", markeredgecolor="k",
            markeredgewidth=1.0, zorder=6, label="Base optimum")
    ax.plot(c_fib_nom, c_alg_nom, "*",
            markersize=14, color="white", markeredgecolor="k",
            markeredgewidth=0.8, zorder=5, label="Nominal formulation")

    ax.set_xlabel("Fibrinogen Concentration [mg/mL]", fontsize=12, fontweight="bold")
    ax.set_ylabel("Alginate Concentration [wt%]", fontsize=12, fontweight="bold")
    ax.set_title(
        "Predicted Optimum Location Under Parameter Perturbations\n"
        "(triangle-up = +20%, triangle-down = -20%; line connects each pair)",
        fontsize=13, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="upper right", framealpha=0.85, ncol=2, columnspacing=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, labelsize=11)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Optimum shift plot saved to {save_path}")
    plt.close(fig)


# =============================================================================
# Console summary
# =============================================================================

def print_summary(results, base):
    """Print a formatted summary table to stdout."""
    W = 100
    print()
    print("=" * W)
    print("SENSITIVITY ANALYSIS SUMMARY — GBM Bioink Model (OAT ±20%)")
    print(f"  Base score at nominal formulation : {base['nom_score']:.4f}")
    print(f"  Base predicted optimum            : "
          f"{base['opt_fib']:.1f} mg/mL fibrinogen, "
          f"{base['opt_alg']:.3f} wt% alginate  "
          f"(score {base['opt_score']:.4f})")
    print("=" * W)
    hdr = (f"{'Parameter':<12} {'Tag':<15} {'Nom. value':<14} "
           f"{'Score +20%':>10} {'Score -20%':>10} {'Swing':>8} "
           f"{'Dopt_fib':>10} {'Dopt_alg':>10}")
    print(hdr)
    print("-" * W)
    for r in results:
        nom_str = (f"{r['nominal_val']:.2f} (+/-{_W_DELTA})" if r["is_weight"]
                   else f"{r['nominal_val']:.3g}")
        d_fib = r["opt_fib_hi"] - r["opt_fib_lo"]
        d_alg = r["opt_alg_hi"] - r["opt_alg_lo"]
        print(f"{r['label']:<12} {r['tag']:<15} {nom_str:<14} "
              f"{r['score_hi']:>10.4f} {r['score_lo']:>10.4f} {r['swing']:>8.4f} "
              f"{d_fib:>10.2f} {d_alg:>10.3f}")
    print("=" * W)

    print()
    print("Key findings:")
    for rank, r in enumerate(results[:3], 1):
        print(f"  #{rank} most influential : {r['label']} "
              f"(swing = {r['swing']:.4f}, tag = {r['tag']})")
    bottom = results[-1]
    print(f"  Least influential  : {bottom['label']} "
          f"(swing = {bottom['swing']:.4f}, tag = {bottom['tag']})")

    # Validation priority message
    est_top = [r for r in results if r["tag"] in ("[EST]", "[PLACEHOLDER]")]
    if est_top:
        print()
        print("Experimental validation priority (highest-swing [EST]/[PLACEHOLDER] first):")
        for r in est_top[:4]:
            print(f"  {r['label']:<12}  swing={r['swing']:.4f}  {r['tag']}")
    print()


# =============================================================================
# Entry point
# =============================================================================

def main(n_fib=_N_FIB, n_alg=_N_ALG):
    """Run full sensitivity analysis and save all outputs."""
    os.makedirs("figures", exist_ok=True)

    print("=" * 60)
    print("SENSITIVITY ANALYSIS — GBM Bioink Model")
    print(f"Grid : {n_fib} x {n_alg}  |  Perturbation : +/-{_DELTA * 100:.0f}%")
    print(f"       (weights shifted by +/-{_W_DELTA})")
    print("=" * 60)

    results, base = run(n_fib=n_fib, n_alg=n_alg)
    print_summary(results, base)
    plot_tornado(results, base)
    plot_optimum_shift(results, base)

    return results, base


if __name__ == "__main__":
    main()
