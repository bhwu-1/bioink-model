"""
main.py
=======
Entry point for the GBM bioink formulation optimisation model.

Runs all seven model layers at a set of representative fibrin–alginate
formulations, then executes the 2-D process window sweep and saves figures.
"""

import os
import numpy as np

os.makedirs("figures", exist_ok=True)

import layer1_viscosity as L1
import layer2_shear as L2
import layer3_kinetics as L3
import layer4_diffusion as L4
import layer5_shape as L5
import layer6_stiffness as L6
import layer7_ecm as L7
import process_window as PW
from parameters import c_Ca_bath, t_culture

# Representative formulation points to evaluate
# (c_alg [wt%], c_fib [mg/mL], label)
TEST_POINTS = [
    (0.25, 10,  "low-low"),
    (0.50, 20,  "nominal"),
    (1.00, 20,  "hi-alg"),
    (1.50, 30,  "hi-hi"),
]


def run_all_layers():
    print("=" * 65)
    print("BIOINK MODEL — LAYER-BY-LAYER SUMMARY  (GBM fibrin–alginate)")
    print("=" * 65)

    for c_alg, c_fib, label in TEST_POINTS:
        print(f"\n--- {label:8s}  c_alg={c_alg} wt%,  c_fib={c_fib} mg/mL ---")

        r1 = L1.run(c_alg, c_fib)
        print(f"  L1 Viscosity : K={r1['K']:.4f} Pa·s^n  n={r1['n']:.3f}  "
              f"gamma_w={r1['gamma_wall']:.1f} 1/s  eta_w={r1['eta_wall']:.4f} Pa·s")

        r2 = L2.run(c_alg, c_fib)
        print(f"  L2 Shear     : tau_wall={r2['tau_wall']:.2f} Pa  "
              f"viability={r2['viability_shear']:.3f}")

        r3 = L3.run()
        print(f"  L3 Kinetics  : k_cross={r3['k_cross']:.3e} 1/s  "
              f"alpha(300s)={r3['alpha_final']:.4f}")

        r4 = L4.run()
        print(f"  L4 Diffusion : Ca2+ centre conc (300 s) = "
              f"{r4['c_centre_final']:.2f} mol/m3  (bath={c_Ca_bath:.0f} mol/m3)")

        r5 = L5.run(c_alg, c_fib)
        print(f"  L5 Shape     : t_gel={r5['t_gel']:.1f} s  "
              f"t_spread={r5['t_spread']:.4f} s  SF={r5['shape_fidelity']:.4f}")

        r6 = L6.run(c_alg, c_fib)
        print(f"  L6 Stiffness : G'={r6['G_prime']:.2f} Pa  "
              f"stiffness_viability={r6['stiffness_viability']:.3f}")

        r7 = L7.run(c_alg, c_fib)
        print(f"  L7 ECM       : G'_eff(30d)={r7['G_eff_final']:.1f} Pa  "
              f"viability(30d)={r7['viability_final']:.3f}")


def run_layer_figures():
    print("\n" + "=" * 65)
    print("GENERATING LAYER FIGURES")
    print("=" * 65)
    L1.plot()
    L2.plot()
    L3.plot()
    L4.plot()
    L5.plot()
    L6.plot()
    L7.plot()


def run_process_window():
    print("\n" + "=" * 65)
    print("PROCESS WINDOW SWEEP  (fibrinogen x alginate)")
    print("=" * 65)
    PW.run()


if __name__ == "__main__":
    run_all_layers()
    run_layer_figures()
    run_process_window()
    print("\nDone. All figures saved to ./figures/")
