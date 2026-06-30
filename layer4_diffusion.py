"""
layer4_diffusion.py
===================
Layer 4 — Calcium ion diffusion into the alginate filament.

Models 1-D radial diffusion of Ca²⁺ from the bath into a cylindrical
alginate filament using the transient diffusion equation (Fick's 2nd law).
Solved numerically with scipy's solve_ivp via the method of lines (MOL).
"""

import numpy as np
from scipy.integrate import solve_ivp
from parameters import D_Ca, R_channel, c_Ca_bath


def build_mol_rhs(r_nodes: np.ndarray, D: float, c_surface: float):
    """
    Build the right-hand side function for the method-of-lines discretisation
    of the radial diffusion equation in cylindrical coordinates:

        dc/dt = D * (1/r) * d/dr (r * dc/dr)

    Boundary conditions:
        - Symmetry at r=0: dc/dr = 0
        - Dirichlet at r=R: c = c_surface (bath concentration)
    """
    dr = r_nodes[1] - r_nodes[0]
    N = len(r_nodes)

    def rhs(t, c):
        dcdt = np.zeros(N)
        # Interior nodes (central differences in cylindrical coords)
        for i in range(1, N - 1):
            r = r_nodes[i]
            d2c = (c[i + 1] - 2 * c[i] + c[i - 1]) / dr**2
            dc_dr = (c[i + 1] - c[i - 1]) / (2 * dr)
            dcdt[i] = D * (d2c + dc_dr / r)
        # Symmetry BC at r = 0 (i=0): use L'Hôpital → dcdt = 2*D*(c[1]-c[0])/dr²
        dcdt[0] = 2 * D * (c[1] - c[0]) / dr**2
        # Dirichlet BC at r = R (i=N-1)
        dcdt[N - 1] = 0.0
        return dcdt

    return rhs


def diffusion_profile(
    R: float = R_channel,
    D: float = D_Ca,
    c_bath: float = c_Ca_bath,
    t_span: tuple = (0, 300),
    N_r: int = 50,
    n_t: int = 200,
):
    """
    Solve 1-D radial Ca²⁺ diffusion into a cylindrical filament.

    Returns
    -------
    r      : radial node positions [m]
    t      : time array [s]
    c      : concentration field, shape (N_r, n_t) [mol/m³]
    """
    r_nodes = np.linspace(0, R, N_r)
    t_eval = np.linspace(t_span[0], t_span[1], n_t)

    c0 = np.zeros(N_r)
    c0[-1] = c_bath  # initial surface concentration

    rhs = build_mol_rhs(r_nodes, D, c_bath)
    sol = solve_ivp(rhs, t_span, c0, t_eval=t_eval, method="BDF", rtol=1e-4, atol=1e-6)

    # Enforce Dirichlet BC at every time step
    sol.y[-1, :] = c_bath

    return r_nodes, sol.t, sol.y  # sol.y shape: (N_r, n_t)


def run(R: float = R_channel, t_end: float = 300.0):
    """Return diffusion profile dict."""
    r, t, c = diffusion_profile(R=R, t_span=(0, t_end))
    c_centre_final = c[0, -1]
    return {"r": r, "t": t, "c": c, "c_centre_final": c_centre_final}


def plot(save_path: str = "figures/layer4_diffusion.png"):
    """
    Radial Ca2+ concentration profiles at t = 10, 30, 60, 120, 300 s.
    Shows the crosslinking front propagating from the filament surface inward.
    """
    import matplotlib.pyplot as plt

    result = run(t_end=300.0)
    r = result["r"]
    t = result["t"]
    c = result["c"]

    t_display = [10, 30, 60, 120, 300]
    cmap = plt.cm.Blues

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, t_target in enumerate(t_display):
        idx = int(np.argmin(np.abs(t - t_target)))
        color = cmap(0.3 + 0.7 * i / (len(t_display) - 1))
        ax.plot(r * 1e6, c[:, idx], color=color, linewidth=2,
                label=f"t = {t_target} s")

    ax.axvline(R_channel * 1e6, color="k", linestyle=":", linewidth=1.2,
               alpha=0.6, label="Filament surface")
    ax.set_xlabel("Radial Position [um]", fontsize=12, fontweight="bold")
    ax.set_ylabel("Ca2+ Concentration [mol/m3]", fontsize=12, fontweight="bold")
    ax.set_title(
        "Layer 4: Ca2+ Radial Diffusion\n"
        "Crosslinking front propagation into the filament (Fick's 2nd Law, cylindrical)",
        fontsize=14, fontweight="bold",
    )
    ax.set_xlim(0, R_channel * 1e6)
    ax.set_ylim(0, c_Ca_bath * 1.05)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, labelsize=11)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    result = run()
    print(f"Ca2+ concentration at filament centre after {result['t'][-1]:.0f} s: "
          f"{result['c_centre_final']:.2f} mol/m3  "
          f"(bath = {c_Ca_bath:.1f} mol/m3)")
    plot()
