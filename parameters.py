"""
parameters.py
=============
Model constants for GBM / brain-tissue bioink formulation optimisation.

Bioink system : Aspect Biosystems fibrin–alginate composite
  - Fibrinogen    : 20 mg/mL  (nominal; sweep 10–30 mg/mL)
  - Sodium alginate : 0.5 wt% (nominal; sweep 0.25–1.5 wt%)
  - Genipin       : 0.3 mg/mL (secondary crosslinker, fixed)
  - Crosslinking  : ionic (Ca²⁺ bath) + thrombin (fibrin) + genipin (amine)

Printer         : Aspect Biosystems RX1 microfluidic bioprinter
Target tissue   : brain / GBM, culture period 30 days

PARAMETER STATUS CODES (in comments):
  [LIT]   — value taken directly from peer-reviewed literature (cite before use)
  [EST]   — estimated from related literature; needs experimental validation
  [CALC]  — derived from hardware specs or exact physical constants
  [PLACEHOLDER] — rough order-of-magnitude; must be replaced before publication
"""

import numpy as np

# =============================================================================
# NOMINAL FORMULATION  (Aspect Biosystems reference point)
# =============================================================================

c_fib_nom    = 20.0    # [mg/mL]  Nominal fibrinogen concentration           [LIT]
c_alg_nom    = 0.5     # [wt%]    Nominal sodium alginate concentration       [LIT]
c_genipin    = 0.3e-3  # [kg/m³]  Genipin concentration (= 0.3 mg/mL)        [LIT]
M_genipin    = 226.23e-3  # [kg/mol]  Genipin molecular weight               [CALC]
c_thrombin   = 0.5     # [NIH U/mL]  Thrombin concentration (fibrin gelation) [EST]
                       #  — REPLACE; range typically 0.5–2 U/mL in bioprinting

# =============================================================================
# BIVARIATE POWER-LAW RHEOLOGY  (tau = K * gamma_dot^n)
#
# Both K and n depend on BOTH fibrinogen and alginate concentration.
# Reference point: c_alg = 0.5 wt%, c_fib = 20 mg/mL (nominal formulation).
#
# K model (multiplicative power-law):
#   K(c_alg, c_fib) = K_ref * (c_alg/c_alg_ref)^alpha_K * (c_fib/c_fib_ref)^beta_K
#
# n model (additive linear perturbation, clamped to [n_min, n_max]):
#   n(c_alg, c_fib) = n_ref + dn_dalg*(c_alg - c_alg_ref)
#                           + dn_dfib*(c_fib  - c_fib_ref)
#
# Use K_func() and n_func() (defined below) in downstream layers.
# =============================================================================

# --- Reference values at nominal formulation ---
K_ref    = 0.18    # [Pa·s^n]  Consistency index at (0.5 wt%, 20 mg/mL)  [EST]
                   #  — REPLACE; literature for 0.5% alginate ~0.1–0.3 Pa·s^n
n_ref    = 0.78    # [-]       Flow behaviour index at nominal             [EST]
                   #  — REPLACE; fibrin-alginate at low conc ~0.75–0.85

# --- K scaling exponents ---
alpha_K  = 2.2     # [-]  Alginate exponent in K (dominant contributor)   [EST]
                   #  — REPLACE; pure alginate: 2.0–2.5 typical
beta_K   = 0.40    # [-]  Fibrinogen exponent in K (minor contributor)     [EST]
                   #  — REPLACE; fibrinogen adds viscosity weakly at <30 mg/mL

# --- n scaling coefficients ---
dn_dalg  = -0.06   # [1/wt%]    dn/d(c_alg): more shear-thinning w/ alginate [EST]
                   #  — REPLACE
dn_dfib  = -0.003  # [1/(mg/mL)] dn/d(c_fib): mild shear-thinning w/ fibrin  [EST]
                   #  — REPLACE
n_min    = 0.30    # [-]  Physical lower bound on n                        [CALC]
n_max    = 0.99    # [-]  Physical upper bound on n                        [CALC]

# Reference concentrations (anchor points for bivariate scaling)
c_alg_ref = c_alg_nom   # [wt%]
c_fib_ref = c_fib_nom   # [mg/mL]


def K_func(c_alg: float, c_fib: float) -> float:
    """
    Consistency index K [Pa·s^n] for the fibrin–alginate bioink.

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    """
    return K_ref * (c_alg / c_alg_ref) ** alpha_K * (c_fib / c_fib_ref) ** beta_K


def n_func(c_alg: float, c_fib: float) -> float:
    """
    Flow behaviour index n [-] (clamped to [n_min, n_max]).

    Parameters
    ----------
    c_alg : alginate concentration [wt%]
    c_fib : fibrinogen concentration [mg/mL]
    """
    n = n_ref + dn_dalg * (c_alg - c_alg_ref) + dn_dfib * (c_fib - c_fib_ref)
    return float(np.clip(n, n_min, n_max))


# Discrete lookup at nominal point — retained for backward-compatibility checks
# NOTE: downstream layer files should call K_func / n_func directly.
# The old POWER_LAW_PARAMS dict (single-component alginate) is no longer valid
# for this bicomponent system; it has been removed.

# =============================================================================
# CELL DAMAGE — SHEAR STRESS THRESHOLD
# Neural cells (GBM, astrocytes, neurons) are more fragile than generic lines.
# =============================================================================

tau_crit = 100.0   # [Pa]  Wall shear stress threshold for significant cell damage [LIT]
                   #  Confirmed for GBM / neural cells; range 50–150 Pa in literature.
                   #  Ref: Nair et al. 2009 (Biotechnol. Bioeng.); Blaeser et al. 2016

# =============================================================================
# CALCIUM ION DIFFUSION  (alginate crosslinking front)
# At 0.5 wt% alginate the gel network is open; D_Ca is higher than in dense gels.
# =============================================================================

D_Ca = 5.0e-10     # [m²/s]  Ca²⁺ diffusion coefficient in 0.5% alginate / fibrin [LIT]
                   #  Free diffusion in water ≈ 7.9e-10 m²/s.
                   #  Dense alginate gel (2%) ≈ 3–5e-10 m²/s; open network → 5e-10.
                   #  Ref: Vreeker et al. 1992; Mørch et al. 2006

# =============================================================================
# CA-ALGINATE CROSSLINKING KINETICS  (Arrhenius)
# Low alginate concentration → lower activation energy, faster relative diffusion.
# =============================================================================

Ea_cross = 35000.0  # [J/mol]  Activation energy, Ca-alginate at low [alg]  [EST]
                    #  — REPLACE; pure alginate ~35–45 kJ/mol
                    #  Ref: Draget et al. 2004 (Food Hydrocolloids)

A_cross  = 8.0e5    # [1/s]    Arrhenius pre-exponential factor              [PLACEHOLDER]
                    #  — REPLACE with value fit to gelation time data

R_gas    = 8.314    # [J/(mol·K)]  Universal gas constant                    [CALC]

T_process = 295.15  # [K]  Print temperature ≈ 22 °C (RX1 cartridge, ambient) [EST]
                    #  — REPLACE with measured cartridge temperature

# =============================================================================
# FIBRIN POLYMERISATION KINETICS  (thrombin-catalysed)
# Treated as pseudo-first-order: d(xi_fib)/dt = k_fibrin * (1 - xi_fib)
# Rate depends strongly on thrombin concentration.
# =============================================================================

k_fibrin  = 5.0e-3  # [1/s]  Effective fibrin polymerisation rate at 0.5 U/mL [EST]
                    #  Corresponds to ~50% gelation in ~2.3 min at 37 °C.
                    #  — REPLACE; measure gelation time at your thrombin conc.
                    #  Ref: Weisel & Litvinov 2017 (Blood)

# =============================================================================
# GENIPIN CROSSLINKING  (secondary, amine-reactive)
# Genipin reacts with primary amines in fibrin (Lys residues), slow reaction.
# =============================================================================

k_gnp    = 5.0e-4   # [1/s]  Pseudo-first-order genipin crosslinking rate   [PLACEHOLDER]
                    #  — REPLACE; depends on [genipin], pH, temperature.
                    #  Reaction typically completes over hours–days at 37 °C.
                    #  Ref: Sundararaghavan et al. 2008 (Biotechnol. Bioeng.)

alpha_gnp = 1.5     # [-]   Stiffness multiplier from complete genipin XL    [PLACEHOLDER]
                    #  G'_gnp = G'_base * (1 + alpha_gnp * xi_gnp)
                    #  — REPLACE with rheometry data (0.3 mg/mL genipin dose)

# =============================================================================
# BRAIN / GBM CELL MECHANOSENSING
# Brain cortex: ~100–800 Pa; GBM tumour: often stiffer, ~500–2000 Pa.
# Bioink starts soft (~30 Pa) and stiffens over 30-day culture via ECM deposition.
# =============================================================================

G_prime_opt  = 300.0  # [Pa]  Optimal G' for GBM / neural cells               [LIT]
                      #  Brain cortex stiffness range ~200–500 Pa.
                      #  Ref: Budday et al. 2017 (Ann. Biomed. Eng.);
                      #       Grundy et al. 2021 (APL Bioeng.)

sigma_stiff  = 200.0  # [Pa]  Gaussian half-width of stiffness viability curve  [EST]
                      #  Controls sharpness of mechanosensing response.
                      #  — REPLACE; calibrate to cell spreading / proliferation data

# =============================================================================
# BIVARIATE GEL MODULUS MODEL
#
# Additive contributions from alginate ionic gel + fibrin network:
#   G'(c_alg, c_fib, alpha_gel) = alpha_gel * (A_alg * c_alg^beta_alg
#                                             + A_fib * c_fib^beta_fib)
#
# Calibrated so that at nominal (0.5 wt%, 20 mg/mL, alpha_gel=1):
#   G'_nominal ≈ 30 Pa  (brain-mimicking initial softness)
#
# Verification:
#   A_alg * 0.5^2.5 + A_fib * 20^1.2
#   = 45.2 * 0.1768  +  0.600 * 36.36
#   = 7.99 + 21.82 = 29.8 Pa  ≈ 30 Pa  ✓
# =============================================================================

A_alg    = 45.2    # [Pa / (wt%)^beta_alg]  Alginate modulus pre-factor       [EST]
                   #  — REPLACE with rheometry on Ca-crosslinked 0.5% alginate
                   #  Ref: Draget et al. 1994; Mørch et al. 2006

beta_alg = 2.5     # [-]  Alginate modulus scaling exponent (c^2 – c^3 range) [LIT]
                   #  Ref: Ouwerx et al. 1998 (Polymer Gels Networks)

A_fib    = 0.600   # [Pa / (mg/mL)^beta_fib]  Fibrin network modulus pre-factor [EST]
                   #  — REPLACE with rheometry on 20 mg/mL fibrin gels
                   #  Ref: Ryan et al. 1999 (Biophys. J.);
                   #       Litvinov & Weisel 2017 (Matrix Biol.)

beta_fib = 1.2     # [-]  Fibrin modulus scaling exponent                      [EST]
                   #  Fibrin: ~c^1 – c^2; lower exponent at dilute concentrations.
                   #  — REPLACE with concentration-series rheometry


def G_prime_func(c_alg: float, c_fib: float, alpha_gel: float = 1.0) -> float:
    """
    Predicted storage modulus G' [Pa] for the crosslinked fibrin–alginate gel.

    Parameters
    ----------
    c_alg     : alginate concentration [wt%]
    c_fib     : fibrinogen concentration [mg/mL]
    alpha_gel : gel conversion fraction (0–1); 1 = fully crosslinked
    """
    return alpha_gel * (A_alg * c_alg ** beta_alg + A_fib * c_fib ** beta_fib)


# =============================================================================
# NEURAL ECM DEPOSITION MODEL  (30-day culture period)
#
# Neural / GBM cells deposit fibronectin, laminin, tenascin-C.
# ECM deposition is slow compared to connective tissue (no fibroblasts).
# Model: first-order approach to steady-state ECM density.
#
# Design intent:
#   G'_eff(t=30d) = G'_initial * (1 + alpha_ecm * rho_ecm(30d))
#                 ≈ 30 * (1 + 9.0 * 0.925) ≈ 279 Pa  → approaches G'_opt = 300 Pa
# =============================================================================

k_dep     = 1.0e-6  # [1/s]  Neural ECM deposition rate constant             [PLACEHOLDER]
                    #  Half-life of ECM accumulation ≈ 8 days.
                    #  — REPLACE; measure ECM density (immunostaining / ELISA)
                    #  over culture time for your specific cell line.

alpha_ecm = 9.0     # [-]    ECM stiffness contribution coefficient           [PLACEHOLDER]
                    #  G'_eff = G'_gel * (1 + alpha_ecm * rho_ecm)
                    #  Calibrated so G'_eff(30d) ≈ G'_opt = 300 Pa for nominal formulation.
                    #  — REPLACE with paired ECM density + rheometry measurements.

t_culture = 30 * 24 * 3600  # [s]  Total culture period = 30 days            [LIT]

# =============================================================================
# RX1 MICROFLUIDIC PRINTHEAD GEOMETRY
# Aspect Biosystems RX1 with 300 µm nozzle (coaxial flow chip).
# =============================================================================

D_nozzle  = 300e-6   # [m]  Nozzle inner diameter                            [LIT]
R_channel = D_nozzle / 2  # [m]  = 150 µm nozzle radius                      [CALC]

L_channel = 5e-3     # [m]  Effective nozzle channel length (~5 mm for RX1 chip) [EST]
                     #  — REPLACE with actual chip geometry from Aspect Biosystems

v_print   = 4.0e-3   # [m/s]  Print head translation speed (RX1 spec)        [LIT]

# Volumetric flow rate derived from print speed and nozzle cross-section.
# This is a geometric approximation valid when the filament diameter ≈ nozzle diameter.
Q_print   = v_print * np.pi * R_channel**2   # [m³/s] ≈ 2.83e-10 m³/s        [CALC]
                                              # — verify against RX1 pump readings

# =============================================================================
# CROSSLINKING BATH CONDITIONS
# =============================================================================

T_bath    = 310.15   # [K]   Bath / incubator temperature = 37 °C             [LIT]
c_Ca_bath = 102.0    # [mol/m³]  CaCl₂ bath concentration ≈ 102 mM           [EST]
                     #  Aspect Biosystems protocols typically use ~100 mM CaCl₂.
                     #  — REPLACE with exact value from your crosslinking protocol.

# =============================================================================
# PROCESS WINDOW SWEEP RANGES
# Axes: fibrinogen concentration vs alginate concentration
# (Replaces the generic single-axis alginate sweep)
# =============================================================================

# Fibrinogen axis
c_fib_min    = 10.0   # [mg/mL]  Lower bound                                  [LIT]
c_fib_max    = 30.0   # [mg/mL]  Upper bound                                  [LIT]
c_fib_points = 20     # [-]      Grid resolution along fibrinogen axis

# Alginate axis
c_alg_min    = 0.25   # [wt%]    Lower bound                                  [LIT]
c_alg_max    = 1.50   # [wt%]    Upper bound                                  [LIT]
c_alg_points = 20     # [-]      Grid resolution along alginate axis

# =============================================================================
# COMPOSITE SCORE WEIGHTS
# Weights for the printability / biocompatibility composite score in
# process_window.py.  Sum must equal 1.0.
# Adjust priorities based on experimental importance ranking.
# =============================================================================

w_shear = 0.30   # Weight: shear-induced cell viability (critical for neural cells) [EST]
w_shape = 0.25   # Weight: filament shape fidelity                                  [EST]
w_stiff = 0.25   # Weight: initial substrate stiffness match to G'_opt              [EST]
w_ecm   = 0.20   # Weight: long-term ECM-driven stiffness trajectory                [EST]
                 #  — REPLACE all weights with values derived from DOE / optimisation

# =============================================================================
# UNIVERSAL PHYSICAL CONSTANTS  (exact / well-established)
# =============================================================================

R_gas = 8.314    # [J/(mol·K)]  Universal gas constant                        [CALC]
