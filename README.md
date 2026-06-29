# GBM Bioink Formulation Optimization Model

A computational model that predicts whether a fibrin–alginate bioink formulation will successfully print a 3D glioblastoma (GBM) tumor model on the **Aspect Biosystems RX1** microfluidic bioprinter — and whether the printed construct will support GBM cell viability over a 30-day culture period.

![Process Window](figures/process_window.png)

*The white star marks the published Aspect Biosystems operating point (Lee et al. 2019: 20 mg/mL fibrinogen, 0.5% alginate). It sits inside the high-scoring (green) region of the process window.*

---

## Background

Glioblastoma multiforme (GBM) is the most lethal primary brain cancer, with a median survival of 14.6 months. Conventional 2D cell culture models fail to predict drug efficacy because they don't replicate the mechanical microenvironment of brain tissue (~200–500 Pa stiffness). Lee et al. (2019) demonstrated that 3D bioprinted fibrin–alginate constructs using the Aspect Biosystems RX1 produced GBM tumor models with significantly higher drug resistance than 2D cultures — more closely matching clinical reality.

Finding a printable, cell-compatible formulation currently requires iterative bench experiments. This model screens 400 formulations in seconds by chaining 7 physical sub-models, each grounded in peer-reviewed literature.

---

## What the Model Does

Seven physical layers are solved in sequence for each formulation point:

```
Formulation (c_alg, c_fib)
        │
        ▼
Layer 1 — Viscosity          Power Law model → apparent viscosity at nozzle shear rate
        │
        ▼
Layer 2 — Shear Stress       Wall shear stress vs τ_crit = 100 Pa (neural cell damage)
        │
        ▼
Layer 3 — Crosslinking       Arrhenius kinetics → gelation time t_gel
        │
        ▼
Layer 4 — Ca²⁺ Diffusion     Fick's 2nd Law (cylindrical) → crosslinking front propagation
        │
        ▼
Layer 5 — Shape Fidelity     φ = t_spread / t_gel → does filament hold shape before gelling?
        │
        ▼
Layer 6 — Gel Stiffness      G'(c_alg, c_fib) vs G'_opt = 300 Pa (brain tissue target)
        │
        ▼
Layer 7 — ECM Deposition     First-order ODE → G'_eff evolution over 30-day culture
        │
        ▼
Composite Score (0–1) → Process Window Map
```

---

## Key Result

The process window identifies the printable, biocompatible formulation space across:
- **Fibrinogen**: 10–30 mg/mL
- **Alginate**: 0.25–1.5% w/v

The published Aspect/UBC operating point (Lee et al. 2019) falls within the high-scoring region, providing direct validation against real experimental data.

### Nominal formulation checks

| Parameter | Value | Target | Status |
|---|---|---|---|
| G' (initial gel stiffness) | 29.8 Pa | ~30 Pa (brain-mimicking) | ✓ |
| Wall shear stress | 7.3 Pa | < 100 Pa (τ_crit) | ✓ |
| G'_eff at 30 days | 278 Pa | ~300 Pa (G'_opt) | ✓ |
| t_gel | 2.3 s | < t_spread | ✓ |

---

## Layer Figures

| Layer | Figure |
|---|---|
| L1 — Viscosity | ![](figures/layer1_viscosity.png) |
| L2 — Shear Stress | ![](figures/layer2_shear.png) |
| L3 — Kinetics | ![](figures/layer3_kinetics.png) |
| L4 — Ca²⁺ Diffusion | ![](figures/layer4_diffusion.png) |
| L5 — Shape Fidelity | ![](figures/layer5_shape.png) |
| L6 — Stiffness | ![](figures/layer6_stiffness.png) |
| L7 — ECM Evolution | ![](figures/layer7_ecm.png) |

---

## How to Run

```bash
git clone https://github.com/bhwu-1/bioink-model.git
cd bioink_model
pip install -r requirements.txt
python main.py
```

All figures are saved to `figures/` at 200 dpi.

---

## File Structure

```
bioink_model/
├── parameters.py          All physical constants (literature-sourced, GBM-tuned)
├── layer1_viscosity.py    Power Law viscosity model
├── layer2_shear.py        Shear-induced cell damage
├── layer3_kinetics.py     Ca-alginate crosslinking kinetics (Arrhenius)
├── layer4_diffusion.py    Ca²⁺ radial diffusion (Fick's 2nd Law, cylindrical)
├── layer5_shape.py        Filament shape fidelity index
├── layer6_stiffness.py    Gel modulus and mechanosensing score
├── layer7_ecm.py          Neural ECM deposition and stiffness evolution
├── process_window.py      2D formulation sweep and composite score map
├── main.py                Entry point
├── figures/               All output figures (200 dpi)
└── requirements.txt
```

---

## Key Parameters

| Parameter | Value | Source |
|---|---|---|
| Fibrinogen (nominal) | 20 mg/mL | Lee et al. 2019 |
| Alginate (nominal) | 0.5% w/v | Lee et al. 2019 |
| Genipin | 0.3 mg/mL | Lee et al. 2019 |
| Nozzle diameter | 300 µm | Aspect Biosystems RX1 |
| Print speed | 4 mm/s | Aspect Biosystems RX1 |
| τ_crit (neural) | 100 Pa | Nair et al. 2009 |
| G'_opt (brain) | 300 Pa | Budday et al. 2017 |
| D_Ca (open gel) | 5×10⁻¹⁰ m²/s | Mørch et al. 2006 |
| Culture period | 30 days | Lee et al. 2019 |

---

## References

1. Lee, J. et al. (2019). 3D bioprinting of human GBM tumor model. *Bioprinting*, 16, e00053.
2. Abelseth, E. et al. (2019). 3D printing of neural tissues from human pluripotent stem cells. *ACS Biomater. Sci. Eng.*, 5(1), 234–244.
3. Budday, S. et al. (2017). Mechanical characterization of human brain tissue. *Acta Biomaterialia*, 48, 319–340.
4. Nair, K. et al. (2009). Characterization of cell viability during bioprinting. *Biotechnol. Bioeng.*, 101(6), 1257–1266.
5. Mørch, Y.A. et al. (2006). Effect of Ca²⁺, Ba²⁺, and Sr²⁺ on alginate microbeads. *Biomacromolecules*, 7(5), 1471–1480.
6. Draget, K.I. et al. (2004). Alginates from algae. *Food Hydrocolloids*, 18(2), 213–228.
7. Ryan, E.A. et al. (1999). Structural origins of fibrin clot rheology. *Biophys. J.*, 77(5), 2813–2826.
8. Sundararaghavan, H.G. et al. (2008). Genipin-induced changes in collagen gels. *Biotechnol. Bioeng.*, 102(2), 632–643.

---

*Ben Wu · UBC Chemical & Biological Engineering · June 2026*
