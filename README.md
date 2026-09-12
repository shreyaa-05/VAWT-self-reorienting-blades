
# Self-Reorienting Blades for VAWTs — Numerical Design Model

A Python simulation and optimization toolkit for the passive, motor-free
"self-reorienting blade" concept for vertical-axis wind turbines (VAWTs).

The blade is pivoted off its aerodynamic center, restrained by a torsional
spring and mechanical stops, and its pitch angle is determined by the balance
between aerodynamic and spring torques as the rotor rotates — without an
active pitch motor.

---

## What this is (and isn't)

This project is a fast engineering-level surrogate model, not a CFD solver.

The baseline model uses single-streamtube momentum theory together with a
2D airfoil polar using Viterna–Corrigan post-stall extrapolation. This makes
large design-space studies and optimization practical on a laptop.

The model is intended to provide an initial computational design and
hypothesis-testing stage before higher-fidelity CFD simulations and
wind-tunnel testing.

The project progressively extends the baseline model through:

- passive self-reorienting blade pitch;
- single-streamtube rotor aerodynamics;
- single-point optimization;
- controlled synthetic turbulence studies;
- DMST-style aerodynamic modeling;
- dynamic pitch simulation;
- multi-TSR optimization; and
- model-to-model comparison.

CFD, URANS/LES simulations, experimental identification of blade inertia and
damping, and wind-tunnel testing remain future validation stages.

---

## Files

| File | Purpose |
|---|---|
| `airfoil.py` | 2D airfoil polar model with linear pre-stall behavior and Viterna–Corrigan post-stall extrapolation. |
| `pitch_model.py` | Core passive blade pitch model. Solves the blade pitch angle from aerodynamic torque and torsional spring torque, including stability and history-aware behavior. |
| `rotor_model.py` | Couples the passive pitch model with single-streamtube momentum theory to calculate rotor performance such as Cp, thrust, and torque. |
| `optimize.py` | Optimizes pivot location, spring stiffness, preload, and mechanical-stop angle using differential evolution. |
| `main.py` | Runs the baseline workflow including optimization, TSR sweep, turbulence study, plots, and summary report generation. |
| `turbulence_model.py` | Generates a reproducible synthetic turbulence field with controllable intensity and spectral shape. |
| `dmst.py` | Engineering-level Double-Multiple-Streamtube (DMST-style) aerodynamic extension with separate upstream and downstream induction factors. |
| `dynamic_pitch.py` | Time-marching 1-DOF dynamic pitch model including inertia, damping, spring torque, aerodynamic torque, and mechanical stops. |
| `advanced_study.py` | Runs the final integrated study including multi-TSR optimization, dynamic pitch, turbulence response, and model comparison. |
| `requirements.txt` | Python package requirements. |
| `summary_report.txt` | Baseline computational summary. |
| `final_advanced_report.txt` | Final integrated research report. |

---

## Turbulence / renewable-energy extension

The project also connects to the research theme:

**"Exploiting Turbulence for Furthering Engineering."**

`turbulence_model.py` generates a controlled and reproducible synthetic
longitudinal turbulence field using correlated Fourier modes.

The turbulence framework allows the study to:

- vary turbulence intensity and examine its effect on energy conversion;
- compare the self-reorienting rotor response under different turbulence
  conditions;
- vary the turbulence spectral exponent at fixed RMS intensity; and
- investigate whether the distribution of turbulent energy across scales
  influences predicted VAWT performance.

The turbulence model is intentionally a controlled engineering surrogate.
It is not intended to reproduce the full physics of atmospheric turbulence or
replace measured wind data, LES, or RANS simulations.

---

## Running the baseline model

Install the required Python packages:

```bash
py -m pip install -r requirements.txt

---

## Final TSR Verification

The final passive blade configuration was evaluated at four target
tip-speed ratios using three aerodynamic approaches:

- single-streamtube quasi-static model;
- DMST-style quasi-static model; and
- dynamic pitch model.

| TSR | Single-streamtube Cp | DMST-style Cp | Dynamic Cp |
|---:|---:|---:|---:|
| 1.5 | 0.05328 | 0.11464 | 0.12808 |
| 2.0 | 0.11399 | 0.20168 | 0.23835 |
| 2.5 | 0.19990 | 0.19843 | 0.35120 |
| 3.0 | 0.30157 | 0.17184 | 0.39292 |

**Verified mean DMST Cp over the reported TSR sweep: 0.17165**

---

## Final Terminal Output

The final computational study completed successfully, including the
DMST-based multi-TSR optimization and final TSR verification.

![Final TSR verification output](final_tsr_output.JPEG)

<img width="741" height="475" alt="image" src="https://github.com/user-attachments/assets/518c12d9-342d-42b0-84eb-247368129992" />





