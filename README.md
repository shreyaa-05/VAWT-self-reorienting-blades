
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

## Final Terminal Output and Graphs

The final computational study completed successfully, including the
DMST-based multi-TSR optimization and final TSR verification.

<img width="728" height="468" alt="final_tsr_output" src="https://github.com/user-attachments/assets/2ed1ec3c-a738-4b47-8f7d-510abda9757d" />

<img width="698" height="498" alt="WhatsApp Image 2026-09-12 at 10 48 17 PM" src="https://github.com/user-attachments/assets/10a3fbbd-f8b7-4127-a716-c85569fc65b3" />

<img width="704" height="502" alt="WhatsApp Image 2026-09-12 at 10 49 23 PM" src="https://github.com/user-attachments/assets/f4bcc8b3-74ee-4a01-966e-7e4db8fe9885" />

<img width="703" height="495" alt="WhatsApp Image 2026-09-12 at 10 49 58 PM" src="https://github.com/user-attachments/assets/38cec65e-d179-4a98-8b6d-bd213d6e26bb" />

<img width="696" height="489" alt="WhatsApp Image 2026-09-12 at 10 50 33 PM" src="https://github.com/user-attachments/assets/93d0dd23-2ab9-4138-b851-efa620c63047" />

---
