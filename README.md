
# Self-Reorienting Blades for VAWTs — Numerical Design Model

This project investigates a passive self-reorienting blade concept for Vertical Axis Wind Turbines (VAWTs). The objective is to develop a computational model in which blade pitch changes naturally in response to aerodynamic loading, reducing the need for active pitch-control mechanisms.

The study progressively develops the aerodynamic and dynamic models from a baseline single-streamtube formulation to a DMST-style aerodynamic model, followed by a dynamic pitch model incorporating blade inertia, damping, and restoring effects. The final design is evaluated across multiple tip-speed ratios (TSRs) and under controlled turbulence perturbations.

The project focuses on computational modeling and design optimization. The results provide a basis for further validation using CFD simulations, wind-tunnel experiments, and experimental identification of mechanical parameters.

---

## Key Concepts
**Vertical Axis Wind Turbine (VAWT):** A wind turbine in which the rotor axis is perpendicular to the incoming wind.
Self-reorienting blade: A blade whose pitch changes passively in response to aerodynamic loading rather than through an actively controlled actuator.

**Tip-Speed Ratio (TSR):** The ratio between the blade's tangential velocity and the incoming wind velocity.

**Power Coefficient (Cp):** A dimensionless measure of the fraction of available wind power extracted by the turbine.

**Single-streamtube model:** A baseline quasi-static aerodynamic model used to estimate blade forces and turbine performance.

**DMST-style model:** A double-multiple-streamtube-style approach that represents the different aerodynamic conditions experienced by the rotor during its upstream and downstream passages.

**Dynamic pitch:** A model in which blade pitch evolves according to aerodynamic torque together with inertia, damping, and restoring effects.

**Passive optimization:** Optimization of mechanical parameters such as pivot location, spring stiffness, and pitch limits without relying on active pitch control.

**Turbulence surrogate:** Controlled perturbations introduced into the computational model to investigate the response of the passive blade system to changing flow conditions.


## Methodology

The project follows a progressive computational modeling and optimization workflow:

**Airfoil characterization**
Aerodynamic lift and drag behavior is represented through an airfoil polar over the required angle-of-attack range.

**Baseline aerodynamic model**
A single-streamtube quasi-static model is developed to calculate blade forces, torque, and power coefficient.

**DMST-style aerodynamic extension**
The baseline model is extended using separate upstream and downstream streamtube treatment to better represent the changing aerodynamic conditions around the rotor.

**Dynamic pitch modeling**
Blade pitch is modeled dynamically using aerodynamic torque together with blade inertia, damping, and restoring effects.

**Turbulence response analysis**
Controlled turbulence perturbations are introduced to study the stability and response of the passive pitch mechanism.

**Multi-TSR optimization**
Mechanical parameters including pivot location, spring stiffness, and pitch limits are optimized across multiple target TSRs.

**Final verification**
The optimized configuration is evaluated using the single-streamtube, DMST-style, and dynamic models, and the resulting performance is compared.

## Implementation / Model Details

The computational framework is implemented in Python using numerical and scientific-computing libraries.

The main components of the model are:

**Airfoil model:** Provides aerodynamic coefficients used for blade force calculations.

**Pitch model:** Defines the passive blade pitch behavior and mechanical restoring characteristics.

**Single-streamtube model:** Provides the baseline aerodynamic performance prediction.

**DMST-style model:** Separates the rotor into upstream and downstream aerodynamic passages and evaluates their contributions to turbine performance.

**Dynamic pitch model:** Incorporates blade inertia, damping, and restoring torque to determine the time-dependent blade pitch response.

**Turbulence model:** Applies controlled flow perturbations to investigate the robustness of the passive response.

**Optimization routine:** Searches for suitable mechanical parameters across multiple TSR conditions.

**Advanced study:** Combines the aerodynamic, dynamic, turbulence, and optimization components for the final evaluation.

The final workflow therefore connects airfoil aerodynamics → rotor aerodynamic modeling → passive pitch dynamics → turbulence response → multi-TSR optimization → final verification.


## Final TSR Verification

The final computational study completed successfully, including the
DMST-based multi-TSR optimization and final TSR verification.

<img width="728" height="468" alt="final_tsr_output" src="https://github.com/user-attachments/assets/2ed1ec3c-a738-4b47-8f7d-510abda9757d" />


## Final Results and Figures

<img width="698" height="498" alt="WhatsApp Image 2026-09-12 at 10 48 17 PM" src="https://github.com/user-attachments/assets/10a3fbbd-f8b7-4127-a716-c85569fc65b3" />

<img width="704" height="502" alt="WhatsApp Image 2026-09-12 at 10 49 23 PM" src="https://github.com/user-attachments/assets/f4bcc8b3-74ee-4a01-966e-7e4db8fe9885" />

<img width="703" height="495" alt="WhatsApp Image 2026-09-12 at 10 49 58 PM" src="https://github.com/user-attachments/assets/38cec65e-d179-4a98-8b6d-bd213d6e26bb" />

<img width="696" height="489" alt="WhatsApp Image 2026-09-12 at 10 50 33 PM" src="https://github.com/user-attachments/assets/93d0dd23-2ab9-4138-b851-efa620c63047" />

---

## Limitations / Future Work
- The aerodynamic predictions are based on reduced-order models rather than full CFD.

- The turbulence treatment is a controlled computational surrogate rather than measured atmospheric turbulence.

- Blade inertia, damping, and other mechanical parameters require experimental identification for physical implementation.

- Three-dimensional effects, dynamic stall, wake interaction, and detailed viscous flow behavior are not fully captured.

- The passive mechanism requires experimental validation to determine whether the predicted pitch response is achievable in a physical rotor.

## Future Work
- CFD validation using higher-fidelity URANS/LES simulations.

- Wind-tunnel testing of the optimized passive blade configuration.

- Experimental identification of inertia and damping parameters.
 
- Comparison between computational and experimental power curves.

- Further optimization of the passive mechanism under realistic turbulent inflow conditions.

## Reproducibility

The project is implemented in Python and includes the scripts,
requirements, and generated figures required to reproduce the
computational study.

