# Self-Reorienting Blades for VAWTs — Numerical Design Model

A Python simulation/optimization toolkit for the passive, motor-free
"self-reorienting blade" concept described in the project brief: a VAWT
blade pivoted off its aerodynamic center, restrained by a torsional spring
and mechanical stops, whose pitch angle is set entirely by the balance of
aerodynamic and spring torques as the rotor spins — no motor required.

## What this is (and isn't)

This is a **fast engineering-level surrogate model**, not a CFD solver. It
uses single-streamtube momentum theory (the classical Templin 1974 VAWT
model) plus a 2D airfoil polar with Viterna–Corrigan post-stall
extrapolation. That combination runs a full design in ~1 second, which is
what makes design-space *optimization* (hundreds of evaluations) tractable
on a laptop. It is meant to sit **upstream** of the CFD simulations and
wind-tunnel testing in the project brief — narrowing down promising pivot
location / spring stiffness / stop-angle combinations cheaply, before
committing those candidates to expensive CFD runs or 3D-printed
wind-tunnel models.

## Files

| File | Purpose |
|---|---|
| `airfoil.py` | 2D airfoil polar (lift/drag) with linear pre-stall + Viterna-Corrigan post-stall extrapolation. |
| `pitch_model.py` | The core physics: solves the passive blade's equilibrium pitch angle from the aero-torque/spring-torque balance. Includes a **stability- and history-aware** solver — important because a softly-sprung self-pitching blade can be genuinely bistable (it snaps between its two mechanical stops rather than settling smoothly), and a naive root-find can return an unphysical, unstable equilibrium. |
| `rotor_model.py` | Couples the pitch solver to single-streamtube momentum theory to get full-rotor performance (Cp, thrust, torque) at a given tip-speed ratio. |
| `optimize.py` | Searches pivot location, spring stiffness, preload, and stop angle (via `scipy.optimize.differential_evolution`) to maximize Cp at a chosen design tip-speed ratio. Also provides two benchmarks: a fixed-pitch (classic Darrieus) baseline, and an "ideal motorized pitch" upper bound. |
| `main.py` | Runs everything end-to-end: optimize → sweep tip-speed ratio → run a controlled turbulence/energy-conversion study → plot → write a summary report. |
| `turbulence_model.py` | Reproducible synthetic turbulence field with controllable intensity and spectral shape for studying turbulence effects on energy conversion. |

## Turbulence / renewable-energy extension

The model now also connects directly to the research theme **“Exploiting
Turbulence for Furthering Engineering.”** `turbulence_model.py` generates a
controlled, reproducible synthetic longitudinal turbulence field using
correlated Fourier modes. `main.py` uses it to:

- sweep turbulence intensity and compare the self-reorienting rotor with the
  fixed-pitch baseline;
- vary the turbulence spectral exponent at fixed RMS intensity to explore
  whether *how* turbulence is distributed across scales changes energy
  conversion; and
- save `turbulence_response.png`, `turbulence_spectral_shaping.png`, and
  the corresponding results in `summary_report.txt`.

This is intentionally a **controlled engineering surrogate**, not an LES/RANS
turbulence model. It provides a reproducible hypothesis-testing layer before
higher-fidelity CFD or wind-tunnel measurements.

## Running it

```bash
py -m pip install numpy scipy matplotlib
py main.py

Run the original turbulence extension:

```bash
py main.py --turbulence
```
```

Outputs (written to the working directory):
- `cp_vs_tsr.png` — power coefficient vs. tip-speed ratio for the optimized self-reorienting design, the fixed-pitch baseline, and the ideal-motorized-pitch upper bound.
- `pitch_schedule_optimized.png` — the optimized design's pitch angle, angle of attack, and tangential-force coefficient over one revolution.
- `summary_report.txt` — plain-text summary of the optimized hardware, performance comparison, and turbulence study.
- `turbulence_response.png` — self-reorienting vs fixed-pitch Cp across turbulence intensity.
- `turbulence_spectral_shaping.png` — Cp response to different synthetic turbulence spectral shapes at fixed TI.

Each module also has a small `if __name__ == "__main__":` self-test you can
run individually, e.g. `python3 airfoil.py` to sanity-check the polar, or
`python3 pitch_model.py` to see the pitch-vs-flow-angle response.

## Key modeling choices worth knowing about

1. **Bistability is real, not a bug.** If the spring is too soft relative
   to the aerodynamic torque at the pivot offset chosen, the "equilibrium"
   the physics predicts is unstable, and the blade actually snaps between
   its two stops as the rotor turns. `pitch_model.solve_beta_continuation`
   detects this via a torque-derivative stability check and tracks the
   correct branch as the blade would, given where it was an instant ago —
   this is what produces the flipping pitch schedule you'll see in
   `pitch_schedule_optimized.png`, which closely resembles the *optimal*
   cyclic pitch schedule used by motor-driven cycloturbines.
2. **Optimizing for one design point doesn't guarantee across-the-board
   superiority.** The optimizer in `main.py` targets a single design tip-
   speed ratio. The TSR sweep in the report is deliberately honest about
   where the optimized passive design does and doesn't beat the
   fixed-pitch baseline — see the "Honest caveat" section of
   `summary_report.txt`.
3. **Single-streamtube momentum theory is a coarse aerodynamic closure.**
   It assumes one induction factor for the whole rotor rather than
   azimuth-varying induction (which a Double-Multiple-Streamtube model or
   CFD would capture), and it has no dynamic-stall or blade-wake-interaction
   physics. Treat absolute Cp values as indicative, not authoritative —
   the tool is meant for comparing *design variants* against each other and
   against the two benchmarks, not for predicting exact turbine output.


## Advanced numerical extension

The project now includes a higher-fidelity computational layer that can be
run without physical hardware.

### Dynamic pitch model

`dynamic_pitch.py` adds a time-marching 1-DOF pitch equation containing blade
rotational inertia, damping, torsional spring torque, aerodynamic torque, and
mechanical stops. This tests whether the passive blade can dynamically follow
the changing aerodynamic environment rather than assuming instantaneous
quasi-static equilibrium.

### DMST-style aerodynamic model

`dmst.py` provides an engineering-level Double-Multiple-Streamtube (DMST-style)
extension. The rotor is separated into upstream and downstream sectors with
separate local induction factors. This is intended to capture more of the
azimuthal variation missing from the original single-streamtube model.

It is still a surrogate model and is not a replacement for validated URANS,
LES, or wind-tunnel measurements.

### Multi-TSR optimization

`advanced_study.py` optimizes the passive hardware against multiple target TSR
values rather than only TSR = 2.0. The resulting design is then evaluated
across the target range and compared with the original single-point optimized
design.

### Improved turbulence study

The advanced study applies the reproducible synthetic turbulence framework to
the DMST-style model through controlled quasi-steady sampling. Turbulence
intensity is varied while the response of the improved aerodynamic model is
measured.

### Model comparison

The advanced study produces direct comparisons between:

- single-streamtube and DMST-style aerodynamics;
- single-point and multi-TSR passive optimization;
- quasi-static and dynamic pitch response; and
- baseline and turbulence-exposed predictions.

### Running the advanced study

On Windows:

```bash
py -m py -m pip install numpy scipy matplotlib
py advanced_study.py
```

The advanced run generates:

- `model_comparison_ss_vs_dmst.png`
- `dynamic_pitch_response.png`
- `dmst_turbulence_response.png`
- `multi_tsr_optimization.png`
- `advanced_summary_report.txt`

The original `main.py` remains available for the baseline study and the
original turbulence experiment.

## Research workflow

The project now follows a staged computational workflow:

```text
Physics-based airfoil model
        ↓
Passive self-reorienting pitch
        ↓
Single-streamtube rotor model
        ↓
Single-point optimization
        ↓
Controlled turbulence study
        ↓
DMST-style aerodynamic extension
        ↓
Dynamic pitch simulation
        ↓
Multi-TSR optimization
        ↓
Model-to-model comparison
        ↓
CFD / wind-tunnel validation (future)
```

The final CFD and experimental stages remain future validation work rather than
claims supported by the current numerical surrogate.

## Suggested next steps (see also `summary_report.txt`)

- Swap in a Double-Multiple-Streamtube (DMST) aerodynamic closure for
  azimuth-resolved induction.
- Replace the quasi-static pitch assumption with a 1-DOF torsional
  dynamic (inertia + damping) time-marching simulation, to check whether
  the snap-through transitions are fast enough for quasi-static to be a
  good approximation.
- Broaden the optimization objective to a TSR-range-weighted average if
  the deployment wind regime spans a range of operating conditions.
- Extend the synthetic turbulence model using measured wind-tunnel spectra.
- Validate the turbulence response against CFD and planned wind-tunnel tests.
- Validate the optimized (or nearby) hardware settings against CFD and
  the planned wind-tunnel tests on 3D-printed blades.

## Final integrated research version

The integrated `advanced_study.py` now performs the main extensions as one reproducible numerical workflow:

1. **DMST-style multi-TSR optimization:** the passive pivot, spring stiffness, preload, and mechanical-stop angle are optimized using the DMST-style aerodynamic model itself at TSR = 1.5, 2.0, 2.5, and 3.0.
2. **Dynamic pitch and dynamic Cp:** `dynamic_pitch.py` time-marches a 1-DOF pitch equation containing inertia, damping, spring torque, aerodynamic torque, and mechanical stops. The resulting pitch trajectory is used directly to calculate aerodynamic forces, torque, power, and dynamic Cp.
3. **Improved turbulence experiments:** the synthetic turbulence field is injected into the local blade velocity used by the DMST-style calculation. Turbulence intensity and spectral exponent can therefore be compared while keeping the experiment reproducible.
4. **Model comparison:** the final study compares single-streamtube quasi-static, DMST-style quasi-static, and dynamic-pitch predictions, and compares the original single-point design with the DMST multi-TSR design.

### Run the final integrated study on Windows

From the folder containing `advanced_study.py`:

```powershell
py -m pip install -r requirements.txt
py advanced_study.py
```

The run creates:

- `model_comparison_final.png`
- `dynamic_pitch_response_final.png`
- `dmst_turbulence_response_final.png`
- `dmst_turbulence_spectrum_final.png`
- `final_design_comparison.png`
- `final_advanced_report.txt`

### Final modeling caution

The DMST implementation is an engineering-level, two-sector surrogate rather than validated CFD. The dynamic inertia and damping values are representative numerical parameters rather than measured blade properties. The turbulence field is a reproducible synthetic surrogate rather than measured atmospheric turbulence. These limitations are intentional and are documented so that the numerical study can be used as a design-space and hypothesis-testing tool before CFD and wind-tunnel validation.


### Optimization and verification note

The DMST-based multi-TSR optimizer uses a deliberately cheaper aerodynamic
evaluation during the search. The selected hardware is then reevaluated at
higher resolution for the reported TSR table. The optimization objective and
the verified sweep average are therefore reported separately. Use the verified
TSR-sweep values when discussing final model performance.

The dynamic pitch calculation uses the time-marched pitch trajectory to compute
aerodynamic torque, power, and dynamic Cp. Its inertia and damping are
representative numerical parameters and are not experimentally identified.
