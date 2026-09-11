"""
main.py
-------
End-to-end driver for the self-reorienting VAWT blade study:

  1. Optimize the passive-pitch hardware (pivot location, spring stiffness,
     preload, mechanical stop angle) to maximize Cp at a design tip-speed
     ratio.
  2. Sweep tip-speed ratio for three designs and plot Cp vs TSR:
       - optimized self-reorienting (passive) blade
       - fixed-pitch baseline (classic straight-bladed Darrieus)
       - ideal motorized pitch (upper-bound benchmark)
  3. Plot the optimized design's pitch and angle-of-attack schedule over
     one revolution at its best operating point.
  4. Print / save a plain-text summary report.

Run:  python3 main.py
Turbulence extension: python3 main.py --turbulence
Outputs (written to the current directory):
  cp_vs_tsr.png
  pitch_schedule_optimized.png
  summary_report.txt
  turbulence_response.png (with --turbulence)
  turbulence_spectral_shaping.png (with --turbulence)
"""

import time
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from airfoil import AirfoilPolar
from pitch_model import PitchHardware
from rotor_model import RotorGeometry, FlowConditions, rotor_performance
from optimize import optimize_design, fixed_pitch_performance, ideal_motorized_pitch_performance
from turbulence_model import SyntheticTurbulence


def sweep_tsr(tsr_values, geo, flow, foil, hw=None, mode="self_pitch",
              beta_fixed_deg=0.0, stop_deg=20.0, n_theta=60):
    cps = []
    for tsr in tsr_values:
        if mode == "self_pitch":
            out = rotor_performance(tsr, geo, flow, hw, foil, n_theta=n_theta)
            cps.append(out["cp"])
        elif mode == "fixed":
            cps.append(fixed_pitch_performance(tsr, geo, flow, foil,
                                                 beta_fixed_deg=beta_fixed_deg, n_theta=n_theta))
        elif mode == "ideal":
            cps.append(ideal_motorized_pitch_performance(tsr, geo, flow, foil,
                                                            stop_deg=stop_deg, n_theta=n_theta))
        else:
            raise ValueError(mode)
    return np.array(cps)



def fixed_pitch_performance_under_turbulence(tsr, geo, flow, foil, turbulence=None,
                                              beta_fixed_deg=0.0, n_theta=36):
    """Fixed-pitch benchmark evaluated with the same turbulence realization."""
    hw_fixed = PitchHardware(
        x_pivot=0.3, k_spring=1e6, beta0=np.radians(beta_fixed_deg),
        beta_min=np.radians(beta_fixed_deg - 1e-3),
        beta_max=np.radians(beta_fixed_deg + 1e-3),
    )
    out = rotor_performance(
        tsr, geo, flow, hw_fixed, foil, n_theta=n_theta, turbulence=turbulence
    )
    return out["cp"]


def run_turbulence_study(tsr, geo, flow, foil, hw, n_theta=36):
    """Controlled turbulence experiment connecting the VAWT model to Project 51491."""
    intensity_values = np.array([0.00, 0.05, 0.10, 0.15, 0.20])
    cp_self, cp_fixed = [], []

    for ti in intensity_values:
        turb = None if ti == 0 else SyntheticTurbulence(
            intensity=float(ti), n_modes=8, spectral_exponent=1.0, seed=7
        )
        self_out = rotor_performance(
            tsr, geo, flow, hw, foil, n_theta=n_theta, turbulence=turb
        )
        fixed_cp = fixed_pitch_performance_under_turbulence(
            tsr, geo, flow, foil, turbulence=turb, n_theta=n_theta
        )
        cp_self.append(self_out["cp"])
        cp_fixed.append(fixed_cp)

    exponents = np.array([0.5, 1.0, 1.5, 2.0])
    cp_shape = []
    for exponent in exponents:
        turb = SyntheticTurbulence(
            intensity=0.10, n_modes=8, spectral_exponent=float(exponent), seed=7
        )
        out = rotor_performance(
            tsr, geo, flow, hw, foil, n_theta=n_theta, turbulence=turb
        )
        cp_shape.append(out["cp"])

    return {
        "intensity_values": intensity_values,
        "cp_self": np.array(cp_self),
        "cp_fixed": np.array(cp_fixed),
        "spectral_exponents": exponents,
        "cp_shape": np.array(cp_shape),
        "shape_intensity": 0.10,
    }

def main():
    parser = argparse.ArgumentParser(description="VAWT self-reorienting blade design study")
    parser.add_argument("--turbulence", action="store_true",
                        help="run the controlled turbulence/energy-conversion extension")
    args = parser.parse_args()

    geo = RotorGeometry(R=1.0, chord=0.15, n_blades=3, span=1.5)
    flow = FlowConditions(U_inf=8.0, rho=1.225)
    foil = AirfoilPolar()   # NACA-0018-like symmetric section

    tsr_design = 2.0

    print("=" * 70)
    print("Optimizing self-reorienting blade hardware "
          f"(design TSR = {tsr_design})...")
    t0 = time.time()
    best_hw, best_cp_design, result = optimize_design(
        geo, flow, foil, tsr_design=tsr_design, n_theta=30, popsize=6, maxiter=8)
    print(f"Done in {time.time()-t0:.1f} s ({result.nfev} evaluations)")
    print(f"  x_pivot   = {best_hw.x_pivot:.3f} (fraction chord from LE, x_ac=0.25)")
    print(f"  k_spring  = {best_hw.k_spring:.3f}")
    print(f"  beta0     = {np.degrees(best_hw.beta0):.2f} deg (preload)")
    print(f"  stops     = +-{np.degrees(best_hw.beta_max):.2f} deg")
    print(f"  Cp at design TSR: {best_cp_design:.4f}")

    print("=" * 70)
    print("Sweeping tip-speed ratio for all three designs...")
    tsr_values = np.linspace(1.0, 3.5, 11)

    t0 = time.time()
    cp_self = sweep_tsr(tsr_values, geo, flow, foil, hw=best_hw, mode="self_pitch", n_theta=45)
    cp_fixed = sweep_tsr(tsr_values, geo, flow, foil, mode="fixed", beta_fixed_deg=0.0, n_theta=60)
    cp_ideal = sweep_tsr(tsr_values, geo, flow, foil, mode="ideal",
                          stop_deg=np.degrees(best_hw.beta_max), n_theta=60)
    print(f"Sweep done in {time.time()-t0:.1f} s")

    # ---- Plot 1: Cp vs TSR for all three designs ----
    plt.figure(figsize=(7, 5))
    plt.plot(tsr_values, cp_self, "o-", label="Self-reorienting (optimized)")
    plt.plot(tsr_values, cp_fixed, "s--", label="Fixed pitch (baseline)")
    plt.plot(tsr_values, cp_ideal, "^:", label="Ideal motorized pitch (upper bound)")
    plt.axhline(16 / 27, color="gray", linewidth=0.8, linestyle=":")
    plt.text(tsr_values[0], 16 / 27 + 0.01, "Betz limit (0.593)", fontsize=8, color="gray")
    plt.xlabel("Tip-speed ratio (TSR)")
    plt.ylabel("Power coefficient Cp")
    plt.title("VAWT performance: passive self-reorienting blade vs. benchmarks")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("cp_vs_tsr.png", dpi=140)
    print("Saved cp_vs_tsr.png")

    # ---- Plot 2: pitch / angle-of-attack schedule at the design TSR ----
    out_design = rotor_performance(tsr_design, geo, flow, best_hw, foil, n_theta=90)
    thetas_deg = np.degrees(out_design["thetas"])
    betas_deg = [np.degrees(r["beta"]) for r in out_design["results"]]
    alphas_deg = [np.degrees(r["alpha"]) for r in out_design["results"]]
    ct_rot = [r["ct_rot"] for r in out_design["results"]]

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    axes[0].plot(thetas_deg, betas_deg, color="tab:blue")
    axes[0].set_ylabel("Pitch beta (deg)")
    axes[0].set_title(f"Optimized self-reorienting blade -- TSR={tsr_design}, "
                       f"Cp={out_design['cp']:.3f}")
    axes[0].grid(alpha=0.3)

    axes[1].plot(thetas_deg, alphas_deg, color="tab:orange")
    axes[1].axhline(12, color="gray", linewidth=0.8, linestyle=":")
    axes[1].axhline(-12, color="gray", linewidth=0.8, linestyle=":")
    axes[1].set_ylabel("Angle of attack (deg)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(thetas_deg, ct_rot, color="tab:green")
    axes[2].axhline(0, color="black", linewidth=0.6)
    axes[2].set_ylabel("Tangential force coeff.")
    axes[2].set_xlabel("Rotor azimuth theta (deg)")
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig("pitch_schedule_optimized.png", dpi=140)
    print("Saved pitch_schedule_optimized.png")

    # ---- Optional turbulence study: Project 51491 connection ----
    turbulence_report = ""
    if args.turbulence:
        print("=" * 70)
        print("Running controlled turbulence study...")
        t0 = time.time()
        turb_study = run_turbulence_study(
            tsr_design, geo, flow, foil, best_hw, n_theta=36
        )
        print(f"Turbulence study done in {time.time()-t0:.1f} s")

        ti_pct = turb_study["intensity_values"] * 100.0
        cp_t_self = turb_study["cp_self"]
        cp_t_fixed = turb_study["cp_fixed"]

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(ti_pct, cp_t_self, "o-", label="Self-reorienting")
        ax.plot(ti_pct, cp_t_fixed, "s--", label="Fixed pitch")
        ax.set_xlabel("Turbulence intensity TI (%)")
        ax.set_ylabel("Power coefficient Cp")
        ax.set_title("VAWT response to controlled turbulence (TSR=2.0)")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig("turbulence_response.png", dpi=140)
        print("Saved turbulence_response.png")

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(turb_study["spectral_exponents"], turb_study["cp_shape"], "o-")
        ax.set_xlabel("Turbulence spectral exponent")
        ax.set_ylabel("Power coefficient Cp")
        ax.set_title("Engineering the turbulence spectrum: VAWT response at TI=10%")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig("turbulence_spectral_shaping.png", dpi=140)
        print("Saved turbulence_spectral_shaping.png")

        turb_rows = "".join(
            f"  {t:>6.1f}        {cs:>10.4f}                {cf:>10.4f}\n"
            for t, cs, cf in zip(ti_pct, cp_t_self, cp_t_fixed)
        )
        shape_rows = "".join(
            f"  {ex:>16.2f} : {cp:.4f}\n"
            for ex, cp in zip(turb_study["spectral_exponents"], turb_study["cp_shape"])
        )
        turbulence_report = f"""
TURBULENCE / ENERGY-CONVERSION STUDY (Project 51491 connection)
------------------------------------------------------------------
The optimized passive rotor was exposed to a controlled, reproducible
synthetic longitudinal turbulence field. This numerical experiment extends
the VAWT study toward the research theme:
"To appreciate flow turbulence and renewable energy. To investigate new ways
to engineer turbulence to improve the performance of energy conversion
systems. A balance between fundamental/curiosity driven and applied
engineering research."

Turbulence-intensity sweep at TSR = {tsr_design}:
  TI (%)        Self-reorienting Cp        Fixed-pitch Cp
{turb_rows}

A second controlled experiment held TI at 10% while changing the spectral
exponent of the synthetic turbulence field. This tests whether the rotor
response depends only on turbulence intensity or also on the distribution
of turbulent energy across scales.

  Spectral exponent : Cp
{shape_rows}

IMPORTANT MODELING NOTE
The turbulence field is a reproducible Fourier-mode surrogate for controlled
numerical experimentation. It is NOT an LES/RANS turbulence model and does
not replace wind-tunnel or CFD validation.
"""
    else:
        print("Turbulence study skipped. Use --turbulence to enable Project 51491 analysis.")

    # ---- Summary report ----
    # Headline comparison at the design point the optimizer actually
    # targeted (TSR = tsr_design), plus a broader-sweep comparison that
    # honestly reports whether the passive design's advantage generalizes.
    idx_design = int(np.argmin(np.abs(tsr_values - tsr_design)))
    cp_self_design = cp_self[idx_design]
    cp_fixed_design = cp_fixed[idx_design]
    cp_ideal_design = cp_ideal[idx_design]
    gap_design = cp_ideal_design - cp_fixed_design
    recovered_design = (100 * (cp_self_design - cp_fixed_design) / gap_design
                         if abs(gap_design) > 1e-9 else float("nan"))

    cp_peak_self = cp_self.max()
    tsr_peak_self = tsr_values[np.argmax(cp_self)]
    cp_peak_fixed = cp_fixed.max()
    tsr_peak_fixed = tsr_values[np.argmax(cp_fixed)]
    cp_peak_ideal = cp_ideal.max()
    tsr_peak_ideal = tsr_values[np.argmax(cp_ideal)]

    beats_fixed_everywhere = bool(np.all(cp_self >= cp_fixed - 1e-9))

    report = f"""
VAWT SELF-REORIENTING BLADE -- DESIGN STUDY SUMMARY
====================================================

MODEL
-----
Aerodynamics : single-streamtube momentum theory (Templin, 1974) with a
               Viterna-Corrigan-extrapolated NACA-0018-like polar
               (cl_alpha={foil.cl_alpha}, alpha_stall={np.degrees(foil.alpha_stall):.1f} deg,
               cd0={foil.cd0}, cd_max={foil.cd_max}).
Rotor        : R={geo.R} m, chord={geo.chord} m, {geo.n_blades} blades, span={geo.span} m.
Flow         : U_inf={flow.U_inf} m/s, rho={flow.rho} kg/m^3.
Pitch physics: quasi-static torque balance between the aerodynamic moment
               about the pitch pivot and a torsional spring, solved with a
               stability- and history-aware continuation solver so that
               genuinely bistable ("snap-through") pitch behaviour is
               captured correctly rather than an unstable equilibrium.

This is a fast engineering-level surrogate for design-space exploration --
NOT a substitute for the CFD simulations and wind-tunnel testing described
in the project brief. Its job is to cheaply narrow the hardware design
space (pivot location, spring stiffness, preload, stop angle) before
committing to expensive CFD / experimental campaigns on the most promising
candidates.

OPTIMIZED HARDWARE (design TSR = {tsr_design})
-----------------------------------------------
  Pivot location (x_pivot)      : {best_hw.x_pivot:.3f}  (fraction of chord from LE; x_ac = 0.25)
  Spring stiffness (k_spring)   : {best_hw.k_spring:.3f}
  Preload (beta0)               : {np.degrees(best_hw.beta0):.2f} deg
  Mechanical stops              : +-{np.degrees(best_hw.beta_max):.2f} deg
  Cp at design TSR              : {best_cp_design:.4f}
  Optimization cost             : {result.nfev} rotor evaluations

PERFORMANCE AT THE DESIGN POINT (TSR = {tsr_design})
-------------------------------------------------------
  Self-reorienting (optimized)  : Cp = {cp_self_design:.4f}
  Fixed pitch (baseline)        : Cp = {cp_fixed_design:.4f}
  Ideal motorized pitch (bound) : Cp = {cp_ideal_design:.4f}

  --> At its design point, the optimized passive mechanism recovers
      approximately {recovered_design:.0f}% of the fixed-pitch-to-ideal-motorized-pitch
      performance gap, with no motor and no external power input.

PERFORMANCE ACROSS THE FULL TSR SWEEP ({tsr_values[0]:.2f} - {tsr_values[-1]:.2f})
-------------------------------------------------------
  Self-reorienting (optimized)  : peak Cp = {cp_peak_self:.4f}  at TSR = {tsr_peak_self:.2f}
  Fixed pitch (baseline)        : peak Cp = {cp_peak_fixed:.4f}  at TSR = {tsr_peak_fixed:.2f}
  Ideal motorized pitch (bound) : peak Cp = {cp_peak_ideal:.4f}  at TSR = {tsr_peak_ideal:.2f}
  Betz limit (reference)        : Cp = {16/27:.4f}

  Honest caveat: this design was optimized for a SINGLE operating point
  (TSR = {tsr_design}). It {"does" if beats_fixed_everywhere else "does NOT"} beat the fixed-pitch baseline at every
  TSR in the sweep -- a passive mechanism tuned for one design point can
  give up performance elsewhere. A multi-point (or TSR-range-weighted)
  objective function is a natural next step if broad-TSR performance
  matters more than peak performance at a single design condition.

{turbulence_report}
NEXT STEPS (toward the full project scope)
-------------------------------------------
  1. Replace the single-streamtube momentum closure with a full
     Double-Multiple-Streamtube (DMST) or 2D/3D URANS CFD model to
     capture azimuth-varying induction, dynamic stall, and blade-wake
     interaction that this surrogate omits.
  2. Replace the quasi-static pitch assumption with a 1-DOF torsional
     dynamic (inertia + damping + spring + aero) time-marching model to
     check whether the snap-through transitions are actually fast enough,
     relative to the rotor's rotation rate, for the quasi-static
     assumption used here to hold.
  3. Broaden the optimization objective to a weighted average of Cp over
     the expected operating TSR range, rather than a single design point,
     if the deployment wind regime spans a range of tip-speed ratios.
  4. Extend the turbulence surrogate to measured wind-tunnel spectra and
     validate the predicted rotor response under controlled turbulence.
  5. Validate against wind-tunnel tests on 3D-printed blades at the
     optimized (or nearby) hardware settings.
"""
    print(report)
    with open("summary_report.txt", "w") as f:
        f.write(report)
    print("Saved summary_report.txt")


if __name__ == "__main__":
    main()
