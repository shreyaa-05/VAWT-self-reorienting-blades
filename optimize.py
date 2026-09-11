"""
optimize.py
-----------
Searches the passive-pitch hardware design space:

    x_pivot   : pivot location, fraction of chord from the leading edge
    k_spring  : torsional spring stiffness (see pitch_model.PitchHardware)
    beta0_deg : spring preload / neutral angle [deg]
    stop_deg  : symmetric mechanical stop limit, i.e. beta in [-stop, +stop] [deg]

to maximize the rotor power coefficient Cp at a chosen design tip-speed
ratio, using differential evolution (a global optimizer -- appropriate here
because the bistable/snap-through pitch behaviour makes Cp(design) a
genuinely non-smooth, potentially multi-modal function of the hardware
parameters).

Also provides two benchmarks to put the optimized self-pitching design in
context:
  - "fixed pitch": beta held at a constant angle (no reorientation at all;
    the classic straight-bladed Darrieus baseline).
  - "ideal motorized pitch": at every azimuth, beta is chosen (within the
    same stop limits) to instantaneously maximize local tangential force.
    This is an upper bound representative of Dr. Vijayaraghavan's
    motor-driven optimal-pitch-function approach -- useful for quantifying
    how much performance the *passive* mechanism gives up in exchange for
    not needing a motor.
"""

import numpy as np
from scipy.optimize import differential_evolution

from airfoil import AirfoilPolar
from pitch_model import PitchHardware
from rotor_model import RotorGeometry, FlowConditions, rotor_performance, blade_azimuth_solution


def _make_hw(params):
    x_pivot, k_spring, beta0_deg, stop_deg = params
    return PitchHardware(
        x_pivot=x_pivot,
        k_spring=k_spring,
        beta0=np.radians(beta0_deg),
        beta_min=np.radians(-stop_deg),
        beta_max=np.radians(stop_deg),
    )


# Search bounds. x_pivot spans from just ahead of the aerodynamic center
# (destabilizing / bistable "snap" regime) to well aft of it (stable
# weathervane regime); k_spring spans soft to stiff relative to the
# aerodynamic torque scale at this rotor's operating conditions (see the
# timing/derivation notes in rotor_model.py's __main__ block).
BOUNDS = [
    (0.15, 0.45),     # x_pivot
    (0.1, 8.0),        # k_spring
    (-15.0, 15.0),     # beta0_deg
    (8.0, 30.0),       # stop_deg
]


def objective(params, geo, flow, foil, tsr_design, n_theta):
    hw = _make_hw(params)
    try:
        out = rotor_performance(tsr_design, geo, flow, hw, foil, n_theta=n_theta)
        cp = out["cp"]
        if not np.isfinite(cp):
            return 1.0
        return -cp   # differential_evolution minimizes
    except Exception:
        return 1.0    # penalize any numerical failure rather than crashing the search


def optimize_design(geo: RotorGeometry, flow: FlowConditions, foil: AirfoilPolar,
                     tsr_design=2.0, n_theta=30, popsize=6, maxiter=8, seed=0):
    """
    Runs differential evolution over BOUNDS to maximize Cp at `tsr_design`.
    Returns (best_hw: PitchHardware, best_cp: float, result: OptimizeResult).
    """
    result = differential_evolution(
        objective, BOUNDS,
        args=(geo, flow, foil, tsr_design, n_theta),
        popsize=popsize, maxiter=maxiter, seed=seed,
        tol=1e-3, mutation=(0.5, 1.2), recombination=0.7,
        polish=True, workers=1, updating="deferred",
    )
    best_hw = _make_hw(result.x)
    best_cp = -result.fun
    return best_hw, best_cp, result


def ideal_motorized_pitch_performance(tsr, geo: RotorGeometry, flow: FlowConditions,
                                        foil: AirfoilPolar, stop_deg=30.0, n_theta=90):
    """
    Upper-bound benchmark: at every azimuth theta, choose beta (within
    +-stop_deg) to instantaneously maximize the local tangential-force
    coefficient, as if a motor could set the pitch perfectly and
    instantaneously (no spring, no inertia, no self-reorientation physics).
    Uses the SAME single-streamtube momentum theory as the passive-pitch
    model, so the comparison to the passive design's Cp is apples-to-apples
    on the aerodynamic side; only the pitch-control mechanism differs.
    """
    U_inf, rho = flow.U_inf, flow.rho
    R, c = geo.R, geo.chord
    omega = tsr * U_inf / R
    beta_grid = np.radians(np.linspace(-stop_deg, stop_deg, 121))

    def blade_forces_ideal(U_local, thetas):
        ct_rot = np.empty_like(thetas)
        f_axial = np.empty_like(thetas)
        for i, th in enumerate(thetas):
            vt = omega * R + U_local * np.cos(th)
            vn = U_local * np.sin(th)
            phi = np.arctan2(vn, vt)
            alpha_grid = phi - beta_grid
            cl_g, cd_g = foil.cl_cd(alpha_grid)
            # Rotate (Cl,Cd) [wind frame] by phi into the rotor (tangential,
            # radial) frame -- same corrected resolution as rotor_model.py.
            ct_rot_g = cl_g * np.sin(phi) - cd_g * np.cos(phi)
            cn_rot_g = cl_g * np.cos(phi) + cd_g * np.sin(phi)
            j = int(np.argmax(ct_rot_g))
            ct_rot[i] = ct_rot_g[j]
            f_axial[i] = ct_rot_g[j] * np.cos(th) + cn_rot_g[j] * np.sin(th)
        return ct_rot, f_axial

    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    A = 2 * R * geo.span

    def ct_from_blades(a):
        U_local = U_inf * (1 - a)
        ct_rot, f_axial = blade_forces_ideal(U_local, thetas)
        q_arr = 0.5 * rho * (np.hypot(omega * R + U_local * np.cos(thetas),
                                       U_local * np.sin(thetas))) ** 2
        F_t = q_arr * c * ct_rot
        torque_per_span = np.mean(F_t) * R * geo.n_blades
        F_ax = q_arr * c * f_axial
        thrust_per_span = np.mean(F_ax) * geo.n_blades
        thrust = thrust_per_span * geo.span
        ct = thrust / (0.5 * rho * U_inf ** 2 * A)
        return ct, torque_per_span

    def ct_momentum(a):
        if a <= 0.4:
            return 4 * a * (1 - a)
        return 0.889 - (0.0203 - (a - 0.143) ** 2) / 0.6427

    from scipy.optimize import brentq
    def residual(a):
        ct_b, _ = ct_from_blades(a)
        return ct_b - ct_momentum(a)

    a_lo, a_hi = 1e-4, 0.6
    f_lo, f_hi = residual(a_lo), residual(a_hi)
    if f_lo * f_hi > 0:
        a_grid = np.linspace(a_lo, a_hi, 40)
        res_grid = [residual(a) for a in a_grid]
        a_star = a_grid[int(np.argmin(np.abs(res_grid)))]
    else:
        a_star = brentq(residual, a_lo, a_hi, xtol=1e-4, maxiter=40)

    _, torque_per_span = ct_from_blades(a_star)
    torque = torque_per_span * geo.span
    power = torque * omega
    cp = power / (0.5 * rho * U_inf ** 3 * A)
    return cp


def fixed_pitch_performance(tsr, geo: RotorGeometry, flow: FlowConditions,
                              foil: AirfoilPolar, beta_fixed_deg=0.0, n_theta=90):
    """
    Baseline benchmark: classic fixed-pitch (straight-bladed) Darrieus rotor.
    Implemented as a self-pitch hardware model with an enormous spring
    stiffness, which pins beta == beta0 everywhere (a clean, exact way to
    reuse the same rotor_performance machinery for the fixed-pitch case).
    """
    hw = PitchHardware(x_pivot=0.3, k_spring=1e6, beta0=np.radians(beta_fixed_deg),
                        beta_min=np.radians(beta_fixed_deg - 1e-3),
                        beta_max=np.radians(beta_fixed_deg + 1e-3))
    out = rotor_performance(tsr, geo, flow, hw, foil, n_theta=n_theta)
    return out["cp"]


if __name__ == "__main__":
    import time
    geo = RotorGeometry()
    flow = FlowConditions()
    foil = AirfoilPolar()

    t0 = time.time()
    best_hw, best_cp, result = optimize_design(geo, flow, foil, tsr_design=2.0,
                                                 n_theta=30, popsize=6, maxiter=8)
    print(f"Optimization took {time.time()-t0:.1f} s, {result.nfev} evaluations")
    print(f"Best Cp (design TSR=2.0): {best_cp:.4f}")
    print(f"  x_pivot={best_hw.x_pivot:.3f}, k_spring={best_hw.k_spring:.3f}, "
          f"beta0={np.degrees(best_hw.beta0):.2f} deg, "
          f"stops=+-{np.degrees(best_hw.beta_max):.2f} deg")

    cp_fixed = fixed_pitch_performance(2.0, geo, flow, foil, beta_fixed_deg=0.0)
    cp_ideal = ideal_motorized_pitch_performance(2.0, geo, flow, foil,
                                                   stop_deg=np.degrees(best_hw.beta_max))
    print(f"Fixed-pitch (beta=0) Cp:        {cp_fixed:.4f}")
    print(f"Ideal motorized-pitch Cp:       {cp_ideal:.4f}")
