"""
rotor_model.py
--------------
VAWT rotor performance model.

Aerodynamics: single-streamtube momentum-theory model (Templin, 1974) --
the historically first, and still commonly taught, engineering-level VAWT
model. The whole rotor is treated as one actuator disk; an induction factor
`a` reduces the free-stream velocity to a single "local" velocity that the
blades see throughout the revolution. This is a coarse approximation of the
real, azimuth-varying induction that a full Double-Multiple-Streamtube (DMST)
model or CFD would capture, but it is fast, robust, and adequate for a first
pass design-space search -- which is exactly the role this script plays
relative to the CFD work described in the project abstract.

At each azimuth angle theta, the LOCAL blade pitch is not prescribed -- it is
solved from the passive torque-balance equilibrium in pitch_model.py. That
coupling (aero force -> pitch equilibrium -> aero force depends on new pitch)
is why this is a genuinely "self-reorienting blade" model rather than a
fixed-pitch VAWT model with the pitch angle set to a constant.

Units: SI throughout.
"""

from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import brentq

from airfoil import AirfoilPolar
from pitch_model import PitchHardware, solve_beta_continuation
try:
    from turbulence_model import SyntheticTurbulence
except ImportError:  # optional module for backward compatibility
    SyntheticTurbulence = object


@dataclass
class RotorGeometry:
    R: float = 1.0          # rotor radius [m]
    chord: float = 0.15     # blade chord [m]
    n_blades: int = 3
    span: float = 1.5       # blade span [m] (for reference power/torque only)


@dataclass
class FlowConditions:
    U_inf: float = 8.0      # free-stream wind speed [m/s]
    rho: float = 1.225      # air density [kg/m^3]


def blade_azimuth_solution(theta, beta_prev, U_local, omega, geo: RotorGeometry,
                             hw: PitchHardware, foil: AirfoilPolar, rho):
    """
    For one blade at azimuth `theta` (0 = pointing directly upwind, increasing
    with rotation), given the LOCAL wind speed U_local already seen at the
    rotor (i.e. after the momentum-theory induction has been applied), solve
    the self-pitch equilibrium -- tracking continuity from `beta_prev`, the
    pitch angle at the previous azimuth step (see pitch_model.solve_beta_continuation
    for why that matters) -- and return the resulting (alpha, beta, phi,
    Vrel, Cn, Ct, Cq_local).
    """
    R, c = geo.R, geo.chord

    # Local relative-wind geometry (standard VAWT blade-element kinematics):
    # tangential blade speed omega*R, wind component along/across it varies
    # sinusoidally with azimuth for this simple (single local velocity) model.
    vt = omega * R + U_local * np.cos(theta)
    vn = U_local * np.sin(theta)
    Vrel = np.hypot(vt, vn)
    phi = np.arctan2(vn, vt)   # local flow angle relative to tangential direction

    q = 0.5 * rho * Vrel ** 2
    q_over_kref = q * c ** 2

    beta = solve_beta_continuation(beta_prev, phi, q_over_kref, hw, foil)
    alpha = phi - beta

    # IMPORTANT: Cl, Cd are defined relative to the RELATIVE WIND (lift
    # perpendicular to Vrel, drag along Vrel). phi is the angle of Vrel
    # relative to the rotor's tangential direction e_t, so rotating (Cl,Cd)
    # by phi gives the force directly in the rotor-fixed (tangential e_t,
    # radial e_r) frame -- this is the standard VAWT blade-element
    # resolution (e.g. Islam et al. 2008; Paraschivoiu's DMST formulation).
    #
    #   ct_rot  = Cl*sin(phi) - Cd*cos(phi)   (drives rotation, along e_t)
    #   cn_rot  = Cl*cos(phi) + Cd*sin(phi)   (radial, along e_r)
    #
    # Rotating the CHORD-frame (Cn, Ct) [which are Cl,Cd rotated by ALPHA,
    # not phi] by phi again would double up/mismatch the rotation angle --
    # that was the earlier bug (it flipped signs for some cases and, more
    # seriously, decoupled the thrust used to close the momentum theory from
    # the actual torque, silently violating the Betz limit at high TSR).
    cl, cd = foil.cl_cd(alpha)
    ct_rot = cl * np.sin(phi) - cd * np.cos(phi)
    cn_rot = cl * np.cos(phi) + cd * np.sin(phi)

    return dict(theta=theta, alpha=alpha, beta=beta, phi=phi, Vrel=Vrel,
                cl=cl, cd=cd, ct_rot=ct_rot, cn_rot=cn_rot, q=q)


def _revolution_average(U_local, omega, geo, hw, foil, rho, n_theta=180, n_rev_settle=2, turbulence=None):
    """
    March the self-pitching blade sequentially around the rotor azimuth
    (continuity-following, see solve_beta_continuation), repeating for
    `n_rev_settle` full revolutions so any history-dependence settles into
    a repeating limit cycle, then average the LAST revolution for
    performance. Returns dict with:
      - torque coefficient integrand info (mean tangential force*R -> torque)
      - a normalized thrust-like coefficient CT_local used to close the
        momentum-theory loop.
    """
    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)

    theta_arr = thetas

    # With no turbulence, keep the original single-blade/axisymmetric path
    # exactly as before. With turbulence enabled, each blade receives a
    # phase-shifted realization of the same correlated synthetic field, so
    # the rotor response is no longer artificially identical blade-to-blade.
    blade_results = []
    torque_per_blade = []
    axial_per_blade = []

    for blade_idx in range(geo.n_blades if turbulence is not None else 1):
        blade_phase = 2 * np.pi * blade_idx / geo.n_blades
        beta_prev = hw.beta0
        results = None
        for _rev in range(n_rev_settle):
            results = []
            for th in thetas:
                if turbulence is None:
                    U_theta = U_local
                else:
                    U_theta = turbulence.local_speed(U_local, th, blade_phase=blade_phase)
                r = blade_azimuth_solution(th, beta_prev, U_theta, omega, geo, hw, foil, rho)
                results.append(r)
                beta_prev = r["beta"]

        ct_rot_arr = np.array([r["ct_rot"] for r in results])
        cn_rot_arr = np.array([r["cn_rot"] for r in results])
        q_arr = np.array([r["q"] for r in results])

        # Mean tangential force per unit span per blade: q * c * ct_rot
        F_t = q_arr * geo.chord * ct_rot_arr
        torque_per_blade.append(np.mean(F_t) * geo.R)

        axial_force_coeff_local = (ct_rot_arr * np.cos(theta_arr) +
                                   cn_rot_arr * np.sin(theta_arr))
        F_axial = q_arr * geo.chord * axial_force_coeff_local
        axial_per_blade.append(np.mean(F_axial))

        if blade_idx == 0:
            blade_results = results

    total_torque_per_span = float(np.sum(torque_per_blade))
    total_axial_per_span = float(np.sum(axial_per_blade))

    return dict(total_torque_per_span=total_torque_per_span,
                total_axial_per_span=total_axial_per_span,
                thetas=theta_arr, results=blade_results)


def solve_induction_factor(omega, geo: RotorGeometry, flow: FlowConditions,
                            hw: PitchHardware, foil: AirfoilPolar, n_theta=180, turbulence=None):
    """
    Single-streamtube momentum theory closure (Templin model):

        CT_blade(a) = CT_momentum(a) = 4*a*(1-a)     [a <= 0.4, standard actuator disk]
        (Glauert empirical correction used above a=0.4 to avoid the
         unphysical CT decrease predicted by simple momentum theory there)

    CT is defined on the *frontal swept area* A = 2*R*span, referenced to
    U_inf: CT = Thrust / (0.5*rho*U_inf^2*A).

    Returns (a, revolution_data_at_converged_a).
    """
    U_inf, rho = flow.U_inf, flow.rho
    A = 2 * geo.R * geo.span

    def ct_from_blades(a):
        U_local = U_inf * (1 - a)
        rev = _revolution_average(U_local, omega, geo, hw, foil, rho, n_theta, turbulence=turbulence)
        thrust = rev["total_axial_per_span"] * geo.span
        ct = thrust / (0.5 * rho * U_inf ** 2 * A)
        return ct, rev

    def ct_momentum(a):
        if a <= 0.4:
            return 4 * a * (1 - a)
        # Glauert high-thrust empirical correction (smooth continuation)
        return 0.889 - (0.0203 - (a - 0.143) ** 2) / 0.6427

    def residual(a):
        ct_b, _ = ct_from_blades(a)
        return ct_b - ct_momentum(a)

    # Bracket and solve. VAWT induction is typically modest (a in ~0.0-0.3
    # for reasonable tip-speed ratios), so search a safe bracket.
    a_lo, a_hi = 1e-4, 0.6
    f_lo, f_hi = residual(a_lo), residual(a_hi)
    if f_lo * f_hi > 0:
        # fall back to a coarse grid search if brentq bracket fails
        a_grid = np.linspace(a_lo, a_hi, 60)
        res_grid = [residual(a) for a in a_grid]
        idx = int(np.argmin(np.abs(res_grid)))
        a_star = a_grid[idx]
    else:
        a_star = brentq(residual, a_lo, a_hi, xtol=1e-4, maxiter=40)

    _, rev = ct_from_blades(a_star)
    return a_star, rev


def rotor_performance(tsr, geo: RotorGeometry, flow: FlowConditions,
                       hw: PitchHardware, foil: AirfoilPolar, n_theta=180, turbulence=None):
    """
    Compute rotor power coefficient Cp at a given tip-speed ratio TSR = omega*R/U_inf.

    Returns dict with Cp, Ct (thrust coeff), induction factor a, mean torque,
    and the full per-azimuth blade solution (useful for plotting pitch schedules).
    """
    omega = tsr * flow.U_inf / geo.R
    a, rev = solve_induction_factor(omega, geo, flow, hw, foil, n_theta, turbulence=turbulence)

    torque = rev["total_torque_per_span"] * geo.span   # N*m
    power = torque * omega                              # W
    A = 2 * geo.R * geo.span
    cp = power / (0.5 * flow.rho * flow.U_inf ** 3 * A)

    return dict(tsr=tsr, cp=cp, a=a, omega=omega, torque=torque, power=power,
                thetas=rev["thetas"], results=rev["results"])


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    geo = RotorGeometry()
    flow = FlowConditions()
    hw = PitchHardware()
    foil = AirfoilPolar()

    out = rotor_performance(tsr=2.0, geo=geo, flow=flow, hw=hw, foil=foil, n_theta=120)
    print(f"TSR=2.0: Cp={out['cp']:.4f}, induction a={out['a']:.4f}, "
          f"torque={out['torque']:.3f} N*m, power={out['power']:.2f} W")

    thetas_deg = np.degrees(out["thetas"])
    betas_deg = [np.degrees(r["beta"]) for r in out["results"]]
    alphas_deg = [np.degrees(r["alpha"]) for r in out["results"]]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(thetas_deg, betas_deg)
    ax[0].set_xlabel("azimuth theta (deg)"); ax[0].set_ylabel("pitch beta (deg)")
    ax[0].set_title("Self-reorienting pitch schedule")
    ax[1].plot(thetas_deg, alphas_deg)
    ax[1].set_xlabel("azimuth theta (deg)"); ax[1].set_ylabel("angle of attack (deg)")
    ax[1].set_title("Resulting angle of attack")
    fig.tight_layout()
    fig.savefig("rotor_single_case_check.png", dpi=120)
    print("Saved rotor_single_case_check.png")
