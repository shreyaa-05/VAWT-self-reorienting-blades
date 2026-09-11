"""
pitch_model.py
--------------
Passive "self-reorienting" blade model.

Physical idea (this is the core mechanism the project abstract asks us to
design/optimize): if the pitch pivot axis is located AFT of the blade's
aerodynamic center, the blade behaves like a wind vane / trailing caster --
any aerodynamic normal force creates a nose-into-the-wind restoring moment
about the pivot. Add a torsional spring (with a preload angle and mechanical
end-stops) and the blade settles, at every azimuthal position, at whatever
pitch angle balances:

    aerodynamic moment about pivot  +  spring moment  =  0

with no motor and no external power input. That equilibrium angle is solved
here as a nonlinear root-find at each azimuth.

Sign / geometry convention
---------------------------
- Chord fraction measured from the leading edge, 0 <= x <= 1.
- x_ac: aerodynamic-center location (0.25 for most thin airfoils).
- x_pivot: pitch-axis location. x_pivot > x_ac  -> self-aligning (stable) blade.
- beta: blade pitch angle relative to the *local relative-wind direction's
  perpendicular* -- practically, alpha = phi - beta, where phi is the local
  flow angle (see dmst.py) and alpha is angle of attack.
- beta0: spring preload / neutral angle.
- (beta_min, beta_max): mechanical end-stop limits.
"""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

from airfoil import AirfoilPolar


@dataclass
class PitchHardware:
    x_pivot: float = 0.32      # pivot location as fraction of chord from LE
    x_ac: float = 0.25         # aerodynamic center, fraction of chord from LE
    k_spring: float = 1.5      # spring "stiffness number": k_spring is compared
                                # directly against q*c^2*(x_pivot-x_ac)*dCn/dalpha,
                                # so it is naturally O(1) at typical operating
                                # dynamic pressures for this rotor scale -- see
                                # optimize.py for the physically sensible search range.
    beta0: float = 0.0         # preload / neutral pitch angle [rad]
    beta_min: float = np.radians(-20.0)
    beta_max: float = np.radians(20.0)


def _residual(beta, phi, q_over_kref, hw: PitchHardware, foil: AirfoilPolar):
    """
    Non-dimensional torque residual about the pivot.

    We non-dimensionalize the spring stiffness by q*c^2 so that k_spring is a
    dimensionless "spring number" independent of operating condition -- this
    is what optimize.py actually searches over. q_over_kref = q*c^2 (dynamic
    pressure * chord^2), representing the aerodynamic torque scale at this
    azimuth/operating point; k_spring is applied directly against
    (beta - beta0) in the same units so both terms are torque-like and
    comparable.
    """
    alpha = phi - beta
    cn, _ = foil.cn_ct(alpha)
    cn = float(cn)
    t_aero = -cn * q_over_kref * (hw.x_pivot - hw.x_ac)   # restoring if x_pivot>x_ac and cn,alpha same sign
    t_spring = -hw.k_spring * (beta - hw.beta0)
    return t_aero + t_spring


def solve_beta(phi, q_over_kref, hw: PitchHardware, foil: AirfoilPolar):
    """
    Solve for AN equilibrium pitch angle beta at one instant, ignoring
    stability/history. Kept for quick sanity checks (see __main__ below).
    For actual rotor performance evaluation use solve_beta_continuation,
    which is stability- and history-aware (see docstring there for why
    that matters).
    """
    lo, hi = hw.beta_min, hw.beta_max
    f_lo = _residual(lo, phi, q_over_kref, hw, foil)
    f_hi = _residual(hi, phi, q_over_kref, hw, foil)

    if f_lo * f_hi > 0:
        return lo if f_lo < 0 else hi

    beta = brentq(_residual, lo, hi, args=(phi, q_over_kref, hw, foil), xtol=1e-6, maxiter=100)
    return beta


def _residual_vec(betas, phi, q_over_kref, hw: PitchHardware, foil: AirfoilPolar):
    """Vectorized torque residual over an array of candidate beta values
    (same physics as `_residual`, but evaluated for the whole grid in one
    numpy call instead of many Python-level scalar calls -- this is the
    dominant cost in solve_beta_continuation, so vectorizing it matters a
    lot for optimization-loop performance)."""
    alpha = phi - betas
    cn, _ = foil.cn_ct(alpha)
    t_aero = -cn * q_over_kref * (hw.x_pivot - hw.x_ac)
    t_spring = -hw.k_spring * (betas - hw.beta0)
    return t_aero + t_spring


def solve_beta_continuation(beta_prev, phi, q_over_kref, hw: PitchHardware,
                             foil: AirfoilPolar, n_grid=41):
    """
    Solve for the equilibrium pitch angle beta, tracking continuity with the
    previous azimuth's solution `beta_prev`.

    WHY THIS MATTERS: a "self-reorienting" blade whose pivot is aft of the
    aerodynamic center is only a *stable* wind-vane if the spring is stiff
    enough. If the spring is soft relative to the aerodynamic torque, the
    torque balance can have an unstable interior root -- physically, the
    blade doesn't sit there, it snaps to one of the two mechanical stops,
    and *which* stop depends on which side of the unstable point it was
    already on (hysteresis). A naive independent root-find at every azimuth
    angle can silently return that unphysical unstable root.

    This function instead finds ALL equilibria (all sign changes of the
    torque residual across the full pitch range, via a coarse grid + bisection
    refinement) and returns whichever one is CLOSEST to beta_prev. Marched
    sequentially around the rotor azimuth, this reproduces the correct
    quasi-static, history-dependent behaviour: the blade tracks a stable
    branch continuously until that branch disappears, then jumps
    (snap-through) to the remaining branch -- exactly the bistable "flapping"
    behaviour a weak-spring self-pitching blade actually exhibits.

    Assumption: pitch dynamics are fast compared to the rotor's rotation
    rate (quasi-static). A full inertia+damping time-domain simulation
    (1-DOF torsional ODE) is the natural next-fidelity extension and is
    noted in the README.

    STABILITY: torque residual f(beta) = T_aero(beta) + T_spring(beta).
    An interior root is only a physically real resting point if f is
    DEcreasing through it (perturb beta up -> net torque goes negative and
    pushes it back down). A root where f INcreases through it is unstable --
    the blade cannot actually sit there; it accelerates away toward whichever
    mechanical stop lies on that side. The two stops themselves are valid
    resting points too, if the torque there still points *into* the stop.
    All of that is accounted for below before picking, among the physically
    real candidates, whichever is closest to beta_prev.
    """
    lo, hi = hw.beta_min, hw.beta_max
    grid = np.linspace(lo, hi, n_grid)
    f = _residual_vec(grid, phi, q_over_kref, hw, foil)

    stable_candidates = []

    # Interior stable roots: f decreasing through the crossing (f[i] > 0 > f[i+1])
    for i in range(len(grid) - 1):
        if f[i] > 0.0 and f[i + 1] < 0.0:
            r = brentq(_residual, grid[i], grid[i + 1],
                       args=(phi, q_over_kref, hw, foil), xtol=1e-7, maxiter=100)
            stable_candidates.append(r)
        elif f[i] == 0.0 and f[i - 1 if i > 0 else 0] >= 0.0:
            stable_candidates.append(grid[i])

    # Mechanical stops are valid resting points if the torque there still
    # points into the stop (nothing else can move it further).
    if f[0] < 0.0:
        stable_candidates.append(lo)
    if f[-1] > 0.0:
        stable_candidates.append(hi)

    if not stable_candidates:
        # Degenerate fallback (shouldn't normally happen): pick whichever
        # stop the torque at the domain edges points toward.
        stable_candidates = [lo if f[0] < 0 else hi]

    stable_candidates = np.array(stable_candidates)
    return float(stable_candidates[np.argmin(np.abs(stable_candidates - beta_prev))])


if __name__ == "__main__":
    # Sanity check: sweep local flow angle phi over a full revolution at fixed
    # q_over_kref and confirm the blade pitch responds smoothly and stays
    # within the stops.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    foil = AirfoilPolar()
    hw = PitchHardware()
    phis = np.linspace(-np.pi, np.pi, 200)
    betas = [np.degrees(solve_beta(p, 1.0, hw, foil)) for p in phis]

    plt.figure(figsize=(6, 4))
    plt.plot(np.degrees(phis), betas)
    plt.xlabel("local flow angle phi (deg)")
    plt.ylabel("equilibrium pitch beta (deg)")
    plt.title("Self-reorienting blade pitch response")
    plt.tight_layout()
    plt.savefig("pitch_response_check.png", dpi=120)
    print("Saved pitch_response_check.png")
