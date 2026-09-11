"""
airfoil.py
----------
Simple 2D airfoil polar model used for the VAWT blade-element calculations.

Pre-stall:  linear lift slope + parabolic drag polar (thin-airfoil-theory-like,
            with an empirical lift-curve-slope reduction to account for
            viscous effects on a real section such as a NACA0018).

Post-stall: Viterna-Corrigan extrapolation, which is the standard technique
            used in wind-turbine BEM codes to extend a 2D polar smoothly out
            to +/-90 degrees of angle of attack.

This is NOT a substitute for the CFD simulations described in the project
abstract -- it is a fast, differentiable surrogate that lets us explore the
self-reorienting-blade design space (spring stiffness, pivot location,
preload, mechanical stops) cheaply before committing to expensive CFD runs.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class AirfoilPolar:
    cl_alpha: float = 5.73        # lift-curve slope [1/rad], ~2*pi reduced for viscosity
    alpha_stall: float = np.radians(12.0)   # static stall angle [rad]
    cd0: float = 0.01             # zero-lift drag coefficient
    k_induced: float = 0.02       # parabolic drag polar coefficient: cd = cd0 + k*cl^2
    cd_max: float = 1.3           # drag coefficient at 90 deg (flat-plate-like), typical 1.0-2.0
    cm_ac: float = 0.0            # pitching moment coefficient about the aerodynamic center
                                   # (~0 for a symmetric section such as NACA0018)

    def __post_init__(self):
        # Pre-compute Viterna-Corrigan coefficients so stall region is continuous
        # with the linear region at alpha_stall.
        a_s = self.alpha_stall
        cl_s = self.cl_alpha * a_s
        cd_s = self.cd0 + self.k_induced * cl_s ** 2

        self.A1 = self.cd_max / 2.0
        self.B1 = self.cd_max
        # Guard against division-by-zero at alpha_stall -> 90 deg edge cases
        self.A2 = (cl_s - self.cd_max * np.sin(a_s) * np.cos(a_s)) * np.sin(a_s) / max(np.cos(a_s) ** 2, 1e-6)
        self.B2 = (cd_s - self.cd_max * np.sin(a_s) ** 2) / max(np.cos(a_s), 1e-6)
        self._cl_stall = cl_s
        self._cd_stall = cd_s

    def cl_cd(self, alpha):
        """Return (Cl, Cd) for angle of attack alpha [rad]. Vectorized over numpy arrays."""
        alpha = np.atleast_1d(np.asarray(alpha, dtype=float))
        sign = np.sign(alpha)
        sign[sign == 0] = 1.0
        aa = np.abs(alpha)

        cl = np.empty_like(aa)
        cd = np.empty_like(aa)

        pre = aa <= self.alpha_stall
        post = ~pre

        # Pre-stall: linear lift, parabolic drag
        cl[pre] = self.cl_alpha * aa[pre]
        cd[pre] = self.cd0 + self.k_induced * cl[pre] ** 2

        # Post-stall: Viterna-Corrigan (formulas are written for positive alpha;
        # apply to |alpha| then re-apply sign to Cl only -- Cd is symmetric)
        a_post = aa[post]
        with np.errstate(divide='ignore', invalid='ignore'):
            cl_post = self.A1 * np.sin(2 * a_post) + self.A2 * np.cos(a_post) ** 2 / np.sin(a_post)
        cl_post = np.nan_to_num(cl_post, nan=self._cl_stall)
        cd_post = self.B1 * np.sin(a_post) ** 2 + self.B2 * np.cos(a_post)

        cl[post] = cl_post
        cd[post] = cd_post

        cl = cl * sign
        if cl.size == 1:
            return float(cl[0]), float(cd[0])
        return cl, cd

    def cn_ct(self, alpha):
        """Normal / tangential (to chord) force coefficients, useful for moment calcs."""
        cl, cd = self.cl_cd(alpha)
        cn = cl * np.cos(alpha) + cd * np.sin(alpha)
        ct = cl * np.sin(alpha) - cd * np.cos(alpha)
        return cn, ct


if __name__ == "__main__":
    # quick self-test / sanity plot of the polar
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    foil = AirfoilPolar()
    alphas = np.radians(np.linspace(-90, 90, 361))
    cl, cd = foil.cl_cd(alphas)

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].plot(np.degrees(alphas), cl)
    ax[0].set_xlabel("alpha (deg)"); ax[0].set_ylabel("Cl"); ax[0].set_title("Lift")
    ax[1].plot(np.degrees(alphas), cd)
    ax[1].set_xlabel("alpha (deg)"); ax[1].set_ylabel("Cd"); ax[1].set_title("Drag")
    fig.tight_layout()
    fig.savefig("airfoil_polar_check.png", dpi=120)
    print("Saved airfoil_polar_check.png")
