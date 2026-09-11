"""dynamic_pitch.py - time-marching passive blade pitch and dynamic Cp."""
from dataclasses import dataclass
import numpy as np
from airfoil import AirfoilPolar

@dataclass
class DynamicPitchParameters:
    inertia: float = 2.0e-4
    damping: float = 2.0e-4
    spring_stiffness: float = 0.15
    beta0: float = 0.0
    beta_min: float = np.radians(-20.0)
    beta_max: float = np.radians(20.0)
    x_ac: float = 0.25

def aerodynamic_pitch_torque(beta, phi, q, chord, x_pivot, x_ac, foil):
    alpha = phi - beta
    cn, _ = foil.cn_ct(alpha)
    return -float(cn) * q * chord**2 * (x_pivot - x_ac)

def simulate_dynamic_pitch(tsr, flow_speed, rho, rotor_radius, chord, x_pivot,
                           params, foil=None, n_steps=720, revolutions=4):
    if foil is None:
        foil = AirfoilPolar()
    omega = tsr * flow_speed / rotor_radius
    dt = (2*np.pi/omega)/n_steps
    total = n_steps*revolutions
    beta = params.beta0
    beta_dot = 0.0
    out = {k: np.empty(total) for k in
           ['theta','beta','beta_dot','phi','alpha','Vrel','q','ct_rot','cn_rot','aero_torque']}
    for i in range(total):
        theta = (omega*i*dt) % (2*np.pi)
        vt = omega*rotor_radius + flow_speed*np.cos(theta)
        vn = flow_speed*np.sin(theta)
        vrel = np.hypot(vt,vn)
        phi = np.arctan2(vn,vt)
        q = 0.5*rho*vrel**2
        ta = aerodynamic_pitch_torque(beta,phi,q,chord,x_pivot,params.x_ac,foil)
        ts = -params.spring_stiffness*(beta-params.beta0)
        beta_dd = (ta+ts-params.damping*beta_dot)/max(params.inertia,1e-12)
        beta_dot += beta_dd*dt
        beta += beta_dot*dt
        if beta < params.beta_min:
            beta=params.beta_min
            if beta_dot<0: beta_dot=0
        elif beta > params.beta_max:
            beta=params.beta_max
            if beta_dot>0: beta_dot=0
        alpha=phi-beta
        cl,cd=foil.cl_cd(alpha)
        ct=cl*np.sin(phi)-cd*np.cos(phi)
        cn=cl*np.cos(phi)+cd*np.sin(phi)
        out['theta'][i]=theta; out['beta'][i]=beta; out['beta_dot'][i]=beta_dot
        out['phi'][i]=phi; out['alpha'][i]=alpha; out['Vrel'][i]=vrel; out['q'][i]=q
        out['ct_rot'][i]=ct; out['cn_rot'][i]=cn; out['aero_torque'][i]=ta
    s=total-n_steps
    result={k:v[s:] for k,v in out.items()}
    result.update(omega=omega,dt=dt,
                  beta_rms_deg=float(np.degrees(np.sqrt(np.mean(result['beta']**2)))),
                  alpha_rms_deg=float(np.degrees(np.sqrt(np.mean(result['alpha']**2)))))
    return result

def _momentum_ct(a):
    if a <= 0.4:
        return 4.0 * a * (1.0 - a)
    return 0.889 - (0.0203 - (a - 0.143) ** 2) / 0.6427


def _dynamic_at_induction(tsr, geo, flow, hw, foil, n_steps, inertia, damping, a):
    local_u = flow.U_inf * (1.0 - a)
    params = DynamicPitchParameters(
        inertia=inertia, damping=damping,
        spring_stiffness=max(0.03, 0.35 * hw.k_spring),
        beta0=hw.beta0, beta_min=hw.beta_min,
        beta_max=hw.beta_max, x_ac=hw.x_ac
    )
    dyn = simulate_dynamic_pitch(
        tsr, local_u, flow.rho, geo.R, geo.chord, hw.x_pivot,
        params, foil, n_steps=n_steps, revolutions=4
    )
    ft = dyn['q'] * geo.chord * dyn['ct_rot']
    fa = dyn['q'] * geo.chord * (
        dyn['ct_rot'] * np.cos(dyn['theta']) +
        dyn['cn_rot'] * np.sin(dyn['theta'])
    )
    torque = np.mean(ft) * geo.R * geo.n_blades * geo.span
    thrust = np.mean(fa) * geo.n_blades * geo.span
    area = 2.0 * geo.R * geo.span
    ct = thrust / (0.5 * flow.rho * flow.U_inf**2 * area)
    omega = dyn['omega']
    power = torque * omega
    cp = power / (0.5 * flow.rho * flow.U_inf**3 * area)
    dyn.update(cp=float(cp), ct=float(ct), torque=float(torque),
               power=float(power), induction=float(a))
    return dyn


def dynamic_rotor_performance(tsr, geo, flow, hw, foil, n_steps=360,
                              inertia=2e-4, damping=2e-4, induction=None):
    """Compute dynamic pitch and Cp with a momentum-consistent induction closure.

    If ``induction`` is None, the axial force from the time-marched dynamic
    solution is iterated against the same actuator-disk/Glauert momentum
    relation used by the baseline rotor model. This prevents the dynamic
    post-processing from reporting an unconstrained Cp above the actuator-disk
    reference simply because induction was omitted.
    """
    if induction is not None:
        return _dynamic_at_induction(tsr, geo, flow, hw, foil, n_steps,
                                     inertia, damping, float(induction))

    from scipy.optimize import brentq

    def residual(a):
        d = _dynamic_at_induction(tsr, geo, flow, hw, foil, n_steps,
                                  inertia, damping, a)
        return d['ct'] - _momentum_ct(a)

    grid = np.linspace(0.0, 0.45, 10)
    vals = [residual(float(a)) for a in grid]
    a_star = None
    for i in range(len(grid)-1):
        if vals[i] == 0 or vals[i] * vals[i+1] < 0:
            a_star = brentq(residual, float(grid[i]), float(grid[i+1]),
                            xtol=2e-3, maxiter=20)
            break
    if a_star is None:
        a_star = float(grid[int(np.argmin(np.abs(vals)))])

    return _dynamic_at_induction(tsr, geo, flow, hw, foil, n_steps,
                                 inertia, damping, a_star)

if __name__ == '__main__':
    from pitch_model import PitchHardware
    class G: R=1.; chord=.15; n_blades=3; span=1.5
    class F: U_inf=8.; rho=1.225
    hw=PitchHardware()
    d=dynamic_rotor_performance(2.0,G,F,hw,AirfoilPolar(),n_steps=240)
    print(f'Dynamic pitch self-test: Cp={d["cp"]:.4f}, pitch RMS={d["beta_rms_deg"]:.3f} deg')
