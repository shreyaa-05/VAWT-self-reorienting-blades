"""dmst.py - engineering-level Double-Multiple-Streamtube-style VAWT model.

This is a DMST-style surrogate: upstream and downstream sectors have separate
induction factors and each blade sees azimuth-dependent local flow. It is not
validated CFD and is intended for model comparison/design exploration.
"""
from dataclasses import dataclass
import numpy as np
from airfoil import AirfoilPolar
from pitch_model import PitchHardware, solve_beta_continuation

@dataclass
class DMSTResult:
    cp: float; induction_upstream: float; induction_downstream: float
    torque: float; power: float; thetas: np.ndarray; results: list

def momentum_ct(a):
    if a <= 0.4: return 4*a*(1-a)
    return 0.889-(0.0203-(a-0.143)**2)/0.6427

def solve_a(ct):
    if ct<=0: return 0.0
    grid=np.linspace(1e-5,0.6,500)
    return float(grid[np.argmin(np.abs(np.array([momentum_ct(a) for a in grid])-ct))])

def _blade(theta,beta_prev,U_local,omega,geo,hw,foil,rho):
    vt=omega*geo.R+U_local*np.cos(theta)
    vn=U_local*np.sin(theta)
    vrel=np.hypot(vt,vn); phi=np.arctan2(vn,vt); q=.5*rho*vrel**2
    beta=solve_beta_continuation(beta_prev,phi,q*geo.chord**2,hw,foil)
    alpha=phi-beta; cl,cd=foil.cl_cd(alpha)
    ct=cl*np.sin(phi)-cd*np.cos(phi); cn=cl*np.cos(phi)+cd*np.sin(phi)
    axial=ct*np.cos(theta)+cn*np.sin(theta)
    return dict(theta=theta,beta=beta,alpha=alpha,phi=phi,Vrel=vrel,q=q,
                cl=cl,cd=cd,ct_rot=ct,cn_rot=cn,axial_coeff=axial)

def dmst_performance(tsr,geo,flow,hw,foil,n_theta=72,max_iter=20,relax=.5,
                     turbulence=None):
    omega=tsr*flow.U_inf/geo.R
    ths=np.linspace(0,2*np.pi,n_theta,endpoint=False)
    up_mask=(ths<=np.pi/2)|(ths>=3*np.pi/2)
    dn_mask=~up_mask
    up=ths[up_mask]; dn=ths[dn_mask]
    def sector_ct(a,sector,phase):
        Umean=flow.U_inf*(1-a); prev=hw.beta0; axial=[]
        for th in sector:
            U=Umean if turbulence is None else turbulence.local_speed(Umean,th,blade_phase=phase)
            r=_blade(th,prev,U,omega,geo,hw,foil,flow.rho); prev=r['beta']
            axial.append(.5*flow.rho*r['Vrel']**2*geo.chord*r['axial_coeff'])
        thrust_per_blade_span=np.mean(axial)
        return thrust_per_blade_span*geo.n_blades*geo.span/(.5*flow.rho*flow.U_inf**2*(geo.R*geo.span))
    au=ad=.15
    for _ in range(max_iter):
        ctu=sector_ct(au,up,0.0); ctd=sector_ct(ad,dn,0.0)
        nu=(1-relax)*au+relax*np.clip(solve_a(ctu),0,.6)
        nd=(1-relax)*ad+relax*np.clip(solve_a(ctd),0,.6)
        if max(abs(nu-au),abs(nd-ad))<1e-4: au,ad=nu,nd; break
        au,ad=nu,nd
    results=[]; torque_sum=[]; thrust_sum=[]
    for b in range(geo.n_blades):
        phase=2*np.pi*b/geo.n_blades; prev=hw.beta0; br=[]
        for th in ths:
            a=au if ((th<=np.pi/2) or (th>=3*np.pi/2)) else ad
            Umean=flow.U_inf*(1-a)
            U=Umean if turbulence is None else turbulence.local_speed(Umean,th,blade_phase=phase)
            r=_blade(th,prev,U,omega,geo,hw,foil,flow.rho); prev=r['beta']; br.append(r)
        if b==0: results=br
        torque_sum.append(np.mean([r['q']*geo.chord*r['ct_rot'] for r in br])*geo.R)
        thrust_sum.append(np.mean([r['q']*geo.chord*r['axial_coeff'] for r in br]))
    torque=float(np.sum(torque_sum)*geo.span); thrust=float(np.sum(thrust_sum)*geo.span)
    power=torque*omega; area=2*geo.R*geo.span
    cp=power/(.5*flow.rho*flow.U_inf**3*area)
    return DMSTResult(float(cp),float(au),float(ad),torque,power,ths,results)

if __name__=='__main__':
    class G: R=1.; chord=.15; n_blades=3; span=1.5
    class F: U_inf=8.; rho=1.225
    out=dmst_performance(2.,G,F,PitchHardware(),AirfoilPolar(),n_theta=36)
    print(f'DMST-style self-test: Cp={out.cp:.4f}, a_up={out.induction_upstream:.3f}, a_down={out.induction_downstream:.3f}')
