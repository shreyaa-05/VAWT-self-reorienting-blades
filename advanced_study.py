"""advanced_study.py

Final computational study for the passive self-reorienting VAWT concept.
Includes DMST-style multi-TSR optimization, dynamic pitch with Cp, turbulence
experiments, and model comparisons. All models are engineering surrogates and
require CFD/wind-tunnel validation.
"""
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
from airfoil import AirfoilPolar
from pitch_model import PitchHardware
from rotor_model import RotorGeometry, FlowConditions, rotor_performance
from dmst import dmst_performance
from dynamic_pitch import dynamic_rotor_performance
from turbulence_model import SyntheticTurbulence


def make_hw(x):
    return PitchHardware(x_pivot=float(x[0]), k_spring=float(x[1]),
        beta0=np.radians(float(x[2])), beta_min=np.radians(-float(x[3])),
        beta_max=np.radians(float(x[3])))


def objective_dmst(x,geo,flow,foil,tsrs):
    hw=make_hw(x)
    try:
        cps=[dmst_performance(float(t),geo,flow,hw,foil,n_theta=8,max_iter=2).cp for t in tsrs]
        if not np.all(np.isfinite(cps)): return 1e3
        return -float(np.mean(cps))
    except Exception:
        return 1e3


def optimize_dmst_multi_tsr(geo,flow,foil,tsrs):
    bounds=[(.15,.45),(.2,8.0),(-15.,15.),(8.,30.)]
    return differential_evolution(objective_dmst,bounds,args=(geo,flow,foil,np.asarray(tsrs)),
        seed=21,popsize=2,maxiter=1,tol=2e-2,polish=True,workers=1,updating='deferred')


def run_dynamic(hw,geo,flow,foil,tsr=2.0):
    # Representative numerical inertia/damping; not measured hardware values.
    return dynamic_rotor_performance(tsr,geo,flow,hw,foil,n_steps=180,
        inertia=2e-4,damping=2e-4)


def turbulence_study(tsr,geo,flow,foil,hw):
    intensities=np.array([0.,.05,.10,.15,.20])
    cps=[]
    for ti in intensities:
        if ti==0:
            cps.append(dmst_performance(tsr,geo,flow,hw,foil,n_theta=18,max_iter=5).cp)
            continue
        vals=[]
        for seed in (7,17):
            turb=SyntheticTurbulence(intensity=float(ti),n_modes=12,spectral_exponent=1.0,seed=seed)
            out=dmst_performance(tsr,geo,flow,hw,foil,n_theta=18,max_iter=5,turbulence=turb)
            vals.append(out.cp)
        cps.append(float(np.mean(vals)))
    spectra=[.5,1.,1.5,2.]
    spec=[]
    for exponent in spectra:
        vals=[]
        for seed in (7,17):
            turb=SyntheticTurbulence(intensity=.10,n_modes=12,spectral_exponent=exponent,seed=seed)
            vals.append(dmst_performance(tsr,geo,flow,hw,foil,n_theta=18,max_iter=5,turbulence=turb).cp)
        spec.append(float(np.mean(vals)))
    return intensities,np.array(cps),np.array(spectra),np.array(spec)


def main():
    geo=RotorGeometry(R=1.0,chord=.15,n_blades=3,span=1.5)
    flow=FlowConditions(U_inf=8.0,rho=1.225)
    foil=AirfoilPolar()
    tsrs=np.array([1.5,2.0,2.5,3.0])

    print('='*72); print('FINAL ADVANCED VAWT STUDY')
    print('1) DMST-based multi-TSR optimization')
    t0=time.time(); opt=optimize_dmst_multi_tsr(geo,flow,foil,tsrs)
    hw=make_hw(opt.x)
    print(f'Optimization: {time.time()-t0:.1f} s, {opt.nfev} evaluations')
    print(f'  x_pivot={hw.x_pivot:.3f}')
    print(f'  k_spring={hw.k_spring:.3f}')
    print(f'  preload={np.degrees(hw.beta0):.2f} deg')
    print(f'  stops=+-{np.degrees(hw.beta_max):.2f} deg')
    print(f'  mean DMST Cp={-opt.fun:.4f}')

    cp_ss=[]; cp_dmst=[]; cp_dynamic=[]
    for t in tsrs:
        cp_ss.append(rotor_performance(float(t),geo,flow,hw,foil,n_theta=24)['cp'])
        cp_dmst.append(dmst_performance(float(t),geo,flow,hw,foil,n_theta=18,max_iter=5).cp)
        cp_dynamic.append(dynamic_rotor_performance(float(t),geo,flow,hw,foil,n_steps=180)['cp'])
    cp_ss=np.array(cp_ss); cp_dmst=np.array(cp_dmst); cp_dynamic=np.array(cp_dynamic)

    verified_mean=float(np.mean(cp_dmst))
    print('\nFinal TSR verification')
    print('TSR    Single-streamtube Cp    DMST-style Cp    Dynamic Cp')
    for t,a,b,c in zip(tsrs,cp_ss,cp_dmst,cp_dynamic):
        print(f'{t:.1f}    {a:.5f}                 {b:.5f}          {c:.5f}')
    print(f'Verified mean DMST Cp over the reported TSR sweep: {verified_mean:.5f}')

    plt.figure(figsize=(7,5)); plt.plot(tsrs,cp_ss,'o-',label='Single-streamtube quasi-static')
    plt.plot(tsrs,cp_dmst,'s--',label='DMST-style quasi-static')
    plt.plot(tsrs,cp_dynamic,'^-.',label='Dynamic pitch (single-streamtube)')
    plt.xlabel('Tip-speed ratio (TSR)'); plt.ylabel('Power coefficient Cp')
    plt.title('Aerodynamic and pitch-model comparison'); plt.grid(alpha=.3); plt.legend(); plt.tight_layout()
    plt.savefig('model_comparison_final.png',dpi=150); plt.close()

    dyn=run_dynamic(hw,geo,flow,foil,2.0)
    plt.figure(figsize=(7,5)); plt.plot(np.degrees(dyn['theta']),np.degrees(dyn['beta']))
    plt.xlabel('Rotor azimuth (deg)'); plt.ylabel('Pitch angle beta (deg)')
    plt.title('Dynamic passive pitch response at TSR=2.0'); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig('dynamic_pitch_response_final.png',dpi=150); plt.close()

    ti,cp_t,spec,spec_cp=turbulence_study(2.0,geo,flow,foil,hw)
    plt.figure(figsize=(7,5)); plt.plot(ti*100,cp_t,'o-')
    plt.xlabel('Turbulence intensity TI (%)'); plt.ylabel('DMST-style Cp')
    plt.title('DMST-style response to controlled turbulence'); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig('dmst_turbulence_response_final.png',dpi=150); plt.close()
    plt.figure(figsize=(7,5)); plt.plot(spec,spec_cp,'o-')
    plt.xlabel('Turbulence spectral exponent'); plt.ylabel('DMST-style Cp')
    plt.title('Effect of turbulence spectral shape at TI=10%'); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig('dmst_turbulence_spectrum_final.png',dpi=150); plt.close()

    # Compare the original single-point hardware against the final DMST-optimized design.
    old=make_hw([.172,.435,7.47,9.73])
    old_ss=np.array([rotor_performance(float(t),geo,flow,old,foil,n_theta=36)['cp'] for t in tsrs])
    new_dmst=np.array([dmst_performance(float(t),geo,flow,hw,foil,n_theta=18,max_iter=5).cp for t in tsrs])
    plt.figure(figsize=(7,5)); plt.plot(tsrs,old_ss,'o--',label='Original single-point design')
    plt.plot(tsrs,new_dmst,'s-',label='Final DMST multi-TSR design')
    plt.xlabel('Tip-speed ratio (TSR)'); plt.ylabel('Power coefficient Cp')
    plt.title('Original vs DMST multi-TSR optimized design'); plt.grid(alpha=.3); plt.legend(); plt.tight_layout()
    plt.savefig('final_design_comparison.png',dpi=150); plt.close()

    lines=['FINAL ADVANCED VAWT RESEARCH STUDY','='*35,'',
      'The final computational layer directly optimizes the passive hardware with the DMST-style model over multiple TSRs. It also computes dynamic pitch trajectories and dynamic Cp, and applies the synthetic turbulence field inside the DMST-style local-flow calculation.','',
      'DMST MULTI-TSR OPTIMIZATION','-'*28,
      f'Target TSRs: {", ".join(f"{x:.1f}" for x in tsrs)}',f'Pivot location: {hw.x_pivot:.4f}',
      f'Spring number: {hw.k_spring:.4f}',f'Preload: {np.degrees(hw.beta0):.3f} deg',
      f'Stops: +- {np.degrees(hw.beta_max):.3f} deg',
      f'Optimization objective mean DMST Cp (search resolution): {-opt.fun:.5f}',
      'The value above is the objective evaluated during the coarse-resolution optimization search.',
      'The selected design is independently verified below using the reported TSR sweep.',
      'TSR    Single-streamtube Cp    DMST-style Cp    Dynamic Cp']
    for t,a,b,c in zip(tsrs,cp_ss,cp_dmst,cp_dynamic): lines.append(f'{t:.1f}    {a:.5f}                 {b:.5f}          {c:.5f}')
    lines += [f'Verified mean DMST Cp over the reported TSR sweep: {verified_mean:.5f}']
    lines += ['', 'DYNAMIC PITCH','-'*14,
      'The pitch equation was time-marched with blade inertia and damping. The resulting beta(t) was then used directly to calculate lift, drag, tangential force, torque, power, and dynamic Cp.',
      f'Pitch RMS at TSR=2: {dyn["beta_rms_deg"]:.4f} deg',f'AoA RMS at TSR=2: {dyn["alpha_rms_deg"]:.4f} deg',f'Dynamic Cp at TSR=2: {dyn["cp"]:.5f}',
      'Dynamic inertia and damping are representative numerical parameters, not measured hardware properties.','',
      'CONTROLLED TURBULENCE / DMST','-'*27,'TI (%)    DMST-style Cp']
    for a,b in zip(ti*100,cp_t): lines.append(f'{a:5.1f}    {b:.5f}')
    lines += ['', 'SPECTRAL SHAPE AT TI=10%','-'*25,'Exponent    DMST-style Cp']
    for a,b in zip(spec,spec_cp): lines.append(f'{a:8.2f}    {b:.5f}')
    lines += ['', 'INTERPRETATION AND LIMITATIONS','-'*31,
      'The DMST-style model captures sector-level azimuthal induction more explicitly than the single-streamtube model, but it remains an engineering surrogate and is not validated CFD.',
      'The dynamic model establishes a computational comparison between instantaneous-equilibrium and time-marching pitch, but quantitative dynamic parameters require blade mass-property and damping measurements.',
      'The turbulence field is a reproducible synthetic surrogate. It is useful for controlled hypothesis testing but does not reproduce all atmospheric turbulence physics and must be validated against measured spectra.',
      '', 'VALIDATION PATH','-'*15,
      '1. Validate DMST predictions against 2D/3D URANS or LES.',
      '2. Identify inertia and damping experimentally.',
      '3. Replace synthetic turbulence with measured wind-tunnel spectra.',
      '4. Validate Cp and loads using the planned 3D-printed blade wind-tunnel tests.']
    with open('final_advanced_report.txt','w',encoding='utf-8') as f: f.write('\n'.join(lines)+'\n')
    print('\nFinal advanced study complete.')
    print('Saved: model_comparison_final.png, dynamic_pitch_response_final.png, dmst_turbulence_response_final.png, dmst_turbulence_spectrum_final.png, final_design_comparison.png, final_advanced_report.txt')

if __name__=='__main__': main()
