"""Actual coupled integrations and independent dense-output continuity check."""
from common import *
from model import *
from scipy.integrate import solve_ivp,quad
import time
CASES=[dict(name='benchmark',eta=1.,gamma=GAMMA,xi=.25),
       dict(name='GR_minimal_limit',eta=1.,gamma=GAMMA,xi=0.,alpha=0.),
       dict(name='Q_stress_1000',eta=1000.,gamma=GAMMA,xi=.25),
       dict(name='transfer_stress',eta=1.,gamma=.05,xi=.25)]
def integrate(config,rtol,coupled=True):
    mod=Model(**{k:v for k,v in config.items() if k!='name'},coupled=coupled)
    z0=mod.initial(); atol=np.full(25,rtol*.01);atol[4]=min(1e-20,rtol*1e-6)
    start=time.monotonic()
    sol=solve_ivp(mod.rhs,(0,350),z0,method='DOP853',rtol=rtol,atol=atol,max_step=.2,dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    grid=np.linspace(0,350,3501);z=sol.sol(grid)
    props=[mod.terms(s) for s in z.T]
    rho=np.array([p['rho'] for p in props]);pressure=np.array([p['pressure'] for p in props])
    C=np.array([p['constraint'] for p in props]);order=np.array([p['order'] for p in props]);qfrac=np.array([(p['rho_Q_bare']+p['rho_Q_nonminimal'])/p['rho'] for p in props])
    # Five-point numerical derivative of interpolated energy, not the RHS formula.
    probes=np.linspace(.01,349.99,700);dh=.001
    def density(t):return mod.terms(sol.sol(t))['rho']
    dr=np.array([(-density(t+2*dh)+8*density(t+dh)-8*density(t-dh)+density(t-2*dh))/(12*dh) for t in probes])
    flux=np.array([3*sol.sol(t)[1]*(mod.terms(sol.sol(t))['rho']+mod.terms(sol.sol(t))['pressure']) for t in probes])
    cont=float(np.max(abs(dr+flux)/(np.abs(dr)+np.abs(flux)+1e-16)))
    hubble_energy=np.array([3*abs(sol.sol(t)[1])*density(t) for t in probes])
    cont_hubble=float(np.max(abs(dr+flux)/(hubble_energy+1e-16)))
    energy_integral=quad(lambda t:3*sol.sol(t)[1]*(mod.terms(sol.sol(t))['rho']+mod.terms(sol.sol(t))['pressure']),
        0,350,epsabs=rho[0]*1e-10,epsrel=1e-9,limit=1500)[0]
    balance=abs(rho[0]-rho[-1]-energy_integral)/rho[0]
    ordered=np.flatnonzero(order>.9)
    out=dict(name=config['name'],rtol=rtol,coupled=coupled,seconds=time.monotonic()-start,nfev=sol.nfev,
        initial_H=z0[1],final_H=z[1,-1],final_N=z[0,-1],final_rho=rho[-1],
        max_constraint_relative=float(np.max(abs(C)/rho)),continuity_finite_difference_relative=cont,
        continuity_finite_difference_over_Hubble_energy=cont_hubble,
        integral_energy_balance_relative=balance,min_QT_margin=min(p['margin'] for p in props),
        min_F_eff=min(p['F_eff'] for p in props),min_radiation=float(z[4].min()),
        max_Q_energy_fraction=float(np.max(qfrac)),first_order_tau=float(grid[ordered[0]]) if len(ordered) else None,
        first_order_N=float(z[0,ordered[0]]) if len(ordered) else None,
        eta=config['eta'],gamma=config['gamma'],xi=config['xi'])
    return out,dict(tau=grid,states=z,rho=rho,pressure=pressure,constraint=C,order=order,Q_fraction=qfrac)
def run():
    outputs=[]
    for config in CASES:
        pair=[]
        for tol in (1e-9,1e-11):
            out,traj=integrate(config,tol);pair.append(out)
            if tol==1e-11:np.savez_compressed(RES/f"{config['name']}.npz",**traj)
            print('COUPLED',out,flush=True)
        diff=abs(pair[0]['final_H']/pair[1]['final_H']-1)
        case=dict(config=config,runs=pair,H_convergence_relative=diff,
            pass_numerics=bool(pair[-1]['max_constraint_relative']<1e-6 and diff<1e-4
                and pair[-1]['min_QT_margin']>0 and pair[-1]['min_F_eff']>0
                and pair[-1]['integral_energy_balance_relative']<1e-7))
        if config['name'] in ('benchmark','Q_stress_1000'):
            spect,traj=integrate(config,1e-11,coupled=False)
            np.savez_compressed(RES/f"{config['name']}_spectator.npz",**traj)
            case['spectator']=spect
            case['backreaction_H_final_relative']=pair[-1]['final_H']/spect['final_H']-1
            case['backreaction_N_final_difference']=pair[-1]['final_N']-spect['final_N']
        outputs.append(case);save('experiments.json',dict(cases=outputs,complete=len(outputs)==len(CASES)))
    derivedF=np.sqrt(2*ZQ)*MUNIT*MPL
    save('units.json',dict(Munit_GeV=MUNIT*MPL,canonical_Q_prefactor_MPl=np.sqrt(C2),
        F_derived_GeV=derivedF,F_previous_GeV=PARAMS['F_GeV'],F_fractional_difference=derivedF/PARAMS['F_GeV']-1,
        U_coefficient_over_V0=UCOEFF,QT_alpha_V0=Model().alpha,gamma_over_sqrtV0=GAMMA,
        note='Small difference in F reflects inherited rounded mass ratios; Q field action and force use unrounded code potential.'))
if __name__=='__main__':run()
