"""Global extrapolation of the soft portal and one declared bounded completion."""
from common import *
from model import PARAMS
def run():
    sigma=PARAMS['m_chi']*PARAMS['sigma_over_m'];V0=PARAMS['V0'];lam=.13;M=10.
    critical=np.sqrt(4*lam*V0)/sigma
    # Toy visible radial h, Vvis = lambda*h^4/4 + sigma f(chi) h^2/2.
    # After minimizing h: -[max(-sigma f,0)]^2/(4 lambda).
    potential=Potential(V0);rows=[]
    for x in (-1e3,-1e4,-1e5,-1e6,0.,1.,1000.):
        bare=float(potential.V(x))
        linear=bare-max(-sigma*x,0)**2/(4*lam)
        bounded=bare-max(-sigma*M*np.tanh(x/M),0)**2/(4*lam)
        rows.append(dict(chi=x,V_linear_minimized=linear,V_bounded_minimized=bounded))
    coefficient=sigma*sigma*M*M/(4*lam)
    if PLATEAU=='sum':
        # Vchi >= V0/Z*tanh^2(chi/f0) >= V0/Z*tanh^2(chi/M), provided M>=f0.
        first_level=V0/potential.Z; proof=bool(M>=potential.sc[0] and coefficient<first_level); pole_bound=None
    else:
        # U26 C4: Vpole(chi) >= coefficient*tanh^2(chi/M) for all chi. inf of Vpole/tanh^2(chi/M) on |chi| in [1e-6,1e6] MPl
        # (both even) and its chi -> 0 limit V0 M^2/(2 C^2); beyond 1e6 the ratio exceeds V0 (1 - C/1e6).
        xs=np.geomspace(1e-6,1e6,200001); ratio=potential.V(xs)/np.tanh(xs/M)**2
        from scipy.optimize import minimize_scalar   # refine the grid minimum (the grid alone is not a lower bound)
        i=int(ratio.argmin()); f=lambda x: float(potential.V(x)/np.tanh(x/M)**2)
        ref=minimize_scalar(f,bounds=(xs[max(i-1,0)],xs[min(i+1,len(xs)-1)]),method='bounded',options=dict(xatol=1e-12))
        limit0=V0*M*M/(2*potential.C**2); first_level=float(min(ratio.min(),ref.fun,limit0,V0*(1-potential.C/1e6)))
        proof=bool(coefficient<first_level)
        pole_bound=dict(grid='|chi| in [1e-6, 1e6] MPl, 200001 log points',grid_min_ratio=float(ratio.min()),
            chi_at_min=float(xs[ratio.argmin()]),refined_min=float(ref.fun),chi_at_refined_min=float(ref.x),
            limit_chi_to_0=limit0,tail_bound=V0*(1-potential.C/1e6))
    save('portal.json',dict(sigma_MPl=sigma,lambda_visible_toy=lam,M_completion_MPl=M,
        linear_portal_unbounded_below=True,linear_instability_scale_asymptotic=critical,
        bounded_completion='f(chi)=M*tanh(chi/M), a new EFT hypothesis, not a UV derivation',
        negative_correction_bound=coefficient,first_level_lower_bound_coefficient=first_level,
        nonnegative_total_toy_potential_proof=proof,plateau=PLATEAU,pole_lower_bound=pole_bound,
        low_field_derivative_identity=True,
        caveats=['toy visible scalar with positive quartic, not an SM Higgs vacuum-stability analysis',
                 'Gamma is retained as phenomenological dissipative transfer; nonlinear preheating is not derived',
                 'This stabilizes the portal extrapolation, not the hierarchy at large chi'],rows=rows))
    print('PORTAL',critical,'completion lower-bound coefficient',first_level-coefficient)
if __name__=='__main__':run()
