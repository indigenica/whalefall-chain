"""Consistent exact scalar endpoint ε_H=1, in V0=1 time units."""
from common import *
from scipy.integrate import solve_ivp
def endpoint(chi_star,rtol=2e-10):
    pot=Potential(1.); v,dv,_=pot.values(chi_star)
    H2=(v+np.sqrt(v*v+2*dv*dv/3))/6
    xd=-dv/(3*np.sqrt(H2))
    def rhs(t,z):
        x,vel,N=z; v,dv,_=pot.values(x); H=np.sqrt((.5*vel*vel+v)/3)
        return [vel,-3*H*vel-dv,H]
    def end(t,z):
        x,vel,N=z; H2=(.5*vel*vel+pot.V(x))/3
        return vel*vel/(2*H2)-1
    end.terminal=True; end.direction=1
    sol=solve_ivp(rhs,(0,300),[chi_star,xd,0],method='DOP853',rtol=rtol,atol=rtol*.01,
                  events=end,max_step=.15)
    if not sol.success or not len(sol.t_events[0]): raise RuntimeError('scalar endpoint absent')
    x,vel,N=sol.y[:,-1]; rho=.5*vel*vel+pot.V(x)
    return dict(chi_end=x,chi_velocity_end=vel,N_actual=N,rho_end_over_V0=rho,H_star_squared_over_V0=H2,
                initial_chi_velocity=xd)
