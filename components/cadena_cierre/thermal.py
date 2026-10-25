"""Continuous condensate-to-radiation transfer; same parameters as C50."""
from common import *
import math
from scipy.integrate import solve_ivp
import c46_chi_reheating_relic as C46
import c35b_pseudo_goldstone_dm_rp3 as C35
from c44i_hierarchy_universality import predict
MPL=C46.MPL
SIGMA_OVER_M=1e-5
def transfer(rho0,gamma,rtol=2e-9):
    def rhs(N,z):
        rc,rr=np.exp(z); H=np.sqrt((rc+rr)/3)
        q=gamma/H
        return [-3-q,-4+q*np.exp(z[0]-z[1])]
    def end(N,z): return z[0]-z[1]-np.log(1e-10)
    end.terminal=True; end.direction=-1
    sol=solve_ivp(rhs,(0,60),[np.log(rho0),np.log(rho0*1e-20)],method='Radau',
                  rtol=rtol,atol=rtol*.1,events=end,dense_output=True,max_step=.08)
    if not sol.success or not len(sol.t_events[0]): raise RuntimeError(sol.message)
    Nf=sol.t[-1]; rc,rr=np.exp(sol.y[:,-1]); Tf=C46.T_of_rho(rr)
    # Differential continuity residual: analytic d(rho_total)/dN +3rc+4rr.
    grid=np.linspace(.01,Nf-.01,1200); z=sol.sol(grid); rcg,rrg=np.exp(z)
    dz=np.array([rhs(n,s) for n,s in zip(grid,z.T)]).T
    res=(rcg*dz[0]+rrg*dz[1]+3*rcg+4*rrg)/(rcg+rrg)
    # Independent integrated balance for comoving matter+rad: integrate derivative of total physical energy.
    from scipy.integrate import quad
    integral=quad(lambda n: float(3*np.exp(sol.sol(n)[0])+4*np.exp(sol.sol(n)[1])),0,Nf,
                  epsabs=rho0*rtol,epsrel=rtol,limit=250)[0]
    balance=abs((rho0*(1+1e-20)-rc-rr-integral)/rho0)
    return dict(N_final=Nf,T_final_GeV=Tf,rho_chi_final=rc,rho_rad_final=rr,
                energy_balance_relative=balance,continuity_residual=float(abs(res).max()),
                N_grid=grid,rho_chi=rcg,rho_rad=rrg)
def normalize_approximate_continuous(rtol):
    last={}
    def callback(v0,m,he,re):
        gamma=SIGMA_OVER_M**2*m/(8*np.pi)
        t=transfer(re,gamma,rtol); last.update(t)
        # Interface supplies matching terms, not an instantaneous rho(T) equation.
        return t['T_final_GeV'],re*np.exp(-3*t['N_final'])
    nr=C46.normalize(callback)
    nr['Gamma_MPl']=SIGMA_OVER_M**2*nr['m_chi']/(8*np.pi)
    nr['Gamma_GeV']=nr['Gamma_MPl']*MPL
    nr['F_GeV']=57.735/np.sqrt(.0809)*.30*nr['H_star_GeV']
    nr['f_eff_GeV']=nr['F_GeV']/2
    nr['mu_GeV']=.30*nr['H_star_GeV']
    nr['xi']=.25; nr['sigma_over_m']=SIGMA_OVER_M
    nr['thermal_history']=last
    nr['T_Gamma_convention_GeV']=C46.T_of_rho(3*nr['Gamma_MPl']**2)
    return nr
def normalize_continuous(rtol):
    from scalar_end import endpoint
    Ntrial=50.; history=[]
    for it in range(20):
        pred=predict(.5,2.,.5,10,Ntrial); xstar=pred['chi_star']; en=endpoint(xstar,rtol*.1)
        V0=1.5*np.pi**2*pred['r']*C46.AS/pred['V_star_over_V0']
        m=np.sqrt(C46.chi_numbers()*V0); gamma=SIGMA_OVER_M**2*m/(8*np.pi)
        rho0=V0*en['rho_end_over_V0']; tr=transfer(rho0,gamma,rtol)
        Hstar=np.sqrt(V0*en['H_star_squared_over_V0'])
        H0=C46.H0_over_c*299792.458/3.0857e19*6.582e-25/MPL
        Nmatch=-np.log(C46.KPIV/C46.H0_over_c)+np.log(Hstar/H0)-tr['N_final']+np.log(C46.T0/tr['T_final_GeV'])+np.log(C46.GS0/C46.GSTAR)/3
        delta=Nmatch-en['N_actual']; history.append(dict(N_SR_trial=Ntrial,N_actual=en['N_actual'],N_match=Nmatch,residual=delta))
        if abs(delta)<1e-7: break
        Ntrial+=delta*.98
    if abs(delta)>=1e-7: raise RuntimeError('endpoint/pivot matching did not converge')
    H_end=np.sqrt(rho0/3)
    nr=dict(N_star=en['N_actual'],N_SR_label=Ntrial,n_s=pred['n_s'],r=pred['r'],chi_star=xstar,
        chi_end=en['chi_end'],chi_velocity_end_over_sqrtV0=en['chi_velocity_end'],
        V0=V0,m_chi=m,m_chi_GeV=m*MPL,H_end=H_end,H_end_GeV=H_end*MPL,
        H_star_GeV=Hstar*MPL,rho_end=rho0,T_rh_GeV=tr['T_final_GeV'],Gamma_MPl=gamma,
        Gamma_GeV=gamma*MPL,thermal_history=tr,matching_history=history,
        T_Gamma_convention_GeV=C46.T_of_rho(3*gamma**2),sigma_over_m=SIGMA_OVER_M,xi=.25)
    nr['F_GeV']=57.735/np.sqrt(.0809)*.30*nr['H_star_GeV']; nr['f_eff_GeV']=nr['F_GeV']/2
    nr['mu_GeV']=.30*nr['H_star_GeV']
    return nr
def run():
    def instant(v0,m,he,re):
        rho=min(3*(SIGMA_OVER_M**2*m/(8*np.pi))**2,re)
        return C46.T_of_rho(rho),rho
    old=C46.normalize(instant)
    a=normalize_continuous(2e-8); b=normalize_continuous(2e-10)
    conv={k:abs(a[k]/b[k]-1) for k in ('N_star','V0','T_final_dummy') if k in a}
    conv['T_final_GeV']=abs(a['thermal_history']['T_final_GeV']/b['thermal_history']['T_final_GeV']-1)
    save('thermal.json',dict(instantaneous_reference=old,continuous=b,convergence=conv,
        pass_energy_and_convergence=bool(b['thermal_history']['energy_balance_relative']<1e-6 and max(conv.values())<1e-3),
        assumption='Exact homogeneous GR scalar endpoint epsilon_H=1; slow-roll CMB spectra; subsequent averaged chi fluid w=0, g*=gs=106.75; no lattice preheating or Q backreaction'))
    print('THERMAL',b['N_star'],b['n_s'],b['r'],b['F_GeV'],b['thermal_history']['energy_balance_relative'],flush=True)
    return b
if __name__=='__main__': run()
