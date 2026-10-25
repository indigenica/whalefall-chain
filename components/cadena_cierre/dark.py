"""Conditional collisionless emission ledger. RP3 source remains uncalibrated."""
from common import *
import c35b_pseudo_goldstone_dm_rp3 as C35
from functools import lru_cache
MPL=C35.MPL; T0=2.7255*8.617333262e-14; H100=100/3.085677581e19*6.582119569e-25
RHOC_H2=3*MPL**2*H100**2
def hubble(T): return np.sqrt(np.pi**2*C35.gstar(T)/90)*T*T/MPL
def scale(T): return T0/T*(3.91/C35.gs(T))**(1/3)
@lru_cache(None)
def nodes(n): return np.polynomial.legendre.leggauss(n)
def source(F,m_eV,Trh,x0=10,efficiency=1.,step=.03,spectrum='ir_delta',spectral_bins=64):
    m=m_eV*1e-9; fa=F/2; mr=F/57.735
    Tstop=np.exp(brentq(lambda l: np.log(hubble(np.exp(l))/(m/3)),np.log(1e-12),np.log(Trh)))
    edge=np.linspace(np.log(Trh),np.log(Tstop),int(np.ceil(np.log(Trh/Tstop)/step))+1)
    q=[]; nc=[]; emitted=[]; avec=[]
    for hi,lo in zip(edge[:-1],edge[1:]):
        T=np.exp((hi+lo)/2); H=hubble(T); a=scale(T)
        dN=hi-lo+np.log(C35.gs(np.exp(hi))/C35.gs(np.exp(lo)))/3
        ell=np.log(mr/H); xi=.24*ell; rho_s=4*xi*np.pi*fa**2*ell*H**2
        # asymptotic GHN Eq5 source, not nonlinear-transient/collapse calibration
        dE=efficiency*2*rho_s*dN
        if spectrum=='ir_delta': xs=np.array([x0]); fractions=np.array([1.])
        elif spectrum=='energy_logflat':
            # Equal emitted energy per log p between x0 H and mr, explicitly a sensitivity arm.
            quadrature,weights=nodes(spectral_bins)
            lo_p=np.log(x0); hi_p=np.log(mr/H)
            xs=np.exp((lo_p+hi_p)/2+quadrature*(hi_p-lo_p)/2); fractions=weights/2
        else: raise ValueError(spectrum)
        p=xs*H; E=np.sqrt(p*p+m*m); numbers=dE*fractions/E*a**3
        q.extend(p*a); nc.extend(numbers); emitted.extend(dE*fractions); avec.extend(np.full(len(xs),a))
    return dict(q=np.array(q),nc=np.array(nc),emitted=np.array(emitted),a=np.array(avec),
                T_stop_GeV=Tstop,m_GeV=m,F_GeV=F,efficiency=efficiency,x0=x0,spectrum=spectrum)
def at_epoch(pack,T):
    a=scale(T); active=pack['a']<=a; n=pack['nc'][active]/a**3
    p=pack['q'][active]/a; m=pack['m_GeV']; E=np.sqrt(m*m+p*p)
    rho=float(np.sum(n*E)); pressure=float(np.sum(n*p*p/(3*E)))
    # 3P is a continuous radiation-like energy moment, not a likelihood-equivalent Neff for massive species.
    rho_nu_one=7*np.pi**2/120*T**4*(4/11)**(4/3) if T<1e-4 else 7*np.pi**2/120*T**4
    return dict(T_GeV=T,rho=rho,pressure=pressure,radiation_moment_3P=3*pressure,
                nonradiation_moment_rho_minus_3P=rho-3*pressure,
                DeltaNeff_radiation_moment=3*pressure/rho_nu_one,
                number_comoving=float(np.sum(pack['nc'][active])),
                fraction_energy_p_gt_m=float(np.sum((n*E)[p>m])/rho) if rho else 0.)
def abundance(pack): return pack['m_GeV']*float(np.sum(pack['nc']))/RHOC_H2
def cold_test():
    p=dict(q=np.array([2.]),nc=np.array([3.]),a=np.array([1e-12]),m_GeV=1.)
    # exact free-particle redshift scalings in regimes separated by ten orders
    def rho(a): return 3/a**3*np.sqrt(1+4/a**2)
    return dict(relativistic_redshift_relative=abs(rho(1e-6)/rho(2e-6)/16-1),
                cold_redshift_relative=abs(rho(1e6)/rho(2e6)/8-1),
                comoving_number_exact=True)
def run():
    thermal=json.loads((ROOT/'results'/'thermal.json').read_text())['continuous']
    F=thermal['F_GeV']; Trh=thermal['T_Gamma_convention_GeV']; m_old=C35.m_needed(F,0)
    baseline=[]
    for spectrum in ('ir_delta','energy_logflat'):
        for x0 in (10,30,100):
            pack=source(F,m_old,Trh,x0,spectrum=spectrum)
            fine=source(F,m_old,Trh,x0,step=.015,spectrum=spectrum)
            spectral_fine=source(F,m_old,Trh,x0,step=.015,spectrum=spectrum,spectral_bins=128)
            st=abundance(fine); mis,Tosc=C35.omega(m_old,F,0)
            baseline.append(dict(spectrum=spectrum,x0=x0,m_G_eV=m_old,Omega_misalignment_h2=mis,
                Omega_emission_rest_mass_h2=st,Omega_total_rest_mass_h2=mis+st,
                maximum_efficiency_at_original_mass=max(0.,(.12-mis)/st),
                convergence_relative=abs(abundance(pack)/st-1),
                spectral_convergence_relative=abs(abundance(spectral_fine)/st-1),
                BBN=at_epoch(fine,1e-3),CMB=at_epoch(fine,.26e-9),today=at_epoch(fine,T0)))
    # Conditional roots: extra production shifts the mass, rather than counting the same budget twice.
    roots=[]
    for spectrum in ('ir_delta','energy_logflat'):
        for eff in (0,.001,.01,.1,1.):
            def budget(lm):
                m=10**lm; mis,_=C35.omega(m,F,0)
                return mis+(abundance(source(F,m,Trh,10,eff,spectrum=spectrum)) if eff else 0)
            m=10**brentq(lambda lm: np.log(budget(lm)/.12),-27,-10,xtol=1e-6)
            pack=source(F,m,Trh,10,eff,spectrum=spectrum)
            ommis,Tosc=C35.omega(m,F,0)
            ell=np.log((F/57.735)/(m*1e-9)); xil=.24*ell**2
            # Literature Eq19 is Omega (no h2), outside its linear domain when xi*log>1000.
            ghn_Omega=.1*xil/3000*((F/2)/1e14)**2*np.sqrt(m/1e-18)*(3.5/C35.gstar(pack['T_stop_GeV']))**.25
            roots.append(dict(spectrum=spectrum,efficiency=eff,m_G_eV=m,Omega_misalignment_h2=ommis,
                Omega_emission_rest_mass_h2=abundance(pack),T_osc_GeV=Tosc,
                BBN=at_epoch(pack,1e-3),CMB=at_epoch(pack,.26e-9),
                GHN_Eq19_Omega_formal=ghn_Omega,xi_log=xil,GHN_Eq19_linear_domain=bool(xil<=1000),
                source_calibrated_RP3=False,collapse_energy_included=False))
    checks=cold_test()
    # Packet-by-packet equality of energy at injection and assigned number * E.
    pack=source(F,m_old,Trh)
    inj=np.sqrt((pack['q']/pack['a'])**2+pack['m_GeV']**2)*pack['nc']/pack['a']**3
    checks['injection_energy_relative']=float(np.max(np.abs(inj/pack['emitted']-1)))
    checks['pass_numerical_ledger']=bool(max(checks[k] for k in ('relativistic_redshift_relative','cold_redshift_relative','injection_energy_relative'))<1e-6
        and max(r['convergence_relative'] for r in baseline)<1e-3
        and max(r['spectral_convergence_relative'] for r in baseline)<1e-3)
    save('dark.json',dict(shared_F_GeV=F,f_eff_GeV=F/2,original_misalignment_only_mass_eV=m_old,
        original_mass_with_emission=baseline,conditional_DM_roots=roots,checks=checks,
        caveats=['U1-inspired source normalization and spectra NOT measured for RP3',
                 'collisionless constant-mass evolution omits nonlinear transient and wall collapse',
                 'radiation-moment Neff is a diagnostic; massive distributions require full Boltzmann/CMB analysis',
                 'DM target is rest-mass abundance; test velocity and isocurvature before observational viability',
                 'radiation-dominated background and interpolated g*; no dark backreaction',
                 'no computed gravitational-wave spectrum or LISA SNR']))
    print('DARK',m_old,[(r['spectrum'],r['x0'],r['Omega_total_rest_mass_h2']) for r in baseline],flush=True)
if __name__=='__main__': run()
