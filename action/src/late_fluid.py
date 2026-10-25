"""Conservative averaged late closure and explicit Q decay-rate envelopes.
All heavy Q energy is approximated as cold dust until it becomes Goldstone
radiation. Relativistic intermediate daughters are not tracked; this is an
averaged effective approximation, not a microscopic Boltzmann solution.
"""
from common import *
from model import *
from scipy.integrate import solve_ivp,quad

def run_one(rate,rtol,matching_tau=350):
    traj=np.load(RES/'benchmark.npz');ix=int(np.argmin(abs(traj['tau']-matching_tau)))
    z=traj['states'][:,ix];term=Model().terms(z)
    scale=term['rho_chi'];Q=(term['rho_Q_bare']+term['rho_Q_nonminimal'])/scale
    y0=np.array([1.,Q,z[4]/scale,0.])
    alpha=Model().alpha
    def props(n,y):
        rho=scale*((y[0]+y[1])*np.exp(-3*n)+(y[2]+y[3])*np.exp(-4*n))
        # Algebraic solution 3h(H^2)=rho. Positive expanding branch.
        psi=rho/np.sqrt(9+(alpha*rho)**2)
        return np.sqrt(psi),rho
    def rhs(n,y):
        H,rho=props(n,y);X,Y,R,G=y
        vx=GAMMA/H*X;vq=rate/H*Y
        return np.array([-vx,-vq,vx*np.exp(n),vq*np.exp(n)])
    def q_gone(n,y):return y[1]-1e-20
    q_gone.terminal=True;q_gone.direction=-1
    atol=np.array([1e-13,1e-23,1e-12,1e-18])*(rtol/1e-9)
    first=solve_ivp(rhs,(0,30),y0,method='Radau',rtol=rtol,atol=atol,
        events=q_gone,max_step=.05,dense_output=True)
    if not first.success or not len(first.t_events[0]):raise RuntimeError(first.message)
    end1=first.t[-1];yi=first.y[:,-1].copy();discard=yi[1];yi[1]=0.
    def last_rhs(n,y):
        H,rho=props(n,y);vx=GAMMA/H*y[0]
        return np.array([-vx,0.,vx*np.exp(n),0.])
    def chi_gone(n,y):
        _,rho=props(n,y)
        return scale*y[0]*np.exp(-3*n)/rho-1e-10
    chi_gone.terminal=True;chi_gone.direction=-1
    second=solve_ivp(last_rhs,(end1,30),yi,method='DOP853',rtol=rtol,atol=atol,
        events=chi_gone,max_step=.03,dense_output=True)
    if not second.success or not len(second.t_events[0]):raise RuntimeError(second.message)
    nf=second.t[-1];grid=np.linspace(0,nf,1500)
    states=np.array([first.sol(n) if n<=end1 else second.sol(n) for n in grid]).T
    densities=scale*np.vstack([(states[0]+states[1])*np.exp(-3*grid),states[2]*np.exp(-4*grid),states[3]*np.exp(-4*grid)])
    H=np.array([props(n,y)[0] for n,y in zip(grid,states.T)])
    rho0=props(0,y0)[1];rhof=props(nf,second.y[:,-1])[1]
    def energy_loss(n,sol):
        y=sol(n);matter=scale*(y[0]+y[1])*np.exp(-3*n);radiation=scale*(y[2]+y[3])*np.exp(-4*n)
        return 3*matter+4*radiation
    loss=quad(lambda n:energy_loss(n,first.sol),0,end1,epsabs=rho0*1e-11,epsrel=1e-10)[0]
    loss+=quad(lambda n:energy_loss(n,second.sol),end1,nf,epsabs=rho0*1e-11,epsrel=1e-10)[0]
    yf=second.y[:,-1];rvis=scale*yf[2]*np.exp(-4*nf)
    temp=(30*(rvis*V0*MPL**4)/(np.pi**2*106.75))**.25
    goldratio=yf[3]/yf[2]
    # Collisionless massless Goldstones, no later production. Entropy dilution
    # relative to neutrinos at their decoupling: g_*s(T_nu)=10.75.
    dneff=goldratio*(4/7)*106.75*(10.75/106.75)**(4/3)
    out=dict(Gamma_Q_GeV=rate*np.sqrt(V0)*MPL,rtol=rtol,
        matching_tau=float(traj['tau'][ix]),matching_Q_fraction=Q/(1+Q+y0[2]),
        H_matching_relative=props(0,y0)[0]/z[1]-1,
        Q_gone_N_since_matching=end1,final_N_since_matching=nf,
        total_N_since_inflation_end=z[0]+nf,final_visible_T_GeV=temp,
        Goldstone_to_visible_radiation_ratio=goldratio,
        Delta_Neff_massless_homogeneous_Goldstones=dneff,
        discarded_comoving_Q_fraction=discard,energy_balance_relative=abs(rho0-rhof-loss)/rho0,
        min_component=float(states.min()),
        approximations=['chi and heavy Q averaged as pressureless fluids',
            'Q cascade compressed into one decay rate; relativistic intermediate daughters omitted',
            'massless collisionless Goldstones, fixed g*=g*s=106.75 during reheating',
            'no strings, misalignment, baryogenesis, dark matter, or late massive-Goldstone evolution'])
    return out,dict(N_since_matching=grid,H=H,densities=densities,comoving_states=states)

def run():
    dec=json.loads((RES/'decays.json').read_text());rates=np.array(dec['widths_GeV'])/(np.sqrt(V0)*MPL)
    results=[]
    for label,rate in [('slow_envelope',rates.min()),('fast_envelope',rates.max())]:
        pair=[]
        for tol in (1e-9,1e-11):
            out,traj=run_one(rate,tol);out['label']=label;pair.append(out)
            if tol==1e-11:np.savez_compressed(RES/f'late_{label}.npz',**traj)
        convergence=abs(pair[0]['final_visible_T_GeV']/pair[1]['final_visible_T_GeV']-1)
        results.append(dict(label=label,runs=pair,T_convergence_relative=convergence,
            pass_numerics=bool(convergence<1e-5 and pair[-1]['energy_balance_relative']<1e-7 and pair[-1]['min_component']>=-1e-20)))
    scan=[]
    for tau in (200,275):
        for label,rate in [('slow_envelope',rates.min()),('fast_envelope',rates.max())]:
            out,_=run_one(rate,1e-11,tau);out['label']=label;scan.append(out)
    save('late_fluid.json',dict(cases=results,matching_scan=scan,scope='Conservative effective envelope; not a rigorous bound or a full Q cascade Boltzmann calculation.'))
    print(json.dumps(results,indent=2))
if __name__=='__main__':run()
