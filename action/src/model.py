"""A declared effective matter action on NPQT homogeneous gravity.
Time tau=sqrt(V0)*t, densities/V0, canonical fields/reduced MPl.
Q is a symmetric 4x4 order parameter, D=Q-X0 I, q=c D.
"""
from common import *
from scipy.optimize import brentq
PARAMS=json.loads(THERMAL.read_text())['continuous']
V0=PARAMS['V0']; MPL=2.435e18; ZQ=400.; LP=192.; LAM=320.
X0=float(np.roots([2*LP,-3*LP,LP+4*LAM,-LAM])[np.abs(np.roots([2*LP,-3*LP,LP+4*LAM,-LAM]).imag)<1e-10].real[0])
MUNIT=(PARAMS['mu_GeV']/MPL)/np.sqrt(.019416)
C2=ZQ*MUNIT**2; UCOEFF=MUNIT**4/V0
GAMMA=PARAMS['Gamma_MPl']/np.sqrt(V0)
POT=Potential(1.)
BASIS=[]
for i in range(4):
    b=np.zeros((4,4)); b[i,i]=1; BASIS.append(b)
for i in range(4):
    for j in range(i+1,4):
        b=np.zeros((4,4)); b[i,j]=b[j,i]=1/np.sqrt(2); BASIS.append(b)
BASIS=np.array(BASIS)
def mat(x): return np.einsum('i,ijk->jk',x,BASIS)
def vec(x): return np.einsum('ijk,jk->i',BASIS,x)
def U_and_gradient(D):
    Q=D+X0*np.eye(4); tr=np.trace(Q)-1; A=Q@Q-Q
    U=LAM/4*tr*tr+LP/4*np.sum(A*A)
    g=LAM/2*tr*np.eye(4)+LP/2*A@(2*Q-np.eye(4))
    return U,g
class Model:
    def __init__(self,alpha=3065566.3114906666*V0,xi=.25,eta=1.,gamma=GAMMA,K0=0.,coupled=True):
        self.alpha=alpha; self.xi=xi; self.eta=eta; self.gamma=gamma; self.K0=K0;self.coupled=coupled
        self.c2=C2*eta; self.u=UCOEFF*eta
    def terms(self,z):
        N,H,x,v,rr=z[:5]; D=mat(z[5:15]); Dv=mat(z[15:25]); K=self.K0*np.exp(-2*N)
        psi=H*H+K; margin=1-(self.alpha*psi)**2
        if margin<=0: raise ValueError('outside QT branch')
        hp=margin**-1.5; h=psi/np.sqrt(margin)
        U,g=U_and_gradient(D);V,dV,_=POT.values(x)
        S=self.c2*np.sum(D*D);Sd=2*self.c2*np.sum(D*Dv);qv2=self.c2*np.sum(Dv*Dv);qUg=self.u*np.sum(D*g)
        rhoc=.5*v*v+V;rhoQ=.5*qv2+self.u*U
        if self.coupled:
            den=2*(hp-self.xi*S)+12*self.xi**2*S
            num=2*(hp-self.xi*S)*K-v*v+(2*self.xi-1)*qv2-4*rr/3-2*self.xi*qUg-4*self.xi*H*Sd-12*self.xi**2*S*(2*H*H+K)
            Hd=num/den
        else: Hd=K-(v*v+4*rr/3)/(2*hp)
        R=6*(Hd+2*H*H+K)
        Dvv=-3*H*Dv-(UCOEFF/C2)*g-self.xi*R*D if self.eta else np.zeros((4,4))
        Sdd=2*self.c2*np.sum(Dv*Dv+D*Dvv)
        rho_nm=self.xi*(3*psi*S+3*H*Sd)
        p_nm=self.xi*(-(2*Hd+3*H*H+K)*S-Sdd-2*H*Sd)
        rho=rhoc+rr+(rhoQ+rho_nm if self.coupled else 0.)
        pressure=.5*v*v-V+rr/3+(.5*qv2-self.u*U+p_nm if self.coupled else 0.)
        return dict(Hd=Hd,R=R,margin=margin,hp=hp,h=h,K=K,Dvv=Dvv,dV=dV,V=V,
            S=S,Sd=Sd,Sdd=Sdd,rho=rho,pressure=pressure,rho_chi=rhoc,rho_Q_bare=rhoQ,
            rho_Q_nonminimal=rho_nm,pressure_Q_nonminimal=p_nm,F_eff=1-self.xi*S,
            constraint=3*h-rho,order=np.sum(D*D)/.75)
    def rhs(self,t,z):
        N,H,x,v,rr=z[:5]; term=self.terms(z)
        return np.r_[H,term['Hd'],v,-(3*H+self.gamma)*v-term['dV'],
            self.gamma*v*v-4*H*rr,z[15:25],vec(term['Dvv'])]
    def initial(self,seed=3203,amplitude=1e-3):
        rng=np.random.default_rng(seed);D=rng.normal(size=(4,4));D=(D+D.T)/2;D-=np.trace(D)/4*np.eye(4)
        D*=amplitude/np.linalg.norm(D)
        z=np.r_[0.,0.,PARAMS['chi_end'],PARAMS['chi_velocity_end_over_sqrtV0'],0.,vec(D),np.zeros(10)]
        def residual(H):
            zz=z.copy();zz[1]=H;return self.terms(zz)['constraint']
        hi=.99/np.sqrt(self.alpha) if self.alpha else 10.
        z[1]=brentq(residual,0,hi,xtol=1e-15)
        return z
