"""Finite-UV-window slow-roll branching diffusion, time = e-folds.
Binary Poisson branching rate 3 is an explicit effective stochastic toy.
"""
from common import *
from scipy.linalg import eigh_tridiagonal, solve_banded
def coeff(p,x):
    v,dv,ddv=p.values(x)
    b=-dv/v; D=v/(24*np.pi**2)
    bp=-ddv/v+(dv/v)**2; Dp=dv/(24*np.pi**2)
    U=3-b*b/(4*D)-bp/2+b*Dp/(2*D)
    return b,D,U
def eigen(p,center,width,n):
    h=width/n; x=np.linspace(center-width/2,center+width/2,n+1)[1:-1]
    b,D,U=coeff(p,x)
    diag=U-2*D/h**2; off=np.sqrt(D[:-1]*D[1:])/h**2
    val=eigh_tridiagonal(diag,off,eigvals_only=True,select='i',select_range=(len(x)-1,len(x)-1))[0]
    return float(val)
def survival(p,center,width,n=1024):
    h=width/n; x=np.linspace(center-width/2,center+width/2,n+1)[1:-1]
    b,D,_=coeff(p,x); left=D/h**2-b/(2*h); right=D/h**2+b/(2*h)
    if min(left.min(),right.min())<=0: raise ValueError('advection underresolved')
    diag=-2*D/h**2
    ev=eigen(p,center,width,n)
    if ev<=0: return dict(lambda_volume=ev,prob_center=0.,residual=0.,n=n)
    s=np.ones(n-1)
    for it in range(80):
        F=diag*s+left*np.r_[0,s[:-1]]+right*np.r_[s[1:],0]+3*s*(1-s)
        band=np.zeros((3,n-1)); band[0,1:]=right[:-1]; band[2,:-1]=left[1:]; band[1]=diag+3*(1-2*s)
        delta=solve_banded((1,1),band,-F)
        fac=1.
        while np.any(s+fac*delta<0) or np.any(s+fac*delta>1+1e-10):
            fac/=2
            if fac<1e-12: raise RuntimeError('survival Newton line search failed')
        s+=fac*delta
        if np.max(abs(F))<1e-9: break
    F=diag*s+left*np.r_[0,s[:-1]]+right*np.r_[s[1:],0]+3*s*(1-s)
    return dict(lambda_volume=ev,prob_center=float(s[len(s)//2]),residual=float(np.max(abs(F))),n=n,
                iterations=it+1,lower_boundary=center-width/2,upper_boundary=center+width/2)
def run():
    thermal=json.loads((ROOT/'results'/'thermal.json').read_text())['continuous']
    p=Potential(thermal['V0']); et=p.eternal()
    critical=brentq(lambda x: np.log(p.pzeta(x)*6),10,et)
    # independent analytic constant-drift eigenvalue and central FD error.
    D=2e-12;b=-2e-6;W=1e-4;n=256;h=W/n
    analytic=3-b*b/(4*D)-D*np.pi**2/W**2
    matrix=eigh_tridiagonal(np.full(n-1,3-b*b/(4*D)-2*D/h**2),np.full(n-2,D/h**2),
        eigvals_only=True,select='i',select_range=(n-2,n-2))[0]
    rows=[]
    for pz in (.1,1/6,.3,1.,3.):
        x=brentq(lambda x: np.log(p.pzeta(x)/pz),10,1000)
        for W in (1e-5,1e-4,1e-3):
            vals=[eigen(p,x,W,n) for n in (128,256,512)]
            row=dict(P_zeta=pz,center=x,width=W,lambda_volume=vals,
                     positive=bool(vals[-1]>0),grid_difference=abs(vals[-1]-vals[-2]))
            if W==1e-4 and pz in (.1,.3,1.):
                row['binary_branching_survival']=[survival(p,x,W,n) for n in (512,1024)]
            rows.append(row)
    save('stochastic.json',dict(V0=p.V0,chi_Pzeta_1=et,chi_local_volume_critical_Pzeta_1_over_6=critical,
        constant_drift_analytic=analytic,constant_drift_numeric=matrix,constant_drift_error=abs(matrix-analytic),
        windows=rows,pass_numerics=bool(abs(matrix-analytic)<1e-5 and max(r['grid_difference'] for r in rows)<1e-3),
        assumptions=['GR slow-roll stochastic description after QT throat; Planck reduced units',
          'Absorbing lower and upper boundaries discard descendants leaving each small window',
          'Survival is for a binary Poisson branching diffusion rate 3, not exact quantum-gravity genealogy',
          'Canonical kinetic term and potential validity assumed at the window field value',
          'A positive eigenvalue demonstrates growth of expected inflating volume for this effective process',
          'Does not prove that the model produces a global closed seed or an eternal controlled cosmology']))
    print('STOCHASTIC',et,critical,'constant check',abs(matrix-analytic),flush=True)
if __name__=='__main__': run()
