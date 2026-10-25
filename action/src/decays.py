"""Tree-level vacuum Q widths, resolved inside degenerate mass subspaces.
No inherited scalaron mass or cosmological shortcut is used.
"""
from common import *
from model import *
import sympy as sp
from c33_q_decays_sympy import build

def run():
    data=build(); fields=data['fields']; zero=dict.fromkeys(fields,0)
    cf=np.array([1/np.sqrt(ZQ)]+[1/np.sqrt(ZQ if j==k else 2*ZQ) for j,k in data['bkeys']])
    m2=np.array([[float(sp.diff(data['U'],a,b).subs(zero))*cf[i]*cf[j]
        for j,b in enumerate(fields)] for i,a in enumerate(fields)])
    w,rot=np.linalg.eigh(m2)
    # L3=phi_a G[a,i,j] d pi_i d pi_j. Width is a Gram matrix, not
    # seven separate widths in an arbitrary basis of a sixfold degeneracy.
    G=np.array([[[float(sp.diff(data['G'][i,j],a).subs(zero))*cf[k]/4
        for j in range(3)] for i in range(3)] for k,a in enumerate(fields)])
    Gmode=np.einsum('an,aij->nij',rot,G)
    gram=np.einsum('aij,bij->ab',Gmode,Gmode)/2
    widthlight=w[0]**1.5/(16*np.pi)*gram[:6,:6]
    rates,decrot=np.linalg.eigh(widthlight)
    rot[:,:6]=rot[:,:6]@decrot
    cubic=np.array([[[float(sp.diff(data['U'],a,b,c).subs(zero))*cf[i]*cf[j]*cf[k]
        for k,c in enumerate(fields)] for j,b in enumerate(fields)] for i,a in enumerate(fields)])
    cubmode=np.einsum('abc,ai,bj,ck->ijk',cubic,rot,rot,rot)
    masses=np.sqrt(w);channels=[];heavywidth=0.
    for i in range(6):
        for j in range(i,6):
            pcm=np.sqrt((masses[6]**2-(masses[i]+masses[j])**2)*(masses[6]**2-(masses[i]-masses[j])**2))/(2*masses[6])
            rate=cubmode[6,i,j]**2*pcm/(8*np.pi*masses[6]**2)*(0.5 if i==j else 1.)
            if rate>1e-16:
                channels.append(dict(daughters=[i,j],width_code=float(rate)));heavywidth+=rate
    widths=np.r_[rates,heavywidth];E=MUNIT*MPL
    result=dict(masses_code=masses,widths_code=widths,masses_GeV=masses*E,widths_GeV=widths*E,
        Gamma_chi_GeV=PARAMS['Gamma_GeV'],min_Q_width_over_Gamma_chi=float(widths.min()*E/PARAMS['Gamma_GeV']),
        light_Goldstone_width_eigenvalues=rates,heavy_Goldstone_width_code=float(masses[6]**3/(16*np.pi)*gram[6,6]),
        heavy_channels=channels,kinetic_normalization=np.eye(3)*.5,
        all_modes_have_channels=bool(np.all(widths>0)),
        scope='Vacuum tree widths. Does not give defect emission, thermal rates, or an inhomogeneous reheating simulation.',
        canonical_mass_eigenvectors=rot,canonical_factors=cf)
    save('decays.json',result);print(json.dumps({k:v for k,v in result.items() if k in ('widths_GeV','masses_GeV','min_Q_width_over_Gamma_chi','all_modes_have_channels')},default=lambda x:x.tolist(),indent=2))
if __name__=='__main__':run()
