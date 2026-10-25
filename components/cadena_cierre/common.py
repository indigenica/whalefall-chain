from pathlib import Path
import os, sys, json
ROOT=Path(__file__).resolve().parents[2]  # adapted: repository root (was campaign root)
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'components'/'cadena_cierre'))  # adapted
sys.path.insert(1,str(ROOT/'components'))  # adapted
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'logs'/'mpl_config'))  # adapted
import numpy as np
from scipy.optimize import brentq
def save(name,data):
    def convert(x):
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,np.generic): return x.item()
        raise TypeError(type(x).__name__)
    (ROOT/'results'/'cadena_rerun'/name).write_text(json.dumps(data,indent=2,default=convert,allow_nan=False),encoding='utf-8')
class Potential:
    def __init__(self,V0=1.,L=10):
        self.V0=V0; self.L=L
        self.w=.5**np.arange(L); self.sc=.5*2.**np.arange(L); self.Z=self.w.sum()
    def values(self,x):
        x=np.asarray(x); u=x[...,None]/self.sc
        t=np.tanh(u); e=np.exp(-2*np.abs(u)); sech=4*e/(1+e)**2
        v=self.V0*np.sum(self.w*t*t,axis=-1)/self.Z
        dv=self.V0*np.sum(self.w*2*t*sech/self.sc,axis=-1)/self.Z
        ddv=self.V0*np.sum(self.w*2*sech*(sech-2*t*t)/self.sc**2,axis=-1)/self.Z
        return v,dv,ddv
    def V(self,x): return self.values(x)[0]
    def dV(self,x): return self.values(x)[1]
    def pzeta(self,x):
        v,dv,_=self.values(x); return v**3/(12*np.pi**2*dv**2)
    def eternal(self): return brentq(lambda x: np.log(self.pzeta(x)),10,2000)
