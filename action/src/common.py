from pathlib import Path
import os, sys, json
ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'runtime'))
sys.path.insert(1,str(ROOT/'vendor'))
os.environ['MPLCONFIGDIR']=str(ROOT/'runtime'/'mpl_config')
import numpy as np
from scipy.optimize import brentq
# U26 (docs/preregistros/ESPEC_U26.md): WHALEFALL_PLATEAU = pole (default) | sum. With sum this folder is the unchanged copy of
# cierre_accion (results/, logs/); with the pole, Potential is the order-4 pole of chain/pole_potential.py, the parameters come from
# vendor/thermal_pole.json and everything is written to results_pole/ and logs_pole/.
PLATEAU=os.environ.get('WHALEFALL_PLATEAU','pole')
if PLATEAU not in ('pole','sum'): raise ValueError(PLATEAU)
RES=ROOT/('results' if PLATEAU=='sum' else 'results_pole'); LOGS=ROOT/('logs' if PLATEAU=='sum' else 'logs_pole')
THERMAL=ROOT/'vendor'/('thermal.json' if PLATEAU=='sum' else 'thermal_pole.json')
def save(name,data):
    def convert(x):
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,np.generic): return x.item()
        raise TypeError(type(x).__name__)
    (RES/name).write_text(json.dumps(data,indent=2,default=convert,allow_nan=False),encoding='utf-8')
class SumPotential:
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
def _pole_class():
    import importlib.util
    spec=importlib.util.spec_from_file_location('pole_potential',ROOT.parent/'chain'/'pole_potential.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod.PolePotential
PolePotential=_pole_class()
Potential=PolePotential if PLATEAU=='pole' else SumPotential
