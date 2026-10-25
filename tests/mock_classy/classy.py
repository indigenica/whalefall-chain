"""MOCK of classy, ONLY for testing the plumbing of lya/u4a_class.py where CLASS cannot be installed.
It answers Class.set/compute/pk/get_background with lya/boltz.py on a coarse grid.  Never used for results."""
import math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lya"))
import boltz as B

__version__ = "MOCK"
KB_TCMB = 8.617333262e-5 * 2.7255
_KG = np.geomspace(0.05, 1000.0, 16)          # h/Mpc


class Class:
    def __init__(self): self.p = {}
    def set(self, p): self.p.update(p)
    def struct_cleanup(self): pass
    def empty(self): self.p = {}

    def _species(self):
        p = self.p
        if not p.get("N_ncdm"): return None
        om = p["omega_ncdm"]; m = p["m_ncdm"]; T = p["T_ncdm"]
        if p.get("use_ncdm_psd_files"):
            q, f0 = np.loadtxt(p["ncdm_psd_filenames"], unpack=True)
            sel = np.linspace(0, len(q) - 1, 24).astype(int); q, f0 = q[sel], np.maximum(f0[sel], 1e-300)
            dl = np.gradient(np.log(q)); W = q ** 3 * f0 * dl; dlnf = np.gradient(np.log(f0)) / dl
            return B.Species(q * T * KB_TCMB / m, W, dlnf, om / B.h ** 2, lmax=8)
        Tnu = (4 / 11) ** (1 / 3)
        return B.thermal_fd(m, om, nq=24, lmax=8)

    def compute(self):
        d, _ = B.solve(_KG * B.h, self._species(), phase_max=0.6)
        self._lp = np.log(_KG ** (B.NS - 4) * d ** 2)

    def pk(self, k_Mpc, z):
        kh = k_Mpc / B.h
        return float(np.exp(np.interp(math.log(kh), np.log(_KG), self._lp)))

    def get_background(self):
        om = self.p.get("omega_ncdm", 0.0)
        return {"(.)rho_ncdm[0]": np.array([om / B.h ** 2]), "(.)rho_crit": np.array([1.0])}
