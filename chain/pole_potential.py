"""The order-4 pole plateau of U22 (docs/preregistros/ESPEC_U22.md), V = V0 [1 - C/sqrt(chi^2 + C^2)], with the exact V0 (1 - C/chi)
tail; C = 0.3516 is the shape fit of the 10-level sum's tail.  Default plateau of the chain since U26 (ESPEC_U26.md)."""
import math
import numpy as np
from scipy.optimize import brentq

C_FIT = 0.3516


class PolePotential:
    def __init__(self, V0=1., C=C_FIT): self.V0, self.C = V0, C
    def values(self, x):
        x = np.asarray(x, float); q = x * x + self.C ** 2
        return (self.V0 * (1 - self.C / np.sqrt(q)), self.V0 * self.C * x / q ** 1.5, self.V0 * self.C * (q - 3 * x * x) / q ** 2.5)
    def V(self, x): return self.values(x)[0]
    def dV(self, x): return self.values(x)[1]
    def pzeta(self, x):
        v, dv, _ = self.values(x); return v ** 3 / (12 * np.pi ** 2 * dv ** 2)
    def eternal(self):
        """chi where P_zeta = 1 (same root as u22_pole.chi_eternal_pole)"""
        return brentq(lambda x: 3 * math.log(float(self.V(x))) - math.log(12 * math.pi ** 2) - 2 * math.log(float(self.dV(x))), 5.0, 1e7, xtol=1e-6)
