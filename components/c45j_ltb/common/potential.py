"""Ballena-4 benchmark potential in C59 units: G = 1 (M_Pl = 1), m_chi = 1.

V(chi) = (V0/Z) sum_{n<L} q^n tanh^2(chi/(f s^n)),  q=1/2, s=2, f = Mbar/2, L=10,
V0 = 1.67e-10 Mbar^4, Mbar^2 = 1/(8 pi).  Lengths/times in 1/m_chi, so the
potential is divided by m_chi^2 (V''(0) = 1 by construction).
"""
import numpy as np

MBAR = 1 / np.sqrt(8 * np.pi)
Q, S, F, L = 0.5, 2.0, 0.5 * MBAR, 10
V0_PLANCK = 1.67e-10 * MBAR ** 4
_Z = sum(Q ** n for n in range(L))
_M2 = 2 * V0_PLANCK / (_Z * F * F) * sum((Q / S ** 2) ** n for n in range(L))
M_CHI = np.sqrt(_M2)                                  # 5.514e-6 (Planck units)
V0 = V0_PLANCK / _M2                                  # plateau in m=1 units: 8.695e-3
H_INF = np.sqrt(8 * np.pi * V0_PLANCK / 3) / M_CHI    # H_inf / m_chi = 0.270
_W = np.array([Q ** n for n in range(L)]) / _Z
_SC = np.array([F * S ** n for n in range(L)])


def V(phi):
    x = np.asarray(phi, float)[..., None] / _SC
    return V0 * np.sum(_W * np.tanh(x) ** 2, axis=-1)


def dV(phi):
    x = np.asarray(phi, float)[..., None] / _SC
    t = np.tanh(x)
    return V0 * np.sum(_W * 2 * t * (1 - t * t) / _SC, axis=-1)


def V_quad(phi):
    return 0.5 * np.asarray(phi, float) ** 2


def dV_quad(phi):
    return np.asarray(phi, float)
