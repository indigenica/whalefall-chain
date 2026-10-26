"""Growing-mode LTB-QT: an overdensity in an Einstein-de Sitter dust parent (uniform bang time),
as in the long-wavelength (gradient-expansion) description of primordial-black-hole formation.

Comoving areal coordinate x (= r-hat): dust mass M(x) = mu x^3 (uniform density at leading order),
and the shell energy fixed by the curvature profile:
    -2E(x) = K(x) x^2 = 1 - (1 + r zeta'(r))^2,   x = r e^{zeta(r)}     (Harada et al. 2015;
    Musco 2019; Escriva, Germani & Sheth 2020: C = (2/3) K r^2 = -(2/3) r zeta' (2 + r zeta')).
Each shell obeys the exact QT shell equation  Rdot^2 = R^2 psi(2M/R^3) - K_j  (ltb/check_ms_mass.py).

Synchronization: at a common early time t_i, when the background has alpha*y_b = 1e-3 (QT negligible
there and later until the bounce), every shell sits on its GR growing-mode cycloid
    R = (M/K)(1 - cos eta),  t = (M/K^{3/2})(eta - sin eta).
From t_i on the motion is integrated with the exact QT equation: expansion to R_ta (turnaround,
R_ta^2 psi(2M/R_ta^3) = K), then the periodic QT orbit of ltb/shells.Shell (time-symmetric about t_ta).
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from ltb.shells import Shell, psi_of_y


def turnaround_radius(M, K, alpha):
    """Large-R root of R^2 psi(2M/R^3) = K (beyond the peak of R^2 psi)."""
    f = lambda R: R * R * psi_of_y(2 * M / R ** 3, alpha) - K
    Rgr = 2 * M / K
    # the peak of R^2 psi lies at R ~ (2 alpha M)^{1/3}; scan from well below it (GR: no peak)
    lo = 1e-3 * (2 * alpha * M) ** (1 / 3) if alpha > 0 else 1e-6 * Rgr
    Rs = np.geomspace(lo, 4 * Rgr, 2000)
    v = f(Rs)
    if v.max() <= 0:
        raise ValueError("no turnaround: shell bound inside the QT core (K above the peak of R^2 psi)")
    k = int(np.argmax(v))
    hi_idx = k + int(np.argmax(v[k:] < 0))
    return brentq(f, Rs[hi_idx - 1], Rs[hi_idx], xtol=1e-14 * Rgr, rtol=1e-14)


def t_turnaround(M, K, alpha, t_i):
    """Proper time of turnaround for a growing-mode shell (bang at t = 0)."""
    # GR cycloid position at t_i
    tau = t_i * K ** 1.5 / M
    eta = brentq(lambda e: e - math.sin(e) - tau, 0.0, math.pi)
    R_i = (M / K) * (1 - math.cos(eta))
    R_ta = turnaround_radius(M, K, alpha)
    U = lambda R: R * R * psi_of_y(2 * M / R ** 3, alpha) - K
    # substitution R = R_ta - (R_ta - R_i) s^2 removes the simple root at R_ta
    D = R_ta - R_i
    g = lambda s: 2 * D * s / math.sqrt(max(U(R_ta - D * s * s), 1e-300)) if s > 0 else 2 * D / math.sqrt(D * max(-dU(R_ta, M, alpha), 1e-300))
    val, err = quad(g, 0.0, 1.0, limit=400, epsabs=0, epsrel=1e-11)
    return t_i + val, R_ta, R_i


def dU(R, M, alpha, h=1e-7):
    f = lambda r: r * r * psi_of_y(2 * M / r ** 3, alpha)
    return (f(R * (1 + h)) - f(R * (1 - h))) / (2 * R * h)


class GrowingShell:
    def __init__(self, M, K, alpha, t_i):
        self.M, self.K = M, K
        self.t_ta, self.R_ta, self.R_i = t_turnaround(M, K, alpha, t_i)
        self.shell = Shell(M, self.R_ta, alpha)
        self.tau_B = self.shell.tau_B
        self.t_B = self.t_ta + self.tau_B

    def R_of_t(self, t):
        return self.shell.R_of_t(np.abs(np.asarray(t, float) - self.t_ta))
