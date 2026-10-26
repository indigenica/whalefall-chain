"""LTB-QT: pressureless shells in spherical NPQT gravity (derivation: ltb/check_ms_mass.py).

Each comoving shell (mass M, starting at rest with areal radius R0) obeys exactly
    Rdot^2 = U(R) = R^2 psi(2 G M / R^3) - R0^2 psi(2 G M / R0^3),   psi(y) = y/sqrt(1+alpha^2 y^2)
until shell crossing.  In QT, R^2 psi -> R^2/alpha as R -> 0, so U has a second root
R_b > 0 (bounce); the motion is periodic: R0 -> R_b in tau_B, back to R0 at 2 tau_B.
GR (alpha = 0): R_b = 0 and tau_B = (pi/2) sqrt(R0^3/(2 G M)) (Oppenheimer-Snyder).

tau(R) on the collapse branch is computed with the substitution
    R(theta) = R_b + (R0 - R_b)(1 - cos theta)/2,   theta in [0, pi]   (theta = 0: bounce),
which removes both inverse-square-root endpoint singularities (simple roots of U).
"""
import numpy as np
from scipy.optimize import brentq


def psi_of_y(y, alpha):
    return y / np.sqrt(1 + (alpha * y) ** 2)


class Shell:
    def __init__(self, M, R0, alpha, G=1.0, n_theta=4001):
        self.M, self.R0, self.alpha, self.G = float(M), float(R0), float(alpha), float(G)
        self.K = R0 ** 2 * psi_of_y(2 * G * M / R0 ** 3, alpha)          # = -2E
        self.E = -0.5 * self.K
        U = lambda R: self.Rpsi(R) - self.K
        # collapse requires U > 0 just below R0 (R0 beyond the peak of R^2 psi)
        self.collapses = U(R0 * (1 - 1e-7)) > 0
        if not self.collapses:
            self.R_b = None; self.tau_B = np.inf; return
        if alpha == 0:
            self.R_b = 0.0
        else:
            # R^2 psi is increasing from 0 up to its peak R_pk, then decreasing; R_b in (0, R_pk)
            Rs = np.geomspace(R0 * 1e-8, R0, 4000)
            vals = U(Rs)
            k = int(np.argmax(vals))
            lo = np.nonzero(vals[:k] < 0)[0]
            self.R_b = brentq(U, Rs[lo[-1]], Rs[lo[-1] + 1], xtol=1e-15, rtol=1e-14) if len(lo) else Rs[0]
        th = np.linspace(0, np.pi, n_theta)
        Rb = self.R_b
        R = Rb + (R0 - Rb) * (1 - np.cos(th)) / 2
        dRdth = (R0 - Rb) * np.sin(th) / 2
        with np.errstate(invalid="ignore", divide="ignore"):
            integrand = dRdth / np.sqrt(np.maximum(U(R), 0))
        # endpoints: limit of dR/dth / sqrt(U) (finite); use one-sided extrapolation
        integrand[0] = 2 * integrand[1] - integrand[2]
        integrand[-1] = 2 * integrand[-2] - integrand[-3]
        # cumulative time from the bounce (theta = 0) outwards
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(th))])
        self.tau_B = float(cum[-1])
        self._t_from_bounce = cum            # time since bounce at radius R (expanding branch)
        self._R = R

    def Rpsi(self, R):
        R = np.asarray(R, float)
        return R ** 2 * psi_of_y(2 * self.G * self.M / R ** 3, self.alpha)

    def R_of_t(self, t):
        """Areal radius at proper time t >= 0 (t = 0 at rest at R0); periodic with 2 tau_B."""
        t = np.asarray(t, float)
        if not self.collapses:
            return np.full_like(t, self.R0)
        if self.alpha == 0:
            tt = np.minimum(t, self.tau_B)       # GR: ends at the singularity
            s = self.tau_B - tt
        else:
            ph = np.mod(t, 2 * self.tau_B)
            s = np.abs(ph - self.tau_B)          # time from the (nearest) bounce
        return np.interp(s, self._t_from_bounce, self._R)

    def Rdot_sign(self, t):
        ph = np.mod(np.asarray(t, float), 2 * self.tau_B)
        return np.where(ph < self.tau_B, -1.0, 1.0)

    def f_min(self):
        """min over the orbit of f = 1 - R^2 psi = 1 + 2E - Rdot^2; trapped if < 0."""
        if not self.collapses:
            return 1 - self.K
        R = self._R
        Ud = self.Rpsi(R) - self.K
        return float(np.min(1 - self.K - Ud))
