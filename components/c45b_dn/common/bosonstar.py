"""C59a-1 - static (complex-proxy) chi-stars in spherical NPQT gravity.

Units G = 1, m_chi = 1.  Phi = phi(r) e^{-i w t}/sqrt(2), potential U = V(phi).
Metric ds^2 = -N^2 f dt^2 + dr^2/f + r^2 dOmega^2,  f = 1 - r^2 psi,
h(psi) = psi/sqrt(1 - alpha^2 psi^2) = y = 2m/r^3  ->  psi = y/sqrt(1+alpha^2 y^2).

  rho   = w^2 phi^2/(2 N^2 f) + f phi'^2/2 + V
  p_r   = w^2 phi^2/(2 N^2 f) + f phi'^2/2 - V
  m'    = 4 pi r^2 rho
  N'/N  = 4 pi r (rho + p_r)/(f h'(psi)),   1/h' = (1 + alpha^2 y^2)^(-3/2)
  (N r^2 f phi')'/(N r^2) + w^2 phi/(N^2 f) = V'(phi)

Same reduction as c58c/regular_center (with a potential added to rho only;
rho + p_r is potential-free).  Shooting in w for the nodeless ground state;
N(0) = 1 during integration, rescaled so that N(inf) = 1 afterwards.
"""
import numpy as np
from scipy.integrate import solve_ivp


def _geom(r, m, alpha):
    y = 2 * m / r ** 3
    inv = 1 / np.sqrt(1 + (alpha * y) ** 2)
    psi = y * inv
    return y, psi, inv ** 3          # y, psi, 1/h'


def rhs_factory(w, alpha, V, dV):
    def rhs(r, z):
        phi, dphi, m, lnN = z
        N = np.exp(lnN)
        y, psi, invhp = _geom(r, m, alpha)
        f = 1 - r * r * psi
        kin = w * w * phi * phi / (2 * N * N * f)
        grad = 0.5 * f * dphi * dphi
        rho = kin + grad + V(phi)
        rpp = 2 * (kin + grad)
        dm = 4 * np.pi * r * r * rho
        dlnN = 4 * np.pi * r * rpp / f * invhp
        dy = 2 * dm / r ** 3 - 3 * y / r
        dpsi = dy * invhp
        df = -2 * r * psi - r * r * dpsi
        ddphi = (-(2 / r + dlnN + df / f) * dphi - w * w * phi / (N * N * f * f) + dV(phi) / f)
        return [dphi, ddphi, dm, dlnN]
    return rhs


def shoot(phi0, w, alpha, V, dV, r_max, r0=1e-4):
    """Integrate outward; return (status, sol). status: +1 diverges up, -1 crosses zero."""
    # regular series at the centre (f=1, N=1): 3 phi2 = V'(phi0) - w^2 phi0
    phi2 = (dV(phi0) - w * w * phi0) / 6.0
    rho0 = w * w * phi0 * phi0 / 2 + V(phi0)
    z0 = [phi0 + phi2 * r0 ** 2, 2 * phi2 * r0, 4 * np.pi / 3 * rho0 * r0 ** 3,
          4 * np.pi * (w * w * phi0 * phi0) * r0 ** 2 / 2]
    rhs = rhs_factory(w, alpha, V, dV)

    def cross(r, z):
        return z[0]
    cross.terminal = True; cross.direction = -1

    def blow(r, z):
        return z[0] - 1.5 * phi0
    blow.terminal = True; blow.direction = 1

    def horizon(r, z):
        y, psi, _ = _geom(r, z[2], alpha)
        return 1 - r * r * psi - 1e-6
    horizon.terminal = True; horizon.direction = -1
    sol = solve_ivp(rhs, (r0, r_max), z0, method="DOP853", rtol=1e-11, atol=1e-14,
                    events=[cross, blow, horizon], dense_output=True)
    if sol.t_events[2].size:
        return 0, sol
    if sol.t_events[0].size:
        return -1, sol
    if sol.t_events[1].size or sol.y[1, -1] > 0:
        return +1, sol
    return (-1 if sol.y[0, -1] < 0 else +1), sol


def ground_state(phi0, alpha=0.0, V=None, dV=None, r_max=60.0, w_lo=0.3, w_hi=1.6, tol=1e-12):
    """Bisection in w (in units of the central lapse) for the nodeless state."""
    from common.potential import V_quad, dV_quad
    V = V or V_quad; dV = dV or dV_quad
    s_lo, _ = shoot(phi0, w_lo, alpha, V, dV, r_max)
    s_hi, _ = shoot(phi0, w_hi, alpha, V, dV, r_max)
    while s_hi != -1 and w_hi < 50:
        w_hi *= 1.5
        s_hi, _ = shoot(phi0, w_hi, alpha, V, dV, r_max)
    while s_lo == -1 and w_lo > 1e-3:
        w_lo *= 0.5
        s_lo, _ = shoot(phi0, w_lo, alpha, V, dV, r_max)
    s_lo = +1 if s_lo == 0 else s_lo
    if not (s_lo == +1 and s_hi == -1):
        raise RuntimeError(f"bracket failed phi0={phi0} ({s_lo},{s_hi})")
    for _ in range(200):
        w = 0.5 * (w_lo + w_hi)
        s, sol = shoot(phi0, w, alpha, V, dV, r_max)
        # a trial that reaches f = 0 has too much field energy: same side as divergence
        if s in (0, +1):
            w_lo = w
        else:
            w_hi = w
        if w_hi - w_lo < tol:
            break
    return 0.5 * (w_lo + w_hi), sol


def measure(phi0, alpha=0.0, V=None, dV=None, r_max=60.0):
    """Mass and radii of the ground state.  The tail (after the field has
    decayed to the bisection limit) is cut where |phi| is minimal."""
    w, sol = ground_state(phi0, alpha, V, dV, r_max)
    r = np.linspace(sol.t[0], sol.t[-1], 20001)
    z = sol.sol(r)
    phi, m, lnN = z[0], z[2], z[3]
    # cut at the minimum of |phi| (beyond it the shooting error grows)
    k = int(np.argmin(np.abs(phi) + (r < 1.0) * 1e9))
    r, phi, m, lnN = r[:k], phi[:k], m[:k], lnN[:k]
    M = m[-1]
    Ninf = np.exp(lnN[-1]) * np.sqrt(max(1 - 2 * M / r[-1], 1e-300))   # Schwarzschild-like tail
    radii = {}
    for frac in (0.5, 0.9, 0.99):
        radii[f"R{int(frac * 100)}"] = float(np.interp(frac * M, m, r))
    rho_bar = {k_: 3 * M * frac / (4 * np.pi * rr ** 3)
               for (k_, rr), frac in zip(radii.items(), (0.5, 0.9, 0.99))}
    y, psi, _ = _geom(r[1:], m[1:], alpha)
    fmin = float(np.min(1 - r[1:] ** 2 * psi))
    return dict(phi0=phi0, alpha=alpha, w=w / Ninf, M=float(M), r_cut=float(r[-1]),
                **radii, rho_bar=rho_bar, f_min=fmin,
                compactness_R99=2 * M / radii["R99"], alpha_y_center=float(alpha * 2 * m[1] / r[1] ** 3))
