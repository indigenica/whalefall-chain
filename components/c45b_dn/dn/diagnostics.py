"""Grid diagnostics: derivatives, psi, Misner-Sharp-like mass m = r^3 h(psi)/(2G),
trapped regions, constraint residuals (independent finite differences)."""
import numpy as np


def derivs(A, h):
    Au = np.full_like(A, np.nan); Av = np.full_like(A, np.nan)
    Au[1:-1, :] = (A[2:, :] - A[:-2, :]) / (2 * h)
    Av[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / (2 * h)
    return Au, Av


def fields(grid, alpha, G=1.0):
    h = grid.h
    R, S, P = grid.r, grid.s, grid.p
    ru, rv = derivs(R, h)
    su, sv = derivs(S, h)
    pu, pv = derivs(P, h)
    e2 = np.exp(2 * S)
    with np.errstate(invalid="ignore", divide="ignore"):
        psi = (1 + 2 * ru * rv / e2) / R ** 2
        hpsi = psi / np.sqrt(1 - (alpha * psi) ** 2)
        mass = R ** 3 * hpsi / (2 * G)
    return dict(ru=ru, rv=rv, su=su, sv=sv, pu=pu, pv=pv, e2=e2, psi=psi, mass=mass)


def constraint_residuals(grid, alpha, G=1.0):
    """C_u = r_uu - 2 s_u r_u + r kappa phi_u^2 m^{3/2}/2 (likewise v), relative to the largest term."""
    h = grid.h; kappa = 8 * np.pi * G
    R, S, P = grid.r, grid.s, grid.p
    f = fields(grid, alpha, G)
    ruu = np.full_like(R, np.nan); rvv = np.full_like(R, np.nan)
    ruu[1:-1, :] = (R[2:, :] - 2 * R[1:-1, :] + R[:-2, :]) / h ** 2
    rvv[:, 1:-1] = (R[:, 2:] - 2 * R[:, 1:-1] + R[:, :-2]) / h ** 2
    with np.errstate(invalid="ignore"):
        m = 1 - (alpha * f["psi"]) ** 2
        m32 = m * np.sqrt(np.clip(m, 0, None))
        Tu = R * kappa * f["pu"] ** 2 * m32 / 2
        Tv = R * kappa * f["pv"] ** 2 * m32 / 2
        Cu = ruu - 2 * f["su"] * f["ru"] + Tu
        Cv = rvv - 2 * f["sv"] * f["rv"] + Tv
        scale_u = np.abs(ruu) + np.abs(2 * f["su"] * f["ru"]) + np.abs(Tu)
        scale_v = np.abs(rvv) + np.abs(2 * f["sv"] * f["rv"]) + np.abs(Tv)
    return Cu, Cv, scale_u, scale_v
