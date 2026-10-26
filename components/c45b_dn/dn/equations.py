"""Double-null NPQT field equations for a real scalar with potential.

Metric  ds^2 = -2 e^{2 sigma} du dv + r^2 dOmega^2,  e2 = e^{2 sigma},
psi = (1 + 2 r_u r_v / e2) / r^2,  margin m = 1 - (alpha psi)^2,  kappa = 8 pi G.
Same reduction and conventions as the team's audited evaluator
(claude_dual/sidequests/c58_spherical_cross_rescaled/rhs.py, eqs. (33)-(34) of
Bueno-Cano-Hennigar-Murcia); re-implemented here independently:

  r_uv / r = kappa T_uv m^{3/2}/2 - e2 psi (3m - 2)/2
  sigma_uv = -r_uv/r + e2/2 [3 alpha^2 psi^3 - 6 alpha^2 psi m J - kappa T_ang m^{3/2}]
  J = b^2 - kappa^2 T_uu T_vv m/(4 e2^2),  b = kappa T_uv sqrt(m)/(2 e2) - 3 psi/2
  constraints: r_uu - 2 sigma_u r_u = -r kappa T_uu m^{3/2}/2,  r_vv - 2 sigma_v r_v = -r kappa T_vv m^{3/2}/2
Scalar (canonical, potential V):
  T_uu = phi_u^2, T_vv = phi_v^2, T_uv = e2 V, T_ang = e^{-2 sigma} phi_u phi_v - V
  phi_uv = -(r_u phi_v + r_v phi_u)/r - e2 V'(phi)/2
"""
import numpy as np


def psi_of_y(y, alpha):
    return y / np.sqrt(1 + (alpha * y) ** 2)


def rhs(r, sigma, phi, ru, rv, su, sv, pu, pv, V, dV, alpha, G=1.0, psi_override=None):
    """Return (r_uv, sigma_uv, phi_uv, psi, margin).  psi_override (array or None)
    replaces the reconstructed psi where finite (regular centre series)."""
    kappa = 8 * np.pi * G
    e2 = np.exp(2 * sigma)
    psi = (1 + 2 * ru * rv / e2) / (r * r)
    if psi_override is not None:
        psi = np.where(np.isfinite(psi_override), psi_override, psi)
    q = alpha * psi
    m = (1 - q) * (1 + q)
    m32 = m * np.sqrt(np.clip(m, 0, None))
    Vp = V(phi)
    Tuu, Tvv, Tuv = pu * pu, pv * pv, e2 * Vp
    Tang = pu * pv / e2 - Vp
    uv_over_r = kappa * Tuv * m32 / 2 - e2 * psi * (3 * m - 2) / 2
    b = kappa * Tuv * np.sqrt(np.clip(m, 0, None)) / (2 * e2) - 1.5 * psi
    J = b * b - kappa ** 2 * Tuu * Tvv * m / (4 * e2 * e2)
    suv = -uv_over_r + e2 / 2 * (3 * q * q * psi - 6 * alpha * q * m * J - kappa * Tang * m32)
    puv = -(ru * pv + rv * pu) / r - e2 * dV(phi) / 2
    return r * uv_over_r, suv, puv, psi, m


def local_density(sigma, pu, pv, phi, V):
    """Energy density seen by n = (d_u + d_v)/sqrt(2 e2): (phi_u^2+phi_v^2)/(2 e2) + V."""
    e2 = np.exp(2 * sigma)
    return (pu * pu + pv * pv) / (2 * e2) + V(phi)
