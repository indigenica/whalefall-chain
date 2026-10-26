"""Characteristic (double-null) evolution with a regular centre, NumPy,
vectorised over anti-diagonals.

Domain: u_i = i h, v_j = j h, 0 <= i <= j <= N.  Centre r = 0 on u = v.
Gauge: r(0, v) = v/2 on the initial outgoing cone; u = v on the centre.
Regularity at r = 0: e2 = -2 r_u r_v = 2 r_v^2 (so that psi stays finite);
sigma at the centre node from that relation with r_v from the ray v = v_k
(second order); phi even in r (quadratic extrapolation along u + v = const).
Cells within n_reg cells of the centre use the regular series
psi = psi(y), y = 8 pi G rho_local / 3 (exact at r = 0, O(r^2) error), instead
of the cancellation-prone (1 + 2 r_u r_v/e2)/r^2.

Diamond scheme (second order, implicit, fixed-point iterated):
  z_N = z_NW + z_SE - z_SW + h^2 F(zbar, z_u, z_v),  zbar = mean of the 4 corners.
No floors, no clipping: a node with margin <= 0, non-finite values or
r < 0 off the centre marks the run as failed at that anti-diagonal.
"""
import numpy as np
from dn.equations import rhs, psi_of_y, local_density


class Grid:
    def __init__(self, N, h):
        self.N, self.h = int(N), float(h)
        shape = (self.N + 1, self.N + 1)
        self.r = np.full(shape, np.nan)
        self.s = np.full(shape, np.nan)
        self.p = np.full(shape, np.nan)
        self.M = np.full(shape, np.nan)      # generalized Misner-Sharp mass r^3 h(psi)/(2G)
        self.ok_until = None


def initial_ray(grid, F, dF, V, dV, alpha, G=1.0, n_sub=4):
    """Data on u = 0: r = v/2, phi = F(r) (F, dF as functions of r).
    Integrate sigma and r_u along v from the centre (RK4 with n_sub substeps)."""
    N, h = grid.N, grid.h
    v = h * np.arange(N + 1)
    grid.r[0] = v / 2
    grid.p[0] = F(v / 2)
    kappa = 8 * np.pi * G

    def deriv(vv, sig, ru):
        r = vv / 2
        rv = 0.5
        pv = 0.5 * dF(r)
        ph = F(r)
        if r < 1e-12:
            # centre: phi_u = phi_v, rho = phi_v^2/e2 + V
            e2 = np.exp(2 * sig)
            rho = pv * pv / e2 + V(ph)
            psi = psi_of_y(8 * np.pi * G * rho / 3, alpha)
            m = 1 - (alpha * psi) ** 2
            return 0.0, r * 0.0
        e2 = np.exp(2 * sig)
        psi = (1 + 2 * ru * rv / e2) / (r * r)
        if r < 4 * h:
            # regular series near the centre: phi_u ~ phi_v there
            rho = pv * pv / e2 + V(ph)
            psi = psi_of_y(8 * np.pi * G * rho / 3, alpha)
        m = 1 - (alpha * psi) ** 2
        m32 = m * np.sqrt(max(m, 0.0))
        dsig = r * kappa * pv * pv * m32 / 2 / (2 * rv)        # from r_vv - 2 sigma_v r_v = C_v with r_vv = 0
        Tuv = e2 * V(ph)
        druv = r * (kappa * Tuv * m32 / 2 - e2 * psi * (3 * m - 2) / 2)
        return dsig, druv

    sig = -0.5 * np.log(2.0)          # e2 = 2 r_v^2 = 1/2 at the centre
    ru = -0.5
    S = np.empty(N + 1); S[0] = sig
    RU = np.empty(N + 1); RU[0] = ru
    dv = h / n_sub
    vv = 0.0
    for j in range(N):
        for _ in range(n_sub):
            k1 = deriv(vv, sig, ru)
            k2 = deriv(vv + dv / 2, sig + dv / 2 * k1[0], ru + dv / 2 * k1[1])
            k3 = deriv(vv + dv / 2, sig + dv / 2 * k2[0], ru + dv / 2 * k2[1])
            k4 = deriv(vv + dv, sig + dv * k3[0], ru + dv * k3[1])
            sig += dv / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            ru += dv / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            vv += dv
        S[j + 1] = sig; RU[j + 1] = ru
    grid.s[0] = S
    # mass along the initial ray: m_v = 4 pi r^2 [r_v V - e^{-2 sigma} r_u phi_v^2]  (trapezoid)
    r0 = v / 2
    pv0 = 0.5 * dF(r0)
    mv = 4 * np.pi * G * 0 + 4 * np.pi * r0 ** 2 * (0.5 * V(F(r0)) - np.exp(-2 * S) * RU * pv0 ** 2)
    Mm = np.concatenate([[0.0], np.cumsum(0.5 * (mv[1:] + mv[:-1]) * h)])
    grid.M[0] = Mm
    return RU


def evolve(grid, V, dV, alpha, G=1.0, n_iter=4, n_reg=2, u_max_index=None, monitor=None,
           stop_on_fail=True):
    """stop_on_fail=False ('causal masking'): a failed node is set to NaN and the
    evolution continues; NaN propagates only into the causal future of the failed
    node (its dependants), so every finite node is still a valid solution of the
    scheme.  grid.ok_until records the first failed anti-diagonal either way."""
    N, h = grid.N, grid.h
    grid.n_failed = 0
    R, Sg, P, Mg = grid.r, grid.s, grid.p, grid.M
    imax = N if u_max_index is None else u_max_index
    status = "COMPLETE"
    for s in range(1, 2 * N + 1):
        i_lo = max(1, s - N)
        i_hi = min(s // 2, imax)
        if i_lo > i_hi:
            continue
        # non-centre nodes: i < j
        ii = np.arange(i_lo, i_hi + 1)
        jj = s - ii
        nc = ii < jj
        i1, j1 = ii[nc], jj[nc]
        if i1.size:
            SW = (i1 - 1, j1 - 1); NW = (i1, j1 - 1); SE = (i1 - 1, j1)
            zs = []
            for A in (R, Sg, P):
                zs.append((A[SW], A[NW], A[SE]))
            Msw, Mnw, Mse = Mg[SW], Mg[NW], Mg[SE]
            Nm = Mnw + Mse - Msw
            # second-row-back values for a second-order r_u on the row u_i (mass flux)
            has2 = i1 >= 2
            i2 = np.where(has2, i1 - 2, 0)
            r2a = np.where(has2, R[i2, j1], np.nan); r2b = np.where(has2, R[i2, j1 - 1], np.nan)
            Nr = zs[0][1] + zs[0][2] - zs[0][0]
            Ns = zs[1][1] + zs[1][2] - zs[1][0]
            Np = zs[2][1] + zs[2][2] - zs[2][0]
            near = (j1 - i1) <= n_reg
            for _ in range(n_iter):
                rb = (zs[0][0] + zs[0][1] + zs[0][2] + Nr) / 4
                sb = (zs[1][0] + zs[1][1] + zs[1][2] + Ns) / 4
                pb = (zs[2][0] + zs[2][1] + zs[2][2] + Np) / 4
                ru = (zs[0][1] - zs[0][0] + Nr - zs[0][2]) / (2 * h)
                rv = (zs[0][2] - zs[0][0] + Nr - zs[0][1]) / (2 * h)
                su = (zs[1][1] - zs[1][0] + Ns - zs[1][2]) / (2 * h)
                sv = (zs[1][2] - zs[1][0] + Ns - zs[1][1]) / (2 * h)
                pu = (zs[2][1] - zs[2][0] + Np - zs[2][2]) / (2 * h)
                pv = (zs[2][2] - zs[2][0] + Np - zs[2][1]) / (2 * h)
                # psi from the generalized mass (keeps |alpha psi| < 1 by construction)
                Mb = (Msw + Mnw + Mse + Nm) / 4
                ov = psi_of_y(2 * G * Mb / rb ** 3, alpha)
                if near.any():
                    rho = local_density(sb[near], pu[near], pv[near], pb[near], V)
                    ov[near] = psi_of_y(8 * np.pi * G * rho / 3, alpha)
                Fr, Fs, Fp, psi, m = rhs(rb, sb, pb, ru, rv, su, sv, pu, pv, V, dV, alpha, G, psi_override=ov)
                # mass update along v on the new row u_i (flux law, midpoint)
                e2b = np.exp(-(zs[1][1] + Ns))          # e^{-2 sigma} at (u_i, v_{j-1/2})
                ru_row = np.where(has2,
                                  (3 * 0.5 * (Nr + zs[0][1]) - 4 * 0.5 * (zs[0][2] + zs[0][0]) + 0.5 * (r2a + r2b)) / (2 * h),
                                  ru)
                rmid = 0.5 * (zs[0][1] + Nr)
                rv_row = (Nr - zs[0][1]) / h
                pv_row = (Np - zs[2][1]) / h
                pmid = 0.5 * (zs[2][1] + Np)
                mv = 4 * np.pi * rmid ** 2 * (rv_row * V(pmid) - e2b * ru_row * pv_row ** 2)
                Nm = Mnw + h * mv
                Nr = zs[0][1] + zs[0][2] - zs[0][0] + h * h * Fr
                Ns = zs[1][1] + zs[1][2] - zs[1][0] + h * h * Fs
                Np = zs[2][1] + zs[2][2] - zs[2][0] + h * h * Fp
            bad = (~np.isfinite(Nr)) | (~np.isfinite(Ns)) | (~np.isfinite(Np)) | (m <= 0) | (Nr <= 0)
            if bad.any():
                if status == "COMPLETE":
                    grid.ok_until = s
                status = "FAILED"
                if stop_on_fail:
                    break
                Nr = np.where(bad, np.nan, Nr); Ns = np.where(bad, np.nan, Ns)
                Np = np.where(bad, np.nan, Np); Nm = np.where(bad, np.nan, Nm)
                grid.n_failed += int(np.sum(bad & np.all(np.isfinite([zs[0][0], zs[0][1], zs[0][2]]), axis=0)))
            R[i1, j1] = Nr; Sg[i1, j1] = Ns; P[i1, j1] = Np; Mg[i1, j1] = Nm
        # centre node
        if s % 2 == 0 and i_lo <= s // 2 <= i_hi:
            k = s // 2
            R[k, k] = 0.0
            Mg[k, k] = 0.0
            # sigma from the regularity condition e2 = 2 r_v^2 at r = 0, with
            # r_v from the two previous nodes on the ray v = v_k (second order):
            # r(k-1,k) = h r_v + h^2 r_uu/2,  r(k-2,k) = 2 h r_v + 2 h^2 r_uu
            if k >= 2:
                rv_c = (4 * R[k - 1, k] - R[k - 2, k]) / (2 * h)
            else:
                rv_c = R[k - 1, k] / h
            Sg[k, k] = 0.5 * np.log(2 * rv_c * rv_c)
            # phi even in r: quadratic extrapolation along u + v = const
            if k >= 2 and k + 2 <= N:
                P[k, k] = (4 * P[k - 1, k + 1] - P[k - 2, k + 2]) / 3
            elif k + 1 <= N:
                P[k, k] = P[k - 1, k + 1]
            else:
                P[k, k] = 2 * P[k - 1, k - 1] - P[k - 2, k - 2]
        if monitor is not None and s % 50 == 0:
            monitor(s)
    if status == "COMPLETE":
        grid.ok_until = 2 * N
    return status
