"""JAX port of the evolution part of network/kern.c (force, KDK leapfrog, drain of the massive modes), for GPUs.
Prototype: same equations, same component order (f[c] with c the Frobenius basis of kern.c), same scale-factor rules
(a = s tau; a_p frozen before TAU_C (U3b) or after TAU_F (U19b); tilt EPS_G with the physical a^2) and the same drain
schedule (every DAMP_EVERY steps, step index from t1).  Fields are stored in float32 as in kern.c; the arithmetic is float64
when x64 is enabled (default here, as kern.c), float32 otherwise.  Measurements stay in kern.c on the CPU (copy at snapshots).
Not a bit-identical port: XLA may reorder sums, so trajectories agree with kern.c to rounding at first and then diverge as
any chaotic evolution does; network/test_jax_kernel.py fixes the tolerance of the short-time agreement.

    import jax_kernel as JK
    f, p = JK.evolve(model, f, p, tau, dtau, s, nsteps, eps_g=0.0, tau_c=0.0, tau_f=0.0, damp_tau=0.0, damp_g=0.0)
    f, p = JK.evolve_lowmem(...)    # the same without the force buffer (U19c: 768^3 on an 80 GB GPU)
"""
from functools import partial
import math, os

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

LAM_Z, LAMP_Z, LAM_U, V2_U = 0.8, 0.48, 0.24, 0.5
ISQ2 = 0.70710678118654752440
OI, OJ = (0, 0, 0, 1, 1, 2), (1, 2, 3, 2, 3, 3)
DAMP_EVERY = 4


def _Q(q):
    """q[10] -> 4x4 nested list of arrays (symmetric)"""
    Q = [[None] * 4 for _ in range(4)]
    for i in range(4): Q[i][i] = q[i]
    for k in range(6):
        v = q[4 + k] * ISQ2; Q[OI[k]][OJ[k]] = v; Q[OJ[k]][OI[k]] = v
    return Q


def _mm(A, B):
    return [[A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] + A[i][3] * B[3][j] for j in range(4)] for i in range(4)]


def rp3_grad(q):
    """dV/dq for the 10 components (kern.c rp3_grad, without the tilt)"""
    Q = _Q(q); Q2 = _mm(Q, Q)
    tr = Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3]
    M = [[Q2[i][j] - Q[i][j] for j in range(4)] for i in range(4)]
    W = [[2 * Q[i][j] - (1.0 if i == j else 0.0) for j in range(4)] for i in range(4)]
    T = _mm(M, W)
    g = [0.5 * LAM_Z * (tr - 1) + 0.5 * LAMP_Z * T[i][i] for i in range(4)]
    g += [0.5 * LAMP_Z * 0.5 * (T[OI[k]][OJ[k]] + T[OJ[k]][OI[k]]) * 2 * ISQ2 for k in range(6)]
    return g


def u1_grad(q):
    d = q[0] * q[0] + q[1] * q[1] - V2_U
    return [LAM_U * d * q[0], LAM_U * d * q[1]]


SLAB_SITES = 1 << 22        # sites per z-slab: the float64 temporaries of force and drain scale with this, not with N^3


def _cz(N):
    """slab thickness: the largest divisor of N with cz N^2 <= SLAB_SITES (at least 1)"""
    return max([d for d in range(1, N + 1) if N % d == 0 and d * N * N <= SLAB_SITES] or [1])


def lap(xh, xd):
    """7-point periodic Laplacian of one float32 component on a z-slab with one halo plane on each side (xh: cz+2, N, N),
    as kern.c: the six neighbours are summed in float32 in the order x+1, x-1, y+1, y-1, z+1, z-1, and 6 x is subtracted in
    the compute precision (xd = the slab in it)"""
    x = xh[1:-1]
    nb = jnp.roll(x, -1, 2) + jnp.roll(x, 1, 2) + jnp.roll(x, -1, 1) + jnp.roll(x, 1, 1) + xh[2:] + xh[:-2]
    return nb.astype(xd.dtype) - 6.0 * xd


def _force_slab(model, fh, a, ap, eps_g, cdt):
    """F on one z-slab; fh is (C, cz+2, N, N) with the halo planes"""
    q = [fh[c, 1:-1].astype(cdt) for c in range(fh.shape[0])]
    g = rp3_grad(q) if model else u1_grad(q)
    a2, ap2 = a * a, ap * ap
    g0 = g[0] - eps_g + eps_g * (1.0 - a2 / ap2)          # kern.c: g[0] -= EPS_G; then += EPS_G (1 - a^2/a_p^2) if a_p != a
    g = [g0] + g[1:]
    return jnp.stack([(a2 * (lap(fh[c], q[c]) - ap2 * g[c])).astype(jnp.float32) for c in range(len(q))])


def force(model, f, a, ap, eps_g, cdt):
    """F = a^2 (lap f - a_p^2 dV), with the tilt on c = 0 multiplied by a^2 (kern.c force); f is (C, N, N, N) float32.
    Computed slab by slab in z and written in place, so the temporaries are those of one slab"""
    C, N = f.shape[0], f.shape[1]; cz = _cz(N)
    def body(k, out):
        idx = (k * cz - 1 + jnp.arange(cz + 2)) % N
        return jax.lax.dynamic_update_slice(out, _force_slab(model, jnp.take(f, idx, axis=1), a, ap, eps_g, cdt), (0, k * cz, 0, 0))
    return jax.lax.fori_loop(0, N // cz, body, jnp.zeros_like(f))


EIGH_BATCH = int(os.environ.get("WHALEFALL_EIGH_BATCH", 1 << 15))   # matrices per eigh call: cuSolver's batched solver refuses
                                                                   # the ~4e6 of one slab ("invalid value"); 2^15 works on an A100


def _jacobi4(A):
    """eigen-decomposition of a batch of symmetric 4x4 matrices (last two axes) with jnp.linalg.eigh, in chunks of EIGH_BATCH
    matrices (each matrix independently, so the chunking does not change the result)"""
    sh = A.shape; M = A.reshape(-1, 4, 4)
    if M.shape[0] <= EIGH_BATCH: return jnp.linalg.eigh(A)
    w, V = jax.lax.map(jnp.linalg.eigh, M, batch_size=EIGH_BATCH)
    return w.reshape(sh[:-1]), V.reshape(sh)


def damp_massive(model, f, p, fac, cdt, do=True):
    """kern.c damp_massive, slab by slab in z (the 4x4 float64 matrices of the RP3 branch exist for one slab at a time).
    do (traced): run the slab loop or not, a loop of 0 or N/cz trips instead of lax.cond, so that p is updated in place
    (on a GPU, a cond can keep a second copy of p: U19c, 768^3)"""
    C, N = f.shape[0], f.shape[1]; cz = _cz(N)
    def body(k, pp):
        fs = jax.lax.dynamic_slice(f, (0, k * cz, 0, 0), (C, cz, N, N)); ps = jax.lax.dynamic_slice(pp, (0, k * cz, 0, 0), (C, cz, N, N))
        return jax.lax.dynamic_update_slice(pp, _damp_local(model, fs, ps, fac, cdt), (0, k * cz, 0, 0))
    return jax.lax.fori_loop(0, jnp.where(do, N // cz, 0), body, p)


def drift(f, p, c, cdt):
    """f += (float)(c p) slab by slab in z, in place (the same arithmetic as the whole-array update of evolve)"""
    C, N = f.shape[0], f.shape[1]; cz = _cz(N)
    def body(k, ff):
        fs = jax.lax.dynamic_slice(ff, (0, k * cz, 0, 0), (C, cz, N, N)); ps = jax.lax.dynamic_slice(p, (0, k * cz, 0, 0), (C, cz, N, N))
        return jax.lax.dynamic_update_slice(ff, fs + (c * ps.astype(cdt)).astype(jnp.float32), (0, k * cz, 0, 0))
    return jax.lax.fori_loop(0, N // cz, body, f)


def _damp_local(model, f, p, fac, cdt):
    """the drain on a block of sites (pointwise, no neighbours)"""
    if model:
        q = jnp.moveaxis(f.astype(cdt), 0, -1); pc = jnp.moveaxis(p.astype(cdt), 0, -1)
        def mat(v):
            Q = jnp.zeros(v.shape[:-1] + (4, 4), v.dtype)
            for i in range(4): Q = Q.at[..., i, i].set(v[..., i])
            for k in range(6):
                Q = Q.at[..., OI[k], OJ[k]].set(v[..., 4 + k] * ISQ2); Q = Q.at[..., OJ[k], OI[k]].set(v[..., 4 + k] * ISQ2)
            return Q
        w, Vv = _jacobi4(mat(q))                           # ascending eigenvalues
        n = Vv[..., :, 3]; gap = w[..., 3] - w[..., 2]
        P = mat(pc); u = jnp.einsum("...rc,...c->...r", P, n); nPn = jnp.einsum("...r,...r->...", n, u)
        wv = u - nPn[..., None] * n
        pt = wv[..., :, None] * n[..., None, :] + n[..., :, None] * wv[..., None, :]
        Pn = pt + fac * (P - pt)
        new = jnp.stack([Pn[..., k, k] for k in range(4)] + [Pn[..., OI[k], OJ[k]] * 2 * ISQ2 for k in range(6)], 0)
        return jnp.where((gap > 0.6)[None], new.astype(jnp.float32), p)
    f0, f1 = f[0].astype(cdt), f[1].astype(cdt); r2 = f0 * f0 + f1 * f1
    rad = (f0 * p[0].astype(cdt) + f1 * p[1].astype(cdt)) / r2
    new = jnp.stack([p[0] - ((1 - fac) * rad * f0).astype(jnp.float32), p[1] - ((1 - fac) * rad * f1).astype(jnp.float32)])
    return jnp.where((r2 > 0.36 * V2_U)[None], new, p)


@partial(jax.jit, static_argnames=("model", "x64"), donate_argnums=(1, 2))
def evolve(model, f, p, tau, dtau, s, nsteps, eps_g=0.0, tau_c=0.0, tau_f=0.0, damp_tau=0.0, damp_g=0.0, x64=True):
    """nsteps of the KDK leapfrog of kern.c evolve; returns (f, p).  s <= 0 means Minkowski (a = 1)"""
    cdt = jnp.float64 if x64 else jnp.float32
    tau, dtau, s, eps_g, tau_c, tau_f, damp_tau, damp_g = (jnp.asarray(v, cdt) for v in (tau, dtau, s, eps_g, tau_c, tau_f, damp_tau, damp_g))
    def aof(t): return jnp.where(s > 0, s * t, 1.0)
    def apf(t):
        return jnp.where((s > 0) & (tau_c > 0) & (t < tau_c), s * tau_c,
                         jnp.where((s > 0) & (tau_f > 0) & (t > tau_f), s * tau_f, aof(t)))
    F0 = force(model, f, aof(tau), apf(tau), eps_g, cdt)
    def step(k, st):
        f, p, F = st
        t0 = tau + k * dtau; th = t0 + 0.5 * dtau; t1 = t0 + dtau
        ih = 1.0 / (aof(th) * aof(th))
        p = p + (0.5 * dtau * F.astype(cdt)).astype(jnp.float32)              # kern.c: p += (float)(0.5 dtau F), float add
        f = f + (dtau * ih * p.astype(cdt)).astype(jnp.float32)
        F = force(model, f, aof(t1), apf(t1), eps_g, cdt)
        p = p + (0.5 * dtau * F.astype(cdt)).astype(jnp.float32)
        do = (damp_tau > 0) & (t1 < damp_tau) & (jnp.round(t1 / dtau).astype(jnp.int64) % DAMP_EVERY == 0)
        p = jax.lax.cond(do, lambda pp: damp_massive(model, f, pp, jnp.exp(-damp_g * aof(t1) * dtau * DAMP_EVERY), cdt), lambda pp: pp, p)
        return f, p, F
    f, p, _ = jax.lax.fori_loop(0, nsteps, step, (f, p, F0))
    return f, p


def kick(model, f, p, a, ap, eps_g, cdt, coef):
    """p += (float)(coef F(f)) slab by slab in z, without storing F (kern.c kick, U19c): F of one slab is that of force()"""
    C, N = f.shape[0], f.shape[1]; cz = _cz(N)
    def body(k, pp):
        idx = (k * cz - 1 + jnp.arange(cz + 2)) % N
        F = _force_slab(model, jnp.take(f, idx, axis=1), a, ap, eps_g, cdt)
        ps = jax.lax.dynamic_slice(pp, (0, k * cz, 0, 0), (C, cz, N, N))
        return jax.lax.dynamic_update_slice(pp, ps + (coef * F.astype(cdt)).astype(jnp.float32), (0, k * cz, 0, 0))
    return jax.lax.fori_loop(0, N // cz, body, p)


@partial(jax.jit, static_argnames=("model", "x64"), donate_argnums=(1, 2))
def evolve_lowmem(model, f, p, tau, dtau, s, nsteps, eps_g=0.0, tau_c=0.0, tau_f=0.0, damp_tau=0.0, damp_g=0.0, x64=True):
    """kern.c evolve_lowmem (U19c): the KDK leapfrog of evolve() without carrying the force between steps (state f and p only,
    twice the force evaluations); the same arithmetic, so the same result as evolve().  Every full-lattice update (kick, drift,
    drain) is a slab loop writing in place, so the device holds f, p and the temporaries of one slab"""
    cdt = jnp.float64 if x64 else jnp.float32
    tau, dtau, s, eps_g, tau_c, tau_f, damp_tau, damp_g = (jnp.asarray(v, cdt) for v in (tau, dtau, s, eps_g, tau_c, tau_f, damp_tau, damp_g))
    def aof(t): return jnp.where(s > 0, s * t, 1.0)
    def apf(t):
        return jnp.where((s > 0) & (tau_c > 0) & (t < tau_c), s * tau_c,
                         jnp.where((s > 0) & (tau_f > 0) & (t > tau_f), s * tau_f, aof(t)))
    def step(k, st):
        f, p, tf = st
        t0 = tau + k * dtau; th = t0 + 0.5 * dtau; t1 = t0 + dtau
        ih = 1.0 / (aof(th) * aof(th))
        p = kick(model, f, p, aof(tf), apf(tf), eps_g, cdt, 0.5 * dtau)
        f = drift(f, p, dtau * ih, cdt)
        p = kick(model, f, p, aof(t1), apf(t1), eps_g, cdt, 0.5 * dtau)
        do = (damp_tau > 0) & (t1 < damp_tau) & (jnp.round(t1 / dtau).astype(jnp.int64) % DAMP_EVERY == 0)
        p = damp_massive(model, f, p, jnp.exp(-damp_g * aof(t1) * dtau * DAMP_EVERY), cdt, do)
        return f, p, t1
    f, p, _ = jax.lax.fori_loop(0, nsteps, step, (f, p, tau))
    return f, p


def bytes_per_network(model, N, keep_force=True):
    """device memory of the state (f, p and the force carried between steps), float32; the slab temporaries come on top
    (of order 100 float64 arrays of SLAB_SITES sites, a few GB)"""
    C = 10 if model else 2
    return (3 if keep_force else 2) * C * N ** 3 * 4


if __name__ == "__main__":
    for N in (512, 640, 768, 1024):
        print(f"N={N}: RP3 state {bytes_per_network(1, N) / 1e9:.1f} GB, U(1) {bytes_per_network(0, N) / 1e9:.1f} GB "
              f"(plus XLA temporaries); host copy for the kern.c measurements {2 * 10 * N ** 3 * 4 / 1e9:.1f} GB")
