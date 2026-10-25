#!/usr/bin/env python3
"""C32b — Do WHALEFALL scars FORM after inflation?  Kibble test with the box and the operator fixed.

Why C32b (see claude/REVISION_C31_C32.md): C32's null (S1: seg_max = 0) is not conclusive.
  (1) Box: the linear analysis shows the correlation wavelength at breaking is L_c ~ 23-35 comoving
      (~ 2-4 Hubble radii physical).  C32's box was 6.4 comoving, a fraction of ONE domain.
  (2) Initial state: white 1e-3 noise per cell.  Here the lattice starts from the LINEAR evolution of the
      adiabatic vacuum spectrum through reheating (Fourier, exact while |D|^2 << 1), switched to the lattice
      when the band variance reaches <|D|^2> = 0.01.
  (3) Operator: the "natural" xi of R + R^2 couples to the JORDAN curvature.  For the Einstein-frame canonical
      field (g_E = F g_J, F = e^{K chi}) the conformal operator identity gives EXACTLY
          m_R^2(t) = R_E/6 + (xi - 1/6) R_J/F,   R_E = 4V - chi'^2,  R_J/F = 3 M^2 (1 - e^{-K chi}).
      Arms 'J' use this; arm 'E' keeps C32's xi R_E for attribution.
  (4) Gate: formation (peak string-length density in units of L_c), not survival in a sub-Hubble box.
Declared simplifications:
  - Q sector potential kept in Jordan form (factor F on the Q potential dropped).  At the switch |F-1| < 2 %
    (reported per arm); the linear stage keeps the exact 1/F on the tachyonic mass.
  - Fluctuation amplitude eps (vacuum normalisation of Q in code units) is not derived: eps = 1e-3 (main) and
    1e-6 (sensitivity arm).
  - The linear stage treats the 9 traceless directions of D = Q - x0 I as independent (exact at quadratic order);
    the heavy trace mode starts at the isotropic stationary point x0 with no fluctuation.
  - Lattice in float32, coarse: dx_phys at the end ~ 6-9 (k_Nyquist ~ 3x the tachyonic scale m_t = 0.139).
    The massive radial modes (m = 0.49, 1.36) are under-resolved; defect cores are ~1 cell.  Topology is read
    with the parent's trusted plaquette flux (same thresholds as v0_72r_topology.compute_flux).
  - Prescribed Starobinsky background a(t), R(t) (Q does not backreact).  Parent RP^3 (m_G = 0), no baths.
Units: code units of the parent (ZQ = 400, LP = 192, lam = 320 [p4_hermite, s = 0]); time in code units,
  M_code = scalaron mass in code units (mu = m_Q/M = 0.49/M_code; 0.98 = tuned C26 mass).

Arms (seed 3202 for all lattice realisations; the control uses the J025 background):
  J025      xi = 0.25, J, M_code 0.98, eps 1e-3, box 3 L_c   (main; heavy during inflation, as C31 requires)
  C0        ordered start (uniform n n^T vacuum) + the J025 fluctuation field, J025 background (detector control)
  J100      xi = 1.0,  J, M_code 0.98, eps 1e-3
  E100      xi = 1.0,  E, M_code 0.98, eps 1e-3   (C32 S2 operator, new box: attribution)
  J025hr    as J025 with box 1.5 L_c (dx / 2): resolution check
  J025e6    xi = 0.25, J, M_code 0.98, eps 1e-6   (Kibble-Zurek amplitude sensitivity)
  J100mu01  xi = 1.0,  J, M_code 4.9,  eps 1e-3   (mu = 0.1: one step towards the light-Q regime of C31)
N = 192 (fallback 160 on GPU OOM, declared), box = 3 L_c (L_c = 2 pi / k_peak of k^2 P(k) at the switch).
Lattice from t_switch to t_switch + 80 (code units), diagnostics every 1.0; dual-graph components on the first
10 ordered samples and every 5 samples.

Smoke tests (development, numpy, not evidence): the tachyonic roll overshoots and the amplitude "breathes"
(order 0.4 <-> 1.3); raw segment counts GROW after formation (Goldstone pumping by the breathing, cf. C18, m_G = 0).
So the gate is taken at the FORMATION epoch, and the late growth is reported separately.

Frozen labels (written before any production data):
  t_ord = first sample with >= 90 % ordered cells; seg_form = max raw segments in [t_ord, t_ord + 5];
  Lambda_form = seg_form * dx * L_c^2 / box^3 (comoving string length per L_c^2);
  persist25 = seg(t_ord + 25) / seg_form; tangle_growth_late = <seg over the last 10> / seg_form.
  C32B_CONTROL_FAIL : C0 seg_max >= 20
  C32B_SCARS_FORM   : J025 Lambda_form >= 0.3 and persist25 >= 0.3 and bianchi_max / seg_max <= 0.01
  C32B_NO_SCARS     : J025 and J100 Lambda_form < 0.03
  C32B_MIXED        : otherwise
  Secondary: C32_NULL_WAS_BOX (E100 Lambda_form >= 0.3); RESOLUTION_OK (J025hr vs J025 Lambda_form within x2);
             TANGLE_RESOLVED (tangle_growth_late within x2 between J025 and J025hr); TANGLE_GROWS_J025 (> 2);
             KZ_WEAK (J025e6 vs J025 Lambda_form within x3); MU01_FORMS (J100mu01 Lambda_form >= 0.3).
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
try:
    import jax
    jax.config.update('jax_enable_x64', False)
    import jax.numpy as xp
    HAVE_JAX = True
except Exception:  # numpy fallback (smoke tests only)
    xp = np; HAVE_JAX = False

# ---------------- parent constants (v0_66 / v0_72 p4_hermite at s = 0) ----------------
ZQ = 400.0; LP = 192.0; LBASE = 64.0; LAM = 320.0
ORDER_D2_THRESHOLD = 0.3850195; ORDER_GAP_THRESHOLD = 0.20
K = math.sqrt(2 / 3)
ARMS = {
    'J025':     dict(xi=0.25, op='J', Mc=0.98, eps=1e-3, start='restored'),
    'J100':     dict(xi=1.0,  op='J', Mc=0.98, eps=1e-3, start='restored'),
    'E100':     dict(xi=1.0,  op='E', Mc=0.98, eps=1e-3, start='restored'),
    'J025e6':   dict(xi=0.25, op='J', Mc=0.98, eps=1e-6, start='restored'),
    'J100mu01': dict(xi=1.0,  op='J', Mc=4.9,  eps=1e-3, start='restored'),
    'C0':       dict(xi=0.25, op='J', Mc=0.98, eps=1e-3, start='ordered'),
    'J025hr':   dict(xi=0.25, op='J', Mc=0.98, eps=1e-3, start='restored', box_lc=1.5),   # resolution: dx / 2
}
SWITCH_VAR = 0.01; BOX_LC = 3.0; BOX_OVERRIDE = None; T_LAT = 80.0; SAMPLE = 1.0; SEED = 3202


def x0_iso():
    r = np.roots([2 * LP, -3 * LP, LP + 4 * LAM, -LAM]); return float(r[np.abs(r.imag) < 1e-10].real[0])


X0 = x0_iso()
MT2 = -0.019416   # tachyonic mass^2 of the 9 traceless directions at x0 I (Hessian / ZQ; checked in selftest)


def V(c): return 0.75 * (1 - np.exp(-K * c)) ** 2
def dV(c): return 1.5 * K * np.exp(-K * c) * (1 - np.exp(-K * c))


def background(Mc, Tcode):
    """Starobinsky from the end of inflation (eps_H = 1), in M = M_Pl = 1 units, returned in code units."""
    c0 = brentq(lambda c: 0.5 * (dV(c) / V(c)) ** 2 - 1, 0.1, 5); cd0 = -dV(c0) / (3 * math.sqrt(V(c0) / 3))
    def rhs(t, y):
        c, cd, la = y; H = math.sqrt((0.5 * cd * cd + V(c)) / 3); return [cd, -3 * H * cd - dV(c), H]
    s = solve_ivp(rhs, (0, Tcode * Mc), [c0, cd0, 0.0], rtol=1e-10, atol=1e-12, dense_output=True, max_step=0.02)
    return s


def bg_parts(s, Mc, tc, xi, op):
    c, cd, la = s.sol(tc * Mc); a = math.exp(la); H = math.sqrt((0.5 * cd * cd + V(c)) / 3) * Mc
    RE = (4 * V(c) - cd * cd) * Mc ** 2; RJF = 3 * (1 - math.exp(-K * c)) * Mc ** 2; F = math.exp(K * c)
    mR2 = RE / 6 + (xi - 1 / 6) * RJF if op == 'J' else xi * RE
    return a, H, mR2, F, RE, RJF


def linear_stage(arm, kmax):
    """Mode functions of one traceless component of D from the adiabatic vacuum at the end of inflation."""
    p = ARMS[arm]; xi, op, Mc, eps = p['xi'], p['op'], p['Mc'], p['eps']
    s = background(Mc, 400.0)
    ks = np.concatenate([np.linspace(0.002, 3.0, 700), np.linspace(3.0, kmax, 200)[1:]])
    n = len(ks)
    a0, H0, m0, F0, _, _ = bg_parts(s, Mc, 0.0, xi, op); m0 += MT2 / F0
    if m0 <= 0: raise RuntimeError(f'{arm}: not restored at the end of inflation (m2 = {m0})')
    w0 = np.sqrt(ks ** 2 + m0)
    y0 = np.concatenate([eps / np.sqrt(2 * w0), np.zeros(n), np.zeros(n), -eps * w0 / np.sqrt(2 * w0)])
    def rhs(t, y):
        a, H, mR2, F, _, _ = bg_parts(s, Mc, t, xi, op); W = ks ** 2 / a ** 2 + mR2 + MT2 / F
        return np.concatenate([y[2 * n:3 * n], y[3 * n:], -3 * H * y[2 * n:3 * n] - W * y[:n], -3 * H * y[3 * n:] - W * y[n:2 * n]])
    def var(y): return 9 * np.trapezoid(ks ** 2 * (y[:n] ** 2 + y[n:2 * n] ** 2), ks) / (2 * np.pi ** 2)
    ev = lambda t, y: var(y) - SWITCH_VAR; ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0, 400.0), y0, events=ev, rtol=1e-7, atol=1e-16, max_step=0.1 / Mc)
    if not len(sol.t_events[0]): raise RuntimeError(f'{arm}: no breaking before t = 400')
    y = sol.y[:, -1]; tsw = float(sol.t[-1])
    u = y[:n] + 1j * y[n:2 * n]; ud = y[2 * n:3 * n] + 1j * y[3 * n:]
    P = np.abs(u) ** 2; kp = float(ks[np.argmax(ks ** 2 * P)])
    a, H, mR2, F, RE, RJF = bg_parts(s, Mc, tsw, xi, op)
    info = dict(t_switch=tsw, a_switch=a, H_switch=H, hubble_radius_phys=1 / H, F_switch=F, F_minus_1=F - 1,
                mR2_switch=mR2, k_peak_com=kp, Lc_com=2 * np.pi / kp, Lc_phys=2 * np.pi / kp * a,
                Lc_over_hubble=2 * np.pi / kp * a * H, m2_start=m0)
    return s, ks, u, ud, info


def traceless_basis():
    B = []
    for i in range(4):
        for j in range(i + 1, 4):
            E = np.zeros((4, 4)); E[i, j] = E[j, i] = 1 / math.sqrt(2); B.append(E)
    for d in ([1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]):
        v = np.array(d, float); B.append(np.diag(v / np.linalg.norm(v)))
    return np.array(B)   # (9,4,4) Frobenius-orthonormal, traceless


def lattice_ic(ks, u, ud, N, box, a_sw, start, rng):
    """Gaussian realisation: D~_k = u_k zeta_k, Pi~_k = a^3 ZQ ud_k zeta_k (same zeta), sqrt2 Re."""
    kx = 2 * np.pi * np.fft.fftfreq(N, d=box / N)
    KK = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kx[None, None, :] ** 2)
    ur = np.interp(KK, ks, u.real, right=0.0); ui = np.interp(KK, ks, u.imag, right=0.0)
    vr = np.interp(KK, ks, ud.real, right=0.0); vi = np.interp(KK, ks, ud.imag, right=0.0)
    norm = N ** 3 / math.sqrt(box ** 3)
    B = traceless_basis(); D = np.zeros((N, N, N, 4, 4)); Dd = np.zeros_like(D)
    for b in range(9):
        z = (rng.normal(size=KK.shape) + 1j * rng.normal(size=KK.shape)) / math.sqrt(2)
        z[0, 0, 0] = 0.0   # no homogeneous component in a finite box (declared)
        d = math.sqrt(2) * np.real(np.fft.ifftn(norm * (ur + 1j * ui) * z))
        dd = math.sqrt(2) * np.real(np.fft.ifftn(norm * (vr + 1j * vi) * z))
        D += d[..., None, None] * B[b]; Dd += dd[..., None, None] * B[b]
    if start == 'restored':
        Q = X0 * np.eye(4) + D
    else:
        n = np.zeros(4); n[0] = 1.0; Q = np.outer(n, n) + D
    return Q, a_sw ** 3 * ZQ * Dd


def grad_U(Q):
    I = xp.eye(4, dtype=Q.dtype); tr = xp.trace(Q, axis1=-2, axis2=-1) - 1.0; A = Q @ Q - Q
    return (LAM / 2.0) * tr[..., None, None] * I + (LP / 2.0) * (A @ (2.0 * Q - I))


def U_density(Q):
    tr = xp.trace(Q, axis1=-2, axis2=-1) - 1.0; A = Q @ Q - Q
    return (LAM / 4.0) * tr * tr + (LP / 4.0) * xp.sum(A * A, axis=(-2, -1))


def lap(Q, dx):
    out = -6.0 * Q
    for ax in range(3): out = out + xp.roll(Q, 1, axis=ax) + xp.roll(Q, -1, axis=ax)
    return out / (dx * dx)


def sym(x): return 0.5 * (x + xp.swapaxes(x, -1, -2))


# GPU eigh fix 2 (main conversation, 26-sep ~17:55): one jitted eigh over chunks of 2^20 matrices asked cuSolver
# syevj for a 4.1 GB workspace (OOM at N192 and N160, which then poisoned the CUDA context).  Now: a small jitted eigh
# per chunk of 2^16 matrices, called from Python, keeping only the top eigenvalue, the gap and the top eigenvector.
# Same algorithm, same thresholds, same result (the flux uses only these three quantities).
_EIG_CHUNK = 1 << 16
_JIT = {}


def _top_eig(Qf, xp_):
    vals, vec = xp_.linalg.eigh(Qf)
    return vals[:, -1], vals[:, -1] - vals[:, -2], vec[:, :, -1]


def top_eig(Q, xp_):
    shp = tuple(Q.shape[:-2]); Qf = Q.reshape((-1, 4, 4)); n = Qf.shape[0]
    if xp_ is np or not HAVE_JAX:
        l1, g, d = _top_eig(Qf, xp_)
    else:
        if 'eig' not in _JIT: _JIT['eig'] = jax.jit(lambda q: _top_eig(q, xp))
        C = _EIG_CHUNK; pad = (-n) % C
        if pad: Qf = xp.concatenate([Qf, xp.broadcast_to(xp.eye(4, dtype=Qf.dtype), (pad, 4, 4))])
        outs = []
        for i in range(0, n + pad, C):
            o = _JIT['eig'](Qf[i:i + C]); jax.block_until_ready(o); outs.append(o)
        l1, g, d = (xp.concatenate([o[k] for o in outs])[:n] for k in range(3))
    return l1.reshape(shp), g.reshape(shp), d.reshape(shp + (4,))


def flux_segments(Q, m=None):
    """Same algorithm and thresholds as tools/v0_72r_topology.compute_flux + v0_75_energy.trusted_defect_support."""
    xp_ = m if m is not None else xp
    l1, gap, d = top_eig(Q, xp_)
    if HAVE_JAX and xp_ is not np:
        if 'flux' not in _JIT: _JIT['flux'] = jax.jit(lambda Q, l1, gap, d: _flux_from(Q, l1, gap, d, xp))
        return _JIT['flux'](Q, l1, gap, d)
    return _flux_from(Q, l1, gap, d, xp_)


def _flux_from(Q, l1, gap, d, xp_):
    d2 = xp_.maximum(xp_.sum(Q * Q, axis=(-2, -1)) + 1.0 - 2.0 * l1, 0.0)
    ordered = (d2 < ORDER_D2_THRESHOLD) & (gap > ORDER_GAP_THRESHOLD)
    s = [xp_.where(xp_.sum(d * xp_.roll(d, -1, axis=ax), axis=-1) >= 0, 1, -1) for ax in range(3)]
    sx, sy, sz = s
    fxy = sx * xp_.roll(sy, -1, 0) * xp_.roll(sx, -1, 1) * sy
    fyz = sy * xp_.roll(sz, -1, 1) * xp_.roll(sy, -1, 2) * sz
    fzx = sz * xp_.roll(sx, -1, 2) * xp_.roll(sz, -1, 0) * sx
    def pm(a, b): return ordered & xp_.roll(ordered, -1, a) & xp_.roll(ordered, -1, b) & xp_.roll(xp_.roll(ordered, -1, a), -1, b)
    seg = xp_.sum((fxy < 0) & pm(0, 1)) + xp_.sum((fyz < 0) & pm(1, 2)) + xp_.sum((fzx < 0) & pm(2, 0))
    co = ordered
    for ax in range(3): co = co & xp_.roll(ordered, -1, axis=ax)
    co = co & xp_.roll(xp_.roll(ordered, -1, 0), -1, 1) & xp_.roll(xp_.roll(ordered, -1, 1), -1, 2) & xp_.roll(xp_.roll(ordered, -1, 2), -1, 0)
    co = co & xp_.roll(xp_.roll(xp_.roll(ordered, -1, 0), -1, 1), -1, 2)
    cube = fxy * xp_.roll(fxy, -1, 2) * fyz * xp_.roll(fyz, -1, 0) * fzx * xp_.roll(fzx, -1, 1)
    bianchi = xp_.sum((cube < 0) & co)
    return seg, bianchi, xp_.mean(ordered.astype(Q.dtype))


def string_components(Q, min_edges):
    """Dual-graph connected components of the trusted negative-flux plaquettes (as tools/v0_72r_topology.dual_graph).
    Returns counts; 'long' = components with >= min_edges edges (one correlation wavelength L_c / dx)."""
    import networkx as nx
    Q = np.asarray(Q, np.float64); N = Q.shape[0]
    vals, vec = np.linalg.eigh(Q)
    d2 = np.maximum(np.sum(Q * Q, axis=(-2, -1)) + 1.0 - 2.0 * vals[..., -1], 0.0); gap = vals[..., -1] - vals[..., -2]
    o = (d2 < ORDER_D2_THRESHOLD) & (gap > ORDER_GAP_THRESHOLD); d = vec[..., -1]
    sx, sy, sz = [np.where(np.sum(d * np.roll(d, -1, axis=ax), axis=-1) >= 0, 1, -1) for ax in range(3)]
    fxy = sx * np.roll(sy, -1, 0) * np.roll(sx, -1, 1) * sy
    fyz = sy * np.roll(sz, -1, 1) * np.roll(sy, -1, 2) * sz
    fzx = sz * np.roll(sx, -1, 2) * np.roll(sz, -1, 0) * sx
    pm = lambda a, b: o & np.roll(o, -1, a) & np.roll(o, -1, b) & np.roll(np.roll(o, -1, a), -1, b)
    G = nx.Graph()
    for i, j, k in zip(*np.where((fxy < 0) & pm(0, 1))): G.add_edge((i, j, (k - 1) % N), (i, j, k))
    for i, j, k in zip(*np.where((fyz < 0) & pm(1, 2))): G.add_edge(((i - 1) % N, j, k), (i, j, k))
    for i, j, k in zip(*np.where((fzx < 0) & pm(2, 0))): G.add_edge((i, (j - 1) % N, k), (i, j, k))
    sizes = sorted((G.subgraph(c).number_of_edges() for c in nx.connected_components(G)), reverse=True)
    tot = int(sum(sizes)); lng = int(sum(x for x in sizes if x >= min_edges))
    return dict(n_components=len(sizes), max_component=int(sizes[0]) if sizes else 0, long_segments=lng,
                long_frac=(lng / tot if tot else 0.0), small_loops_le8=int(sum(1 for x in sizes if x <= 8)))


def spectrum_peak(Q, box):
    """k of the maximum of k^2 P(k) of D = Q - <Q> (all components), comoving."""
    Qn = np.asarray(Q, np.float64); D = Qn - Qn.mean(axis=(0, 1, 2))
    N = D.shape[0]; P = np.zeros((N, N, N))
    for i in range(4):
        for j in range(i, 4):
            P += np.abs(np.fft.fftn(D[..., i, j])) ** 2 * (1 if i == j else 2)
    kx = 2 * np.pi * np.fft.fftfreq(N, d=box / N)
    KK = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kx[None, None, :] ** 2)
    dk = 2 * np.pi / box; idx = np.rint(KK / dk).astype(int); nb = idx.max() + 1
    pk = np.bincount(idx.ravel(), P.ravel(), nb) / np.maximum(np.bincount(idx.ravel(), None, nb), 1)
    kb = np.arange(nb) * dk; w = kb ** 2 * pk; w[0] = 0
    return float(kb[np.argmax(w)])


def run_arm(arm, N, out, T=T_LAT, lin_cache=None):
    p = ARMS[arm]; bgarm = 'J025' if arm == 'C0' else arm
    t0w = time.perf_counter()
    if bgarm == 'J025hr': bgarm_lin = 'J025'
    else: bgarm_lin = bgarm
    blc = BOX_OVERRIDE if BOX_OVERRIDE is not None else p.get('box_lc', BOX_LC)
    s, ks, u, ud, info = linear_stage(bgarm_lin, kmax=4.0)
    box = blc * info['Lc_com']; kmax = math.sqrt(3) * math.pi * N / box * 1.01
    if kmax > 4.0: s, ks, u, ud, info = linear_stage(bgarm_lin, kmax=kmax)
    bgarm = bgarm_lin
    dx = box / N; Mc = ARMS[bgarm]['Mc']; xi, op = ARMS[bgarm]['xi'], ARMS[bgarm]['op']
    tsw = info['t_switch']; T1 = tsw + T
    dt = 0.05 if Mc < 2 else 0.02
    # background tables for the lattice
    tt = np.arange(tsw - 1.0, T1 + 1.0, 0.005)
    tab = np.array([bg_parts(s, Mc, t, xi, op) for t in tt])
    a_tab, H_tab, mR2_tab = tab[:, 0], tab[:, 1], tab[:, 2]
    rng = np.random.default_rng(SEED)
    Q0, Pi0 = lattice_ic(ks, u, ud, N, box, info['a_switch'], p['start'], rng)
    dtype = xp.float32
    Q = xp.asarray(Q0, dtype=dtype); Pi = xp.asarray(Pi0, dtype=dtype); del Q0, Pi0
    I0 = xp.asarray(X0 * np.eye(4), dtype=dtype)
    tj = xp.asarray(tt); aj = xp.asarray(a_tab); mj = xp.asarray(mR2_tab)
    a_of = lambda t: xp.interp(t, tj, aj).astype(dtype); m_of = lambda t: xp.interp(t, tj, mj).astype(dtype)

    def force(Q, t):
        a = a_of(t)
        return a * ZQ * lap(Q, dx) - a ** 3 * (grad_U(Q) + ZQ * m_of(t) * (Q - I0))

    def step(Q, Pi, t):
        Pi = Pi + 0.5 * dt * force(Q, t)
        Q = sym(Q + dt * Pi / (a_of(t + 0.5 * dt) ** 3 * ZQ))
        return Q, Pi + 0.5 * dt * force(Q, t + dt)

    nsub = int(round(SAMPLE / dt))
    if HAVE_JAX:
        from jax import lax
        @jax.jit
        def advance(Q, Pi, t0):
            def body(i, st): return step(st[0], st[1], t0 + i * dt)
            return lax.fori_loop(0, nsub, body, (Q, Pi))
        _mf = flux_segments   # chunked eigh from Python + jitted flux (fix 2); not jit-able as a whole
        def meas_flux(Q):
            try:
                r = _mf(Q); jax.block_until_ready(r); return r   # materialise inside the try (async errors)
            except Exception as e:   # e.g. GPU OOM in batched eigh: fall back to host numpy (same algorithm)
                print('flux on host:', type(e).__name__, flush=True); return flux_segments(np.asarray(Q, np.float64), np)
    else:
        def advance(Q, Pi, t0):
            for i in range(nsub): Q, Pi = step(Q, Pi, t0 + i * dt)
            return Q, Pi
        meas_flux = flux_segments

    def measure(Q, Pi, t):
        a = float(np.interp(t, tt, a_tab))
        Vd = Pi / (a ** 3 * ZQ); Kd = 0.5 * ZQ * xp.sum(Vd * Vd, axis=(-2, -1))
        Gd = sum(0.5 * ZQ * xp.sum((xp.roll(Q, -1, axis=ax) - Q) ** 2, axis=(-2, -1)) for ax in range(3)) / (a * a * dx * dx)
        Ud = U_density(Q); D = Q - I0
        seg, bia, ordf = meas_flux(Q)
        return dict(time=t, a=a, mR2=float(np.interp(t, tt, mR2_tab)), rho=float(xp.mean(Kd + Gd + Ud)),
                    p=float(xp.mean(Kd - Gd / 3 - Ud)), order=float(xp.mean(xp.sum(D * D, axis=(-2, -1)))) / 0.75,
                    ordered_frac=float(ordf), segments=int(seg), bianchi=int(bia))

    od = Path(out) / arm; od.mkdir(parents=True, exist_ok=True)
    rows = []; t = tsw; kpk = {}; n_ord = 0; minE = max(8, int(round(info['Lc_com'] / dx)))
    kpk['switch'] = spectrum_peak(np.asarray(Q), box)
    nsamp = int(round(T / SAMPLE))
    for i in range(nsamp + 1):
        rows.append(measure(Q, Pi, t))
        if rows[-1]['ordered_frac'] >= 0.9: n_ord += 1
        if (1 <= n_ord <= 10) or i % 5 == 0 or i == nsamp:
            rows[-1].update(string_components(Q, min_edges=minE) if rows[-1]['segments'] > 0
                            else dict(n_components=0, max_component=0, long_segments=0, long_frac=0.0, small_loops_le8=0))
        if i % 5 == 0 or i == nsamp:
            r = rows[-1]
            print(f"[{arm}] t={t:.1f} a={r['a']:.2f} order={r['order']:.3f} ordf={r['ordered_frac']:.3f} seg={r['segments']} "
                  f"bia={r['bianchi']} long={r.get('long_segments')} maxc={r.get('max_component')} el={time.perf_counter()-t0w:.0f}s", flush=True)
        if i == nsamp: break
        Q, Pi = advance(Q, Pi, xp.asarray(t)); t = tsw + (i + 1) * SAMPLE
        if HAVE_JAX: jax.block_until_ready(Q)
        if not np.isfinite(rows[-1]['rho']): raise RuntimeError('non-finite energy')
        if i + 1 == int(round(20 / SAMPLE)): kpk['switch_plus_20'] = spectrum_peak(np.asarray(Q), box)
    import pandas as pd
    df = pd.DataFrame(rows); df.to_csv(od / 'trace.csv', index=False)
    seg = df.segments.values; tm = df.time.values; Lc = info['Lc_com']; V3 = box ** 3
    lam = lambda n: n * dx * Lc ** 2 / V3          # comoving string length per L_c^2
    iord = np.where(df.ordered_frac.values >= 0.9)[0]
    if len(iord):
        t_ord = float(tm[iord[0]]); win = (tm >= t_ord - 1e-9) & (tm <= t_ord + 5 + 1e-9)
        jf = np.where(win)[0][np.argmax(seg[win])]; seg_form = int(seg[jf]); t_form = float(tm[jf])
        j25 = np.searchsorted(tm, t_ord + 25.0 - 1e-9)
        persist = float(seg[j25] / seg_form) if (seg_form > 0 and j25 < len(seg)) else (0.0 if seg_form == 0 else float('nan'))
        late = tm >= tm[-1] - 10 - 1e-9
        tangle = float(seg[late].mean() / seg_form) if seg_form > 0 else float('nan')
        comp_f = df.iloc[jf]; long_frac_form = float(comp_f.get('long_frac', float('nan')))
    else:
        t_ord = t_form = float('nan'); seg_form = 0; persist = 0.0; tangle = float('nan'); long_frac_form = float('nan')
    res = dict(arm=arm, **p, N=N, box_com=box, box_over_Lc=box / Lc, dx_com=dx, dx_phys_switch=dx * info['a_switch'],
               dx_phys_end=dx * float(df.a.iloc[-1]), dt=dt, linear=info, k_peak_lattice=kpk,
               t_ordered=t_ord, t_form=t_form, seg_form=seg_form, Lambda_form=lam(seg_form), persist25=persist,
               long_frac_form=long_frac_form, tangle_growth_late=tangle,
               seg_max=int(seg.max()), Lambda_max=lam(int(seg.max())), bianchi_max=int(df.bianchi.max()),
               bianchi_frac=float(df.bianchi.max() / max(int(seg.max()), 1)), seg_end=int(seg[-1]),
               order_end=float(df.order.iloc[-1]), ordered_frac_end=float(df.ordered_frac.iloc[-1]),
               w_late=float(df.p.iloc[-20:].sum() / df.rho.iloc[-20:].sum()), a_end=float(df.a.iloc[-1]),
               runtime_s=time.perf_counter() - t0w, backend='jax' if HAVE_JAX else 'numpy')
    (od / 'arm_summary.json').write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != 'linear'}, indent=1, default=float), flush=True)
    return res


def verdict(R):
    """Gate at the FORMATION epoch: t_ord = first sample with >= 90 % ordered cells; seg_form = max raw segments
    in [t_ord, t_ord + 5].  Lambda_form = seg_form * dx * L_c^2 / box^3 (string length per L_c^2, comoving)."""
    need = ['J025', 'J100', 'C0']
    if not all(k in R for k in need): return None
    c0, j, j1 = R['C0'], R['J025'], R['J100']
    if c0['seg_max'] >= 20: v = 'C32B_CONTROL_FAIL'
    elif (j['Lambda_form'] >= 0.3 and j['persist25'] >= 0.3 and j['bianchi_frac'] <= 0.01): v = 'C32B_SCARS_FORM'
    elif j['Lambda_form'] < 0.03 and j1['Lambda_form'] < 0.03: v = 'C32B_NO_SCARS'
    else: v = 'C32B_MIXED'
    sec = {}
    if 'E100' in R: sec['C32_NULL_WAS_BOX'] = bool(R['E100']['Lambda_form'] >= 0.3)
    if 'J025e6' in R:
        a, b = j['Lambda_form'], R['J025e6']['Lambda_form']
        sec['KZ_WEAK'] = bool(a > 0 and b > 0 and max(a, b) / min(a, b) <= 3.0)
    if 'J100mu01' in R: sec['MU01_FORMS'] = bool(R['J100mu01']['Lambda_form'] >= 0.3)
    if 'J025hr' in R:
        a, b = j['Lambda_form'], R['J025hr']['Lambda_form']
        sec['RESOLUTION_OK'] = bool(a > 0 and b > 0 and max(a, b) / min(a, b) <= 2.0)
        sec['TANGLE_RESOLVED'] = bool(np.isfinite(j['tangle_growth_late']) and np.isfinite(R['J025hr']['tangle_growth_late'])
                                      and max(j['tangle_growth_late'], R['J025hr']['tangle_growth_late']) /
                                      max(min(j['tangle_growth_late'], R['J025hr']['tangle_growth_late']), 1e-9) <= 2.0)
    sec['TANGLE_GROWS_J025'] = bool(np.isfinite(j['tangle_growth_late']) and j['tangle_growth_late'] > 2.0)
    return dict(verdict=v, secondary=sec)


def selftest():
    """Parity with the parent (only when the parent modules import, i.e. in WSL) + Hessian check."""
    rep = {}
    try:
        import v0_66_3d_autonomous as H
        import v0_72_polynomial_trigger as M
        lam = float(M.coefficient_np(0.0, 'p4_hermite')[0]); assert abs(lam - LAM) < 1e-12
        assert H.ZQ == ZQ and H.LP == LP and H.LBASE == LBASE
        rng = np.random.default_rng(1); Qr = rng.normal(size=(3, 3, 3, 4, 4)); Qr = 0.5 * (Qr + np.swapaxes(Qr, -1, -2))
        g_par = np.asarray(H.q_terms_j(H.jnp.asarray(Qr, dtype=H.jnp.float32), lam)[4], np.float64); g_me = np.asarray(grad_U(xp.asarray(Qr, dtype=xp.float32)), np.float64)
        rep['grad_parity_rel'] = float(np.max(np.abs(g_par - g_me)) / np.max(np.abs(g_par)))
        import v0_72r_topology as T
        Qt = rng.normal(size=(12, 12, 12, 4, 4)) * 0.05; Qt = 0.5 * (Qt + np.swapaxes(Qt, -1, -2))
        n = rng.normal(size=(12, 12, 12, 4)); n /= np.linalg.norm(n, axis=-1, keepdims=True)
        Qt = Qt + n[..., :, None] * n[..., None, :]
        r = T.compute_flux(Qt)
        seg_par = int(((r['fxy'] < 0) & r['mxy']).sum() + ((r['fyz'] < 0) & r['myz']).sum() + ((r['fzx'] < 0) & r['mzx']).sum())
        seg_me, bia_me, _ = flux_segments(Qt, np)
        rep['flux_parity'] = dict(parent=(seg_par, int(r['bianchi'])), mine=(int(seg_me), int(bia_me)))
        rep['parity_ok'] = bool(rep['grad_parity_rel'] < 1e-5 and seg_par == int(seg_me) and int(r['bianchi']) == int(bia_me))
    except ImportError as e:
        rep['parity_skipped'] = str(e)
    # Hessian of the traceless directions at x0 I
    # float64 numpy copy of U_density: in WSL xp is jax with x64 disabled (float32), and a finite-difference Hessian with
    # h = 1e-4 in float32 is pure round-off (fix by Claude, main conversation, 26-sep ~17:35; physics and gates unchanged)
    B = traceless_basis(); h = 1e-4
    Ub = lambda Q: float((LAM / 4.0) * (np.trace(Q) - 1.0) ** 2 + (LP / 4.0) * np.sum((Q @ Q - Q) ** 2))
    Q0 = X0 * np.eye(4); m2 = []
    for b in B:
        m2.append((Ub(Q0 + h * b) - 2 * Ub(Q0) + Ub(Q0 - h * b)) / (h * h) / ZQ)
    rep['traceless_m2_at_x0'] = [round(x, 6) for x in m2]; rep['MT2_ok'] = bool(max(abs(x - MT2) for x in m2) < 1e-4)
    print('SELFTEST', rep, flush=True)
    return rep


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--N', type=int, default=192); ap.add_argument('--out', required=True)
    ap.add_argument('--T', type=float, default=T_LAT); ap.add_argument('--arms', default='J025,C0,J100,E100,J025hr,J025e6,J100mu01')
    ap.add_argument('--linear-only', action='store_true')
    ap.add_argument('--box-lc', type=float, default=None, help='smoke tests only; the preregistered value is 3')
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.box_lc is not None: globals()['BOX_OVERRIDE'] = a.box_lc
    st = selftest()
    if st.get('parity_ok') is False or not st['MT2_ok']:
        (out / 'SELFTEST_FAILED.json').write_text(json.dumps(st, indent=2)); raise SystemExit('selftest failed')
    if a.linear_only:
        rows = {}
        for arm in a.arms.split(','):
            if arm == 'C0': continue
            _, _, _, _, info = linear_stage(arm, 4.0); info['box_com'] = ARMS[arm].get('box_lc', BOX_LC) * info['Lc_com']
            info['dx_com'] = info['box_com'] / a.N; rows[arm] = info; print(arm, json.dumps(info, default=float), flush=True)
        (out / 'C32B_LINEAR.json').write_text(json.dumps(rows, indent=2, default=float)); return
    R = dict(selftest=st); cache = {}
    for arm in a.arms.split(','):
        R[arm] = run_arm(arm, a.N, out, a.T)
        R['verdict'] = verdict(R)
        (out / 'C32B_RESULTS.json').write_text(json.dumps(R, indent=2, default=float))
    print('VERDICT', R.get('verdict'), flush=True)


if __name__ == '__main__':
    main()
