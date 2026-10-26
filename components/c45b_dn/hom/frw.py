"""C45b-0 - homogeneous closed QT patch: dust (optional) + chi, 8 pi G = 1.

Port of the audited C45 v2 integrator (WHALEFALL/claude/tools/c45_fold_launches_chi.py),
generalised to a chi-only patch (rho_dust = 0) whose turnaround state is set by
the C59a chi-star density.  Same equations, same declared constraint
relaxation, same potential:
   3 h(psi) = rho_d + chidot^2/2 + V0 Vn(chi),  psi = H^2 + 1/a^2
   Hdot = 1/a^2 - (rho_d + chidot^2)(1 - z^2)^{3/2}/2,   z = alpha psi
   chiddot = -3 H chidot - V0 Vn'(chi)
"""
import math
import numpy as np
from scipy.integrate import solve_ivp, quad

Lh, qh, sh, fh = 10, 0.5, 2.0, 0.5
W = qh ** np.arange(Lh); Z = W.sum(); SC = fh * sh ** np.arange(Lh)
V0_BENCH = 1.67e-10


def Vn(x):
    return float((W * np.tanh(x / SC) ** 2).sum() / Z)


def dVn(x):
    t = np.tanh(x / SC); return float((W * 2 * t * (1 - t * t) / SC).sum() / Z)


def _chi_end():
    xs = np.geomspace(1e-2, 1e6, 3000); eps = [0.5 * (dVn(v) / Vn(v)) ** 2 for v in xs]
    return float(xs[next(k for k in range(len(xs) - 1) if eps[k] >= 1 and eps[k + 1] < 1)])


CHI_END = _chi_end()


def N_remaining(chi):
    x = abs(chi); xs = np.geomspace(1e-2, 1e6, 3000)
    eps = [0.5 * (dVn(v) / Vn(v)) ** 2 for v in xs]
    i = next((k for k in range(len(xs) - 1) if eps[k] >= 1 and eps[k + 1] < 1), None)
    if i is None or x <= xs[i]:
        return 0.0
    xe = xs[i]
    if dVn(x) <= 0 or x > 20 * SC[-1]:
        return float('inf')
    pts = np.geomspace(xe, x, 60)
    return float(sum(quad(lambda y: Vn(y) / max(dVn(y), 1e-300), pts[k], pts[k + 1], limit=200)[0] for k in range(59)))


def run(a0, al, V0, rho_d0, chi0, chid0, max_cycles=60, rtol=1e-11, wall=300.0, efolds_stop=60.0, t_extra=None,
        return_traj=False):
    """Start at turnaround-like state a = a0, H = 0 (constraint fixes the dust
    constant so that 3 h(1/a0^2) = E exactly).  rho_d0 = initial dust density
    requested; the chi part is (chi0, chid0).  Returns summary dict."""
    import time
    h = lambda p: p / math.sqrt(1 - al * al * p * p)
    E_chi = 0.5 * chid0 ** 2 + V0 * Vn(chi0)
    # exact constraint at H = 0: dust absorbs the residual (as C45) -> Aeff
    Aeff = (3 * h(1 / a0 ** 2) - E_chi) * a0 ** 3 / 3
    if Aeff < -1e-12 * abs(E_chi) * a0 ** 3:
        raise ValueError("chi energy exceeds the closed-patch budget at a0 (choose a0 from the constraint)")
    Aeff = max(Aeff, 0.0)
    t0w = time.time()

    def rhs(t, y):
        if time.time() - t0w > wall:
            raise TimeoutError
        a, H, x, xd = y
        z = al * (H * H + 1 / (a * a)); z = min(z, 1 - 1e-14)
        Hd = 1 / (a * a) - (3 * Aeff / a ** 3 + xd * xd) * (1 - z * z) ** 1.5 / 2
        E = 3 * Aeff / a ** 3 + 0.5 * xd * xd + V0 * Vn(x); yv = E / 3
        psic = yv / math.sqrt(1 + al * al * yv * yv); H2c = psic - 1 / (a * a)
        if H2c > 0.05 * psic:
            Hc = math.copysign(math.sqrt(H2c), H); Hd -= 10 * abs(Hc) * (H - Hc); H = Hc
        return [a * H, Hd, xd, -3 * H * xd - V0 * dVn(x)]
    ev_small = lambda t, y: y[0] - 1e-3 * math.sqrt(al); ev_small.terminal = True
    ev_big = lambda t, y: math.log(y[0] / a0) - efolds_stop; ev_big.terminal = True
    # chi-only bank: room for the e-folds themselves (t_extra=None -> (efolds_stop + 10)/H_inf);
    # the C45 reference keeps C45's window exactly (t_extra = 0)
    T = max_cycles * 4.0 * a0 + ((efolds_stop + 10) / math.sqrt(V0 / 3) if t_extra is None else t_extra)
    try:
        s = solve_ivp(rhs, (0, T), [a0, 0.0, chi0, chid0], method='LSODA', rtol=rtol, atol=1e-16,
                      events=[ev_small, ev_big])
    except TimeoutError:
        return dict(outcome='timeout')
    a, H, x, xd = s.y
    V = V0 * np.array([Vn(v) for v in x]); rhod = 3 * Aeff / a ** 3
    psi = H * H + 1 / a ** 2
    E = rhod + 0.5 * xd * xd + V
    zz = np.minimum(al * psi, 1 - 1e-14)
    cons = np.abs(3 * psi / np.sqrt(1 - zz ** 2) - E) / E
    Hd = 1 / a ** 2 - (rhod + xd * xd) * (1 - zz * zz) ** 1.5 / 2
    infl = (H > 0) & (Hd + H * H > 0) & (V > xd * xd) & (np.abs(x) > CHI_END)
    N_inf = float(np.sum(np.diff(np.log(a))[infl[:-1]]))
    nb = int(np.sum((H[:-1] < 0) & (H[1:] >= 0)))
    out = dict(outcome=('launched60' if len(s.t_events[1]) else ('crunch' if len(s.t_events[0]) else 'cycling')),
               n_bounces=nb, a_min=float(a.min()), chi_abs_max=float(np.max(np.abs(x))), chi_end=float(x[-1]),
               N_inf=N_inf, constraint_max=float(np.nanmax(cons)), Aeff=float(Aeff),
               dust_fraction_start=float(3 * Aeff / a0 ** 3 / E[0]), z_max=float(zz.max()), t_end=float(s.t[-1]))
    if out['outcome'] == 'launched60':
        out['N_remaining_attractor'] = N_remaining(x[-1])
    out['fertile'] = bool(out['outcome'] == 'launched60' and out['chi_abs_max'] > 5 and out['constraint_max'] < 1e-3)
    out['efolds_total'] = float(np.log(a[-1] / a0))
    out['direct_inflation'] = bool(nb == 0 and out['outcome'] == 'launched60')
    if return_traj:
        out['traj'] = dict(t=s.t, a=a, H=H, x=x, xd=xd)
    return out


def c45_reference(A, al, lam, sigma, sign):
    """Exactly the C45 initial state: a0 = A, chi = 0, chid = sign sigma sqrt(3A/A^3)."""
    V0 = lam * 4.0 / (9.0 * A * A)
    cd0 = sign * sigma * math.sqrt(3 * A / A ** 3)
    return run(A, al, V0, None, 0.0, cd0, max_cycles=200, t_extra=0.0)
