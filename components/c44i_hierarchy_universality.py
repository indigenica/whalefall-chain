#!/usr/bin/env python3
"""C44i — Independent check of the C44E 'two-sheet recursive balance' inflaton, and WHY it gives n_s ~ 0.974 (CPU, seconds).

C44E (external conversation, folder claudio_astro/WHALEFALL_C44E_TWO_SHEET_ESCAPE):
   V_L(chi) = V0/Z_L sum_{n<L} q^n tanh^2( chi / (f s^n) ),  f = 1/2, q = 1/2, s = 2  ->  n_s ~ 0.9737, r ~ 0.0048 at N = 51.8.
Claim tested here (written before running): the result is NOT a coincidence of (1/2, 2, 1/2) but a structural law.
For f s^k << chi << f s^(k+1) the levels n <= k are saturated and the missing plateau energy is
   1 - V/V0 ~ sum_{n>k} q^n ~ q^k ~ (chi/f)^(-p),   p = ln(1/q) / ln(s),
i.e. the hierarchy is an inverse-power plateau V ~ V0 (1 - (mu/chi)^p), whose slow-roll attractor is
   n_s = 1 - 2 (p+1) / ((p+2) N),    r = 16 eps(N) -> depends on mu (~ f).
q s = 1 (energy halves when the scale doubles: a scale-invariant, 'marginal' hierarchy) gives p = 1 and n_s = 1 - 4/(3N) = 0.97426 at
N = 51.8, independently of the particular q and of f (f only sets r).  Labels (frozen here):
  C44I_REPRODUCED      : our L = 6 and L = 10 numbers agree with C44E within 1e-4 (n_s) and 2 % (r).
  C44I_UNIVERSAL_QS1   : for q s = 1 with q in {0.3, 0.5, 0.7} and f in {0.25, 0.5, 1}, L large, |n_s - (1 - 4/(3N))| < 0.002 for all.
  C44I_P_LAW           : for q s != 1 (p = 0.5, 2, 3) |n_s - formula(p)| < 0.003.
"""
import json, math, sys
from pathlib import Path
import numpy as np

NSTAR = 51.8


def model(q, s, f, L):
    w = q ** np.arange(L); Z = w.sum(); sc = f * s ** np.arange(L)
    def V(x):
        return (w[:, None] * np.tanh(x[None, :] / sc[:, None]) ** 2).sum(0) / Z
    def dV(x):
        t = np.tanh(x[None, :] / sc[:, None]); sech2 = 1 - t * t
        return (w[:, None] * 2 * t * sech2 / sc[:, None]).sum(0) / Z
    def d2V(x):
        t = np.tanh(x[None, :] / sc[:, None]); sech2 = 1 - t * t
        return (w[:, None] * 2 * sech2 * (sech2 - 2 * t * t) / sc[:, None] ** 2).sum(0) / Z
    return V, dV, d2V


def predict(q, s, f, L, N=NSTAR):
    """Adaptive: chi_end from eps = 1 (brentq), N(chi) = int V/V' (quad), chi_* from N(chi_*) = N (brentq)."""
    from scipy.optimize import brentq
    from scipy.integrate import quad
    V, dV, d2V = model(q, s, f, L)
    v = lambda x: float(V(np.array([x]))[0]); dv = lambda x: float(dV(np.array([x]))[0]); d2v = lambda x: float(d2V(np.array([x]))[0])
    eps = lambda x: 0.5 * (dv(x) / v(x)) ** 2
    xs = np.geomspace(1e-3, 1e7, 4000); e = [eps(x) for x in xs]
    i = next(k for k in range(len(xs) - 1) if e[k] >= 1 and e[k + 1] < 1)
    xend = brentq(lambda x: eps(x) - 1, xs[i], xs[i + 1])
    def Nof(x):
        pts = [xend * (x / xend) ** (k / 40) for k in range(41)]
        return sum(quad(lambda y: v(y) / dv(y), pts[k], pts[k + 1], limit=200)[0] for k in range(40))
    hi = xend * 2
    while Nof(hi) < N: hi *= 2
    xs_ = brentq(lambda x: Nof(x) - N, xend, hi, xtol=1e-10)
    return dict(n_s=1 - 6 * eps(xs_) + 2 * d2v(xs_) / v(xs_), r=16 * eps(xs_), chi_star=xs_, chi_end=xend, V_star_over_V0=v(xs_))


def ns_formula(p, N=NSTAR): return 1 - 2 * (p + 1) / ((p + 2) * N)


def main(root):
    out = Path(root) / 'c44i_hierarchy_universality'; out.mkdir(parents=True, exist_ok=True)
    res = dict(N=NSTAR)
    ref = {6: (0.97363390, 0.00480837), 10: (0.97370295, 0.00482692)}
    rep = {L: predict(0.5, 2.0, 0.5, L) for L in (1, 2, 4, 6, 10, 16)}
    res['reproduction'] = rep
    res['C44I_REPRODUCED'] = bool(all(abs(rep[L]['n_s'] - ref[L][0]) < 1e-4 and abs(rep[L]['r'] / ref[L][1] - 1) < 0.02 for L in ref))
    uni = []
    for q in (0.3, 0.5, 0.7):
        s = 1 / q
        L = int(math.ceil(math.log(1e4) / math.log(s))) + 2        # levels up to ~1e4 f
        for f in (0.25, 0.5, 1.0):
            p_ = predict(q, s, f, L); p_.update(q=q, s=s, f=f, L=L); uni.append(p_)
    res['qs1_family'] = uni; res['ns_formula_p1'] = ns_formula(1.0)
    res['C44I_UNIVERSAL_QS1'] = bool(all(abs(u['n_s'] - ns_formula(1.0)) < 0.002 for u in uni))
    plaw = []
    for q, s in ((0.5, 4.0), (0.25, 2.0), (0.125, 2.0)):          # p = 0.5, 2, 3
        p = math.log(1 / q) / math.log(s); L = int(math.ceil(math.log(1e4) / math.log(s))) + 2
        pr = predict(q, s, 0.5, L); pr.update(q=q, s=s, p=p, formula=ns_formula(p)); plaw.append(pr)
    res['p_law'] = plaw
    res['C44I_P_LAW'] = bool(all(abs(x['n_s'] - x['formula']) < 0.003 for x in plaw))
    # r vs f for q s = 1 (mu ~ f): analytic p = 1 attractor r = 8 mu^(2/3) / (3N)^(4/3) -> fit mu/f
    res['r_vs_f_qs1'] = [dict(f=u['f'], q=u['q'], r=u['r'], mu_eff=(u['r'] * (3 * NSTAR) ** (4 / 3) / 8) ** 1.5) for u in uni]
    (out / 'C44I_RESULTS.json').write_text(json.dumps(res, indent=2, default=float))
    for L, v in rep.items(): print(f"L={L:2d} n_s={v['n_s']:.6f} r={v['r']:.6f} chi*={v['chi_star']:.3f}")
    print('formula p=1:', round(ns_formula(1.0), 5))
    for u in uni: print(f"q={u['q']} s={u['s']:.3f} f={u['f']} L={u['L']} n_s={u['n_s']:.5f} r={u['r']:.5f}")
    for x in plaw: print(f"p={x['p']:.2f} n_s={x['n_s']:.5f} formula={x['formula']:.5f} r={x['r']:.5f}")
    print({k: res[k] for k in ('C44I_REPRODUCED', 'C44I_UNIVERSAL_QS1', 'C44I_P_LAW')})


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'claude/results')
