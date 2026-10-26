"""Coherent-region measures with intrusion of outer shells (amendment after audit_pbh/AUDIT_PBH.md).

For a region [0..k] of inner shells (shells uniform in x, masses w_j):
  internal validity: no adjacent crossing among 0..k, no local caustic (twin shells) at 0..k,
                     all shells bounced (t >= max t_B) and still expanding (t < min(t_ta + 2 tau_B));
  intrusion:  m_in(k, t) = mass of outer shells j > k with R_j(t) < R_k(t)
      strict    -> region ends at the first t with m_in > 0
      tolerant  -> region ends when m_in > f * M_k, f in FRACS (the intruders are then a minor perturbation
                   of the region's own mass; f = 0.1 is the reference value, 0.01 and 0.3 sensitivity).
R_coh(criterion) = max_k R_k(T_k) over regions with B_k < T_k.  The strict criterion also kills the
homogeneous control (a top-hat's edge layer falls onto the re-expanding core), so the tolerant one is the
discriminating measure; both are reported.
"""
import math
import numpy as np
from common.potential import H_INF

FRACS = (0.01, 0.1, 0.3)


def analyse_intrusion(shells, twins, xs, x_m, alpha, dt=None, t_max_factor=1.0):
    n = len(shells)
    tB = np.array([s.t_B for s in shells]); tta = np.array([s.t_ta for s in shells])
    tauB = np.array([s.tau_B for s in shells])
    dx = xs[1] - xs[0]
    w = (xs + dx / 2) ** 3 - np.maximum(xs - dx / 2, 0) ** 3          # mass weights (mu factored out)
    Wlt = np.concatenate([[0.0], np.cumsum(w)[:-1]])                   # mass of shells j < k
    Mk = Wlt + w
    core = xs <= x_m
    T0 = tta.min() * 0.999
    T_end = float((tta[core] + 2 * tauB[core]).max()) * t_max_factor
    dt = dt or math.sqrt(alpha) / 100
    n_total = int(math.ceil((T_end - T0) / dt)) + 1
    cpair = np.full(n, np.inf); cl = np.full(n, np.inf); Istr = np.full(n, np.inf)
    Itol = {f: np.full(n, np.inf) for f in FRACS}
    chunk = max(1, int(4e6 // n))
    tend_all = tta + 2 * tauB
    cumB_all = np.maximum.accumulate(tB)
    for a in range(0, n_total, chunk):
        if a > 0:
            # exact early termination: Tint_k (internal validity limit) can only decrease as more crossings are
            # found, and T_k <= Tint_k; once the scan has passed max Tint_k over the regions that can still be
            # valid (B_k < Tint_k), nothing else can change any R_coh.
            t_now = T0 + dt * (a - 1)
            cp = np.concatenate([[np.inf], cpair[:-1]])
            Tint = np.minimum.accumulate(np.minimum(np.minimum(cl, tend_all), cp))
            cand = cumB_all < Tint
            if not cand.any() or t_now >= float(Tint[cand].max()):
                break
        tt = T0 + dt * np.arange(a, min(a + chunk, n_total))
        Rk = np.vstack([s.R_of_t(tt) for s in shells])
        Tw = np.vstack([s.R_of_t(tt) for s in twins])
        nt = len(tt)

        def first(mask, arr):
            hit = mask.any(axis=1) & np.isinf(arr)
            idx = np.nonzero(hit)[0]
            arr[idx] = tt[np.argmax(mask[idx], axis=1)]
        pm = np.zeros((n, nt), bool); pm[:-1] = Rk[1:] <= Rk[:-1]       # pair (k, k+1) crossing, stored at k
        first(pm, cpair)
        first(Tw <= Rk, cl)
        suf = np.minimum.accumulate(Rk[::-1], axis=0)[::-1]
        sm = np.zeros((n, nt), bool); sm[:-1] = suf[1:] <= Rk[:-1]
        first(sm, Istr)
        order = np.argsort(Rk, axis=0)
        cws = np.cumsum(w[order], axis=0)
        inv = np.empty_like(order)
        np.put_along_axis(inv, order, np.broadcast_to(np.arange(n)[:, None], order.shape), axis=0)
        Mbelow = np.take_along_axis(cws, inv, axis=0)                  # mass with R <= R_k (incl. k)
        m_in = Mbelow - w[:, None] - Wlt[:, None]                      # outer mass inside R_k
        for f in FRACS:
            first(m_in > f * Mk[:, None] * (1 + 1e-12), Itol[f])
    out = {}
    tend = tta + 2 * tauB
    for name, I in [("strict", Istr)] + [(f"tol{f}", Itol[f]) for f in FRACS]:
        best = dict(R=0.0, k=-1)
        Bk = -np.inf; Tint = np.inf
        for k in range(n):
            Bk = max(Bk, tB[k])
            Tint = min(Tint, cl[k], tend[k], cpair[k - 1] if k > 0 else np.inf)
            Tk = min(Tint, I[k])
            if Tint <= Bk:
                break
            if Tk > Bk:
                Rk_ = float(shells[k].R_of_t(Tk - 1e-9))
                if Rk_ > best["R"]:
                    best = dict(R=Rk_, k=k, t=float(Tk))
        out[name] = dict(R_coh_H=best["R"] * H_INF, x_coh_over_xm=float(xs[best["k"]] / x_m) if best["k"] >= 0 else 0.0,
                         M_coh_over_Mm=float(Mk[best["k"]] / Mk[np.searchsorted(xs, x_m) - 1]) if best["k"] >= 0 else 0.0)
    jm = int(np.argmin(np.abs(xs - x_m)))
    out["ceiling_H"] = float(shells[jm].R_ta * H_INF)
    out["n_t"] = n_total
    return out
