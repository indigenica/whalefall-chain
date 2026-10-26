"""Initial-cone profiles F(r) for the double-null runs.

gaussian(A, w):  F = A exp(-(r/w)^2)
chistar(s, alpha): F = sqrt2 * s * phi_star(r), the static complex-proxy ground
  state of C59a (phi0 = 0.12, ballena potential), tail cut at its first minimum
  and smoothly tapered over 2 length units (C59a-2 used a hard cut on a t = 0 slice;
  on a null cone a hard cut would put a delta in phi_v).  Same family as C59a-2.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from common.bosonstar import ground_state
from common.potential import V, dV


def gaussian(A, w):
    F = lambda r: A * np.exp(-(np.asarray(r) / w) ** 2)
    dF = lambda r: -2 * np.asarray(r) / w ** 2 * F(r)
    return F, dF


_cache = {}


def star_shape(alpha, phi0=0.12):
    key = (alpha, phi0)
    if key not in _cache:
        w, sol = ground_state(phi0, alpha, V, dV, r_max=50.0)
        r = np.linspace(0, sol.t[-1], 4001)
        rr = np.clip(r, sol.t[0], sol.t[-1])
        prof = sol.sol(rr)[0]
        prof[r < sol.t[0]] = phi0
        k = int(np.argmin(np.abs(prof) + (r < 1.0) * 9))
        rc = r[k]
        taper = np.where(r < rc - 2, 1.0, np.where(r < rc, 0.5 * (1 + np.cos(np.pi * (r - rc + 2) / 2)), 0.0))
        prof = prof * taper
        _cache[key] = (r, prof, float(rc), float(w))
    return _cache[key]


def chistar(s, alpha, phi0=0.12):
    r, prof, rc, w = star_shape(alpha, phi0)
    cs = CubicSpline(r, np.sqrt(2) * s * prof, bc_type=((1, 0.0), (1, 0.0)))
    F = lambda x: np.where(np.asarray(x) <= r[-1], cs(np.minimum(x, r[-1])), 0.0)
    dF = lambda x: np.where(np.asarray(x) <= r[-1], cs(np.minimum(x, r[-1]), 1), 0.0)
    return F, dF, rc
