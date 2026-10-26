"""Pointwise self-convergence between runs at N and 2N on the common nodes."""
import numpy as np


def compare(lo, hi):
    """lo, hi: dicts/npz with R, S, P, M. Returns |differences| on the coarse grid."""
    out = {}
    for k in ("R", "S", "P", "M"):
        a = lo[k]; b = hi[k][::2, ::2]
        out[k] = np.abs(a - b)
    return out


def trusted_mask(lo, hi, tolP=1e-2, tolS=5e-2, tolR=1e-2):
    d = compare(lo, hi)
    R = lo["R"]
    ok = (d["P"] < tolP) & (d["S"] < tolS) & (d["R"] < tolR * np.maximum(1.0, np.abs(R)))
    # trusted = ok and all causal past nodes ok (cumulative AND along u and v)
    t = ok.copy()
    n = t.shape[0]
    for i in range(n):
        for j in range(i, n):
            if i > 0 and not t[i - 1, j]:
                t[i, j] = False
            if j > i and not t[i, j - 1]:
                t[i, j] = False
    return t, d
