import numpy as np
from dn.solver import Grid, initial_ray, evolve

ZERO = lambda x: 0.0 * np.asarray(x)


def flat_err(N, L=20.0, A=1e-3):
    h = L / N
    F = lambda r: A * np.exp(-((np.asarray(r) - 3.0) / 1.0) ** 2)
    dF = lambda r: -2 * (np.asarray(r) - 3.0) * F(r)
    g = Grid(N, h)
    initial_ray(g, F, dF, ZERO, ZERO, alpha=0.0, G=0.0)
    st = evolve(g, ZERO, ZERO, alpha=0.0, G=0.0)
    u = h * np.arange(N + 1)
    U, Vv = np.meshgrid(u, u, indexing="ij")
    with np.errstate(invalid="ignore", divide="ignore"):
        ex = (Vv * F(Vv / 2) - U * F(U / 2)) / (Vv - U)
    mask = (Vv > U + 1e-12) & np.isfinite(g.p)
    err = np.nanmax(np.abs(g.p - ex)[mask]) / A
    rerr = np.nanmax(np.abs(g.r - (Vv - U) / 2)[mask])
    return st, err, rerr


def test_flat_exact():
    """Amended after the independent audit: in flat space the diamond scheme is EXACT for
    the radial wave equation (r phi)_uv = 0 with r = (v-u)/2, so the error is round-off
    and an 'order' is meaningless (the original order test failed for that reason).
    Second-order convergence is tested in curved space by audit/check_convergence.py."""
    for N in (400, 800):
        st, e, re = flat_err(N)
        assert st == "COMPLETE" and e < 1e-9 and re < 1e-12, (N, e, re)
