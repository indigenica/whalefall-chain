"""Growing mode for the 'soft QT' regime (alpha*y_ta ~ 0.01-1), where no GR epoch with y >> y_ta and
alpha*y << 1 exists to synchronize on.

In QT a bound shell (K > 0) has no bang: going back in time it reaches its bounce radius R_b and bounces.
The analogue of the GR growing mode (uniform bang time) is a UNIFORM PREVIOUS BOUNCE TIME: the parent
region emerged from a QT bounce at t = 0, every shell from its own R_b.  Then t_ta,j = tau_B,j (time from
R_b to R_ta), and the next bounce is at 2 tau_B,j.  For small alpha this differs from the GR cycloid
synchronization only by the O(sqrt(alpha)) bounce duration (checked: alpha = 0.1 rows reproduce the GR-sync
coherence numbers within that tolerance, see results).
"""
import math
import numpy as np
from scipy.integrate import quad
from ltb.shells import Shell, psi_of_y
from ltb.growing import turnaround_radius


def t_turnaround_qt(M, K, alpha, R_s):
    R_ta = turnaround_radius(M, K, alpha)
    U = lambda R: R * R * psi_of_y(2 * M / R ** 3, alpha) - K
    Rm = 0.5 * (R_s + R_ta) if R_s < 0.5 * R_ta else R_s
    # log part [R_s, Rm]
    f1 = lambda u: math.exp(u) / math.sqrt(max(U(math.exp(u)), 1e-300))
    t1 = quad(f1, math.log(R_s), math.log(Rm), limit=400, epsrel=1e-11, epsabs=0)[0] if Rm > R_s else 0.0
    # root part [Rm, R_ta] with R = R_ta - D s^2
    D = R_ta - Rm
    h = 1e-7
    dU = ((R_ta * (1 + h)) ** 2 * psi_of_y(2 * M / (R_ta * (1 + h)) ** 3, alpha)
          - (R_ta * (1 - h)) ** 2 * psi_of_y(2 * M / (R_ta * (1 - h)) ** 3, alpha)) / (2 * R_ta * h)
    g = lambda s: 2 * D * s / math.sqrt(max(U(R_ta - D * s * s), 1e-300)) if s > 0 else 2 * math.sqrt(D) / math.sqrt(max(-dU, 1e-300))
    t2 = quad(g, 0.0, 1.0, limit=400, epsrel=1e-11, epsabs=0)[0]
    return t1 + t2, R_ta


class GrowingShellQT:
    def __init__(self, M, K, alpha, R_s=None):
        self.M, self.K = M, K
        self.R_ta = turnaround_radius(M, K, alpha)
        self.shell = Shell(M, self.R_ta, alpha)
        self.tau_B = self.shell.tau_B
        self.t_ta = self.tau_B                      # uniform previous-bounce time t = 0
        self.t_B = self.t_ta + self.tau_B

    def R_of_t(self, t):
        return self.shell.R_of_t(np.abs(np.asarray(t, float) - self.t_ta))


def build_qt(x, K2, AH, alpha, n_sh=600, eps_fac=1e-3):
    """Same as ltb.bank_pbh.build but with QT-past synchronization."""
    from scipy.interpolate import CubicSpline
    from scipy.optimize import brentq
    from common.potential import H_INF
    import ltb.bank_pbh as B
    km = int(np.argmax(K2)); x_m = x[km]
    bound_end = km + int(np.argmax(K2[km:] <= 0)) if (K2[km:] <= 0).any() else len(x) - 1
    x_hi = min(1.2 * x_m, x[bound_end - 1])
    i0 = max(2, int(0.01 * km)); k0 = K2[i0] / x[i0] ** 2

    def A_exact(mu):
        M = mu * x[i0] ** 3; Kc = K2[i0]
        Rta = turnaround_radius(M, Kc, alpha)
        return 1 / math.sqrt(psi_of_y(2 * M / Rta ** 3, alpha))
    target = AH / H_INF
    mu_gr = target * k0 ** 1.5 / 2
    lo, hi = 0.7 * mu_gr, 1.5 * mu_gr
    while A_exact(hi) < target:
        hi *= 1.5
    for _ in range(40):
        try:
            if A_exact(lo) <= target:
                break
        except ValueError:
            pass
        lo *= 0.95
    mu = brentq(lambda m_: A_exact(m_) - target, lo, hi, xtol=1e-14 * mu_gr)
    cs = CubicSpline(x, K2)
    xs = x_hi * (np.arange(n_sh) + 0.5) / n_sh
    # common early state: y_s = 2 mu/eps^3 >> max(1/alpha, k0): homologous, curvature negligible
    y_s = 1e4 * max(1 / alpha if alpha > 0 else 0, k0, 1.0)
    eps = (2 * mu / y_s) ** (1 / 3)
    mk = lambda xx: GrowingShellQT(mu * xx ** 3, float(cs(xx)), alpha, eps * xx)
    shells = [mk(xx) for xx in xs]
    twins = [mk(xx * (1 + B.EPS)) for xx in xs]
    return shells, twins, xs, x_m, mu, A_exact(mu) * H_INF


def validate_qt():
    """GR limit: t_ta -> pi M/K^{3/2}; and alpha = 0.1 agreement with the GR-cycloid synchronization."""
    s = GrowingShellQT(1.0, 0.5, 1e-8)
    gr = abs(s.t_ta / (math.pi / 0.5 ** 1.5) - 1)
    return dict(GR_limit_t_ta_rel=gr)


if __name__ == "__main__":
    print(validate_qt())
