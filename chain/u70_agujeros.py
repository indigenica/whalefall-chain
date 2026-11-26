"""U70 (docs/preregistros/ESPEC_U70.md): O(4)-symmetric Euclidean wormholes of Q (charge of a U(1) in SO(4), radial valley of the
10-component potential of the common action) in Einstein gravity and in the NPQT realisation 3h(psi) = rho, h = psi/sqrt(1-alpha^2 psi^2),
alpha = l_QT^2.  Reduced Planck units (8 pi G = 1, lengths in l_P).  Written independently of the external napkin code.

Model (per unit 2 pi^2):
  kinetic (f_a^2/2)[g(r) r'^2 + r^2 vartheta'^2],  g = 1 + (ds/dr)^2 + 4 (du/dr)^2 along the valley min_{s,u} U(s, r, u);
  V = M^4 U_v(r), M = f_a / sqrt(Z_Q/2), f_a = F/2;  charge q = n/(2 pi^2), flux energy Phi = q^2/(2 f_a^2 r^2 a^6).
  rho_E = -K r'^2/2 + V + Phi,  psi_E = (1 - a'^2)/a^2 = h^{-1}(rho_E/3),
  a''  = -a psi + a (2 Phi - K r'^2) / (2 h'(psi)),   K r'' + K_r r'^2/2 + 3 (a'/a) K r' = V_r - 2 Phi / r.
Action of the half wormhole (one asymptotic region), reduced first-order Lagrangian L = -a^3 F(-a'^2/a^2, 1/a^2) + a^3 (K r'^2/2 + V + Phi),
  F(y, k) = 3 h(k) - 3 int_0^1 dt t^-2 [h(k + y t^2) - h(k)],  minus flat space with the same outer boundary (pure GR boundary term):
  S/(2 pi^2) = 3 a0^2 + int dtau [ -3 a (1 - a')^2 - a^3 dF + a^3 (K r'^2/2 + V + Phi) ],  dF = F - F_GR.
    python3 -B chain/u70_agujeros.py [--quick]  ->  results/u70/U70_RESULT.json"""
import json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "u70"
MPL_GEV, F0_GEV = 2.435e18, 1.0691347310310972e15
LAM, LP, ZQ = 320.0, 192.0, 400.0
EDGES = {"upper": (0.5, 1.3e-17), "lower": (0.0625, 9.2e-14)}          # (F/F0, m_G in eV), the band edges of the preprint
F_ROWS = (0.5, 0.2, 0.0625)
LQT_ROWS = (0.0, 300.0, 1000.0, 2000.0)
N_GRID = np.unique(np.round(np.geomspace(1.0, 1e6, 40)))
GL_T, GL_W = np.polynomial.legendre.leggauss(64); GL_T = 0.5 * (GL_T + 1); GL_W = 0.5 * GL_W


# ---------------------------------------------------------------- band and the required action
def m_of_F(ff):
    (fu, mu), (fl, ml) = EDGES["upper"], EDGES["lower"]
    s = math.log(ml / mu) / math.log(fu / fl)
    return mu * (fu / ff) ** s


def S_required(ff):
    F = ff * F0_GEV; m = m_of_F(ff) * 1e-9
    eps = 0.5 * m * m * F * F
    return math.log(MPL_GEV ** 4 / eps)


# ---------------------------------------------------------------- the radial valley of Q
def U_parts(s, r, u):
    lp, lm = 0.5 * (s + r), 0.5 * (s - r)
    t = s + 2 * u - 1
    U = 0.25 * LAM * t * t + 0.25 * LP * (lp * lp * (lp - 1) ** 2 + lm * lm * (lm - 1) ** 2 + 2 * u * u * (u - 1) ** 2)
    dlp = 0.5 * LP * lp * (lp - 1) * (2 * lp - 1); dlm = 0.5 * LP * lm * (lm - 1) * (2 * lm - 1)
    return U, 0.5 * LAM * t + 0.5 * (dlp + dlm), 0.5 * (dlp - dlm), LAM * t + LP * u * (u - 1) * (2 * u - 1)


def build_valley(xmax=16.0, nx=1601):
    xs = np.linspace(0.0, xmax, nx)
    S = np.empty(nx); Uu = np.empty(nx); Uv = np.empty(nx); dUr = np.empty(nx)
    guess = np.array([1.0, 0.0])
    for i, x in enumerate(xs):
        r = math.exp(x); sc = 1.0 / max(r ** 4, 1.0)
        fun = lambda z: U_parts(z[0], r, z[1])[0] * sc
        jac = lambda z: np.array([U_parts(z[0], r, z[1])[1], U_parts(z[0], r, z[1])[3]]) * sc
        res = minimize(fun, guess, jac=jac, method="BFGS", options=dict(gtol=1e-14, maxiter=2000))
        guess = res.x
        S[i], Uu[i] = res.x
        U, _, dr, _ = U_parts(res.x[0], r, res.x[1])
        Uv[i], dUr[i] = U, dr
    Uv[0] = 0.0; dUr[0] = 0.0
    return dict(xs=xs, s=CubicSpline(xs, S), u=CubicSpline(xs, Uu), lnU=None, Uv=CubicSpline(xs, Uv), dUr=CubicSpline(xs, dUr),
                s_raw=S, u_raw=Uu, Uv_raw=Uv)


VAL = None


def valley(r):
    """U_v, dU_v/dr (envelope), g(r), dg/dr at r >= 1."""
    x = math.log(max(r, 1.0)); x = min(x, VAL["xs"][-1])
    sx, ux = float(VAL["s"](x, 1)), float(VAL["u"](x, 1)); sxx, uxx = float(VAL["s"](x, 2)), float(VAL["u"](x, 2))
    ds, du = sx / r, ux / r                                  # d/dr = (1/r) d/dx
    d2s, d2u = (sxx - sx) / r ** 2, (uxx - ux) / r ** 2
    g = 1.0 + ds * ds + 4.0 * du * du
    gr = 2.0 * ds * d2s + 8.0 * du * d2u
    return float(VAL["Uv"](x)), float(VAL["dUr"](x)), g, gr


# ---------------------------------------------------------------- NPQT
def h(p, al):
    return p / np.sqrt(1.0 - (al * p) ** 2) if al else p


def hprime(p, al):
    return (1.0 - (al * p) ** 2) ** -1.5 if al else 1.0


def hinv(R, al):
    return R / math.sqrt(1.0 + (al * R) ** 2) if al else R


def gfun(x, al):
    """h(x) - x without cancellation."""
    if not al:
        return 0.0 * x
    w = np.sqrt(1.0 - (al * x) ** 2)
    return x * (al * x) ** 2 / (w * (1.0 + w))


def dF(y, k, al):
    """F - F_GR = 3 [g(k)] - 3 int_0^1 t^-2 [g(k + y t^2) - g(k)] dt."""
    if not al:
        return 0.0
    t2 = GL_T ** 2
    return 3.0 * gfun(k, al) - 3.0 * float(np.sum(GL_W * (gfun(k + y * t2, al) - gfun(k, al)) / t2))


# ---------------------------------------------------------------- one wormhole
class Problem:
    def __init__(self, ff, n, lqt, dynamical):
        self.fa = 0.5 * ff * F0_GEV / MPL_GEV
        self.M4 = (self.fa / math.sqrt(0.5 * ZQ)) ** 4
        self.q = n / (2.0 * math.pi ** 2)
        self.al = lqt * lqt
        self.dyn = dynamical

    def Phi(self, r, a):
        return self.q ** 2 / (2.0 * self.fa ** 2 * r * r * a ** 6)

    def rho(self, a, r, rp):
        if not self.dyn:
            return self.Phi(1.0, a)
        Uv, _, g, _ = valley(r)
        return -0.5 * self.fa ** 2 * g * rp * rp + self.M4 * Uv + self.Phi(r, a)

    def throat(self, r0):
        """Outer root of a^2 psi(rho(a, r0, 0)) = 1 (a' = 0)."""
        G = lambda la: math.exp(2 * la) * hinv(self.rho(math.exp(la), r0, 0.0) / 3.0, self.al) - 1.0
        lo = math.log(max(math.sqrt(self.al), 1e-8)) if self.al else math.log(1e-8)
        las = np.linspace(lo, lo + 40.0, 801)
        vals = np.array([G(x) for x in las])
        # the throat is the first crossing from G >= 0 (flux-dominated) to G < 0 going outwards; a later rise of G is the
        # de Sitter region of the potential, not a throat (fix of 25 Nov: the first version took the outermost root)
        down = np.where((vals[:-1] >= 0) & (vals[1:] < 0))[0]
        if len(down) == 0:
            return None
        i = down[0]
        return math.exp(brentq(G, las[i], las[i + 1], xtol=1e-14, rtol=1e-14))

    def rhs(self, tau, y):
        a, ap, r, rp = y[0], y[1], y[2], y[3]
        if self.dyn:
            Uv, dUr, g, gr = valley(r)
            K, Kr = self.fa ** 2 * g, self.fa ** 2 * gr
            Ph = self.Phi(r, a)
            rho = -0.5 * K * rp * rp + self.M4 * Uv + Ph
        else:
            K = Kr = 0.0; Ph = self.Phi(1.0, a); rho = Ph; Uv = dUr = 0.0
        ps = hinv(rho / 3.0, self.al)
        app = -a * ps + a * (2.0 * Ph - K * rp * rp) / (2.0 * hprime(ps, self.al))
        if self.dyn:
            rpp = (self.M4 * dUr - 2.0 * Ph / r - 0.5 * Kr * rp * rp - 3.0 * (ap / a) * K * rp) / K
        else:
            rpp = 0.0
        one_m_ap = a * a * ps / (1.0 + ap) if ap > -1 else 1.0 - ap
        y_e = -(1.0 - a * a * ps) / (a * a)                     # -a'^2/a^2 from the constraint
        dS = (-3.0 * a * one_m_ap ** 2 - a ** 3 * dF(y_e, 1.0 / (a * a), self.al)
              + a ** 3 * ((0.5 * K * rp * rp if self.dyn else 0.0) + (self.M4 * Uv if self.dyn else 0.0) + Ph))
        return [ap, app, rp, rpp, dS]

    def integrate(self, r0, a0, tol=1.0, dense=False):
        y0 = [a0, 0.0, r0, 0.0, 0.0]
        span = (0.0, 1e4 * a0 if not self.dyn else 1e9)
        ev = []
        if self.dyn:
            def over(t, y): return y[2] - 1.0
            over.terminal = True; over.direction = -1
            def under(t, y): return y[3] if t > 1e-3 * a0 else -1.0
            under.terminal = True; under.direction = 1
            def far(t, y): return (y[0] - 1e6 * a0) if abs(y[2] - 1) < 1e-6 else -1.0
            far.terminal = True
            ev = [over, under, far]
        else:
            def far(t, y): return y[0] - 1e4 * a0
            far.terminal = True
            ev = [far]
        return solve_ivp(self.rhs, span, y0, method="DOP853", rtol=1e-10 * tol, atol=1e-14 * tol * a0, events=ev, dense_output=dense)

    def solve(self, tol=1.0):
        if not self.dyn:
            a0 = self.throat(1.0)
            if a0 is None:
                return None
            sol = self.integrate(1.0, a0, tol)
            return self.finish(sol, a0, 1.0)
        # shooting in ln r0: small r0 overshoots (r crosses 1), large r0 undershoots (r turns back up)
        def outcome(lr0):
            r0 = math.exp(lr0); a0 = self.throat(r0)
            if a0 is None:
                return None, None, None
            sol = self.integrate(r0, a0, tol)
            if len(sol.t_events[0]):
                return "OVER", sol, a0
            if len(sol.t_events[1]):
                return "UNDER", sol, a0
            return "OVER" if sol.y[2, -1] < 1.0 else "UNDER", sol, a0
        lo, hi = 1e-9, 16.0
        o_hi = outcome(hi)[0]
        o_lo = outcome(lo)[0]
        if o_lo is None and o_hi is None:
            return None
        # bracket: scan for a sign change among points with a throat
        grid = np.linspace(lo, hi, 65); oc = [outcome(x)[0] for x in grid]
        pairs = [(grid[i], grid[i + 1]) for i in range(64) if oc[i] == "OVER" and oc[i + 1] == "UNDER"]
        if not pairs:
            return dict(throat=False, note="no overshoot/undershoot bracket", outcomes=oc)
        lo, hi = pairs[0]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            o, _, _ = outcome(mid)
            if o == "OVER":
                lo = mid
            else:
                hi = mid
        o, sol, a0 = outcome(lo)
        return self.finish(sol, a0, math.exp(lo))

    def finish(self, sol, a0, r0):
        a, ap, r, rp, S = sol.y
        # tail beyond the end: flux term int a^3 Phi dtau ~ q^2/(4 f^2 r^2 A^2)
        A = a[-1]
        tail = self.q ** 2 / (4.0 * self.fa ** 2 * A * A)
        S_half = 2.0 * math.pi ** 2 * (3.0 * a0 * a0 + S[-1] + tail)
        # constraint residual: 1 - a'^2 - a^2 psi(rho)
        res = 0.0
        for k in range(0, len(sol.t), max(1, len(sol.t) // 400)):
            ps = hinv(self.rho(a[k], r[k], rp[k]) / 3.0, self.al)
            res = max(res, abs(1.0 - ap[k] ** 2 - a[k] ** 2 * ps) / max(1.0, abs(a[k] ** 2 * ps)))
        # validity of the tail (fix of 25 Nov, after the Codex review of PR #106): the tail assumes r = 1 and a nearly flat
        # outer region; a row whose integration stopped before that regime (flux term >= 1e-3 in a'^2, or a tail >= 1% of S)
        # is reported as numerically unresolved and left out of S_min
        flux_end = self.q ** 2 / (6.0 * self.fa ** 2 * r[-1] ** 2 * A ** 4)
        tail_frac = 2.0 * math.pi ** 2 * tail / S_half if S_half > 0 else float("inf")
        resolved = bool(sol.status >= 0 and flux_end < 1e-3 and tail_frac < 0.01)
        return dict(throat=True, a0=a0, r0=r0, S=S_half, constraint=res, a_end=float(A), steps=len(sol.t),
                    flux_end=flux_end, tail_fraction=tail_frac, resolved=resolved)


def gs_analytic(ff, n):
    fa = 0.5 * ff * F0_GEV / MPL_GEV; q = n / (2 * math.pi ** 2)
    a0 = (q * q / (6.0 * fa * fa)) ** 0.25
    return a0, 3.0 * math.pi ** 3 * a0 * a0


def n_min_throat(ff, lqt):
    """smallest (real) n with a throat in (C), by bisection on ln n."""
    has = lambda ln: Problem(ff, math.exp(ln), lqt, False).throat(1.0) is not None
    if has(0.0):
        return 1.0
    lo, hi = 0.0, math.log(1e8)
    if not has(hi):
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if has(mid) else (mid, hi)
    return math.exp(hi)


def task(args):
    ff, n, lqt, dyn, tol = args
    t0 = time.time()
    try:
        r = Problem(ff, n, lqt, dyn).solve(tol)
    except Exception as e:               # reported, not hidden
        r = dict(throat=False, error=repr(e))
    if r is None:
        r = dict(throat=False)
    r.update(F_over_F0=ff, n=float(n), l_QT=lqt, version="D" if dyn else "C", tol=tol, runtime_s=time.time() - t0)
    return r


def init_worker():
    global VAL
    VAL = build_valley()


def gate_G3(rng):
    errs = []
    for _ in range(20):
        lqt = rng.choice([300.0, 1000.0, 2000.0]); al = lqt * lqt
        psi = rng.uniform(0, 0.999) / al
        k = rng.uniform(0, 1) * (1.0 / al)
        y = psi - k                                   # psi = k + y
        Fv = 3 * k - 3 * y + dF(y, k, al)
        FY = -3.0 * float(np.sum(GL_W * hprime(k + y * GL_T ** 2, al)))      # dF/dY = -3 int_0^1 h'(k + Y t^2) dt (fix of 25 Nov:
        E = Fv - 2 * y * FY                                                  # the first version used finite differences)
        errs.append(abs(E / (3 * h(psi, al)) - 1) if psi > 0 else abs(E))
    return max(errs)


def main(argv):
    global VAL
    quick = "--quick" in argv
    OUT.mkdir(parents=True, exist_ok=True)
    VAL = build_valley()
    Sreq = {ff: S_required(ff) for ff in F_ROWS}
    rng = np.random.default_rng(70)
    g3 = gate_G3(rng)
    log = open(OUT / "run.log", "a")
    def lg(m): print(m, flush=True); log.write(m + "\n"); log.flush()
    lg(f"U70 start {time.strftime('%Y-%m-%d %H:%M')}  S_req={Sreq}  G3 max rel err={g3:.2e}")

    # G1 and G2 (alpha = 0)
    g1 = []
    for ff in F_ROWS:
        a0, Sa = gs_analytic(ff, 1.0)
        r = task((ff, 1.0, 0.0, False, 1.0))
        g1.append(dict(F_over_F0=ff, a0_num=r.get("a0"), a0_analytic=a0, S_num=r.get("S"), S_analytic=Sa,
                       rel_err=abs(r["S"] / Sa - 1) if r.get("throat") else None))
        lg(f"G1 F={ff}: S={r.get('S')}, analytic={Sa}, a0={r.get('a0')} vs {a0}")
    g2rows = [task((ff, 1.0, 0.0, True, 1.0)) for ff in F_ROWS]
    for r in g2rows:
        lg(f"G2 F={r['F_over_F0']}: {r}")
    if quick:
        return
    jobs = []
    nmin = {}
    for lqt in LQT_ROWS:
        for ff in F_ROWS:
            nm = n_min_throat(ff, lqt) if lqt else 1.0
            nmin[f"{lqt:g}_{ff}"] = nm
            ns = list(N_GRID) + ([math.ceil(nm)] if nm and nm > 1 else [])
            for n in sorted(set(ns)):
                for dyn in (False, True):
                    jobs.append((ff, float(n), lqt, dyn, 1.0))
    lg(f"n_min(C) = {nmin};  {len(jobs)} jobs")
    rows = []
    with Pool(4, initializer=init_worker) as pool:
        for r in pool.imap_unordered(task, jobs, chunksize=2):
            rows.append(r)
            if r.get("throat"):
                lg(f"row l={r['l_QT']:g} F={r['F_over_F0']} n={r['n']:g} {r['version']}: a0={r['a0']:.4g} r0={r['r0']:.4g} S={r['S']:.6g} res={r['constraint']:.1e}")
    (OUT / "rows.json").write_text(json.dumps(rows, default=float))
    lg("rows done")


if __name__ == "__main__":
    main(sys.argv[1:])
