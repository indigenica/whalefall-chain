"""U59 (docs/preregistros/ESPEC_U59.md): the scalar perturbations through the cuscuton bounce, independent audit.
The quadratic ADM action of Einstein + cuscuton + chi + radiation (P = X^2) on S^3 is derived here with sympy from the
covariant action, with an explicit l = 0 harmonic Y(r) = sin((n+1) r)/sin r, expanded to second order point by point and
integrated over S^3 (no code of the explorer is used; only its variable conventions and published numbers are compared).

Conventions (those of the explorer's C2, for comparison): h_ij = a^2 e^{2 zeta Y} gamma_ij, N = 1 + alpha Y, N_i = D_i(beta Y),
(delta chi, delta theta, delta phi) = (x, w, f) Y;  L2 / a^3 = 1/2 v K0 v + v B0 q + 1/2 q C0 q + v T u + u R q + 1/2 u S u,
q = (zeta, x, w), v = qdot, u = (alpha, beta, f), after integrating the zeta fdot term by parts; Y normalised on S^3.

    python3 -B chain/u59_escalares_rebote.py      ->  results/u59_escalares_rebote.json
Internal, not exported."""
import json, math, pickle, time
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from ucommon import ROOT, save
import explora_cuscuton as EC

CACHE = ROOT / "results" / "u59_runs" / "L2_symbolic.pkl"
ADOPTED = dict(mu=0.16, m=0.02, pinf=0.25)
C1 = dict(mu=0.3, m=0.05, pinf=0.25)
NS = (2, 3, 10, 100, 1000, 10000, 30000, 100000, 300000)

# ---------------------------------------------------------------- symbols
r, th = sp.symbols("r theta", positive=True)
eps = sp.Symbol("eps")
a, H, chid, phid, q, mu2 = sp.symbols("a H chid phid q mu2", positive=True)
Vv, Vp, Vpp, Uv, Up, Upp = sp.symbols("Vv Vp Vpp Uv Up Upp")
zeta, x, w, al, be, f = sp.symbols("zeta x w al be f")
zd, xd, wd, fd = sp.symbols("zd xd wd fd")
Ys, Y1, Y2 = sp.symbols("Ys Y1 Y2")
QV = (zeta, x, w); VV = (zd, xd, wd); UV = (al, be, f)


def derive():
    """L2 (the eps^2 coefficient of the ADM Lagrangian density, with sqrt(gamma) = sin^2 r, angles integrated out) as a
    polynomial in (Ys, Y1, Y2) with r-dependent coefficients"""
    Y = sp.Function("Y")(r)
    wc = sp.log(a) + eps * zeta * Y                                       # h = e^{2 wc} gamma
    g = sp.diag(1, sp.sin(r) ** 2, sp.sin(r) ** 2 * sp.sin(th) ** 2)
    h = sp.exp(2 * wc) * g; hinv = h.inv()
    X = (r, th, sp.Symbol("ph"))
    Gam = [[[sum(hinv[k, l] * (sp.diff(h[l, i], X[j]) + sp.diff(h[l, j], X[i]) - sp.diff(h[i, j], X[l])) / 2 for l in range(3))
             for j in range(3)] for i in range(3)] for k in range(3)]
    N = 1 + eps * al * Y
    Ncov = [eps * be * sp.diff(Y, r), 0, 0]
    DN = sp.Matrix(3, 3, lambda i, j: sp.diff(Ncov[j], X[i]) - sum(Gam[k][i][j] * Ncov[k] for k in range(3)))
    hdot = 2 * (H + eps * zd * Y) * h
    K = (hdot - DN - DN.T) / (2 * N)
    Kup = hinv * K * hinv
    KK = sum(K[i, j] * Kup[i, j] for i in range(3) for j in range(3)); Ktr = sum(hinv[i, j] * K[i, j] for i in range(3) for j in range(3))
    wr = sp.diff(wc, r)
    R3 = sp.exp(-2 * wc) * (6 - 4 * sp.diff(sp.sin(r) ** 2 * wr, r) / sp.sin(r) ** 2 - 2 * wr ** 2)
    sqrth = sp.exp(3 * wc) * sp.sin(r) ** 2                              # angles: sin(theta) integrated with the 4 pi below
    Nr = hinv[0, 0] * Ncov[0]
    def kin(dot, rad):                                                  # (field_dot - N^r field_r)^2 / N^2 - h^rr field_r^2
        return (dot - Nr * rad) ** 2 / N ** 2 - hinv[0, 0] * rad ** 2
    xx = eps * x * Y
    Lchi = N * sqrth * (sp.Rational(1, 2) * kin(chid + eps * xd * Y, sp.diff(xx, r)) - (Vv + Vp * xx + Vpp * xx ** 2 / 2))
    ff = eps * f * Y
    twoXphi = kin(phid + eps * fd * Y, sp.diff(ff, r))
    Lphi = N * sqrth * (-mu2 * sp.sqrt(twoXphi) - (Uv + Up * ff + Upp * ff ** 2 / 2))
    Xth = kin(q + eps * wd * Y, sp.diff(eps * w * Y, r)) / 2
    Lrad = N * sqrth * Xth ** 2
    Lg = N * sqrth * (R3 / 2 + (KK - Ktr ** 2) / 2)
    Ltot = Lg.subs(th, sp.pi / 2) + Lchi + Lphi + Lrad                  # th only enters through h_phph: sin(th)=1 after the ratios
    L2 = sp.diff(Ltot, eps, 2).subs(eps, 0) / 2
    L2 = L2.subs(sp.Derivative(Y, (r, 2)), Y2).subs(sp.Derivative(Y, r), Y1).subs(Y, Ys)
    L2 = sp.expand(sp.simplify(L2))
    return L2


def load_L2():
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    if CACHE.exists():
        return pickle.loads(CACHE.read_bytes())
    t0 = time.time(); L2 = derive(); CACHE.write_bytes(pickle.dumps(L2))
    print(f"derived L2 in {time.time() - t0:.0f}s", flush=True)
    return L2


def split(L2):
    """L2 = sum over (pair of perturbation symbols) coef_pair(bg) * integrand_pair(r, Ys, Y1, Y2)"""
    pert = list(QV + VV + (fd,) + UV)
    terms = {}
    P = sp.Poly(L2, *pert)
    for mon, c in P.terms():
        terms[mon] = c
    return pert, terms


def harmonic_integrals(expr_r, n, npts=None):
    """integral over S^3 (4 pi * int_0^pi dr) of expr(r, Ys, Y1, Y2) with Y = sin((n+1) r)/(pi sqrt(2) sin r) (normalised)"""
    npts = npts or max(400, 8 * (n + 2))
    panels = max(20, npts // 16); xg, wg = leggauss(16)                # composite Gauss-Legendre, 16 points per panel
    edges = np.linspace(0, math.pi, panels + 1); hw = 0.5 * (edges[1:] - edges[:-1]); mid = 0.5 * (edges[1:] + edges[:-1])
    rr = (mid[:, None] + hw[:, None] * xg[None, :]).ravel(); ww = (hw[:, None] * wg[None, :]).ravel()
    k = n + 1; s, c = np.sin(rr), np.cos(rr); nrm = 1 / (math.pi * math.sqrt(2))
    Y = nrm * np.sin(k * rr) / s
    Yr = nrm * (k * np.cos(k * rr) * s - np.sin(k * rr) * c) / s ** 2
    # Y'' from the eigen-equation of the radial Laplacian on S^3: Y'' + 2 cot r Y' = -n(n+2) Y
    Yrr = -n * (n + 2) * Y - 2 * c / s * Yr
    fn = expr_r if callable(expr_r) else sp.lambdify((r, Ys, Y1, Y2), expr_r, "numpy")
    vals = fn(rr, Y, Yr, Yrr) * np.ones_like(rr)
    return float(4 * math.pi * np.sum(ww * vals))


class Action:
    """matrices K0, B0, C0, T, R, S (per a^3) as functions of the background, for any n via a polynomial in k^2 = n(n+2)"""
    def __init__(self):
        L2 = load_L2(); self.pert, terms = split(L2)
        bgs = (a, H, chid, phid, q, mu2, Vv, Vp, Vpp, Uv, Up, Upp)
        self.items = []
        for mon, c in terms.items():
            # c = sum_j coefbg_j * rfun_j(r, Ys, Y1, Y2): separate the r part from the background part
            c = sp.expand(c)
            parts = {}
            for term in sp.Add.make_args(c):
                bgpart = sp.Mul(*[fct for fct in sp.Mul.make_args(term) if not (fct.free_symbols & {r, Ys, Y1, Y2})])
                rpart = sp.powsimp(term / bgpart)
                parts.setdefault(rpart, 0); parts[rpart] += bgpart
            for rp, bp in parts.items():
                self.items.append((mon, sp.lambdify(bgs, bp, "numpy"), sp.lambdify((r, Ys, Y1, Y2), rp, "numpy")))
        self.cache = {}

    def integrals(self, n):
        if n not in self.cache:
            self.cache[n] = [harmonic_integrals(rp, n) for (_, _, rp) in self.items]
        return self.cache[n]

    def fit(self, ns=(2, 3, 4, 5, 6, 8)):
        """each integral as a polynomial of degree <= 3 in k^2 (least squares on ns), checked at n = 37 (gate G1)"""
        K2 = np.array([n * (n + 2) for n in ns], float)
        M = np.vstack([np.ones_like(K2), K2, K2 ** 2, K2 ** 3]).T
        I = np.array([self.integrals(n) for n in ns])
        self.coef = np.linalg.lstsq(M, I, rcond=None)[0]
        I37 = np.array(self.integrals(37)); k2 = 37 * 39.0
        P37 = self.coef[0] + self.coef[1] * k2 + self.coef[2] * k2 ** 2 + self.coef[3] * k2 ** 3
        scale = np.maximum(np.abs(I37), 1e-12 * np.max(np.abs(I37)))
        return float(np.max(np.abs(P37 - I37) / np.max(np.abs(I37)))), float(np.max(np.abs(P37 - I37) / scale))

    def quad_check(self, n, factor=2):
        """amendment 1 (gate G1): the integrals with factor x more Gauss points change by <= 1e-10 relative"""
        I1 = np.array(self.integrals(n)); I2 = np.array([harmonic_integrals(rp, n, factor * max(400, 8 * (n + 2))) for (_, _, rp) in self.items])
        return float(np.max(np.abs(I2 - I1)) / np.max(np.abs(I1)))

    def matrices(self, bg, n):
        I = self.integrals(n)                                           # amendment 1: direct quadrature at every n
        idx = {s: i for i, s in enumerate(self.pert)}
        M = np.zeros((10, 10))
        args = tuple(bg[k] for k in ("a", "H", "chid", "phid", "q", "mu2", "Vv", "Vp", "Vpp", "Uv", "Up", "Upp"))
        for (mon, fb, _), Ii in zip(self.items, I):
            val = float(fb(*args)) * Ii
            nz = [i for i, e_ in enumerate(mon) for _ in range(e_)]
            if len(nz) == 2:
                i, j = nz
                if i == j: M[i, i] += 2 * val
                else: M[i, j] += val; M[j, i] += val
        M /= bg["a"] ** 3
        iq = [idx[s] for s in QV]; iv = [idx[s] for s in VV]; iu = [idx[s] for s in UV]; ifd = idx[fd]
        # the zeta fdot term, c * zeta * fdot with c = M[zeta, fd] (per a^3): by parts -> -(3 H c) zeta f - c zetadot f
        c = M[idx[zeta], ifd]
        K0 = M[np.ix_(iv, iv)]; B0 = M[np.ix_(iv, iq)]; C0 = M[np.ix_(iq, iq)]
        T = M[np.ix_(iv, iu)].copy(); R = M[np.ix_(iu, iq)].copy(); S = M[np.ix_(iu, iu)]
        T[0, 2] += -c                                                    # zetadot * f
        R[2, 0] += -3 * bg["H"] * c                                      # f * zeta
        others = [M[ifd, j] for j in range(10) if j not in (idx[zeta],)]
        return dict(K0=K0, B0=B0, C0=C0, T=T, R=R, S=S, c_zeta_fdot=c, fd_other_max=float(np.max(np.abs(others))))


def reduce(m):
    S = m["S"]; Si_T = np.linalg.solve(S, m["T"].T); Si_R = np.linalg.solve(S, m["R"])
    K = m["K0"] - m["T"] @ Si_T; B = m["B0"] - m["T"] @ Si_R; C = m["C0"] - m["R"].T @ Si_R
    return dict(K=(K + K.T) / 2, B=B, C=(C + C.T) / 2, detS_rel=abs(np.linalg.det(S)) / np.prod(np.linalg.norm(S, axis=1)))


# ---------------------------------------------------------------- backgrounds
def bg_from_state(cu, Crad, y):
    """y = (lna, chi, chidot, phi) in mode 0 of explora_cuscuton; the background dict the matrices need"""
    Q, _ = EC.deriv(cu, 0, Crad, 0.0)
    qq = Q(y); phi, U, dU, ddU = cu.all(0, y[3])
    chi = y[1]; Cp = EC.CP; V0 = EC.V0
    Vpp_ = V0 * Cp * (Cp * Cp - 2 * chi * chi) / (Cp * Cp + chi * chi) ** 2.5
    return dict(a=math.exp(y[0]), H=qq["H"], chid=y[2], phid=qq["phidot"], q=(4 * qq["rho"] / 3) ** 0.25, mu2=cu.mu2,
                Vv=EC.Vp(chi), Vp=EC.dVp(chi), Vpp=Vpp_, Uv=U, Up=dU, Upp=ddU, D=qq["D"], W=qq["W"], phi=phi)


def window(cu, Crad, y_b, phis):
    """integrate in t from the bounce state y_b (phi = 0) to phi = -0.2 and +0.2; returns dense solutions both ways"""
    _, rhs = EC.deriv(cu, 0, Crad, 0.0)
    out = {}
    for sgn, ph in ((-1, -0.2), (1, 0.2)):
        ev = lambda t, y: y[3] - ph; ev.terminal = True
        s = solve_ivp(rhs, (0, sgn * 1e7), y_b, method="DOP853", rtol=1e-12, atol=[1e-14, 1e-15, 1e-24, 1e-24], events=ev, dense_output=True)
        out[sgn] = s
    return out


# ---------------------------------------------------------------- production
PUB_C2 = {2: (0.7175965997875235, 2.6095240811535536), 10: (None, 3.4793120085846123), 1000: (None, 11.563681480021323)}
PUB_C3E = {2: 12.3027, 3: 10.7995, 10: 9.8412, 100: 9.5579, 1000: 2.5548, 3000: 2.0493, 10000: 1.7114, 20000: 1.9949, 30000: 1.6026}
HOST = dict(ADOPTED, AH=100.0)


def host_state(par):
    """the bounce state of a mother (explora_cuscuton.run) -> cuscuton, C_rad, y_b = (ln a, chi, chidot, phi = 0)"""
    cu = EC.Cusc(par["mu"], par["m"], par["pinf"])
    o = EC.run(par); b = o["bounce"]
    return cu, b["rho_r"] * b["a"] ** 4, [math.log(b["a"]), b["chi"], b["chidot"], 0.0], o


class Trajectory:
    def __init__(self, cu, Crad, y_b):
        self.cu, self.Crad = cu, Crad
        self.sol = window(cu, Crad, y_b, None)
        self.tm, self.tp = float(self.sol[-1].t[-1]), float(self.sol[1].t[-1])

    def y(self, t):
        return self.sol[-1].sol(t) if t < 0 else self.sol[1].sol(t)

    def bg(self, t):
        return bg_from_state(self.cu, self.Crad, list(self.y(t)))

    def t_of_phi(self, ph):
        from scipy.optimize import brentq
        s = self.sol[-1] if ph < 0 else self.sol[1]
        lo, hi = (self.tm, 0.0) if ph < 0 else (0.0, self.tp)
        return brentq(lambda t: s.sol(t)[3] - ph, lo, hi, xtol=1e-14 * max(abs(lo), abs(hi)))


def red2(A, bg, n):
    m = A.matrices(bg, n); R = reduce(m)
    return R["K"][:2, :2], R["B"][:2, :2], R["C"][:2, :2], R["detS_rel"]


def energy_matrix(K, B, C):
    Ki = np.linalg.inv(K)
    return np.block([[B.T @ Ki @ B - C, -B.T @ Ki], [-Ki @ B, Ki]])


def transfer(A, tr, n, rtol):
    def rhs(t, y):
        bg = tr.bg(t); K, B, C, _ = red2(A, bg, n); a3 = bg["a"] ** 3
        Y = y.reshape(4, 4); out = np.empty_like(Y)
        for j in range(4):
            qv, p = Y[j, :2], Y[j, 2:]
            qd = np.linalg.solve(K, p / a3 - B @ qv)
            out[j, :2] = qd; out[j, 2:] = a3 * (B.T @ qd + C @ qv)
        return out.ravel()
    bgi, bgf = tr.bg(tr.tm), tr.bg(tr.tp); a3i, a3f = bgi["a"] ** 3, bgf["a"] ** 3
    # initial basis in z = (q, P): unit vectors; p = a^3 P
    Z0 = np.eye(4); Y0 = Z0.copy(); Y0[:, 2:] *= a3i
    s = solve_ivp(rhs, (tr.tm, tr.tp), Y0.ravel(), method="DOP853", rtol=rtol, atol=1e-300)
    Yf = s.y[:, -1].reshape(4, 4); Zf = Yf.copy(); Zf[:, 2:] /= a3f
    Tz = Zf.T                                                            # columns: images of the basis
    Ki, Bi, Ci, _ = red2(A, bgi, n); Kf, Bf, Cf, _ = red2(A, bgf, n)
    Ei, Ef = energy_matrix(Ki, Bi, Ci), energy_matrix(Kf, Bf, Cf)
    Ei, Ef = 0.5 * (Ei + Ei.T), 0.5 * (Ef + Ef.T)
    from scipy.linalg import sqrtm
    pos = bool(np.min(np.linalg.eigvalsh(Ei)) > 0 and np.min(np.linalg.eigvalsh(Ef)) > 0)
    G = float(np.max(np.linalg.svd(np.real(sqrtm(Ef)) @ Tz @ np.linalg.inv(np.real(sqrtm(Ei))), compute_uv=False))) if pos else None
    # symplectic in (q, p)
    Tqp = np.diag([1, 1, a3f, a3f]) @ Tz @ np.diag([1, 1, 1 / a3i, 1 / a3i])
    J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    symp = float(np.max(np.abs(Tqp.T @ J @ Tqp - J)) / np.max(np.abs(Tqp)) ** 2)
    a_b = math.exp(tr.y(0.0)[0])
    return dict(n=n, G=G, energy_positive=pos, symp=symp, nfev=int(s.nfev), k_phys_bounce=math.sqrt(n * (n + 2.0)) / a_b)


def main():
    t0 = time.time(); A = Action(); out = dict(spec="docs/preregistros/ESPEC_U59.md")
    # G0: C1 bounce against the published C2 eigenvalues
    ref = json.loads((ROOT / "data" / "explorador" / "cuscuton_c1" / "result.json").read_text())["bounce"]
    cu1 = EC.Cusc(**C1); _, Uv0, dU0, ddU0 = cu1.all(0, 0.0); chi = ref["chi"]; Cp, V0 = EC.CP, EC.V0
    bg1 = dict(a=ref["a"], H=0.0, chid=ref["chi_dot"], phid=ref["phidot"], q=(4 * ref["rho_r"] / 3) ** 0.25, mu2=cu1.mu2, Vv=EC.Vp(chi),
               Vp=EC.dVp(chi), Vpp=V0 * Cp * (Cp * Cp - 2 * chi * chi) / (Cp * Cp + chi * chi) ** 2.5, Uv=Uv0, Up=dU0, Upp=ddU0)
    g0 = {}
    for n, (lo, hi) in PUB_C2.items():
        ev = np.linalg.eigvalsh(red2(A, bg1, n)[0])
        g0[n] = dict(eig=[float(e) for e in ev], rel=max([abs(ev[1] / hi - 1)] + ([abs(ev[0] / lo - 1)] if lo else [])))
    G0 = all(v["rel"] <= 1e-8 for v in g0.values())
    G1v = dict(n37=A.quad_check(37), n3e5=A.quad_check(300000)); G1 = all(v <= 1e-10 for v in G1v.values())
    out.update(G0=g0, G1=G1v)
    print("G0", g0, "G1", G1v, flush=True)
    res = {}
    for label, par in (("adoptado", HOST | dict(f=None)), ("C1", None)):
        if label == "adoptado":
            rows = json.loads((ROOT / "results" / "u57b_runs" / "rows_rtol11.json").read_text())["rows"]
            k = sorted([r for r in rows if r["set"] == "eps_pequeno" and r["par"]["AH"] == 100.0 and r["outcome"] == "HIJO"], key=lambda r: r["par"]["f"])[0]
            cu, Crad, yb, o = host_state(dict(ADOPTED, AH=100.0, f=k["par"]["f"]))
        else:
            cu = cu1; o = EC.run(EC.c1_row()); b = o["bounce"]; Crad = b["rho_r"] * b["a"] ** 4; yb = [math.log(b["a"]), b["chi"], b["chidot"], 0.0]
        tr = Trajectory(cu, Crad, yb)
        # positivity scan
        phis = np.linspace(-0.2, 0.2, 81); kmin = math.inf; dmin = math.inf; Dmin = math.inf; Wmin = math.inf
        for ph in phis:
            t = tr.t_of_phi(ph) if abs(ph) > 1e-12 else 0.0; bg = tr.bg(t)
            Dmin = min(Dmin, bg["D"]); Wmin = min(Wmin, bg["W"])
            for n in NS:
                K, _, _, dS = red2(A, bg, n); ev = np.linalg.eigvalsh(K)
                kmin = min(kmin, ev[0] / ev[1]); dmin = min(dmin, dS)
        tf = []
        for n in NS:
            tf.append(transfer(A, tr, n, 1e-10)); print(label, tf[-1], flush=True)
        conv = {n: abs(transfer(A, tr, n, 1e-9)["G"] / next(x["G"] for x in tf if x["n"] == n) - 1) for n in (2, 1000, 100000)}
        res[label] = dict(min_K_ratio=kmin, min_detS_rel=dmin, D_min=Dmin, W_min=Wmin, a_b=math.exp(yb[0]), t_window=[tr.tm, tr.tp],
                          transfers=tf, conv_rtol=conv)
    ad = res["adoptado"]; Gs = {x["n"]: x["G"] for x in ad["transfers"]}
    G2 = all(v <= 1e-4 for v in ad["conv_rtol"].values()) and all(v <= 1e-4 for v in res["C1"]["conv_rtol"].values())
    G3 = all(x["symp"] <= 1e-6 for r_ in res.values() for x in r_["transfers"])
    G4 = all(x["energy_positive"] for r_ in res.values() for x in r_["transfers"])
    ratio = Gs[300000] / Gs[10000] if Gs[300000] and Gs[10000] else None
    sane = bool(ad["min_K_ratio"] > 0 and ad["min_detS_rel"] > 1e-12 and ratio is not None and ratio <= 2)
    unstable = bool(ad["min_K_ratio"] < -1e-9 or (ratio is not None and ratio > 2 and Gs[30000] < Gs[100000] < Gs[300000]))
    gates = dict(G0=G0, G1=G1, G2=G2, G3=G3, G4=G4)
    verdict = "INCONCLUSO" if not all(gates.values()) else ("ESCALARES_SANOS" if sane else ("INESTABILIDAD" if unstable else "INCONCLUSO"))
    out.update(gates=gates, ratio_G_3e5_over_1e4=ratio, would_be=("ESCALARES_SANOS" if sane else ("INESTABILIDAD" if unstable else "INCONCLUSO")),
               verdict=verdict, results=res, C3e_published=PUB_C3E, seconds=time.time() - t0)
    save("u59_escalares_rebote.json", out)
    print(json.dumps(dict(gates=gates, ratio=ratio, would_be=out["would_be"], verdict=verdict,
                          adopted=dict((k_, ad[k_]) for k_ in ("min_K_ratio", "min_detS_rel", "D_min", "W_min", "a_b", "conv_rtol")),
                          G_adopted=Gs, G_C1={x["n"]: x["G"] for x in res["C1"]["transfers"]}), indent=1, default=float))


if __name__ == "__main__":
    main()
