"""U58 (docs/preregistros/ESPEC_U58.md): the cuscuton in the pocket's inflation and the late universe.
(1) the observational bound on eps_G (Afshordi 2009, eq. 13-14; Robbers, Afshordi & Doran 2008): -0.025 <= eps_G <= 0.017;
(2) n_s and r from the quadratic action of Boruah, Kim & Geshnizjani (2017, eqs. 3.16-3.24) on the pole plateau of the chain,
    flat FLRW, with the Boruah cuscuton at its late attractor, for GR, the small set and the C1 set.
    python3 -B chain/u58_cuscuton_inflacion.py   ->  results/u58_cuscuton_inflacion.json
Internal, not exported."""
import json, math, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, benchmark, save
import explora_cuscuton as EC

B = benchmark(); V0 = B["V0"]; CP = EC.CP
SETS = {"GR": None, "pequeno": dict(mu=0.16, m=0.02, pinf=0.25), "C1": dict(mu=0.3, m=0.05, pinf=0.25)}
BOUND = (0.003 - 0.014, 0.003 + 0.014)                          # lambda - 1 at 95 % (Afshordi 2009, eq. 14)
EPS_BOUND = (-1.5 * BOUND[1], -1.5 * BOUND[0])                  # eps_G = (3/2)(1 - lambda)


def V(x): return V0 * (1 - CP / math.hypot(CP, x))
def dV(x): return V0 * CP * x / (CP * CP + x * x) ** 1.5


class Model:
    def __init__(self, par):
        self.cu = EC.Cusc(**par) if par else None

    def d_of_H(self, H):
        """phi = phi_inf - d with U'(phi) = 3 mu^2 H; Newton from the quadratic guess (H is tiny, d ~ 3 mu^2 H / (4 phi_inf^2))"""
        cu = self.cu; t = 3 * cu.mu2 * H; d = t / (4 * cu.pinf ** 2)
        for _ in range(60):
            _, U, dU, ddU = cu.all(1, d)
            step = (dU - t) / (-ddU)                              # d(U')/dd = -U''(phi)
            d -= step
            if abs(step) <= 1e-16 * d: break
        return d

    def H_of(self, x, xN):
        if self.cu is None:
            return math.sqrt(V(x) / (3 - 0.5 * xN * xN))
        H = math.sqrt(V(x) / (3 - 0.5 * xN * xN))
        for _ in range(80):                                       # 3 H^2 (1 - xN^2/6) = V + U(d(H)); U ~ -3 eps_G H^2
            Hn = math.sqrt((V(x) + self.cu.all(1, self.d_of_H(H))[1]) / (3 - 0.5 * xN * xN))
            if abs(Hn / H - 1) <= 1e-16: H = Hn; break
            H = Hn
        return H

    def quantities(self, y):
        x, xN = y; H = self.H_of(x, xN); xd = H * xN
        if self.cu is None:
            phid = 0.0; D = None; U = 0.0
        else:
            d = self.d_of_H(H); _, U, dU, ddU = self.cu.all(1, d); D = 1.5 * self.cu.mu4 - ddU
            phid = 3 * self.cu.mu2 * 0.5 * xd * xd / D
        Hd = -0.5 * xd * xd + (0.5 * self.cu.mu2 * phid if self.cu else 0.0)
        eps = -Hd / H ** 2; alpha = 0.5 * xN * xN
        res = abs(3 * H * H - 0.5 * xd * xd - V(x) - U) / (3 * H * H)
        return dict(H=H, eps=eps, alpha=alpha, sigma=eps - alpha, D=D, phid=phid, res=res)

    def rhs(self, N, y):
        q = self.quantities(y)
        return [y[1], -(3 - q["eps"]) * y[1] - dV(y[0]) / q["H"] ** 2]


def background(model, rtol):
    x0 = 6.0; H0 = math.sqrt(V(x0) / 3); xN0 = -dV(x0) / (3 * H0 * H0)
    ev = lambda N, y: model.quantities(y)["eps"] - 1; ev.terminal = True; ev.direction = 1
    s = solve_ivp(model.rhs, (0, 400), [x0, xN0], method="DOP853", rtol=rtol, atol=1e-14, events=ev, dense_output=True)
    return s, float(s.t_events[0][0])


def make_bg(model, s, Nend, n=40001):
    """the background tabulated on a fine grid of N and splined: a H, H, alpha, sigma, eps, eta = eps_N/eps, beta = alpha_N/alpha"""
    from scipy.interpolate import CubicSpline
    Ns = np.linspace(0, Nend, n); Q = [model.quantities(s.sol(N)) for N in Ns]
    g = lambda k: np.array([q[k] for q in Q])
    sp = {k: CubicSpline(Ns, g(k)) for k in ("H", "alpha", "sigma", "eps")}
    lnH = CubicSpline(Ns, np.log(g("H")))
    D = g("D") if model.cu else None; phid = g("phid"); res = g("res")
    def f(N):
        ep, al = float(sp["eps"](N)), float(sp["alpha"](N))
        return dict(aH=math.exp(N + float(lnH(N))), H=math.exp(float(lnH(N))), alpha=al, sigma=float(sp["sigma"](N)), eps=ep,
                    eta=float(sp["eps"](N, 1)) / ep, beta=float(sp["alpha"](N, 1)) / al)
    f.res_max = float(res.max()); f.D_min = float(D.min()) if D is not None else None
    f.phid_min = float(phid.min()) if model.cu else None
    return f


def z2_cs2(b, k, N, cs_one):
    a = math.exp(N); cH = b["aH"]; al, sg, ep, et, be = b["alpha"], b["sigma"], b["eps"], b["eta"], b["beta"]
    z2 = a * a * al * (k * k + 3 * al * cH * cH) / (k * k + al * cH * cH * (3 - sg))
    if cs_one:
        return z2, 1.0
    A1 = 6 * al - al * sg; A2 = 9 * al * al - 3 * al * al * sg
    B1 = A1 + sg * (6 + et + be - 2 * ep) + al * (et - be)
    B2 = A2 + sg * al * (12 - 4 * sg + 3 * et) + 3 * al * al * (et - be)
    c2 = (k ** 4 + k * k * cH * cH * B1 + cH ** 4 * B2) / (k ** 4 + k * k * cH * cH * A1 + cH ** 4 * A2)
    return z2, c2


def mode(bg, k, Ncross, rtol, cs_one=False, tensor=False):
    """zeta (or h) of comoving mode k from k = 100 aH to k = aH/1000; returns P and the minimum of z^2"""
    Ni = Ncross - math.log(100.0); Nf = Ncross + math.log(1000.0)
    zmin = [math.inf]

    def coef(N):
        b = bg(N)
        if tensor:
            return b, math.exp(2 * N), 1.0                       # "z^2" = a^2 for h with Q = a^2 h'
        z2, c2 = z2_cs2(b, k, N, cs_one); zmin[0] = min(zmin[0], z2); return b, z2, c2

    def rhs(N, y):
        b, z2, c2 = coef(N); cH = b["aH"]
        zr, zi, pr, pi_ = y
        return [pr / (z2 * cH), pi_ / (z2 * cH), -z2 * c2 * k * k * zr / cH, -z2 * c2 * k * k * zi / cH]
    b, z2, c2 = coef(Ni); cH = b["aH"]
    h = 1e-5; _, z2p, _ = coef(Ni + h); _, z2m, _ = coef(Ni - h)
    dlnz_dtau = cH * (math.log(z2p) - math.log(z2m)) / (4 * h)  # z'/z, conformal
    v = 1 / math.sqrt(2 * k); vd = -1j * k * v
    if tensor:                                                   # h = 2 v / a, Q = a^2 h'
        a = math.exp(Ni); hh = 2 * v / a; hd = 2 * (vd - cH * v) / a; Q = a * a * hd; y0 = hh
    else:                                                        # v = sqrt(2) z zeta, Pi = z^2 zeta'
        z = math.sqrt(z2); y0 = v / (math.sqrt(2) * z); zd = (vd - dlnz_dtau * v) / (math.sqrt(2) * z); Q = z2 * zd
    s = solve_ivp(rhs, (Ni, Nf), [y0.real, y0.imag, Q.real, Q.imag], method="DOP853", rtol=rtol, atol=1e-300)
    amp2 = s.y[0, -1] ** 2 + s.y[1, -1] ** 2
    P = (k ** 3 * amp2 / math.pi ** 2) if tensor else (k ** 3 * amp2 / (2 * math.pi ** 2))
    return P, zmin[0]


def run_set(name, rtol, cs_one=False):
    model = Model(SETS[name]); s, Nend = background(model, rtol); bg = make_bg(model, s, Nend)
    Npiv = Nend - B["N_star"]; bp = bg(Npiv); kp = bp["aH"]
    out = dict(set=name, N_end=Nend, N_piv=Npiv, chi_piv=float(s.sol(Npiv)[0]), rtol=rtol, cs_one=cs_one)
    P = {}; zmins = []
    for dl in (-0.5, 0.0, 0.5):
        k = kp * math.exp(dl)
        Nc = brentq(lambda N: bg(N)["aH"] - k, Npiv - 2, Npiv + 2)
        P[dl], zm = mode(bg, k, Nc, rtol, cs_one); zmins.append(zm)
    Pt, _ = mode(bg, kp, Npiv, rtol, tensor=True)
    out.update(P_zeta=P[0.0], n_s=1 + (math.log(P[0.5]) - math.log(P[-0.5])) / 1.0, r=Pt / P[0.0], P_t=Pt,
               P_sr=bp["H"] ** 2 / (8 * math.pi ** 2 * bp["alpha"]), z2_min=min(zmins))
    out.update(res_max=bg.res_max, D_min=bg.D_min, phid_min=bg.phid_min,
               sigma_over_eps_piv=bp["sigma"] / bp["eps"], eps_G=(model.cu.eps_G() if model.cu else 0.0))
    return out


def main():
    t0 = time.time(); R = {}
    for name in SETS:
        R[name] = dict(main=run_set(name, 1e-10), rtol9=run_set(name, 1e-9))
        if SETS[name]:
            R[name]["cs1"] = run_set(name, 1e-10, cs_one=True)
        print(name, json.dumps({k: R[name]["main"][k] for k in ("n_s", "r", "P_zeta", "P_sr", "N_end", "chi_piv", "res_max", "D_min", "phid_min", "z2_min", "sigma_over_eps_piv", "eps_G")}, default=float), flush=True)
    g = R["GR"]["main"]
    G0 = bool(abs(g["n_s"] - B["n_s"]) <= 0.002 and abs(g["r"] / B["r"] - 1) <= 0.10 and abs(g["P_zeta"] / g["P_sr"] - 1) <= 0.05)
    G1 = all(abs(R[n]["main"][q] / R[n]["rtol9"][q] - 1) <= 1e-4 for n in SETS for q in ("n_s", "r"))
    G2 = all(abs(R[n]["main"]["P_zeta"] / R[n]["cs1"]["P_zeta"] - 1) <= 1e-4 for n in SETS if SETS[n])
    G3 = all(R[n]["main"]["res_max"] <= 1e-10 and R[n]["main"]["z2_min"] > 0 and (SETS[n] is None or (R[n]["main"]["D_min"] > 0 and R[n]["main"]["phid_min"] > 0)) for n in SETS)
    sets = {}
    for n in SETS:
        if not SETS[n]:
            continue
        m = R[n]["main"]; e = m["eps_G"]
        sets[n] = dict(eps_G=e, lambda_minus_1=-2 * e / 3, allowed=bool(EPS_BOUND[0] <= e <= EPS_BOUND[1]),
                       dn_s=m["n_s"] - g["n_s"], dr_rel=m["r"] / g["r"] - 1, dr_rel_napkin=(1 + e) ** (2 / 3) - 1,
                       dNeff_BBN=-43 / 7 * e, children_U57b=True)
    if not (G0 and G1 and G2 and G3):
        verdict = "INCONCLUSO"
    else:
        ok = [n for n, v in sets.items() if v["allowed"] and v["children_U57b"]]
        if not ok:
            verdict = "CUSCUTON_EXCLUIDO"
        elif any(abs(sets[n]["dn_s"]) <= 0.001 for n in ok):
            verdict = "CUSCUTON_COMPATIBLE"
        else:
            verdict = "CUSCUTON_TENSO"
    out = dict(spec="docs/preregistros/ESPEC_U58.md", eps_bound_95=EPS_BOUND, chain=dict(n_s=B["n_s"], r=B["r"], N_star=B["N_star"]),
               gates=dict(G0=G0, G1=G1, G2=G2, G3=G3), runs=R, sets=sets, verdict=verdict, seconds=time.time() - t0)
    save("u58_cuscuton_inflacion.json", out)
    print(json.dumps(dict(gates=out["gates"], sets=sets, verdict=verdict, GR=dict(n_s=g["n_s"], r=g["r"], P_over_Psr=g["P_zeta"] / g["P_sr"])), indent=1, default=float))


if __name__ == "__main__":
    main()
