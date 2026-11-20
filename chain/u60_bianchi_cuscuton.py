"""U60 (docs/preregistros/ESPEC_U60.md): Bianchi IX with Einstein + cuscuton + chi (pole plateau) + isotropic radiation.  The
cuscuton fixes the mean expansion, H = U'(phi)/(3 mu^2); phidot = 3 mu^2 W / D with W = 2/3 rho_r + chidot^2/2 + sigma^2 - R3/6.
Shape deformations q_i = ln abar + eps d_i at the mother's turnaround (matter conserved, abar from the constraint), three Misner
directions x 12 amplitudes, for the 52 hosts of U61b.  Legs and projections as u61b_cuscuton_huella.host.  explora_cuscuton is
imported unchanged (Cusc.all, Vp, dVp, H_PLATEAU).
    python3 -B chain/u60_bianchi_cuscuton.py   ->  results/u60_runs/runs.json, results/u60_bianchi_cuscuton.json
Internal, not exported."""
import json, math, sys, time
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, save
import explora_cuscuton as EC
import u61_cuscuton_huella as U61

RUNS = ROOT / "results" / "u60_runs"
ADOPTED = dict(mu=0.16, m=0.02, pinf=0.25)
S3 = math.sqrt(3.0)
DIRS = {"+b+": (1.0, 1.0, -2.0), "-b+": (-1.0, -1.0, 2.0), "+b-": (S3, -S3, 0.0)}
EPS = tuple(float(x) for x in np.logspace(-6, -0.5, 12))
ZETA0 = 4.6e-5


def ricci(q):
    c = [2 * math.exp(q[i] - q[(i + 1) % 3] - q[(i + 2) % 3]) for i in range(3)]
    R = [0.5 * (c[i] ** 2 - (c[(i + 1) % 3] - c[(i + 2) % 3]) ** 2) for i in range(3)]
    return R, sum(R)


class Model:
    def __init__(self, par):
        self.cu = EC.Cusc(par["mu"], par["m"], par["pinf"]); self.par = par
        a_t = par["AH"] / EC.H_PLATEAU; self.a_t = a_t; self.Crad = 3 * a_t ** 2 / (1 + par["f"])

    def Q(self, mode, y):
        q = y[0:3]; s1, s2 = y[3], y[4]; s = (s1, s2, -s1 - s2); chi, chid, d = y[5], y[6], y[7]
        cu = self.cu; alpha = sum(q) / 3; rho = self.Crad * math.exp(-4 * alpha)
        phi, U, dU, ddU = cu.all(mode, d)
        R, R3 = ricci(q); sig2 = 0.5 * sum(x * x for x in s)
        H = dU / (3 * cu.mu2); D = 1.5 * cu.mu4 - ddU
        W = 2 / 3 * rho + 0.5 * chid * chid + sig2 - R3 / 6
        phid = 3 * cu.mu2 * W / D; Hd = ddU * phid / (3 * cu.mu2)
        rhs_ = rho + 0.5 * chid * chid + EC.Vp(chi) + sig2 - R3 / 2
        res = abs(3 * H * H - U - rhs_) / (rho + 0.5 * chid * chid + EC.Vp(chi) + abs(U) + sig2 + abs(R3) / 2)
        return dict(alpha=alpha, rho=rho, phi=phi, U=U, H=H, Hdot=Hd, W=W, D=D, phid=phid, R=R, R3=R3, sig2=sig2, s=s, res=res,
                    Wscale=2 / 3 * rho + 0.5 * chid * chid + sig2 + abs(R3) / 6, rhs_=rhs_)

    def rhs(self, mode):
        def f(t, y):
            Q = self.Q(mode, y); H = Q["H"]; R, R3 = Q["R"], Q["R3"]; s = Q["s"]
            return [H + s[0], H + s[1], H + s[2],
                    -3 * H * s[0] - (R[0] - R3 / 3), -3 * H * s[1] - (R[1] - R3 / 3),
                    y[6], -3 * H * y[6] - EC.dVp(y[5]), (1 if mode <= 0 else -1) * Q["phid"]]
        return f

    def project(self, y):
        """mode +1: reset d so that the Hamiltonian constraint holds (as u61b.project, with shear and anisotropic curvature)"""
        Q = self.Q(1, y); rhs_ = Q["rhs_"]; cu = self.cu
        g = lambda d: (lambda r: 3 * (r[2] / (3 * cu.mu2)) ** 2 - r[1])(cu.all(1, d)) - rhs_
        lo, hi = y[7] * 0.9, y[7] * 1.1
        if g(lo) * g(hi) < 0:
            y = list(y); y[7] = brentq(g, lo, hi, xtol=1e-300, rtol=1e-15)
        return y


def slice_state(model, eps, dvec):
    q0 = [eps * x for x in dvec]; _, r3 = ricci(q0)
    if r3 <= 0: return None
    abar = model.a_t * math.sqrt(6 / r3); rho = model.Crad / abar ** 4
    return [math.log(abar) + x for x in q0] + [0.0, 0.0, 0.0, math.sqrt(2 * model.par["f"] * rho), 0.0]


ATOL = [1e-13] * 3 + [1e-24, 1e-24, 1e-14, 1e-24, 1e-24]


def run(par, eps, dname, rtol=1e-12, keep=False):
    M = Model(par); cu = M.cu; half = cu.pinf / 2
    y = slice_state(M, eps, DIRS[dname] if dname else (0, 0, 0))
    if y is None: return dict(par=par, eps=eps, dir=dname, outcome="SIN_RODAJA")
    legs = []; res = [0.0]; t = 0.0; a_b_iso = None
    def leg(mode, y, t, events, name):
        f = M.rhs(mode); evs = []
        for kind, direction, term in events:
            fn = {"acc": lambda tt, yy: (lambda q: q["Hdot"] + q["H"] ** 2)(M.Q(mode, yy)),
                  "d_half": lambda tt, yy: yy[7] - half, "d0": lambda tt, yy: yy[7], "turn": lambda tt, yy: yy[7],
                  "chi_end": lambda tt, yy: abs(yy[5]) - U61.XEND}[kind]
            g = (lambda h: (lambda tt, yy: h(tt, yy)))(fn); g.direction = direction; g.terminal = term; evs.append(g)
        wev = lambda tt, yy: (lambda q: q["W"] / q["Wscale"])(M.Q(mode, yy)) - 1e-12; wev.terminal = True; wev.direction = -1
        aev = lambda tt, yy: 10 - (max(yy[0:3]) - min(yy[0:3])); aev.terminal = True; aev.direction = -1
        evs += [wev, aev]
        s = solve_ivp(f, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=ATOL, events=evs, dense_output=True)
        res[0] = max(res[0], max(M.Q(mode, s.y[:, i])["res"] for i in range(s.y.shape[1])))
        legs.append((name, mode, s))
        if len(s.t_events[-2]): return s, "RAMA_ROTA"
        if len(s.t_events[-1]): return s, "ANISOTROPO"
        if s.status != 1: return s, "SIN_EVENTO"
        hit = [kind for (kind, _, term), te in zip(events, s.t_events[:-2]) if term and len(te)]
        return s, (hit[0] if hit else None)
    def fail(why):
        return dict(par=par, eps=eps, dir=dname, outcome=why, constraint_max=res[0])
    s, hit = leg(-1, y, t, [("d_half", 1, True)], "half-")
    if hit != "d_half": return fail(hit)
    t = s.t[-1]; y = list(s.y[:, -1]); y[7] = -half
    s, hit = leg(0, y, t, [("d0", 1, False), ("d_half", 1, True)], "half+")
    if hit != "d_half": return fail(hit)
    tb = float(s.t_events[0][0]); yb = s.sol(tb); Qb = M.Q(0, yb)
    bounce = dict(t=tb, a=math.exp(Qb["alpha"]), sig2_over_rho=Qb["sig2"] / Qb["rho"], aniso=max(yb[0:3]) - min(yb[0:3]),
                  chid=yb[6], W_rel=Qb["W"] / Qb["Wscale"])
    t = s.t[-1]; y = list(s.y[:, -1]); y[7] = half
    s, hit = leg(1, y, t, [("acc", -1, True), ("turn", -1, True)], "decel")
    if hit != "acc": return fail("RECOLAPSO" if hit == "turn" else hit)
    t = s.t[-1]; y = M.project(list(s.y[:, -1]))
    s, hit = leg(1, y, t, [("acc", 1, True), ("turn", -1, True)], "onset")
    if hit != "acc": return fail("RECOLAPSO" if hit == "turn" else hit)
    t_on = s.t[-1]; al_on = sum(s.y[0:3, -1]) / 3; y = M.project(list(s.y[:, -1]))
    s, hit = leg(1, y, t_on, [("acc", -1, True)], "end")
    if hit != "acc": return fail(hit)
    al_end = sum(s.y[0:3, -1]) / 3; t = s.t[-1]; y = M.project(list(s.y[:, -1]))
    s, hit = leg(1, y, t, [("chi_end", -1, True)], "tail")
    if hit != "chi_end": return fail(hit)
    # landing: first chidot = 0 with |chi| > 1 after the bounce
    land = None
    for name, mode, sl in legs[1:]:
        lo = max(sl.t[0], tb); hi = sl.t[-1]
        if hi <= lo: continue
        ts = np.linspace(lo, hi, 4001); Y = sl.sol(ts)
        for i in range(len(ts) - 1):
            if abs(Y[5, i]) > 1 and Y[6, i] * Y[6, i + 1] < 0:
                tl = brentq(lambda tt: sl.sol(tt)[6], ts[i], ts[i + 1]); land = abs(sl.sol(tl)[5]); break
        if land is not None: break
    if land is None: return fail("SIN_ATERRIZAJE")
    Qe = M.Q(1, s.y[:, -1])
    out = dict(par=par, eps=eps, dir=dname, outcome="HIJO", chi_land=land, N_acc=float(al_end - al_on),
               Omega_sigma_end=Qe["sig2"] / (3 * Qe["H"] ** 2), constraint_max=res[0], bounce=bounce)
    if keep: out["_legs"] = legs; out["_model"] = M
    return out


def job(a):
    try:
        return run(*a)
    except Exception as ex:
        return dict(par=a[0], eps=a[1], dir=a[2], outcome="ERROR", error=repr(ex))


def gate_B1(hosts):
    """eps = 1e-7 along +b-: delta = (q1 - q2)/2 against the linear equation on the isotropic run, turnaround -> bounce"""
    out = []
    for h in hosts[:3]:
        par = h["par"]; r0 = run(par, 0.0, None, keep=True); r1 = run(par, 1e-7, "+b-", keep=True)
        # isotropic background functions on [0, t_b]: H(t), a(t)
        segs = [(m, s) for _, m, s in r0["_legs"][:2]]; M0 = r0["_model"]; tb = r0["bounce"]["t"]
        def bg(t):
            for m, s in segs:
                if s.t[0] <= t <= s.t[-1]: y = s.sol(t); Q = M0.Q(m, y); return Q["H"], math.exp(Q["alpha"])
            m, s = segs[-1]; y = s.sol(t); Q = M0.Q(m, y); return Q["H"], math.exp(Q["alpha"])
        d0 = S3 * 1e-7
        lin = solve_ivp(lambda t, y: (lambda Ha: [y[1], -3 * Ha[0] * y[1] - 8 * y[0] / Ha[1] ** 2])(bg(t)), (0.0, tb), [d0, 0.0],
                        method="DOP853", rtol=1e-12, atol=1e-30, dense_output=True)
        segs1 = [s for _, m, s in r1["_legs"][:2]]
        ts = np.linspace(0, min(tb, r1["bounce"]["t"]), 2001); dn = []
        for t in ts:
            for s in segs1:
                if s.t[0] <= t <= s.t[-1]: y = s.sol(t); dn.append(0.5 * (y[0] - y[1])); break
            else: dn.append(np.nan)
        dl = lin.sol(ts)[0]; dn = np.array(dn)
        out.append(dict(AH=par["AH"], f=par["f"], rel=float(np.nanmax(np.abs(dn - dl)) / np.max(np.abs(dl))),
                        growth=float(np.max(np.abs(dl)) / d0)))
    return out


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    u61 = json.loads((ROOT / "results" / "u61b_cuscuton_huella.json").read_text())
    hosts = [h for h in u61["hosts"] if "chi_land" in h]
    print(len(hosts), "hosts", flush=True)
    jobs = [(dict(ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]), 0.0, None) for h in hosts]
    jobs += [(dict(ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]), e, d) for h in hosts for d in DIRS for e in EPS]
    path = RUNS / "runs.json"
    if path.exists():
        R = json.loads(path.read_text())["runs"]
    else:
        with Pool(4) as p:
            R = p.map(job, jobs, chunksize=1)
        path.write_text(json.dumps(dict(runs=R), default=float))
    print("runs done %.0fs" % (time.time() - t0), flush=True)
    key = lambda r: (r["par"]["AH"], r["par"]["f"])
    iso = {key(r): r for r in R if r["dir"] is None}
    # B0
    b0 = []
    for h in hosts:
        r = iso[(h["par"]["AH"], h["par"]["f"])]
        b0.append(dict(AH=h["par"]["AH"], f=h["par"]["f"], outcome=r["outcome"],
                       chi_rel=abs(r["chi_land"] / h["chi_land"] - 1) if r.get("chi_land") else None,
                       dN=abs(r["N_acc"] - h["N_acc"]) if r.get("N_acc") is not None else None))
    B0 = all(x["chi_rel"] is not None and x["chi_rel"] <= 1e-5 and x["dN"] <= 1e-4 for x in b0)
    print("B0", B0, max((x["chi_rel"] or 9) for x in b0), max((x["dN"] or 9) for x in b0), flush=True)
    b1 = gate_B1([dict(par=dict(ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"])) for h in hosts])
    B1 = all(x["rel"] <= 1e-4 for x in b1); print("B1", B1, b1, flush=True)
    # B3: two hosts (lowest and highest A.H), eps = 1e-3 and 1e-1.5, directions +b+ and +b-
    hs = [min(hosts, key=lambda h: h["par"]["AH"]), max(hosts, key=lambda h: h["par"]["AH"])]
    b3 = []
    for h in hs:
        par = dict(ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"])
        for e in (EPS[5], EPS[9]):
            for d in ("+b+", "+b-"):
                a = run(par, e, d, rtol=1e-11); b = run(par, e, d, rtol=1e-12)
                b3.append(dict(AH=par["AH"], eps=e, dir=d, same=a["outcome"] == b["outcome"],
                               dN=abs(a.get("N_acc", 0) - b.get("N_acc", 0)), dO=abs(a.get("Omega_sigma_end", 0) - b.get("Omega_sigma_end", 0))))
    B3 = all(x["same"] and x["dN"] <= 1e-3 and x["dO"] <= 1e-4 for x in b3); print("B3", B3, b3, flush=True)
    # robustness and eps*
    def robust(r, r0):
        return (r["outcome"] == "HIJO" and abs(r["N_acc"] - r0["N_acc"]) <= 1 and r["Omega_sigma_end"] < 1e-2 and r["constraint_max"] <= 1e-7)
    valid = lambda r: r.get("constraint_max", 0) <= 1e-7
    small = [r for r in R if r["dir"] and r["eps"] <= 1e-2 and r["outcome"] != "SIN_RODAJA"]
    inv_frac = sum(not valid(r) for r in small) / max(1, len(small))
    children = []
    for h in hosts:
        k_ = (h["par"]["AH"], h["par"]["f"]); r0 = iso[k_]; es = 0.0; first_fail = None
        for e in EPS:
            rs = [r for r in R if r["dir"] and key(r) == k_ and r["eps"] == e]
            if all(r["outcome"] != "SIN_RODAJA" and robust(r, r0) for r in rs): es = e
            else:
                first_fail = [dict(dir=r["dir"], outcome=r["outcome"], dN=(r["N_acc"] - r0["N_acc"]) if r.get("N_acc") is not None else None,
                                   Om=r.get("Omega_sigma_end"), res=r.get("constraint_max")) for r in rs]; break
        bn = {d: next((r["bounce"]["sig2_over_rho"] for r in R if r["dir"] == d and key(r) == k_ and r["eps"] == EPS[4] and r.get("bounce")), None) for d in DIRS}
        children.append(dict(AH=k_[0], f=k_[1], N_acc_iso=r0.get("N_acc"), eps_star=es, robust=es >= ZETA0, first_fail=first_fail,
                             sig2_over_rho_bounce_eps1e4=bn))
    f_rob = sum(c["robust"] for c in children) / len(children)
    gates = dict(B0=B0, B1=B1, B2=inv_frac <= 0.2, B3=B3)
    verdict = "INCONCLUSO" if not all(gates.values()) else ("ANISOTROPIA_INOCUA" if f_rob >= 0.9 else ("ANISOTROPIA_FATAL" if f_rob <= 0.1 else "ANISOTROPIA_ACOTADA"))
    would = "ANISOTROPIA_INOCUA" if f_rob >= 0.9 else ("ANISOTROPIA_FATAL" if f_rob <= 0.1 else "ANISOTROPIA_ACOTADA")
    # informative: initial shear velocity (s_i != 0) on the first host is not run here (Sec. 7 lists it as an extra row)
    out = dict(spec="docs/preregistros/ESPEC_U60.md", gates=gates, B0=b0, B1=b1, B3=b3, B2_invalid_fraction=inv_frac, children=children,
               f_rob=f_rob, would_be=would, verdict=verdict, seconds=time.time() - t0)
    save("u60_bianchi_cuscuton.json", out)
    es = [c["eps_star"] for c in children]
    print(json.dumps(dict(gates=gates, B2_invalid_fraction=inv_frac, f_rob=f_rob, verdict=verdict, would_be=would,
                          eps_star_min=min(es), eps_star_median=float(np.median(es)), eps_star_max=max(es)), indent=1), flush=True)


if __name__ == "__main__":
    main()
