"""U64 (docs/preregistros/ESPEC_U64.md): the full cubic order in the cuscuton unitary gauge.  The linear solution of the U59
(uniform-radiation) gauge on the U63b vacuum (A2) is moved to the unitary gauge (delta phi = 0) with xi0 = f/phidot, and the
eps^3 / eps^2 densities of the full ADM Lagrangian of U59 are evaluated there (no 1/phidot terms: the cuscuton is
-sqrt(h) (mu^2 phidot + N U)).  U59, U62d, U63 and U63b are imported unchanged.
    python3 -B chain/u64_cubico_unitario.py   ->  results/u64_cubico_unitario.json
Internal, not exported."""
import json, math, time
from multiprocessing import Pool
import numpy as np
import sympy as sp
from ucommon import ROOT, save
import explora_cuscuton as EC
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62d_cl_rebote as D
import u63_cubico_cuscuton as Q
import u63b_gauge_unitario as G

RUNS = ROOT / "results" / "u64_runs"
NS = (10, 30, 100, 300, 1000, 3000, 10000)
NRM = Q.NRM


def build_sectors():
    """L2, L3 of the total and L3 of each sector (gravity, chi, cuscuton, radiation)"""
    src = open(ROOT / "chain" / "u59_escalares_rebote.py").read()
    body = src[src.index("def derive():"):src.index("    L2 = sp.diff(Ltot, eps, 2)")]
    ns = dict(vars(U)); exec(body.replace("def derive():", "def build():") + "    return Lg.subs(th, sp.pi / 2), Lchi, Lphi, Lrad\n", ns)
    parts = ns["build"](); Y = sp.Function("Y")(U.r)
    def coef(L, k, fac):
        L = sp.diff(L, U.eps, k).subs(U.eps, 0) / fac
        L = L.subs(sp.Derivative(Y, (U.r, 2)), U.Y2).subs(sp.Derivative(Y, U.r), U.Y1).subs(Y, U.Ys)
        return sp.lambdify(Q.ARGS, L, "numpy", cse=True)
    tot = sum(parts)
    return dict(L2=coef(tot, 2, 2), L3=coef(tot, 3, 6), sectors=[coef(p, 3, 6) for p in parts])


SECTORS = ("gravedad", "chi", "cuscuton", "radiacion")


def unitary(A, bg, n, at, t, j, c):
    """fields of the U59 gauge -> unitary gauge (linear time shift xi0 = f/phidot)"""
    b, qv, qd, u = Q.fields(A, bg, n, at, t, j, c); h = Q.step(b, n)
    up = Q.fields(A, bg, n, at, t + h, j, c)[3]; um = Q.fields(A, bg, n, at, t - h, j, c)[3]
    fd = (up[2] - um[2]) / (2 * h)
    bp, bm = bg.bg(t + h), bg.bg(t - h)
    phid = b["phid"]; phidd = (bp["phid"] - bm["phid"]) / (2 * h); Hd = (bp["H"] - bm["H"]) / (2 * h)
    H, chid, qrad = b["H"], b["chid"], b["q"]
    m_, y_, Qf = bg.state(t); chidd = -3 * H * chid - EC.dVp(y_[1]); qd_rad = -H * qrad
    al, be, f = u
    xi = f / phid; xid = fd / phid - f * phidd / phid ** 2
    zeta_u = qv[0] - H * xi; x_u = qv[1] - chid * xi; w_u = -qrad * xi
    zd_u = qd[0] - Hd * xi - H * xid; xd_u = qd[1] - chidd * xi - chid * xid; wd_u = -qd_rad * xi - qrad * xid
    al_u = al - xid; be_u = be + xi
    return b, dict(zeta=zeta_u, x=x_u, w=w_u, zd=zd_u, xd=xd_u, wd=wd_u, fd=0.0, al=al_u, be=be_u, f=0.0), dict(xi=xi, xid=xid, al_u=al_u)


ORDER = ("zeta", "x", "w", "zd", "xd", "wd", "fd", "al", "be", "f")


def measures(A, bg, n, at, t, L, npts=64):
    rr = Q.local_r(n, npts); Yv, Yr, Yrr = Q.Yprof(n, rr); best = None
    for j in range(2):
        for c in range(2):
            b, v, aux = unitary(A, bg, n, at, t, j, c)
            args = [rr, Yv, Yr, Yrr] + [b[k] for k in Q.BGS] + [v[k] for k in ORDER]
            l2 = np.mean(np.abs(L["L2"](*args) * np.ones_like(rr))); l3 = np.mean(np.abs(L["L3"](*args) * np.ones_like(rr)))
            rho3 = n ** 1.5 * float(l3 / l2)
            if best is None or rho3 > best["rho3"]:
                sec = [float(np.mean(np.abs(Ls(*args) * np.ones_like(rr)))) for Ls in L["sectors"]]
                best = dict(rho3=rho3, al_u=n ** 1.5 * abs(aux["al_u"]) * NRM, sectors=dict(zip(SECTORS, [s / max(sum(sec), 1e-300) for s in sec])))
    best.update(phid_over_H=b["phid"] / max(abs(b["H"]), 1e-300), lna=math.log(b["a"]))
    return best


def row(args):
    par, n, rtol, npts, dens = args
    t0 = time.time(); A = U.Action(); L = build_sectors(); bg = D.Background(par)
    try:
        at, info, chk = G.integrate(A, bg, n, "A2", rtol)
    except (RuntimeError, ValueError) as ex:
        return dict(par=par, n=n, status="NO_DEFINIDA", error=str(ex))
    pts = [dict(phase=ph, **measures(A, bg, n, at, t, L, npts)) for ph, t in G.times(bg, chk["t0"], dens)]
    best = max(pts, key=lambda p: p["rho3"])
    return dict(par=par, n=n, rtol=rtol, npts=npts, dens=dens, status="OK", rho3_max=best["rho3"], where=dict(phase=best["phase"], lna=best["lna"],
                phid_over_H=best["phid_over_H"], al_u=best["al_u"], sectors=best["sectors"]),
                phase_max={ph: max([p["rho3"] for p in pts if p["phase"] == ph] or [0.0]) for ph in ("contraccion", "ventana", "expansion")},
                al_u_max=max(p["al_u"] for p in pts), a_b=math.exp(bg.wl.sol(0.0)[0]),
                Hmax=max(abs(bg.bg(t)["H"]) for ph, t in G.times(bg, chk["t0"]) if ph == "ventana"), checks=chk, seconds=time.time() - t0)


def cached_row(args):
    par, n, rtol, npts, dens = args
    f = RUNS / "rows" / ("AH%g_f%.6f_n%d_rt%g_p%d_d%d.json" % (par["AH"], par["f"], n, rtol, npts, dens)); f.parent.mkdir(parents=True, exist_ok=True)
    if f.exists(): return json.loads(f.read_text())
    r = row(args); f.write_text(json.dumps(r, default=float))
    with open(RUNS / "produce.out", "a") as log:
        log.write("row AH %g f %.4f n %d rt %g p %d d %d %s %.0fs\n" % (par["AH"], par["f"], n, rtol, npts, dens, r["status"], r.get("seconds", 0)))
    return r


def gate_G1(pars, A):
    """the alpha and beta Euler-Lagrange equations of the integrated quadratic action (U59 items) on the unitary-gauge solution"""
    worst = 0.0; ia, ib = ORDER.index("al"), ORDER.index("be")
    for par in pars[:3]:
        bg = D.Background(par)
        for n in (10, 300, 10000):
            try:
                at, info, chk = G.integrate(A, bg, n, "A2", 1e-10)
            except (RuntimeError, ValueError):
                continue                                                  # a vacuum not defined at A2 (U63b): no row, no check
            Ts = [t for _, t in G.times(bg, chk["t0"])][::80]
            I = A.integrals(n)
            for t in Ts:
                b, v, _ = unitary(A, bg, n, at, t, 0, 0); vals = [v[k] for k in ORDER]; bgv = [b[k] for k in Q.BGS]
                for idx in (ia, ib):
                    s = 0.0; sc = 0.0
                    for (mon, fb, _), Ii in zip(A.items, I):
                        e = mon[idx]
                        if e == 0: continue
                        pv = e * np.prod([val ** (ee - (1 if k_ == idx else 0)) for k_, (val, ee) in enumerate(zip(vals, mon))])
                        term = float(fb(*bgv)) * Ii * pv; s += term; sc += abs(term)
                    worst = max(worst, abs(s) / max(sc, 1e-300))
    return worst


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    pars, nm, hd = Q.pick_hosts()
    L = build_sectors(); L2c = sp.lambdify(Q.ARGS, U.load_L2(), "numpy")
    rng = np.random.default_rng(64); g0 = 0.0
    for _ in range(20):
        v = [rng.uniform(0.2, 2.9)] + list(rng.normal(size=3)) + list(rng.uniform(0.5, 2, size=12)) + list(rng.normal(size=10))
        g0 = max(g0, abs(L["L2"](*v) - L2c(*v)) / max(abs(L2c(*v)), 1e-300))
    print("G0", g0, flush=True)
    A = U.Action(); g1 = gate_G1(pars, A); print("G1", g1, flush=True)
    jobs = [(par, n, 1e-10, 64, 1) for par in pars for n in NS]
    with Pool(4) as p:
        R = p.map(cached_row, jobs, chunksize=1)
    (RUNS / "rows.json").write_text(json.dumps(dict(rows=R), default=float))
    print("rows done %.0fs" % (time.time() - t0), flush=True)
    g2, g3p, g3t = {}, {}, {}
    extra = [(pars[0], n, rt, p_, d_) for n in (100, 3000, 10000) for rt, p_, d_ in ((1e-9, 64, 1), (1e-10, 128, 1), (1e-10, 64, 2))]
    with Pool(4) as p:
        X = p.map(cached_row, extra, chunksize=1)
    for n in (100, 3000, 10000):
        base = next(r for r in R if r["par"] == pars[0] and r["n"] == n)["rho3_max"]
        pick = lambda rt, p_, d_: next(x for x in X if x["n"] == n and x["rtol"] == rt and x["npts"] == p_ and x["dens"] == d_)["rho3_max"]
        g2[n] = abs(pick(1e-9, 64, 1) / base - 1); g3p[n] = abs(pick(1e-10, 128, 1) / base - 1); g3t[n] = abs(pick(1e-10, 64, 2) / base - 1)
        print("G2/G3", n, g2[n], g3p[n], g3t[n], flush=True)
    gates = dict(G0=bool(g0 <= 1e-12), G1=bool(g1 <= 1e-6), G2=bool(max(g2.values()) <= 1e-3),
                 G3=bool(max(g3p.values()) <= 1e-3 and max(g3t.values()) <= 0.05))
    hosts = []
    u63b = {(h["par"]["AH"], h["par"]["f"]): h for h in json.loads((ROOT / "results" / "u63b_gauge_unitario.json").read_text())["hosts"]}
    for par in pars:
        k_ = (par["AH"], par["f"]); h = hd[k_]
        s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))) if h.get("OmegaK0") else None; n30 = 30 / s if s else None
        rows = sorted([r for r in R if r["par"] == par and r["status"] == "OK"], key=lambda r: r["n"])
        nsc = next((r["n"] for r in rows if r["rho3_max"] >= 1), None)
        ir = max([r["rho3_max"] for r in rows if n30 and r["n"] <= n30] or [0.0])
        a_b = rows[0]["a_b"]
        hosts.append(dict(par=par, n_mem=nm[k_], n30=n30, n_sc=nsc, k_sc=(nsc / a_b) if nsc else None, k_mem=nm[k_] / a_b, Hmax=rows[0]["Hmax"],
                          DGQ=rows[0]["Hmax"] ** (2 / 3), n_sc_U63b=u63b[k_]["n_sc"], IR_max=ir,
                          per_n=[dict(n=r["n"], rho3=r["rho3_max"], al_u=r["al_u_max"], where=r["where"], phase_max=r["phase_max"]) for r in rows]))
    strong = any((H_["n_sc"] is not None and H_["n_sc"] <= H_["n_mem"]) or H_["IR_max"] >= 1 for H_ in hosts)
    pert = all((H_["n_sc"] is None or H_["n_sc"] >= 10 * H_["n_mem"]) and H_["IR_max"] <= 0.1 for H_ in hosts)
    would = "ACOPLAMIENTO_FUERTE" if strong else ("PERTURBATIVO" if pert else "LIMITE")
    verdict = would if all(gates.values()) else "INCONCLUSO"
    out = dict(spec="docs/preregistros/ESPEC_U64.md", gates=gates, G0=g0, G1=g1, G2=g2, G3_points=g3p, G3_times=g3t, hosts=hosts,
               would_be=would, verdict=verdict, seconds=time.time() - t0)
    save("u64_cubico_unitario.json", out)
    print(json.dumps(dict(gates=gates, G0=g0, G1=g1, G2=g2, G3_points=g3p, G3_times=g3t, would_be=would, verdict=verdict), indent=1, default=float))
    for H_ in hosts:
        print(H_["par"]["AH"], "%.3f" % H_["par"]["f"], "n_mem %.0f n_sc %s (U63b %s) k_sc %s IR %.2f" % (H_["n_mem"], H_["n_sc"], H_["n_sc_U63b"],
              ("%.2e" % H_["k_sc"]) if H_["k_sc"] else "-", H_["IR_max"]),
              [(p["n"], "%.1e" % p["rho3"], "%.1e" % p["al_u"], p["where"]["phase"], {k: round(v, 2) for k, v in p["where"]["sectors"].items()}) for p in H_["per_n"]], flush=True)


if __name__ == "__main__":
    main()
