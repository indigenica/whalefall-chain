"""U63 (docs/preregistros/ESPEC_U63.md, amendment 1): cubic order and strong coupling of the cuscuton near phidot -> 0.  On the
linear vacuum modes (B_s-free frame, set at A1 near the mother's turnaround, background of U62d integrated forward), from A1 to
rho_r/rho = 1e-5 in early inflation: rho3 = <|L3|>/<|L2|> of the full ADM Lagrangian density of U59 (eps^3 against eps^2, local
wavelength average at r ~ pi/2), Q_u (perturbation of 2X on constant-phi slices) and Q_t (tilt of constant-phi slices), with log-band
vacuum amplitudes.  U59, U62b, U62c and U62d are imported unchanged.
    python3 -B chain/u63_cubico_cuscuton.py   ->  results/u63_cubico_cuscuton.json
Internal, not exported."""
import json, math, pickle, sys, time
from multiprocessing import Pool
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from ucommon import ROOT, save
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62b_bogoliubov_fisico as B62
import u62c_espectro_rebote as C
import u62d_cl_rebote as D

RUNS = ROOT / "results" / "u63_runs"
NS = (2, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
BGS = ("a", "H", "chid", "phid", "q", "mu2", "Vv", "Vp", "Vpp", "Uv", "Up", "Upp")
ARGS = (U.r, U.Ys, U.Y1, U.Y2, U.a, U.H, U.chid, U.phid, U.q, U.mu2, U.Vv, U.Vp, U.Vpp, U.Uv, U.Up, U.Upp) + U.QV + U.VV + (U.fd,) + U.UV
NRM = 1 / (math.pi * math.sqrt(2))


def build_L23():
    src = open(ROOT / "chain" / "u59_escalares_rebote.py").read()
    body = src[src.index("def derive():"):src.index("    L2 = sp.diff(Ltot, eps, 2)")]
    ns = dict(vars(U)); exec(body.replace("def derive():", "def build():") + "    return Ltot\n", ns)
    Ltot = ns["build"](); Y = sp.Function("Y")(U.r)
    out = []
    for k, fac in ((2, 2), (3, 6)):
        L = sp.diff(Ltot, U.eps, k).subs(U.eps, 0) / fac
        L = L.subs(sp.Derivative(Y, (U.r, 2)), U.Y2).subs(sp.Derivative(Y, U.r), U.Y1).subs(Y, U.Ys)
        out.append(sp.lambdify(ARGS, L, "numpy", cse=True))
    return out


def Yprof(n, rr):
    k = n + 1; s, c = np.sin(rr), np.cos(rr)
    Y = NRM * np.sin(k * rr) / s; Yr = NRM * (k * np.cos(k * rr) * s - np.sin(k * rr) * c) / s ** 2
    return Y, Yr, -n * (n + 2) * Y - 2 * c / s * Yr


def local_r(n, npts):
    half = math.pi / (n + 1); return np.linspace(math.pi / 2 - half, math.pi / 2 + half, npts)


def integrate_modes(A, bg, n, rtol):
    """the two in-modes (vacuum at A1, B_s-free frame) as dense solutions of (q, p = a^3 P), per leg"""
    t0, Z, info = C.in_modes(A, bg, n, bg.t["A1"], 1)
    a30 = bg.bg(t0)["a"] ** 3
    Y0 = np.concatenate([np.concatenate([Z[:2, j].real, Z[:2, j].imag, a30 * Z[2:, j].real, a30 * Z[2:, j].imag]) for j in range(2)])
    def rhs(t, y):
        b = bg.bg(t); K, Bm, Cm, _ = U.red2(A, b, n); a3 = b["a"] ** 3; out = np.empty_like(y)
        for j in range(2):
            o = 8 * j
            for c in range(2):                                      # real, imaginary parts
                qv = y[o + 2 * c:o + 2 * c + 2]; p = y[o + 4 + 2 * c:o + 4 + 2 * c + 2]
                qd = np.linalg.solve(K, p / a3 - Bm @ qv)
                out[o + 2 * c:o + 2 * c + 2] = qd; out[o + 4 + 2 * c:o + 4 + 2 * c + 2] = a3 * (Bm.T @ qd + Cm @ qv)
        return out
    sols = []; y = Y0
    for a_, b_ in ((t0, bg.tm), (bg.tm, bg.tp), (bg.tp, bg.t["S5"])):
        s = solve_ivp(rhs, (a_, b_), y, method="DOP853", rtol=rtol, atol=1e-300, dense_output=True); sols.append(s); y = s.y[:, -1]
    def at(t):
        for s in sols:
            if s.t[0] <= t <= s.t[-1]: return s.sol(t)
        return sols[-1].sol(t)
    return at, info


def fields(A, bg, n, at, t, j, c):
    """(bg dict, q, qdot, u = (alpha, beta, f)) for mode j, quadrature c, at time t"""
    b = bg.bg(t); m = A.matrices(b, n); K, Bm, _, _ = U.red2(A, b, n); y = at(t); o = 8 * j
    qv = y[o + 2 * c:o + 2 * c + 2]; P = y[o + 4 + 2 * c:o + 4 + 2 * c + 2] / b["a"] ** 3
    qd = np.linalg.solve(K, P - Bm @ qv)
    q3 = np.array([qv[0], qv[1], 0.0]); qd3 = np.array([qd[0], qd[1], 0.0])
    u = -np.linalg.solve(m["S"], m["T"].T @ qd3 + m["R"] @ q3)
    return b, qv, qd, u


def step(b, n):
    k = math.sqrt(n * (n + 2.0)) / b["a"]
    return 1e-4 * min(1 / k, 1 / max(abs(b["H"]), 1e-300), 46.0)


def measures(A, bg, n, at, t, L2f, L3f, npts=64):
    out = dict(rho3=0.0, Qu=0.0, Qt=0.0)
    rr = local_r(n, npts); Yv, Yr, Yrr = Yprof(n, rr)
    for j in range(2):
        for c in range(2):
            b, qv, qd, u = fields(A, bg, n, at, t, j, c); h = step(b, n)
            bp, _, _, up = fields(A, bg, n, at, t + h, j, c); bm, _, _, um = fields(A, bg, n, at, t - h, j, c)
            fd = (up[2] - um[2]) / (2 * h); phid = b["phid"]; phidd = (bp["phid"] - bm["phid"]) / (2 * h)
            al, be, f = u
            bgv = [b[k_] for k_ in BGS]
            args = [rr, Yv, Yr, Yrr] + bgv + [qv[0], qv[1], 0.0, qd[0], qd[1], 0.0, fd, al, be, f]
            l2 = np.abs(L2f(*args) * np.ones_like(rr)); l3 = np.abs(L3f(*args) * np.ones_like(rr))
            n15 = n ** 1.5
            rho3 = n15 * float(np.mean(l3) / np.mean(l2))
            Qu = n15 * abs(2 * fd / phid - 2 * al - 2 * f * phidd / phid ** 2) * NRM
            Qt = n15 * math.sqrt(n * (n + 2.0)) * abs(f) / (b["a"] * phid) * NRM
            out["rho3"] = max(out["rho3"], rho3); out["Qu"] = max(out["Qu"], Qu); out["Qt"] = max(out["Qt"], Qt)
    b = bg.bg(t); out.update(t=t, lna=math.log(b["a"]), phid_over_H=b["phid"] / max(abs(b["H"]), 1e-300), H=b["H"], phid=b["phid"])
    return out


def times(bg):
    """120 in the contraction (uniform in ln a), 120 in the window (uniform in t), 160 in the expansion (uniform in ln a)"""
    from scipy.optimize import brentq
    T = []
    s = bg.cl; l0, l1 = s.sol(bg.t["A1"])[0], s.sol(bg.tm)[0]
    for l in np.linspace(l0, l1, 121)[:-1]:
        T.append(("contraccion", brentq(lambda tt: s.sol(tt)[0] - l, bg.t["A1"], bg.tm, xtol=1e-12 * abs(bg.t["A1"]), maxiter=500) if l != l0 else bg.t["A1"] * (1 - 1e-12)))
    T += [("ventana", float(t)) for t in np.linspace(bg.tm, bg.tp, 122)[1:-1]]
    s = bg.el; l0, l1 = s.sol(bg.tp)[0], s.sol(bg.t["S5"])[0]
    for l in np.linspace(l0, l1, 161)[1:]:
        T.append(("expansion", brentq(lambda tt: s.sol(tt)[0] - l, bg.tp, bg.t["S5"], xtol=1e-12 * abs(bg.t["S5"]), maxiter=500) if l != l1 else bg.t["S5"] * (1 - 1e-12)))
    return T


def row(args):
    par, n, rtol, npts = args
    t0 = time.time(); A = U.Action(); L2f, L3f = build_L23(); bg = D.Background(par)
    at, info = integrate_modes(A, bg, n, rtol)
    pts = [dict(phase=ph, **measures(A, bg, n, at, t, L2f, L3f, npts)) for ph, t in times(bg) if t >= info["t"]]
    best = {k: max(pts, key=lambda p: p[k]) for k in ("rho3", "Qu", "Qt")}
    phase_max = {ph: {k: max(p[k] for p in pts if p["phase"] == ph) for k in ("rho3", "Qu", "Qt")} for ph in ("contraccion", "ventana", "expansion")}
    return dict(par=par, n=n, rtol=rtol, npts=npts, A1_state=dict(shifts=info["shifts"], adiab=info["adiab"]),
                max=dict(rho3=best["rho3"]["rho3"], Qu=best["Qu"]["Qu"], Qt=best["Qt"]["Qt"]),
                where={k: dict(phase=best[k]["phase"], lna=best[k]["lna"], phid_over_H=best[k]["phid_over_H"]) for k in best},
                phase_max=phase_max, a_b=math.exp(bg.wl.sol(0.0)[0]), Hmax=max(abs(p["H"]) for p in pts if p["phase"] == "ventana"),
                phid_over_H_end=pts[-1]["phid_over_H"], phid_over_H_min=min(p["phid_over_H"] for p in pts if p["phase"] == "expansion"),
                trend=[(p["phase"], p["lna"], p["phid_over_H"], p["rho3"], p["Qu"], p["Qt"]) for p in pts[::4]], seconds=time.time() - t0)


def pick_hosts():
    d = json.loads((ROOT / "results" / "u61b_cuscuton_huella.json").read_text()); H = [h for h in d["hosts"] if "chi_land" in h]
    vis = M.visible_hosts(); sel = [h["par"] for h in vis]
    for A_, mode in ((48.0, min), (300.0, min), (100.0, max)):
        c = [h for h in H if h["par"]["AH"] == A_]; sel.append(mode(c, key=lambda h: h["par"]["f"])["par"])
    nm = {(h["par"]["AH"], h["par"]["f"]): h["n_mem"] for h in H}
    return [dict(U.ADOPTED, AH=p["AH"], f=p["f"]) for p in sel], nm, {(h["par"]["AH"], h["par"]["f"]): h for h in H}


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    pars, nm, hd = pick_hosts()
    # G0: L2 of this code against the cached U59 L2
    L2f, L3f = build_L23(); L2c = sp.lambdify(ARGS, U.load_L2(), "numpy")
    rng = np.random.default_rng(63); g0 = 0.0
    for _ in range(20):
        v = [rng.uniform(0.2, 2.9)] + list(rng.normal(size=3)) + list(rng.uniform(0.5, 2, size=12)) + list(rng.normal(size=10))
        x1, x2 = L2f(*v), L2c(*v); g0 = max(g0, abs(x1 - x2) / max(abs(x2), 1e-300))
    print("G0", g0, flush=True)
    # G1: int L2 dr (full composite GL, 4 pi) against the U59 items, 3 hosts x n = 2, 10, 100 x 5 times
    A = U.Action(); g1 = 0.0
    from numpy.polynomial.legendre import leggauss
    for par in pars[:3]:
        bg = D.Background(par)
        for n in (2, 10, 100):
            at, inf_ = integrate_modes(A, bg, n, 1e-10); Ts = [t for _, t in times(bg) if t >= inf_["t"]][::80]
            npts = max(400, 8 * (n + 2)); panels = max(20, npts // 16); xg, wg = leggauss(16)
            edges = np.linspace(0, math.pi, panels + 1); hw = 0.5 * (edges[1:] - edges[:-1]); mid = 0.5 * (edges[1:] + edges[:-1])
            rr = (mid[:, None] + hw[:, None] * xg[None, :]).ravel(); ww = (hw[:, None] * wg[None, :]).ravel(); Yv, Yr, Yrr = Yprof(n, rr)
            I = A.integrals(n)
            for t in Ts:
                b, qv, qd, u = fields(A, bg, n, at, t, 0, 0); h = step(b, n)
                fd = (fields(A, bg, n, at, t + h, 0, 0)[3][2] - fields(A, bg, n, at, t - h, 0, 0)[3][2]) / (2 * h)
                vals = dict(zip([str(s_) for s_ in A.pert], [qv[0], qv[1], 0.0, qd[0], qd[1], 0.0, fd, u[0], u[1], u[2]]))
                bgv = [b[k_] for k_ in BGS]
                direct = 4 * math.pi * float(np.sum(ww * L2f(rr, Yv, Yr, Yrr, *bgv, qv[0], qv[1], 0.0, qd[0], qd[1], 0.0, fd, *u)))
                items = 0.0
                for (mon, fb, _), Ii in zip(A.items, I):
                    pv = 1.0
                    for s_, e_ in zip(A.pert, mon): pv *= vals[str(s_)] ** e_
                    items += float(fb(*bgv)) * Ii * pv
                g1 = max(g1, abs(direct - items) / max(abs(items), 1e-300))
    print("G1", g1, flush=True)
    jobs = [(par, n, 1e-10, 64) for par in pars for n in NS]
    path = RUNS / "rows.json"
    if path.exists():
        R = json.loads(path.read_text())["rows"]
    else:
        with Pool(4) as p:
            R = p.map(row, jobs, chunksize=1)
        path.write_text(json.dumps(dict(rows=R), default=float))
    print("rows done %.0fs" % (time.time() - t0), flush=True)
    # G2 (rtol) and G3 (64 -> 128 points) on the first host, n = 100, 3000, 30000
    g2 = {}; g3 = {}
    for n in (100, 3000, 30000):
        base = next(r for r in R if r["par"] == pars[0] and r["n"] == n)
        mx = lambda r: max(r["max"].values())
        a9 = row((pars[0], n, 1e-9, 64)); a128 = row((pars[0], n, 1e-10, 128))
        g2[n] = abs(mx(a9) / mx(base) - 1); g3[n] = abs(a128["max"]["rho3"] / base["max"]["rho3"] - 1)
        print("G2/G3", n, g2[n], g3[n], flush=True)
    gates = dict(G0=bool(g0 <= 1e-12), G1=bool(g1 <= 1e-8), G2=bool(max(g2.values()) <= 1e-3), G3=bool(max(g3.values()) <= 1e-3))
    hosts = []
    for par in pars:
        rows = sorted([r for r in R if r["par"] == par], key=lambda r: r["n"]); k_ = (par["AH"], par["f"]); h = hd[k_]
        s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))) if h.get("OmegaK0") else None; n30 = 30 / s if s else None
        mx = [max(r["max"].values()) for r in rows]
        nsc = next((r["n"] for r, m in zip(rows, mx) if m >= 1), None)
        ir = max([m for r, m in zip(rows, mx) if n30 and r["n"] <= n30] or [0.0])
        a_b = rows[0]["a_b"]; Hmax = rows[0]["Hmax"]
        hosts.append(dict(par=par, n_mem=nm[k_], n30=n30, n_sc=nsc, k_sc=(nsc / a_b) if nsc else None, k_mem=nm[k_] / a_b, Hmax=Hmax,
                          DGQ=Hmax ** (2 / 3), a_b=a_b, IR_max=ir, per_n=[dict(n=r["n"], max=r["max"], where=r["where"], phase_max=r["phase_max"]) for r in rows],
                          phid_over_H_min=min(r["phid_over_H_min"] for r in rows)))
    strong = any((H_["n_sc"] is not None and H_["n_sc"] <= H_["n_mem"]) or H_["IR_max"] >= 1 for H_ in hosts)
    pert = all((H_["n_sc"] is None or H_["n_sc"] >= 10 * H_["n_mem"]) and H_["IR_max"] <= 0.1 for H_ in hosts)
    would = "ACOPLAMIENTO_FUERTE" if strong else ("PERTURBATIVO" if pert else "LIMITE")
    verdict = would if all(gates.values()) else "INCONCLUSO"
    out = dict(spec="docs/preregistros/ESPEC_U63.md", gates=gates, G0=g0, G1=g1, G2=g2, G3=g3, hosts=hosts, would_be=would, verdict=verdict,
               seconds=time.time() - t0)
    save("u63_cubico_cuscuton.json", out)
    print(json.dumps(dict(gates=gates, G0=g0, G1=g1, G2=g2, G3=g3, would_be=would, verdict=verdict), indent=1, default=float))
    for H_ in hosts:
        print(H_["par"]["AH"], "%.3f" % H_["par"]["f"], "n_mem %.0f n_sc %s k_sc %s k_mem %.2e Hmax %.2e DGQ %.2e IR %.2e phid/H_min %.1e" %
              (H_["n_mem"], H_["n_sc"], ("%.2e" % H_["k_sc"]) if H_["k_sc"] else "-", H_["k_mem"], H_["Hmax"], H_["DGQ"], H_["IR_max"], H_["phid_over_H_min"]),
              [(p["n"], "%.1e/%.1e/%.1e" % (p["max"]["rho3"], p["max"]["Qu"], p["max"]["Qt"]), p["where"]["Qt"]["phase"]) for p in H_["per_n"]], flush=True)


if __name__ == "__main__":
    main()
