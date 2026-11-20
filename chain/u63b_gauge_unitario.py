"""U63b (docs/preregistros/ESPEC_U63b.md): strong coupling of the cuscuton in unitary-gauge variables.  On vacuum modes (B_s-free
frame) set at A2 (main path) or A1 (informative), propagated with the U59 matrices on the U62d background up to rho_r/rho = 1e-5:
M1 = |zeta - H f/phidot| (curvature on constant-phi slices), M2 = k |f|/phidot (tilt of constant-phi slices), M3 = |f/phidot| |phiddot/
phidot| (nonlinearity of the time shift to the unitary gauge), with log-band vacuum amplitudes.  U59, U62b, U62c, U62d and U63 are
imported unchanged.
    python3 -B chain/u63b_gauge_unitario.py   ->  results/u63b_gauge_unitario.json
Internal, not exported."""
import json, math, time
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, save
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62b_bogoliubov_fisico as B62
import u62c_espectro_rebote as C
import u62d_cl_rebote as D
import u63_cubico_cuscuton as Q

RUNS = ROOT / "results" / "u63b_runs"
NS = Q.NS
NRM = Q.NRM


def integrate(A, bg, n, key, rtol):
    t0, Z, info = C.in_modes(A, bg, n, bg.t[key], 1)
    a30 = bg.bg(t0)["a"] ** 3
    Y0 = np.concatenate([np.concatenate([Z[:2, j].real, Z[:2, j].imag, a30 * Z[2:, j].real, a30 * Z[2:, j].imag]) for j in range(2)])
    def rhs(t, y):
        b = bg.bg(t); K, Bm, Cm, _ = U.red2(A, b, n); a3 = b["a"] ** 3; out = np.empty_like(y)
        for j in range(2):
            o = 8 * j
            for c in range(2):
                qv = y[o + 2 * c:o + 2 * c + 2]; p = y[o + 4 + 2 * c:o + 4 + 2 * c + 2]
                qd = np.linalg.solve(K, p / a3 - Bm @ qv)
                out[o + 2 * c:o + 2 * c + 2] = qd; out[o + 4 + 2 * c:o + 4 + 2 * c + 2] = a3 * (Bm.T @ qd + Cm @ qv)
        return out
    segs = [(a_, b_) for a_, b_ in ((t0, bg.tm), (bg.tm, bg.tp), (bg.tp, bg.t["S5"])) if b_ > a_]
    sols = []; y = Y0
    for a_, b_ in segs:
        s = solve_ivp(rhs, (a_, b_), y, method="DOP853", rtol=rtol, atol=1e-300, dense_output=True); sols.append(s); y = s.y[:, -1]
    def at(t):
        for s in sols:
            if s.t[0] <= t <= s.t[-1]: return s.sol(t)
        return sols[-1].sol(t)
    def wr(y):
        out = []
        for j in range(2):
            o = 8 * j; q = y[o:o + 2] + 1j * y[o + 2:o + 4]; p = y[o + 4:o + 6] + 1j * y[o + 6:o + 8]
            out.append(float(2 * np.imag(q @ np.conj(p))))
        return out
    g0 = max(max(abs(M.Om(info_F[:, j], info_F[:, k].conj()) - (1j if j == k else 0)) for j in range(2) for k in range(2))
             for info_F in [B62.end_info(A, bg, n, t0, 1)["F"]])
    return at, info, dict(t0=t0, G0=float(g0), W_start=wr(Y0), W_end=wr(y))


def times(bg, t0, dens=1):
    T = []
    if t0 < bg.tm:
        s = bg.cl; l0, l1 = s.sol(t0)[0], s.sol(bg.tm)[0]
        for l in np.linspace(l0, l1, 120 * dens + 1)[1:-1]:
            T.append(("contraccion", brentq(lambda tt: s.sol(tt)[0] - l, t0, bg.tm, xtol=1e-12 * abs(t0), maxiter=500)))
    T += [("ventana", float(t)) for t in np.linspace(bg.tm, bg.tp, 120 * dens + 2)[1:-1]]
    s = bg.el; l0, l1 = s.sol(bg.tp)[0], s.sol(bg.t["S5"])[0]
    for l in np.linspace(l0, l1, 160 * dens + 1)[1:-1]:
        T.append(("expansion", brentq(lambda tt: s.sol(tt)[0] - l, bg.tp, bg.t["S5"], xtol=1e-12 * abs(bg.t["S5"]), maxiter=500)))
    return T


def measures(A, bg, n, at, t):
    out = dict(M1=0.0, M2=0.0, M3=0.0)
    band = n ** 1.5 * NRM
    for j in range(2):
        for c in range(2):
            b, qv, qd, u = Q.fields(A, bg, n, at, t, j, c); h = Q.step(b, n)
            bp = bg.bg(t + h); bm = bg.bg(t - h)
            phid = b["phid"]; phidd = (bp["phid"] - bm["phid"]) / (2 * h); xi = u[2] / phid
            out["M1"] = max(out["M1"], band * abs(qv[0] - b["H"] * xi))
            out["M2"] = max(out["M2"], band * math.sqrt(n * (n + 2.0)) / b["a"] * abs(xi))
            out["M3"] = max(out["M3"], band * abs(xi * phidd / phid))
    b = bg.bg(t); out.update(phid_over_H=b["phid"] / max(abs(b["H"]), 1e-300), lna=math.log(b["a"]), H=b["H"])
    return out


def row(args):
    par, n, key, rtol, dens = args
    t0_ = time.time(); A = U.Action(); bg = D.Background(par)
    try:
        at, info, chk = integrate(A, bg, n, key, rtol)
    except (RuntimeError, ValueError) as ex:                       # amendment 1: an in-state that is never defined before the window
        return dict(par=par, n=n, key=key, status="NO_DEFINIDA", error=str(ex))
    pts = [dict(phase=ph, **measures(A, bg, n, at, t)) for ph, t in times(bg, chk["t0"], dens)]
    mx = {k: max(p[k] for p in pts) for k in ("M1", "M2", "M3")}
    where = {k: (lambda p: dict(phase=p["phase"], lna=p["lna"], phid_over_H=p["phid_over_H"]))(max(pts, key=lambda p: p[k])) for k in mx}
    phase_max = {ph: {k: max([p[k] for p in pts if p["phase"] == ph] or [0.0]) for k in mx} for ph in ("contraccion", "ventana", "expansion")}
    return dict(par=par, n=n, key=key, rtol=rtol, dens=dens, status="OK", max=mx, where=where, phase_max=phase_max, checks=chk,
                shifts=info["shifts"], a_b=math.exp(bg.wl.sol(0.0)[0]), Hmax=max(abs(p["H"]) for p in pts if p["phase"] == "ventana"),
                phid_over_H_min={ph: min([p["phid_over_H"] for p in pts if p["phase"] == ph] or [math.inf]) for ph in ("contraccion", "ventana", "expansion")},
                seconds=time.time() - t0_)


def cached_row(args):
    par, n, key, rtol, dens = args
    f = RUNS / "rows" / ("%s_AH%g_f%.6f_n%d.json" % (key, par["AH"], par["f"], n)); f.parent.mkdir(parents=True, exist_ok=True)
    if f.exists(): return json.loads(f.read_text())
    r = row(args); f.write_text(json.dumps(r, default=float))
    with open(RUNS / "produce.out", "a") as log:
        log.write("row %s AH %g f %.4f n %d %s %.0fs\n" % (key, par["AH"], par["f"], n, r["status"], r.get("seconds", 0)))
    return r


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    pars, nm, hd = Q.pick_hosts()
    jobs = [(par, n, key, 1e-10, 1) for key in ("A2", "A1") for par in pars for n in NS]
    # each row is saved as soon as it ends and skipped on a relaunch (the container can restart without warning)
    with Pool(4) as p:
        R = p.map(cached_row, jobs, chunksize=1)
    (RUNS / "rows.json").write_text(json.dumps(dict(rows=R), default=float))
    print("rows done %.0fs" % (time.time() - t0), flush=True)
    main_rows = [r for r in R if r["key"] == "A2"]
    G0v = max(r["checks"]["G0"] for r in main_rows if r["status"] == "OK")
    G1v = max(max(abs(w - 1) for w in r["checks"]["W_end"]) for r in main_rows if r["status"] == "OK")
    g2, g3 = {}, {}
    mxf = lambda r: max(r["max"].values())
    for n in (100, 3000, 30000):
        base = next(r for r in main_rows if r["par"] == pars[0] and r["n"] == n)
        a9 = row((pars[0], n, "A2", 1e-9, 1)); a2x = row((pars[0], n, "A2", 1e-10, 2))
        g2[n] = abs(mxf(a9) / mxf(base) - 1); g3[n] = abs(mxf(a2x) / mxf(base) - 1)
        print("G2/G3", n, g2[n], g3[n], flush=True)
    gates = dict(G0=bool(G0v <= 1e-10), G1=bool(G1v <= 1e-6), G2=bool(max(g2.values()) <= 1e-3), G3=bool(max(g3.values()) <= 0.05))                                   # amendment 1: no extra gate
    hosts = []
    for par in pars:
        k_ = (par["AH"], par["f"]); h = hd[k_]
        s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))) if h.get("OmegaK0") else None; n30 = 30 / s if s else None
        rows = sorted([r for r in main_rows if r["par"] == par and r["status"] == "OK"], key=lambda r: r["n"])
        mx = [mxf(r) for r in rows]
        nsc = next((r["n"] for r, m in zip(rows, mx) if m >= 1), None)
        ir = max([m for r, m in zip(rows, mx) if n30 and r["n"] <= n30] or [0.0])
        a_b, Hmax = rows[0]["a_b"], rows[0]["Hmax"]
        turn = sorted([r for r in R if r["key"] == "A1" and r["par"] == par], key=lambda r: r["n"])
        hosts.append(dict(par=par, n_mem=nm[k_], n30=n30, n_sc=nsc, k_sc=(nsc / a_b) if nsc else None, k_mem=nm[k_] / a_b, Hmax=Hmax,
                          DGQ=Hmax ** (2 / 3), a_b=a_b, IR_max=ir,
                          per_n=[dict(n=r["n"], max=r["max"], where=r["where"], phase_max=r["phase_max"], phid_over_H_min=r["phid_over_H_min"]) for r in rows],
                          turnaround=[dict(n=r["n"], status=r["status"], max=r.get("max"), where=r.get("where"), shifts=r.get("shifts"),
                                           phase_max=r.get("phase_max")) for r in turn]))
    strong = any((H_["n_sc"] is not None and H_["n_sc"] <= H_["n_mem"]) or H_["IR_max"] >= 1 for H_ in hosts)
    pert = all((H_["n_sc"] is None or H_["n_sc"] >= 10 * H_["n_mem"]) and H_["IR_max"] <= 0.1 for H_ in hosts)
    would = "ACOPLAMIENTO_FUERTE" if strong else ("PERTURBATIVO" if pert else "LIMITE")
    verdict = would if all(gates.values()) else "INCONCLUSO"
    undefined_main = [(r["par"]["AH"], r["par"]["f"], r["n"]) for r in main_rows if r["status"] != "OK"]
    out = dict(spec="docs/preregistros/ESPEC_U63b.md", undefined_main=undefined_main, gates=gates, G0=G0v, G1=G1v, G2=g2, G3=g3, hosts=hosts, would_be=would, verdict=verdict,
               seconds=time.time() - t0)
    save("u63b_gauge_unitario.json", out)
    print(json.dumps(dict(gates=gates, G0=G0v, G1=G1v, G2=g2, G3=g3, would_be=would, verdict=verdict), indent=1, default=float))
    for H_ in hosts:
        print(H_["par"]["AH"], "%.3f" % H_["par"]["f"], "n_mem %.0f n_sc %s k_sc %s k_mem %.2e Hmax %.2e DGQ %.2e IR %.2e" %
              (H_["n_mem"], H_["n_sc"], ("%.2e" % H_["k_sc"]) if H_["k_sc"] else "-", H_["k_mem"], H_["Hmax"], H_["DGQ"], H_["IR_max"]))
        print("   main", [(p["n"], "%.1e/%.1e/%.1e" % (p["max"]["M1"], p["max"]["M2"], p["max"]["M3"]), p["where"]["M2"]["phase"]) for p in H_["per_n"]])
        print("   turn", [(p["n"], p["status"], ("%.1e/%.1e/%.1e" % (p["max"]["M1"], p["max"]["M2"], p["max"]["M3"])) if p["max"] else "-",
                          p["where"]["M2"]["phase"] if p["where"] else "-") for p in H_["turnaround"]], flush=True)


if __name__ == "__main__":
    main()
