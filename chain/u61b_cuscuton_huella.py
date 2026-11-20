"""U61b (docs/preregistros/ESPEC_U61b.md): U61 with new data (f at the quarter points of the U57/U57b grid) and the constraint
projection of explora_cuscuton.run (reimplemented: the mode +1 part is split into its decel / onset / end legs, plus a tail to
chi_end, each of the last three starting with the projection).  explora_cuscuton and u61_cuscuton_huella are imported unchanged.
    python3 -B chain/u61b_cuscuton_huella.py   ->  results/u61b_runs/rows.json, results/u61b_cuscuton_huella.json
Internal, not exported."""
import json, math, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, save
import explora_cuscuton as EC
import u61_cuscuton_huella as U61

RUNS = ROOT / "results" / "u61b_runs"
AHS = (3.0, 10.0, 20.0, 30.0, 48.0, 100.0, 300.0, 1000.0)
FS = tuple(sorted(10 ** ((k + o) / 16) for k in range(-64, 0) for o in (0.25, 0.75)))
ATOL = [1e-13, 1e-14, 1e-24, 1e-24]


def project(cu, Crad, y):
    ia2 = math.exp(-2 * y[0]); rhs_ = Crad * ia2 * ia2 + 0.5 * y[2] ** 2 + EC.Vp(y[1]) - 3 * ia2
    g = lambda d: (lambda r: 3 * (r[2] / (3 * cu.mu2)) ** 2 - r[1])(cu.all(1, d)) - rhs_
    lo, hi = y[3] * 0.9, y[3] * 1.1
    if g(lo) * g(hi) < 0:
        y[3] = brentq(g, lo, hi, xtol=1e-300, rtol=1e-15)
    return y


def host(par, rtol=1e-12):
    cu = EC.Cusc(par["mu"], par["m"], par["pinf"]); half = cu.pinf / 2
    a_t = par["AH"] / EC.H_PLATEAU; Crad = 3 * a_t ** 2 / (1 + par["f"]); rho_t = Crad / a_t ** 4
    y = [math.log(a_t), 0.0, math.sqrt(2 * par["f"] * rho_t), 0.0]
    legs = []; res = [0.0]; t = 0.0; landing = [None]; tb = [None]

    def resid(Q, yy):
        q = Q(yy); tot = q["rho"] + q["K"] + q["V"]
        return abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"])

    def leg(mode, y, t, events, name):
        Q, f = EC.deriv(cu, mode, Crad, 0.0)
        acc = lambda tt, yy: (lambda q: q["Hdot"] + q["H"] ** 2)(Q(yy))
        land = lambda tt, yy: yy[2] if (abs(yy[1]) > 1 and tb[0] is not None and tt > tb[0]) else 1.0
        evs = []
        for kind, direction, term in events:
            fn = {"acc": acc, "d_half": lambda tt, yy: yy[3] - half, "d0": lambda tt, yy: yy[3], "turn": lambda tt, yy: yy[3],
                  "chi_end": lambda tt, yy: abs(yy[1]) - U61.XEND}[kind]
            g = (lambda h: (lambda tt, yy: h(tt, yy)))(fn); g.direction = direction; g.terminal = term; evs.append(g)
        evs.append(land)
        s = solve_ivp(f, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=ATOL, events=evs, dense_output=True)
        res[0] = max(res[0], max(resid(Q, s.y[:, i]) for i in range(s.y.shape[1])))
        legs.append((name, s, Q))
        if landing[0] is None and len(s.t_events[-1]):
            landing[0] = (len(legs) - 1, float(s.t_events[-1][0]))
        hit = [kind for (kind, _, term), te in zip(events, s.t_events[:-1]) if term and len(te)]
        return s, Q, (hit[0] if hit else None)

    s, Q, _ = leg(-1, y, t, [("d_half", 1, True)], "half-"); t = s.t[-1]; y = list(s.y[:, -1]); y[3] = -half
    s, Q, _ = leg(0, y, t, [("d0", 1, False), ("d_half", 1, True)], "half+")
    tb[0] = float(s.t_events[0][0])
    # the landing event of this leg was evaluated before tb was known: redo the landing search on this leg
    ts = np.linspace(tb[0], s.t[-1], 4001); xd = s.sol(ts)[2]; xx = s.sol(ts)[1]
    for i in range(len(ts) - 1):
        if abs(xx[i]) > 1 and xd[i] * xd[i + 1] < 0:
            landing[0] = (len(legs) - 1, float(brentq(lambda tt: s.sol(tt)[2], ts[i], ts[i + 1]))); break
    else:
        landing[0] = None
    t = s.t[-1]; y = list(s.y[:, -1]); y[3] = half
    s, Q, hit = leg(1, y, t, [("acc", -1, True), ("turn", -1, True)], "decel")
    if hit == "turn":
        return dict(par=par, error="turn in decel")
    t = s.t[-1]; y = project(cu, Crad, list(s.y[:, -1]))
    s, Q, hit = leg(1, y, t, [("acc", 1, True), ("turn", -1, True)], "onset")
    if hit == "turn":
        return dict(par=par, error="turn in onset")
    t_on = s.t[-1]; lna_on = s.y[0, -1]; y = project(cu, Crad, list(s.y[:, -1]))
    s, Q, hit = leg(1, y, t_on, [("acc", -1, True)], "end")
    lna_endacc = s.y[0, -1]; t = s.t[-1]; y = project(cu, Crad, list(s.y[:, -1]))
    s, Q, hit = leg(1, y, t, [("chi_end", -1, True)], "tail")
    if landing[0] is None:
        return dict(par=par, error="no landing")
    li, tl = landing[0]; sL, QL = legs[li][1], legs[li][2]
    yl = sL.sol(tl); a_land = math.exp(yl[0]); chi_land = abs(yl[1])
    yend = s.y[:, -1]; N_tot = float(yend[0] - math.log(a_land))
    out = dict(par=par, chi_land=chi_land, N_tot=N_tot, N_acc=float(lna_endacc - lna_on), constraint_max=res[0])
    if N_tot >= U61.NSTAR:
        lnp = yend[0] - U61.NSTAR
        # the pivot lies in the 'end' or 'tail' leg (or earlier): search the legs after the landing
        for (name, sl, Ql) in legs[li:]:
            if sl.y[0, 0] <= lnp <= sl.y[0, -1]:
                tp = brentq(lambda tt: sl.sol(tt)[0] - lnp, sl.t[0], sl.t[-1], xtol=1e-14 * abs(sl.t[-1]), rtol=1e-15)
                yy = sl.sol(tp); aHs = math.exp(yy[0]) * Ql(yy)["H"]
                out.update(aH_star=aHs, OmegaK0=-(U61.KSTAR / aHs) ** 2); break
    vals = []
    for k_, (name, sl, Ql) in enumerate(legs):
        lo = max(sl.t[0], tb[0]); hi = min(sl.t[-1], tl) if k_ <= li else None
        if hi is None or hi <= lo:
            continue
        for tt in np.linspace(lo, hi, 3000):
            yy = sl.sol(tt); q = Ql(yy); a = math.exp(yy[0])
            vals.append((a * q["H"], a * math.sqrt(abs(q["Hdot"] + 2 * q["H"] ** 2 + q["ia2"]))))
    v = np.array(vals); out.update(aH_max=float(v[:, 0].max()), aR_max=float(v[:, 1].max()), n_mem=float(v.max()))
    return out


def job(par):
    try:
        return host(par)
    except Exception as ex:
        return dict(par=par, error=repr(ex))


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    path = RUNS / "rows.json"
    if path.exists():
        rows = json.loads(path.read_text())["rows"]
    else:
        jobs = [("eps_pequeno", dict(U61.ADOPTED, AH=A, f=f)) for A in AHS for f in FS]
        with Pool(4) as p:
            rows = p.map(EC.job, jobs, chunksize=1)
        path.write_text(json.dumps(dict(rows=rows), default=float))
    hosts = [r for r in rows if r.get("N_inf") and r["N_inf"] >= U61.NSTAR and r["outcome"] in ("HIJO", "INFLACION_CORTA")]
    print(f"rows {len(rows)}, hosts {len(hosts)} ({time.time() - t0:.0f}s)", flush=True)
    pars = [dict(U61.ADOPTED, AH=r["par"]["AH"], f=r["par"]["f"]) for r in hosts]
    with Pool(4) as p:
        R = p.map(job, pars, chunksize=1)
    by = {(r["par"]["AH"], r["par"]["f"]): r for r in hosts}
    g0 = []
    for o in R:
        ref = by[(o["par"]["AH"], o["par"]["f"])]
        g0.append(dict(AH=o["par"]["AH"], f=o["par"]["f"], err=o.get("error"),
                       chi_rel=abs(o["chi_land"] / ref["chi_land"] - 1) if "chi_land" in o else None,
                       dN=(o["N_acc"] - ref["N_inf"]) if "N_acc" in o else None))
    G0 = all(x["err"] is None and x["chi_rel"] <= 1e-5 and abs(x["dN"]) <= 0.01 for x in g0)
    G1 = all(o.get("constraint_max", 1) <= 1e-7 for o in R)
    low = sorted([o for o in R if "OmegaK0" in o], key=lambda o: o["N_tot"])[:5]
    g2 = []
    for o in low:
        o2 = host(o["par"], rtol=1e-11)
        g2.append(dict(par=o["par"], dOm=abs(o2["OmegaK0"] / o["OmegaK0"] - 1), dn=abs(o2["n_mem"] / o["n_mem"] - 1)))
    G2 = all(x["dOm"] <= 0.01 and x["dn"] <= 1e-4 for x in g2)
    per = defaultdict(lambda: dict(n=0, vis=0, curv=0, excl=0, vis_desi=0, excl_desi=0))
    for o in R:
        if "OmegaK0" not in o:
            continue
        A = o["par"]["AH"]; om = o["OmegaK0"]; x = U61.D_LSS * math.sqrt(abs(om))
        ell = math.sqrt(max(o["n_mem"] ** 2 - 1, 0)) * math.sin(x) if x < math.pi / 2 else float("nan")
        o["ell_mem"] = ell
        al, ad = om >= U61.OMK_PL, om >= U61.OMK_DESI
        d = per[A]; d["n"] += 1; d["vis"] += bool(al and ell >= 2); d["curv"] += bool(al and abs(om) >= 1e-3); d["excl"] += (not al)
        d["vis_desi"] += bool(ad and ell >= 2); d["excl_desi"] += (not ad)
    frac = {k: float(np.mean([d[k] / d["n"] for d in per.values() if d["n"]])) for k in ("vis", "curv", "excl", "vis_desi", "excl_desi")}
    nh = sum(1 for o in R if "OmegaK0" in o)
    gates = dict(G0=G0, G1=G1, G2=G2)
    verdict = "INCONCLUSO" if (not all(gates.values()) or nh < 20) else \
        ("HUELLA_OCULTA" if frac["vis"] == 0 else ("HUELLA_RARA" if frac["vis"] < 1e-2 else "HUELLA_POSIBLE"))
    out = dict(spec="docs/preregistros/ESPEC_U61b.md", n_rows=len(rows), n_hosts=nh, gates=gates, G0=g0, G2=g2, fractions=frac,
               per_AH=dict(per), verdict=verdict, hosts=R, seconds=time.time() - t0,
               n_mem_range=[min(o["n_mem"] for o in R if "n_mem" in o), max(o["n_mem"] for o in R if "n_mem" in o)],
               constraint_max=max(o.get("constraint_max", 0) for o in R),
               G0_worst=dict(chi=max((x["chi_rel"] or 0) for x in g0), dN=max(abs(x["dN"] or 0) for x in g0), errors=[x for x in g0 if x["err"]]))
    save("u61b_cuscuton_huella.json", out)
    print(json.dumps({k: out[k] for k in ("n_rows", "n_hosts", "gates", "G0_worst", "constraint_max", "fractions", "per_AH", "verdict", "n_mem_range")}, indent=1, default=float))
    for o in sorted([o for o in R if "OmegaK0" in o], key=lambda o: o["N_tot"])[:14]:
        print(o["par"]["AH"], round(o["par"]["f"], 4), "N %.2f chi %.3f OmK %.2e n_mem %.0f ell %.2f" % (o["N_tot"], o["chi_land"], o["OmegaK0"], o["n_mem"], o["ell_mem"]))


if __name__ == "__main__":
    main()
