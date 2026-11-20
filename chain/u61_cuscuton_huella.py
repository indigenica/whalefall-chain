"""U61 (docs/preregistros/ESPEC_U61.md): curvature and bounce imprint with the cuscuton children (adopted set).
The background is integrated here with explora_cuscuton.deriv (unchanged), in legs: mode -1 to phi = -phi_inf/2, mode 0 through
the bounce to +phi_inf/2, mode +1 to the end of inflation at the chain's chi_end.
    python3 -B chain/u61_cuscuton_huella.py   ->  results/u61_cuscuton_huella.json
Internal, not exported."""
import json, math, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, benchmark, save
import explora_cuscuton as EC

ADOPTED = dict(mu=0.16, m=0.02, pinf=0.25)
B = benchmark(); NSTAR, XEND = B["N_star"], B["chi_end"]
KSTAR = 0.05 * 299792.458 / (100 * 0.6736)
D_LSS = 3.1158200301196635                                           # U55 (G3 of U55: flat LCDM, Planck 2018)
OMK_PL, OMK_DESI = 0.0007 - 2 * 0.0019, 0.0023 - 2 * 0.0011


def host(par, rtol=1e-12):
    cu = EC.Cusc(par["mu"], par["m"], par["pinf"]); half = cu.pinf / 2
    a_t = par["AH"] / EC.H_PLATEAU; Crad = 3 * a_t ** 2 / (1 + par["f"]); rho_t = Crad / a_t ** 4
    y = [math.log(a_t), 0.0, math.sqrt(2 * par["f"] * rho_t), 0.0]
    atol = [1e-13, 1e-14, 1e-24, 1e-24]; res = [0.0]; acc_ev = []; landing = [None]
    def resid(Q, yy):
        q = Q(yy); tot = q["rho"] + q["K"] + q["V"]
        return abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"])
    # leg 1: mode -1 to phi = -phi_inf/2
    Q1, f1 = EC.deriv(cu, -1, Crad, 0.0)
    e = lambda t, yy: yy[3] - half; e.terminal = True; e.direction = 1
    s1 = solve_ivp(f1, (0, 1e13), y, method="DOP853", rtol=rtol, atol=atol, events=e)
    res[0] = max(res[0], max(resid(Q1, s1.y[:, i]) for i in range(s1.y.shape[1])))
    t = s1.t[-1]; y = list(s1.y[:, -1]); y[3] = -half
    # leg 2: mode 0 through the bounce
    Q2, f2 = EC.deriv(cu, 0, Crad, 0.0)
    eb = lambda t, yy: yy[3]; eb.direction = 1
    ee = lambda t, yy: yy[3] - half; ee.terminal = True; ee.direction = 1
    el = lambda t, yy: yy[2] if abs(yy[1]) > 1 else 1.0
    ea = lambda t, yy: (lambda q: q["Hdot"] + q["H"] ** 2)(Q2(yy))
    s2 = solve_ivp(f2, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=atol, events=[eb, ee, el, ea], dense_output=True)
    res[0] = max(res[0], max(resid(Q2, s2.y[:, i]) for i in range(s2.y.shape[1])))
    tb = float(s2.t_events[0][0])
    acc_ev += [(float(tt), "2") for tt in s2.t_events[3] if tt > tb]
    if len(s2.t_events[2]): landing[0] = ("2", float(s2.t_events[2][0]))
    t = s2.t[-1]; y = list(s2.y[:, -1]); y[3] = half
    # leg 3: mode +1 to the landing (if not yet) and to |chi| = chi_end
    Q3, f3 = EC.deriv(cu, 1, Crad, 0.0)
    el3 = lambda t, yy: yy[2] if abs(yy[1]) > 1 else 1.0
    landed = [landing[0] is not None]
    def eend(t, yy): return abs(yy[1]) - XEND if landed[0] else 1.0
    eend.terminal = True; eend.direction = -1
    ea3 = lambda t, yy: (lambda q: q["Hdot"] + q["H"] ** 2)(Q3(yy))
    if not landed[0]:
        el3.terminal = True
        s3a = solve_ivp(f3, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=atol, events=[el3, ea3], dense_output=True)
        res[0] = max(res[0], max(resid(Q3, s3a.y[:, i]) for i in range(s3a.y.shape[1])))
        landing[0] = ("3a", float(s3a.t_events[0][0])); acc_ev += [(float(tt), "3a") for tt in s3a.t_events[1]]
        t = s3a.t[-1]; y = list(s3a.y[:, -1]); landed[0] = True
    else:
        s3a = None
    s3 = solve_ivp(f3, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=atol, events=[eend, ea3], dense_output=True)
    res[0] = max(res[0], max(resid(Q3, s3.y[:, i]) for i in range(s3.y.shape[1])))
    acc_ev += [(float(tt), "3") for tt in s3.t_events[1]]
    sols = {"2": (s2, Q2), "3a": (s3a, Q3), "3": (s3, Q3)}
    def state(leg, tt): s, Q = sols[leg]; yy = s.sol(tt); return yy, Q(yy)
    # landing, end, pivot
    yl, _ = state(*landing[0]); a_land = math.exp(yl[0]); chi_land = abs(yl[1])
    yend = s3.y[:, -1]; a_end = math.exp(yend[0]); N_tot = math.log(a_end / a_land)
    out = dict(par=par, chi_land=chi_land, N_tot=N_tot, constraint_max=res[0])
    if N_tot >= NSTAR:
        lnp = yend[0] - NSTAR
        tp = brentq(lambda tt: s3.sol(tt)[0] - lnp, s3.t[0], s3.t[-1], xtol=1e-14 * abs(s3.t[-1]), rtol=1e-15)
        yy, q = state("3", tp); aHs = math.exp(yy[0]) * q["H"]
        out.update(aH_star=aHs, OmegaK0=-(KSTAR / aHs) ** 2)
    # n_mem between the bounce and the landing
    vals = []
    legs = [("2", tb, landing[0][1] if landing[0][0] == "2" else s2.t[-1])]
    if landing[0][0] == "3a": legs.append(("3a", s3a.t[0], landing[0][1]))
    for leg, ta, tz in legs:
        for tt in np.linspace(ta, tz, 4000):
            yy, q = state(leg, tt); a = math.exp(yy[0])
            vals.append((a * q["H"], a * math.sqrt(abs(q["Hdot"] + 2 * q["H"] ** 2 + q["ia2"]))))
    v = np.array(vals); out.update(aH_max=float(v[:, 0].max()), aR_max=float(v[:, 1].max()), n_mem=float(v.max()))
    # e-folds of the accelerated phase, for G0 (onset: last up-crossing of a-ddot before the end; end: first down-crossing after it)
    ups, downs = [], []
    for tt, leg in sorted(acc_ev):
        yy, q = state(leg, tt); h = 1e-6 * max(1.0, abs(tt))
        y2, q2 = state(leg, tt + h) if (leg != "3" or tt + h <= s3.t[-1]) else state(leg, tt - h)
        (ups if (q2["Hdot"] + q2["H"] ** 2) > 0 else downs).append((tt, leg, yy[0]))
    end = next((d for d in downs if any(u[0] < d[0] for u in ups) and d[0] > landing[0][1]), None)
    if end:
        on = max((u for u in ups if u[0] < end[0]), key=lambda u: u[0])
        out["N_acc"] = end[2] - on[2]
    return out


def job(par):
    try:
        return host(par)
    except Exception as ex:
        return dict(par=par, error=repr(ex))


def main():
    t0 = time.time()
    rows = json.loads((ROOT / "results" / "u57_runs" / "rows_rtol11.json").read_text())["rows"] + \
           json.loads((ROOT / "results" / "u57b_runs" / "rows_rtol11.json").read_text())["rows"]
    rows = [r for r in rows if r["set"] == "eps_pequeno"]
    hosts = [r for r in rows if r.get("N_inf") and r["N_inf"] >= NSTAR and r["outcome"] in ("HIJO", "INFLACION_CORTA")]
    pars = [dict(ADOPTED, AH=r["par"]["AH"], f=r["par"]["f"]) for r in hosts]
    with Pool(4) as p:
        R = p.map(job, pars, chunksize=1)
    by = {(r["par"]["AH"], r["par"]["f"]): r for r in hosts}
    g0 = []
    for o in R:
        ref = by[(o["par"]["AH"], o["par"]["f"])]
        g0.append(dict(AH=o["par"]["AH"], f=o["par"]["f"], err=o.get("error"),
                       chi_rel=abs(o["chi_land"] / ref["chi_land"] - 1) if "chi_land" in o else None,
                       dN=(o.get("N_acc", float("nan")) - ref["N_inf"]) if "N_acc" in o else None))
    G0 = all(x["err"] is None and x["chi_rel"] is not None and x["chi_rel"] <= 1e-6 and x["dN"] is not None and abs(x["dN"]) <= 0.05 for x in g0)
    G1 = all(o.get("constraint_max", 1) <= 1e-8 for o in R)
    low = sorted([o for o in R if "OmegaK0" in o], key=lambda o: o["N_tot"])[:5]
    g2 = []
    for o in low:
        o2 = host(o["par"], rtol=1e-11)
        g2.append(dict(par=o["par"], dOm=abs(o2["OmegaK0"] / o["OmegaK0"] - 1), dn=abs(o2["n_mem"] / o["n_mem"] - 1)))
    G2 = all(x["dOm"] <= 0.01 and x["dn"] <= 1e-4 for x in g2)
    # fractions (per A.H, then averaged)
    per = defaultdict(lambda: dict(n=0, vis=0, curv=0, excl=0, vis_desi=0, excl_desi=0))
    for o in R:
        if "OmegaK0" not in o:
            continue
        A = o["par"]["AH"]; om = o["OmegaK0"]; ell = math.sqrt(max(o["n_mem"] ** 2 - 1, 0)) * math.sin(D_LSS * math.sqrt(abs(om)))
        o["ell_mem"] = ell
        al, ad = om >= OMK_PL, om >= OMK_DESI
        d = per[A]; d["n"] += 1; d["vis"] += bool(al and ell >= 2); d["curv"] += bool(al and abs(om) >= 1e-3); d["excl"] += (not al)
        d["vis_desi"] += bool(ad and ell >= 2); d["excl_desi"] += (not ad)
    frac = {k: float(np.mean([d[k] / d["n"] for d in per.values() if d["n"]])) for k in ("vis", "curv", "excl", "vis_desi", "excl_desi")}
    nh = sum(1 for o in R if "OmegaK0" in o)
    gates = dict(G0=G0, G1=G1, G2=G2)
    if not all(gates.values()) or nh < 20:
        verdict = "INCONCLUSO"
    else:
        verdict = "HUELLA_OCULTA" if frac["vis"] == 0 else ("HUELLA_RARA" if frac["vis"] < 1e-2 else "HUELLA_POSIBLE")
    out = dict(spec="docs/preregistros/ESPEC_U61.md", n_hosts=nh, gates=gates, G0=g0, G2=g2, fractions=frac, per_AH=dict(per),
               verdict=verdict, hosts=R, seconds=time.time() - t0,
               n_mem_range=[min(o["n_mem"] for o in R if "n_mem" in o), max(o["n_mem"] for o in R if "n_mem" in o)],
               OmK_range=[min(o["OmegaK0"] for o in R if "OmegaK0" in o), max(o["OmegaK0"] for o in R if "OmegaK0" in o)],
               constraint_max=max(o.get("constraint_max", 0) for o in R))
    save("u61_cuscuton_huella.json", out)
    print(json.dumps({k: out[k] for k in ("n_hosts", "gates", "fractions", "verdict", "n_mem_range", "OmK_range", "constraint_max", "G2")}, indent=1, default=float))
    print("G0 worst", max((x["chi_rel"] or 0) for x in g0), max(abs(x["dN"] or 0) for x in g0), [x for x in g0 if x["err"]])
    for o in sorted([o for o in R if "OmegaK0" in o], key=lambda o: o["N_tot"])[:12]:
        print(o["par"]["AH"], round(o["par"]["f"], 4), "N %.2f chi %.3f OmK %.2e n_mem %.1f ell %.2f" % (o["N_tot"], o["chi_land"], o["OmegaK0"], o["n_mem"], o["ell_mem"]))


if __name__ == "__main__":
    main()
