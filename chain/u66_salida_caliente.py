"""U66 (docs/preregistros/ESPEC_U66.md): the hot exit of the cuscuton birth.

The 57 children of U57b (19 of the adopted set, 38 of C1) are evolved again from the mother's turnaround with the chain's
effective reheating switched on all along,

    chi'' + (3H + Gamma_chi) chi' + V'(chi) = 0,     rho_r' = -4 H rho_r + Gamma_chi chi'^2,

with rho_r a dynamical variable (in U57b it was C_rad/a^4).  Everything else is the system of chain/explora_cuscuton.py
(the class Cusc and the plan of legs of U57b).  After the end of inflation:
  phase 2: exact integration through the first 25 turning points of chi (chi' = 0), watching the branch W > 0;
  phase 3: chi as a w = 0 fluid decaying into radiation, until rho_chi/rho_r < 1e-10, requiring (2/3) rho_r a^2 > 1
           (the worst case of W, at a turning point of chi).

    python3 -B chain/u66_salida_caliente.py   ->  results/u66_salida_caliente.json   (a few minutes with 4 processes)"""
import json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from explora_cuscuton import Cusc, Vp, dVp, H_PLATEAU, N_HORIZON          # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
GAMMA = json.loads((ROOT / "data" / "cadena_cierre" / "thermal_pole.json").read_text())["continuous"]["Gamma_MPl"]
M_CHI = math.sqrt(json.loads((ROOT / "data" / "cadena_cierre" / "thermal_pole.json").read_text())["continuous"]["V0"]) / 0.3516
N_TURNS = 25
BREAK2 = 1e-12            # phase-2 branch threshold on W / W_scale (ESPEC_U66 §2)
OM_M = 0.315              # informative late-time reading only
LOG10 = math.log(10.0)


def deriv(cu, mode, gamma):
    def Q(y):
        lna, chi, p, d, rho = y
        rho = max(rho, 0.0)
        ia2 = math.exp(-2 * lna)
        phi, U, dU, ddU = cu.all(mode, d)
        H = dU / (3 * cu.mu2); D = 1.5 * cu.mu4 - ddU
        W = 2 / 3 * rho + 0.5 * p * p - ia2
        phidot = 3 * cu.mu2 * W / D
        return dict(H=H, D=D, W=W, Wscale=2 / 3 * rho + 0.5 * p * p + ia2, phidot=phidot, Hdot=ddU * phidot / (3 * cu.mu2),
                    rho=rho, ia2=ia2, U=U, K=0.5 * p * p, V=Vp(chi))

    def rhs(t, y):
        q = Q(y)
        return [q["H"], y[2], -(3 * q["H"] + gamma) * y[2] - dVp(y[1]), (1 if mode <= 0 else -1) * q["phidot"],
                -4 * q["H"] * q["rho"] + gamma * y[2] * y[2]]
    return Q, rhs


def resid(q):
    tot = q["rho"] + q["K"] + q["V"]
    return abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"])


ATOL = [1e-13, 1e-14, 1e-24, 1e-24, 1e-60]


def run(par, gamma, rtol, follow=True):
    """one U57b mother; returns phase-1 classification and, for children, phases 2-3"""
    cu = Cusc(par["mu"], par["m"], par["pinf"])
    a_t = par["AH"] / H_PLATEAU
    rho_t = 3 / (a_t ** 2 * (1 + par["f"]))                 # Crad / a_t^4 with Crad = 3 a_t^2 / (1 + f), as U57b
    y = [math.log(a_t), 0.0, math.sqrt(2 * par["f"] * rho_t), 0.0, rho_t]
    t, res1 = 0.0, 0.0
    out = dict(AH=par["AH"], f=par["f"])

    def leg(mode, y, t, target, direction, extra=None, keep=None):
        nonlocal res1
        Q, rhs = deriv(cu, mode, gamma)
        funcs = {"d": lambda t, y: y[3] - target[1], "acc": lambda t, y: (lambda q: q["Hdot"] + q["H"] ** 2)(Q(y)),
                 "turn": lambda t, y: y[3], "break": lambda t, y: Q(y)["W"] / Q(y)["Wscale"] - 1e-9}
        names = [target[0], "break"] + ([extra] if extra else [])
        evs = []
        for n in names:
            f = (lambda g: (lambda t, y: g(t, y)))(funcs[n]); f.terminal = True
            f.direction = direction if n == target[0] else -1; evs.append(f)
        q0 = Q(y)
        sol = solve_ivp(rhs, (t, t + 1e13), y, method="DOP853", rtol=rtol, atol=ATOL, events=evs,
                        first_step=min(1.0, 1e-3 / max(abs(q0["H"]), 1e-9)))
        chimax, rl = 0.0, 0.0
        for i in range(sol.y.shape[1]):
            yy = sol.y[:, i]; q = Q(yy); rl = max(rl, resid(q)); chimax = max(chimax, abs(yy[1]))
            if keep is not None: keep.append((yy[0], yy[0] + math.log(abs(q["H"]) + 1e-300)))
        hit = "numerico"
        if sol.status == 1: hit = [n for n, te in zip(names, sol.t_events) if len(te)][0]
        if hit != "break": res1 = max(res1, rl)
        return hit, list(sol.y[:, -1]), float(sol.t[-1]), chimax

    def project(y):
        ia2 = math.exp(-2 * y[0]); r = max(y[4], 0.0) + 0.5 * y[2] ** 2 + Vp(y[1]) - 3 * ia2
        g = lambda d: (lambda c: 3 * (c[2] / (3 * cu.mu2)) ** 2 - c[1])(cu.all(1, d)) - r
        lo, hi = y[3] * 0.9, y[3] * 1.1
        if g(lo) * g(hi) < 0: y[3] = brentq(g, lo, hi, xtol=1e-300, rtol=1e-15)
        return y

    half = cu.pinf / 2
    plan = [(-1, ("d", half), 1, None, "half-"), (0, ("d", 0.0), 1, None, "bounce"), (0, ("d", half), 1, None, "half+"),
            (1, ("acc", 0), -1, "turn", "decel"), (1, ("acc", 0), 1, "turn", "onset"), (1, ("acc", 0), -1, None, "end")]
    chi_land, lna_onset, samples, outcome = 0.0, None, [], None
    for mode, target, direction, extra, stage in plan:
        if stage == "bounce": y[3] = -half
        if stage == "decel": y[3] = half
        if stage in ("onset", "end"): y = project(y)
        hit, y, t, cm = leg(mode, y, t, target, direction, extra, keep=samples if stage == "end" else None)
        if stage in ("half+", "decel", "onset", "end"): chi_land = max(chi_land, cm)
        if hit == "numerico": outcome = "NUMERICO"; break
        if hit == "break": outcome = "RAMA_ROTA"; break
        if hit == "turn" and stage in ("decel", "onset"): outcome = "SIN_INFLACION"; break
        if stage == "onset": lna_onset = y[0]
        if stage == "end":
            N = y[0] - lna_onset; out["N_inf"] = N
            outcome = "HIJO" if N >= 60 else "INFLACION_CORTA"
            if N > N_HORIZON:                                   # log10 |Omega_K,0| = -2 log10(a H) at N_HORIZON before the end
                s = np.array(samples)
                out["log10_OmegaK0"] = float(-2 * np.interp(y[0] - N_HORIZON, s[:, 0], s[:, 1]) / LOG10)
    out.update(outcome=outcome, chi_land=chi_land, res_phase1=res1)
    if outcome != "HIJO" or not follow: return out

    # phase 2: exact oscillation through the first N_TURNS turning points of chi
    y = project(y)
    Q, rhs = deriv(cu, 1, gamma)
    ev_turn = lambda t, y: y[2]; ev_turn.terminal = False
    ev_brk = lambda t, y: Q(y)["W"] / Q(y)["Wscale"] - BREAK2; ev_brk.terminal = True; ev_brk.direction = -1
    T = 2 * math.pi / M_CHI
    sol = solve_ivp(rhs, (t, t + (N_TURNS + 3) * T), y, method="DOP853", rtol=rtol, atol=ATOL, events=[ev_turn, ev_brk],
                    max_step=T / 40)
    res2, minW = 0.0, 1.0
    for i in range(sol.y.shape[1]):
        q = Q(sol.y[:, i]); res2 = max(res2, resid(q)); minW = min(minW, q["W"] / q["Wscale"])
    turns = [yt for yt in sol.y_events[0]][:N_TURNS]
    margins = [math.log10(2 / 3 * max(yt[4], 1e-300)) + 2 * yt[0] / LOG10 for yt in turns]
    out.update(phase2_broken=bool(len(sol.t_events[1])), n_turns=len(turns), res_phase2=res2, minWrel_phase2=minW,
               M_turn1=margins[0] if margins else None, M_min_phase2=min(margins) if margins else None,
               M_turn_last=margins[-1] if margins else None)
    if out["phase2_broken"] or len(turns) < N_TURNS: return out

    # phase 3: averaged fluid from the last turning point (K = 0 there, so rho_chi = V)
    yl = turns[-1]; lna0 = yl[0]; epsG = cu.eps_G()

    def frhs(N, z):
        rc, rr = math.exp(z[0]), math.exp(z[1])
        H = math.sqrt(max((rc + rr) / (3 * (1 + epsG)) - math.exp(-2 * (lna0 + N)), 1e-300))
        g = gamma / H
        return [-3 - g, -4 + g * math.exp(z[0] - z[1])]
    ev_end = lambda N, z: (z[0] - z[1]) - math.log(1e-10); ev_end.terminal = True; ev_end.direction = -1
    fs = solve_ivp(frhs, (0, 60), [math.log(Vp(yl[1])), math.log(max(yl[4], 1e-300))], method="Radau", rtol=1e-10, atol=1e-12,
                   events=ev_end, max_step=0.02)
    M3 = [(math.log(2 / 3) + fs.y[1, i] + 2 * (lna0 + fs.t[i])) / LOG10 for i in range(fs.y.shape[1])]
    rr_end = math.exp(fs.y[1, -1])
    out.update(phase3_done=bool(len(fs.t_events[0])), M_min_phase3=min(M3), M_end_phase3=M3[-1], N_phase3=float(fs.t[-1]),
               T_end_GeV=(30 * rr_end / (math.pi ** 2 * 106.75)) ** 0.25 * 2.435e18)
    if out.get("log10_OmegaK0") is not None:                   # informative: margin today and e-folds of future (ESPEC §5)
        out["ln_margin_today"] = math.log(1.5 * OM_M) - out["log10_OmegaK0"] * LOG10
    return out


def is_hot(x):
    """the complete hot-exit classification of ESPEC_U66 §4 (inflation, phase 2 and phase 3)"""
    return (x["outcome"] == "HIJO" and not x.get("phase2_broken", True) and x.get("n_turns") == N_TURNS
            and bool(x.get("phase3_done")) and x.get("M_min_phase3", -1) > 0)


def job(args):
    row, gamma, rtol, follow = args
    p = dict(mu=row["par"]["mu"], m=row["par"]["m"], pinf=row["par"]["pinf"], AH=row["par"]["AH"], f=row["par"]["f"])
    r = run(p, gamma, rtol, follow); r["set"] = row["set"]; r["N_inf_U57b"] = row["N_inf"]; r["chi_land_U57b"] = row["chi_land"]
    return r


def main_u66b():
    """ESPEC_U66b: Gamma_chi x {0.1, 1, 10}; G0 on N_inf <= 1e-6 and chi_land <= 1e-5 (declared after U66)"""
    t0 = time.time()
    rows = [r for r in json.loads((ROOT / "results" / "u57b_runs" / "rows_rtol11.json").read_text())["rows"] if r["outcome"] == "HIJO"]
    facs = (0.1, 1.0, 10.0)
    with Pool(4) as pool:
        g0 = pool.map(job, [(r, 0.0, 1e-11, False) for r in rows])
        prod = {k: pool.map(job, [(r, k * GAMMA, 1e-11, True) for r in rows]) for k in facs}
        conv12 = {k: pool.map(job, [(r, k * GAMMA, 1e-12, True) for r in rows]) for k in facs}
    rel = lambda a, b: abs(a - b) / abs(b)
    G0_N = max(rel(x["N_inf"], x["N_inf_U57b"]) if x["outcome"] == "HIJO" else 1.0 for x in g0)
    G0_chi = max(rel(x["chi_land"], x["chi_land_U57b"]) if x["outcome"] == "HIJO" else 1.0 for x in g0)
    G1_res = max(max(x["res_phase1"], x.get("res_phase2", 0.0)) for k in facs for x in prod[k])
    G2_dM = max(abs(a["M_turn1"] - b["M_turn1"]) if a.get("M_turn1") is not None and b.get("M_turn1") is not None else 1.0
                for k in facs for a, b in zip(prod[k], conv12[k]))
    G2_class = all(a["outcome"] == b["outcome"] and is_hot(a) == is_hot(b) for k in facs for a, b in zip(prod[k], conv12[k]))
    gates = dict(G0=bool(G0_N <= 1e-6 and G0_chi <= 1e-5), G1=bool(G1_res <= 1e-6), G2=bool(G2_dM <= 0.01 and G2_class))
    ishot = is_hot
    allhot = all(ishot(x) for k in facs for x in prod[k])
    verdict = "INCONCLUSO" if not all(gates.values()) else ("SALIDA_CALIENTE" if allhot else "SALIDA_FRIA")
    def pick(k, s, key): return [x[key] for x in prod[k] if x["set"] == s and x.get(key) is not None]
    summ = {f"{k:g}": {s: dict(n=len([x for x in prod[k] if x["set"] == s]), n_hot=sum(ishot(x) for x in prod[k] if x["set"] == s),
                              incomplete=[dict(AH=x["AH"], f=x["f"]) for x in prod[k] if x["set"] == s and not ishot(x)],
                              M_turn1_min=min(pick(k, s, "M_turn1"), default=None),
                              M_min_phase3=min(pick(k, s, "M_min_phase3"), default=None),
                              dN_max=max((abs(x["N_inf"] - x["N_inf_U57b"]) for x in prod[k] if x["set"] == s and x.get("N_inf") is not None), default=None),
                              T_end_GeV_median=float(np.median(pick(k, s, "T_end_GeV"))) if pick(k, s, "T_end_GeV") else None)
                       for s in ("eps_pequeno", "C1")} for k in facs}
    out = dict(spec="docs/preregistros/ESPEC_U66b.md", Gamma_MPl=GAMMA, factors=facs, gates=gates, G0_max_rel_N=G0_N,
               G0_max_rel_chi_land=G0_chi, G1_max_res=G1_res, G2_max_dM=G2_dM, G2_same_class=G2_class, verdict=verdict,
               summary=summ, rows={f"{k:g}": prod[k] for k in facs}, seconds=time.time() - t0)
    conv = lambda o: o.item() if hasattr(o, "item") else str(o)
    (ROOT / "results" / "u66b_salida_caliente.json").write_text(json.dumps(out, indent=1, default=conv))
    print(json.dumps(dict(gates=gates, G0_N=G0_N, G0_chi=G0_chi, G1=G1_res, G2=G2_dM, verdict=verdict, summary=summ), indent=1, default=conv))


def main():
    t0 = time.time()
    rows = [r for r in json.loads((ROOT / "results" / "u57b_runs" / "rows_rtol11.json").read_text())["rows"] if r["outcome"] == "HIJO"]
    with Pool(4) as pool:
        g0 = pool.map(job, [(r, 0.0, 1e-11, False) for r in rows])
        main_ = pool.map(job, [(r, GAMMA, 1e-11, True) for r in rows])
        g2 = pool.map(job, [(r, GAMMA, 1e-12, True) for r in rows])
    rel = lambda a, b: abs(a - b) / abs(b)
    G0_err = max(max(rel(x["N_inf"], x["N_inf_U57b"]), rel(x["chi_land"], x["chi_land_U57b"])) if x["outcome"] == "HIJO" else 1.0
                 for x in g0)
    G1_res = max(max(x["res_phase1"], x.get("res_phase2", 0.0)) for x in main_)
    G2_dM = max(abs(a["M_turn1"] - b["M_turn1"]) if a.get("M_turn1") is not None and b.get("M_turn1") is not None else 1.0
                for a, b in zip(main_, g2))
    G2_class = all(a["outcome"] == b["outcome"] and is_hot(a) == is_hot(b) for a, b in zip(main_, g2))     # Codex review: full class
    gates = dict(G0=bool(G0_err <= 1e-6), G1=bool(G1_res <= 1e-6), G2=bool(G2_dM <= 0.01 and G2_class))
    hot = [is_hot(x) for x in main_]
    verdict = "INCONCLUSO" if not all(gates.values()) else ("SALIDA_CALIENTE" if all(hot) else "SALIDA_FRIA")
    summ = {}
    for s in ("eps_pequeno", "C1"):
        xs = [x for x in main_ if x["set"] == s]
        pick = lambda k: [x[k] for x in xs if x.get(k) is not None]            # rows that reached the phase (Codex review)
        summ[s] = dict(n=len(xs), n_hot=sum(h for h, x in zip(hot, main_) if x["set"] == s),
                       incomplete=[dict(AH=x["AH"], f=x["f"]) for x in xs if not is_hot(x)],
                       M_turn1_min=min(pick("M_turn1"), default=None), M_min_phase2=min(pick("M_min_phase2"), default=None),
                       M_min_phase3=min(pick("M_min_phase3"), default=None),
                       dN_max=max((abs(x["N_inf"] - x["N_inf_U57b"]) for x in xs if x.get("N_inf") is not None), default=None),
                       T_end_GeV_median=float(np.median(pick("T_end_GeV"))) if pick("T_end_GeV") else None,
                       future_efolds_min=min(pick("ln_margin_today"), default=None))
    out = dict(spec="docs/preregistros/ESPEC_U66.md", Gamma_MPl=GAMMA, gates=gates, G0_max_rel=G0_err, G1_max_res=G1_res,
               G2_max_dM=G2_dM, G2_same_class=G2_class, verdict=verdict, summary=summ, rows=main_,
               rows_rtol12=[dict(AH=x["AH"], f=x["f"], set=x["set"], M_turn1=x.get("M_turn1")) for x in g2],
               seconds=time.time() - t0)
    conv = lambda o: o.item() if hasattr(o, "item") else str(o)        # numpy scalars (first run crashed here, before writing)
    (ROOT / "results" / "u66_salida_caliente.json").write_text(json.dumps(out, indent=1, default=conv))
    print(json.dumps(dict(gates=gates, G0=G0_err, G1=G1_res, G2=G2_dM, verdict=verdict, summary=summ), indent=1, default=conv))


if __name__ == "__main__":
    main_u66b() if "--u66b" in sys.argv[1:] else main()
