"""Exploratory, no verdict (docs/informes/NOTA_CUSCUTON_V4.md): the homogeneous birth with a cuscuton instead of the NPQT
tower.  Independent re-implementation, in cosmic time, of the background of REVENGE C1 (whalefall-birth-explorer, branch
research/cuscuton-C1-C3-audit, revenge/C1/c1_cuscuton.py, which integrates in the cuscuton clock phi); no code is shared.

Action (M_Pl = 1, closed FLRW, k = +1), as in C1 and in Boruah, Kim, Rouben & Geshnizjani, JCAP 08 (2018) 031:

    R/2 - (d chi)^2/2 - V_pole(chi) - mu^2 sqrt(2X) - U(phi) + radiation,
    U = m^2 (phi^2 - phi_inf^2) - m^4 [exp((phi^2 - phi_inf^2)/m^2) - 1]

Branch phi_dot > 0:   U'(phi) = 3 mu^2 H      (H is slaved to phi)
                      phi_dot = 3 mu^2 W / D,  W = (2/3) rho_r + chi_dot^2/2 + sigma^2 - 1/a^2,  D = (3/2) mu^4 - U''
                      3 (H^2 + 1/a^2) = rho_r + chi_dot^2/2 + V + sigma^2 + U    (checked, never imposed)
sigma^2 = Sigma2/a^6 is an optional Bianchi-I-like shear energy (napkin only: no anisotropic curvature).
The branch exists only while W > 0 and D > 0; W -> 0 is reported as RAMA_ROTA.

A mother is a closed radiation universe at its turnaround (phi = -phi_inf, H = 0), of size A.H in plateau Hubble radii
(as the parents of explora_padres.py), with chi at its minimum and a kinetic fraction f = K_chi/rho_r.

    python3 -B chain/explora_cuscuton.py  ->  results/explora_cuscuton.json   (a few minutes with 4 processes)"""
import json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
V0, CP = 1.71311981e-10, 0.3516                      # pole plateau, as in C1 (the chain's V0 and C)
C1PAR = dict(mu=0.3, m=0.05, pinf=0.25)              # cuscuton parameters of C1 (those of the example of Boruah et al.)
H_PLATEAU = math.sqrt(V0 / 3)
N_HORIZON = 57.0                                     # e-folds before the end at which a0 H0 left the horizon (N_* = 51.6 at k = 0.05/Mpc, + ln 222)
RTOL = 1e-11


def Vp(x): return V0 * (1 - CP / math.hypot(CP, x))
def dVp(x): return V0 * CP * x / (CP * CP + x * x) ** 1.5


class Cusc:
    """U and its derivatives from (mode, d): phi = -pinf + d (mode -1), d (mode 0), pinf - d (mode +1); near +-pinf the
    exponent z = (phi^2 - pinf^2)/m^2 is computed from d, so H = U'/(3 mu^2) keeps its relative precision when H -> 0"""
    def __init__(self, mu, m, pinf):
        self.mu2, self.m2, self.pinf = mu * mu, m * m, pinf
        self.mu4 = self.mu2 ** 2

    def phi(self, mode, d): return d if mode == 0 else mode * (self.pinf - d)

    def all(self, mode, d):
        if mode != 0 and d < 0: d = 0.0             # a trial stage of the integrator beyond +-phi_inf (the events use the unclamped d)
        p = self.phi(mode, d)
        z = (p * p - self.pinf ** 2 if mode == 0 else d * d - 2 * self.pinf * d) / self.m2
        em = math.expm1(z)
        U = self.m2 * self.m2 * (z - em)
        dU = -2 * self.m2 * p * em
        ddU = -2 * self.m2 * em - 4 * p * p * (em + 1)
        return p, U, dU, ddU

    def eps_G(self): return 3 * self.mu4 / (8 * self.pinf ** 2)       # late times: 3 H^2 (1 + eps_G) = rho
    def U0(self): return self.all(0, 0.0)[1]


def deriv(cu, mode, Crad, Sig2):
    def quantities(y):
        lna, chi, p, d = y
        ia2 = math.exp(-2 * lna)
        rho = Crad * ia2 * ia2; sig = Sig2 * ia2 ** 3
        phi, U, dU, ddU = cu.all(mode, d)
        H = dU / (3 * cu.mu2); D = 1.5 * cu.mu4 - ddU
        W = 2 / 3 * rho + 0.5 * p * p + sig - ia2
        phidot = 3 * cu.mu2 * W / D
        Hdot = ddU * phidot / (3 * cu.mu2)
        return dict(H=H, D=D, W=W, phidot=phidot, Hdot=Hdot, rho=rho, sig=sig, ia2=ia2, U=U, phi=phi, K=0.5 * p * p, V=Vp(chi),
                    Wscale=2 / 3 * rho + 0.5 * p * p + sig + ia2)

    def rhs(t, y):
        q = quantities(y)
        return [q["H"], y[2], -3 * q["H"] * y[2] - dVp(y[1]), (1 if mode <= 0 else -1) * q["phidot"]]
    return quantities, rhs


def run(par):
    """one mother; par: AH, f (K_chi/rho_r at turnaround, chi at its minimum), s (shear/rho_r at turnaround), cuscuton
    parameters, cycles (> 1: U repeated periodically, so that the mother bounces again if it did not inflate), or the C1 data"""
    cu = Cusc(par.get("mu", C1PAR["mu"]), par.get("m", C1PAR["m"]), par.get("pinf", C1PAR["pinf"]))
    rtol = par.get("rtol", RTOL)
    out = dict(par=par, eps_G=cu.eps_G(), U0=cu.U0(), cycles=[])
    if "c1" in par:                                  # the initial data of C1, for the comparison row
        Crad, Sig2, y = par["c1"]["Crad"], 0.0, [math.log(par["c1"]["a"]), 0.0, par["c1"]["p"], par["c1"]["d"]]
    else:
        a_t = par["AH"] / H_PLATEAU                  # a_t = A
        Crad = 3 * a_t ** 2 / (1 + par["f"] + par.get("s", 0.0))
        rho_t = Crad / a_t ** 4
        Sig2 = par.get("s", 0.0) * rho_t * a_t ** 6
        y = [math.log(a_t), 0.0, math.sqrt(2 * par["f"] * rho_t), 0.0]
    t, res_max = 0.0, 0.0

    def leg(mode, y, t, target, direction, extra=None, keep=None, span=1e13):
        """integrates until the event `target` (name); returns (name of the event that stopped it, y, t, snapshot)"""
        nonlocal res_max
        Q, rhs = deriv(cu, mode, Crad, Sig2)
        acc = lambda t, y: (lambda q: q["Hdot"] + q["H"] ** 2)(Q(y))
        funcs = {"d": lambda t, y: y[3] - target[1], "acc": acc, "turn": lambda t, y: y[3],
                 "break": lambda t, y: Q(y)["W"] / Q(y)["Wscale"] - 1e-9}
        names = [target[0], "break"] + ([extra] if extra else [])
        evs = []
        for n in names:
            f = (lambda g: (lambda t, y: g(t, y)))(funcs[n])
            f.terminal = True; f.direction = direction if n == target[0] else -1
            evs.append(f)
        q0 = Q(y)
        sol = solve_ivp(rhs, (t, t + span), y, method="DOP853", rtol=rtol, atol=[1e-13, 1e-14, 1e-24, 1e-24], events=evs,
                        first_step=min(1.0, 1e-3 / max(abs(q0["H"]), 1e-9)))
        mW, chimax, res_leg = 1.0, 0.0, 0.0
        for i in range(sol.y.shape[1]):
            yy = sol.y[:, i]; q = Q(yy)
            tot = q["rho"] + q["K"] + q["V"] + q["sig"]
            res_leg = max(res_leg, abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"]))
            mW = min(mW, q["W"] / q["Wscale"]); chimax = max(chimax, abs(yy[1]))
            if keep is not None: keep.append((yy[0], math.exp(yy[0]) * q["H"]))
            if q["D"] <= 0: out["D_negative"] = True
        y = list(sol.y[:, -1]); q = Q(y)
        snap = dict(t=float(sol.t[-1]), a=math.exp(y[0]), chi=float(y[1]), chidot=float(y[2]), phi=q["phi"], H=q["H"], Hdot=q["Hdot"],
                    rho_r=q["rho"], K=q["K"], V=q["V"], U=q["U"], phidot=q["phidot"], D=q["D"], minWrel=mW, chimax=chimax)
        hit = "numerico"
        if sol.status == 1:
            hit = [n for n, te in zip(names, sol.t_events) if len(te)][0]
        if hit != "break": res_max = max(res_max, res_leg)       # the absolute error of the constraint is frozen while the terms dilute
        return hit, y, float(sol.t[-1]), snap

    def project(y):
        """late legs (mode +1): d is reset so that the Friedmann constraint holds; the drift it removes is in constraint_max"""
        ia2 = math.exp(-2 * y[0]); rhs_ = Crad * ia2 * ia2 + 0.5 * y[2] ** 2 + Vp(y[1]) + Sig2 * ia2 ** 3 - 3 * ia2
        g = lambda d: (lambda r: 3 * (r[2] / (3 * cu.mu2)) ** 2 - r[1])(cu.all(1, d)) - rhs_
        lo, hi = y[3] * 0.9, y[3] * 1.1
        if g(lo) * g(hi) < 0: y[3] = brentq(g, lo, hi, xtol=1e-300, rtol=1e-15)
        return y

    half = cu.pinf / 2
    for cyc in range(par.get("cycles", 1)):
        c = dict(cycle=cyc); out["cycles"].append(c); samples = []
        plan = [(-1, ("d", half), 1, None, "half-"), (0, ("d", 0.0), 1, None, "bounce"), (0, ("d", half), 1, None, "half+"),
                (1, ("acc", 0), -1, "turn", "decel"), (1, ("acc", 0), 1, "turn", "onset"), (1, ("acc", 0), -1, None, "end"),
                (1, ("turn", 0), -1, None, "post")]
        stop = None
        for mode, target, direction, extra, stage in plan:
            if stage == "post" and par.get("cycles", 1) == 1: break            # what a short inflation does next only matters for the cycles
            if stage == "bounce": y[3] = -half
            if stage == "decel": y[3] = half
            if stage == "half-" and cyc > 0: y[3] = max(-y[3], 1e-300)       # phi crosses phi_inf: the next well starts
            if stage in ("onset", "end", "post"): y = project(y)
            hit, y, t, snap = leg(mode, y, t, target, direction, extra, keep=samples if stage == "end" else None,
                                  span=3e8 if stage == "post" else 1e13)
            c[stage] = snap
            if stage in ("half+", "decel", "onset", "end"): c["chi_land"] = max(c.get("chi_land", 0.0), snap["chimax"])
            if stage in ("decel", "onset", "end"): c["minWrel_after_bounce"] = min(c.get("minWrel_after_bounce", 1.0), snap["minWrel"])
            if hit == "numerico": stop = "NUMERICO" if stage != "post" else "INFLACION_CORTA"
            elif hit == "break": stop = "RAMA_ROTA" if stage != "post" else "INFLACION_CORTA"      # post: first chi_dot = 0, no reheating here
            if stop: c["where"] = stage; break
            if hit == "turn" and stage in ("decel", "onset"): c["outcome"] = "SIN_INFLACION"; break       # turnaround: next cycle
            if stage == "end":
                N = math.log(c["end"]["a"] / c["onset"]["a"]); c["N_inf"] = N
                c["minus_U_over_3H2_end"] = -snap["U"] / (3 * snap["H"] ** 2)  # late-time shift of Friedmann, against eps_G
                if N >= 60:
                    c["outcome"] = "HIJO"
                    if N > N_HORIZON:
                        s = np.array(samples)
                        c["OmegaK0"] = float(np.exp(-2 * np.interp(math.log(c["end"]["a"]) - N_HORIZON, s[:, 0], np.log(s[:, 1]))))
                    break
                c["outcome"] = "INFLACION_CORTA"
            if stage == "post" and hit == "turn": c["outcome"] = "INFLACION_CORTA_RECOLAPSA"
        if "bounce" in c:
            b = c["bounce"]; b["rho_matter"] = b["rho_r"] + b["K"] + b["V"]; b["K_over_rho_r"] = b["K"] / b["rho_r"]
        if stop: c["outcome"] = stop
        if c["outcome"] not in ("SIN_INFLACION", "INFLACION_CORTA_RECOLAPSA"): break
    last = out["cycles"][-1]
    for k in ("outcome", "where", "N_inf", "chi_land", "OmegaK0", "minus_U_over_3H2_end", "minWrel_after_bounce", "bounce", "onset"):
        if k in last: out[k] = last[k]
    out["n_cycles"] = len(out["cycles"]); out["constraint_max"] = res_max
    if par.get("cycles", 1) == 1: del out["cycles"]
    else: out["cycles"] = [dict(cycle=c["cycle"], outcome=c["outcome"], N_inf=c.get("N_inf"), chi_land=c.get("chi_land"),
                                K_over_rho_r_bounce=c.get("bounce", {}).get("K_over_rho_r")) for c in out["cycles"]]
    return out


def c1_row():
    """C1's own initial data (phi_i = -phi_inf + 1e-6, C_rad from a_ref = 3e4, chi_dot_i = 4.3e-7), rebuilt here from its formulas"""
    cu = Cusc(**C1PAR); d = 1e-6
    Crad = (3 / 3e4 ** 2 - cu.U0()) * 3e4 ** 4
    _, U, dU, _ = cu.all(-1, d); Hi = dU / (3 * cu.mu2); p = 4.3e-7
    A = 3 * Hi ** 2 - U - 0.5 * p * p - Vp(0.0)
    a = math.sqrt(2 * Crad / (3 + math.sqrt(9 + 4 * A * Crad)))
    return dict(c1=dict(Crad=Crad, a=a, p=p, d=d))


def job(j):
    try: r = run(j[1])
    except (OverflowError, FloatingPointError, ValueError, ZeroDivisionError) as e:      # a row that leaves the domain is a result, not a crash
        r = dict(par=j[1], outcome="ERROR", error=repr(e), n_cycles=0, constraint_max=float("nan"))
    return dict(set=j[0], **r)


def main():
    t0 = time.time()
    AHS = (3.0, 10.0, 20.0, 30.0, 48.0, 100.0, 300.0, 1000.0)
    FS = tuple(10 ** (k / 16) for k in range(-64, 1))                         # f = 1e-4 ... 1, 16 per decade
    jobs = [("c1", c1_row()), ("c1_rtol", dict(c1_row(), rtol=1e-10))]
    jobs += [("grid", dict(AH=A, f=f)) for A in AHS for f in FS]
    jobs += [("shear48", dict(AH=48.0, f=0.076, s=s)) for s in (0.0, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1)]
    alt = dict(mu=0.16, m=0.02, pinf=0.25)                                     # a smaller late-time shift of G (eps_G = 0.4 %)
    jobs += [("alt", dict(alt, AH=A, f=f)) for A in (48.0, 100.0, 300.0) for f in FS]
    jobs += [("cycles", dict(AH=48.0, f=10 ** (k / 2), cycles=8)) for k in range(-32, -7)]   # f = 1e-16 ... 1e-4, periodic U
    rows_path = ROOT / "results" / "explora_cuscuton_rows.jsonl"              # one line per finished row: a restart loses nothing
    rows = []
    with Pool(4) as pool, open(rows_path, "w") as fh:
        for r in pool.imap_unordered(job, jobs, chunksize=1):
            rows.append(r); fh.write(json.dumps(r, default=float) + "\n"); fh.flush()
    order = {json.dumps([j[0], j[1]], sort_keys=True, default=float): i for i, j in enumerate(jobs)}
    rows.sort(key=lambda r: order[json.dumps([r["set"], r["par"]], sort_keys=True, default=float)])
    out = dict(kind="explora_cuscuton", status="EXPLORATORIO_SIN_VEREDICTO", V0=V0, C=CP, c1_parameters=C1PAR, H_plateau=H_PLATEAU,
               N_horizon=N_HORIZON, rtol=RTOL, rows=rows, seconds=time.time() - t0)
    (ROOT / "results" / "explora_cuscuton.json").write_text(json.dumps(out, indent=1, default=float))
    rows_path.unlink()
    summary(out)


def summary(out):
    from collections import Counter
    for r in out["rows"]:
        if r["set"] in ("c1", "c1_rtol"):
            b = r["bounce"]
            print("%s: a_b %.6f chi_b %.10f chidot_b %.9e Hdot_b %.8e rho_b %.6e N %.3f chi_land %.5f -U/3H^2 %.5f (eps_G %.5f) OmK0 %.1e minW %.3f constraint %.1e" % (
                r["set"], b["a"], b["chi"], b["chidot"], b["Hdot"], b["rho_matter"], r["N_inf"], r["chi_land"], r["minus_U_over_3H2_end"], r["eps_G"],
                r["OmegaK0"], r["minWrel_after_bounce"], r["constraint_max"]))
    for name in ("grid", "alt"):
        print(name, "(eps_G = %.4f, U0 = %.3e)" % next((r["eps_G"], r["U0"]) for r in out["rows"] if r["set"] == name and "eps_G" in r))
        for A in sorted({r["par"]["AH"] for r in out["rows"] if r["set"] == name}):
            rs = [r for r in out["rows"] if r["set"] == name and r["par"]["AH"] == A]
            kids = [r for r in rs if r["outcome"] == "HIJO"]
            line = "  AH %7.0f  %s" % (A, dict(Counter(r["outcome"] + ("@" + r["where"] if r["outcome"] == "RAMA_ROTA" else "") for r in rs)))
            if kids:
                line += "  HIJO: f %.3g-%.3g, chi_land %.2f-%.2f, N %.0f-%.0f, K/rho_r(b) %.2f-%.2f, rho_b %.3e, |OmK0| <= %.0e, minW >= %.2f" % (
                    min(r["par"]["f"] for r in kids), max(r["par"]["f"] for r in kids), min(r["chi_land"] for r in kids), max(r["chi_land"] for r in kids),
                    min(r["N_inf"] for r in kids), max(r["N_inf"] for r in kids), min(r["bounce"]["K_over_rho_r"] for r in kids),
                    max(r["bounce"]["K_over_rho_r"] for r in kids), max(r["bounce"]["rho_matter"] for r in kids),
                    max(r.get("OmegaK0", 0) for r in kids), min(r["minWrel_after_bounce"] for r in kids))
            br = [r for r in rs if r["outcome"] == "RAMA_ROTA"]
            if br: line += "  RAMA_ROTA: chi_land >= %.2f" % min(r.get("chi_land", float("inf")) for r in br)
            print(line)
        print("   max constraint %.1e, max rho_b %.3e, D<0: %s" % (np.nanmax([r["constraint_max"] for r in out["rows"] if r["set"] == name]),
              max(r["bounce"]["rho_matter"] for r in out["rows"] if r["set"] == name and "bounce" in r), any(r.get("D_negative") for r in out["rows"] if r["set"] == name)))
        print("   max constraint among HIJO %.1e" % max([r["constraint_max"] for r in out["rows"] if r["set"] == name and r["outcome"] == "HIJO"] or [float("nan")]))
    print("shear48")
    for r in out["rows"]:
        if r["set"] == "shear48": print("  s %.0e  %-16s N %6.2f chi_land %.3f" % (r["par"]["s"], r["outcome"], r.get("N_inf", float("nan")), r.get("chi_land", float("nan"))))
    print("cycles (periodic U)")
    for A in (48.0,):
        rs = [r for r in out["rows"] if r["set"] == "cycles" and r["par"]["AH"] == A and "cycles" in r]
        print("  AH %.0f  final: %s  cycles used: %s" % (A, dict(Counter(r["outcome"] for r in rs)), dict(Counter(r["n_cycles"] for r in rs))))
        amp = [c2["K_over_rho_r_bounce"] / c1_["K_over_rho_r_bounce"] for r in rs for c1_, c2 in zip(r["cycles"], r["cycles"][1:])
               if c1_["outcome"] == "SIN_INFLACION" and c2["K_over_rho_r_bounce"] and c2["K_over_rho_r_bounce"] < 1e-3]
        if amp: print("    gain of K/rho_r at the bounce per cycle (linear regime): median %.3g, range %.3g-%.3g (n = %d)" % (np.median(amp), min(amp), max(amp), len(amp)))
    print("seconds %.0f" % out["seconds"])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "summary": summary(json.loads((ROOT / "results" / "explora_cuscuton.json").read_text()))
    else: main()
