"""U62d (docs/preregistros/ESPEC_U62d.md): the primordial spectrum of the four visible cuscuton children of U61b on the chi > 0 sheet,
with the background integrated forward from the mother's turnaround (no constraint loss), projected to C_l with small-angle
Sachs-Wolfe on S^3 and compared with LCDM in l = 2-30 (Wishart Delta chi^2).  The propagation is that of U62c (imported unchanged):
two-field stage (U59 matrices) until rho_r/rho = 1e-4, then the Boruah-Kim-Geshnizjani action until k/aH = 1e-3, against the U58
attractor at |chi_x|.
    python3 -B chain/u62d_cl_rebote.py   ->  results/u62d_cl_rebote.json
Internal, not exported."""
import json, math, time
import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import spherical_jn
from ucommon import ROOT, save
import explora_cuscuton as EC
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62_diagnostico as DG
import u62b_bogoliubov_fisico as B62
import u62c_espectro_rebote as C

AT = B62.ATOL_BG
RUNS = ROOT / "results" / "u62d_runs"
ELLS = np.arange(2, 31)
W_VERDICT, W_ALT, FSKY = 0.6, 1.0, 0.8
NGRID = 100


class Shifted:
    """a solve_ivp solution re-labelled in time: t_new = t_old - t0"""
    def __init__(self, s, t0): self.s, self.t0 = s, t0; self.t = s.t - t0; self.y = s.y; self.status = s.status
    def sol(self, t): return self.s.sol(t + self.t0)


class Background:
    """forward from the turnaround: mode -1 to phi = -0.2, mode 0 to phi = +0.2, mode +1 to rho_r/rho = 0.9e-5"""

    def __init__(self, par):
        self.par = par; cu = self.cu = EC.Cusc(par["mu"], par["m"], par["pinf"])
        a_t = par["AH"] / EC.H_PLATEAU; Crad = self.Crad = 3 * a_t ** 2 / (1 + par["f"]); rho_t = Crad / a_t ** 4
        self.lna_t = math.log(a_t)
        y = [math.log(a_t), 0.0, math.sqrt(2 * par["f"] * rho_t), 0.0]
        self.Qm, rhs = EC.deriv(cu, -1, Crad, 0.0)
        ev = lambda t, y: y[3] - (cu.pinf - 0.2); ev.terminal = True
        self.cl = solve_ivp(rhs, (0.0, 1e15), y, method="DOP853", rtol=1e-12, atol=AT, dense_output=True, events=ev)
        self.tm = float(self.cl.t[-1]); y = list(self.cl.y[:, -1]); y[3] = -cu.pinf + y[3]
        self.Q0, rhs = EC.deriv(cu, 0, Crad, 0.0)
        ev = lambda t, y: y[3] - 0.2; ev.terminal = True
        self.wl = solve_ivp(rhs, (self.tm, self.tm + 1e15), y, method="DOP853", rtol=1e-12, atol=AT, dense_output=True, events=ev)
        self.tp = float(self.wl.t[-1]); y = list(self.wl.y[:, -1]); y[3] = cu.pinf - y[3]
        self.Qp, rhs = EC.deriv(cu, 1, Crad, 0.0)
        ev = lambda t, y: math.log(self.rr(y)) - math.log(0.9 * C.SWITCH_G4); ev.terminal = True
        self.el = solve_ivp(rhs, (self.tp, self.tp + 1e15), y, method="DOP853", rtol=1e-12, atol=AT, dense_output=True, events=ev)
        if not (self.cl.status == 1 and self.wl.status == 1 and self.el.status == 1): raise RuntimeError("background legs")
        # time origin at the bounce (phi = 0), as in U62b/U62c: the finite-difference steps of U62b.frame/end_info are relative to |t|
        tb = brentq(lambda t: self.wl.sol(t)[3], self.tm, self.tp, xtol=1e-12 * self.tp)
        self.cl, self.wl, self.el = Shifted(self.cl, tb), Shifted(self.wl, tb), Shifted(self.el, tb)
        self.tm, self.tp, self.t_turn = self.tm - tb, self.tp - tb, -tb
        self.lm, self.lp = float(self.cl.y[0, -1]), float(self.wl.y[0, -1])       # ln a at phi = -0.2, +0.2
        self.t = {}
        s = self.cl
        self.t["A1"] = brentq(lambda t: s.sol(t)[0] - (self.lna_t - 0.1), s.t[0], s.t[-1], xtol=1e-13 * abs(s.t[0]))
        self.t["A2"] = brentq(lambda t: s.sol(t)[0] - (self.lm + 2.0), s.t[0], s.t[-1], xtol=1e-13 * abs(s.t[0]))
        s = self.el
        self.t["B"] = brentq(lambda t: s.sol(t)[0] - (self.lp + 2.0), s.t[0], s.t[-1], xtol=1e-13 * s.t[-1])
        for name, target in (("S4", C.SWITCH), ("S5", C.SWITCH_G4)):
            self.t[name] = brentq(lambda t: math.log(self.rr(s.sol(t))) - math.log(target), self.t["B"], s.t[-1], xtol=1e-13 * s.t[-1])
        self.L = self.lna_t - self.lm

    def rr(self, y):
        q = self.Qp(list(y)); return q["rho"] / (q["rho"] + q["K"] + q["V"])

    def state(self, t):
        if t < self.tm: return -1, list(self.cl.sol(t)), self.Qm
        if t > self.tp: return 1, list(self.el.sol(t)), self.Qp
        return 0, list(self.wl.sol(t)), self.Q0

    def bg(self, t):
        m, y, _ = self.state(t)
        return DG.bg_mode(self.cu, self.Crad, y, m)

    def constraint_max(self):
        r = 0.0
        for s, Q in ((self.cl, self.Qm), (self.wl, self.Q0), (self.el, self.Qp)):
            for i in range(s.y.shape[1]):
                q = Q(list(s.y[:, i])); tot = q["rho"] + q["K"] + q["V"]
                r = max(r, abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"]))
        return r


class Attractor(C.Attractor):
    """U62c's attractor, returning (None, None) when |chi_x| is outside the attractor background (a mode that exits before landing)"""
    def P(self, chi_x, rtol=1e-10):
        x = abs(chi_x); lo, hi = float(self.s.sol(self.Nend - 1)[0]), float(self.s.sol(0.0)[0])
        if not (lo < x < hi): return None, None
        return C.Attractor.P(self, chi_x, rtol)


def grid_for(h):
    s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"])))
    n_lo = 0.25 * 2 / s; n_hi = min(max(2.5 * 30 / s, 3 * h["n_mem"]), 3e4)
    return sorted(set(int(round(x)) for x in np.logspace(math.log10(n_lo), math.log10(n_hi), NGRID))), s


def host_job(h):
    t0 = time.time(); A = U.Action(); att = Attractor(); par = h["par"]; bg = Background(par)
    grid, s = grid_for(h); nc = int(round(10 * h["n_mem"]))
    tag = "AH%g_f%.4f" % (par["AH"], par["f"]); log = open(RUNS / ("produce_%s.out" % tag), "a"); rows = []
    for n in grid:
        r = C.spectrum_row(A, bg, att, n); r["ell"] = n * s; rows.append(r)
        S = r["S4"]
        log.write(json.dumps(dict(n=n, ell=n * s, F_A1=S.get("F_tot_A1"), F_A2=S.get("F_tot_A2"), F_B=S.get("F_tot_B"), chi_x=S["chi_exit"],
                                  W=[S[X]["W_end"] for X in ("A1", "A2", "B")], A1_shifts=S["A1"]["state"]["shifts"], sec=time.time() - t0)) + "\n"); log.flush()
    rc = C.spectrum_row(A, bg, att, nc, states=("A2", "B")); rc["ell"] = nc * s
    log.write("control " + json.dumps(dict(n=nc, F_A2=rc["S4"].get("F_tot_A2"), sec=time.time() - t0)) + "\n"); log.close()
    return dict(par=par, OmegaK0=h["OmegaK0"], n_mem=h["n_mem"], ell_mem=h["ell_mem"], L_contraction=bg.L, t=bg.t, sin_chiLS=s,
                chi_LS=M.D_LSS * math.sqrt(abs(h["OmegaK0"])), constraint_max=bg.constraint_max(), grid=grid, rows=rows, control=rc,
                seconds=time.time() - t0)


def project(H, X, every=1):
    """R_l = C_l[F P_att] / C_l[P_att] (small-angle Sachs-Wolfe on S^3); F from the rows (every-th point), 1 outside the grid"""
    rows = [r for r in H["rows"] if r["S4"].get("F_tot_" + X) and r["S4"].get("P_att")][::every]
    ln_n = np.log([r["n"] for r in rows]); lnF = np.log([r["S4"]["F_tot_" + X] for r in rows])
    allr = [r for r in H["rows"] if r["S4"].get("P_att")]
    ln_na = np.log([r["n"] for r in allr]); lnPa = np.log([r["S4"]["P_att"] for r in allr])
    slope, icpt = np.polyfit(ln_na, lnPa, 1)
    chi = H["chi_LS"]; nmax = int(50 * 30 / H["sin_chiLS"])
    n = np.arange(2, nmax + 1, dtype=float); ln = np.log(n)
    lnP = np.where((ln >= ln_na[0]) & (ln <= ln_na[-1]), np.interp(ln, ln_na, lnPa), icpt + slope * ln)
    F = np.where((ln >= ln_n[0]) & (ln <= ln_n[-1]), np.exp(np.interp(ln, ln_n, lnF)), 1.0)
    Pa = np.exp(lnP); R = []; low = None
    for l in ELLS:
        K = spherical_jn(int(l), n * chi) ** 2 / n
        R.append(float(np.sum(F * Pa * K) / np.sum(Pa * K)))
        if l == 2: low = float(np.sum((Pa * K)[n < math.exp(ln_n[0])]) / np.sum(Pa * K))
    return np.array(R), low


def dchi2(R, w):
    Rt = 1 + w * (R - 1)
    return float(FSKY * np.sum((2 * ELLS + 1) * (1 / Rt + np.log(Rt) - 1)))


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    hosts = [h for h in M.visible_hosts() if not (h["par"]["AH"] == 1000.0 and abs(h["par"]["f"] - 0.8977) < 1e-3)]
    with Pool(min(4, len(hosts))) as pool:
        res = pool.map(host_job, hosts, chunksize=1)
    out = dict(spec="docs/preregistros/ESPEC_U62d.md", hosts=res, excluded_host="AH=1000, f=0.898 (hoja chi<0)")
    save("u62d_cl_rebote.json", out)
    # G2, G4 on the first host
    A = U.Action(); att = Attractor(); h0 = hosts[0]; bg = Background(h0["par"]); r0 = res[0]; s = r0["sin_chiLS"]
    pick = [min(r0["grid"], key=lambda n: abs(n - l / s)) for l in (2, 10, 30)]
    g2, g4 = {}, {}
    for n in pick:
        base = next(r for r in r0["rows"] if r["n"] == n)["S4"]["A1"]["P"]
        g2[n] = abs(C.spectrum_row(A, bg, att, n, rtol=1e-9, states=("A1",))["S4"]["A1"]["P"] / base - 1)
        g4[n] = abs(C.spectrum_row(A, bg, att, n, switches=("S5",), states=("A1",))["S5"]["A1"]["P"] / base - 1)
        print("G2/G4", n, g2[n], g4[n], flush=True)
    G0v = []
    for n in r0["grid"][:3]:
        for X, tow in (("A1", 1), ("A2", 1), ("B", -1)):
            e = B62.end_info(A, bg, n, bg.t[X], tow)
            if e["pos"]:
                F = e["F"]; G0v.append(max(abs(M.Om(F[:, j], F[:, k].conj()) - (1j if j == k else 0)) for j in range(2) for k in range(2)))
    allrows = [r for H in res for r in H["rows"]] + [H["control"] for H in res]
    W = [abs(r["S4"][X]["W_end"] / 0.25 - 1) for r in allrows for X in ("A1", "A2", "B") if X in r["S4"]]
    ctrl = [abs(H["control"]["S4"]["F_tot_A2"] - 1) if H["control"]["S4"].get("F_tot_A2") else math.inf for H in res]
    proj = []; g6 = []
    for H in res:
        p = dict(par=H["par"])
        for X in ("A1", "A2"):
            R, low = project(H, X); R2, _ = project(H, X, every=2)
            d, d2 = dchi2(R, W_VERDICT), dchi2(R2, W_VERDICT)
            p[X] = dict(R_l=R.tolist(), dchi2_w06=d, dchi2_w1=dchi2(R, W_ALT), dchi2_w06_half=d2, low_n_fraction_C2=low)
            g6.append(abs(d2 - d) <= max(2.0, 0.2 * d))
        p["dep_A1_A2_median"] = float(np.median([abs(r["S4"]["A1"]["P"] / r["S4"]["A2"]["P"] - 1) for r in H["rows"]]))
        p["excluded"] = bool(p["A1"]["dchi2_w06"] >= 25 and p["A2"]["dchi2_w06"] >= 25)
        p["compatible"] = bool(p["A1"]["dchi2_w06"] < 9 and p["A2"]["dchi2_w06"] < 9)
        proj.append(p)
    gates = dict(G0=bool(max(G0v) <= 1e-10), G1=bool(max(W) <= 1e-3), G2=bool(max(g2.values()) <= 1e-3), G3=bool(max(ctrl) <= 0.03),
                 G4=bool(max(g4.values()) <= 1e-2), G5=bool(max(H["constraint_max"] for H in res) <= 3e-5), G6=bool(all(g6)))
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif all(p["excluded"] for p in proj): verdict = "VISIBLES_EXCLUIDOS"
    elif all(p["compatible"] for p in proj): verdict = "VISIBLES_COMPATIBLES"
    else: verdict = "MIXTO"
    would = ("VISIBLES_EXCLUIDOS" if all(p["excluded"] for p in proj) else ("VISIBLES_COMPATIBLES" if all(p["compatible"] for p in proj) else "MIXTO"))
    out.update(gates=gates, G0_max=max(G0v), G1_max=max(W), G2=g2, G3_control=ctrl, G4=g4, constraint_max=max(H["constraint_max"] for H in res),
               projection=proj, would_be=would, verdict=verdict, seconds=time.time() - t0)
    save("u62d_cl_rebote.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_max", "G1_max", "G2", "G3_control", "G4", "constraint_max", "would_be", "verdict", "seconds")},
                     indent=1, default=float))
    for p in proj:
        print(p["par"]["AH"], "%.3f" % p["par"]["f"], {X: (round(p[X]["dchi2_w06"], 2), round(p[X]["dchi2_w1"], 2), round(p[X]["dchi2_w06_half"], 2),
                                                          [round(x, 2) for x in p[X]["R_l"][:6]]) for X in ("A1", "A2")})


if __name__ == "__main__":
    main()
