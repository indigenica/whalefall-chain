"""U62c (docs/preregistros/ESPEC_U62c.md): the primordial spectrum of the five visible cuscuton children of U61b, propagated from
the contracting mother through the bounce to horizon exit.  Stage 1: the U59 two-field equations (uniform-radiation gauge) on a
contraction leg (mode -1, backwards), the U59 window and an expansion leg (mode +1) until rho_r/rho = 1e-4 (1e-5 for G4); in-states
are the instantaneous vacua of the B_s-free frame of U62b at A1 (0.5 e-folds from the mother's turnaround), A2 (Delta ln a = 2 before
the window) and B (Delta ln a = 2 after).  Stage 2: zeta_C = zeta_R - (H/chidot) x_R, then the Boruah-Kim-Geshnizjani action
(u58_cuscuton_inflacion.z2_cs2) on a flat U58 background started from the child's (chi, chi_N), until k/aH = 1e-3.  Reference: the U58
Bunch-Davies spectrum on the attractor at the chi where the mode exits.  U58, U59, U62 and U62b are imported unchanged.
    python3 -B chain/u62c_espectro_rebote.py   ->  results/u62c_espectro_rebote.json
Internal, not exported."""
import json, math, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from ucommon import ROOT, save
import explora_cuscuton as EC
import u58_cuscuton_inflacion as U58
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62_diagnostico as DG
import u62b_bogoliubov_fisico as B62

ATOL_BG = B62.ATOL_BG
RUNS = ROOT / "results" / "u62c_runs"
SWITCH = 1e-4
SWITCH_G4 = 1e-5


class Background:
    def __init__(self, par):
        self.par = par
        cu, Crad, yb, _ = U.host_state(par); self.cu, self.Crad = cu, Crad
        self.tr = tr = U.Trajectory(cu, Crad, yb); self.tm, self.tp = tr.tm, tr.tp
        self.L = math.log(par["AH"] / EC.H_PLATEAU) - tr.y(tr.tm)[0]          # contraction length, phi = -0.2 -> turnaround
        # contraction leg, backwards from phi = -0.2 to 0.45 e-folds before the turnaround
        y0 = list(tr.y(tr.tm)); self.l0m = y0[0]; y0[3] = cu.pinf + y0[3]
        self.Qm, rhs = EC.deriv(cu, -1, Crad, 0.0)
        ev = lambda t, y: y[0] - self.l0m - (self.L - 0.45); ev.terminal = True
        self.cl = solve_ivp(rhs, (tr.tm, tr.tm - 1e15), y0, method="DOP853", rtol=1e-12, atol=ATOL_BG, dense_output=True, events=ev)
        if self.cl.status != 1: raise RuntimeError("contraction leg")
        # expansion leg, from phi = +0.2 until rho_r/rho = 0.9e-5
        y0 = list(tr.y(tr.tp)); self.l0p = y0[0]; y0[3] = cu.pinf - y0[3]
        self.Qp, rhs = EC.deriv(cu, 1, Crad, 0.0)
        ev = lambda t, y: math.log(self.rr(y)) - math.log(0.9 * SWITCH_G4); ev.terminal = True
        self.el = solve_ivp(rhs, (tr.tp, tr.tp + 1e15), y0, method="DOP853", rtol=1e-12, atol=ATOL_BG, dense_output=True, events=ev)
        if self.el.status != 1: raise RuntimeError("expansion leg")
        self.t = {}
        s = self.cl
        self.t["A1"] = brentq(lambda t: s.sol(t)[0] - self.l0m - (self.L - 0.5), s.t[-1], s.t[0], xtol=1e-13 * abs(s.t[-1]))
        self.t["A2"] = brentq(lambda t: s.sol(t)[0] - self.l0m - 2.0, s.t[-1], s.t[0], xtol=1e-13 * abs(s.t[-1]))
        s = self.el
        self.t["B"] = brentq(lambda t: s.sol(t)[0] - self.l0p - 2.0, s.t[0], s.t[-1], xtol=1e-13 * abs(s.t[-1]))
        for name, target in (("S4", SWITCH), ("S5", SWITCH_G4)):
            self.t[name] = brentq(lambda t: math.log(self.rr(s.sol(t))) - math.log(target), self.t["B"], s.t[-1], xtol=1e-13 * abs(s.t[-1]))

    def rr(self, y):
        q = self.Qp(list(y)); return q["rho"] / (q["rho"] + q["K"] + q["V"])

    def state(self, t):
        if t < self.tm: return -1, list(self.cl.sol(t)), self.Qm
        if t > self.tp: return 1, list(self.el.sol(t)), self.Qp
        return 0, list(self.tr.y(t)), None

    def bg(self, t):
        m, y, _ = self.state(t)
        return self.tr.bg(t) if m == 0 else DG.bg_mode(self.cu, self.Crad, y, m)

    def constraint_max(self):
        r = 0.0
        for s, Q in ((self.cl, self.Qm), (self.el, self.Qp)):
            for i in range(s.y.shape[1]):
                q = Q(list(s.y[:, i])); tot = q["rho"] + q["K"] + q["V"]
                r = max(r, abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"]))
        return r


def in_modes(A, bg, n, t0, toward):
    """in-modes in z = (q, P) at t0 (moved toward the bounce in 0.1 e-fold steps if E' is not positive); returns (t, Z, info)"""
    t = t0; shifts = 0
    while True:
        e = B62.end_info(A, bg, n, t, toward)
        if e["pos"]: break
        shifts += 1
        if shifts > 30: raise RuntimeError("no positive E' for the in-state")
        m, y, _ = bg.state(t); s = bg.cl if t < bg.tm else bg.el
        tgt = y[0] - 0.1                                                   # toward the bounce: smaller a on both legs
        lo, hi = sorted((t, bg.tm if t < bg.tm else bg.tp))
        t = brentq(lambda tt: s.sol(tt)[0] - tgt, lo, hi, xtol=1e-13 * abs(t))
    Z = np.linalg.solve(e["Sh"], e["F"])                                  # xi' -> z
    return t, Z, dict(t=t, shifts=shifts, adiab=e["adiab"], w_over_k=e["w_over_k"], k_over_aH=e["k_over_aH"])


def to_zetaC(A, bg, n, t, Zs):
    """two-field z = (q, P) at t (columns: modes) -> (zeta_C, d zeta_C / dt) for each mode"""
    b = bg.bg(t); K, Bm, _, _ = U.red2(A, b, n)
    m, y, Q = bg.state(t); q = Q(y); H, Hd, p = q["H"], q["Hdot"], y[2]
    pd = -3 * H * p - EC.dVp(y[1])
    out = []
    for j in range(Zs.shape[1]):
        qv, P = Zs[:2, j], Zs[2:, j]; qd = np.linalg.solve(K, P - Bm @ qv)
        zc = qv[0] - H / p * qv[1]
        zcd = qd[0] - H / p * qd[1] - (Hd / p - H * pd / p ** 2) * qv[1]
        out.append((zc, zcd))
    return out, dict(a=b["a"], H=H, chi=y[1], chiN=p / H, chidot=p)


def stage2_bg(model, chi, chiN, N0, span=25.0, npts=20001):
    s = solve_ivp(model.rhs, (N0, N0 + span), [chi, chiN], method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    Ns = np.linspace(N0, s.t[-1], npts); Qs = [model.quantities(s.sol(N)) for N in Ns]
    g = lambda k: np.array([q[k] for q in Qs])
    sp = {k: CubicSpline(Ns, g(k)) for k in ("alpha", "sigma", "eps")}; lnH = CubicSpline(Ns, np.log(g("H")))

    def f(N):
        ep, al = float(sp["eps"](N)), float(sp["alpha"](N))
        return dict(aH=math.exp(N + float(lnH(N))), H=math.exp(float(lnH(N))), alpha=al, sigma=float(sp["sigma"](N)), eps=ep,
                    eta=float(sp["eps"](N, 1)) / ep, beta=float(sp["alpha"](N, 1)) / al)
    f.N0, f.N1, f.sol = N0, float(s.t[-1]), s.sol
    f.lnaH = lambda N: N + float(lnH(N))
    return f


def stage2(bg2, k, N0, modes, rtol):
    """modes: list of (zeta_C, d zeta_C/dt) at N0 (complex); BKG equations in N until k/aH = 1e-3"""
    Nf = brentq(lambda N: bg2.lnaH(N) - math.log(1000 * k), N0, bg2.N1)
    def rhs(N, y):
        b = bg2(N); z2, c2 = U58.z2_cs2(b, k, N, False); cH = b["aH"]
        zr, zi, pr, pi_ = y
        return [pr / (z2 * cH), pi_ / (z2 * cH), -z2 * c2 * k * k * zr / cH, -z2 * c2 * k * k * zi / cH]
    b0 = bg2(N0); z20, _ = U58.z2_cs2(b0, k, N0, False); a0 = math.exp(N0)
    res = []; W0 = 0.0; W1 = 0.0
    for zc, zcd in modes:
        Pi = z20 * a0 * zcd                                                  # Pi = z^2 zeta' (conformal)
        W0 += (zc * np.conj(Pi)).imag
        s = solve_ivp(rhs, (N0, Nf), [zc.real, zc.imag, Pi.real, Pi.imag], method="DOP853", rtol=rtol, atol=1e-300)
        zf = s.y[0, -1] + 1j * s.y[1, -1]; Pf = s.y[2, -1] + 1j * s.y[3, -1]
        W1 += (zf * np.conj(Pf)).imag; res.append(abs(zf) ** 2)
    return float(sum(res)), dict(W_start=float(W0), W_end=float(W1), N_end=Nf)


def last_exit(bg, bg2, k, ts):
    """chi where k = aH for the last time (stage 2 if the mode is inside at the switch, else on the expansion leg)"""
    if bg2.lnaH(bg2.N0) < math.log(k):
        N = brentq(lambda N: bg2.lnaH(N) - math.log(k), bg2.N0, bg2.N1); return float(bg2.sol(N)[0]), "etapa2"
    ts_ = np.linspace(bg.tp, ts, 20001)[::-1]
    f = lambda t: (lambda m_y_Q: math.log(math.exp(m_y_Q[1][0]) * abs(m_y_Q[2](m_y_Q[1])["H"])) - math.log(k))(bg.state(t))
    prev = f(ts_[0])
    for t0, t1 in zip(ts_[:-1], ts_[1:]):
        cur = f(t1)
        if cur * prev <= 0 or cur < 0 <= prev:
            t = brentq(f, t1, t0); return float(bg.state(t)[1][1]), "etapa1"
        prev = cur
    return None, "sin_salida"


class Attractor:
    def __init__(self):
        self.model = U58.Model(U.ADOPTED); self.s, self.Nend = U58.background(self.model, 1e-10)
        self.bg = U58.make_bg(self.model, self.s, self.Nend)

    def P(self, chi_x, rtol=1e-10):
        chi_x = abs(chi_x)                                                   # the plateau is even: a child on the chi < 0 sheet uses |chi|
        N = brentq(lambda N: self.s.sol(N)[0] - chi_x, 0.0, self.Nend - 1)
        if N - math.log(100.0) < 0: return None, N
        k = self.bg(N)["aH"]; return U58.mode(self.bg, k, N, rtol)[0], N


def spectrum_row(A, bg, att, n, rtol=1e-10, switches=("S4",), states=("A1", "A2", "B")):
    k = math.sqrt(n * (n + 2.0)); out = dict(n=n)
    ins = {}
    for X in states:
        toward = 1 if X in ("A1", "A2") else -1
        ins[X] = in_modes(A, bg, n, bg.t[X], toward)
    # stage-1 segments from the earliest in-time to each switch
    marks = sorted({ins[X][0] for X in ins} | {bg.tm, bg.tp} | {bg.t[s] for s in switches})
    T = {}; Tcum = np.eye(4); cum = {marks[0]: np.eye(4)}
    for a_, b_ in zip(marks[:-1], marks[1:]):
        Tz, _ = B62.segment(A, bg, n, a_, b_, rtol); Tcum = Tz @ Tcum; cum[b_] = Tcum.copy()
    model = U58.Model(U.ADOPTED)
    for sw in switches:
        ts = bg.t[sw]
        res_sw = {}
        for X, (t0, Z, info) in ins.items():
            Tx = cum[ts] @ np.linalg.inv(cum[t0])
            modes, st = to_zetaC(A, bg, n, ts, Tx @ Z)
            N0 = math.log(st["a"])
            if "bg2" not in res_sw:
                res_sw["bg2"] = stage2_bg(model, st["chi"], st["chiN"], N0)
            amp2, w = stage2(res_sw["bg2"], k, N0, modes, rtol)
            res_sw[X] = dict(P=n ** 3 * amp2 / (2 * math.pi ** 2), **w, state=info if sw == switches[0] else None)
        chi_x, where = last_exit(bg, res_sw["bg2"], k, ts)
        Patt, Natt = att.P(chi_x) if chi_x is not None else (None, None)
        res_sw.pop("bg2")
        out[sw] = dict(chi_exit=chi_x, exit_in=where, P_att=Patt, N_att=Natt, **{X: res_sw[X] for X in ins})
        if Patt:
            for X in ins: out[sw]["F_tot_" + X] = res_sw[X]["P"] / Patt
        if "B" in ins:
            for X in ins:
                if X != "B": out[sw]["F_reb_" + X] = res_sw[X]["P"] / res_sw["B"]["P"]
    return out


def host_job(h):
    t0 = time.time(); A = U.Action(); att = Attractor(); par = h["par"]; bg = Background(par)
    s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))); n2, n40 = l2n(2, s), l2n(40, s)
    grid = sorted(set(int(round(x)) for x in np.logspace(math.log10(n2), math.log10(n40), 16)))
    nc = int(round(10 * h["n_mem"]))
    tag = "AH%g_f%.4f" % (par["AH"], par["f"]); log = open(RUNS / ("produce_%s.out" % tag), "a")
    rows = []
    for n in grid:
        r = spectrum_row(A, bg, att, n); r["ell"] = n * s; rows.append(r)
        log.write(json.dumps(dict(n=n, ell=n * s, F_A1=r["S4"].get("F_tot_A1"), F_A2=r["S4"].get("F_tot_A2"), F_B=r["S4"].get("F_tot_B"),
                                  chi_x=r["S4"]["chi_exit"], W=[r["S4"][X]["W_end"] for X in ("A1", "A2", "B")], sec=time.time() - t0)) + "\n"); log.flush()
    rc = spectrum_row(A, bg, att, nc, states=("A2", "B")); rc["ell"] = nc * s
    log.write("control " + json.dumps(dict(n=nc, F_A2=rc["S4"].get("F_tot_A2"), sec=time.time() - t0)) + "\n"); log.close()
    return dict(par=par, OmegaK0=h["OmegaK0"], n_mem=h["n_mem"], ell_mem=h["ell_mem"], L_contraction=bg.L, t=bg.t,
                constraint_max=bg.constraint_max(), grid=grid, rows=rows, control=rc, seconds=time.time() - t0)


def l2n(l, s): return l / s


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True); hosts = M.visible_hosts()
    with Pool(min(4, len(hosts))) as pool:
        res = pool.map(host_job, hosts, chunksize=1)
    out = dict(spec="docs/preregistros/ESPEC_U62c.md", hosts=res)
    save("u62c_espectro_rebote.json", out)
    # G2 and G4 on the first host
    A = U.Action(); att = Attractor(); h0 = hosts[0]; bg = Background(h0["par"]); r0 = res[0]
    s = math.sin(M.D_LSS * math.sqrt(abs(h0["OmegaK0"])))
    pick = [min(r0["grid"], key=lambda n: abs(n - l / s)) for l in (2, 10, 40)]
    g2, g4 = {}, {}
    for n in pick:
        base = next(r for r in r0["rows"] if r["n"] == n)["S4"]["A1"]["P"]
        a9 = spectrum_row(A, bg, att, n, rtol=1e-9, states=("A1",))["S4"]["A1"]["P"]
        s5 = spectrum_row(A, bg, att, n, switches=("S5",), states=("A1",))["S5"]["A1"]["P"]
        g2[n] = abs(a9 / base - 1); g4[n] = abs(s5 / base - 1)
        print("G2/G4", n, g2[n], g4[n], flush=True)
    allrows = [r for H in res for r in H["rows"]] + [H["control"] for H in res]
    W = [abs(r["S4"][X]["W_end"] / 0.25 - 1) for r in allrows for X in ("A1", "A2", "B") if X in r["S4"]]
    g0 = None  # the in-mode normalisation is checked inside B62.end_info / M.normal_modes; recomputed here for the record
    G0v = []
    for H in res[:1]:
        for n in H["grid"][:3]:
            for X, tow in (("A1", 1), ("A2", 1), ("B", -1)):
                e = B62.end_info(A, bg, n, bg.t[X], tow)
                if e["pos"]:
                    F = e["F"]; G0v.append(max(abs(M.Om(F[:, j], F[:, kk].conj()) - (1j if j == kk else 0)) for j in range(2) for kk in range(2)))
    ctrl = [abs(H["control"]["S4"]["F_tot_A2"] - 1) if H["control"]["S4"].get("F_tot_A2") else math.inf for H in res]
    gates = dict(G0=bool(max(G0v) <= 1e-10), G1=bool(max(W) <= 1e-3), G2=bool(max(g2.values()) <= 1e-3), G3=bool(max(ctrl) <= 0.03),
                 G4=bool(max(g4.values()) <= 1e-2), G5=bool(max(H["constraint_max"] for H in res) <= 3e-5))
    D = []
    for H in res:
        d = {}
        for X in ("A1", "A2"):
            v = [abs(r["S4"]["F_tot_" + X] - 1) for r in H["rows"] if r["S4"].get("F_tot_" + X)]
            v10 = [abs(r["S4"]["F_tot_" + X] - 1) for r in H["rows"] if r["S4"].get("F_tot_" + X) and r["ell"] >= 10]
            d[X] = max(v) if v else math.inf; d[X + "_l10"] = max(v10) if v10 else None
        d["dep_A1_A2_median"] = float(np.median([abs(r["S4"]["A1"]["P"] / r["S4"]["A2"]["P"] - 1) for r in H["rows"]]))
        D.append(d)
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif all(d["A1"] <= 0.1 and d["A2"] <= 0.1 for d in D): verdict = "ESPECTRO_LISO"
    elif sum(min(d["A1"], d["A2"]) >= 0.5 for d in D) >= 3: verdict = "RASGOS_FUERTES"
    else: verdict = "RASGOS_MODERADOS"
    out.update(gates=gates, G0_max=max(G0v), G1_max=max(W), G2=g2, G3_control=ctrl, G4=g4, constraint_max=max(H["constraint_max"] for H in res),
               D=D, verdict=verdict, seconds=time.time() - t0)
    save("u62c_espectro_rebote.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_max", "G1_max", "G2", "G3_control", "G4", "constraint_max", "D", "verdict", "seconds")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
