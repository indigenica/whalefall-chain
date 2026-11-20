"""U65 (docs/preregistros/ESPEC_U65.md): tensors (and, informatively, vectors) through the cuscuton bounce.  Tensor modes on S^3
(TT harmonics, k^2 + 2K = n(n+2)) obey the GR equation on the cuscuton background; canonical pair (h, pi = a^3 hdot / 4), vacuum
omega = k/a at A1 / A2 / B; stage 2 on the flat U58 background until k/aH = 1e-3; reference: U58 tensor Bunch-Davies on the
attractor at |chi_x|.  Background: U62d (forward from the turnaround).  U58, U62c, U62d imported unchanged.
    python3 -B chain/u65_tensores_vectores.py   ->  results/u65_tensores_vectores.json
Internal, not exported."""
import json, math, time
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, save
import explora_cuscuton as EC
import u58_cuscuton_inflacion as U58
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62c_espectro_rebote as C
import u62d_cl_rebote as D

RUNS = ROOT / "results" / "u65_runs"
NB = (100, 300, 1000, 3000, 10000)


def k2(n): return n * (n + 2.0)


def transfer(bg, n, t0, t1, rtol):
    """2x2 real transfer of (h, pi) from t0 to t1 on the child's background"""
    def rhs(t, y):
        a = bg.bg(t)["a"]; return [4 * y[1] / a ** 3, -a * k2(n) * y[0] / 4, 4 * y[3] / a ** 3, -a * k2(n) * y[2] / 4]
    s = solve_ivp(rhs, (t0, t1), [1.0, 0.0, 0.0, 1.0], method="DOP853", rtol=rtol, atol=1e-300)
    y = s.y[:, -1]; return np.array([[y[0], y[2]], [y[1], y[3]]])


def vacuum(bg, n, t):
    a = bg.bg(t)["a"]; w = math.sqrt(k2(n)) / a; A = math.sqrt(2 / (w * a ** 3))
    return np.array([A, -1j * w * a ** 3 / 4 * A])                        # (h, pi): Omega(f, f*) = i


Om = lambda f, g: f[0] * g[1] - f[1] * g[0]


def stage2(bg2, k, N0, h, hd, a, rtol):
    Nf = brentq(lambda N: bg2.lnaH(N) - math.log(1000 * k), N0, bg2.N1)
    def rhs(N, y):
        b = bg2(N); z2 = math.exp(2 * N); cH = b["aH"]
        return [y[2] / (z2 * cH), y[3] / (z2 * cH), -z2 * k * k * y[0] / cH, -z2 * k * k * y[1] / cH]
    Q = a ** 3 * hd
    s = solve_ivp(rhs, (N0, Nf), [h.real, h.imag, Q.real, Q.imag], method="DOP853", rtol=rtol, atol=1e-300)
    hf = s.y[0, -1] + 1j * s.y[1, -1]; Qf = s.y[2, -1] + 1j * s.y[3, -1]
    return abs(hf) ** 2, float((hf * np.conj(Qf)).imag), float((h * np.conj(Q)).imag)


class Attractor(D.Attractor):
    def Pt(self, chi_x, rtol=1e-10):
        x = abs(chi_x); lo, hi = float(self.s.sol(self.Nend - 1)[0]), float(self.s.sol(0.0)[0])
        if not (lo < x < hi): return None
        N = brentq(lambda N: self.s.sol(N)[0] - x, 0.0, self.Nend - 1)
        if N - math.log(100.0) < 0: return None
        return U58.mode(self.bg, self.bg(N)["aH"], N, rtol, tensor=True)[0]


def row(args):
    par, n, rtol, switches, states = args
    t_0 = time.time(); bg = D.Background(par); att = Attractor(); model = U58.Model(U.ADOPTED); k = math.sqrt(k2(n))
    out = dict(par=par, n=n, rtol=rtol)
    marks = sorted({bg.t[X] for X in states} | {bg.tm, bg.tp} | {bg.t[s] for s in switches})
    cum = {marks[0]: np.eye(2)}; T = np.eye(2)
    for a_, b_ in zip(marks[:-1], marks[1:]):
        T = transfer(bg, n, a_, b_, rtol) @ T; cum[b_] = T.copy()
    out["G0"] = max(abs(Om(vacuum(bg, n, bg.t[X]), np.conj(vacuum(bg, n, bg.t[X]))) - 1j) for X in states)
    for sw in switches:
        ts = bg.t[sw]; m, y, Q = bg.state(ts); q = Q(y); a = math.exp(y[0]); N0 = y[0]
        bg2 = C.stage2_bg(model, y[1], y[2] / q["H"], N0)
        res = {}
        for X in states:
            f0 = vacuum(bg, n, bg.t[X]); f = (cum[ts] @ np.linalg.inv(cum[bg.t[X]])) @ f0
            h, pi_ = f[0], f[1]; hd = 4 * pi_ / a ** 3
            amp2, W1, W0 = stage2(bg2, k, N0, h, hd, a, rtol)
            res[X] = dict(P=n ** 3 * amp2 / math.pi ** 2, W_end=W1, W_start=W0)
        chi_x, where = C.last_exit(bg, bg2, k, ts)
        Pa = att.Pt(chi_x) if chi_x is not None else None
        out[sw] = dict(chi_exit=chi_x, exit_in=where, P_att=Pa, **res)
        if Pa:
            for X in states: out[sw]["F_" + X] = res[X]["P"] / Pa
    out["seconds"] = time.time() - t_0
    return out


def bogo_row(args):
    par, n = args; bg = D.Background(par)
    f = transfer(bg, n, bg.t["A2"], bg.t["B"], 1e-10) @ vacuum(bg, n, bg.t["A2"]); g = vacuum(bg, n, bg.t["B"])
    beta = -1j * Om(g, f); alpha = 1j * Om(np.conj(g), f)
    return dict(par=par, n=n, N_t=float(abs(beta) ** 2), norm=float(abs(alpha) ** 2 - abs(beta) ** 2))


def cached(args):
    par, n, rtol, switches, states = args
    f = RUNS / "rows" / ("AH%g_f%.6f_n%d_rt%g_%s_%s.json" % (par["AH"], par["f"], n, rtol, "-".join(switches), "-".join(states)))
    f.parent.mkdir(parents=True, exist_ok=True)
    if f.exists(): return json.loads(f.read_text())
    r = row(args); f.write_text(json.dumps(r, default=float))
    with open(RUNS / "produce.out", "a") as log: log.write("row AH %g f %.4f n %d %.0fs\n" % (par["AH"], par["f"], n, r["seconds"]))
    return r


def vectors():
    d = json.loads((ROOT / "results" / "u61b_cuscuton_huella.json").read_text()); out = []
    for h in [h for h in d["hosts"] if "chi_land" in h]:
        par = dict(U.ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]); bg = D.Background(par)
        a_t = par["AH"] / EC.H_PLATEAU; a_b = math.exp(bg.wl.sol(0.0)[0]); g = (a_t / a_b) ** 2
        out.append(dict(AH=par["AH"], f=par["f"], a_t=a_t, a_b=a_b, growth=g, max_initial_PhiV=0.1 / g))
    return out


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    hosts = M.visible_hosts(); jobs = []; grids = {}
    for h in hosts:
        par = dict(U.ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]); s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"])))
        grid = sorted(set(int(round(x)) for x in np.logspace(math.log10(2 / s), math.log10(40 / s), 16)))
        nc = int(round(10 * h["n_mem"])); grids[(par["AH"], par["f"])] = (grid, nc, s)
        jobs += [(par, n, 1e-10, ("S4",), ("A1", "A2", "B")) for n in grid] + [(par, nc, 1e-10, ("S4",), ("A2", "B"))]
    with Pool(4) as p:
        R = p.map(cached, jobs, chunksize=1)
    print("rows done %.0fs" % (time.time() - t0), flush=True)
    par0 = jobs[0][0]; g0grid, _, s0 = grids[(par0["AH"], par0["f"])]
    pick = [min(g0grid, key=lambda n: abs(n - l / s0)) for l in (2, 10, 40)]
    extra = [(par0, n, 1e-9, ("S4",), ("A2",)) for n in pick] + [(par0, n, 1e-10, ("S5",), ("A2",)) for n in pick]
    with Pool(4) as p:
        X = p.map(cached, extra, chunksize=1)
    g2, g4 = {}, {}
    for n in pick:
        base = next(r for r in R if r["par"] == par0 and r["n"] == n)["S4"]["A2"]["P"]
        g2[n] = abs(next(x for x in X if x["n"] == n and x["rtol"] == 1e-9)["S4"]["A2"]["P"] / base - 1)
        g4[n] = abs(next(x for x in X if x["n"] == n and "S5" in x)["S5"]["A2"]["P"] / base - 1)
    with Pool(4) as p:
        BG = p.map(bogo_row, [(dict(U.ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]), n) for h in hosts for n in NB], chunksize=1)
    allr = R
    G0 = max(r["G0"] for r in allr); G1 = max(abs(r["S4"][X_]["W_end"] / 2 - 1) for r in allr for X_ in ("A1", "A2", "B") if X_ in r["S4"])
    ctrl = []; hostsout = []
    for h in hosts:
        par = dict(U.ADOPTED, AH=h["par"]["AH"], f=h["par"]["f"]); grid, nc, s = grids[(par["AH"], par["f"])]
        rows = sorted([r for r in R if r["par"] == par and r["n"] in grid], key=lambda r: r["n"]); rc = next(r for r in R if r["par"] == par and r["n"] == nc)
        ctrl.append(abs(rc["S4"]["F_A2"] - 1) if rc["S4"].get("F_A2") else math.inf)
        Dx = {X_: max([abs(r["S4"]["F_" + X_] - 1) for r in rows if r["S4"].get("F_" + X_)] or [math.inf]) for X_ in ("A1", "A2")}
        hostsout.append(dict(par=par, ell_mem=h["ell_mem"], n_mem=h["n_mem"], D=Dx, control_F=rc["S4"].get("F_A2"),
                             rows=[dict(n=r["n"], ell=r["n"] * s, F_A1=r["S4"].get("F_A1"), F_A2=r["S4"].get("F_A2"), F_B=r["S4"].get("F_B"),
                                        chi_exit=r["S4"]["chi_exit"]) for r in rows]))
    gates = dict(G0=bool(G0 <= 1e-10), G1=bool(G1 <= 1e-3), G2=bool(max(g2.values()) <= 1e-3), G3=bool(max(ctrl) <= 0.03), G4=bool(max(g4.values()) <= 1e-2))
    smooth = all(H_["D"]["A1"] <= 0.1 and H_["D"]["A2"] <= 0.1 for H_ in hostsout)
    feat = sum(min(H_["D"]["A1"], H_["D"]["A2"]) >= 0.5 for H_ in hostsout) >= 3
    would = "TENSORES_LISOS" if smooth else ("TENSORES_CON_RASGOS" if feat else "TENSORES_MODERADOS")
    verdict = would if all(gates.values()) else "INCONCLUSO"
    vec = vectors()
    out = dict(spec="docs/preregistros/ESPEC_U65.md", gates=gates, G0=G0, G1=G1, G2=g2, G3_control=ctrl, G4=g4, hosts=hostsout, bogoliubov=BG,
               vectors=vec, would_be=would, verdict=verdict, seconds=time.time() - t0)
    save("u65_tensores_vectores.json", out)
    print(json.dumps(dict(gates=gates, G0=G0, G1=G1, G2=g2, G3_control=ctrl, G4=g4, would_be=would, verdict=verdict), indent=1, default=float))
    for H_ in hostsout:
        print(H_["par"]["AH"], "%.3f" % H_["par"]["f"], "D", {k_: round(v, 3) for k_, v in H_["D"].items()}, "ctrl %.4f" % (H_["control_F"] or 0),
              [(round(r["ell"], 1), r["F_A1"] and round(r["F_A1"], 3), r["F_A2"] and round(r["F_A2"], 3), r["F_B"] and round(r["F_B"], 3)) for r in H_["rows"][::3]])
    print("Bogoliubov", [(b["par"]["AH"], b["n"], "%.2e" % b["N_t"]) for b in BG])
    gv = [v["growth"] for v in vec]; print("vectors growth min %.2e median %.2e max %.2e" % (min(gv), float(np.median(gv)), max(gv)))


if __name__ == "__main__":
    main()
