"""U62b (docs/preregistros/ESPEC_U62b.md): Bogoliubov coefficients of the physical scalar modes through the cuscuton bounce.
The in/out states are the instantaneous vacua of the frozen Hamiltonian in the frame where the symmetric part of B is integrated by
parts (p' = p - a^3 B_s q), set Delta ln a = 2 (and 1.5, for the stability check) away from the window phi in [-0.2, 0.2], on a
contraction leg (mode -1 of explora_cuscuton, integrated backwards) and an expansion leg (mode +1).  The evolution is that of U59
(uniform-radiation gauge, q = (zeta, x)); U59, U62 and the U62 diagnostics are imported unchanged.
    python3 -B chain/u62b_bogoliubov_fisico.py   ->  results/u62b_bogoliubov_fisico.json
Internal, not exported."""
import json, math, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT, save
import explora_cuscuton as EC
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M
import u62_diagnostico as DG

NS = M.NS
ELLS = M.ELLS
DELTAS = (2.0, 1.5)
ATOL_BG = [1e-14, 1e-15, 1e-24, 1e-24]
Om = M.Om


class Background:
    """window (U59.Trajectory) plus a contraction leg (mode -1, backwards) and an expansion leg (mode +1), Delta ln a up to 2.05"""

    def __init__(self, par):
        self.par = par
        cu, Crad, yb, _ = U.host_state(par); self.cu, self.Crad = cu, Crad
        self.tr = tr = U.Trajectory(cu, Crad, yb); self.tm, self.tp = tr.tm, tr.tp
        self.legs = {}
        for sgn, mode, t0 in ((-1, -1, tr.tm), (1, 1, tr.tp)):
            y0 = list(tr.y(t0)); phi0 = y0[3]
            y0[3] = cu.pinf + phi0 if mode == -1 else cu.pinf - phi0           # d for modes -1 / +1
            _, rhs = EC.deriv(cu, mode, Crad, 0.0)
            ev = (lambda l0: (lambda t, y: y[0] - l0 - 2.05))(y0[0]); ev.terminal = True
            s = solve_ivp(rhs, (t0, t0 + sgn * 1e13), y0, method="DOP853", rtol=1e-12, atol=ATOL_BG, dense_output=True, events=ev)
            if s.status != 1:
                raise RuntimeError("leg %d did not reach Delta ln a = 2.05" % sgn)
            self.legs[sgn] = (mode, s, y0[0])
        self.t_at = {}
        for d in DELTAS:
            for sgn in (-1, 1):
                mode, s, l0 = self.legs[sgn]
                lo, hi = sorted((s.t[0], s.t[-1]))
                self.t_at[(sgn, d)] = brentq(lambda t: s.sol(t)[0] - l0 - d, lo, hi, xtol=1e-13 * abs(hi))

    def y_mode(self, t):
        if t < self.tm: m, s, _ = self.legs[-1]; return m, list(s.sol(t))
        if t > self.tp: m, s, _ = self.legs[1]; return m, list(s.sol(t))
        return 0, list(self.tr.y(t))

    def bg(self, t):
        m, y = self.y_mode(t)
        return self.tr.bg(t) if m == 0 else DG.bg_mode(self.cu, self.Crad, y, m)

    def constraint_max(self):
        r = 0.0
        for sgn in (-1, 1):
            mode, s, _ = self.legs[sgn]; Q, _ = EC.deriv(self.cu, mode, self.Crad, 0.0)
            for i in range(s.y.shape[1]):
                q = Q(list(s.y[:, i])); tot = q["rho"] + q["K"] + q["V"]
                r = max(r, abs(3 * (q["H"] ** 2 + q["ia2"]) - tot - q["U"]) / (tot + abs(q["U"]) + 3 * q["ia2"]))
        return r


def segment(A, bg, n, t0, t1, rtol):
    """transfer in z = (q, P = p / a^3) from t0 to t1 (U59 equations)"""
    def rhs(t, y):
        b = bg.bg(t); K, B, C, _ = U.red2(A, b, n); a3 = b["a"] ** 3
        Y = y.reshape(4, 4); out = np.empty_like(Y)
        for j in range(4):
            qv, p = Y[j, :2], Y[j, 2:]
            qd = np.linalg.solve(K, p / a3 - B @ qv)
            out[j, :2] = qd; out[j, 2:] = a3 * (B.T @ qd + C @ qv)
        return out.ravel()
    a3i, a3f = bg.bg(t0)["a"] ** 3, bg.bg(t1)["a"] ** 3
    Y0 = np.eye(4); Y0[:, 2:] *= a3i
    s = solve_ivp(rhs, (t0, t1), Y0.ravel(), method="DOP853", rtol=rtol, atol=1e-300)
    Zf = s.y[:, -1].reshape(4, 4).copy(); Zf[:, 2:] /= a3f
    return Zf.T, int(s.nfev)


def frame(A, bg, n, t, inward):
    """E' and the map z -> xi' = a^1.5 (q, P - B_s q) at t; Bdot by central differences (relative step 1e-5)"""
    b = bg.bg(t); K, B, C, _ = U.red2(A, b, n)
    h = 1e-5 * abs(t)
    Bd = (U.red2(A, bg.bg(t + h), n)[1] - U.red2(A, bg.bg(t - h), n)[1]) / (2 * h)
    Bs, Ba = 0.5 * (B + B.T), 0.5 * (B - B.T)
    E = U.energy_matrix(K, Ba, C - 0.5 * (Bd + Bd.T) - 3 * b["H"] * Bs); E = 0.5 * (E + E.T)
    Sh = b["a"] ** 1.5 * np.block([[np.eye(2), np.zeros((2, 2))], [-Bs, np.eye(2)]])
    return E, Sh, b


def end_info(A, bg, n, t, inward):
    E, Sh, b = frame(A, bg, n, t, inward)
    k = math.sqrt(n * (n + 2.0)) / b["a"]
    if np.min(np.linalg.eigvalsh(E)) <= 0:
        return dict(pos=False, k_over_aH=k / abs(b["H"]))
    w, F = M.normal_modes(E)
    t2 = t + inward * 1e-3 * abs(t); E2 = frame(A, bg, n, t2, inward)[0]
    ad = (np.abs(M.normal_modes(E2)[0] - w) / abs(t2 - t) / w ** 2).tolist() if np.min(np.linalg.eigvalsh(E2)) > 0 else [math.inf, math.inf]
    return dict(pos=True, E=E, Sh=Sh, w=w, F=F, w_over_k=(w / k).tolist(), adiab=ad, k_over_aH=k / abs(b["H"]))


def bogo(ei, eo, Tz):
    Tp = eo["Sh"] @ Tz @ np.linalg.inv(ei["Sh"])
    Fi, Ff = ei["F"], eo["F"]
    g0 = max(max(abs(Om(Fi[:, j], Fi[:, k].conj()) - (1j if j == k else 0)) for j in range(2) for k in range(2)), abs(Om(Fi[:, 0], Fi[:, 1])))
    out = Tp @ Fi
    al = np.array([[1j * Om(Ff[:, k].conj(), out[:, j]) for k in range(2)] for j in range(2)])
    be = np.array([[-1j * Om(Ff[:, k], out[:, j]) for k in range(2)] for j in range(2)])
    norm = [float(np.sum(np.abs(al[j]) ** 2 - np.abs(be[j]) ** 2)) for j in range(2)]
    ichi = int(np.argmin([abs(x - 1) for x in eo["w_over_k"]])); irad = 1 - ichi
    N = [float(np.sum(np.abs(be[:, k]) ** 2)) for k in range(2)]
    return dict(N_chi=N[ichi], N_rad=N[irad], i_chi=ichi,
                R_min=float(np.sum((np.abs(al[:, ichi]) - np.abs(be[:, ichi])) ** 2)),
                R_max=float(np.sum((np.abs(al[:, ichi]) + np.abs(be[:, ichi])) ** 2)),
                G0=float(g0), G1=float(max(abs(x - 1) for x in norm)), abs_alpha=np.abs(al).tolist(), abs_beta=np.abs(be).tolist())


def row(A, bg, n, rtol):
    ta2, ta15 = bg.t_at[(-1, 2.0)], bg.t_at[(-1, 1.5)]; tb15, tb2 = bg.t_at[(1, 1.5)], bg.t_at[(1, 2.0)]
    ends = {d: (end_info(A, bg, n, bg.t_at[(-1, d)], 1), end_info(A, bg, n, bg.t_at[(1, d)], -1)) for d in DELTAS}
    strip = lambda e: {k: v for k, v in e.items() if k in ("pos", "w_over_k", "adiab", "k_over_aH")}
    r = dict(n=n, ends={str(d): [strip(e) for e in ends[d]] for d in DELTAS})
    if not all(e["pos"] for d in DELTAS for e in ends[d]):
        r["status"] = "NO_DEFINIDO: E' no positiva en un extremo"; return r
    nf = 0; T = {}
    for a_, b_ in ((ta2, ta15), (ta15, bg.tm), (bg.tm, bg.tp), (bg.tp, tb15), (tb15, tb2)):
        T[(a_, b_)], k_ = segment(A, bg, n, a_, b_, rtol); nf += k_
    core = T[(bg.tp, tb15)] @ T[(bg.tm, bg.tp)] @ T[(ta15, bg.tm)]
    full = T[(tb15, tb2)] @ core @ T[(ta2, ta15)]
    r["b2"] = bogo(ends[2.0][0], ends[2.0][1], full); r["b15"] = bogo(ends[1.5][0], ends[1.5][1], core); r["nfev"] = nf
    ad_ok = all(max(e["adiab"]) <= 0.1 for e in ends[2.0])
    st_ok = abs(r["b2"]["N_chi"] - r["b15"]["N_chi"]) <= max(0.005, 0.25 * r["b2"]["N_chi"])
    r["status"] = "DEFINIDO" if (ad_ok and st_ok) else ("NO_DEFINIDO: no adiabatico" if not ad_ok else "NO_DEFINIDO: N inestable entre 1.5 y 2")
    return r


def main():
    t0 = time.time(); A = U.Action(); hosts = M.visible_hosts()
    out = dict(spec="docs/preregistros/ESPEC_U62b.md", hosts=[]); bgs = []
    for h in hosts:
        par = h["par"]; bg = Background(par); bgs.append(bg)
        s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))); nobs = {l: int(round(l / s)) for l in ELLS}
        H = dict(par=par, OmegaK0=h["OmegaK0"], n_obs=nobs, constraint_max_legs=bg.constraint_max(),
                 t_ends={str(d): [bg.t_at[(-1, d)], bg.t_at[(1, d)]] for d in DELTAS}, window_t=[bg.tm, bg.tp], rows=[])
        for n in sorted(set(NS) | set(nobs.values())):
            r = row(A, bg, n, 1e-10); H["rows"].append(r)
            b = r.get("b2")
            print(par["AH"], "%.3f" % par["f"], n, r["status"], ("N_chi %.4g (1.5: %.4g) N_rad %.4g G1 %.1e" %
                  (b["N_chi"], r["b15"]["N_chi"], b["N_rad"], b["G1"])) if b else "", "w/k", [np.round(e["w_over_k"], 3).tolist() if e["pos"] else None for e in r["ends"]["2.0"]],
                  "ad", [np.round(e["adiab"], 3).tolist() if e["pos"] else None for e in r["ends"]["2.0"]], flush=True)
        out["hosts"].append(H)
        save("u62b_bogoliubov_fisico.json", out)                                # partial, overwritten at the end
    # G2 on the first host
    conv = {}
    for n in (1000, 10000, 100000):
        a_ = row(A, bgs[0], n, 1e-9); b_ = next(r for r in out["hosts"][0]["rows"] if r["n"] == n)
        if "b2" in a_ and "b2" in b_:
            conv[n] = max(abs(a_["b2"][k] - b_["b2"][k]) / (1e-4 + 1e-4 * abs(b_["b2"][k])) for k in ("N_chi", "N_rad"))
    rows = [r for H in out["hosts"] for r in H["rows"] if "b2" in r]
    defined = [(H, r) for H in out["hosts"] for r in H["rows"] if r["status"] == "DEFINIDO"]
    lab_ok = all(all(0.95 <= e["w_over_k"][r["b2"]["i_chi"]] <= 1.05 and 0.95 <= e["w_over_k"][1 - r["b2"]["i_chi"]] * math.sqrt(3) <= 1.05
                     for e in r["ends"]["2.0"]) for _, r in defined)
    gates = dict(G0=bool(max(max(r["b2"]["G0"], r["b15"]["G0"]) for r in rows) <= 1e-10),
                 G1=bool(max(max(r["b2"]["G1"], r["b15"]["G1"]) for r in rows) <= 1e-6),
                 G2=bool(len(conv) > 0 and max(conv.values()) <= 1.0), G3=bool(lab_ok),
                 G4=bool(max(H["constraint_max_legs"] for H in out["hosts"]) <= 1e-6))
    obs = []
    for H in out["hosts"]:
        for l, n in H["n_obs"].items():
            r = next(r for r in H["rows"] if r["n"] == n)
            obs.append(dict(AH=H["par"]["AH"], f=H["par"]["f"], ell=int(l), n=n, status=r["status"],
                            N_chi=r["b2"]["N_chi"] if "b2" in r else None, N_rad=r["b2"]["N_rad"] if "b2" in r else None))
    d_obs = [o for o in obs if o["status"] == "DEFINIDO"]
    if not all(gates.values()) or len(d_obs) < 3: verdict = "INCONCLUSO"
    elif all(o["N_chi"] <= 0.01 for o in d_obs): verdict = "REBOTE_ADIABATICO"
    elif sum(o["N_chi"] >= 0.1 for o in d_obs) >= len(d_obs) / 3: verdict = "REBOTE_EXCITA"
    else: verdict = "INTERMEDIO"
    out.update(gates=gates, G0_max=max(max(r["b2"]["G0"], r["b15"]["G0"]) for r in rows), G1_max=max(max(r["b2"]["G1"], r["b15"]["G1"]) for r in rows),
               G2_conv_over_tol=conv, constraint_max=max(H["constraint_max_legs"] for H in out["hosts"]), obs=obs, n_defined_obs=len(d_obs),
               verdict=verdict, seconds=time.time() - t0)
    save("u62b_bogoliubov_fisico.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_max", "G1_max", "G2_conv_over_tol", "constraint_max", "obs", "n_defined_obs", "verdict", "seconds")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
