"""U62 (docs/preregistros/ESPEC_U62.md, amendment 1): Bogoliubov coefficients of the two scalar modes through the cuscuton
bounce window phi in [-0.2, 0.2], adopted set, for the five visible children of U61b.  In-state: the instantaneous vacuum of the
frozen Hamiltonian at phi = -0.2; out: decomposition on the instantaneous normal modes at phi = +0.2.  The U59 machinery
(action, reduced matrices in the uniform-radiation gauge, trajectory) is imported unchanged.
    python3 -B chain/u62_bogoliubov_rebote.py   ->  results/u62_bogoliubov_rebote.json
Internal, not exported."""
import json, math, time
import numpy as np
from scipy.integrate import solve_ivp
from ucommon import ROOT, save
import u59_escalares_rebote as U
from u61_cuscuton_huella import D_LSS

NS = (100, 200, 500, 1000, 2000, 5000, 10000, 30000, 100000)
ELLS = (2, 10, 30)
J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
Om = lambda f, g: f @ J @ g                                             # symplectic product, no conjugation


def visible_hosts():
    d = json.loads((ROOT / "results" / "u61b_cuscuton_huella.json").read_text())
    return [h for h in d["hosts"] if h.get("OmegaK0") is not None and h["OmegaK0"] >= -0.0031 and h.get("ell_mem") and h["ell_mem"] >= 2]


def normal_modes(E):
    """positive-frequency normal modes of the frozen Hamiltonian H = 1/2 xi^T E xi (xi canonical): columns f_j with J E f = -i w f,
    Om(f_j, f_k^*) = i delta_jk; sorted by decreasing frequency"""
    ev, V = np.linalg.eig(J @ E)
    sel = [i for i in range(4) if ev[i].imag < 0]
    w = -ev[sel].imag; F = V[:, sel]
    for j in range(2):
        c = (Om(F[:, j], F[:, j].conj()) / 1j).real
        F[:, j] /= math.sqrt(c)
    o = np.argsort(-w)
    return w[o], F[:, o]


def transfer_z(A, tr, n, rtol):
    """the U59 transfer in z = (q, P = p / a^3), returned as a matrix (U59.transfer only returns the gain)"""
    def rhs(t, y):
        bg = tr.bg(t); K, B, C, _ = U.red2(A, bg, n); a3 = bg["a"] ** 3
        Y = y.reshape(4, 4); out = np.empty_like(Y)
        for j in range(4):
            qv, p = Y[j, :2], Y[j, 2:]
            qd = np.linalg.solve(K, p / a3 - B @ qv)
            out[j, :2] = qd; out[j, 2:] = a3 * (B.T @ qd + C @ qv)
        return out.ravel()
    a3i, a3f = tr.bg(tr.tm)["a"] ** 3, tr.bg(tr.tp)["a"] ** 3
    Y0 = np.eye(4); Y0[:, 2:] *= a3i
    s = solve_ivp(rhs, (tr.tm, tr.tp), Y0.ravel(), method="DOP853", rtol=rtol, atol=1e-300)
    Zf = s.y[:, -1].reshape(4, 4).copy(); Zf[:, 2:] /= a3f
    return Zf.T, int(s.nfev)


def bogoliubov(A, tr, n, rtol):
    bgi, bgf = tr.bg(tr.tm), tr.bg(tr.tp); ai, af = bgi["a"], bgf["a"]
    Ei = U.energy_matrix(*U.red2(A, bgi, n)[:3]); Ef = U.energy_matrix(*U.red2(A, bgf, n)[:3])
    Ei, Ef = 0.5 * (Ei + Ei.T), 0.5 * (Ef + Ef.T)
    wi, Fi = normal_modes(Ei); wf, Ff = normal_modes(Ef)
    g0 = max(max(abs(Om(Fi[:, j], Fi[:, k].conj()) - (1j if j == k else 0)) for j in range(2) for k in range(2)),
             abs(Om(Fi[:, 0], Fi[:, 1])))
    Tz, nfev = transfer_z(A, tr, n, rtol)
    Tp = (af / ai) ** 1.5 * Tz                                           # map in the rescaled canonical variables (a^1.5 q, a^-1.5 p)
    out = Tp @ Fi                                                        # evolved in-modes (columns j)
    al = np.array([[1j * Om(Ff[:, k].conj(), out[:, j]) for k in range(2)] for j in range(2)])
    be = np.array([[-1j * Om(Ff[:, k], out[:, j]) for k in range(2)] for j in range(2)])
    norm = [float(np.sum(np.abs(al[j]) ** 2 - np.abs(be[j]) ** 2)) for j in range(2)]
    N = [float(np.sum(np.abs(be[:, k]) ** 2)) for k in range(2)]
    Rmin = float(np.sum((np.abs(al[:, 0]) - np.abs(be[:, 0])) ** 2)); Rmax = float(np.sum((np.abs(al[:, 0]) + np.abs(be[:, 0])) ** 2))
    kc = lambda a_: math.sqrt(n * (n + 2.0)) / a_
    return dict(n=n, N_chi=N[0], N_other=N[1], R_min=Rmin, R_max=Rmax, G0=float(g0), G1=float(max(abs(x - 1) for x in norm)),
                w_in_over_k=[float(x / kc(ai)) for x in wi], w_out_over_k=[float(x / kc(af)) for x in wf],
                k_phys_bounce=kc(math.exp(tr.y(0.0)[0])), nfev=nfev,
                abs_alpha=np.abs(al).tolist(), abs_beta=np.abs(be).tolist())


def main():
    t0 = time.time(); A = U.Action(); hosts = visible_hosts()
    out = dict(spec="docs/preregistros/ESPEC_U62.md", hosts=[])
    print(len(hosts), "visible hosts", flush=True)
    for h in hosts:
        par = h["par"]; cu, Crad, yb, _ = U.host_state(par); tr = U.Trajectory(cu, Crad, yb)
        s = math.sin(D_LSS * math.sqrt(abs(h["OmegaK0"])))
        nobs = {l: int(round(l / s)) for l in ELLS}
        rows = []
        for n in sorted(set(NS) | set(nobs.values())):
            r = bogoliubov(A, tr, n, 1e-10); rows.append(r)
            print(par["AH"], "%.3f" % par["f"], n, "N_chi %.4g N_o %.4g R %.3g-%.3g G1 %.1e w_in %s w_out %s" %
                  (r["N_chi"], r["N_other"], r["R_min"], r["R_max"], r["G1"], np.round(r["w_in_over_k"], 4), np.round(r["w_out_over_k"], 4)),
                  flush=True)
        out["hosts"].append(dict(par=par, OmegaK0=h["OmegaK0"], n_mem=h["n_mem"], ell_mem=h["ell_mem"], N_tot=h["N_tot"],
                                 a_b=math.exp(yb[0]), window_t=[tr.tm, tr.tp], n_obs=nobs, rows=rows))
        out["hosts"][-1]["_tr"] = tr
    # G2: rtol 1e-9 against 1e-10 at three n, first host
    h0 = out["hosts"][0]; conv = {}
    for n in (100, 10000, 100000):
        a_ = bogoliubov(A, h0["_tr"], n, 1e-9); b_ = next(r for r in h0["rows"] if r["n"] == n)
        conv[n] = max(abs(a_[k] - b_[k]) / (1e-4 + 1e-4 * abs(b_[k])) for k in ("N_chi", "N_other"))
    for h in out["hosts"]: h.pop("_tr")
    rowsall = [r for h in out["hosts"] for r in h["rows"]]
    gates = dict(G0=bool(max(r["G0"] for r in rowsall) <= 1e-10), G1=bool(max(r["G1"] for r in rowsall) <= 1e-6),
                 G2=bool(max(conv.values()) <= 1.0))
    obs = [[next(r for r in h["rows"] if r["n"] == h["n_obs"][l])["N_chi"] for l in ELLS] for h in out["hosts"]]
    excited = sum(max(v) >= 0.1 for v in obs); strict = sum(min(v) >= 0.1 for v in obs)
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif all(x <= 0.01 for v in obs for x in v): verdict = "REBOTE_ADIABATICO"
    elif excited >= 3: verdict = "REBOTE_EXCITA"
    else: verdict = "INTERMEDIO"
    out.update(gates=gates, G0_max=max(r["G0"] for r in rowsall), G1_max=max(r["G1"] for r in rowsall), G2_conv_over_tol=conv,
               N_chi_obs=obs, n_excited=excited, n_excited_strict=strict, verdict=verdict, seconds=time.time() - t0)
    save("u62_bogoliubov_rebote.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_max", "G1_max", "G2_conv_over_tol", "N_chi_obs", "n_excited", "n_excited_strict",
                                          "verdict", "seconds")}, indent=1, default=float))


if __name__ == "__main__":
    main()
