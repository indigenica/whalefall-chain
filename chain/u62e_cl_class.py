"""U62e (docs/preregistros/ESPEC_U62e.md): the U62d primordial spectra of the four visible cuscuton children (chi > 0 sheet; states
A1, A2) through the full CLASS C_l (flat, Planck 2018 base, external primordial table), against the Planck 2018 Commander TT
D_l at l = 2-29 with a full-sky Wishart likelihood (f_sky = 0.8).
    python3 -B chain/u62e_cl_class.py   ->  results/u62e_cl_class.json
Internal, not exported."""
import json, math, time, tempfile, os
import numpy as np
from classy import Class
from ucommon import ROOT, save

BASE = dict(h=0.6736, omega_b=0.02237, omega_cdm=0.1200, tau_reio=0.0544, N_ur=3.046,
            output="tCl", l_max_scalars=40, l_logstep=1.001, l_linstep=1)
AS, NS_, KP = 2.100e-9, 0.9649, 0.05
QLIN, QLIN_G1 = 0.02, 0.01
FSKY = 0.8
L = np.arange(2, 30)
TCMB_UK = 2.7255e6
H0_C = 0.6736 / 2997.92458                                              # H0/c [1/Mpc]
RUNS = ROOT / "results" / "u62e_runs"


def planck():
    d = np.loadtxt(ROOT / "data" / "planck2018" / "COM_PowerSpect_CMB-TT-full_R3.01.txt")
    return {int(r[0]): r[1] for r in d if 2 <= r[0] <= 29}


def run_class(table=None, qlin=QLIN):
    p = dict(BASE); p["q_linstep"] = qlin
    tmp = None
    if table is None:
        p.update(A_s=AS, n_s=NS_)
    else:
        tmp = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False, dir=RUNS)
        for k, P in table: tmp.write("%.12e %.12e\n" % (k, P))
        tmp.close(); p.update({"P_k_ini type": "external_Pk", "command": "cat " + tmp.name})
    c = Class(); c.set(p); c.compute(); cl = c.raw_cl(40)["tt"]; c.struct_cleanup(); c.empty()
    if tmp: os.unlink(tmp.name)
    return np.array([l * (l + 1) * cl[l] / (2 * math.pi) * TCMB_UK ** 2 for l in L])


def F_of_n(rows, X, every=1):
    rr = [r for r in rows if r["S4"].get("F_tot_" + X)][::every]
    ln_n = np.log([r["n"] for r in rr]); lnF = np.log([r["S4"]["F_tot_" + X] for r in rr])
    def F(n):
        ln = np.log(n)
        return np.where((ln >= ln_n[0]) & (ln <= ln_n[-1]), np.exp(np.interp(ln, ln_n, lnF)), 1.0)
    return F, int(round(math.exp(ln_n[0]))), int(round(math.exp(ln_n[-1])))


def table(H, F, nlo, nhi):
    sK = math.sqrt(abs(H["OmegaK0"])) * H0_C
    n_int = np.arange(max(nlo, 2), nhi + 1, dtype=float)
    k = np.unique(np.concatenate([np.logspace(-8, 1, 4000), np.sqrt(n_int * (n_int + 2)) * sK]))
    n_of_k = -1 + np.sqrt(1 + (k / sK) ** 2)
    P = AS * (k / KP) ** (NS_ - 1) * (F(np.maximum(n_of_k, 1.0)) if F else 1.0)
    return list(zip(k, P))


def m2lnL(D, Dhat):
    return float(FSKY * np.sum((2 * L + 1) * (Dhat / D + np.log(D))))


def main():
    t0 = time.time(); RUNS.mkdir(parents=True, exist_ok=True)
    d62 = json.loads((ROOT / "results" / "u62d_cl_rebote.json").read_text())
    Dp = np.array([planck()[l] for l in L])
    analytic = run_class()
    out = dict(spec="docs/preregistros/ESPEC_U62e.md", ell=L.tolist(), D_planck=Dp.tolist(), D_analytic=analytic.tolist(), hosts=[])
    g0 = []; g1 = []
    for H in d62["hosts"]:
        rows = H["rows"]; _, nlo, nhi = F_of_n(rows, "A2")
        ref = run_class(table(H, None, nlo, nhi)); ref01 = run_class(table(H, None, nlo, nhi), QLIN_G1)
        g0.append(float(np.max(np.abs(ref / analytic - 1))))
        res = dict(par=H["par"], OmegaK0=H["OmegaK0"], ell_mem=H["ell_mem"], D_ref=ref.tolist(), states={})
        for X in ("A1", "A2"):
            F, a, b = F_of_n(rows, X); Fh, ah, bh = F_of_n(rows, X, every=2)
            D = run_class(table(H, F, a, b)); Dh = run_class(table(H, Fh, ah, bh)); D01 = run_class(table(H, F, a, b), QLIN_G1)
            dP = m2lnL(D, Dp) - m2lnL(ref, Dp); dPh = m2lnL(Dh, Dp) - m2lnL(ref, Dp); dP01 = m2lnL(D01, Dp) - m2lnL(ref01, Dp)
            dL = m2lnL(D, ref) - m2lnL(ref, ref)
            low = lambda DD, RR, Dh_: float(FSKY * np.sum(((2 * L + 1) * (Dh_ / DD + np.log(DD)) - (2 * L + 1) * (Dh_ / RR + np.log(RR)))[:4]))
            defined = bool(abs(dPh - dP) <= max(2.0, 0.2 * abs(dP)))
            g1.append(abs(dP01 - dP) <= max(0.5, 0.05 * abs(dP)))
            res["states"][X] = dict(ratio=(D / ref).tolist(), dchi2_planck=dP, dchi2_planck_half=dPh, dchi2_planck_q001=dP01, dchi2_lcdm=dL,
                                    dchi2_planck_l2to5=low(D, ref, Dp), defined=defined)
            print(H["par"]["AH"], "%.3f" % H["par"]["f"], X, "dchi2 Planck %.2f (half %.2f, q0.01 %.2f) LCDM %.2f  l2-5 %.2f  ratio l2..5 %s" %
                  (dP, dPh, dP01, dL, res["states"][X]["dchi2_planck_l2to5"], np.round(D[:4] / ref[:4], 2)), flush=True)
        out["hosts"].append(res)
    gates = dict(G0=bool(max(g0) <= 5e-3), G1=bool(all(g1)), G2=bool(all(any(s["defined"] for s in H["states"].values()) for H in out["hosts"])))
    cls = []
    for H in out["hosts"]:
        v = [s["dchi2_planck"] for s in H["states"].values() if s["defined"]]
        H["excluded"] = bool(v and all(x >= 25 for x in v)); H["compatible"] = bool(v and all(x < 9 for x in v)); H["preferred"] = bool(v and all(x <= -4 for x in v))
        cls.append(H)
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif all(H["excluded"] for H in cls): verdict = "VISIBLES_EXCLUIDOS"
    elif all(H["preferred"] for H in cls): verdict = "VISIBLES_PREFERIDOS"
    elif all(H["compatible"] for H in cls): verdict = "VISIBLES_COMPATIBLES"
    else: verdict = "MIXTO"
    out.update(gates=gates, G0_max=max(g0), verdict=verdict, seconds=time.time() - t0)
    save("u62e_cl_class.json", out)
    print(json.dumps(dict(gates=gates, G0_max=max(g0), verdict=verdict,
                          per_host=[(H["par"]["AH"], round(H["par"]["f"], 3), H["excluded"], H["compatible"], H["preferred"]) for H in cls]), indent=1))


if __name__ == "__main__":
    main()
