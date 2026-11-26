"""U67c analysis (docs/preregistros/ESPEC_U67c.md): deficit D = 1 - <xi_A>/<xi_B> of a very sparse RP3 network (xi_i = 0.05)
against the RD_SPARSE reference (xi_i = 0.5) at 256^3, its closing rate (exp and power-law fits over 2.5 <= ell <= ell_end,
jackknife over seed pairs) and the declared extrapolation to ell = 80 from a birth at ell_b = 3.4 with xi = 0.05.
    python3 -B network/u67c_analyze.py  ->  network/results/u67c/U67C_RESULT.json"""
import glob, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE / "results" / "u67c" / "runs"
SEEDS = (601, 602, 603, 604)
ELL_LO, ELL_B, XI_B, ELL_C, ELL_G2 = 2.5, 3.4, 0.05, 80.0, 4.56


def load(row):
    out = []
    for sd in SEEDS:
        f = RUNS / f"RP3_N256_{row}_s{sd}.json"
        out.append(json.loads(f.read_text()) if f.exists() else None)
    return out


def fits(ell, xa, xb, drop_last=False):
    """xa, xb: arrays (nseed, nsnap).  Returns slope_exp, tau, k on the mean over the given seeds."""
    R = xa.mean(0) / xb.mean(0)
    D = 1.0 - R
    m = (ell >= ELL_LO) & (D > 0)
    if drop_last:
        m[np.where(m)[0][-1]] = False
    se, _ = np.polyfit(ell[m], np.log(D[m]), 1)
    sp, _ = np.polyfit(np.log(ell[m]), np.log(D[m]), 1)
    return dict(slope_exp=float(se), tau_ell=float(-1.0 / se) if se < 0 else float("inf"), k=float(-sp), n_points=int(m.sum()))


def jack(ell, xa, xb, key, **kw):
    n = xa.shape[0]
    th = np.array([fits(ell, np.delete(xa, i, 0), np.delete(xb, i, 0), **kw)[key] for i in range(n)])
    return float(math.sqrt((n - 1) / n * np.sum((th - th.mean()) ** 2)))


def u67_reference():
    vals = []
    for f in sorted(glob.glob(str(HERE / "results" / "u67" / "runs" / "RP3_N192_RD_SPARSE_s*.json"))):
        rows = json.loads(Path(f).read_text())["rows"]
        e = np.array([r["ell"] for r in rows]); x = np.array([r["xi"] for r in rows])
        vals.append(float(x[np.argmin(abs(e - ELL_G2))]))
    return float(np.mean(vals)), vals


def main():
    A, B = load("A"), load("B")
    res = dict(espec="docs/preregistros/ESPEC_U67c.md", n_done=dict(A=sum(r is not None for r in A), B=sum(r is not None for r in B)))
    if None in A or None in B:
        res["verdict"] = "PENDIENTE"; print(json.dumps(res)); return res
    ell = np.array([r["ell"] for r in B[0]["rows"]])
    for r in A + B:
        assert np.allclose([x["ell"] for x in r["rows"]], ell, atol=1e-9), "snapshot grids differ"
    xa = np.array([[x["xi"] for x in r["rows"]] for r in A])
    xb = np.array([[x["xi"] for x in r["rows"]] for r in B])
    ell_end = float(ell[-1])

    # gates
    finite = bool(np.all(np.isfinite(xa)) and np.all(np.isfinite(xb)))
    xib_end = float(xb[:, -1].mean())
    g1 = finite and 0.5 <= xib_end <= 1.2
    ref, ref_vals = u67_reference()
    i_g2 = int(np.argmin(abs(ell - ELL_G2)))
    xib_g2 = float(xb[:, i_g2].mean())
    g2 = abs(xib_g2 / ref - 1.0) < 0.15

    f = fits(ell, xa, xb)
    s_slope = jack(ell, xa, xb, "slope_exp")
    s_tau = jack(ell, xa, xb, "tau_ell")
    s_k = jack(ell, xa, xb, "k")
    xib_b = float(np.interp(ELL_B, ell, xb.mean(0)))
    D_b = 1.0 - XI_B / xib_b
    tau_p, k_p = f["tau_ell"] + s_tau, f["k"] - s_k
    D80_exp = D_b * math.exp(-(ELL_C - ELL_B) / tau_p) if tau_p > 0 else float("nan")
    D80_pot = D_b * (ELL_B / ELL_C) ** k_p
    D80_exp_c = D_b * math.exp(-(ELL_C - ELL_B) / f["tau_ell"]) if f["tau_ell"] > 0 else float("nan")
    D80_pot_c = D_b * (ELL_B / ELL_C) ** f["k"]

    if not (g1 and g2):
        verdict = "INCONCLUSO"
    elif f["slope_exp"] + 2 * s_slope < 0 and D80_exp < 0.10 and D80_pot < 0.10:
        verdict = "OLVIDA"
    elif f["slope_exp"] + s_slope >= 0 or D80_exp >= 0.30 or D80_pot >= 0.30:
        verdict = "NO OLVIDA"
    else:
        verdict = "INTERMEDIA"

    R = xa.mean(0) / xb.mean(0)
    fd = fits(ell, xa, xb, drop_last=True)
    # 192^3 exploratory comparison (explora_disperso, seeds 101-102, against U67 RD_SPARSE)
    cmp192 = None
    ex = sorted(glob.glob(str(HERE / "results" / "explora_disperso" / "RP3_N192_XI005_s*.json")))
    if ex:
        ea = np.array([[r["xi"] for r in json.loads(Path(p).read_text())["rows"]] for p in ex])
        e192 = np.array([r["ell"] for r in json.loads(Path(ex[0]).read_text())["rows"]])
        eb = np.array([[r["xi"] for r in json.loads(Path(p).read_text())["rows"]] for p in
                       sorted(glob.glob(str(HERE / "results" / "u67" / "runs" / "RP3_N192_RD_SPARSE_s*.json")))])
        n = min(ea.shape[1], eb.shape[1])
        cmp192 = dict(ell=e192[:n].tolist(), R=(ea[:, :n].mean(0) / eb[:, :n].mean(0)).tolist(), fit=fits(e192[:n], ea[:, :n], eb[:, :n]))

    res.update(dict(
        ell_end=ell_end,
        gates=dict(G1=dict(passed=g1, finite=finite, xi_B_end=xib_end),
                   G2=dict(passed=g2, ell=float(ell[i_g2]), xi_B=xib_g2, u67_192_ref=ref, u67_192_seeds=ref_vals)),
        fit=dict(**f, sigma_slope_exp=s_slope, sigma_tau_ell=s_tau, sigma_k=s_k),
        extrapolation=dict(ell_b=ELL_B, xi_B_at_ell_b=xib_b, D_b=D_b, tau_pess=tau_p, k_pess=k_p,
                           D80_exp_pess=D80_exp, D80_pot_pess=D80_pot, D80_exp_central=D80_exp_c, D80_pot_central=D80_pot_c),
        verdict=verdict,
        informative=dict(table=[dict(ell=float(e), xi_A=float(a), xi_B=float(b), R=float(r)) for e, a, b, r in
                                zip(ell, xa.mean(0), xb.mean(0), R)],
                         fit_without_last=fd, xi_A_end_seeds=xa[:, -1].tolist(), xi_B_end_seeds=xb[:, -1].tolist(),
                         cmp_192=cmp192)))
    out = HERE / "results" / "u67c" / "U67C_RESULT.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("gates", "fit", "extrapolation", "verdict")}, indent=1, default=float))
    return res


if __name__ == "__main__":
    main()
