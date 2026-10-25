"""U19d (docs/preregistros/ESPEC_U19d.md): the U19c analysis (u19c_analyze) on the U19c networks plus the U19d seeds 93-96.
Only seeds complete at all four masses (both protocols) enter, besides the U19c seeds 91 and 92.  Also reported, outside the
verdict: the verdict with the new seeds only, the power law of the mean X_lat per mass (chosen after seeing U19c), and p
leaving each mass out.

    python -B network/u19d_analyze.py            ->  network/results/u19d/U19D_RESULT.json"""
import glob, json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u19b_analyze as A                                                      # noqa: E402
import u19c_analyze as C                                                      # noqa: E402
import u19d_produce as D                                                      # noqa: E402

HERE = Path(__file__).resolve().parent; OUT = HERE / "results" / "u19d"; UC = HERE / "results" / "u19c"; UB = HERE / "results" / "u19b"
L_UNIV, HIGH = 90.0, 576


def runs(N, seeds):
    """per-run rows (u19b_analyze.per_run) of mass N for the given seeds, from U19c/U19b and U19d"""
    lg = D.LG[N]; rows = []
    if N == 512: val = json.loads((UB / "U19B_VALIDATION.json").read_text()); src = {91: UB, 92: UB}
    else: val = C.val_u19b_format(json.loads((UC / "U19C_VALIDATION.json").read_text())["masses"][f"{lg:g}"], lg); src = {91: UC, 92: UC}
    for s in seeds:
        for prot in ("sparse", "dense"):
            f = (src.get(s, OUT) / "runs" / f"RP3_N{N}_lG{lg:g}_{prot}_s{s}.json")
            rows.append(A.per_run(json.loads(f.read_text()), val))
    return rows


def fit(per_mass):
    """the U19c power law, models and verdict from {N: rows}"""
    out = dict(masses={})
    for N, rows in sorted(per_mass.items()):
        a = C.affine([dict(r, ell_1=D.LG[N] + math.log(3.0)) for r in rows])
        out["masses"][str(N)] = dict(ell_1=D.LG[N] + math.log(3.0), n_runs=len(rows), valid=all(r["valid"] for r in rows), affine=a,
                                     X_lat_mean=float(np.mean([r["X_lat"] for r in rows])), X_lat_sd=float(np.std([r["X_lat"] for r in rows], ddof=1)))
    good = {int(k): v for k, v in out["masses"].items() if v["valid"]}
    l1 = np.array([good[N]["ell_1"] for N in sorted(good)]); X0 = np.array([good[N]["affine"]["X0"] for N in sorted(good)])
    sX = np.array([good[N]["affine"]["X0_err"] for N in sorted(good)])
    ys = 162 * max(max(good[N]["affine"]["y_s"] + 2 * good[N]["affine"]["y_s_err"] for N in good), 0.0)
    w = 1 / (sX / X0) ** 2; lp = float(np.exp(np.sum(w * np.log(l1)) / np.sum(w)))
    beta, cov, chi2, nu = C.wls(np.log(l1 / lp), np.log(X0), sX / X0, lambda x: np.column_stack([np.ones_like(x), x]))
    p, sp = float(beta[1]), float(math.sqrt(cov[1, 1]) * (math.sqrt(chi2 / nu) if nu > 0 and chi2 / nu > 1 else 1.0))
    Xu = lambda pp: float(math.exp(beta[0] + pp * math.log(L_UNIV / lp)) + ys)
    out["power_law"] = dict(p=p, sigma_p=sp, pivot_ell_1=lp, chi2=chi2, nu=nu, strings_term=ys, X_univ=Xu(p), X_lo=Xu(p - sp), X_hi=Xu(p + sp))
    models = {}
    b, c, x2, n2 = C.wls(l1, X0, sX, lambda x: np.ones((len(x), 1))); models["X0_constante"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] + ys))
    b, c, x2, n2 = C.wls(l1, X0, sX, lambda x: x[:, None]); models["X0_proporcional_a_l1"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] * L_UNIV + ys))
    b, c, x2, n2 = C.wls(l1, X0, sX, lambda x: np.column_stack([np.ones_like(x), x])); models["X0_lineal_en_l1"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] + b[1] * L_UNIV + ys))
    models["potencia_libre"] = dict(chi2=chi2, nu=nu, X_univ=Xu(p)); out["models"] = models
    if len(good) < 3 or HIGH not in good or not sp <= 0.5:
        v = "INDETERMINADO (" + ("menos de 3 masas válidas" if len(good) < 3 else "falta la masa alta" if HIGH not in good else "sigma_p > 0.5") + ")"
    else:
        lo, hi = Xu(p - sp), Xu(p + sp)
        v = ("X_UNIV_BAJO" if hi <= C.X_F0 else "X_UNIV_ALTO" if lo > C.X_085 else "X_UNIV_MEDIO" if (lo > C.X_F0 and hi <= C.X_085)
             else f"X_UNIV_INDECISO (entre {lo:.0f} y {hi:.0f})")
    out["verdict"] = v
    return out


def mean_xlat_power(res):
    """informative (chosen after seeing U19c): power law of the mean X_lat per mass, errors sd/sqrt(n)"""
    m = res["masses"]; Ns = sorted(m, key=int)
    l1 = np.array([m[N]["ell_1"] for N in Ns]); x = np.array([m[N]["X_lat_mean"] for N in Ns])
    se = np.array([m[N]["X_lat_sd"] / math.sqrt(m[N]["n_runs"]) for N in Ns]); lp = 4.0
    b, cov, chi2, nu = C.wls(np.log(l1 / lp), np.log(x), se / x, lambda z: np.column_stack([np.ones_like(z), z]))
    sp = math.sqrt(cov[1, 1]) * (math.sqrt(chi2 / nu) if nu > 0 and chi2 / nu > 1 else 1.0)
    return dict(p=float(b[1]), sigma_p=float(sp), chi2=float(chi2), nu=int(nu), X_lat_at_90=float(math.exp(b[0] + b[1] * math.log(L_UNIV / lp))))


def main():
    new = [s for s in D.SEEDS if all((OUT / "runs" / f"RP3_N{N}_lG{D.LG[N]:g}_{prot}_s{s}.json").exists() for N in D.LG for prot in ("sparse", "dense"))]
    seeds = [91, 92] + new
    res = fit({N: runs(N, seeds) for N in D.LG}); res.update(spec="docs/preregistros/ESPEC_U19d.md", seeds=seeds, new_seeds_complete=new)
    out = dict(res, verdict_U19d=res.pop("verdict"))
    if len(new) >= 2:
        o = fit({N: runs(N, new) for N in D.LG}); out["only_new_seeds"] = dict(verdict=o["verdict"], power_law=o["power_law"])
    out["mean_X_lat_power_informative"] = mean_xlat_power(res)
    out["leave_one_mass_out"] = {str(N): (lambda o: dict(p=o["power_law"]["p"], sigma_p=o["power_law"]["sigma_p"]))(
        fit({M: runs(M, seeds) for M in D.LG if M != N})) for N in D.LG}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "U19D_RESULT.json").write_text(json.dumps(A.finite(out), indent=1, default=float, allow_nan=False))
    print("U19d:", out["verdict_U19d"], "seeds", seeds, json.dumps(out["power_law"], default=float), flush=True)


if __name__ == "__main__":
    main()
