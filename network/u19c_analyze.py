"""U19c (docs/preregistros/ESPEC_U19c.md): X0(l_1) in four masses with the same relative protocol, its power law and the
X of the universe against the U17f thresholds.  Per-run quantities are u19b_analyze.per_run (definitions of U19b, enmiendas
1-3); per mass the affine split X = X0 + y_s Z of U19b enmienda 5 (Z = (8 pi/27) xi_1 l_1, 2 protocols x 2 seeds).
Masses: the new U19c networks (network/results/u19c/runs) and the four primary U19b networks with l_G = 3.2 at 512^3.
Fit: ln X0 = ln X0(l_p) + p (ln l_1 - ln l_p), weighted with sigma(X0)/X0, pivot l_p at the weighted mean of ln l_1 (so that
p - sigma_p and p + sigma_p move the extrapolation only); sigma_p scaled by sqrt(chi2/nu) if chi2/nu > 1.
X_univ = X0(90) + 162 max(y_s + 2 sigma, 0) with the least favourable y_s of the masses (U19b enmienda 5).

    python -B network/u19c_analyze.py            ->  network/results/u19c/U19C_RESULT.json"""
import glob, json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u19b_analyze as A                                                      # noqa: E402
import u19c_produce as P                                                      # noqa: E402

HERE = Path(__file__).resolve().parent; OUT = P.OUT; L_UNIV = 90.0
X_F0, X_085, X_07 = 194.0, 275.0, 320.0                                       # U17f: F0 with X <= 194, 0.85 F0 with X <= 275, 0.7 F0 at ~320
HIGH = (768, 640, 576)                                                        # the high mass (enmienda 1: 640^3, or 576^3; enmienda 2: 768^3 if valid)


def val_u19b_format(v, lg):
    """the U19c validation of one mass in the key format that u19b_analyze.per_run/block read"""
    k = f"{lg:g}"
    return dict(V2={f"{nm}_{k}_{rel}": v["V2"][f"{nm}_{rel}"] for nm in ("RP3", "U1") for rel in ("early", "dense")},
                V5={f"{nm}_{k}": v["V5"][nm] for nm in ("RP3", "U1")},
                anharmonic={f"{nm}_{k}_{x}": v["anharmonic"][f"{nm}_{x}"] for nm in ("RP3", "U1") for x in ("early", "dense", "delta")})


def affine(rr):
    Z = np.array([A.C827 * r["xi_1"] * r["ell_1"] for r in rr]); X = np.array([r["X_lat"] for r in rr])
    M = np.column_stack([np.ones_like(Z), Z]); coef, *_ = np.linalg.lstsq(M, X, rcond=None)
    dof = len(X) - 2; res = X - M @ coef
    cov = (float(res @ res) / dof) * np.linalg.inv(M.T @ M) if dof > 0 else np.full((2, 2), np.nan)
    return dict(X0=float(coef[0]), X0_err=float(np.sqrt(cov[0, 0])), y_s=float(coef[1]), y_s_err=float(np.sqrt(cov[1, 1])), n=len(X))


def wls(x, y, s, design):
    """weighted least squares y ~ design(x) @ beta with errors s: beta, cov, chi2, nu"""
    M = design(x) / s[:, None]; b = y / s
    beta, *_ = np.linalg.lstsq(M, b, rcond=None); cov = np.linalg.inv(M.T @ M); chi2 = float(np.sum((b - M @ beta) ** 2))
    return beta, cov, chi2, len(x) - M.shape[1]


def main():
    out = dict(spec="docs/preregistros/ESPEC_U19c.md", masses={})
    vpath = OUT / "U19C_VALIDATION.json"; vals = json.loads(vpath.read_text())["masses"] if vpath.exists() else {}
    vb = json.loads((A.D / "U19B_VALIDATION.json").read_text())
    sets = {512: (3.2, sorted(glob.glob(str(A.D / "runs" / "RP3_N512_lG3.2_*.json"))), vb, bool(vb.get("PASS")))}
    for N, lg in P.MASSES.items():
        v = vals.get(f"{lg:g}")
        sets[N] = (lg, sorted(glob.glob(str(OUT / "runs" / f"RP3_N{N}_lG{lg:g}_*.json"))), val_u19b_format(v, lg) if v else None, bool(v and v.get("PASS")))
    for N, (lg, files, val, vpass) in sorted(sets.items()):
        m = dict(N=N, ell_G=lg, ell_1=lg + math.log(3.0), validation_pass=vpass, n_runs=len(files), reused_from_U19b=(N == 512))
        if files and val:
            rows = [A.per_run(json.loads(Path(f).read_text()), val) for f in files]
            m["runs"] = [{k: r[k] for k in ("tag", "protocol", "seed", "xi_1", "X_lat", "Y_lat", "r_W", "r_N", "c_p", "valid", "checks")} for r in rows]
            full = all(sum(1 for r in rows if r["protocol"] == pr) >= 2 for pr in ("sparse", "dense"))
            m.update(complete=full, valid=bool(vpass and full and all(r["valid"] for r in rows)))
            if full: m["affine"] = affine(rows)
            b = A.block(rows, val, "RP3", lg)
            if b: m["Y_block"] = {k: b[k] for k in ("Y", "sigma", "r_W", "r_N", "c_p")}
        else:
            m.update(complete=False, valid=False)
        out["masses"][str(N)] = m
    good = {int(k): v for k, v in out["masses"].items() if v.get("valid") and "affine" in v}
    high = next((N for N in HIGH if N in good), None)
    out["masses_in_fit"] = sorted(good); out["high_mass"] = high
    if len(good) >= 2:
        l1 = np.array([good[N]["ell_1"] for N in sorted(good)]); X0 = np.array([good[N]["affine"]["X0"] for N in sorted(good)])
        sX = np.array([good[N]["affine"]["X0_err"] for N in sorted(good)])
        ys = 162 * max(max(good[N]["affine"]["y_s"] + 2 * good[N]["affine"]["y_s_err"] for N in good), 0.0)
        w = 1 / (sX / X0) ** 2; lp = float(np.exp(np.sum(w * np.log(l1)) / np.sum(w)))
        beta, cov, chi2, nu = wls(np.log(l1 / lp), np.log(X0), sX / X0, lambda x: np.column_stack([np.ones_like(x), x]))
        p, sp = float(beta[1]), float(math.sqrt(cov[1, 1]) * (math.sqrt(chi2 / nu) if nu > 0 and chi2 / nu > 1 else 1.0))
        Xu = lambda pp: float(math.exp(beta[0] + pp * math.log(L_UNIV / lp)) + ys)
        out["power_law"] = dict(p=p, sigma_p=sp, pivot_ell_1=lp, X0_at_pivot=float(math.exp(beta[0])), chi2=chi2, nu=nu, strings_term=ys,
                                X_univ=Xu(p), X_lo=Xu(p - sp), X_hi=Xu(p + sp))
        models = {}
        b, c, x2, n2 = wls(l1, X0, sX, lambda x: np.ones((len(x), 1))); models["X0_constante"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] + ys))
        b, c, x2, n2 = wls(l1, X0, sX, lambda x: x[:, None]); models["X0_proporcional_a_l1"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] * L_UNIV + ys))
        if len(good) >= 3:
            b, c, x2, n2 = wls(l1, X0, sX, lambda x: np.column_stack([np.ones_like(x), x])); models["X0_lineal_en_l1"] = dict(chi2=x2, nu=n2, X_univ=float(b[0] + b[1] * L_UNIV + ys))
        models["potencia_libre"] = dict(chi2=chi2, nu=nu, X_univ=Xu(p))
        out["models"] = models
    pl = out.get("power_law")
    if len(good) < 3 or high is None or pl is None or not (pl["sigma_p"] <= 0.5):
        v = "INDETERMINADO (" + ("menos de 3 masas válidas" if len(good) < 3 else "falta la masa alta" if high is None else "sigma_p > 0.5") + ")"
    else:
        lo, hi = pl["X_lo"], pl["X_hi"]
        v = ("X_UNIV_BAJO" if hi <= X_F0 else "X_UNIV_ALTO" if lo > X_085 else "X_UNIV_MEDIO" if (lo > X_F0 and hi <= X_085)
             else f"X_UNIV_INDECISO (entre {lo:.0f} y {hi:.0f})")
    out["verdict_U19c"] = v
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "U19C_RESULT.json").write_text(json.dumps(A.finite(out), indent=1, default=float, allow_nan=False))
    print("U19c:", v, json.dumps(pl, default=float) if pl else "", flush=True)


if __name__ == "__main__":
    main()
