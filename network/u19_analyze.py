"""U19 analysis (docs/preregistros/ESPEC_U19.md, amendments 1-2): K_lat = R / R_V2 - 1 per network, fit K_lat = A + B xi_1 ell_1,
c/eta = (8 pi / 27) / B.  Dense networks carry the anharmonic systematic delta (amendment 2) in quadrature.
Output: network/results/u19/U19_RESULT.json"""
import glob, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; D = HERE / "results" / "u19"
MR = math.sqrt(0.24)


def finite(x):
    if isinstance(x, dict): return {k: finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [finite(v) for v in x]
    if isinstance(x, (float, np.floating)): return float(x) if math.isfinite(x) else None
    if isinstance(x, np.bool_): return bool(x)
    return x


def per_run(run, val):
    R = run["rows"]; m = run["m_G"]; s = run["s"]; l1 = run["ell_1"]
    ell = np.array([x["ell"] for x in R])
    meas = [x for x in R if not x.get("initial")]
    em = np.array([x["ell"] for x in meas]); xi = np.array([x["xi"] for x in meas])
    xi1 = float(np.interp(l1, em, xi))
    kin = [(x["ell"], x["N_kin"]) for x in meas if "N_kin" in x]
    k_e, k_n = np.array(kin).T
    Nk1 = float(np.interp(l1, k_e, k_n))
    last = meas[-1]; Nlin = last["N_lin"]
    if run["protocol"] == "sparse":
        early = [x for x in meas if x["ell"] <= 2.5]; Vearly = early[-1]["V_eps"]
    else:
        Vearly = R[0]["V_eps"]
    tau1 = math.sqrt(math.exp(l1) / (MR * s)); a1 = s * tau1
    Rv = m * (Nlin - Nk1) / (Vearly * a1 ** 3)
    rel = "early" if run["protocol"] == "sparse" else "dense"
    Rv2 = val["V2"][f"{run['model']}_{rel}"]["R"]
    wall = np.array([x["wall_frac"] for x in meas])
    return dict(tag=f"{run['model']}_{run['protocol']}_s{run['seed']}", model=run["model"], protocol=run["protocol"], seed=run["seed"],
                xi_1=xi1, ell_1=l1, x=xi1 * l1, N_kin_1=Nk1, N_lin_end=Nlin, V_early=Vearly, R=Rv, R_V2=Rv2, K_lat=Rv / Rv2 - 1,
                wall_end_over_max=float(wall[-1] / max(wall.max(), 1e-30)), far_frac_end=last["far_frac"], ell_end=float(em[-1]),
                N_kin_end=last.get("N_kin"))


def fit(rows, delta):
    x = np.array([r["x"] for r in rows]); K = np.array([r["K_lat"] for r in rows])
    dense = np.array([r["protocol"] == "dense" for r in rows])
    sig = []
    for prot in ("sparse", "dense"):
        k = K[[r["protocol"] == prot for r in rows]]
        sig.append(float(np.std(k, ddof=1)) if len(k) > 1 else float("nan"))
    s_seed = np.where(dense, sig[1], sig[0]); s_seed = np.where(np.isfinite(s_seed), s_seed, np.nanmax(sig))
    s_tot = np.sqrt(s_seed ** 2 + np.where(dense, delta ** 2, 0.0))
    w = 1 / np.maximum(s_tot, 1e-6) ** 2
    Xm = np.vstack([np.ones_like(x), x]).T; C = np.linalg.inv(Xm.T @ (w[:, None] * Xm)); A, B = C @ (Xm.T @ (w * K))
    sB = math.sqrt(C[1, 1]); ce = 8 * math.pi / 27 / B if B != 0 else float("inf")
    sce = abs(ce) * sB / abs(B) if B != 0 else float("inf")
    return dict(A=float(A), B=float(B), sigma_B=sB, c_over_eta=ce, sigma_c_over_eta=sce, sigma_seed=sig, delta_anh=delta,
                x_mean_sparse=float(x[~dense].mean()), x_mean_dense=float(x[dense].mean()))


def verdict(f, rows, validation_pass=True):
    # vetoes before any classification of the fit (codex audit, 14 Oct 2026): validation, linear-counter validity, complete runs
    if not validation_pass: return "INDETERMINADO (validación fallida)"
    if not all(np.isfinite(r.get("far_frac_end", np.nan)) and r["far_frac_end"] < 0.01 for r in rows): return "INDETERMINADO (contador no válido)"
    if not all(r.get("N_kin_end") is not None and np.isfinite(r["K_lat"]) for r in rows): return "INDETERMINADO (ejecución incompleta)"
    if not all(r["wall_end_over_max"] <= 0.05 for r in rows): return "INDETERMINADO (colapso sin terminar)"
    if f["B"] < 2 * f["sigma_B"]: return "INDETERMINADO (B < 2 sigma)"
    if f["sigma_c_over_eta"] / f["c_over_eta"] > 0.5: return "INDETERMINADO (error relativo > 50 %)"
    ce, s = f["c_over_eta"], f["sigma_c_over_eta"]
    if ce + s <= 6: return "C_ETA_BAJO"
    if ce - s >= 20: return "C_ETA_ALTO"
    return "C_ETA_MEDIO"


def main():
    val = json.loads((D / "U19_VALIDATION.json").read_text())
    runs = [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(D / "runs" / "*.json")))]
    out = dict(spec="docs/preregistros/ESPEC_U19.md", validation_pass=val["PASS"], models={})
    for model in ("RP3", "U1"):
        rows = [per_run(r, val) for r in runs if r["model"] == model]
        if len({r["protocol"] for r in rows}) < 2: out["models"][model] = dict(rows=rows, note="incomplete"); continue
        f = fit(rows, abs(val["anharmonic"][f"{model}_delta"]))
        out["models"][model] = dict(rows=rows, fit=f, verdict=verdict(f, rows, bool(val["PASS"])),
                                    far_frac_end_ok=bool(all(r["far_frac_end"] < 0.01 for r in rows)))
    if "fit" in out["models"].get("RP3", {}): out["verdict_U19"] = out["models"]["RP3"]["verdict"]
    else:   # amendment 3: one RP3 confirmation run; stop if the linear counter is outside its validity
        rp = out["models"].get("RP3", {}).get("rows", [])
        if rp and any(r["far_frac_end"] > 0.01 for r in rp): out["verdict_U19"] = "INDETERMINADO (contador no válido)"
    (D / "U19_RESULT.json").write_text(json.dumps(finite(out), indent=1, allow_nan=False))
    for mdl, v in out["models"].items():
        for r in v["rows"]: print(mdl, r["tag"], f"x={r['x']:.2f} K_lat={r['K_lat']:.3f} wall_end={r['wall_end_over_max']:.3f} far={r['far_frac_end']:.4f}")
        if "fit" in v: print(mdl, {k: v["fit"][k] for k in ("A", "B", "sigma_B", "c_over_eta", "sigma_c_over_eta")}, v["verdict"])


if __name__ == "__main__":
    main()
