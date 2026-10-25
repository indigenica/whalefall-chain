"""U19b analysis (docs/preregistros/ESPEC_U19b.md, enmienda 1): per network the misalignment number N_mis (homogeneous control of
its release), K_lat, the soft momentum X_lat = sum_{k < 10 m_G a_1} k n_k / (m_G a_1 N_mis) - kappa, Y_lat = X_lat / [(8 pi/27) xi_1 ell_1],
the epoch (r_W, r_N), c_p, c_E, eta_lat and the hard energy fraction; validity checks before any classification.
    python -B network/u19b_analyze.py [runs|pilot]      ->  network/results/u19b/U19B_RESULT.json (or U19B_PILOT.json)"""
import glob, json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; D = HERE / "results" / "u19b"
C827 = 8 * math.pi / 27
Y_STRICT, Y_PRIMARY, Y_B4A = 0.59, 1.14, 0.10                  # ESPEC_U19b: X = 162 Y (B2), X <= 95 / 185; B4a X = 1800 Y


def finite(x):
    if isinstance(x, dict): return {k: finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [finite(v) for v in x]
    if isinstance(x, (float, np.floating)): return float(x) if math.isfinite(x) else None
    if isinstance(x, np.bool_): return bool(x)
    return x


def per_run(run, val):
    R = run["rows"]; m = run["m_G"]; l1 = run["ell_1"]; a1 = run["a_1"]; lg = run["ell_G"]; nm = run["model"]
    meas = [x for x in R if not x.get("initial")]
    em = np.array([x["ell"] for x in meas]); xi = np.array([x["xi"] for x in meas]); wall = np.array([x["wall_frac"] for x in meas])
    a = np.array([x["a"] for x in meas])
    xi1 = float(np.interp(l1, em, xi))
    kin = np.array([(x["ell"], x["N_kin"]) for x in meas if "N_kin" in x]); Nk1 = float(np.interp(l1, kin[:, 0], kin[:, 1]))
    lin = [x for x in meas if "N_lin" in x]; last = meas[-1]; Nend = last["N_lin"]
    if run["protocol"] == "sparse":
        V_early = [x for x in meas if x["ell"] <= lg - 0.3][-1]["V_eps"]; rel = "early"
    else:
        V_early = R[0]["V_eps"]; rel = "dense"
    ctl = val["V2"][f"{nm}_{lg:g}_{rel}"]; N_mis = ctl["R"] * V_early * a1 ** 3 / m
    kappa = val["V5"][f"{nm}_{lg:g}"]["kappa"]
    e = np.array(run["k_edges"]); kc = 0.5 * (e[1:] + e[:-1]); sp = last["lin_spectrum"]; n = np.array(sp["n"]); kn = np.array(sp["kn"])
    soft = kc < 10 * m * a1
    # enmienda 3: subtract the quanta already radiated at ell_1 (kinetic spectrum of measure(), interpolated in ell), as K_lat does
    sk = [x for x in meas if "spectrum" in x]; ks = np.array([x["ell"] for x in sk])
    def kin(x):
        ee = np.array(x["spectrum"]["k_edges"]); kk = 0.5 * (ee[1:] + ee[:-1]); w1 = np.sqrt(kk ** 2 / x["a"] ** 2 + m ** 2)
        nn = x["a"] ** 3 * 2 * np.array(x["spectrum"]["dE"]) / w1; sel = kk < 10 * m * a1
        return float(nn[sel].sum()), float((kk * nn)[sel].sum())
    KN = np.array([kin(x) for x in sk]); N1_soft = float(np.interp(l1, ks, KN[:, 0])); P1_soft = float(np.interp(l1, ks, KN[:, 1]))
    X = float((kn[soft].sum() - P1_soft) / (m * a1 * N_mis)) - kappa
    Y = X / (C827 * xi1 * l1)
    K = (Nend - Nk1) / N_mis - 1
    # epoch
    imax = int(np.argmax(wall)); after = np.flatnonzero((np.arange(len(wall)) > imax) & (wall < 0.1 * wall[imax]))
    r_W = float(a[after[0]] / a1) if len(after) else float("nan")
    aN = np.array([x["a"] for x in lin]); NN = np.array([x["N_lin"] for x in lin])
    out_ = np.flatnonzero(np.abs(NN / Nend - 1) >= 0.1)                   # enmienda 3: from where N stays within 10 % of its final value
    r_N = float(aN[out_[-1] + 1] / a1) if len(out_) and out_[-1] + 1 < len(aN) else (float(aN[0] / a1) if not len(out_) else float("nan"))
    # excess quanta (misalignment taken as N_mis quanta with mean comoving momentum kappa m a_1)
    n_exc = n[soft].sum() - N_mis - N1_soft; k_exc = kn[soft].sum() - P1_soft - kappa * m * a1 * N_mis
    c_p = float(k_exc / n_exc / (m * a1 * r_W)) if n_exc > 0 and math.isfinite(r_W) else float("nan")
    cE_all = np.sum(n[soft] * np.sqrt(1 + (kc[soft] / (m * a1 * r_W)) ** 2)) if math.isfinite(r_W) else float("nan")
    c_E = float((cE_all - N_mis - N1_soft) / n_exc) if n_exc > 0 else float("nan")             # misalignment and pre-ell_1 quanta ~ at rest
    eta = (27 / (8 * math.pi)) * K * c_E / (xi1 * l1)
    aend = last["a"]; w = np.sqrt(kc ** 2 / aend ** 2 + m ** 2); f_hard = float(np.sum((n * w)[~soft]) / np.sum(n * w))
    ell_B = em[-1]; Nhalf = float(np.interp(ell_B - 0.5, [x["ell"] for x in lin], NN))
    checks = dict(strings=bool(xi[-1] < 0.01 * xi1), walls=bool(wall[-1] < 0.05 * wall.max()), counter=bool(last["far_frac"] < 0.01),
                  plateau=bool(abs(Nend / Nhalf - 1) < 0.03))
    spec = dict(k_over_mGa1=(kc / (m * a1)).tolist(), dN_dlnk_over_Nmis=(n / np.diff(np.log(np.maximum(e, 1e-12))) / N_mis).tolist())
    return dict(tag=f"{nm}_lG{lg:g}_{run['protocol']}_s{run['seed']}", model=nm, ell_G=lg, protocol=run["protocol"], seed=run["seed"],
                xi_1=xi1, ell_1=l1, N_mis=N_mis, N_kin_1=Nk1, N_kin_1_soft=N1_soft, P_kin_1_soft=P1_soft, N_lin_end=Nend, K_lat=K, kappa=kappa, X_lat=X, Y_lat=Y, Y_kappa=kappa / (C827 * xi1 * l1),
                r_W=r_W, r_N=r_N, c_p=c_p, c_E=c_E, eta_lat=eta, f_E_hard=f_hard, far_frac_end=last["far_frac"], xi_end=float(xi[-1]),
                wall_end_over_max=float(wall[-1] / max(wall.max(), 1e-30)), plateau_rel=Nend / Nhalf - 1, checks=checks, valid=all(checks.values()),
                spectrum_end=spec)


def block(rows, val, nm, lg):
    rr = [r for r in rows if r["model"] == nm and r["ell_G"] == lg]
    if not rr: return None
    Y = np.array([r["Y_lat"] for r in rr]); n = len(Y)
    stat = float(np.std(Y, ddof=1) / math.sqrt(n)) if n > 1 else float("inf")
    kap = float(np.mean([r["Y_kappa"] for r in rr]))
    # U19's delta is an absolute change (in N_mis units) of the anharmonic misalignment yield between releases; those quanta carry
    # momentum ~ kappa m_G a_1, so on Y it is |delta| * Y_kappa, for the dense share of the block
    dense = [r for r in rr if r["protocol"] == "dense"]
    anh = abs(val["anharmonic"].get(f"{nm}_{lg:g}_delta", 0.0)) * (float(np.mean([r["Y_kappa"] for r in dense])) if dense else 0.0) * len(dense) / n
    anh_bias = (val["anharmonic"][f"{nm}_{lg:g}_early"]["ratio"] - 1) * kap             # informative: large-angle misalignment under-subtracted
    sig = math.sqrt(stat ** 2 + kap ** 2 + anh ** 2)
    return dict(n=n, Y=float(np.mean(Y)), sigma=sig, sigma_parts=dict(stat=stat, kappa=kap, anharmonic=anh), Y_anharmonic_bias_upper=anh_bias,
                valid=all(r["valid"] for r in rr),
                r_W=float(np.nanmean([r["r_W"] for r in rr])), r_N=float(np.nanmean([r["r_N"] for r in rr])),
                c_p=float(np.nanmean([r["c_p"] for r in rr])), c_E=float(np.nanmean([r["c_E"] for r in rr])),
                eta_lat=float(np.nanmean([r["eta_lat"] for r in rr])), f_E_hard=float(np.mean([r["f_E_hard"] for r in rr])))


def verdict(b, val_pass):
    if not val_pass: return "INVALIDO (validación)"
    if b is None or b["n"] < 4: return "INDETERMINADO (redes primarias incompletas)"
    if not b["valid"]: return "INDETERMINADO (condición de validez)"
    if not (b["sigma"] / abs(b["Y"]) <= 0.5): return "INDETERMINADO (error relativo > 50 %)"
    if b["Y"] + b["sigma"] <= Y_STRICT: return "COLAPSO_DURO"
    if b["Y"] - b["sigma"] > Y_PRIMARY: return "COLAPSO_BLANDO"
    return "COLAPSO_INTERMEDIO"


def diagnostic_e4(rows):
    """ESPEC_U19b enmienda 4 (informative, no verdict): does the soft yield X follow the string energy xi_1 ell_1?"""
    def mean(nm, lg, pr, key):
        v = [r[key] if key != "xiell" else r["xi_1"] * r["ell_1"] for r in rows if r["model"] == nm and r["ell_G"] == lg and r["protocol"] == pr]
        return (float(np.mean(v)), len(v)) if v else (None, 0)
    out = dict(D1={}, D2={})
    for lg in (3.2, 2.4):
        (xs, ns), (xd, nd) = mean("RP3", lg, "sparse", "X_lat"), mean("RP3", lg, "dense", "X_lat")
        (es, _), (ed, _) = mean("RP3", lg, "sparse", "xiell"), mean("RP3", lg, "dense", "xiell")
        if ns and nd: out["D1"][f"{lg:g}"] = dict(X_dense_over_sparse=xd / xs, xiell_dense_over_sparse=ed / es, n=(ns, nd))
    full = True
    for pr in ("sparse", "dense"):
        (x2, n2), (x3, n3) = mean("RP3", 2.4, pr, "X_lat"), mean("RP3", 3.2, pr, "X_lat")
        (e2, _), (e3, _) = mean("RP3", 2.4, pr, "xiell"), mean("RP3", 3.2, pr, "xiell")
        if n2 and n3: out["D2"][pr] = dict(rho_X=x2 / x3, rho_xiell=e2 / e3, n=(n2, n3))
        full = full and n2 == 2 and n3 == 2
    # enmienda 5 (informative): affine split per l_G, X = X0 + y_s Z with Z = (8 pi/27) xi_1 ell_1, least squares over the runs
    out["affine"] = {}
    for lg in (3.2, 2.4):
        rr = [r for r in rows if r["model"] == "RP3" and r["ell_G"] == lg]
        if len({r["protocol"] for r in rr}) < 2: continue
        Z = np.array([C827 * r["xi_1"] * r["ell_1"] for r in rr]); X = np.array([r["X_lat"] for r in rr])
        A = np.column_stack([np.ones_like(Z), Z]); coef, *_ = np.linalg.lstsq(A, X, rcond=None)
        dof = len(X) - 2; res = X - A @ coef
        cov = (float(res @ res) / dof) * np.linalg.inv(A.T @ A) if dof > 0 else np.full((2, 2), np.nan)
        full_lg = all(sum(1 for r in rr if r["protocol"] == p) >= 2 for p in ("sparse", "dense"))
        out["affine"][f"{lg:g}"] = dict(X0=float(coef[0]), X0_err=float(np.sqrt(cov[0, 0])), y_s=float(coef[1]), y_s_err=float(np.sqrt(cov[1, 1])),
                                        n=len(rr), complete=bool(full_lg), X_universe_B2=float(coef[0] + 162 * coef[1]),
                                        X_universe_B2_upper=float(coef[0] + 162 * max(coef[1] + 2 * np.sqrt(cov[1, 1]), 0.0)) if dof > 0 else None)
    af = out["affine"]
    if "3.2" in af and "2.4" in af and af["3.2"]["complete"] and af["2.4"]["complete"]:   # enmienda 5: three X0 extrapolations to l ~ 90
        a3, a2 = af["3.2"], af["2.4"]; l3, l2, lu = 3.2 + math.log(3), 2.4 + math.log(3), 90.0
        ys = 162 * max(max(a3["y_s"] + 2 * a3["y_s_err"], a2["y_s"] + 2 * a2["y_s_err"]), 0.0)
        ext = {"X0_constante": a3["X0"], "X0_lineal_en_l1": a3["X0"] + (a3["X0"] - a2["X0"]) / (l3 - l2) * (lu - l3),
               "X0_proporcional_a_l": a3["X0"] * lu / l3}
        out["X_universe"] = {k: dict(X=v + ys, strict_dA=bool(v + ys <= 95), primary=bool(v + ys <= 185)) for k, v in ext.items()}
    if full:
        d = out["D2"].values()
        a = all(abs(v["rho_X"] / v["rho_xiell"] - 1) <= 0.15 for v in d)
        b = all(abs(v["rho_X"] - 1) <= 0.15 and abs(v["rho_xiell"] - 1) > 0.3 for v in d)
        out["reading"] = "(c) mixta" if a == b else ("(a) escalado de U17 (Y universal)" if a else "(b) rendimiento blando independiente de las cuerdas")
    else:
        out["reading"] = "pendiente (faltan redes de l_G = 2.4)"
    return out


def main(sub="runs"):
    val = json.loads((D / "U19B_VALIDATION.json").read_text())
    rows = [per_run(json.loads(Path(p).read_text()), val) for p in sorted(glob.glob(str(D / sub / "*.json")))]
    out = dict(spec="docs/preregistros/ESPEC_U19b.md", validation_pass=val["PASS"], runs=rows,
               blocks={f"{nm}_{lg:g}": block(rows, val, nm, lg) for nm in ("RP3", "U1") for lg in (3.2, 2.4)})
    b = out["blocks"]["RP3_3.2"]
    if sub == "runs":
        out["verdict_U19b"] = verdict(b, val["PASS"])
        if b and b["valid"]:
            out["B4a_viable"] = bool(b["Y"] + b["sigma"] <= Y_B4A)
            b2 = out["blocks"]["RP3_2.4"]
            if b2 and b2["Y"] > 0 and b["Y"] > 0:
                g1 = math.log(b["Y"] / b2["Y"]) / math.log((3.2 + math.log(3)) / (2.4 + math.log(3)))
                sg = math.sqrt((b["sigma"] / b["Y"]) ** 2 + (b2["sigma"] / b2["Y"]) ** 2) / math.log((3.2 + math.log(3)) / (2.4 + math.log(3)))
                out["gamma_minus_1"] = dict(value=g1, sigma=sg)
        out["diagnostic_enmienda4"] = diagnostic_e4(rows)
    name = "U19B_RESULT.json" if sub == "runs" else f"U19B_{sub.upper()}.json"
    (D / name).write_text(json.dumps(finite(out), indent=1, allow_nan=False))
    for r in rows:
        print(f"{r['tag']:<28} xi1={r['xi_1']:.3f} K={r['K_lat']:.3f} X={r['X_lat']:.3f} Y={r['Y_lat']:.3f} r_W={r['r_W']:.2f} r_N={r['r_N']:.2f} "
              f"c_p={r['c_p']:.2f} eta={r['eta_lat']:.3f} f_hard={r['f_E_hard']:.3f} far={r['far_frac_end']:.4f} valid={r['valid']} {r['checks']}")
    if "diagnostic_enmienda4" in out: print("enmienda 4:", out["diagnostic_enmienda4"])
    print("U19b:", out.get("verdict_U19b", "(pilot: no verdict)"), {k: (v and {kk: v[kk] for kk in ("Y", "sigma", "n", "valid")}) for k, v in out["blocks"].items()})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs")
