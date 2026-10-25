"""U7 (docs/preregistros/ESPEC_U7.md): calibration of the S4/S6 emission against arXiv:2512.13653 (Correia, Hindmarsh et al.):
nu_ax = n_ax/(f_a^2 H) = 1.66 +- 0.17 for the scaling U(1) network; zeta_* ~ 1.49.

    python -B lya/u7_class.py --nu-only           # U7a only (no CLASS, seconds)
    python -B lya/u7_class.py --workers 16        # U7a + U7b (24 CLASS rows, ~1-2 h)
    python -B lya/u7_class.py --workers 24 --u7c  # U7c (ESPEC_U7c): same rows with R_Gamma in {1.39, 1.69, 2.09}

Output: lya/results/u7/U7_RESULT.json"""
import itertools, json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
OUT7 = HERE / "results" / "u7"
NU_AX, NU_AX_ERR, ZETA = 1.66, 0.17, 1.49
H_EVALS = (1e-30, 1e-33, 1e-36)         # GeV


def nu_emission(F, Trh, eps_fn, x0, H_eval, step=0.03, bins=64):
    """number density of massless Goldstones emitted (log-flat spectrum on [x0 H, m_r]) up to H = H_eval, over f_a^2 H.
    Same emission as chain/s6_measured.source_eps (rho_s = xi mu / t^2 with xi = 0.24 ell, times eps), with E = p."""
    from s4_dark import hubble, scale, gs
    from scipy.optimize import brentq
    fa = F / 2; mr = F / 57.735
    Tev = math.exp(brentq(lambda l: math.log(hubble(math.exp(l)) / H_eval), math.log(1e-12), math.log(Trh)))
    edge = np.linspace(math.log(Trh), math.log(Tev), int(math.ceil(math.log(Trh / Tev) / step)) + 1)
    hi, lo = edge[:-1], edge[1:]
    T = np.exp((hi + lo) / 2); H = hubble(T); a = scale(T)
    dN = hi - lo + np.log(gs(np.exp(hi)) / gs(np.exp(lo))) / 3
    ell = np.log(mr / H); rho_s = 4 * (0.24 * ell) * np.pi * fa ** 2 * ell * H ** 2
    dE = np.clip(eps_fn(ell), 0, None) * 2 * rho_s * dN
    qd, wt = np.polynomial.legendre.leggauss(bins)
    lo_p = math.log(x0); hi_p = np.log(mr / H)[:, None]
    X = np.exp((lo_p + hi_p) / 2 + qd[None, :] * (hi_p - lo_p) / 2); W = wt / 2
    Ncom = float(np.sum(dE[:, None] * W[None, :] / (X * H[:, None]) * a[:, None] ** 3))
    a_ev = float(scale(Tev))
    return Ncom / a_ev ** 3 / (fa ** 2 * H_eval), float(np.log(mr / H_eval))


def u7a(F, Trh):
    eps = lambda ell: 1.20 / (0.24 * np.asarray(ell))            # B2, U(1) (R = 1)
    tab = {}
    for x0 in (10, 30):
        vals = [nu_emission(F, Trh, eps, x0, h) for h in H_EVALS]
        nus = [v[0] for v in vals]; mean = float(np.mean(nus))
        tab[str(x0)] = dict(nu_S4=nus, ell=[v[1] for v in vals], nu_mean=mean, spread=float(np.std(nus) / mean),
                            kappa=NU_AX / mean, kappa_err=NU_AX_ERR / mean)
    return tab


def main(argv):
    from ucommon import benchmark
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    OUT7.mkdir(parents=True, exist_ok=True); t0 = time.time()
    out = dict(spec="docs/preregistros/ESPEC_U7.md", nu_ax_reference=dict(value=NU_AX, err=NU_AX_ERR, zeta=ZETA, source="arXiv:2512.13653"))
    out["U7a"] = u7a(F, Trh)
    for x0, v in out["U7a"].items():
        print(f"U7a x0={x0}: nu_S4 = {[round(n, 3) for n in v['nu_S4']]} (ell {[round(l, 1) for l in v['ell']]}), "
              f"kappa = {v['kappa']:.3f} +- {v['kappa_err']:.3f}", flush=True)
    if "--nu-only" in argv:
        (OUT7 / "U7A_RESULT.json").write_text(json.dumps(out, indent=1, default=float)); return out

    import u4a_class as U4
    import u5_class as U5
    import s6_measured as S6
    from s4_dark import lambda_wdm_reference
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else U4.default_workers()
    import classy
    out["classy_version"] = getattr(classy, "__version__", "?")
    P0 = U4.p_cdm()
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    Rv, Rs, Nr = U5.r_gamma(); Rset = {"-1s": max(Rv - Rs, 1e-6), "central": Rv, "+1s": Rv + Rs}
    out.update(R_Gamma=dict(value=Rv, total=Rs, source=f"U4b N={Nr}"), deltaA_strict=strict, deltaA_wide=wide)
    U7C = "--u7c" in argv
    if U7C:                                  # ESPEC_U7c: R_Gamma protocol spread from U6 (sparse) and U9 (dense, window scan)
        Rset = {"R1.39": 1.39, "R1.69": 1.69, "R2.09": 2.09}
        out.update(spec="docs/preregistros/ESPEC_U7c.md", R_Gamma=dict(values=Rset, source="U6 (xi_i=0.5), U9 (xi_i=2; window, max of window scan)"))
    s3 = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text()); dneff_late = s3["DeltaNeff_late_homogeneous_max"]
    lam_wdm, _ = lambda_wdm_reference()
    rows, params, idx = [], [], []
    for fam_name, (rl, R), x0, K in itertools.product(("B2k", "B2z"), Rset.items(), (10, 30), (0, 1)):
        k = out["U7a"][str(x0)]["kappa"]
        eps = (lambda ell, R=R, k=k: k * R * 1.20 / (0.24 * np.asarray(ell))) if fam_name == "B2k" else (lambda ell, R=R: R * ZETA / (0.24 * np.asarray(ell)))
        r = S6.node(F, Trh, dneff_late, lam_wdm, eps, "energy_logflat", x0, K)
        r.update(family=fam_name, R_label=rl, R=R, x0=x0, K_sw=K, kappa=(k if fam_name == "B2k" else None), tag=f"U7_{fam_name}_{rl}_x{x0}_K{K}")
        if not (r.get("P1") and all(r.get(q) for q in ("P2", "P3", "P4"))):
            r.update(LYA_primary=False, excluded_by="S6 P1-P4"); rows.append(r); continue
        pk = S6.source_eps(F, r["m_G_eV"], Trh, x0, eps, spectrum="energy_logflat")
        path, T, tail = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], r["tag"])
        idx.append(len(rows)); params.append(U4.packet_params(path, T, r["f_em"])); rows.append(r)
    print(f"{len(params)} CLASS row runs", flush=True)
    for i, P in zip(idx, U4.pmap(params, workers, [rows[i]["tag"] for i in idx], "row")):
        if P is None: rows[i].update(class_error=True, LYA_primary=False); continue
        rows[i].update(U5.gates(P / P0, P0, fam, strict, wide))
        r = rows[i]; print(f"{r['tag']:<26} f_em={r['f_em']:.3f} dA={r['deltaA']:.3f} primary={r['LYA_primary']}", flush=True)
    out["rows"] = rows
    ok = lambda f, rl: any(r.get("LYA_primary") for r in rows if r["family"] == f and r["R_label"] == rl)
    out["table"] = {f"{f}|{rl}": ok(f, rl) for f in ("B2k", "B2z") for rl in Rset}
    out["runtime_h"] = (time.time() - t0) / 3600
    if U7C:
        okx = lambda rl, x0: any(r.get("LYA_primary") for r in rows if r["family"] == "B2k" and r["R_label"] == rl and r["x0"] == x0)
        out["table_x0"] = {f"B2k|{rl}|x0={x0}": okx(rl, x0) for rl in Rset for x0 in (10, 30)}
        out["verdict_U7c"] = ("CALIBRADA_ROBUSTA" if okx("R2.09", 10) and okx("R2.09", 30) else
                              "CALIBRADA_SENSIBLE_A_R" if ok("B2k", "R1.69") else "CALIBRADA_EXCLUIDA_CON_REDES_DENSAS")
        (OUT7 / "U7C_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
        print("U7c:", out["verdict_U7c"], out["table_x0"], flush=True); return out
    out["verdict_U7"] = ("CALIBRADA_PERMITIDA" if ok("B2k", "+1s") else "CALIBRADA_EXCLUIDA" if not ok("B2k", "-1s") else "MARGINAL")
    (OUT7 / "U7_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print("U7:", out["verdict_U7"], out["table"], flush=True)
    return out


if __name__ == "__main__":
    main(sys.argv)
