"""U17e (docs/preregistros/ESPEC_U17e.md): the collapse population measured in the RP3 lattice (U19b, l_G = 3.2, four primary
networks) passed through CLASS on the U21 F grid.  Soft excess (k < 10 m_G a_1) minus the quanta radiated before ell_1, minus the
misalignment taken from the lowest-k bins (conservative), number-weighted shape w(x), x = k/(m_G a_1); in the universe q = x m_G a_osc.
Rows Y (Y = 1.70, U17 scaling K = (8 pi/27) R_rho xi ell Y/xbar), Y+ (Y = 2.24), X (K = measured soft number, reading (b) of
ESPEC_U19b enmienda 4); informative A (B4a, Y), Yp (misalignment removed in proportion), Xd (X with the lattice hard quanta).
Gates: U12 Lya, own isocurvature (U20 estimator, C = 1, r = r_N = 2.40), DeltaNeff <= 0.3, CMB strings.

    python -B lya/u17e_class.py --dry           # shape and rows, no CLASS
    python -B lya/u17e_class.py --workers 4     # cloud (CLASS cache of U12/U21 for references and family)
Output: lya/results/u17e/U17E_RESULT.json"""
import glob, json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain")); sys.path.insert(0, str(ROOT / "network"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402

OUT = HERE / "results" / "u17e"
FGRID = (1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03)
Q, X0, R_RHO = 0.87, 10, 1.56
Y_C, Y_HI = 1.70, 2.24                                                         # U19b block: 1.70 +- 0.55 (rounded as preregistered)
R_ISO, R_ISO_ALT = 2.40, 1.78                                                  # r_N (primary, latest) and r_W of U19b
ISO_C, ISO_MAX, DNEFF_MAX, F_CMB, LISA_FA, COLD = 1.0, 0.10, 0.3, 1.26e15, 5e14, 0.05
INFO_F = {"Yp": (0.5, 0.35), "Xd": (1.0, 0.5)}


def measured_shape(mode="low", hard=False):
    """ESPEC_U17e: number-weighted shape of the collapse excess of the four primary U19b networks (equal weight).
    mode 'low': misalignment (N_mis quanta) removed from the lowest-k bins; 'prop': removed in proportion to the soft spectrum.
    hard: keep the lattice's k >= 10 m_G a_1 quanta (minus their pre-ell_1 part) in the population."""
    import u19b_analyze as A
    val = json.loads((A.D / "U19B_VALIDATION.json").read_text())
    W, per = [], []
    for p in sorted(glob.glob(str(A.D / "runs" / "RP3_N512_lG3.2_*.json"))):
        run = json.loads(Path(p).read_text()); r = A.per_run(run, val); m = run["m_G"]; a1 = run["a_1"]
        meas = [x for x in run["rows"] if not x.get("initial")]; last = meas[-1]
        e = np.array(run["k_edges"]); n = np.array(last["lin_spectrum"]["n"])
        sk = [x for x in meas if "spectrum" in x]; ls = np.array([x["ell"] for x in sk])
        def cum(x):                                       # cumulative kinetic number of measure(), on the counter's edges
            ee = np.array(x["spectrum"]["k_edges"]); kk = 0.5 * (ee[1:] + ee[:-1]); w1 = np.sqrt(kk ** 2 / x["a"] ** 2 + m ** 2)
            nn = x["a"] ** 3 * 2 * np.array(x["spectrum"]["dE"]) / w1
            return np.interp(e, ee, np.concatenate([[0.0], np.cumsum(nn)]))
        Cm = np.array([cum(x) for x in sk]); j = int(np.searchsorted(ls, r["ell_1"])); t = (r["ell_1"] - ls[j - 1]) / (ls[j] - ls[j - 1])
        exc = n - np.diff((1 - t) * Cm[j - 1] + t * Cm[j])
        x = 0.5 * (e[1:] + e[:-1]) / (m * a1); keep = np.ones_like(x, bool) if hard else x < 10
        ex = exc.copy(); soft = x < 10
        if mode == "low":
            left = r["N_mis"]
            for i in range(len(ex)):
                if left <= 0: break
                take = min(max(ex[i], 0.0), left); ex[i] -= take; left -= take
        else:
            pos = np.clip(ex, 0, None) * soft; ex = ex - r["N_mis"] * pos / pos.sum()
        neg = float(-ex[keep & (ex < 0)].sum() / r["N_mis"]); ex = np.clip(ex, 0, None) * keep
        K = float(ex.sum() / r["N_mis"]); w = ex / ex.sum()
        W.append(w); per.append(dict(tag=r["tag"], K=K, xbar=float((w * x).sum()), negative_clipped=neg, X_lat=r["X_lat"]))
    w = np.mean(W, axis=0); sel = w > 0
    return dict(mode=mode, hard=hard, x=x[sel].tolist(), w=w[sel].tolist(), xbar=float((w * x).sum()),
                K=float(np.mean([q["K"] for q in per])), runs=per)


def row(F, Trh, branch, shape, kind, Y=None):
    import c35b_pseudo_goldstone_dm_rp3 as C35
    import u17_ksw as U17
    from s4_dark import abundance, scale, RHOC_H2
    from s6_measured import source_eps
    xb = shape["xbar"]
    s = U17.solve(F, Trh, branch, Q, X0, c_eta=xb / Y, R_rho=R_RHO) if Y is not None else U17.solve(F, Trh, branch, Q, X0, K_fixed=shape["K"])
    m_eV = s["m_G_eV"]; m = m_eV * 1e-9; K = s["K"]
    om_mis = C35.omega(m_eV, F, 0.0)[0]; om_coll = om_mis * K
    pk = source_eps(F, m_eV, Trh, X0, U17.eps_fn(branch), spectrum=f"powerlaw_q{Q:g}"); om_em = abundance(pk)
    a_osc = float(scale(s["T_osc_MeV"] * 1e-3))
    q_coll = np.array(shape["x"]) * m * a_osc; n_coll = om_coll * RHOC_H2 / m * np.array(shape["w"])
    q_all = np.concatenate([pk["q"], q_coll]); n_all = np.concatenate([pk["nc"], n_coll])
    return dict(kind=kind, branch=branch, Y=Y, K=K, xbar=xb, m_G_eV=m_eV, T_osc_MeV=s["T_osc_MeV"], a_osc=a_osc, Omega_mis_bare=om_mis,
                Omega_collapse=om_coll, Omega_em=om_em, f_nc=(om_em + om_coll) / 0.12, DeltaNeff_emission=s["DeltaNeff_CMB_emission"],
                ell_osc=s["ell_osc"]), q_all / m, n_all


def isocurvature(rows, P0):
    import u20_iso as U20
    k = U4.KH * U20.H_LITTLE; D2ad = float(np.interp(math.log(U20.KH), np.log(U4.KH), k ** 3 * P0 / (2 * math.pi ** 2)))
    D_iso = 1 + 1.5 * U20.G_B * (1 + U20.Z_EQ) / (1 + U20.Z); kk = U20.KH * U20.H_LITTLE
    for r in rows:
        rat = {f"r{rr:g}": ISO_C * (kk / U20.kstar(r["T_osc_MeV"] * 1e-3, rr)) ** 3 * D_iso ** 2 / D2ad for rr in (R_ISO, R_ISO_ALT)}
        r.update(iso_ratio=rat, iso_ok=bool(rat[f"r{R_ISO:g}"] <= ISO_MAX), iso_ok_rW=bool(rat[f"r{R_ISO_ALT:g}"] <= ISO_MAX))
    return dict(Delta2_ad_k10=D2ad, D_iso=D_iso)


def specs(shapes):
    out = [(k, "B2", shapes["low"], Y) for k, Y in (("Y", Y_C), ("Y+", Y_HI), ("X", None))] + [("A", "B4a", shapes["low"], Y_C)]
    return [(k, b, sh, Y, fr) for k, b, sh, Y in out for fr in FGRID] + \
           [("Yp", "B2", shapes["prop"], Y_C, fr) for fr in INFO_F["Yp"]] + [("Xd", "B2", shapes["hard"], None, fr) for fr in INFO_F["Xd"]]


def verdict(fmax):
    if fmax is None: return "SIN_BANDA"
    return "BANDA_U21_CONFIRMADA" if fmax >= 0.5 else "BANDA_RECORTADA" if fmax >= 0.18 else "BANDA_MINIMA"


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    from ucommon import benchmark
    from s4_dark import lambda_fs, lambda_wdm_reference
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; lam_wdm, _ = lambda_wdm_reference()
    shapes = {"low": measured_shape("low"), "prop": measured_shape("prop"), "hard": measured_shape("low", hard=True)}
    for k, s in shapes.items(): print(f"shape {k}: xbar={s['xbar']:.2f} K={s['K']:.3f} runs={[(q['tag'], round(q['K'], 3), round(q['xbar'], 2)) for q in s['runs']]}")
    rows = []
    for kind, br, sh, Y, fr in specs(shapes):
        F = F0 * fr; r, x, n = row(F, Trh, br, sh, kind, Y)
        lam = float(lambda_fs(np.array([r["xbar"] * r["m_G_eV"] * r["a_osc"]]), r["m_G_eV"], np.array([r["a_osc"]]))[0]) / lam_wdm
        r.update(F_over_F0=fr, F_GeV=F, lambda_proxy_over_WDM=lam, cold_trivial=bool(lam < COLD), f_a_GeV=F / 2, LISA_ok=bool(F / 2 >= LISA_FA),
                 dneff_ok=bool(r["DeltaNeff_emission"] <= DNEFF_MAX), cmb_strings_ok=bool(F <= F_CMB), tag=f"U17e_{kind}_F{fr:g}")
        r["_xn"] = (x, n); rows.append(r)
        print(f"{r['tag']:<16} K={r['K']:8.3f} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.3g} MeV lam/WDM={lam:.3g} cold={r['cold_trivial']} "
              f"f_nc={r['f_nc']:.3f}", flush=True)
    ctl = next((r for r in rows if r["kind"] == "Y" and r["cold_trivial"]), None)
    need = [r for r in rows if not r["cold_trivial"]] + ([ctl] if ctl else [])
    print(f"{len(need)} CLASS rows (control {ctl['tag'] if ctl else None})", flush=True)
    if dry: return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U17e.md", classy_version=getattr(classy, "__version__", "?"), F0_GeV=F0, grid=FGRID,
               shapes={k: {kk: v for kk, v in s.items()} for k, s in shapes.items()})
    P0 = U4.p_cdm()
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict); out["isocurvature"] = isocurvature(rows, P0)
    params = []
    for r in need:
        x, n = r["_xn"]; path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
        params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
    for r, P in zip(need, U4.pmap(params, workers, [r["tag"] for r in need], "row")):
        if P is None: r.update(class_error=True, LYA_primary=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
        print(f"{r['tag']:<16} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    control_ok = bool(ctl is None or ("deltaA" in ctl and ctl["deltaA"] < COLD))
    out["cold_control"] = dict(tag=ctl["tag"] if ctl else None, deltaA=ctl.get("deltaA") if ctl else None, ok=control_ok)
    if not control_ok:
        (OUT / "U17E_PARTIAL.json").write_text(json.dumps(dict(out, rows=[{k: v for k, v in r.items() if k != "_xn"} for r in rows]), indent=1, default=float))
        raise SystemExit("cold-trivial control failed: re-run with every row through CLASS (ESPEC_U17e)")
    for r in rows:
        r.pop("_xn", None)
        if r["cold_trivial"] and "LYA_primary" not in r: r["LYA_primary"] = True; r["LYA_by"] = "frio_trivial"
        r["passes"] = bool(r.get("LYA_primary") and r["iso_ok"] and r["dneff_ok"] and r["cmb_strings_ok"])
    out["rows"] = rows
    bands = {}
    for k in ("Y", "Y+", "X", "A", "Yp", "Xd"):
        rk = [r for r in rows if r["kind"] == k]; ok = [r["F_over_F0"] for r in rk if r["passes"]]
        bands[k] = dict(F_max_over_F0=max(ok) if ok else None, F_min_over_F0=min(ok) if ok else None, n_pass=len(ok), n_rows=len(rk),
                        contiguous=bool(ok and all(r["passes"] for r in rk if min(ok) <= r["F_over_F0"] <= max(ok))),
                        F_max_lya_only=max([r["F_over_F0"] for r in rk if r.get("LYA_primary")], default=None))
    out["bands"] = bands
    out["verdict_U17e"] = verdict(bands["Y"]["F_max_over_F0"])
    out["X_rows_note"] = ("BENCHMARK_PERMITIDO_SI_(b)" if bands["X"]["F_max_over_F0"] == 1.0 else "el benchmark no pasa ni con (b)") + \
                         " (solo aplica si la enmienda 4 de U19b da la lectura (b) con l_G = 2.4)"
    out["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U17E_RESULT.json").write_text(json.dumps(out, indent=1, default=float, allow_nan=False))
    print("U17e:", out["verdict_U17e"], out["X_rows_note"], bands, flush=True)


if __name__ == "__main__":
    main(sys.argv)
