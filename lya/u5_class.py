"""U5 (docs/preregistros/ESPEC_U5.md, with amendment 1): does the dark sector survive a growing xi (U5a, branch B4) and the
measured spectral slope q (U5f)?  Plus the extended shape family (U5b), the white-noise isocurvature bound (U5c) and the
extra U(1) seeds (U5d, read if present).  Same CLASS set-up and gates as U4a (lya/u4a_class.py).

    python -B lya/u5_class.py --workers 16          # WSL, CLASS with env/class_qm_trapz_log.patch (~3 h with 16-32 cores)
    python -B lya/u5_class.py --quick               # x0 = 10, K_sw = 0, central R only (pipeline check)

Output: lya/results/u5/U5_RESULT.json (CLASS cache shared with U4a: lya/results/u4a/cache)."""
import itertools, json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402

OUT5 = HERE / "results" / "u5"
FAM_M_EXT = sorted([0.2, 0.3, 0.4] + U4.FAM_M); FAM_F_EXT = sorted([0.01, 0.02, 0.03] + U4.FAM_F)
XI_REF, ELL_REF, XI_STAR = 0.616, 5.65, 1.20          # U4b phase 1 (ESPEC_U5)
SLOPES = (0.062, 0.105, 0.181); QS = (0.9, 1.1, 1.2)
GEV_TO_INV_MPC = 1.5637e38                            # 1 GeV = 1.5637e38 Mpc^-1 (hbar c = 1.9733e-16 GeV m, 1 Mpc = 3.0857e22 m)


def eps_B2(R): return lambda ell: R * XI_STAR / (0.24 * ell)


def eps_B4(R, s, norm):
    return lambda ell: R * norm * np.maximum(XI_REF + s * (np.asarray(ell) - ELL_REF), 0.3) / (0.24 * np.asarray(ell))


def r_gamma():
    """R_Gamma at the attractor: largest-N U4b block with GATES_PASSED and RP3 runs (ESPEC_U5 amendment 1)"""
    b = json.loads((ROOT / "network" / "results" / "u4" / "U4B_RESULT.json").read_text())["blocks"]
    ok = [x for x in b if x.get("status") == "GATES_PASSED" and "R_Gamma_summary" in x]
    if not ok: raise RuntimeError("no U4b block with GATES_PASSED and R_Gamma: run network/u4b_analyze.py first")
    x = max(ok, key=lambda y: y["N"]); return x["R_Gamma_summary"]["value"], x["R_Gamma_summary"]["total"], x["N"]


def row_specs(Rv, Rs, quick):
    Rset = {"-1s": max(Rv - Rs, 1e-6), "central": Rv, "+1s": Rv + Rs}
    x0s, Ks = ((10,), (0,)) if quick else ((10, 30), (0, 1))
    Rlab = ("central",) if quick else tuple(Rset)
    out = []
    for nm, norm in (("B4a", 1.0), ("B4b", XI_STAR / XI_REF)):              # U5a
        for s, rl, x0, K in itertools.product(SLOPES, Rlab, x0s, Ks):
            out.append(dict(part="U5a", branch=nm, slope=s, R_label=rl, R=Rset[rl], spectrum="energy_logflat", q=1.0, x0=x0, K_sw=K,
                            eps=eps_B4(Rset[rl], s, norm)))
    for (nm, mk), q, rl, x0, K in itertools.product((("B2", None), ("B4a", 0.105)), QS, Rlab, x0s, Ks):     # U5f
        out.append(dict(part="U5f", branch=nm, slope=mk, R_label=rl, R=Rset[rl], spectrum=f"powerlaw_q{q:g}", q=q, x0=x0, K_sw=K,
                        eps=eps_B2(Rset[rl]) if nm == "B2" else eps_B4(Rset[rl], mk, 1.0)))
    for rl, x0, K in itertools.product(Rlab, x0s, Ks):                                                      # reference
        out.append(dict(part="ref", branch="B2", slope=None, R_label=rl, R=Rset[rl], spectrum="energy_logflat", q=1.0, x0=x0, K_sw=K,
                        eps=eps_B2(Rset[rl])))
    return out


def tag_of(r):
    sl = "" if r["slope"] is None else f"_s{r['slope']:g}"
    return f"U5_{r['branch']}{sl}_{r['R_label']}_q{r['q']:g}_x{r['x0']}_K{r['K_sw']}"


def gates(T2, P0, fam, strict, wide):
    dA = U4.deltaA(T2, P0); shp = U4.fit_shape(T2, fam)
    g = dict(deltaA=dA, shape=shp, LYA_dA_strict=bool(dA <= strict), LYA_dA_wide=bool(dA <= wide),
             LYA_shape_strict=bool(shp["f"] <= U4.f_max_strict(shp["m_keV"])), LYA_shape_wide=bool(shp["f"] <= U4.f_max_wide(shp["m_keV"])),
             shape_representable=bool(shp["rms"] <= 0.1), shape_extrapolated_below_1keV=bool(shp["m_keV"] < 1.0),
             shape_at_grid_edge=bool(shp["f"] <= fam["f"][0] + 1e-9 or shp["m_keV"] <= fam["m"][0] * (1 + 1e-9)))
    g["LYA_primary"] = g["LYA_shape_strict"] if g["shape_representable"] else g["LYA_dA_strict"]
    return g


def isocurvature(rows, P0):
    """U5c: white-noise isocurvature below k_* = a H at the end of emission (H = m/3), Delta^2 ~ f^2 (k/k_*)^3, frozen until
    equality and growing ~ (1 + 1.5 a/a_eq) after; compared with the adiabatic CDM power at z = 3, k = 20 h/Mpc."""
    from s4_dark import hubble, scale
    kh = 20.0; k = kh * U4.H; P = float(np.exp(np.interp(math.log(kh), np.log(U4.KH), np.log(P0)))); D2ad = k ** 3 * P / (2 * math.pi ** 2)
    grow = (1 + 1.5 * (1 / (1 + U4.Z)) * 3403.0) ** 2
    out = []
    for r in rows:
        T = r["T_stop_GeV"]; kstar = scale(T) * hubble(T) * GEV_TO_INV_MPC
        D2iso = 1.0 * (k / kstar) ** 3 * grow                                   # f_DM = 1 (conservative)
        out.append(dict(tag=r["tag"], k_star_invMpc=kstar, D2_iso=D2iso, D2_ad=D2ad, ratio=D2iso / D2ad))
    return dict(k_h_Mpc=kh, rows=out, max_ratio=max((x["ratio"] for x in out), default=0.0),
                passed=bool(all(x["ratio"] < 0.01 for x in out)))


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else U4.default_workers()
    quick = "--quick" in argv; t0 = time.time(); OUT5.mkdir(parents=True, exist_ok=True)
    import classy
    import s6_measured as S6
    from s4_dark import lambda_wdm_reference
    from ucommon import benchmark
    out = dict(spec="docs/preregistros/ESPEC_U5.md (+ amendment 1)", classy_version=getattr(classy, "__version__", "?"), quick=quick)
    print("classy", out["classy_version"], "workers", workers, flush=True)
    P0 = U4.p_cdm()

    # C0b again (the patch must be present: rows use strategy 4)
    q = np.geomspace(1e-3, 40, 800); U4.PSD.mkdir(parents=True, exist_ok=True); fd = U4.PSD / "C0_fermi_dirac.dat"
    np.savetxt(fd, np.column_stack([q, 1 / (np.exp(q) + 1)]), fmt="%.10e")
    pa = U4.thermal_params(3.5, 1.0); pq = dict(pa, use_ncdm_psd_files=1, ncdm_psd_filenames=str(fd), ncdm_quadrature_strategy=4, ncdm_N_momentum_bins=20)
    Pa, Pq = U4.pmap([pa, pq], min(workers, 2), ["C0 analytic", "C0b"]); s20 = (U4.KH >= 0.5) & (U4.KH <= 20)
    c0b = abs(U4.deltaA(Pq / P0, P0) / U4.deltaA(Pa / P0, P0) - 1) < 5e-3 and float(np.max(np.abs(Pq[s20] / Pa[s20] - 1))) < 5e-3
    if not c0b and out["classy_version"] != "MOCK":
        raise RuntimeError("C0b failed: CLASS without env/class_qm_trapz_log.patch?")

    # references (as U4a) and the extended family (U5b)
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = (FAM_M_EXT, FAM_F_EXT)
    labels = list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f]
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers, labels, "thermal")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run: see lya/results/u4a/class_errors.log")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}
    strict, wide = min(ref_dA.values()), max(ref_dA.values())
    T2f = np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1)
    fam_ext = dict(m=fam_m, f=fam_f, T2=T2f.tolist())
    im = [fam_m.index(m) for m in U4.FAM_M]; jf = [fam_f.index(f) for f in U4.FAM_F]
    fam_u4 = dict(m=U4.FAM_M, f=U4.FAM_F, T2=T2f[np.ix_(im, jf)].tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict, deltaA_wide=wide)
    print("refs", ref_dA, "strict", strict, "wide", wide, flush=True)

    # S6 rows for U5a / U5f / reference with the attractor R_Gamma
    Rv, Rs, Nr = r_gamma(); out["R_Gamma"] = dict(value=Rv, total=Rs, source=f"U4b N={Nr}")
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    s3 = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text()); dneff_late = s3["DeltaNeff_late_homogeneous_max"]
    lam_wdm, _ = lambda_wdm_reference()
    specs = row_specs(Rv, Rs, quick); rows, params, meta = [], [], []
    for sp in specs:
        r = S6.node(F, Trh, dneff_late, lam_wdm, sp["eps"], sp["spectrum"], sp["x0"], sp["K_sw"])
        r.update({k: v for k, v in sp.items() if k != "eps"}); r["tag"] = tag_of(sp)
        if not (r.get("P1") and all(r.get(k) for k in ("P2", "P3", "P4"))):
            r.update(LYA_primary=False, excluded_by="S6 P1-P4"); rows.append(r); continue
        if r["f_em"] < 1e-4:
            r.update(gates(np.ones_like(U4.KH), P0, fam_ext, strict, wide)); rows.append(r); continue
        pk = S6.source_eps(F, r["m_G_eV"], Trh, sp["x0"], sp["eps"], spectrum=sp["spectrum"])
        path, T, tail = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], r["tag"])
        r["tail_number_fraction"] = tail; meta.append(len(rows)); params.append(U4.packet_params(path, T, r["f_em"])); rows.append(r)
    print(f"{len(params)} CLASS row runs (of {len(specs)} rows)", flush=True)
    Prows = U4.pmap(params, workers, [rows[i]["tag"] for i in meta], "row")
    for i, P in zip(meta, Prows):
        r = rows[i]
        if P is None:
            r.update(class_error=True, LYA_primary=False); continue
        r.update(gates(P / P0, P0, fam_ext, strict, wide))
        r["shape_U4a_family"] = U4.fit_shape(P / P0, fam_u4)
        print(f"{r['tag']:<36} f_em={r['f_em']:.3f} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} "
              f"primary={r['LYA_primary']}", flush=True)
    out["rows"] = rows

    # U5b: refit the U4a rows (old S6 file, cached CLASS spectra) with the extended family
    old = ROOT / "results" / "s6_measured_closure_R_U2.json"; u4a_res = HERE / "results" / "u4a" / "U4A_RESULT.json"
    if old.exists() and u4a_res.exists() and not quick:
        from run_u3a import eps_fn_for
        U2 = json.loads((ROOT / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
        prev = {x["tag"]: x for x in json.loads(u4a_res.read_text())["rows"]}
        s6o = json.loads(old.read_text()); rb, pb = [], []
        for r in s6o["rows"]:
            if not (r.get("P1") and r["K_sw"] in (0, 1) and all(r.get(k) for k in ("P2", "P3", "P4"))) or r["f_em"] < 1e-4: continue
            tag = f"{r['branch'][:2]}_{r['R_label']}_{r['spectrum'][:6]}_x{r['x0']}_K{r['K_sw']}"
            pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
            path, T, _ = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], tag)
            rb.append(tag); pb.append(U4.packet_params(path, T, r["f_em"]))
        Pb = U4.pmap(pb, workers, rb, "U4a row (cache)")
        refit = []
        for tag, P in zip(rb, Pb):
            if P is None or tag not in prev: continue
            g = gates(P / P0, P0, fam_ext, strict, wide)
            refit.append(dict(tag=tag, U4a_primary=prev[tag]["LYA_primary"], U5b_primary=g["LYA_primary"], shape_ext=g["shape"],
                              shape_U4a=prev[tag]["shape"], extrapolated_below_1keV=g["shape_extrapolated_below_1keV"], at_edge=g["shape_at_grid_edge"]))
        out["U5b"] = dict(rows=refit, n=len(refit), n_changed=sum(x["U4a_primary"] != x["U5b_primary"] for x in refit),
                          n_still_at_edge=sum(x["at_edge"] for x in refit))
        print("U5b", {k: out["U5b"][k] for k in ("n", "n_changed", "n_still_at_edge")}, flush=True)

    # U5c
    out["U5c"] = isocurvature([r for r in rows if r.get("LYA_primary")], P0)
    print("U5c max ratio", out["U5c"]["max_ratio"], "passed", out["U5c"]["passed"], flush=True)

    # verdicts (ESPEC_U5 and amendment 1)
    def allowed(part, branch, slope=None, q=None, R_label=None):
        return any(r.get("LYA_primary") for r in rows if r["part"] in part and r["branch"] == branch and (slope is None or r["slope"] == slope)
                   and (q is None or r["q"] == q) and (R_label is None or r["R_label"] == R_label))
    def v3(ok_hi, ok_lo, yes, no):
        return yes if ok_hi else (no if not ok_lo else "MARGINAL")
    out["verdict_U5a"] = v3(allowed(("U5a",), "B4a", 0.105, R_label="+1s"), allowed(("U5a",), "B4a", 0.105, R_label="-1s"),
                            "SOBREVIVE_A_LA_PENDIENTE", "DEPENDE_DEL_ESCALADO")
    out["verdict_U5f"] = v3(allowed(("U5f",), "B2", q=1.2, R_label="central"), allowed(("U5f",), "B2", q=1.1, R_label="-1s"),
                            "SOBREVIVE_A_q", "DEPENDE_DE_q")
    out["U5a_table"] = {f"{b}|s={s:g}|{rl}": allowed(("U5a",), b, s, R_label=rl) for b in ("B4a", "B4b") for s in SLOPES for rl in ("-1s", "central", "+1s")}
    out["U5f_table"] = {f"{b}|q={q:g}|{rl}": allowed(("U5f", "ref", "U5a"), b, None if b == "B2" else 0.105, q, rl)
                        for b in ("B2", "B4a") for q in (0.9, 1.0, 1.1, 1.2) for rl in ("-1s", "central", "+1s")}

    def crit(part, branch, key, vals, fixed):
        pts = [(r[key], r["deltaA"]) for r in rows if r["part"] in part and r["branch"] == branch and r["R_label"] == "central" and r["x0"] == 10
               and r["K_sw"] == 0 and "deltaA" in r and all(r.get(k) == v for k, v in fixed.items())]
        pts = sorted(set(pts))
        for (a, da), (b, db) in zip(pts[:-1], pts[1:]):
            if (da - strict) * (db - strict) <= 0 and da != db: return a + (strict - da) * (b - a) / (db - da)
        return None
    out["s_crit_B4a"] = crit(("U5a",), "B4a", "slope", SLOPES, {})
    out["q_crit_B2"] = crit(("U5f", "ref"), "B2", "q", QS, {})
    if quick: out["verdict_U5a"] = out["verdict_U5f"] = "QUICK (comprobación de la tubería)"

    # U5d (if the extra seeds exist)
    runs = sorted((ROOT / "network" / "results" / "u4" / "runs").glob("U1_N768_s6*_xi0.5.json"))
    sl = []
    for p in runs:
        R = json.loads(p.read_text()); l = np.array([x["ell"] for x in R["rows"]]); x = np.array([x["xi"] for x in R["rows"]])
        m = (l >= ELL_REF - 0.7) & (l <= ELL_REF); sl.append(dict(run=p.name, slope=float(np.polyfit(l[m], x[m], 1)[0]), xi_ref=float(np.interp(ELL_REF, l, x))))
    out["U5d"] = dict(runs=sl, mean_slope=float(np.mean([s["slope"] for s in sl])) if sl else None,
                      sem=float(np.std([s["slope"] for s in sl], ddof=1) / math.sqrt(len(sl))) if len(sl) > 1 else None)
    out["runtime_h"] = (time.time() - t0) / 3600
    (OUT5 / "U5_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print("U5a:", out["verdict_U5a"], "s_crit", out["s_crit_B4a"], "| U5f:", out["verdict_U5f"], "q_crit", out["q_crit_B2"], "| U5c:", out["U5c"]["passed"], flush=True)
    return out


if __name__ == "__main__":
    main(sys.argv)
