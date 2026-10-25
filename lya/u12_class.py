"""U12 (docs/preregistros/ESPEC_U12.md): does the cold branch survive a growing string density?  Cross of the xi branches
B2 (control, xi* = 1.20), B4a / B4b (measured growth, s = 0.133 from U5d) and B1 (xi = 0.24 ell, GHV) with the spectral slope
q in {0.6, 0.75, 0.87, 1.0}, at the clean R_Gamma = 1.56 of U11.  Same CLASS set-up, references, shape family and gates as U5.

    python -B lya/u12_class.py --workers 32     # WSL, CLASS with the quadrature patch (64 rows; rows with f_em < 1e-4 skip CLASS)
    python -B lya/u12_class.py --quick          # plumbing check: x0 = 10, K_sw = 0 only -> U12_QUICK.json (no verdict)

Output: lya/results/u12/U12_RESULT.json (CLASS cache shared with U4a/U5: lya/results/u4a/cache)."""
import itertools, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402

OUT12 = HERE / "results" / "u12"
R_GAMMA, R_SOURCE = 1.56, "U11 (RP3 drained, R_Gamma = 1.562 +- 0.094)"
SLOPE = 0.133                                       # U5d: four U(1) seeds at 768^3
QS = (0.6, 0.75, 0.87, 1.0)
BRANCHES = ("B2", "B4a", "B4b", "B1")
GROWING = ("B4a", "B4b", "B1")


def eps_of(branch, R):
    if branch == "B2": return U5.eps_B2(R)
    if branch == "B4a": return U5.eps_B4(R, SLOPE, 1.0)
    if branch == "B4b": return U5.eps_B4(R, SLOPE, U5.XI_STAR / U5.XI_REF)
    if branch == "B1": return lambda ell: np.full_like(np.asarray(ell, dtype=float), R)
    raise ValueError(branch)


def row_specs(quick):
    x0s, Ks = ((10,), (0,)) if quick else ((10, 30), (0, 1))
    return [dict(branch=b, q=q, x0=x0, K_sw=K, R=R_GAMMA, spectrum=f"powerlaw_q{q:g}", eps=eps_of(b, R_GAMMA))
            for b, q, x0, K in itertools.product(BRANCHES, QS, x0s, Ks)]


def tag_of(r): return f"U12_{r['branch']}_q{r['q']:g}_x{r['x0']}_K{r['K_sw']}"


def status(rows, branch, q, key="LYA_primary"):
    n = sum(bool(r.get(key)) for r in rows if r["branch"] == branch and r["q"] == q)
    tot = sum(1 for r in rows if r["branch"] == branch and r["q"] == q)
    return dict(n=n, of=tot, state="PASA" if n == tot and tot else ("FALLA" if n == 0 else "PARCIAL"))


def verdict(tab):
    st = lambda b, q: tab[f"{b}|q={q:g}"]["state"] == "PASA"
    if all(st(b, 0.87) for b in GROWING): return "FRIA_SOBREVIVE_AL_CRECIMIENTO"
    if st("B4a", 0.87) and st("B4b", 0.87): return "FRIA_SOBREVIVE_AL_CRECIMIENTO_MEDIDO"
    if any(st("B4a", q) and st("B4b", q) for q in (0.6, 0.75)): return "SOLO_CON_q_BAJA"
    return "DEPENDE_DE_LA_ESTABILIZACION"


def q_crit(rows, branch, strict):
    pts = sorted((r["q"], r["deltaA"]) for r in rows if r["branch"] == branch and r["x0"] == 10 and r["K_sw"] == 0 and "deltaA" in r)
    for (a, da), (b, db) in zip(pts[:-1], pts[1:]):
        if (da - strict) * (db - strict) <= 0 and da != db: return a + (strict - da) * (b - a) / (db - da)
    return None


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else U4.default_workers()
    quick = "--quick" in argv; t0 = time.time(); OUT12.mkdir(parents=True, exist_ok=True)
    import classy
    import s6_measured as S6
    from s4_dark import lambda_wdm_reference
    from ucommon import benchmark
    out = dict(spec="docs/preregistros/ESPEC_U12.md", classy_version=getattr(classy, "__version__", "?"), quick=quick,
               R_Gamma=dict(value=R_GAMMA, source=R_SOURCE), slope_B4=SLOPE, qs=list(QS))
    print("classy", out["classy_version"], "workers", workers, flush=True)
    P0 = U4.p_cdm()

    # C0b: the quadrature patch must be present (rows use strategy 4)
    qg = np.geomspace(1e-3, 40, 800); U4.PSD.mkdir(parents=True, exist_ok=True); fd = U4.PSD / "C0_fermi_dirac.dat"
    np.savetxt(fd, np.column_stack([qg, 1 / (np.exp(qg) + 1)]), fmt="%.10e")
    pa = U4.thermal_params(3.5, 1.0); pq = dict(pa, use_ncdm_psd_files=1, ncdm_psd_filenames=str(fd), ncdm_quadrature_strategy=4, ncdm_N_momentum_bins=20)
    Pa, Pq = U4.pmap([pa, pq], min(workers, 2), ["C0 analytic", "C0b"]); s20 = (U4.KH >= 0.5) & (U4.KH <= 20)
    c0b = abs(U4.deltaA(Pq / P0, P0) / U4.deltaA(Pa / P0, P0) - 1) < 5e-3 and float(np.max(np.abs(Pq[s20] / Pa[s20] - 1))) < 5e-3
    out["C0b_passed"] = bool(c0b)
    if not c0b and out["classy_version"] != "MOCK":
        raise RuntimeError("C0b failed: CLASS without the quadrature patch (env/cloud_setup.sh or env/class_qm_trapz_log.patch)?")

    # references and extended shape family (as U5; cached)
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run: see lya/results/u4a/class_errors.log")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict, deltaA_wide=wide)
    print("refs", ref_dA, "strict", strict, flush=True)

    # S6 rows
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    s3 = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text()); dneff_late = s3["DeltaNeff_late_homogeneous_max"]
    lam_wdm, _ = lambda_wdm_reference()
    rows, params, meta = [], [], []
    for sp in row_specs(quick):
        r = S6.node(F, Trh, dneff_late, lam_wdm, sp["eps"], sp["spectrum"], sp["x0"], sp["K_sw"])
        r.update({k: v for k, v in sp.items() if k != "eps"}); r["tag"] = tag_of(sp)
        if not (r.get("P1") and all(r.get(k) for k in ("P2", "P3", "P4"))):
            r.update(LYA_primary=False, LYA_dA_strict=False, excluded_by="S6 P1-P4"); rows.append(r); continue
        if r["f_em"] < 1e-4:
            r.update(U5.gates(np.ones_like(U4.KH), P0, fam, strict, wide), class_skipped="f_em < 1e-4"); rows.append(r); continue
        pk = S6.source_eps(F, r["m_G_eV"], Trh, sp["x0"], sp["eps"], spectrum=sp["spectrum"])
        path, T, tail = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], r["tag"])
        r["tail_number_fraction"] = tail; meta.append(len(rows)); params.append(U4.packet_params(path, T, r["f_em"])); rows.append(r)
    print(f"{len(params)} CLASS row runs (of {len(rows)} rows)", flush=True)
    for i, P in zip(meta, U4.pmap(params, workers, [rows[i]["tag"] for i in meta], "row")):
        r = rows[i]
        if P is None:
            r.update(class_error=True, LYA_primary=False, LYA_dA_strict=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
    for r in rows:
        if "deltaA" in r:
            print(f"{r['tag']:<24} f_em={r['f_em']:.4f} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} "
                  f"primary={r['LYA_primary']} dA_only={r['LYA_dA_strict']}", flush=True)
    out["rows"] = rows

    out["table_primary"] = {f"{b}|q={q:g}": status(rows, b, q) for b in BRANCHES for q in QS}
    out["table_dA_only"] = {f"{b}|q={q:g}": status(rows, b, q, "LYA_dA_strict") for b in BRANCHES for q in QS}
    out["q_crit"] = {b: q_crit(rows, b, strict) for b in BRANCHES}
    out["cold"] = {f"{b}|q={q:g}": bool(all(r["f_em"] < 0.01 for r in rows if r["branch"] == b and r["q"] == q and "f_em" in r))
                   for b in BRANCHES for q in QS}
    out["verdict_U12"] = "QUICK (comprobación de la tubería, sin veredicto)" if quick else verdict(out["table_primary"])
    if out["classy_version"] == "MOCK": out["verdict_U12"] = "MOCK (sin validez)"
    out["runtime_h"] = (time.time() - t0) / 3600
    name = "U12_QUICK.json" if (quick or out["classy_version"] == "MOCK") else "U12_RESULT.json"
    (OUT12 / name).write_text(json.dumps(out, indent=1, default=float))
    print("U12:", out["verdict_U12"], flush=True)
    for b in BRANCHES:
        print(f"  {b:<4}", "  ".join(f"q={q:g}: {out['table_primary'][f'{b}|q={q:g}']['n']}/{out['table_primary'][f'{b}|q={q:g}']['of']}"
                                    f" (dA {out['table_dA_only'][f'{b}|q={q:g}']['n']})" for q in QS), "q_crit", out["q_crit"][b], flush=True)
    return out


if __name__ == "__main__":
    main(sys.argv)
