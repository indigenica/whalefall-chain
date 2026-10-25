"""U17c (docs/preregistros/ESPEC_U17c.md): eta scan of the string-wall collapse population for the growing-xi branches.
Rows: B4a, B1 x c in {1, 3} x eta in {1, 0.3, 0.1, 0.03, 0.01} x sigma_ln in {0.5, 1.0}, plus B2 (c = 1, eta = 1); q = 0.87,
x0 = 10, R_Gamma = 1.56.  Row construction of U17b (lya/u17b_class.row), U12 gates.

    python -B lya/u17c_class.py --workers 16      # WSL (reuses the U12 CLASS cache for references and family)
    python -B lya/u17c_class.py --dry             # rows only, no CLASS
Output: lya/results/u17c/U17C_RESULT.json"""
import itertools, json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402
import u17b_class as U17B                                                     # noqa: E402

OUT = HERE / "results" / "u17c"
ETAS = (1.0, 0.3, 0.1, 0.03, 0.01); CS = (1.0, 3.0); SIGS = (0.5, 1.0)


def specs():
    out = [("B4a", c, e, s) for c, e, s in itertools.product(CS, ETAS, SIGS)]
    out += [("B1", c, e, s) for c, e, s in itertools.product(CS, ETAS, SIGS)]
    return out + [("B2", 1.0, 1.0, 0.5)]


def eta_crit(rows, branch, c, sig, strict):
    rr = sorted([r for r in rows if r["branch"] == branch and r["c"] == c and r["sigma_ln"] == sig], key=lambda r: -r["eta"])
    ok = [r for r in rr if r.get("LYA_primary")]
    if not ok: return None
    e_ok = max(r["eta"] for r in ok)
    above = [r for r in rr if r["eta"] > e_ok and "deltaA" in r]
    if not above: return e_ok
    hi = min(above, key=lambda r: r["eta"]); lo = [r for r in rr if r["eta"] == e_ok][0]
    if hi["deltaA"] == lo["deltaA"]: return e_ok
    t = (strict - lo["deltaA"]) / (hi["deltaA"] - lo["deltaA"])          # interpolation of deltaA in ln eta
    return float(math.exp(math.log(e_ok) + min(max(t, 0.0), 1.0) * (math.log(hi["eta"]) - math.log(e_ok))))


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    from ucommon import benchmark
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    rows, params = [], []
    for b, c, e, s in specs():
        U17B.SIG_LN = s
        r, x, n = U17B.row(F, Trh, b, c, e, legacy_energy=True); r["sigma_ln"] = s     # as run (energy/momentum mixed: U17d); r["tag"] = f"U17c_{b}_c{c:g}_eta{e:g}_s{s:g}"
        print(f"{r['tag']:<28} K={r['K_eff']:8.2f} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.4f} MeV f_nc={r['f_nc']:.4f} "
              f"M_mc={r['M_minicluster_Msun']:.3g}", flush=True)
        if not dry:
            path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
            params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
        rows.append(r)
    if dry: return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U17c.md", classy_version=getattr(classy, "__version__", "?"))
    P0 = U4.p_cdm()
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict)
    for r, P in zip(rows, U4.pmap(params, workers, [r["tag"] for r in rows], "row")):
        if P is None: r.update(class_error=True, LYA_primary=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
        print(f"{r['tag']:<28} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    out["rows"] = rows
    out["eta_crit"] = {f"{b}|c={c:g}|sigma={s:g}": eta_crit(rows, b, c, s, strict) for b in ("B4a", "B1") for c in CS for s in SIGS}
    ec = out["eta_crit"]["B4a|c=3|sigma=0.5"]
    out["verdict_U17c"] = ("RAMA_COMODA_MUERTA" if ec is None or ec < 0.03 else "RAMA_COMODA_FORZADA" if ec < 0.3 else "RAMA_COMODA_RAZONABLE")
    out["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U17C_RESULT.json").write_text(json.dumps(out, indent=1, default=float, allow_nan=False))
    print("U17c:", out["verdict_U17c"], out["eta_crit"], flush=True)


if __name__ == "__main__":
    main(sys.argv)
