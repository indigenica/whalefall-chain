"""C45j (ESPEC_C45IJ.md): do the eternal sub-Planckian children survive the multi-bounce desynchronisation with realistic
perturbations?  For each eternal-sub child of C45h (dust, window separations): eps_est at its (A.H, alpha H^2) from the LTB-QT
bounce-time spread of a recollapsing closed parent (as C45d), eps*_n = eps_est / n_b, zeta0_max = eps*_n / [(2/3) pi^2 (AH/LH)^2]."""
import glob, json, math, time
import numpy as np
from common.provenance import new_run_dir, write_json, write_manifest, source_hashes, environment
from common.potential import H_INF
from ltb.growing_qt import build_qt

SEPS = [1e-4, 1e-5, 1e-6]
LHS = [0.3, 1.0, 3.0]
ZETA0 = 4.6e-5
EPS = 1e-4


def eps_est(AH, aH2, LH):
    al = aH2 / H_INF ** 2                      # LTB units (G = 1, m_chi = 1): alpha H^2 is dimensionless and unit independent
    x = np.linspace(0, 0.95, 20001)
    K2 = x ** 2 * (1 + EPS * np.sinc(x / (LH / AH)))
    sh, tw, xs, xm, mu, AHx = build_qt(x, K2, AH, al, n_sh=600)
    tB = np.array([s.t_B for s in sh])
    return math.sqrt(al) * EPS / (tB.max() - tB.min()), AHx


def main():
    d = new_run_dir("c45j_multibounce")
    src = sorted(glob.glob("results/c45h_content_2*"))[-1]
    H = json.load(open(src + "/results.json"))["rows"]
    rows = [r for r in H if r["frac_rad"] == 0.0 and r["sep"] in SEPS]
    write_json(d / "PROTOCOL.json", dict(spec="ESPEC_C45IJ.md (C45j)", input=src, sources=source_hashes(), env=environment()))
    cache, out = {}, []
    t0 = time.time()
    for r in rows:
        if not r["eternal_sub"]:
            continue
        AH = 0.385 * math.sqrt(r["lam"]); aH2 = r["sep"] * 4 * r["lam"] / 27
        rec = dict(sep=r["sep"], lam=r["lam"], sigma=r["sigma"], AH=AH, alphaH2=aH2, n_b=r["n_bounces"], chi_land=r["chi_land"])
        for LH in LHS:
            key = (round(AH, 6), round(aH2, 12), LH)
            if key not in cache:
                cache[key] = eps_est(AH, aH2, LH)
            e, _ = cache[key]
            en = e / max(r["n_bounces"], 1)
            rec[f"eps_est_LH{LH:g}"] = e
            rec[f"zeta0_max_LH{LH:g}"] = en / ((2 / 3) * math.pi ** 2 * (AH / LH) ** 2)
        out.append(rec)
    print(f"{len(out)} children, {len(cache)} LTB evaluations, {time.time()-t0:.0f}s", flush=True)
    summ = {}
    for sep in SEPS:
        fe = [r for r in rows if r["sep"] == sep and r["fertile"]]
        ch = [o for o in out if o["sep"] == sep]
        n = max(len(fe), 1)
        summ[f"{sep:g}"] = dict(n_fertile=len(fe), n_eternal_sub=len(ch),
                                **{f"f_viable_LH{LH:g}": sum(o[f"zeta0_max_LH{LH:g}"] >= ZETA0 for o in ch) / n for LH in LHS},
                                **{f"median_zeta0_max_LH{LH:g}": float(np.median([o[f"zeta0_max_LH{LH:g}"] for o in ch])) if ch else None for LH in LHS},
                                n_b=sorted(o["n_b"] for o in ch))
    nf = sum(v["n_fertile"] for v in summ.values())
    fv = sum(o["zeta0_max_LH1"] >= ZETA0 for o in out) / max(nf, 1)
    verdict = "VIABLE" if fv >= 0.05 else ("MARGINAL" if fv >= 0.01 else "NO_VIABLE")
    allz = {f"LH{LH:g}": dict(p10=float(np.percentile([o[f"zeta0_max_LH{LH:g}"] for o in out], 10)),
                              median=float(np.median([o[f"zeta0_max_LH{LH:g}"] for o in out])),
                              p90=float(np.percentile([o[f"zeta0_max_LH{LH:g}"] for o in out], 90)),
                              max=float(max(o[f"zeta0_max_LH{LH:g}"] for o in out))) for LH in LHS}
    write_json(d / "results.json", dict(verdict=verdict, frac_viable_LH1_over_fertile=fv, summary=summ, zeta0_max_distribution=allz, children=out))
    write_manifest(d)
    print("VERDICT", verdict, round(fv, 4)); print(json.dumps(allz, indent=1))
    for k, v in summ.items():
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "n_b"}, "n_b:", v["n_b"])


if __name__ == "__main__":
    main()
