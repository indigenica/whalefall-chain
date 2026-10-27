"""U28, informative reading decided after the data (INFORME_U28, note of 27 October): why the sampled cells of the runs that
pass the bounce-time test but fail f_loc land as 'none' (separate-universe classification of amendment 6, unchanged).

For every valid E-v and D-v run at L.H = 1 with pass_dt and not coherent, the run is repeated with the production code
(u28_c45k.one_run's steps) and each sampled cell gets a reason: 'open' (psi <= H^2), 'infertile' (no 60 e-folds or
|chi| <= 5), 'super_planck' (fertile but max(E_max, rho_max)/c2 >= 1), or its class E / D. The labels must reproduce the
stored su_labels exactly. Hypothesis to test (paper/zenodo_whale/README.md, item 3): eternal children lose cells mainly
above the Planck density, because K_max alpha = (2/pi^2) chi_land^2.

    python3 -B chain/u28_motivos.py  ->  results/u28_motivos.json   (~10-20 min)"""
import json, math, time
from collections import Counter
import numpy as np
import u28_c45k as K
import u28_comoving as C
from ucommon import ROOT, benchmark


def reasons(M, S, st, sel, rho_max):
    """the loop of u28_c45k.classify_cells, returning the reason of every sampled cell"""
    B = benchmark(); c2 = st["V0"] / B["V0_BENCH_C45"]
    r, U, G, phi, P, Y = S
    rp, eL, php, q, Pi, rho_d, rho, pr = M.derived(S)
    idx = np.nonzero(sel)[0]
    pick = idx[np.unique(np.linspace(0, len(idx) - 1, K.N_SU).round().astype(int))]
    out = []
    for i in pick:
        H = U[i] / r[i]; psi = rho[i] / 3 / math.sqrt(1 + (st["al"] * rho[i] / 3) ** 2)
        inv_a2 = psi - H * H
        if inv_a2 <= 0:
            out.append(dict(reason="open", label="open")); continue
        o = C.hom_continue(1 / math.sqrt(inv_a2), H, phi[i], Pi[i], rho_d[i], st)
        fert = o["outcome"] == "launched60" and o["chi_abs_max"] > 5
        e_over = max(o["E_max"], float(rho_max[i])) / c2
        if not fert:
            out.append(dict(reason="infertile", label="none", outcome=o["outcome"], chi_abs_max=float(o["chi_abs_max"]), E_over_planck=e_over))
        elif e_over >= 1:
            out.append(dict(reason="super_planck", label="none", E_over_planck=e_over, chi_end=float(o["chi_end"])))
        else:
            lab = "E" if abs(o["chi_end"]) > K.U27_CHI_ET else "D"
            out.append(dict(reason=lab, label=lab, E_over_planck=e_over, chi_end=float(o["chi_end"])))
    return out


def rerun(r, o):
    st = K.state_of(r)
    M, S, info = C.parent_data(st, K.N_PROD, o["family"], o["eps"], o["LH"] / r["AH"])
    ev = C.evolve(M, S, mode="bounce", t_max=4 * (info["t_half"] * 2), track_rho=True)
    return reasons(M, ev["S"], st, M.x <= K.X_HI, ev["rho_max"])


def main():
    ch, nfert, et = K.children(); K.U27_CHI_ET = et
    sel = json.loads((K.RUNS / "selection.json").read_text())
    rows = {x["tag"]: x for g in ("E-v", "D-v") for x in sel[g]}
    group = {x["tag"]: g for g in ("E-v", "D-v") for x in sel[g]}
    runs = [json.loads(f.read_text()) for f in sorted(K.RUNS.glob("sep*.json"))]
    todo = [o for o in runs if o["child"] in rows and o["LH"] == 1.0 and o.get("valid") and o.get("pass_dt") and not o["coherent"]]
    res, mismatches, t0 = [], 0, time.time()
    for o in todo:
        cells = rerun(rows[o["child"]], o)
        same = [c["label"] for c in cells] == o["su_labels"]
        mismatches += not same
        res.append(dict(child=o["child"], group=group[o["child"]], cls=o["cls"], family=o["family"], eps=o["eps"], f_loc=o["f_loc"],
                        labels_reproduced=same, cells=cells))
        print(group[o["child"]], o["child"], o["family"], o["eps"], Counter(c["reason"] for c in cells), "reproduced" if same else "MISMATCH", flush=True)
    summ = {}
    for g in ("E-v", "D-v"):
        cnt = Counter(c["reason"] for x in res if x["group"] == g for c in x["cells"])
        sp = [c["E_over_planck"] for x in res if x["group"] == g for c in x["cells"] if c["reason"] == "super_planck"]
        summ[g] = dict(n_runs=sum(x["group"] == g for x in res), reasons=dict(cnt),
                       super_planck_E_over_planck_median=float(np.median(sp)) if sp else None)
    out = dict(note="informative, decided after the data; INFORME_U28 dated note of 27 October", n_runs=len(res),
               labels_mismatches=mismatches, summary=summ, runs=res, runtime_s=time.time() - t0)
    (ROOT / "results" / "u28_motivos.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(dict(n_runs=len(res), labels_mismatches=mismatches, summary=summ), indent=1))


if __name__ == "__main__":
    main()
