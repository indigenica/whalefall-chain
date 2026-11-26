"""U70 analysis (docs/preregistros/ESPEC_U70.md): gates G1-G5, S_min(l_QT) over n, F and the two versions among the rows with a
throat, and the verdict decided at l_QT = 300 l_P.  G4 re-solves three rows chosen in the specification with tolerances 10x tighter.
    python3 -B chain/u70_analyze.py  ->  results/u70/U70_RESULT.json"""
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import u70_agujeros as U

OUT = U.OUT


def main():
    U.VAL = U.build_valley()
    rows = json.loads((OUT / "rows.json").read_text())
    Sreq = {ff: U.S_required(ff) for ff in U.F_ROWS}
    Sreq_max = max(Sreq.values())

    # G1: alpha = 0, (C), against 3 pi^3 a0^2
    g1 = []
    for ff in U.F_ROWS:
        a0, Sa = U.gs_analytic(ff, 1.0)
        r = U.task((ff, 1.0, 0.0, False, 1.0))
        g1.append(dict(F_over_F0=ff, S=r["S"], S_analytic=Sa, rel_err=abs(r["S"] / Sa - 1)))
    G1 = all(x["rel_err"] < 0.01 for x in g1)
    # G2: alpha = 0, (D), n = 1
    d = {r["F_over_F0"]: r for r in rows if r["l_QT"] == 0 and r["version"] == "D" and r["n"] == 1.0 and r.get("throat") and r.get("resolved")}
    c = {r["F_over_F0"]: r for r in rows if r["l_QT"] == 0 and r["version"] == "C" and r["n"] == 1.0 and r.get("throat") and r.get("resolved")}
    ok2 = len(d) == 3 and len(c) == 3
    g2 = dict(rows=len(d))
    if ok2:
        lnif = np.array([math.log(1.0 / (0.5 * ff)) for ff in U.F_ROWS]); S = np.array([d[ff]["S"] for ff in U.F_ROWS])
        slope = float(np.polyfit(lnif, S, 1)[0])
        ratio = max(d[ff]["S"] / c[ff]["S"] for ff in U.F_ROWS)
        g2.update(S_D=S.tolist(), slope=slope, max_ratio_D_over_C=ratio)
        ok2 = ratio < 0.10 and 0.3 <= slope <= 3.0 and bool(np.all(np.diff(S) > 0))
    # G3
    G3err = U.gate_G3(np.random.default_rng(70)); G3 = G3err < 1e-8
    # G4: constraint along every solution, and three rows re-solved with tol/10
    found = [r for r in rows if r.get("throat") and r.get("resolved")]
    unresolved = [dict(l_QT=r["l_QT"], F_over_F0=r["F_over_F0"], n=r["n"], version=r["version"], tail_fraction=r.get("tail_fraction"))
                  for r in rows if r.get("throat") and not r.get("resolved")]
    max_res = max(r["constraint"] for r in found)
    by_l = {}
    for l in U.LQT_ROWS:
        cand = [r for r in found if r["l_QT"] == l]
        by_l[l] = min(cand, key=lambda r: r["S"]) if cand else None
    picks = [by_l[300.0],
             next(r for r in found if r["l_QT"] == 0 and r["version"] == "C" and r["n"] == 1.0 and r["F_over_F0"] == 0.5),
             next(r for r in found if r["l_QT"] == 0 and r["version"] == "D" and r["n"] == 1.0 and r["F_over_F0"] == 0.5)]
    g4 = []
    for p in picks:
        t = U.task((p["F_over_F0"], p["n"], p["l_QT"], p["version"] == "D", 0.1))
        g4.append(dict(l_QT=p["l_QT"], F_over_F0=p["F_over_F0"], n=p["n"], version=p["version"], S=p["S"], S_tight=t.get("S"),
                       rel_change=abs(t["S"] / p["S"] - 1) if t.get("throat") else None))
    G4 = max_res < 1e-6 and all(x["rel_change"] is not None and x["rel_change"] < 0.01 for x in g4)
    # G5
    viol = [r for r in found if r["l_QT"] > 0 and r["a0"] < r["l_QT"]]
    G5 = len(viol) == 0

    gates = dict(G1=dict(passed=G1, rows=g1), G2=dict(passed=ok2, **g2), G3=dict(passed=G3, max_rel_err=G3err),
                 G4=dict(passed=G4, max_constraint=max_res, rows=g4), G5=dict(passed=G5, violations=len(viol)))
    S300 = by_l[300.0]["S"] if by_l[300.0] else None
    if not all(g["passed"] for g in gates.values()):
        verdict = "INCONCLUSO"
    elif S300 is None:
        verdict = "SIN GARGANTA (ver la regla de la cota en ESPEC §5)"
    elif S300 >= 2 * Sreq_max:
        verdict = "PROTEGE"
    elif S300 < Sreq_max:
        verdict = "NO PROTEGE"
    else:
        verdict = "INTERMEDIA"

    # informative: existence of dynamical (D) solutions with NPQT, minimal charge with a throat, S_min per version
    info = {}
    for l in U.LQT_ROWS:
        for v in ("C", "D"):
            cand = [r for r in found if r["l_QT"] == l and r["version"] == v]
            info[f"{l:g}_{v}"] = dict(n_rows_with_throat=len(cand),
                                     S_min=min((r["S"] for r in cand), default=None),
                                     at=min(cand, key=lambda r: r["S"]) if cand else None)
    res = dict(espec="docs/preregistros/ESPEC_U70.md", S_required=Sreq, S_required_max=Sreq_max, gates=gates,
               S_min={f"{l:g}": (by_l[l]["S"] if by_l[l] else None) for l in U.LQT_ROWS},
               S_min_row={f"{l:g}": by_l[l] for l in U.LQT_ROWS}, verdict=verdict, informative=info, n_rows=len(rows), unresolved=unresolved)
    (OUT / "U70_RESULT.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(dict(gates={k: v["passed"] for k, v in gates.items()}, S_min=res["S_min"], S_req_max=Sreq_max, verdict=verdict,
                          D_rows={k: v["n_rows_with_throat"] for k, v in info.items()}, n_unresolved=len(unresolved)), indent=1, default=float))


if __name__ == "__main__":
    main()
