"""U57b (docs/preregistros/ESPEC_U57b.md): U57 repeated with new data (the 64 log-midpoints of the U57 f grid), G1 only on
the children, and the children at the edge of the branch (min W_rel < 1e-2) set apart.  explora_cuscuton is imported unchanged.
    python3 -B chain/u57b_cuscuton_fondo.py   ->  results/u57b_runs/rows_rtol{11,12}.json, results/u57b_cuscuton_fondo.json
Internal, not exported."""
import json, math, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
from ucommon import ROOT, save
import explora_cuscuton as EC
import u57_cuscuton_fondo as U57

RUNS = ROOT / "results" / "u57b_runs"
FS = tuple(10 ** ((k + 0.5) / 16) for k in range(-64, 0))
EDGE = 1e-2


def produce(rtol):
    path = RUNS / f"rows_rtol{-int(round(math.log10(rtol)))}.json"
    if path.exists():
        return json.loads(path.read_text())["rows"]
    jobs = [(name, dict(p, AH=A, f=f, rtol=rtol)) for name, p in U57.SETS.items() for A in U57.AHS for f in FS]
    t0 = time.time()
    with Pool(4) as pool:
        rows = pool.map(EC.job, jobs, chunksize=1)
    path.write_text(json.dumps(dict(rtol=rtol, rows=rows, seconds=time.time() - t0), default=float))
    print(f"rtol {rtol:g}: {len(rows)} rows in {time.time() - t0:.0f}s", flush=True)
    return rows


def main():
    RUNS.mkdir(parents=True, exist_ok=True)
    ref = json.loads((ROOT / "data" / "explorador" / "cuscuton_c1" / "result.json").read_text())["bounce"]
    b = EC.run(EC.c1_row())["bounce"]
    g0 = {k: abs(b[k2] / ref[k] - 1) for k, k2 in (("a", "a"), ("chi", "chi"), ("chi_dot", "chidot"))}
    G0 = bool(max(g0.values()) <= 1e-8)
    A11 = produce(1e-11); A12 = produce(1e-12)
    key = lambda r: (r["set"], r["par"]["AH"], r["par"]["f"])
    R12 = {key(r): r for r in A12}
    valid = lambda r: r.get("constraint_max") is not None and np.isfinite(r["constraint_max"]) and r["constraint_max"] <= 1e-6
    kids = [r for r in A11 if r["outcome"] == "HIJO"]
    edge = lambda r: r["minWrel_after_bounce"] < EDGE
    G1 = all(valid(r) for r in kids)
    bad = [dict(key=key(r), N11=r["N_inf"], out12=R12[key(r)]["outcome"], N12=R12[key(r)].get("N_inf")) for r in kids
           if not edge(r) and not (R12[key(r)]["outcome"] == "HIJO" and abs(R12[key(r)]["N_inf"] - r["N_inf"]) <= 0.1)]
    G2 = not bad
    G3 = not any(r.get("D_negative") for r in A11 + A12)
    pct = {s: sum(valid(r) for r in A11 if r["set"] == s) / sum(1 for r in A11 if r["set"] == s) for s in U57.SETS}
    table = {}
    for s in U57.SETS:
        table[s] = {}
        for A in U57.AHS:
            rs = [r for r in A11 if r["set"] == s and r["par"]["AH"] == A]
            k_ = [r for r in rs if r["outcome"] == "HIJO" and valid(r) and not edge(r)]
            e_ = [r for r in rs if r["outcome"] == "HIJO" and edge(r)]
            span = math.log10(max(r["par"]["f"] for r in k_) / min(r["par"]["f"] for r in k_)) if len(k_) >= 2 else 0.0
            table[s][f"{A:g}"] = dict(classes=dict(Counter(r["outcome"] for r in rs)), n_kids=len(k_), n_edge=len(e_), f_span_dex=span,
                                      interval=bool(len(k_) >= 2 and span >= 0.1),
                                      chi_land=[min(r["chi_land"] for r in k_), max(r["chi_land"] for r in k_)] if k_ else None,
                                      N=[min(r["N_inf"] for r in k_), max(r["N_inf"] for r in k_)] if k_ else None)
    fert = {s: any(v["interval"] for v in table[s].values()) for s in U57.SETS}
    if not (G0 and G1 and G2 and G3):
        verdict = "INCONCLUSO"
    elif fert["C1"] and fert["eps_pequeno"]:
        verdict = "CUSCUTON_FERTIL"
    elif fert["C1"]:
        verdict = "SOLO_C1"
    elif fert["eps_pequeno"]:
        verdict = "SOLO_EPS_PEQUENO"
    else:
        verdict = "ESTERIL" if all(p >= 0.9 for p in pct.values()) else "INCONCLUSO"
    out = dict(spec="docs/preregistros/ESPEC_U57b.md", gates=dict(G0=G0, G1=G1, G2=G2, G3=G3), G0_rel=g0, G2_changed=bad,
               valid_fraction=pct, n_kids=len(kids), n_edge=sum(edge(r) for r in kids),
               AH_min_with_kids={s: min([float(A) for A, v in table[s].items() if v["n_kids"]], default=None) for s in U57.SETS},
               verdict=verdict, table=table)
    save("u57b_cuscuton_fondo.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_rel", "G2_changed", "valid_fraction", "n_kids", "n_edge", "AH_min_with_kids", "verdict")}, indent=1, default=float))
    for s in U57.SETS:
        for A, v in table[s].items():
            print(s, A, v)


if __name__ == "__main__":
    main()
