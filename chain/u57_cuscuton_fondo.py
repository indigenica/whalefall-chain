"""U57 (docs/preregistros/ESPEC_U57.md): the cuscuton birth, homogeneous background.  Imports chain/explora_cuscuton.py
unchanged (its `run`); adds the preregistered grid (two cuscuton parameter sets), the gates and the verdict.

    python3 -B chain/u57_cuscuton_fondo.py   ->  results/u57_runs/rows_rtol{11,12}.json, results/u57_cuscuton_fondo.json

Internal, not exported."""
import json, math, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
from ucommon import ROOT, save
import explora_cuscuton as EC

RUNS = ROOT / "results" / "u57_runs"
SETS = {"C1": dict(mu=0.3, m=0.05, pinf=0.25), "eps_pequeno": dict(mu=0.16, m=0.02, pinf=0.25)}
AHS = (3.0, 10.0, 20.0, 30.0, 48.0, 100.0, 300.0, 1000.0)
FS = tuple(10 ** (k / 16) for k in range(-64, 1))


def jobs(rtol):
    return [(name, dict(p, AH=A, f=f, rtol=rtol)) for name, p in SETS.items() for A in AHS for f in FS]


def produce(rtol):
    path = RUNS / f"rows_rtol{-int(round(math.log10(rtol)))}.json"
    if path.exists():
        return json.loads(path.read_text())
    t0 = time.time()
    with Pool(4) as pool:
        rows = pool.map(EC.job, jobs(rtol), chunksize=1)
    path.write_text(json.dumps(dict(rtol=rtol, rows=rows, seconds=time.time() - t0), default=float))
    print(f"rtol {rtol:g}: {len(rows)} rows in {time.time() - t0:.0f}s", flush=True)
    return json.loads(path.read_text())


def main():
    RUNS.mkdir(parents=True, exist_ok=True)
    # G0
    ref = json.loads((ROOT / "data" / "explorador" / "cuscuton_c1" / "result.json").read_text())["bounce"]
    c1 = EC.run(EC.c1_row()); b = c1["bounce"]
    g0 = {k: abs(b[k2] / ref[k] - 1) for k, k2 in (("a", "a"), ("chi", "chi"), ("chi_dot", "chidot"))}
    G0 = bool(max(g0.values()) <= 1e-8)
    A11 = produce(1e-11)["rows"]; A12 = produce(1e-12)["rows"]
    key = lambda r: (r["set"], r["par"]["AH"], r["par"]["f"])
    R12 = {key(r): r for r in A12}
    valid = lambda r: r.get("constraint_max") is not None and np.isfinite(r["constraint_max"]) and r["constraint_max"] <= 1e-6
    # G1
    g1 = {}
    for s in SETS:
        rs = [r for r in A11 if r["set"] == s]
        g1[s] = dict(n=len(rs), n_valid=sum(valid(r) for r in rs), kids_invalid=sum(r["outcome"] == "HIJO" and not valid(r) for r in rs))
    G1 = all(v["kids_invalid"] == 0 and v["n_valid"] >= 0.9 * v["n"] for v in g1.values())
    # G2
    kid_bad = []; same = 0; nval = 0
    for r in A11:
        q = R12[key(r)]
        if r["outcome"] == "HIJO" and not (q["outcome"] == "HIJO" and abs(q["N_inf"] - r["N_inf"]) <= 0.1):
            kid_bad.append(dict(key=key(r), N11=r.get("N_inf"), out12=q["outcome"], N12=q.get("N_inf")))
        if valid(r):
            nval += 1; same += (q["outcome"] == r["outcome"])
    G2 = bool(not kid_bad and same >= 0.95 * nval)
    # G3
    G3 = not any(r.get("D_negative") for r in A11 + A12)
    # per set and A.H
    table = {}
    for s in SETS:
        table[s] = {}
        for A in AHS:
            rs = [r for r in A11 if r["set"] == s and r["par"]["AH"] == A and valid(r)]
            kids = [r for r in rs if r["outcome"] == "HIJO"]
            span = math.log10(max(r["par"]["f"] for r in kids) / min(r["par"]["f"] for r in kids)) if len(kids) >= 2 else 0.0
            table[s][f"{A:g}"] = dict(
                classes=dict(Counter(r["outcome"] for r in rs)), n_invalid=sum(1 for r in A11 if r["set"] == s and r["par"]["AH"] == A and not valid(r)),
                n_kids=len(kids), f_span_dex=span, interval=bool(len(kids) >= 2 and span >= 0.1),
                f_kids=[min(r["par"]["f"] for r in kids), max(r["par"]["f"] for r in kids)] if kids else None,
                chi_land=[min(r["chi_land"] for r in kids), max(r["chi_land"] for r in kids)] if kids else None,
                N=[min(r["N_inf"] for r in kids), max(r["N_inf"] for r in kids)] if kids else None,
                rho_bounce=max((r["bounce"]["rho_matter"] for r in rs if "bounce" in r), default=None),
                minus_U_over_3H2_end=[min(r["minus_U_over_3H2_end"] for r in kids), max(r["minus_U_over_3H2_end"] for r in kids)] if kids else None,
                constraint_max_kids=max((r["constraint_max"] for r in kids), default=None))
    fert = {s: any(v["interval"] for v in table[s].values()) for s in SETS}
    if not (G0 and G1 and G2 and G3):
        verdict = "INCONCLUSO"
    elif fert["C1"] and fert["eps_pequeno"]:
        verdict = "CUSCUTON_FERTIL"
    elif fert["C1"]:
        verdict = "SOLO_C1"
    elif fert["eps_pequeno"]:
        verdict = "SOLO_EPS_PEQUENO"
    else:
        verdict = "ESTERIL"
    out = dict(spec="docs/preregistros/ESPEC_U57.md", gates=dict(G0=G0, G1=G1, G2=G2, G3=G3), G0_rel=g0, G1=g1,
               G2_kids_changed=kid_bad, G2_same_class=[same, nval], eps_G={s: EC.Cusc(**p).eps_G() for s, p in SETS.items()},
               AH_min_with_kids={s: min([float(A) for A, v in table[s].items() if v["n_kids"]], default=None) for s in SETS},
               verdict=verdict, table=table)
    save("u57_cuscuton_fondo.json", out)
    print(json.dumps({k: out[k] for k in ("gates", "G0_rel", "G1", "G2_same_class", "eps_G", "AH_min_with_kids", "verdict")}, indent=1, default=float))
    print("G2 kids changed:", kid_bad)
    for s in SETS:
        for A, v in table[s].items():
            print(s, A, v["classes"], "inv", v["n_invalid"], "kids", v["n_kids"], "span %.2f" % v["f_span_dex"], v["f_kids"], v["chi_land"], v["N"])


if __name__ == "__main__":
    main()
