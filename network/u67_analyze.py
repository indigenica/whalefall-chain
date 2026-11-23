"""U67 analysis (docs/preregistros/ESPEC_U67.md): gates G2-G3 and the verdict OLVIDA / RECUERDA / INTERMEDIO / INCONCLUSO.
G0 and G1 come from network/results/u67/U67_CHECKS.json (u67_produce.py --checks).

    python -B network/u67_analyze.py   ->  network/results/u67/U67_RESULT.json"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "u67"
ELL_C, NBOOT = 4.10, 10000
PROTS = ("RD_SPARSE", "RD_DENSE", "MD2", "MD3")


def xi_at(rows, ell):
    e = np.array([r["ell"] for r in rows]); x = np.array([r["xi"] for r in rows])
    if ell < e.min() or ell > e.max(): return float("nan")
    return float(np.interp(ell, e, x))


def main(u67b=False):
    chk = json.loads((OUT / "U67_CHECKS.json").read_text())
    if u67b:        # ESPEC_U67b: MDA2/MDA3 decide, MDO2 informative; RD controls from U67
        out_dir = HERE / "results" / "u67b"
        runs = {p: [json.loads(f.read_text()) for f in sorted((OUT / "runs").glob(f"RP3_N192_{p}_s*.json"))] for p in ("RD_SPARSE", "RD_DENSE")}
        for p in ("MDA2", "MDA3", "MDO2"):
            runs[p] = [json.loads(f.read_text()) for f in sorted((out_dir / "runs").glob(f"RP3_N192_{p}_s*.json"))]
        return analyse(runs, chk, ("MDA2", "MDA3"), ("RD_DENSE", "MDO2"), out_dir / "U67B_RESULT.json", "docs/preregistros/ESPEC_U67b.md",
                       g2_mode="vacuum")
    runs = {p: [json.loads(f.read_text()) for f in sorted((OUT / "runs").glob(f"RP3_N192_{p}_s*.json"))] for p in PROTS}
    return analyse(runs, chk, ("MD2", "MD3"), ("RD_DENSE",), OUT / "U67_RESULT.json", "docs/preregistros/ESPEC_U67.md", g2_mode="ordered")


def analyse(runs, chk, decide, extra, path, spec, g2_mode):
    PROTS_ = tuple(runs)
    rng = np.random.default_rng(67)
    G2_rows, G2 = [], True
    for p in decide:
        for r in runs[p]:
            done = next((x["ell"] for x in r["rows"] if x.get("gap", 0) >= 0.5), None)
            if g2_mode == "ordered":         # U67: ordered (gap >= 0.5) at the first snapshot with ell >= 3.5
                first = next((x for x in r["rows"] if x["ell"] >= 3.5), None); g = first["gap"] if first else float("nan")
            else:                             # U67b: still in the vacuum (gap >= 0.5) at every snapshot with ell >= 2
                g = min(x["gap"] for x in r["rows"] if x["ell"] >= 2.0)
            G2_rows.append(dict(protocol=p, seed=r["seed"], gap=g, ell_ordered=done)); G2 = G2 and bool(g >= 0.5)
    for p in extra:
        for r in runs[p]:
            if "gap" in r["rows"][0]:
                G2_rows.append(dict(protocol=p, seed=r["seed"], informative=True,
                                    ell_ordered=next((x["ell"] for x in r["rows"] if x.get("gap", 0) >= 0.5), None)))
    gated = tuple(decide) + ("RD_SPARSE", "RD_DENSE")
    plaq_end = {f"{p}_s{r['seed']}": r["rows"][-1]["n_plaq"] for p in PROTS_ for r in runs[p]}
    G3 = all(plaq_end[f"{p}_s{r['seed']}"] >= 200 for p in gated for r in runs[p]) and all(len(runs[p]) == 3 for p in PROTS_)
    xi_c = {p: [xi_at(r["rows"], ELL_C) for r in runs[p]] for p in PROTS_}
    mean = {p: float(np.mean(v)) for p, v in xi_c.items()}

    def boot_ratio(p):
        a = np.array(xi_c[p]); b = np.array(xi_c["RD_SPARSE"])
        ia = rng.integers(0, len(a), (NBOOT, len(a))); ib = rng.integers(0, len(b), (NBOOT, len(b)))
        return float(np.std(a[ia].mean(1) / b[ib].mean(1) - 1))
    others = tuple(decide) + tuple(extra)
    D = {p: mean[p] / mean["RD_SPARSE"] - 1 for p in others}
    sig = {p: boot_ratio(p) for p in others}
    gates = dict(G0=chk["G0"], G1=chk["G1"], G2=G2, G3=G3)
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif all(abs(D[p]) <= 0.15 for p in decide): verdict = "OLVIDA"
    elif any(abs(D[p]) > 0.30 and abs(D[p]) > 3 * sig[p] for p in decide): verdict = "RECUERDA"
    else: verdict = "INTERMEDIO"
    # informative: D(ell) along the window
    ells = [2.5, 3.0, 3.5, 3.8, 4.0, ELL_C]
    D_ell = {p: {f"{e:g}": float(np.nanmean([xi_at(r["rows"], e) for r in runs[p]]) /
                             np.nanmean([xi_at(r["rows"], e) for r in runs["RD_SPARSE"]]) - 1) for e in ells}
             for p in others}
    xi_ell = {p: {f"{e:g}": [float(np.nanmean([xi_at(r["rows"], e) for r in runs[p]])),
                             float(np.nanstd([xi_at(r["rows"], e) for r in runs[p]]))] for e in ells} for p in PROTS_}
    km = {p: float(np.mean([r["rows"][-1]["KM_all"] / max(r["rows"][-1]["KG_all"], 1e-300) for r in runs[p]])) for p in PROTS_}
    res = dict(spec=spec, ell_c=ELL_C, gates=gates, checks=chk, G2_rows=G2_rows, plaq_end=plaq_end,
               xi_at_ell_c=xi_c, mean_xi=mean, D=D, sigma_D=sig, verdict=verdict, D_of_ell=D_ell, xi_of_ell=xi_ell,
               KM_over_KG_end=km, runtime_s={p: [r["runtime_s"] for r in runs[p]] for p in PROTS_})
    path.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(dict(gates=gates, D=D, sigma=sig, mean_xi=mean, verdict=verdict, G2=G2_rows, D_of_ell=D_ell), indent=1, default=float))


if __name__ == "__main__":
    import sys
    main(u67b="--u67b" in sys.argv[1:])
