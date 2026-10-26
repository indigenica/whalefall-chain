"""U27 (docs/preregistros/ESPEC_U27.md, amendments 1-2): multi-bounce coherence of the direct (non-eternal) children, in the
C45j model, with the pole plateau.  V1: the C45j code (components/c45j_ltb, unchanged) reproduces data/c45b/c45j_results.json;
V2: the pole grid (same integration as chain/u22b_births.grid, plus chi_abs_max) reproduces the stored U22b pole summary.
Rows R1-R5 and the verdicts of ESPEC_U27 sections 4-5.  eps_est runs in a subprocess inside components/c45j_ltb, whose
package `common` would clash with components/cadena_cierre/common.py.
    python -B chain/u27_direct_children.py      ->  results/u27_direct_children.json (~3 min)"""
import glob, json, math, subprocess, sys, time
import numpy as np
from ucommon import ROOT, benchmark, load, save
import u22_pole as U22
import u22b_births as U22B

C45J = ROOT / "components" / "c45j_ltb"
SEPS = (1e-4, 1e-5, 1e-6)
LHS = (0.3, 1.0, 3.0)
ZETA0 = 4.6e-5
BAND = (10.0, 50.0)
CHI_ET_SUM = 417.28352731287106                  # C45h / C45j (amendment 2)


def v1():
    """run C45j as is and compare with the archived result (amendment 1: full-precision values and the same integers)"""
    subprocess.run([sys.executable, "-B", "-m", "ltb.multibounce_viability"], cwd=C45J, check=True, capture_output=True)
    run = sorted(glob.glob(str(C45J / "results" / "c45j_multibounce_*")))[-1]
    new = json.loads(open(run + "/results.json").read()); old = load("c45b/c45j_results.json")
    k = lambda c: (c["sep"], c["lam"], c["sigma"])
    N, O = {k(c): c for c in new["children"]}, {k(c): c for c in old["children"]}
    rel = max(abs(N[key][f] - o[f]) / abs(o[f]) for key, o in O.items() for f in o if f.startswith(("eps_est", "zeta0")) and key in N) if set(N) == set(O) else math.inf
    nv = sum(c["zeta0_max_LH1"] >= ZETA0 for c in new["children"]); nf = sum(v["n_fertile"] for v in new["summary"].values())
    ok = bool(set(N) == set(O) and len(N) == 91 and rel <= 1e-6 and nv == 20 and nf == 557
              and abs(new["frac_viable_LH1_over_fertile"] - old["frac_viable_LH1_over_fertile"]) <= 1e-9)
    return dict(passed=ok, run=str(Path_rel(run)), n_children=len(N), max_rel_diff=rel, n_viable=nv, n_fertile=nf,
                frac=new["frac_viable_LH1_over_fertile"], frac_archived=old["frac_viable_LH1_over_fertile"])


def Path_rel(p):
    from pathlib import Path
    return Path(p).relative_to(ROOT)


def pole_grid():
    """chain/u22b_births.grid in mode 4, also recording chi_abs_max (out[1])"""
    from frw_c import lib, c45_state
    B = benchmark(); H = load("c45b/c45h_results.json"); rows = [r for r in H["rows"] if r["frac_rad"] == 0.0]
    nr = U22.s2_pole(U22.C_FIT); chi_et = U22.chi_eternal_pole(nr["V0"], U22.C_FIT)
    P = U22.PolePotential(1., U22.C_FIT)
    from scipy.optimize import brentq
    xend = brentq(lambda x: 0.5 * (float(P.dV(x)) / float(P.V(x))) ** 2 - 1, 1e-3, 50)
    lib().set_deform(4, U22.C_FIT, 0.0); out = []
    try:
        for r in rows:
            st = c45_state(B["A_C45"], r["lam"], r["sigma"], 1, r["sep"]); o = np.zeros(12)
            code = lib().run(st["a0"], st["al"], st["V0"], st["chi0"], st["chid0"], 60 * 4.0 * st["a0"], 60.0, 1e-11, 400_000_000, xend, o, 0.0)
            c2 = st["V0"] / B["V0_BENCH_C45"]; land = abs(o[2]); fert = bool(code == 1 and o[1] > 5 and o[4] < 1e-3)
            out.append(dict(sep=r["sep"], lam=r["lam"], sigma=r["sigma"], fertile=fert, chi_land=land, chi_abs_max=float(o[1]),
                            E_phys_max=o[10] / c2, n_b=int(o[0]), eternal=bool(fert and land > chi_et and o[10] / c2 < 1)))
    finally:
        lib().set_deform(0, 1e30, 0.0)
    return out, chi_et


def eps_table(keys):
    """eps_est(AH, alphaH2, LH) from the unchanged C45j code, in a subprocess"""
    code = ("import json, sys\nfrom ltb.multibounce_viability import eps_est\n"
            "ks = json.load(sys.stdin)\nprint(json.dumps([eps_est(a, b, c)[0] for a, b, c in ks]))")
    r = subprocess.run([sys.executable, "-B", "-c", code], cwd=C45J, input=json.dumps(keys), capture_output=True, text=True, check=True)
    return dict(zip(map(tuple, keys), json.loads(r.stdout.strip().splitlines()[-1])))


def parent(lam, sep):                               # as in ltb/multibounce_viability.py
    return 0.385 * math.sqrt(lam), sep * 4 * lam / 27


def rows_eval(children, E):
    out = []
    for c in children:
        AH, aH2 = parent(c["lam"], c["sep"]); rec = dict(c, AH=AH, alphaH2=aH2)
        for LH in LHS:
            e = E[(round(AH, 6), round(aH2, 12), LH)]
            rec[f"zeta0_max_LH{LH:g}"] = e / max(c["n_b"], 1) / ((2 / 3) * math.pi ** 2 * (AH / LH) ** 2)
        out.append(rec)
    return out


def row_summary(children, n_fertile):
    s = dict(n=len(children))
    for LH in LHS:
        v = [c for c in children if c[f"zeta0_max_LH{LH:g}"] >= ZETA0]
        s[f"LH{LH:g}"] = dict(n_viable=len(v), f_v=len(v) / n_fertile, f_cond=(len(v) / len(children)) if children else None,
                              median_zeta0_max=float(np.median([c[f"zeta0_max_LH{LH:g}"] for c in children])) if children else None,
                              by_sep={f"{sp:g}": sum(c["sep"] == sp for c in v) for sp in SEPS})
    s["by_sep"] = {f"{sp:g}": sum(c["sep"] == sp for c in children) for sp in SEPS}
    s["n_b"] = {str(k): sum(c["n_b"] == k for c in children) for k in sorted({c["n_b"] for c in children})}
    s["AH_median"] = float(np.median([c["AH"] for c in children])) if children else None
    return s


def classes(rows, chi_et, eternal_key):
    w = [r for r in rows if r["sep"] in SEPS and r["fertile"]]
    E = [r for r in w if r[eternal_key]]
    D = [r for r in w if r["E_phys_max"] < 1 and r["chi_land"] <= chi_et]
    Db = [r for r in D if BAND[0] <= r["chi_abs_max"] <= BAND[1]]
    return w, dict(E=E, D=D, D10_50=Db)


def main():
    t0 = time.time()
    V1 = v1()
    if not V1["passed"]:
        save("u27_direct_children.json", dict(spec="docs/preregistros/ESPEC_U27.md", V1=V1, verdict="INVALIDO (V1)")); sys.exit("V1 failed")
    pole, et_pole = pole_grid()
    stored = load_result("u22b_births.json")["pole"]; S = U22B.summary(pole)
    v2 = bool(S["n"] == stored["n"] and S["n_fertile"] == stored["n_fertile"] and S["n_window_fertile"] == stored["n_window_fertile"]
              and S["by_sep"] == stored["by_sep"] and abs(S["f_ES"] - stored["f_ES"]) <= 1e-9)
    V2 = dict(passed=v2, chi_et_pole=et_pole, summary=S, stored=stored)
    if not v2:
        save("u27_direct_children.json", dict(spec="docs/preregistros/ESPEC_U27.md", V1=V1, V2=V2, verdict="INVALIDO (V2)")); sys.exit("V2 failed")
    H = [dict(r, n_b=r["n_bounces"]) for r in load("c45b/c45h_results.json")["rows"] if r["frac_rad"] == 0.0]
    wP, CP = classes(pole, et_pole, "eternal"); wS, CS = classes(H, CHI_ET_SUM, "eternal_sub")
    exact = {}                                        # as the C45j cache: keyed by the rounded values, evaluated at the exact ones
    for c in wP + wS:
        a, b = parent(c["lam"], c["sep"])
        for LH in LHS:
            exact.setdefault((round(a, 6), round(b, 12), LH), [a, b, LH])
    keys = sorted(exact)
    vals = eps_table([exact[k] for k in keys]); E = {k: vals[tuple(exact[k])] for k in keys}
    P = {k: row_summary(rows_eval(v, E), len(wP)) for k, v in CP.items()}
    Sm = {k: row_summary(rows_eval(v, E), len(wS)) for k, v in CS.items()}
    fE, fD, fB, nB = P["E"]["LH1"]["f_v"], P["D"]["LH1"]["f_v"], P["D10_50"]["LH1"]["f_v"], P["D10_50"]["LH1"]["n_viable"]
    main_v = ("RAMA_DIRECTA_VIABLE" if fD >= 0.05 else "RAMA_DIRECTA_MARGINAL_COMO_ETERNA" if fD >= fE else
              "RAMA_DIRECTA_MARGINAL_PEOR" if fD >= 0.01 else "RAMA_DIRECTA_NO_COHERENTE")
    band_v = "BANDA_10_50_POBLADA" if (nB >= 5 and fB >= 0.01) else "BANDA_10_50_VACIA"
    out = dict(spec="docs/preregistros/ESPEC_U27.md", V1=V1, V2=dict(passed=v2, chi_et_pole=et_pole), n_window_fertile_pole=len(wP),
               n_window_fertile_sum=len(wS), n_eps_evaluations=len(keys), zeta0=ZETA0, band_chi_abs_max=BAND,
               R1_pole_E=P["E"], R2_pole_D=P["D"], R3_pole_D10_50=P["D10_50"], R5_sum=Sm,
               verdict_U27=main_v, verdict_band=band_v, runtime_s=time.time() - t0)
    save("u27_direct_children.json", out)
    print(f"V1 {V1['passed']} ({V1['n_viable']}/{V1['n_fertile']}), V2 {v2}; fertile window pole {len(wP)}, sum {len(wS)}")
    for lab, d in (("pole", P), ("sum", Sm)):
        for k, s in d.items():
            print(lab, k, "n", s["n"], {LH: (s[f'LH{LH:g}']['n_viable'], round(s[f'LH{LH:g}']['f_v'], 4)) for LH in LHS}, "n_b", s["n_b"], "AH_med", s["AH_median"])
    print("U27:", main_v, "| band:", band_v, f"({out['runtime_s']:.0f}s)")


def load_result(name):
    return json.loads((ROOT / "results" / name).read_text())


if __name__ == "__main__":
    main()
