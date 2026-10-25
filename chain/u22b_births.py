"""U22b (docs/preregistros/ESPEC_U22b.md): the whole C45h grid (frac_rad = 0) with the sum of 10 levels (mode 0, must reproduce) and
the order-4 pole (mode 4, C = 0.3516): fertile fraction, chi_et and f_ES (fertile in sep 1e-4..1e-6 with chi_land > chi_et, E < 1).
    python -B chain/u22b_births.py      ->  results/u22b_births.json"""
import json, math, sys, time
import numpy as np
from ucommon import ROOT, benchmark, load, save
sys.path.insert(0, str(ROOT / "components" / "c45b_hom"))
import u22_pole as U22

WINDOW = (1e-4, 1e-5, 1e-6)


def grid(mode, C, chi_et, chi_end):
    from frw_c import lib, c45_state
    B = benchmark(); H = load("c45b/c45h_results.json"); rows = [r for r in H["rows"] if r["frac_rad"] == 0.0]
    lib().set_deform(mode, C if mode == 4 else 1e30, 0.0); out = []
    try:
        for r in rows:
            st = c45_state(B["A_C45"], r["lam"], r["sigma"], 1, r["sep"]); o = np.zeros(12)
            code = lib().run(st["a0"], st["al"], st["V0"], st["chi0"], st["chid0"], 60 * 4.0 * st["a0"], 60.0, 1e-11, 400_000_000, chi_end, o, 0.0)
            c2 = st["V0"] / B["V0_BENCH_C45"]; land = abs(o[2]); fert = bool(code == 1 and o[1] > 5 and o[4] < 1e-3)
            out.append(dict(sep=r["sep"], lam=r["lam"], sigma=r["sigma"], fertile=fert, chi_land=land, E_phys_max=o[10] / c2, n_b=int(o[0]),
                            eternal=bool(fert and land > chi_et and o[10] / c2 < 1), rel_diff_stored=abs(land - r["chi_land"]) / max(r["chi_land"], 1e-300),
                            fertile_stored=r["fertile"]))
    finally:
        lib().set_deform(0, 1e30, 0.0)
    return out


def summary(rows):
    w = [r for r in rows if r["fertile"] and r["sep"] in WINDOW]
    by = {f"{s:g}": dict(n_fertile=sum(r["fertile"] for r in rows if r["sep"] == s),
                         n_eternal=sum(r["eternal"] for r in rows if r["sep"] == s)) for s in sorted({r["sep"] for r in rows})}
    return dict(n=len(rows), n_fertile=sum(r["fertile"] for r in rows), f_ES=(sum(r["eternal"] for r in w) / len(w)) if w else 0.0,
                n_window_fertile=len(w), by_sep=by, median_chi_land_fertile=float(np.median([r["chi_land"] for r in rows if r["fertile"]])))


def main():
    t0 = time.time(); B = benchmark("sum")                     # the sum's reference V0 (U26: the chain default is now the pole)
    from frw import CHI_END
    from s1_birth import chi_eternal
    et_sum = chi_eternal(B["V0"])
    nr = U22.s2_pole(U22.C_FIT); et_pole = U22.chi_eternal_pole(nr["V0"], U22.C_FIT)
    P = U22.PolePotential(1., U22.C_FIT)
    from scipy.optimize import brentq
    xend = brentq(lambda x: 0.5 * (float(P.dV(x)) / float(P.V(x))) ** 2 - 1, 1e-3, 50)
    s = grid(0, None, et_sum, CHI_END); p = grid(4, U22.C_FIT, et_pole, xend)
    maxdiff = max(r["rel_diff_stored"] for r in s if r["fertile_stored"])
    S, Pp = summary(s), summary(p)
    R = Pp["f_ES"] / S["f_ES"] if S["f_ES"] > 0 else float("nan")
    v = ("INVALIDO" if maxdiff > 1e-8 else "NACIMIENTO_ROBUSTO" if R >= 0.5 else "NACIMIENTO_REDUCIDO" if R >= 0.1 else "NACIMIENTO_FRAGIL")
    out = dict(spec="docs/preregistros/ESPEC_U22b.md", chi_et_sum=et_sum, chi_et_pole=et_pole, mode0_max_rel_diff=maxdiff, sum=S, pole=Pp,
               R=R, verdict_U22b=v, runtime_s=time.time() - t0)
    save("u22b_births.json", out)
    print("sum", S); print("pole", Pp); print("R", R, "U22b:", v, f"({out['runtime_s']:.0f}s)")


if __name__ == "__main__":
    main()
