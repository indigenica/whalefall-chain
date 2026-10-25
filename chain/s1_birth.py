"""S1 · birth (ballena 5): contracting closed seed -> NPQT bounce(s) -> chi pump + throat shot -> landing chi_land.

Re-runs, with the vendored C port of the C45 equations (components/c45b_hom/frw_core.c, unchanged), every C45j-viable child
(eternal, sub-Planckian and tolerant to zeta0 = 4.6e-5 at L.H = 1) plus a deterministic control sample of C45h rows, and
checks bit-level reproduction.  Then re-classifies 'eternal' with the plateau height V0 of the pivot matching (S2), which
differs by 2.8 % from the V0_BENCH of the C45 scans.

U26 (docs/preregistros/ESPEC_U26.md): with the pole plateau (the default) the reproduction above stays as the control of the C port
(mode 0, the sum, with the sum's V0), and the primary row is the U22b statistic: the whole C45h grid (frac_rad = 0) in mode 4 (pole)
and mode 0, f_ES for each with its own chi_et, and R = f_ES(pole)/f_ES(sum) >= 0.5.  With WHALEFALL_PLATEAU=sum, S1 is as before."""
import math, shutil, time
from ucommon import load, save, benchmark, provenance, PLATEAU
import numpy as np
from scipy.optimize import brentq

ZETA0 = 4.6e-5


def logPz(x, V0):
    from frw import Vn, dVn
    return 3 * math.log(V0 * Vn(x)) - math.log(12 * math.pi ** 2) - 2 * math.log(max(V0 * dVn(x), 1e-300))


def chi_eternal(V0):
    return brentq(lambda x: logPz(x, V0), 5.0, 1e5, xtol=1e-6)


def main():
    B = benchmark("sum")                      # the reproduction runs are mode 0 (the sum), with the sum's common V0
    H = load("c45b/c45h_results.json"); J = load("c45b/c45j_results.json")
    rows = [r for r in H["rows"] if r["frac_rad"] == 0.0]
    viable = [c for c in J["children"] if c["zeta0_max_LH1"] >= ZETA0]
    out = dict(stage="S1_birth", n_viable_C45j=len(viable), provenance=provenance())
    chi_et_bench, chi_et_common = chi_eternal(B["V0_BENCH_C45"]), chi_eternal(B["V0"])
    out.update(chi_et_V0_BENCH=chi_et_bench, chi_et_V0_common=chi_et_common, chi_et_C45h_stored=H["chi_et"])
    if shutil.which("gcc") is None:
        out.update(status="SKIPPED_NO_GCC", reproduction=None)
        save("s1_birth.json", out); print("S1 SKIPPED (no gcc)"); return out
    from frw_c import run_c, c45_state
    t0 = time.time()
    # deterministic control: every 13th fertile row of C45h with frac_rad = 0 (covers all separations)
    ctrl = [r for r in rows if r["fertile"]][::13]
    targets = [("viable", next(r for r in rows if r["sep"] == c["sep"] and r["lam"] == c["lam"] and r["sigma"] == c["sigma"]), c)
               for c in viable] + [("control", r, None) for r in ctrl]
    rep = []
    for kind, r, c in targets:
        st = c45_state(B["A_C45"], r["lam"], r["sigma"], 1, r["sep"])
        x = run_c(**st)
        c2 = st["V0"] / B["V0_BENCH_C45"]
        rep.append(dict(kind=kind, sep=r["sep"], lam=r["lam"], sigma=r["sigma"], n_b=x["n_bounces"], n_b_stored=r["n_bounces"],
                        chi_land=abs(x["chi_end"]), chi_land_stored=r["chi_land"],
                        rel_diff=abs(abs(x["chi_end"]) - r["chi_land"]) / max(r["chi_land"], 1e-300),
                        fertile=x["fertile"], E_phys_max=x["E_max"] / c2, alpha_phys=r["sep"] * B["A_C45"] ** 2 * c2,
                        AH=(c or {}).get("AH"), zeta0_max_LH1=(c or {}).get("zeta0_max_LH1")))
    maxdiff = max(p["rel_diff"] for p in rep); same_nb = all(p["n_b"] == p["n_b_stored"] for p in rep)
    v = [p for p in rep if p["kind"] == "viable"]
    for p in v:
        p["eternal_common_V0"] = bool(p["chi_land"] > chi_et_common and p["E_phys_max"] < 1)
    fertile_window = [r for r in rows if r["fertile"] and r["sep"] in (1e-4, 1e-5, 1e-6)]
    es_common = [r for r in fertile_window if r["chi_land"] > chi_et_common and r["E_phys_max"] < 1]
    out.update(status="PASS" if (maxdiff < 1e-8 and same_nb) else "FAIL", max_rel_diff_chi_land=maxdiff, same_bounce_counts=same_nb,
               runtime_s=time.time() - t0, reproduction=rep,
               viable_all_eternal_with_common_V0=all(p["eternal_common_V0"] for p in v),
               f_ES_window_common_V0=len(es_common) / len(fertile_window),
               f_viable_C45j=J["frac_viable_LH1_over_fertile"],
               viable_ranges=dict(AH=[min(p["AH"] for p in v), max(p["AH"] for p in v)],
                                  chi_land=[min(p["chi_land"] for p in v), max(p["chi_land"] for p in v)],
                                  n_bounces=sorted(p["n_b"] for p in v),
                                  alpha_phys=[min(p["alpha_phys"] for p in v), max(p["alpha_phys"] for p in v)],
                                  ell_QT_over_lP=[math.sqrt(min(p["alpha_phys"] for p in v)), math.sqrt(max(p["alpha_phys"] for p in v))]))
    if PLATEAU == "pole":
        out["control_status_sum_mode0"] = out["status"]; out.update(pole_births())
        out["status"] = "PASS" if (out["control_status_sum_mode0"] == "PASS" and out["grid_mode0_max_rel_diff"] <= 1e-8
                                   and out["R_f_ES_pole_over_sum"] >= 0.5) else "FAIL"          # mode-0 grid gate: ESPEC_U26 enmienda 2
        print(f"S1 (pole) f_ES {out['f_ES_pole']:.3f} vs sum {out['f_ES_sum']:.3f}: R = {out['R_f_ES_pole_over_sum']:.3f}; "
              f"viable C45j children eternal with the pole: {out['viable_eternal_pole']}/{out['n_viable_C45j']}")
    save("s1_birth.json", out)
    print(f"S1 {out['status']}: {len(rep)} runs, max rel diff {maxdiff:.1e}, chi_et {chi_et_bench:.2f} -> {chi_et_common:.2f}, "
          f"viable eternal with common V0: {out['viable_all_eternal_with_common_V0']}, {out['runtime_s']:.0f}s")
    return out


def pole_births():
    """U26 row C2: the U22b grid (chain/u22b_births.py), both plateaus, each with its own chi_et"""
    import u22_pole as U22, u22b_births as U22B
    from frw import CHI_END
    from scipy.optimize import brentq as _brentq
    Bp = benchmark("pole"); et_sum = chi_eternal(benchmark("sum")["V0"])
    P = U22.PolePotential(Bp["V0"], U22.C_FIT); et_pole = P.eternal()
    xend = _brentq(lambda x: 0.5 * (float(P.dV(x)) / float(P.V(x))) ** 2 - 1, 1e-3, 50)
    t0 = time.time(); g_sum = U22B.grid(0, None, et_sum, CHI_END); g_pole = U22B.grid(4, U22.C_FIT, et_pole, xend)
    S, Pp = U22B.summary(g_sum), U22B.summary(g_pole)
    J = load("c45b/c45j_results.json"); keys = {(c["sep"], c["lam"], c["sigma"]) for c in J["children"] if c["zeta0_max_LH1"] >= ZETA0}
    return dict(plateau_row="C2 of ESPEC_U26: U22b statistic", chi_et_sum=et_sum, chi_et_pole=et_pole, chi_end_pole=xend,
                f_ES_sum=S["f_ES"], f_ES_pole=Pp["f_ES"], R_f_ES_pole_over_sum=Pp["f_ES"] / S["f_ES"], grid_sum=S, grid_pole=Pp,
                grid_mode0_max_rel_diff=max(r["rel_diff_stored"] for r in g_sum if r["fertile_stored"]),
                viable_eternal_pole=sum(r["eternal"] for r in g_pole if (r["sep"], r["lam"], r["sigma"]) in keys),
                note_C45j="the C45j perturbation filter (zeta0, f_viable) was computed with the sum and is not redone", grid_runtime_s=time.time() - t0)


if __name__ == "__main__":
    main()
