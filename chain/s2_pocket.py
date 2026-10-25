"""S2 · pocket: the SAME hierarchical plateau gives (i) the eternal phase where chi lands (P_zeta > 1) and (ii) the
observable pocket inflation, its end, reheating by sigma and the pivot matching.

(ii) re-runs the continuous pivot matching of cadena_cierre (components/cadena_cierre/thermal.py, unchanged) and compares
with the archived thermal.json.  (i) computes the eternal threshold, the volume-weighted threshold and the D3 flatness
requirement for the undeformed plateau with the common V0, and carries the C45i deformation verdict.

U26 (docs/preregistros/ESPEC_U26.md), rows C1 and C3: with the pole plateau (the default) the matching is the U22 one (chain/u22_pole.py,
s2_pole) against data/cadena_cierre/thermal_pole.json, the CMB gate |dn_s| <= 0.002 and |dr/r| <= 0.10 against the sum's thermal.json,
and the eternal phase with the pole potential.  The C45i deformation verdict was computed with the sum and is carried, declared."""
import math, time
from ucommon import load, save, benchmark, provenance, PLATEAU, THERMAL
import numpy as np
from scipy.optimize import brentq

PLANCK18 = dict(n_s=0.9649, sig=0.0042)          # Planck 2018 TT,TE,EE+lowE+lensing
ACT_DR6 = dict(n_s=0.9743, sig=0.0034)           # ACT DR6 + Planck (P-ACT), 2025
R_MAX = 0.036                                    # BICEP/Keck 2021, 95 %


def main():
    B = benchmark(); arch = load(THERMAL[PLATEAU])["continuous"]
    out = dict(stage="S2_pocket", plateau=PLATEAU, provenance=provenance())
    t0 = time.time()
    import thermal
    if PLATEAU == "pole":
        import u22_pole as U22
        nr = U22.s2_pole(U22.C_FIT, 2e-8)
    else:
        nr = thermal.normalize_continuous(2e-8)
    keys = ("N_star", "n_s", "r", "V0", "chi_star", "chi_end", "T_rh_GeV", "F_GeV", "H_star_GeV")
    cmp = {k: dict(rerun=nr[k], archived=arch[k], rel=abs(nr[k] / arch[k] - 1)) for k in keys}
    worst = max(v["rel"] for v in cmp.values())
    out.update(pivot_rerun=cmp, pivot_max_rel=worst, pivot_status="PASS" if worst < 1e-5 else "FAIL", pivot_runtime_s=time.time() - t0)
    if PLATEAU == "pole":                     # U26 C1: same CMB as the sum (the U22 criteria)
        sm = load(THERMAL["sum"])["continuous"]; dn, dr = nr["n_s"] - sm["n_s"], nr["r"] / sm["r"] - 1
        out["pole_vs_sum"] = dict(dn_s=dn, dr_rel=dr, cmb_equivalent=bool(abs(dn) <= 0.002 and abs(dr) <= 0.10),
                                  rel_change={k: nr[k] / sm[k] - 1 for k in ("V0", "T_rh_GeV", "F_GeV", "H_star_GeV", "m_chi_GeV", "Gamma_GeV", "chi_end")})
        if not out["pole_vs_sum"]["cmb_equivalent"]: out["pivot_status"] = "FAIL"
    # observables against data (informative)
    out["observables"] = dict(n_s=nr["n_s"], r=nr["r"], pull_Planck18=(nr["n_s"] - PLANCK18["n_s"]) / PLANCK18["sig"],
                              pull_ACT_DR6=(nr["n_s"] - ACT_DR6["n_s"]) / ACT_DR6["sig"], r_below_BK18=bool(nr["r"] < R_MAX))
    # eternal phase on the same plateau, common V0

    from vendor_potential import Potential
    pot = Potential(B["V0"])
    chi_et = pot.eternal()
    chi_vol = brentq(lambda x: math.log(pot.pzeta(x) * 6), 10, 2000 if PLATEAU == "sum" else 1e5)
    flat = math.sqrt(B["V0"] / 3) / (2 * math.pi)
    xs = np.geomspace(chi_et, 1e4, 400); v, dv, _ = pot.values(xs); ratio = dv / v
    out["eternal"] = dict(chi_Pzeta_1=chi_et, chi_volume_Pzeta_1_over_6=chi_vol, flatness_bound_per_MPl=flat,
                          max_ratio_above_threshold=float(ratio.max()), undeformed_plateau_meets_bound=bool(ratio.max() <= flat * 1.0001),
                          ratio_over_bound_at=dict(zip(("chi_et", "1e3", "1e4"), [float(np.interp(c, xs, ratio) / flat) for c in (chi_et, 1e3, 1e4)])))
    di = load("c45b/c45i_deform_summary.json")
    out["D3"] = dict(C45i_verdict=di.get("verdict_after_obs_fix", di.get("verdict")),
                     requirement=(f"|V'/V| < sqrt(V0/3)/(2 pi) ~ 1.2e-6 per M_Pl from chi ~ {chi_et:.0f} to ~1e3 M_Pl" if PLATEAU == "pole" else
                                  "|V'/V| < sqrt(V0/3)/(2 pi) ~ 1.2e-6 per M_Pl from chi ~ 4e2 to ~1e3 M_Pl"),
                     UV_status="(b) incompatible with the swampland distance and dS conjectures; (c) open only if they fail (docs/c45b_ballena5/NOTA_D3_UV.md)")
    save("s2_pocket.json", out)
    print(f"S2 pivot {out['pivot_status']} (max rel {worst:.1e}, {out['pivot_runtime_s']:.0f}s); n_s={nr['n_s']:.5f} r={nr['r']:.5f}; "
          f"chi_et={chi_et:.2f} chi_vol={chi_vol:.2f} flat={flat:.2e} plateau ok={out['eternal']['undeformed_plateau_meets_bound']}")
    return out


if __name__ == "__main__":
    main()
