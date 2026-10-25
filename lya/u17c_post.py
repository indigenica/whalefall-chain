"""U17c, post-hoc reading (no verdict; docs/informes/INFORME_U17c.md §3-§4).  Reads lya/results/u17c/U17C_RESULT.json and
writes lya/results/u17c/U17C_POST.json:
  - A*eta = K_eff*c/eta*eta = (8 pi/27) xi*ell*eta, the variable that sets the free-streaming of the collapse population
    when K_eff >> 1 (m_G ~ K^-2, so lambda_fs ~ c/sqrt(m_G) ~ c*K ~ A*eta, independent of c);
  - the value of A*eta where deltaA crosses the strict threshold (linear interpolation in ln(A*eta), sigma_ln = 0.5 rows);
  - an upper estimate of Delta N_eff from the (1-eta) part of the network energy, if it is dark radiation:
        rho_DR/rho_coll = (1-eta) c/eta * a_osc/a,   a_osc/a_eq <= T_eq/T_osc  (g*s(T_osc) >= g*s(T_eq)),
        Delta N_eff = rho_DR/rho_nu1 at a_eq = (1-eta) c/eta * (T_eq/T_osc) * (omega_coll/omega_m) * 7.45;
  - the hard-quantum momentum needed for that component to stay below 1 % of the dark matter once it is non-relativistic:
        f_h = (1-eta) c/eta * (m_G/p_h) * f_nc  <  0.01   =>   p_h/m_G > 100 (1-eta) c/eta * f_nc.
    python -B lya/u17c_post.py"""
import json, math
from pathlib import Path

D = Path(__file__).resolve().parent / "results" / "u17c"
T_EQ_EV, OMEGA_M, RHO_R_OVER_NU1 = 0.80, 0.1424, (1 + 0.2271 * 3.044) / 0.2271


def main():
    R = json.loads((D / "U17C_RESULT.json").read_text()); strict = R["deltaA_strict"]
    rows = []
    for r in R["rows"]:
        A_eta = r["K_eff"] * r["c"]                                   # = (8 pi/27) xi ell eta
        ratio = (1 - r["eta"]) * r["c"] / r["eta"]
        dneff = ratio * T_EQ_EV / (r["T_osc_MeV"] * 1e6) * r["Omega_collapse"] / OMEGA_M * RHO_R_OVER_NU1
        rows.append(dict(tag=r["tag"], branch=r["branch"], c=r["c"], eta=r["eta"], sigma_ln=r["sigma_ln"], A_eta=A_eta,
                         deltaA=r["deltaA"], LYA_primary=r["LYA_primary"], dNeff_upper=dneff,
                         p_hard_over_mG_min=(100 * ratio * r["f_nc"] if r["eta"] < 1 else None)))
    s5 = sorted([x for x in rows if x["sigma_ln"] == 0.5], key=lambda x: x["A_eta"])
    cross = None
    for lo, hi in zip(s5, s5[1:]):
        if lo["deltaA"] <= strict < hi["deltaA"]:
            t = (strict - lo["deltaA"]) / (hi["deltaA"] - lo["deltaA"])
            cross = math.exp(math.log(lo["A_eta"]) + t * (math.log(hi["A_eta"]) - math.log(lo["A_eta"])))
    out = dict(source="lya/results/u17c/U17C_RESULT.json", post_hoc=True, deltaA_strict=strict, A_eta_at_strict=cross,
               xi_ell_eta_at_strict=(cross / (8 * math.pi / 27) if cross else None), rows=rows)
    (D / "U17C_POST.json").write_text(json.dumps(out, indent=1, allow_nan=False))
    for x in sorted(rows, key=lambda x: x["A_eta"]):
        print(f"{x['tag']:<28} A*eta={x['A_eta']:8.1f} dA={x['deltaA']:.3f} primary={x['LYA_primary']!s:<5} "
              f"dNeff<={x['dNeff_upper']:.1e} p_h/m_G>{x['p_hard_over_mG_min'] or 0:.0f}")
    print(f"A*eta at deltaA = {strict:.4f}: {cross:.1f}  (xi*ell*eta = {out['xi_ell_eta_at_strict']:.1f})")


if __name__ == "__main__":
    main()
