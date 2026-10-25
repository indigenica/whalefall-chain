"""U14 (docs/preregistros/ESPEC_U14.md): can the whale's miniclusters microlens towards the LMC?

Hard-sphere radius (Fairbairn et al. 2017, arXiv:1707.03310, eq. 43) against the point-lens Einstein radius (eq. 33, x = 1/2,
D_S = 50 kpc); delta distribution F(delta > d0) = (1 + d0/1.023)^-0.462 (their fit to Kolb & Tkachev); point-lens bounds
f_max = 0.08 (EROS-2, Tisserand et al. 2007) and 0.20 (EROS-2 + MACHO, Blaineau et al. 2022).

    python3 lya/u14_microlens.py  ->  lya/results/u14/U14_RESULT.json"""
import json, math
from pathlib import Path
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent / "results" / "u14"
G, C, MSUN, KPC = 6.674e-8, 2.998e10, 1.989e33, 3.0857e21          # cgs
OMEGA_M, D_S, X = 0.32, 50 * KPC, 0.5
MASSES = (0.23, 0.5, 1.0, 2.0, 5.0, 10.0)
BOUNDS = {"estricta (EROS-2, 0.08)": 0.08, "amplia (EROS-2+MACHO, 0.20)": 0.20}
RADIUS_FACTORS = {"esfera dura": 1.0, "radio x3": 3.0, "radio /3": 1 / 3}


def r_mc(delta, M, fac=1.0): return fac * 4e16 / (delta * ((delta + 1) * OMEGA_M) ** (1 / 3)) * M ** (1 / 3)
def R_E(M): return 2 * math.sqrt(G * M * MSUN * X * (1 - X) * D_S) / C
def F_tail(d0): return (1 + d0 / 1.023) ** -0.462
def M0(m_eV): return 2.3e-7 * (m_eV / 1e-10) ** -1.5


def delta_lens(M, fac=1.0):
    return brentq(lambda d: math.log(r_mc(d, M, fac) / R_E(M)), 1e-3, 1e8)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = dict(spec="docs/preregistros/ESPEC_U14.md", M0_range_Msun=dict(U7c=[M0(4.9e-15), M0(0.8e-15)], cold_U11=[M0(10e-15), M0(2e-15)]), rows=[])
    for fn, fac in RADIUS_FACTORS.items():
        for M in MASSES:
            d = delta_lens(M, fac); F = F_tail(d)
            row = dict(radius=fn, M_Msun=M, R_E_AU=R_E(M) / 1.496e13, r_mc_delta1_AU=r_mc(1, M, fac) / 1.496e13, delta_lens=d, F_above=F,
                       f_MC_max={b: fm / F for b, fm in BOUNDS.items()})
            out["rows"].append(row)
            print(f"{fn:<12} M={M:5.2f}  R_E={row['R_E_AU']:6.1f} AU  r_mc(δ=1)={row['r_mc_delta1_AU']:7.0f} AU  δ_lens={d:7.1f}  "
                  f"F(δ>δ_lens)={F:.3f}  f_MC,max = {fm_str(row)}")
    hs = [r for r in out["rows"] if r["radius"] == "esfera dura"]
    mn = {b: min(r["f_MC_max"][b] for r in hs) for b in BOUNDS}
    strict, wide = list(BOUNDS)
    out["min_f_MC_max"] = mn
    out["verdict_U14"] = ("MICROLENTE_SIN_EFECTO" if mn[strict] >= 1 else "MICROLENTE_ACOTA_f_MC" if mn[wide] >= 1 else "MICROLENTE_EXIGE_f_MC_PEQUEÑA")
    print("U14:", out["verdict_U14"], {k: round(v, 3) for k, v in mn.items()})
    (OUT / "U14_RESULT.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


def fm_str(row): return ", ".join(f"{v:.2f}" for v in row["f_MC_max"].values())


if __name__ == "__main__":
    main()
