"""U58b (docs/preregistros/ESPEC_U58b.md): U58 with the adopted cuscuton set only (mu = 0.16, m = 0.02, phi_inf = 0.25), plus
the amplitude control (V0 rescaled by (1 + eps_G)^3).  u58_cuscuton_inflacion is imported unchanged.
    python3 -B chain/u58b_cuscuton_inflacion.py   ->  results/u58b_cuscuton_inflacion.json
Internal, not exported."""
import json, time
from ucommon import save
import u58_cuscuton_inflacion as U58

ADOPTED = dict(mu=0.16, m=0.02, pinf=0.25)


def main():
    t0 = time.time()
    U58.SETS = {"GR": None, "adoptado": ADOPTED}
    R = {}
    for name in U58.SETS:
        R[name] = dict(main=U58.run_set(name, 1e-10), rtol9=U58.run_set(name, 1e-9))
        if U58.SETS[name]:
            R[name]["cs1"] = U58.run_set(name, 1e-10, cs_one=True)
    g = R["GR"]["main"]; m = R["adoptado"]["main"]; e = m["eps_G"]
    G0 = bool(abs(g["n_s"] - U58.B["n_s"]) <= 0.002 and abs(g["r"] / U58.B["r"] - 1) <= 0.10 and abs(g["P_zeta"] / g["P_sr"] - 1) <= 0.05)
    G1 = all(abs(R[n]["main"][q] / R[n]["rtol9"][q] - 1) <= 1e-4 for n in R for q in ("n_s", "r"))
    G2 = abs(m["P_zeta"] / R["adoptado"]["cs1"]["P_zeta"] - 1) <= 1e-4
    G3 = all(R[n]["main"]["res_max"] <= 1e-10 and R[n]["main"]["z2_min"] > 0 for n in R) and m["D_min"] > 0 and m["phid_min"] > 0
    # amplitude control: V0 * (1 + eps_G)^3
    V0 = U58.V0
    U58.V0 = V0 * (1 + e) ** 3
    amp = U58.run_set("adoptado", 1e-10)
    U58.V0 = V0
    allowed = U58.EPS_BOUND[0] <= e <= U58.EPS_BOUND[1]
    dn, dr = m["n_s"] - g["n_s"], m["r"] / g["r"] - 1
    if not (G0 and G1 and G2 and G3):
        verdict = "INCONCLUSO"
    elif not allowed:
        verdict = "CUSCUTON_EXCLUIDO"
    else:
        verdict = "CUSCUTON_COMPATIBLE" if abs(dn) <= 0.001 else "CUSCUTON_TENSO"
    out = dict(spec="docs/preregistros/ESPEC_U58b.md", adopted=ADOPTED, eps_G=e, lambda_minus_1=-2 * e / 3, allowed=allowed,
               gates=dict(G0=G0, G1=G1, G2=G2, G3=G3), G2_rel=m["P_zeta"] / R["adoptado"]["cs1"]["P_zeta"] - 1,
               dn_s=dn, dr_rel=dr, dr_rel_napkin=(1 + e) ** (2 / 3) - 1, dNeff_BBN=-43 / 7 * e,
               amplitude_control=dict(P_zeta_over_GR=amp["P_zeta"] / g["P_zeta"], dn_s=amp["n_s"] - g["n_s"], dr_rel=amp["r"] / g["r"] - 1),
               GR=dict(n_s=g["n_s"], r=g["r"]), runs=R, verdict=verdict, seconds=time.time() - t0)
    save("u58b_cuscuton_inflacion.json", out)
    print(json.dumps({k: out[k] for k in ("eps_G", "lambda_minus_1", "allowed", "gates", "G2_rel", "dn_s", "dr_rel", "dr_rel_napkin", "dNeff_BBN", "amplitude_control", "GR", "verdict")}, indent=1, default=float))


if __name__ == "__main__":
    main()
