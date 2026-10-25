"""U20 (docs/preregistros/ESPEC_U20.md): own order-of-magnitude white-noise isocurvature estimate on the U21 rows.
Delta2_iso(k) = C (k/k*)^3 D_iso^2, k* = a_osc H_osc / r, D_iso = 1 + 1.5 g_b (1 + z_eq)/(1 + z) at z = 3; criterion Delta2_iso/Delta2_ad <= 0.1
at k = 10 h/Mpc (Delta2_ad from the CLASS LCDM spectrum cached by U4a).  No free-streaming (conservative).
    python -B chain/u20_iso.py      ->  results/u20_iso.json"""
import json, math, sys
import numpy as np
from ucommon import ROOT, save
from s4_dark import scale, gstar
sys.path.insert(0, str(ROOT / "lya"))

GEV_TO_MPCINV = 1.5637e38; MPL_FULL = 1.2209e19; H_LITTLE = 0.674; Z_EQ, Z, G_B, KH = 3400.0, 3.0, 0.9, 10.0
CS = (0.3, 1.0, 3.0, 10.0); RATIO_MAX = 0.10; R_DECAY = {"P": 1.0, "H": 3.0}


def kstar(T_GeV, r):
    H = 1.66 * math.sqrt(gstar(T_GeV)) * T_GeV ** 2 / MPL_FULL
    return float(scale(T_GeV)) * H * GEV_TO_MPCINV / r


def main():
    import u4a_class as U4
    P0 = U4.p_cdm(); k = U4.KH * H_LITTLE; D2ad = float(np.interp(math.log(KH), np.log(U4.KH), k ** 3 * P0 / (2 * math.pi ** 2)))
    D_iso = 1 + 1.5 * G_B * (1 + Z_EQ) / (1 + Z); kk = KH * H_LITTLE
    u21 = json.loads((ROOT / "lya" / "results" / "u21" / "U21_RESULT.json").read_text())
    rows = []
    for r in u21["rows"]:
        if r["variant"] not in R_DECAY: continue
        ks = kstar(r["T_osc_MeV"] * 1e-3, R_DECAY[r["variant"]])
        ratio = {f"{C:g}": C * (kk / ks) ** 3 * D_iso ** 2 / D2ad for C in CS}
        rows.append(dict(variant=r["variant"], F_over_F0=r["F_over_F0"], m_G_eV=r["m_G_eV"], T_osc_MeV=r["T_osc_MeV"], k_star_Mpc=ks,
                         ratio=ratio, passes_C1=bool(ratio["1"] <= RATIO_MAX), U21_passes=r["passes"]))
    out = dict(spec="docs/preregistros/ESPEC_U20.md", Delta2_ad_k10=D2ad, D_iso=D_iso, rows=rows)
    for v in R_DECAY:
        rv = sorted([x for x in rows if x["variant"] == v], key=lambda x: x["m_G_eV"])
        mins = {}
        for C in CS:
            ok = [x for x in rv if x["ratio"][f"{C:g}"] <= RATIO_MAX]
            if ok:   # interpolate in ln m between the last failing and the first passing row (ratio ~ m^-3/2)
                j = rv.index(ok[0]); a = rv[j - 1] if j > 0 else ok[0]; b = ok[0]
                if a is b: mins[f"{C:g}"] = b["m_G_eV"]
                else:
                    ra, rb = math.log(a["ratio"][f"{C:g}"]), math.log(b["ratio"][f"{C:g}"]); t = (math.log(RATIO_MAX) - ra) / (rb - ra)
                    mins[f"{C:g}"] = math.exp(math.log(a["m_G_eV"]) + t * (math.log(b["m_G_eV"]) - math.log(a["m_G_eV"])))
            else: mins[f"{C:g}"] = None
        band = [x["F_over_F0"] for x in rv if x["passes_C1"] and x["U21_passes"]]
        out[f"variant_{v}"] = dict(m_min_eV=mins, F_max_over_F0_C1=max(band) if band else None)
    p = out["variant_P"]["F_max_over_F0_C1"]
    out["verdict_U20"] = ("ISOCURVATURA_CIERRA_BANDA" if p is None else "ISOCURVATURA_CONFIRMA_BANDA" if p >= 0.7 else "ISOCURVATURA_RECORTA_BANDA")
    save("u20_iso.json", out)
    for x in rows: print(f"{x['variant']} F/F0={x['F_over_F0']:<6g} m={x['m_G_eV']:.2e} k*={x['k_star_Mpc']:.3g}/Mpc ratio(C=1)={x['ratio']['1']:.3g} pass={x['passes_C1']}")
    print("Delta2_ad(k=10h/Mpc,z=3) =", round(D2ad, 3), "D_iso =", round(D_iso, 1))
    print({v: out[f"variant_{v}"] for v in R_DECAY}); print("U20:", out["verdict_U20"])


if __name__ == "__main__":
    main()
