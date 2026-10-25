"""U23 (docs/preregistros/ESPEC_U23.md): black-hole superradiance with self-interactions on the F band.
Primary gate: stellar BHs, the 2 sigma excluded region of Witte & Mummery 2025 (arXiv:2412.03655, Fig. 1), digitized in
data/literatura/witte_mummery_2412.03655_fig1.json; f_eff in [F/4, F] (f_a = F/2 with an O(1) RP3 normalization uncertainty).
Informative: SMBH window of Stott 2020 (no self-interaction) with the BGLS quenching threshold f_CD (Eq. 62, the SMBH mass giving alpha = 0.1
in [1e6, 1e10] Msun, a* = 0.9, T = t_Edd = 4e8 yr; enmienda 1).  Rows: U17e (Y, Y+, X) and U21 (P, H).
    python -B chain/u23_superradiance.py      ->  results/u23_superradiance.json"""
import json, math
import numpy as np
from ucommon import ROOT, save

FIG = ROOT / "data" / "literatura" / "witte_mummery_2412.03655_fig1.json"
YCLIP, MED = -11.4, 9
SMBH_WIN = (5.5e-20, 1.3e-16); A_STAR, T_EDD = 0.9, 4e8
CS_SMBH_MIN = 1.7e-17      # enmienda 2: Chathirathas & Schwetz (2511.15790) with Unal, Pacucci & Loeb (2012.12790), self-interactions included
KINDS = {"U17e": ("Y", "Y+", "X"), "U21": ("P", "H")}


def boundary(key):
    pts = np.array(json.loads(FIG.read_text())[key], dtype=float); x, y = pts[:, 0], pts[:, 1]
    ys = np.array([np.median(y[max(0, i - MED // 2): i + MED // 2 + 1]) for i in range(len(y))])
    return x, np.minimum(ys, YCLIP)


def inside(lm, linvf, region):
    x, y = region
    return bool(x[0] <= lm <= x[-1] and linvf < float(np.interp(lm, x, y)))


def f_cd(m_eV):
    """ESPEC_U23 enmienda 1: the SMBH mass is the one giving alpha = 0.1 (the fastest 211 growth), within [1e6, 1e10] Msun"""
    M = min(max(0.1 / (7.49e9 * m_eV), 1e6), 1e10); alpha = 7.49e9 * M * m_eV
    if alpha > 0.45: return alpha, M, float("inf")                                 # 211 not superradiant / outside Eq. 62
    return alpha, M, 3e14 * (1e10 / T_EDD) ** 0.5 * (1e-13 / m_eV) ** 0.5 * (0.01 / alpha) ** 2.5 * (0.9 / A_STAR) ** 0.75


def classify(m_eV, F_GeV, dark, light):
    lm = math.log10(m_eV); fs = np.geomspace(F_GeV / 4, F_GeV, 9)
    ins = [inside(lm, -math.log10(f), dark) for f in fs]
    if all(ins): st = "excluida"
    elif any(ins): st = "frontera"
    elif any(inside(lm, -math.log10(f), light) for f in fs): st = "probable"
    else: st = "libre"
    out = dict(stellar=st, f_eff_range_GeV=[F_GeV / 4, F_GeV], SMBH_CS_excluded=bool(m_eV < CS_SMBH_MIN))
    if SMBH_WIN[0] <= m_eV <= SMBH_WIN[1]:
        a, M, fc = f_cd(m_eV); out.update(SMBH_alpha=a, SMBH_M_Msun=M, SMBH_f_CD_GeV=fc, SMBH_possible=bool(F_GeV / 2 > fc), SMBH_possible_fF=bool(F_GeV > fc))
    return out


def verdict(rows, smbh=False):
    band = sorted([r for r in rows if r["passes"]], key=lambda r: r["F_over_F0"])
    if not band: return "SIN_BANDA_PREVIA", None
    surv = [r for r in band if r["stellar"] != "excluida" and not (smbh and r["SMBH_CS_excluded"])]
    if len(surv) == len(band): return "SUPERRADIANCIA_NO_RECORTA", [band[0]["F_over_F0"], band[-1]["F_over_F0"]]
    if not surv: return "SUPERRADIANCIA_CIERRA", None
    idx = [band.index(r) for r in surv]; contiguous = idx == list(range(idx[0], idx[-1] + 1))
    rng = [surv[0]["F_over_F0"], surv[-1]["F_over_F0"]]
    if contiguous and surv[-1] is band[-1]: return "SUPERRADIANCIA_RECORTA_EL_BORDE_INFERIOR", rng
    if contiguous: return "SUPERRADIANCIA_RECORTA (borde superior; caso no previsto en el ESPEC)", rng
    return "SUPERRADIANCIA_PARTE_LA_BANDA", rng


def main():
    dark, light = boundary("dark"), boundary("light")
    src = {"U17e": ROOT / "lya" / "results" / "u17e" / "U17E_RESULT.json", "U21": ROOT / "lya" / "results" / "u21" / "U21_RESULT.json"}
    out = dict(spec="docs/preregistros/ESPEC_U23.md", figure=str(FIG.relative_to(ROOT)), dark_mass_range_eV=[10 ** dark[0][0], 10 ** dark[0][-1]],
               light_mass_range_eV=[10 ** light[0][0], 10 ** light[0][-1]], rows=[], verdicts={})
    for s, kinds in KINDS.items():
        if not src[s].exists(): print("missing", src[s]); continue
        res = json.loads(src[s].read_text())
        for k in kinds:
            rk = [r for r in res["rows"] if r.get("kind", r.get("variant")) == k]
            rows = []
            for r in rk:
                c = classify(r["m_G_eV"], r["F_GeV"], dark, light)
                rows.append(dict(source=s, kind=k, F_over_F0=r["F_over_F0"], F_GeV=r["F_GeV"], m_G_eV=r["m_G_eV"], passes=bool(r["passes"]), **c))
            v, rng = verdict(rows); v2, rng2 = verdict(rows, smbh=True)              # enmienda 3: with the SMBH (Unal et al.) gate
            out["verdicts"][f"{s}_{k}"] = dict(verdict=v, band_after=rng, verdict_with_SMBH=v2, band_after_with_SMBH=rng2); out["rows"] += rows
            for r in rows:
                print(f"{s}_{k:<3} F/F0={r['F_over_F0']:<7g} m={r['m_G_eV']:.2e} f_a={r['F_GeV'] / 2:.2e} pasa_antes={r['passes']!s:<5} "
                      f"estelar={r['stellar']:<9} SMBH(C&S)={'excl' if r['SMBH_CS_excluded'] else 'ok':<4} {'SMBH f_CD=%.1e posible=%s' % (r['SMBH_f_CD_GeV'], r['SMBH_possible']) if 'SMBH_f_CD_GeV' in r else ''}")
            print(f"  -> {s}_{k}: {v} {rng} | con supermasivos: {v2} {rng2}")
    if "U17e_Y" in out["verdicts"]:
        out["verdict_U23"] = out["verdicts"]["U17e_Y"]["verdict"]; out["verdict_U23_with_SMBH"] = out["verdicts"]["U17e_Y"]["verdict_with_SMBH"]
    save("u23_superradiance.json", out)
    print("U23:", out.get("verdict_U23"))


if __name__ == "__main__":
    main()
