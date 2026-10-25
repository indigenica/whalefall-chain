"""U20b (docs/preregistros/ESPEC_U20b.md): white-noise isocurvature against the observational f_iso bounds.
f_iso = C (k_CMB/k*)^3 / A_s, k* = a_osc H_osc / r (as U20), k_CMB = 0.05/Mpc; primary gate Garcia-Gallego et al. 2026
(arXiv:2603.04401) Lya f_iso < 0.0084 (95 %); informative: Lya 4e-3 and UFD 6e-4 (Chathirathas & Schwetz, arXiv:2511.15790,
Table 1), GG tentative detection 0.0064 (+0.0012 -0.0014).  Rows: U17e (Y, Y+, X; r = r_N = 2.40) and U21 (P r = 1, H r = 3).
    python -B chain/u20b_iso.py      ->  results/u20b_iso.json"""
import json
from ucommon import ROOT, save
import u20_iso as U20

K_CMB, A_S = 0.05, 2.10e-9
CS = (0.3, 1.0, 3.0, 10.0); R_N, R_W = 2.40, 1.78
GG, LYA_OLD, UFD, DET = 0.0084, 4e-3, 6e-4, (0.0050, 0.0076)
FGRID = (1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03)


def fiso(T_MeV, r, C=1.0):
    return C * (K_CMB / U20.kstar(T_MeV * 1e-3, r)) ** 3 / A_S


def main():
    e = json.loads((ROOT / "lya" / "results" / "u17e" / "U17E_RESULT.json").read_text())
    u21 = json.loads((ROOT / "lya" / "results" / "u21" / "U21_RESULT.json").read_text())
    p23 = ROOT / "results" / "u23_superradiance.json"; u23 = json.loads(p23.read_text()) if p23.exists() else None
    st = {(r["source"], r["kind"], r["F_over_F0"]): r for r in (u23 or {}).get("rows", [])}
    rows = []
    for src, res, kinds in (("U17e", e, {"Y": R_N, "Y+": R_N, "X": R_N}), ("U21", u21, {"P": 1.0, "H": 3.0})):
        for r in res["rows"]:
            k = r.get("kind", r.get("variant"))
            if k not in kinds: continue
            rr = kinds[k]; ks = U20.kstar(r["T_osc_MeV"] * 1e-3, rr)
            f = {f"C{C:g}": fiso(r["T_osc_MeV"], rr, C) for C in CS}; f1 = f["C1"]
            s = st.get((src, k, r["F_over_F0"]), {})
            row = dict(source=src, kind=k, F_over_F0=r["F_over_F0"], m_G_eV=r["m_G_eV"], T_osc_MeV=r["T_osc_MeV"], r=rr, k_star_Mpc=ks,
                       f_iso=f, f_iso_rW=fiso(r["T_osc_MeV"], R_W) if src == "U17e" else None, valid_lya=bool(ks > 10), valid_ufd=bool(ks > 1000),
                       iso_GG=bool(f1 < GG), iso_lya_old=bool(f1 < LYA_OLD), iso_UFD=bool(f1 < UFD), GG_detection_compatible=bool(DET[0] <= f1 <= DET[1]),
                       LYA_primary=bool(r.get("LYA_primary")), dneff_ok=bool(r.get("dneff_ok", True)), cmb_ok=bool(r.get("cmb_strings_ok", True)),
                       stellar=s.get("stellar"), SMBH_CS_excluded=s.get("SMBH_CS_excluded"))
            row["combined"] = bool(row["LYA_primary"] and row["iso_GG"] and row["dneff_ok"] and row["cmb_ok"] and row["stellar"] != "excluida")
            rows.append(row)
    out = dict(spec="docs/preregistros/ESPEC_U20b.md", rows=rows, bands={})
    for src, k in (("U17e", "Y"), ("U17e", "Y+"), ("U17e", "X"), ("U21", "P"), ("U21", "H")):
        rk = sorted([r for r in rows if r["source"] == src and r["kind"] == k], key=lambda r: -r["F_over_F0"])
        iso = [r["F_over_F0"] for r in rk if r["iso_GG"]]; comb = [r["F_over_F0"] for r in rk if r["combined"]]
        out["bands"][f"{src}_{k}"] = dict(F_max_iso=max(iso) if iso else None, combined=sorted(comb, reverse=True),
                                          combined_SMBH_flag=[r["F_over_F0"] for r in rk if r["combined"] and r["SMBH_CS_excluded"]],
                                          combined_UFD_fail=[r["F_over_F0"] for r in rk if r["combined"] and not r["iso_UFD"]],
                                          GG_detection_compatible=[r["F_over_F0"] for r in rk if r["GG_detection_compatible"]])
    y = out["bands"]["U17e_Y"]["F_max_iso"]
    out["verdict_U20b"] = "ISOCURVATURA_NO_LIMITA" if y == 1.0 else f"ISOCURVATURA_LIMITA (F_max = {y} F0)"
    save("u20b_iso.json", out)
    for r in rows:
        print(f"{r['source']}_{r['kind']:<3} F/F0={r['F_over_F0']:<7g} m={r['m_G_eV']:.2e} k*={r['k_star_Mpc']:.3g} f_iso(C=1)={r['f_iso']['C1']:.2e} "
              f"GG={r['iso_GG']!s:<5} UFD={r['iso_UFD']!s:<5} Lya={r['LYA_primary']!s:<5} estelar={r['stellar']} SMBH={r['SMBH_CS_excluded']} -> {r['combined']}")
    for k, b in out["bands"].items(): print(k, b)
    print("U20b:", out["verdict_U20b"])


if __name__ == "__main__":
    main()
