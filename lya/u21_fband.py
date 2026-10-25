"""U21 (docs/preregistros/ESPEC_U21.md): band of the Q scale F allowed by Lya (own CLASS gate), isocurvature (literature bound),
dark radiation and the CMB string bound, for the U17d collapse population.  F/F0 grid x variants P (B2 eta 1), H (B2, r = 3),
D (B2, eta 0.3), A (B4a eta 1); c = 3, sigma_ln = 0.5, R_rho = 1.56.

    python -B lya/u21_fband.py --workers 4      # cloud (reuses the CLASS cache for references and family)
    python -B lya/u21_fband.py --dry            # rows only
Output: lya/results/u21/U21_RESULT.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402
import u17b_class as U17B                                                     # noqa: E402

OUT = HERE / "results" / "u21"
FGRID = (1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03)
VARIANTS = {"P": ("B2", 1.0, 1.0), "H": ("B2", 1.0, 3.0), "D": ("B2", 0.3, 1.0), "A": ("B4a", 1.0, 1.0)}   # branch, eta, r_decay
C, R_RHO = 3.0, 1.56
ISO_STRICT, ISO_LAX, DNEFF_MAX, F_CMB, LISA_FA = 1e-17, 1.73e-18, 0.3, 1.26e15, 5e14
COLD = 0.05


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    import u17_ksw as U17
    from ucommon import benchmark
    from s4_dark import lambda_wdm_reference
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; lam_wdm, _ = lambda_wdm_reference(); U17B.SIG_LN = 0.5
    rows = []
    for v, (br, eta, rd) in VARIANTS.items():
        for fr in FGRID:
            F = F0 * fr
            r, x, n = U17B.row(F, Trh, br, C, eta, R_rho=R_RHO, r_decay=rd)
            s = U17.solve(F, Trh, br, 0.87, 10, c_eta=r["c_E"] / eta, R_rho=R_RHO)
            lam = s["lambda_collapse_Mpc"] * rd / lam_wdm
            r.update(variant=v, F_over_F0=fr, F_GeV=F, mu_over_Hstar=0.30 * fr, lambda_proxy_over_WDM=lam, cold_trivial=bool(lam < COLD),
                     DeltaNeff_emission=s["DeltaNeff_CMB_emission"], f_a_GeV=F / 2, LISA_ok=bool(F / 2 >= LISA_FA),
                     Lambda_H_keV=math.sqrt(r["m_G_eV"] * F / 2 * 1e9) / 1e3, iso_strict=bool(r["m_G_eV"] >= ISO_STRICT),
                     iso_lax=bool(r["m_G_eV"] >= ISO_LAX), dneff_ok=bool(s["DeltaNeff_CMB_emission"] <= DNEFF_MAX), cmb_strings_ok=bool(F <= F_CMB),
                     superradiance_band_no_selfint=("estelar 6e-13-1e-11" if 6e-13 <= r["m_G_eV"] <= 1e-11 else
                                                    "SMBH indirecta 1e-19-1e-13" if 1e-19 <= r["m_G_eV"] <= 1e-13 else "fuera"),
                     tag=f"U21_{v}_F{fr:g}")
            print(f"{r['tag']:<14} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.3g} MeV lam/WDM={lam:.3g} cold={r['cold_trivial']} f_nc={r['f_nc']:.3f} "
                  f"iso={r['iso_strict']} dNeff={r['DeltaNeff_emission']:.1e} Lambda_H={r['Lambda_H_keV']:.3g} keV", flush=True)
            r["_xn"] = (x, n); rows.append(r)
    ctl = next((r for r in rows if r["variant"] == "P" and r["cold_trivial"]), None)
    need = [r for r in rows if not r["cold_trivial"]] + ([ctl] if ctl else [])
    if dry:
        print(f"{len(need)} CLASS rows (incl. control {ctl['tag'] if ctl else None})"); return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U21.md", classy_version=getattr(classy, "__version__", "?"), F0_GeV=F0, grid=FGRID)
    P0 = U4.p_cdm()
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict)
    params = []
    for r in need:
        x, n = r["_xn"]; path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
        params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
    for r, P in zip(need, U4.pmap(params, workers, [r["tag"] for r in need], "row")):
        if P is None: r.update(class_error=True, LYA_primary=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
        print(f"{r['tag']:<14} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    control_ok = bool(ctl is None or ("deltaA" in ctl and ctl["deltaA"] < COLD))
    out["cold_control"] = dict(tag=ctl["tag"] if ctl else None, deltaA=ctl.get("deltaA") if ctl else None, ok=control_ok)
    if not control_ok:
        raise SystemExit("cold-trivial control failed: re-run with every row through CLASS (ESPEC_U21)")
    for r in rows:
        r.pop("_xn", None)
        if r["cold_trivial"] and "LYA_primary" not in r: r["LYA_primary"] = True; r["LYA_by"] = "frio_trivial"
        r["passes"] = bool(r.get("LYA_primary") and r["iso_strict"] and r["dneff_ok"] and r["cmb_strings_ok"])
    out["rows"] = rows
    bands = {}
    for v in VARIANTS:
        ok = [r["F_over_F0"] for r in rows if r["variant"] == v and r["passes"]]
        bands[v] = dict(F_max_over_F0=max(ok) if ok else None, F_min_over_F0=min(ok) if ok else None, n_pass=len(ok),
                        contiguous=bool(ok and all(r["passes"] for r in rows if r["variant"] == v and min(ok) <= r["F_over_F0"] <= max(ok))))
    out["bands"] = bands
    p = bands["P"]
    out["verdict_U21"] = ("BENCHMARK_PERMITIDO" if p["F_max_over_F0"] == 1.0 else "BANDA_BAJO_EL_BENCHMARK" if p["n_pass"] else "SIN_BANDA")
    out["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U21_RESULT.json").write_text(json.dumps(out, indent=1, default=float, allow_nan=False))
    print("U21:", out["verdict_U21"], bands, flush=True)


if __name__ == "__main__":
    main(sys.argv)
