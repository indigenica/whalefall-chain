"""U17f (docs/preregistros/ESPEC_U17f.md): the upper edge of the F band with the U19b collapse extrapolated to the universe,
over the whole interval X_univ ~ 170-320 (Y = 1.05, 1.20, 1.36, 1.70, 1.98; X ~ 162 Y at F0) and, informative, reading (ii)
of U19b section 8 (Y = 1.0, 0.5).  Population, number and branch as the U17e Y rows (measured shape 'low', B2, U17 scaling).
F grid of U21 plus 0.85.  CLASS for F in {1, 0.85, 0.7, 0.5}; F <= 0.35 by dominance over the U17e Y+ rows, with one CLASS
control (Y = 1.98, F = 0.35); U21 cold-trivial rule.  Gates: U12 Lya, U20b f_iso < 0.0084 (C = 1, r = 2.40), DeltaNeff <= 0.3,
CMB strings, U23 stellar superradiance not 'excluida'.
    python -B lya/u17f_class.py --dry           # rows, masses and which go to CLASS
    python -B lya/u17f_class.py --workers 4     # cloud (CLASS cache of U12/U21/U17e for references and family)
Output: lya/results/u17f/U17F_RESULT.json"""
import json, os, sys, time
from pathlib import Path
import numpy as np

os.environ["WHALEFALL_PLATEAU"] = "sum"                                       # enmienda 1: the plateau of U17e and U21

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain")); sys.path.insert(0, str(ROOT / "network"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402
import u17e_class as E                                                        # noqa: E402

OUT = HERE / "results" / "u17f"
FGRID = (1.0, 0.85, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03)
F_CLASS = (1.0, 0.85, 0.7, 0.5)
Y_PRIMARY = (1.05, 1.20, 1.36, 1.70, 1.98)                                    # X_univ(F0) ~ 170, 195, 221, 275, 320
Y_INFO = (1.0, 0.5)                                                            # reading (ii) of U19b section 8
Y_DOMINANCE = 2.24                                                             # the U17e Y+ rows
CONTROL = (1.98, 0.35)
FISO_MAX, R_ISO, R_ISO_ALT = 0.0084, 2.40, 1.78


def fam(Y):
    return f"Y{Y:.2f}"


def gates_cheap(r, F, fr, lam_wdm, dark, light):
    import u20b_iso as U20B
    import u23_superradiance as U23
    from s4_dark import lambda_fs
    lam = float(lambda_fs(np.array([r["xbar"] * r["m_G_eV"] * r["a_osc"]]), r["m_G_eV"], np.array([r["a_osc"]]))[0]) / lam_wdm
    fi = {f"r{rr:g}": U20B.fiso(r["T_osc_MeV"], rr) for rr in (R_ISO, R_ISO_ALT)}
    r.update(F_over_F0=fr, F_GeV=F, lambda_proxy_over_WDM=lam, cold_trivial=bool(lam < E.COLD), f_a_GeV=F / 2, LISA_ok=bool(F / 2 >= E.LISA_FA),
             dneff_ok=bool(r["DeltaNeff_emission"] <= E.DNEFF_MAX), cmb_strings_ok=bool(F <= E.F_CMB), f_iso=fi,
             iso_ok=bool(fi[f"r{R_ISO:g}"] < FISO_MAX), iso_ok_rW=bool(fi[f"r{R_ISO_ALT:g}"] < FISO_MAX), **U23.classify(r["m_G_eV"], F, dark, light))
    r["superradiance_ok"] = r["stellar"] != "excluida"


def band(rk, smbh=False):
    ok = [r for r in rk if r["passes"] and not (smbh and r["SMBH_CS_excluded"])]
    fs = sorted(r["F_over_F0"] for r in ok)
    if not fs: return dict(F_max_over_F0=None, F_min_over_F0=None, n_pass=0, n_rows=len(rk), contiguous=False, F_max_lya_only=None)
    top = max(fs); seg = []
    for r in sorted(rk, key=lambda r: -r["F_over_F0"]):
        if r["F_over_F0"] > top: continue
        if not (r["passes"] and not (smbh and r["SMBH_CS_excluded"])): break
        seg.append(r["F_over_F0"])
    return dict(F_max_over_F0=top, F_min_over_F0=min(fs), F_min_contiguous_over_F0=min(seg), n_pass=len(fs), n_rows=len(rk),
                contiguous=bool(len(seg) == len(fs)), passing=fs,
                F_max_lya_only=max([r["F_over_F0"] for r in rk if r.get("LYA_primary")], default=None))


def verdict(bands):
    fm = [bands[fam(Y)]["F_max_over_F0"] for Y in Y_PRIMARY]
    if any(f is None for f in fm): return "SIN_BANDA", None, None
    lo, hi = min(fm), max(fm)
    v = "BENCHMARK_PERMITIDO" if lo >= 1.0 else "BORDE_SUBE" if lo >= 0.85 else "BORDE_0.7_CONFIRMADO" if lo >= 0.7 else "BORDE_BAJA"
    if hi >= 1.0 > lo:
        xs = max(round(162 * Y) for Y in Y_PRIMARY if bands[fam(Y)]["F_max_over_F0"] >= 1.0)
        v += f"; F0 PERMITIDO CON X <= {xs}"
    return v, lo, hi


def _finite(o):
    """JSON without inf: U23's f_cd is inf where the 211 mode of the SMBH is not superradiant (alpha > 0.45); written as null"""
    if isinstance(o, dict): return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_finite(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o): return None
    return o


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    import u23_superradiance as U23
    from ucommon import benchmark
    from s4_dark import lambda_wdm_reference
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; lam_wdm, _ = lambda_wdm_reference()
    dark, light = U23.boundary("dark"), U23.boundary("light")
    shape = E.measured_shape("low"); print(f"shape: xbar={shape['xbar']:.3f} K_lattice={shape['K']:.3f}")
    u17e = json.loads((HERE / "results" / "u17e" / "U17E_RESULT.json").read_text())
    dom = {r["F_over_F0"]: bool(r.get("LYA_primary")) for r in u17e["rows"] if r["kind"] == "Y+"}
    rows = []
    for Y, primary in [(Y, True) for Y in Y_PRIMARY] + [(Y, False) for Y in Y_INFO]:
        for fr in FGRID:
            F = F0 * fr; r, x, n = E.row(F, Trh, "B2", shape, fam(Y), Y)
            gates_cheap(r, F, fr, lam_wdm, dark, light)
            r.update(family=fam(Y), primary=primary, X_univ=r["K"] * r["xbar"], tag=f"U17f_{fam(Y)}_F{fr:g}")
            if fr in F_CLASS: r["lya_route"] = "CLASS"
            elif r["cold_trivial"]: r["lya_route"] = "frio_trivial"
            elif Y <= Y_DOMINANCE and dom.get(fr): r["lya_route"] = "dominancia_U17e_Y+"
            else: r["lya_route"] = "CLASS"
            if (Y, fr) == CONTROL: r["lya_route"] = "CLASS"; r["dominance_control"] = True
            r["_xn"] = (x, n); rows.append(r)
            print(f"{r['tag']:<18} K={r['K']:7.2f} X={r['X_univ']:6.1f} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.3g} MeV lam/WDM={r['lambda_proxy_over_WDM']:.3g} "
                  f"f_iso={r['f_iso']['r2.4']:.2e} estelar={r['stellar']:<9} SMBH={'excl' if r['SMBH_CS_excluded'] else 'ok':<4} -> {r['lya_route']}", flush=True)
    need = [r for r in rows if r["lya_route"] == "CLASS"]
    print(f"{len(need)} CLASS rows", flush=True)
    if dry: return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U17f.md", classy_version=getattr(classy, "__version__", "?"), F0_GeV=F0, grid=FGRID,
               Y_primary=Y_PRIMARY, Y_informative=Y_INFO, shape=shape)
    P0 = U4.p_cdm()
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    famT = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict)
    params = []
    for r in need:
        x, n = r["_xn"]; path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
        params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
    for r, P in zip(need, U4.pmap(params, workers, [r["tag"] for r in need], "row")):
        if P is None: r.update(class_error=True); continue                     # an execution failure, not a physical exclusion
        r.update(U5.gates(P / P0, P0, famT, strict, wide))
        print(f"{r['tag']:<18} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    failed = [r["tag"] for r in need if r.get("class_error")]
    if failed:                                    # no verdict with missing rows: a failed CLASS run must not lower F_max
        for r in rows: r.pop("_xn", None)
        (OUT / "U17F_PARTIAL.json").write_text(json.dumps(_finite(dict(out, rows=rows, class_failed=failed)), indent=1, default=float, allow_nan=False))
        raise SystemExit(f"CLASS failed on {failed}: re-run (finished rows come from the cache); no verdict written")
    e_dA = {r["F_over_F0"]: r.get("deltaA") for r in u17e["rows"] if r["kind"] == "Y"}
    out["pipeline_control"] = {f"F{fr:g}": dict(U17e=e_dA[fr], U17f=r.get("deltaA"), ok=bool(r.get("deltaA") is not None and abs(r["deltaA"] - e_dA[fr]) <= 1e-3))
                               for r in rows for fr in (1.0, 0.7, 0.5) if r["family"] == fam(1.70) and r["F_over_F0"] == fr}
    print("pipeline control (Y = 1.70 vs U17e):", out["pipeline_control"], flush=True)
    ctl = next(r for r in rows if r.get("dominance_control"))
    control_ok = bool(ctl.get("LYA_primary"))
    out["dominance_control"] = dict(tag=ctl["tag"], deltaA=ctl.get("deltaA"), ok=control_ok)
    for r in rows:
        r.pop("_xn", None)
        if r["lya_route"] != "CLASS": r["LYA_primary"] = True
        r["passes"] = bool(r.get("LYA_primary") and r["iso_ok"] and r["dneff_ok"] and r["cmb_strings_ok"] and r["superradiance_ok"])
    out["rows"] = rows
    if not control_ok:
        (OUT / "U17F_PARTIAL.json").write_text(json.dumps(out, indent=1, default=float))
        raise SystemExit("dominance control failed: every F <= 0.35 row that is not cold-trivial goes through CLASS (ESPEC_U17f)")
    fams = sorted({r["family"] for r in rows}, key=lambda f: -float(f[1:]))
    out["bands"] = {f: band([r for r in rows if r["family"] == f]) for f in fams}
    out["bands_with_SMBH"] = {f: band([r for r in rows if r["family"] == f], smbh=True) for f in fams}
    v, lo, hi = verdict(out["bands"])
    out.update(verdict_U17f=v, F_max_low=lo, F_max_high=hi, runtime_h=(time.time() - t0) / 3600)
    (OUT / "U17F_RESULT.json").write_text(json.dumps(_finite(out), indent=1, default=float, allow_nan=False))
    print("U17f:", v, {f: (b["F_min_contiguous_over_F0"] if b["n_pass"] else None, b["F_max_over_F0"]) for f, b in out["bands"].items()}, flush=True)


if __name__ == "__main__":
    main(sys.argv)
