"""U17g (docs/preregistros/ESPEC_U17g.md): the upper edge of the F band over the X interval of U19d (X_univ 244-1367), with
the U17f population, number and gates, for X ~ 162 Y at F0 up to Y = 8.44.  F grid of U17f plus 0.6; CLASS for every row that
is not cold-trivial (no dominance: U17e has no rows with Y > 2.24).  Control: the Y = 1.98 family reproduces U17f at 1e-3.
    python -B lya/u17g_class.py --dry           # rows, masses and which go to CLASS
    python -B lya/u17g_class.py --workers 4     # cloud (CLASS cache of U12/U21/U17e/U17f)
Output: lya/results/u17g/U17G_RESULT.json"""
import json, os, sys, time
from pathlib import Path
import numpy as np

os.environ["WHALEFALL_PLATEAU"] = "sum"                                       # as U17f (its enmienda 1)

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain")); sys.path.insert(0, str(ROOT / "network"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402
import u17e_class as E                                                        # noqa: E402
import u17f_class as G                                                        # noqa: E402

OUT = HERE / "results" / "u17g"
FGRID = (1.0, 0.85, 0.7, 0.6, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03)
Y_PRIMARY = (2.5, 3.0, 3.43, 4.25, 5.5, 8.44)                                 # X_univ(F0) ~ 405, 486, 556, 689, 891, 1367
Y_CONTROL = 1.98                                                              # U17f family: pipeline control
X_LO, X_MID, X_HI = 244.0, 556.0, 1367.0                                      # U19d: X_lo, X_univ, X_hi


def fam(Y):
    return f"Y{Y:.2f}"


def quadrature_patch():
    """the WHALEFALL quadrature patch (env/class_qm_trapz_log.patch) that the rows use (ncdm_quadrature_strategy = 4): an unpatched
    CLASS does not complain but leaves the momentum weights unset.  As env/check_env.py: the ncdm background with strategy 4 must
    match the default sampling to 2e-3 (two background-only runs, seconds; never cached, so it checks the build that runs now)"""
    import tempfile, classy
    with tempfile.TemporaryDirectory() as d:
        q = np.geomspace(1e-3, 40, 800); f = Path(d) / "fd.dat"; np.savetxt(f, np.column_stack([q, 1 / (np.exp(q) + 1)]))
        base = {"h": 0.674, "omega_b": 0.0224, "omega_cdm": 0.0, "A_s": 2.1e-9, "n_s": 0.9649, "N_ur": 3.044, "N_ncdm": 1, "m_ncdm": 3500.0,
                "T_ncdm": 0.1477 * (4 / 11) ** (1 / 3), "omega_ncdm": 0.12, "ncdm_fluid_approximation": 3, "use_ncdm_psd_files": 1,
                "ncdm_psd_filenames": str(f)}
        bg = []
        for extra in ({}, {"ncdm_quadrature_strategy": 4, "ncdm_N_momentum_bins": 45}):
            c = classy.Class(); c.set(dict(base, **extra)); c.compute(); b = c.get_background(); c.struct_cleanup(); c.empty()
            bg.append((np.array(b["z"]), np.array(b["(.)rho_ncdm[0]"]), np.array(b["(.)p_ncdm[0]"])))
    zz = np.geomspace(1e-2, 1e9, 60)
    rel = max(float(np.max(np.abs(np.interp(np.log(zz), np.log(z1[::-1] + 1e-30), x1[::-1]) / np.interp(np.log(zz), np.log(z0[::-1] + 1e-30), x0[::-1]) - 1)))
              for (z0, r0, p0), (z1, r1, p1) in [bg] for x0, x1 in ((r0, r1), (p0, p1)))
    return dict(max_rel_rho_p=rel, passed=bool(rel < 2e-3))


def verdict(bands):
    """F_max,bajo over the primary families (the hottest decides), and the edge per family"""
    fm = {Y: bands[fam(Y)]["F_max_over_F0"] for Y in Y_PRIMARY}
    if any(f is None for f in fm.values()): return "SIN_BANDA", None
    lo = min(fm.values())
    v = ("BORDE_0.7_SOSTENIDO" if lo >= 0.7 else "BORDE_0.6" if lo >= 0.6 else "BORDE_0.5" if lo >= 0.5 else
         f"BORDE_BAJA_A_{lo:g}")
    return v, lo


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    import u23_superradiance as U23
    from ucommon import benchmark
    from s4_dark import lambda_wdm_reference
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; lam_wdm, _ = lambda_wdm_reference()
    dark, light = U23.boundary("dark"), U23.boundary("light")
    shape = E.measured_shape("low"); print(f"shape: xbar={shape['xbar']:.3f} K_lattice={shape['K']:.3f}")
    rows = []
    for Y, primary in [(Y, True) for Y in Y_PRIMARY] + [(Y_CONTROL, False)]:
        for fr in FGRID:
            if not primary and fr not in (0.85, 0.7, 0.5): continue           # the control: three U17f rows
            F = F0 * fr; r, x, n = E.row(F, Trh, "B2", shape, fam(Y), Y)
            G.gates_cheap(r, F, fr, lam_wdm, dark, light)
            r.update(family=fam(Y), primary=primary, X_univ=r["K"] * r["xbar"], tag=f"U17g_{fam(Y)}_F{fr:g}",
                     lya_route="frio_trivial" if (r["cold_trivial"] and primary) else "CLASS")
            r["_xn"] = (x, n); rows.append(r)
            print(f"{r['tag']:<18} K={r['K']:7.2f} X={r['X_univ']:6.1f} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.3g} MeV lam/WDM={r['lambda_proxy_over_WDM']:.3g} "
                  f"f_iso={r['f_iso']['r2.4']:.2e} estelar={r['stellar']:<9} SMBH={'excl' if r['SMBH_CS_excluded'] else 'ok':<4} -> {r['lya_route']}", flush=True)
    need = [r for r in rows if r["lya_route"] == "CLASS"]
    print(f"{len(need)} CLASS rows", flush=True)
    if dry: return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U17g.md", classy_version=getattr(classy, "__version__", "?"), F0_GeV=F0, grid=FGRID,
               Y_primary=Y_PRIMARY, Y_control=Y_CONTROL, U19d=dict(X_lo=X_LO, X_univ=X_MID, X_hi=X_HI), shape=shape)
    out["quadrature_patch"] = quadrature_patch(); print("quadrature patch:", out["quadrature_patch"], flush=True)
    if not out["quadrature_patch"]["passed"] and out["classy_version"] != "MOCK":
        raise SystemExit("CLASS without the quadrature patch (env/class_qm_trapz_log.patch): no rows are computed or read from the cache")
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
        x, n = r["_xn"]; tag = r["tag"].replace("U17g", "U17f") if not r["primary"] else r["tag"]   # control: the U17f packets
        path, T, tail = U4.psd_from_packets(x, n, tag); r["tail_number_fraction"] = tail
        params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
    for r, P in zip(need, U4.pmap(params, workers, [r["tag"] for r in need], "row")):
        if P is None: r.update(class_error=True); continue
        r.update(U5.gates(P / P0, P0, famT, strict, wide))
        print(f"{r['tag']:<18} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    failed = [r["tag"] for r in need if r.get("class_error")]
    if failed:
        for r in rows: r.pop("_xn", None)
        (OUT / "U17G_PARTIAL.json").write_text(json.dumps(G._finite(dict(out, rows=rows, class_failed=failed)), indent=1, default=float, allow_nan=False))
        raise SystemExit(f"CLASS failed on {failed}: re-run (finished rows come from the cache); no verdict written")
    f_dA = {r["F_over_F0"]: r.get("deltaA") for r in json.loads((HERE / "results" / "u17f" / "U17F_RESULT.json").read_text())["rows"]
            if r["family"] == fam(Y_CONTROL)}
    out["pipeline_control"] = {f"F{r['F_over_F0']:g}": dict(U17f=f_dA[r["F_over_F0"]], U17g=r.get("deltaA"),
                                                           ok=bool(r.get("deltaA") is not None and abs(r["deltaA"] - f_dA[r["F_over_F0"]]) <= 1e-3))
                               for r in rows if not r["primary"]}
    print("pipeline control (Y = 1.98 vs U17f):", out["pipeline_control"], flush=True)
    for r in rows:
        r.pop("_xn", None)
        if r["lya_route"] != "CLASS": r["LYA_primary"] = True
        r["passes"] = bool(r.get("LYA_primary") and r["iso_ok"] and r["dneff_ok"] and r["cmb_strings_ok"] and r["superradiance_ok"])
    out["rows"] = rows
    prim = [f for f in (fam(Y) for Y in Y_PRIMARY)]
    out["bands"] = {f: G.band([r for r in rows if r["family"] == f]) for f in prim}
    out["bands_with_SMBH"] = {f: G.band([r for r in rows if r["family"] == f], smbh=True) for f in prim}
    v, lo = verdict(out["bands"])
    out["pipeline_ok"] = bool(all(c["ok"] for c in out["pipeline_control"].values()))
    if not out["pipeline_ok"]: v = "INVALIDO (control de la tubería)"
    out.update(verdict_U17g=v, F_max_low=lo, F_max_at_U19d_central=out["bands"][fam(3.43)]["F_max_over_F0"], runtime_h=(time.time() - t0) / 3600)
    (OUT / "U17G_RESULT.json").write_text(json.dumps(G._finite(out), indent=1, default=float, allow_nan=False))
    print("U17g:", v, {f: (b["F_min_contiguous_over_F0"] if b["n_pass"] else None, b["F_max_over_F0"]) for f, b in out["bands"].items()}, flush=True)


if __name__ == "__main__":
    main(sys.argv)
