"""U17d (docs/preregistros/ESPEC_U17d.md): the decisive collapse rows of U17b/U17c with the energy balance corrected
(number = eta rho_net / <E>, <E> of the log-normal whose mean momentum is c m_G; codex audit 14 Oct 2026) and the two bridge
hypotheses declared: R_rho = rho_RP3/rho_U1 of the network at the collapse, and r = a_decay/a_osc (momentum c m_G at a_decay,
energy frozen at H = m_G/3: the Hiramatsu et al. prescription).  c = 3, sigma_ln = 0.5, q = 0.87, x0 = 10, R_Gamma = 1.56.

    python -B lya/u17d_class.py --workers 4      # cloud or WSL (reuses the U12 CLASS cache for references and family)
    python -B lya/u17d_class.py --dry            # rows only, no CLASS
Output: lya/results/u17d/U17D_RESULT.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402
import u17b_class as U17B                                                     # noqa: E402

OUT = HERE / "results" / "u17d"
C = 3.0; R_PRIMARY = 1.56


def specs():
    """(branch, eta, R_rho, r_decay, legacy)"""
    out = [("B2", e, R, 1.0, False) for R in (R_PRIMARY, 1.0) for e in (1.0, 0.3, 0.1)]
    out += [("B2", e, R_PRIMARY, r, False) for r in (2.0, 4.0) for e in (1.0, 0.3, 0.1)]
    out += [("B4a", e, R_PRIMARY, 1.0, False) for e in (0.1, 0.03)] + [("B4a", 0.1, 1.0, 1.0, False)]
    return out + [("B2", 1.0, 1.0, 1.0, True)]                              # control: U17b row as run (deltaA 0.247)


def x_var(r):
    """the U17c reading variable, with the corrections: (8 pi/27) R_rho xi ell eta (c/c_E) r = K_eff c r"""
    return r["K_eff"] * r["c"] * r["r_decay"]


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    dry = "--dry" in argv; t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    from ucommon import benchmark
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; U17B.SIG_LN = 0.5
    rows, params = [], []
    for b, e, R, rd, leg in specs():
        r, x, n = U17B.row(F, Trh, b, C, e, legacy_energy=leg, R_rho=R, r_decay=rd)
        r["tag"] = f"U17d_{b}_eta{e:g}_R{R:g}_r{rd:g}" + ("_legacy" if leg else ""); r["X"] = x_var(r)
        print(f"{r['tag']:<30} c_E={r['c_E']:.3f} K={r['K_eff']:8.2f} X={r['X']:7.1f} m_G={r['m_G_eV']:.3e} T_osc={r['T_osc_MeV'] * 1e3:.1f} keV "
              f"f_nc={r['f_nc']:.4f}", flush=True)
        if not dry:
            path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
            params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0)))
        rows.append(r)
    if dry: return
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U17d.md", classy_version=getattr(classy, "__version__", "?"), c=C, sigma_ln=0.5)
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
    for r, P in zip(rows, U4.pmap(params, workers, [r["tag"] for r in rows], "row")):
        if P is None: r.update(class_error=True, LYA_primary=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
        print(f"{r['tag']:<30} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    out["rows"] = rows
    ctl = [r for r in rows if r["legacy_energy"]][0]
    out["control_ok"] = bool("deltaA" in ctl and abs(ctl["deltaA"] - 0.2470) < 0.002)
    prim = {r["eta"]: r.get("LYA_primary", False) for r in rows if r["branch"] == "B2" and r["R_rho"] == R_PRIMARY and r["r_decay"] == 1.0}
    eta_ok = [e for e, ok in prim.items() if ok]
    out["eta_crit_B2_primary_grid"] = max(eta_ok) if eta_ok else None
    if not out["control_ok"] or any(r.get("class_error") for r in rows): v = "INCONCLUSO (control o CLASS)"
    elif prim.get(1.0): v = "RAMA_VIVA_CON_COLAPSO_TIPO_U1"
    elif prim.get(0.3) or prim.get(0.1): v = "RAMA_VIVA_REQUIERE_COLAPSO_DURO"
    else: v = "RAMA_VIVA_EN_PELIGRO"
    out["verdict_U17d"] = v; out["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U17D_RESULT.json").write_text(json.dumps(out, indent=1, default=float, allow_nan=False))
    print("U17d:", v, "eta_crit(B2, primary grid) =", out["eta_crit_B2_primary_grid"], flush=True)


if __name__ == "__main__":
    main(sys.argv)
