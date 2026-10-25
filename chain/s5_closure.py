"""S5 · conditional cosmological closure (ESPEC_CIERRE.md §2, preregistered 2026-10-06 06:55 UTC).

For every node (spectrum, x0, K_sw, eps) find m_G with Omega_mis + Omega_em = 0.12 and apply gates P1-P6."""
import itertools, json, math, time
from ucommon import load, save, benchmark, provenance
import numpy as np
import c35b_pseudo_goldstone_dm_rp3 as C35
from s4_dark import source, abundance, at_epoch, dm_root, emitted_lambda, lambda_wdm_reference

SPECTRA = ("ir_delta", "energy_logflat"); X0S = (10, 30, 100); KSW = (0, 1, 10)
EPS = [0.0] + list(10 ** np.linspace(-4, 0, 17))
DNEFF_CMB, DNEFF_BBN = 0.29, 0.38
M_LYA = 2e-20
SR_WINDOWS = [(2.5e-21, 1.2e-20), (5.5e-20, 1.3e-16), (3.8e-14, 3.4e-11)]   # Stott 2020, arXiv:2009.07206 (68 %); ESPEC miscited it as Unal et al.
SR_FA_MIN = 1.5e15   # preregistered single threshold (kept for the gate); Stott gives per-population thresholds below
SR_FA_STOTT = dict(xray_binaries=1.2e15, gw_mergers=2.5e14, smbh=2.0e16)  # audit correction (audit_u1/AUDIT_U1.md 4f)
F_EM_FREE = 0.1
LISA_FA = 5e14; CMBS4_DETECT = 0.06


def node(F, Trh, dneff_late, lam_wdm, spec, x0, K, eps, step=0.03, bins=64):
    try:
        m = dm_root(F, Trh, eps, spec, x0, K, step)
    except ValueError:
        return dict(spectrum=spec, x0=x0, K_sw=K, eps=eps, P1=False)
    mis, Tosc = C35.omega(m, F, K)
    pk = source(F, m, Trh, x0, eps, step, spec, bins) if eps else None
    em = abundance(pk) if eps else 0.0
    cmb = at_epoch(pk, 0.26e-9)["DeltaNeff_radiation_moment"] if eps else 0.0
    bbn = at_epoch(pk, 1e-3)["DeltaNeff_radiation_moment"] if eps else 0.0
    lam = emitted_lambda(pk) if eps else 0.0
    f_em = em / 0.12
    r = dict(spectrum=spec, x0=x0, K_sw=K, eps=eps, m_G_eV=m, T_osc_GeV=Tosc, T_stop_GeV=(pk["T_stop_GeV"] if pk else None),
             Omega_mis=mis, Omega_em=em, f_em=f_em, DeltaNeff_CMB=cmb + dneff_late, DeltaNeff_BBN=bbn + dneff_late,
             lambda_em_Mpc=lam, lambda_em_over_WDM=lam / lam_wdm, M_minicluster_Msun=C35.minicluster_Msun(m), P1=True)
    r["P2"] = r["DeltaNeff_CMB"] < DNEFF_CMB
    r["P3"] = r["DeltaNeff_BBN"] < DNEFF_BBN
    r["P4"] = m >= M_LYA
    r["P5"] = bool(f_em <= F_EM_FREE or lam < lam_wdm)
    r["P6_window_hit"] = any(lo <= m <= hi for lo, hi in SR_WINDOWS)
    r["conservative"] = all(r[k] for k in ("P1", "P2", "P3", "P4", "P5"))
    r["strict"] = r["conservative"] and not r["P6_window_hit"]
    r["CMBS4_detectable"] = r["DeltaNeff_CMB"] >= CMBS4_DETECT
    return r


def _contig(sub):
    """largest eps such that every grid node from eps = 0 up to it is allowed (P1-P5)"""
    srt = sorted(sub, key=lambda r: r["eps"]); best = None
    for r in srt:
        if not r.get("conservative"): break
        best = r["eps"]
    return best


def _preds(rows):
    if not rows: return None
    g = lambda k: [min(r[k] for r in rows), max(r[k] for r in rows)]
    return dict(m_G_eV=g("m_G_eV"), DeltaNeff_CMB=g("DeltaNeff_CMB"), f_em=g("f_em"), lambda_em_Mpc=g("lambda_em_Mpc"),
                M_minicluster_Msun=g("M_minicluster_Msun"), T_osc_GeV=g("T_osc_GeV"))


def main():
    B = benchmark(); F = B["F_GeV"]; Trh = B["T_Gamma_GeV"]; fa = B["f_a_GeV"]
    s3 = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text())
    dneff_late = s3["DeltaNeff_late_homogeneous_max"]
    lam_wdm, _ = lambda_wdm_reference()
    t0 = time.time(); rows = []
    for K in KSW:
        rows.append(node(F, Trh, dneff_late, lam_wdm, "none", 0, K, 0.0))
        for spec, x0, eps in itertools.product(SPECTRA, X0S, EPS[1:]):
            rows.append(node(F, Trh, dneff_late, lam_wdm, spec, x0, K, eps))
        print(f"K_sw={K} done {time.time()-t0:.0f}s", flush=True)
    # numerical convergence at representative nodes (step 0.015, 128 bins)
    conv = []
    for spec, x0, eps in (("ir_delta", 10, 0.01), ("ir_delta", 100, 0.1), ("energy_logflat", 10, 0.1), ("energy_logflat", 100, 1.0)):
        a = node(F, Trh, dneff_late, lam_wdm, spec, x0, 0, eps); b = node(F, Trh, dneff_late, lam_wdm, spec, x0, 0, eps, 0.015, 128)
        conv.append(dict(spectrum=spec, x0=x0, eps=eps, dlog10_m=abs(math.log10(a["m_G_eV"] / b["m_G_eV"])),
                         rel_DeltaNeff=abs(a["DeltaNeff_CMB"] / b["DeltaNeff_CMB"] - 1),
                         rel_lambda=abs(a["lambda_em_Mpc"] / b["lambda_em_Mpc"] - 1)))
    # regions
    regions = {}
    for spec, x0, K in itertools.product(SPECTRA, X0S, KSW):
        sub = [r for r in rows if r.get("spectrum") in (spec, "none") and r.get("x0") in (x0, 0) and r["K_sw"] == K]
        def rng(key, flag):
            v = [r[key] for r in sub if r.get(flag)]
            return [min(v), max(v)] if v else None
        regions[f"{spec}|x0={x0}|K_sw={K}"] = dict(
            eps_conservative=rng("eps", "conservative"), eps_strict=rng("eps", "strict"),
            m_G_conservative=rng("m_G_eV", "conservative"), m_G_strict=rng("m_G_eV", "strict"),
            DeltaNeff_CMB_conservative=rng("DeltaNeff_CMB", "conservative"),
            f_em_conservative=rng("f_em", "conservative"),
            failed_by={g: [r["eps"] for r in sub if r.get("P1") and not r.get(g)] for g in ("P2", "P3", "P4", "P5")},
            CMBS4_detectable_and_allowed=[r["eps"] for r in sub if r.get("conservative") and r["CMBS4_detectable"]],
            eps_max_contiguous=_contig(sub))
    k0 = [v for k, v in regions.items() if k.endswith("K_sw=0")]
    closure = any(v["eps_conservative"] for v in k0)
    narrow = closure and any(v["eps_strict"] for v in regions.values())
    verdict = "CIERRE_ESTRECHO" if narrow else ("CIERRE_CONDICIONAL" if closure else "SIN_CIERRE")
    allowed = [r for r in rows if r.get("conservative")]
    preds = dict(m_G_eV=[min(r["m_G_eV"] for r in allowed), max(r["m_G_eV"] for r in allowed)],
                 DeltaNeff_CMB=[min(r["DeltaNeff_CMB"] for r in allowed), max(r["DeltaNeff_CMB"] for r in allowed)],
                 f_em=[min(r["f_em"] for r in allowed), max(r["f_em"] for r in allowed)],
                 eps=[min(r["eps"] for r in allowed), max(r["eps"] for r in allowed)],
                 M_minicluster_Msun=[min(r["M_minicluster_Msun"] for r in allowed), max(r["M_minicluster_Msun"] for r in allowed)],
                 f_a_GeV=fa, LISA_threshold_GeV=LISA_FA, f_a_over_LISA=fa / LISA_FA,
                 superradiance_applies=bool(fa >= SR_FA_MIN), f_a_over_superradiance_min=fa / SR_FA_MIN,
                 superradiance_applies_Stott={k: bool(fa >= v) for k, v in SR_FA_STOTT.items()}) if allowed else None
    em = {k: v["eps_max_contiguous"] for k, v in regions.items()}
    base = [r for r in allowed if r["K_sw"] in (0, 1)]
    preds_K01 = _preds(base)
    req = dict(eps_max_K0=[min(v for k, v in em.items() if k.endswith("K_sw=0")), max(v for k, v in em.items() if k.endswith("K_sw=0"))],
               eps_max_by_config=em)
    out = dict(stage="S5_closure", requirement_on_source=req, predictions_K_sw_0_1=preds_K01, spec="ESPEC_CIERRE.md", verdict=verdict, provenance=provenance(),
               gates=dict(DNEFF_CMB=DNEFF_CMB, DNEFF_BBN=DNEFF_BBN, M_LYA=M_LYA, SR_WINDOWS=SR_WINDOWS, SR_FA_MIN=SR_FA_MIN,
                          F_EM_FREE=F_EM_FREE, lambda_WDM_Mpc=lam_wdm),
               n_nodes=len(rows), n_conservative=len(allowed), n_strict=sum(r.get("strict", False) for r in rows),
               regions=regions, predictions=preds, convergence=conv, runtime_s=time.time() - t0, rows=rows)
    save("s5_closure.json", out)
    print("VERDICT", verdict, f"{len(allowed)}/{len(rows)} conservative, {out['n_strict']} strict, {out['runtime_s']:.0f}s")
    print(json.dumps(preds, indent=1, default=float))
    for c in conv: print("conv", c)
    for k, v in regions.items():
        print(k, "eps_cons", v["eps_conservative"], "eps_strict", v["eps_strict"], "m", v["m_G_conservative"])
    return out


if __name__ == "__main__":
    main()
