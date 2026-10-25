"""S6 · physical closure with the MEASURED network efficiency (network/ESPEC_U2.md §5-§6).
eps(ell) per branch: B1 (GHN, xi ~ 0.24 ell): eps = R_Gamma;  B2 (standard scaling, xi* = 1.20): eps = R_Gamma * 1.20/(0.24 ell);
B3 only if the in-window trend of R_Gamma is significant (> 2 sigma): R(ell) = c0 + c1/ell fitted in the window, eps = R(ell)."""
import itertools, json, math
from ucommon import ROOT, RESULTS, save, benchmark, provenance, load
import numpy as np
from scipy.optimize import brentq
import c35b_pseudo_goldstone_dm_rp3 as C35
from s4_dark import gstar, gs, hubble, scale, abundance, at_epoch, emitted_lambda, lambda_wdm_reference
from s5_closure import DNEFF_CMB, DNEFF_BBN, M_LYA, SR_WINDOWS, F_EM_FREE

XI_STAR = 1.20


def source_eps(F, m_eV, Trh, x0, eps_fn, step=0.03, spectrum="ir_delta", spectral_bins=64):
    """= s4_dark.source with efficiency -> eps_fn(ell) evaluated per emission interval"""
    m = m_eV * 1e-9; fa = F / 2; mr = F / 57.735
    Tstop = math.exp(brentq(lambda l: math.log(hubble(math.exp(l)) / (m / 3)), math.log(1e-12), math.log(Trh)))
    edge = np.linspace(math.log(Trh), math.log(Tstop), int(math.ceil(math.log(Trh / Tstop) / step)) + 1)
    hi, lo = edge[:-1], edge[1:]
    T = np.exp((hi + lo) / 2); H = hubble(T); a = scale(T)
    dN = hi - lo + np.log(gs(np.exp(hi)) / gs(np.exp(lo))) / 3
    ell = np.log(mr / H); xi = 0.24 * ell; rho_s = 4 * xi * np.pi * fa ** 2 * ell * H ** 2
    dE = np.clip(eps_fn(ell), 0, None) * 2 * rho_s * dN
    if spectrum == "ir_delta":
        X = np.full((len(T), 1), float(x0)); W = np.ones((len(T), 1))
    else:
        qd, wt = np.polynomial.legendre.leggauss(spectral_bins)
        lo_p = math.log(x0); hi_p = np.log(mr / H)[:, None]
        X = np.exp((lo_p + hi_p) / 2 + qd[None, :] * (hi_p - lo_p) / 2); W = np.broadcast_to(wt / 2, X.shape)
        if spectrum.startswith("powerlaw_q"):        # ESPEC_U5 amendment 1: dE/dln x ~ x^(1-q) on the same support (q = 1: log-flat)
            qq = float(spectrum[len("powerlaw_q"):]); w = (wt / 2)[None, :] * np.exp((1 - qq) * (np.log(X) - math.log(x0)))
            W = w / w.sum(1, keepdims=True)
    p = X * H[:, None]; E = np.sqrt(p * p + m * m); nc = dE[:, None] * W / E * a[:, None] ** 3
    A = np.broadcast_to(a[:, None], X.shape)
    return dict(q=(p * A).ravel(), nc=nc.ravel(), emitted=(dE[:, None] * W).ravel(), a=A.ravel().copy(), T_stop_GeV=Tstop, m_GeV=m,
                F_GeV=F, efficiency=1.0, x0=x0, spectrum=spectrum, ell_range=[float(ell.min()), float(ell.max())],
                eps_range=[float(np.clip(eps_fn(ell), 0, None).min()), float(np.clip(eps_fn(ell), 0, None).max())])


def node(F, Trh, dneff_late, lam_wdm, eps_fn, spec, x0, K):
    def budget(lm):
        mm = 10 ** lm
        return C35.omega(mm, F, K)[0] + abundance(source_eps(F, mm, Trh, x0, eps_fn, spectrum=spec))
    try:
        m = 10 ** brentq(lambda lm: math.log(budget(lm) / 0.12), -27, -10, xtol=1e-7)
    except ValueError:
        return dict(P1=False)
    pk = source_eps(F, m, Trh, x0, eps_fn, spectrum=spec)
    em = abundance(pk); f_em = em / 0.12; lam = emitted_lambda(pk)
    r = dict(m_G_eV=m, Omega_em=em, f_em=f_em, lambda_em_Mpc=lam, lambda_em_over_WDM=lam / lam_wdm,
             DeltaNeff_CMB=at_epoch(pk, 0.26e-9)["DeltaNeff_radiation_moment"] + dneff_late,
             DeltaNeff_BBN=at_epoch(pk, 1e-3)["DeltaNeff_radiation_moment"] + dneff_late,
             T_stop_GeV=pk["T_stop_GeV"], ell_range=pk["ell_range"], eps_range=pk["eps_range"],
             M_minicluster_Msun=C35.minicluster_Msun(m), P1=True)
    r.update(P2=r["DeltaNeff_CMB"] < DNEFF_CMB, P3=r["DeltaNeff_BBN"] < DNEFF_BBN, P4=m >= M_LYA,
             P5=bool(f_em <= F_EM_FREE or lam < lam_wdm), P6_window_hit=any(lo <= m <= hi for lo, hi in SR_WINDOWS))
    r["allowed"] = all(r[k] for k in ("P1", "P2", "P3", "P4", "P5"))
    r["allowed_no_franchise"] = r["allowed"] and bool(lam < lam_wdm or f_em <= 1e-3)
    return r


def main():
    B = benchmark(); F = B["F_GeV"]; Trh = B["T_Gamma_GeV"]
    U2 = json.loads((ROOT / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
    U3B = ROOT / "network" / "results" / "u3" / "U3B_RESULT.json"
    if U3B.exists() and json.loads(U3B.read_text()).get("status") == "GATES_PASSED":      # ESPEC_U3: use the attractor runs
        u3 = json.loads(U3B.read_text()); U2 = dict(U2, R_Gamma=u3["R_Gamma"], U1_control_gate=dict(passed=True, source="U3b"))
    U4B = ROOT / "network" / "results" / "u4" / "U4B_RESULT.json"                         # ESPEC_U4: larger lattices
    if U4B.exists():
        ok = [b for b in json.loads(U4B.read_text()).get("blocks", []) if b.get("status") == "GATES_PASSED" and "R_Gamma_summary" in b]
        if ok:
            b = max(ok, key=lambda x: x["N"]); rs = b["R_Gamma_summary"]
            U2 = dict(U2, R_Gamma=dict(value=rs["value"], total=rs["total"], source=f"U4b N={b['N']}"),
                      U1_control_gate=dict(passed=True, source=f"U4b N={b['N']}"))
    s3 = json.loads((RESULTS / "s3_sectors.json").read_text())
    dneff_late = s3["DeltaNeff_late_homogeneous_max"]; lam_wdm, _ = lambda_wdm_reference()
    Rv, Rs = U2["R_Gamma"]["value"], U2["R_Gamma"]["total"]
    Rset = {"-1s": max(Rv - Rs, 1e-6), "central": Rv, "+1s": Rv + Rs}
    branches = {"B1_GHN": lambda R: (lambda ell: np.full_like(ell, R)),
                "B2_standard_scaling": lambda R: (lambda ell: R * XI_STAR / (0.24 * ell))}
    if U2["R_Gamma_trend"]["significant"]:
        g = np.array(U2["R_Gamma_A"]["grid"]); pe = np.array(U2["R_Gamma_A"]["per_ell"]); mk = np.isfinite(pe)
        c1, c0 = np.polyfit(1 / g[mk], pe[mk], 1)
        branches["B3_trend_1_over_ell"] = lambda R, c0=c0, c1=c1, Rv=Rv: (lambda ell: (c0 + c1 / ell) * (R / Rv))
    rows = []
    for (bn, bf), (rn, R), spec, x0, K in itertools.product(branches.items(), Rset.items(), ("ir_delta", "energy_logflat"), (10, 30), (0, 1)):
        r = node(F, Trh, dneff_late, lam_wdm, bf(R), spec, x0, K)
        r.update(branch=bn, R_label=rn, R=R, spectrum=spec, x0=x0, K_sw=K); rows.append(r)
        print(bn, rn, f"R={R:.3f}", spec, x0, K, "allowed" if r.get("allowed") else "EXCLUDED",
              {k: (round(r[k], 4) if isinstance(r.get(k), float) else r.get(k)) for k in ("m_G_eV", "f_em", "lambda_em_over_WDM", "DeltaNeff_CMB")}, flush=True)

    def branch_ok(bn, rl, K=0): return any(r.get("allowed") for r in rows if r["branch"] == bn and r["R_label"] == rl and r["K_sw"] == K)
    names = list(branches)
    ok_plus = {b: branch_ok(b, "+1s") for b in names}; ok_minus = {b: branch_ok(b, "-1s") for b in names}; ok_c = {b: branch_ok(b, "central") for b in names}
    if all(ok_plus.values()): verdict = "CIERRE_FISICO"
    elif not any(ok_minus.values()): verdict = "FALSADA_EN_FORMA_ACTUAL"
    elif any(ok_plus.values()): verdict = "CIERRE_SEGUN_RAMA"
    else: verdict = "INCONCLUSO"
    if not U2["U1_control_gate"]["passed"]: verdict = "INCONCLUSO (control U(1) fuera de rango)"
    # EXPLORATORY (not preregistered): massless-Goldstone variant (no explicit breaking, no walls, network never collapses):
    # the emitted Goldstones are dark radiation; DeltaNeff at BBN and CMB per branch (m = 1e-29 eV keeps every packet relativistic)
    expl = {}
    for (bn, bf), (rn, R) in itertools.product(branches.items(), Rset.items()):
        for spec in ("ir_delta", "energy_logflat"):
            pk = source_eps(F, 1e-29, Trh, 10, bf(R), spectrum=spec)
            expl[f"{bn}|{rn}|{spec}"] = dict(DeltaNeff_BBN=at_epoch(pk, 1e-3)["DeltaNeff_radiation_moment"],
                                             DeltaNeff_CMB=at_epoch(pk, 0.26e-9)["DeltaNeff_radiation_moment"])
    out_expl = dict(note="exploratory, not preregistered: massless Goldstones (no m_G), DM from another source", DeltaNeff=expl)
    # EXPLORATORY 2: mixed cold+warm Lyman-alpha gate (Garcia-Gallego et al. arXiv:2504.06367; audit correction of the attribution, 2 sigma: f_WDM < 0.14 (m/keV)^1.1, m in 1-4 keV)
    # with m_eq = thermal-relic mass giving the same mass-weighted lambda_fs (same integrator); replaces the fixed 10 % franchise of P5
    from s4_dark import lambda_wdm_reference as lw
    mgrid = np.geomspace(0.3, 20, 60); lgrid = np.array([lw(mk)[0] for mk in mgrid])
    def m_eq(lam): return float(np.exp(np.interp(-math.log(lam), -np.log(lgrid), np.log(mgrid)))) if lam > 0 else 1e9
    for r in rows:
        if not r.get("P1"): continue
        me = m_eq(r["lambda_em_Mpc"]); fmax = 1.0 if me >= 5.3 else min(1.0, 0.14 * me ** 1.1)
        r["m_eq_keV"] = me; r["f_max_mixed"] = fmax; r["P5_mixed"] = bool(r["f_em"] <= fmax)
        r["allowed_mixed_gate"] = all(r[k] for k in ("P1", "P2", "P3", "P4")) and r["P5_mixed"]
    out_expl["mixed_gate_allowed"] = {f"{b}|{rl}": any(r.get("allowed_mixed_gate") for r in rows if r["branch"] == b and r["R_label"] == rl and r["K_sw"] == 0)
                                      for b in names for rl in Rset}
    # EXPLORATORY 3: suppression factor mu (eps -> mu eps) needed to pass the preregistered gates, central R, K_sw = 0
    supp = {}
    for (bn, bf), spec, x0 in itertools.product(branches.items(), ("ir_delta", "energy_logflat"), (10, 30)):
        lo, hi = -4.0, 0.0
        if node(F, Trh, dneff_late, lam_wdm, bf(Rv), spec, x0, 0).get("allowed"): supp[f"{bn}|{spec}|x0={x0}"] = 1.0; continue
        for _ in range(14):
            mid = 0.5 * (lo + hi)
            ok = node(F, Trh, dneff_late, lam_wdm, (lambda e, f=bf(Rv), mu=10 ** mid: mu * f(e)), spec, x0, 0).get("allowed")
            lo, hi = (mid, hi) if ok else (lo, mid)
        supp[f"{bn}|{spec}|x0={x0}"] = 10 ** lo
    out_expl["suppression_needed_central_K0"] = supp
    out = dict(stage="S6_measured_closure", exploratory_massless_variant=out_expl, spec="docs/preregistros/ESPEC_U2.md", verdict=verdict, R_Gamma=U2["R_Gamma"],
               allowed_at_plus1sigma=ok_plus, allowed_at_central=ok_c, allowed_at_minus1sigma=ok_minus,
               U1_control_gate=U2["U1_control_gate"], provenance=provenance(), rows=rows)
    save("s6_measured_closure.json", out)
    print("VERDICT", verdict, ok_minus, ok_c, ok_plus)


if __name__ == "__main__":
    main()
