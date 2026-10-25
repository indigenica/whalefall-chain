"""U17b (docs/preregistros/ESPEC_U17b.md): the string-wall collapse population of U17 through CLASS, with the U12 gates.
Rows: B2, B4a, B1 at q = 0.87, x0 = 10, R_Gamma = 1.56; (c, eta) in {(3, 1), (30, 1), (3, 0.1)}.  m_G solved as in U17
(chain/u17_ksw.py); the non-cold species = emission packets (S6) + collapse packets (log-normal in q, sigma_ln = 0.5, mean
momentum c m_G at H = m_G/3, number m_G n = K_eff Omega_mis_bare rho_c); bare misalignment stays cold.

    python -B lya/u17b_class.py --workers 2      ->  lya/results/u17b/U17B_RESULT.json   (CLASS cache shared with U4a/U5/U12)"""
import json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
import u4a_class as U4                                                        # noqa: E402
import u5_class as U5                                                         # noqa: E402

OUT = HERE / "results" / "u17b"
BRANCHES = ("B2", "B4a", "B1"); Q, X0 = 0.87, 10
VARIANTS = ((3.0, 1.0), (30.0, 1.0), (3.0, 0.1))     # (c, eta)
SIG_LN = 0.5


def c_energy(c, sig):
    """mean energy / m of the log-normal collapse packet whose mean momentum is c m"""
    zh, wh = np.polynomial.hermite_e.hermegauss(48); wh = wh / wh.sum()
    return float(np.sum(wh * np.sqrt(1 + (c * np.exp(sig * zh - sig ** 2 / 2)) ** 2)))


def row(F, Trh, branch, c, eta, legacy_energy=False, R_rho=1.0, r_decay=1.0):
    """c = mean MOMENTUM / m_G of the collapse packet.  The number is the energy eta*rho_net divided by the packet's mean
    energy c_E m_G (codex audit 14 Oct 2026).  legacy_energy=True reproduces U17b/U17c as run (number from c m_G, i.e. energy
    and momentum mixed); kept only so the stored results can be regenerated.  R_rho: RP3/U(1) network energy ratio;
    r_decay = a_decay/a_osc: momentum c m_G at a_decay, energy frozen at H = m_G/3 (ESPEC_U17d)."""
    import c35b_pseudo_goldstone_dm_rp3 as C35
    import u17_ksw as U17
    from s4_dark import abundance, scale, RHOC_H2
    from s6_measured import source_eps
    cE = c if legacy_energy else c_energy(c, SIG_LN)
    r = U17.solve(F, Trh, branch, Q, X0, c_eta=cE / eta, R_rho=R_rho)
    m_eV = r["m_G_eV"]; m = m_eV * 1e-9; K = r["K"]
    om_mis_bare = C35.omega(m_eV, F, 0.0)[0]; om_coll = om_mis_bare * K
    pk = source_eps(F, m_eV, Trh, X0, U17.eps_fn(branch), spectrum=f"powerlaw_q{Q:g}")
    om_em = abundance(pk)
    a_osc = float(scale(r["T_osc_MeV"] * 1e-3))
    zh, wh = np.polynomial.hermite_e.hermegauss(48); wh = wh / wh.sum()
    q_coll = c * m * a_osc * r_decay * np.exp(SIG_LN * zh - SIG_LN ** 2 / 2)    # comoving momentum [GeV], mean = c m a_decay
    n_coll = om_coll * RHOC_H2 / m * wh                                          # comoving number with m n = Omega_coll rho_c
    q_all = np.concatenate([pk["q"], q_coll]); n_all = np.concatenate([pk["nc"], n_coll])
    f_nc = (om_em + om_coll) / 0.12
    return dict(branch=branch, q=Q, x0=X0, c=c, eta=eta, c_over_eta=c / eta, c_E=cE, R_rho=R_rho, r_decay=r_decay,
                legacy_energy=legacy_energy, K_eff=K, m_G_eV=m_eV, T_osc_MeV=r["T_osc_MeV"],
                Omega_mis_bare=om_mis_bare, Omega_collapse=om_coll, Omega_em=om_em, f_nc=f_nc, f_em=om_em / 0.12,
                M_minicluster_Msun=C35.minicluster_Msun(m_eV), below_ULDM_bound=bool(m_eV < 2e-20),
                tag=f"U17b_{branch}_c{c:g}_eta{eta:g}"), q_all / m, n_all


def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    import classy
    from ucommon import benchmark
    out = dict(spec="docs/preregistros/ESPEC_U17b.md", classy_version=getattr(classy, "__version__", "?"), q=Q, x0=X0, R_Gamma=1.56,
               variants=[dict(c=c, eta=e) for c, e in VARIANTS], sigma_ln=SIG_LN)
    P0 = U4.p_cdm()
    # C0b: the quadrature patch must be present (as U12)
    qg = np.geomspace(1e-3, 40, 800); U4.PSD.mkdir(parents=True, exist_ok=True); fd = U4.PSD / "C0_fermi_dirac.dat"
    np.savetxt(fd, np.column_stack([qg, 1 / (np.exp(qg) + 1)]), fmt="%.10e")
    pa = U4.thermal_params(3.5, 1.0); pq = dict(pa, use_ncdm_psd_files=1, ncdm_psd_filenames=str(fd), ncdm_quadrature_strategy=4, ncdm_N_momentum_bins=20)
    Pa, Pq = U4.pmap([pa, pq], min(workers, 2), ["C0 analytic", "C0b"]); s20 = (U4.KH >= 0.5) & (U4.KH <= 20)
    c0b = abs(U4.deltaA(Pq / P0, P0) / U4.deltaA(Pa / P0, P0) - 1) < 5e-3 and float(np.max(np.abs(Pq[s20] / Pa[s20] - 1))) < 5e-3
    out["C0b_passed"] = bool(c0b)
    if not c0b: raise RuntimeError("C0b failed: CLASS without the quadrature patch")
    # references and the extended shape family (as U5/U12)
    refs = {"thermal_5.3keV": U4.thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": U4.thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    fam_m, fam_f = U5.FAM_M_EXT, U5.FAM_F_EXT
    Ps = U4.pmap(list(refs.values()) + [U4.thermal_params(m, f) for m in fam_m for f in fam_f], workers,
                 list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f], "thermal (cache)")
    if any(P is None for P in Ps): raise RuntimeError("CLASS failed on a thermal run: see lya/results/u4a/class_errors.log")
    ref_dA = {k: U4.deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}; strict, wide = min(ref_dA.values()), max(ref_dA.values())
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    out.update(references_deltaA=ref_dA, deltaA_strict=strict, deltaA_wide=wide); print("refs", ref_dA, "strict", strict, flush=True)
    # rows
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    rows, params = [], []
    for b in BRANCHES:
        for c, e in VARIANTS:
            r, x, n = row(F, Trh, b, c, e, legacy_energy=True)      # as run on 14 Oct (see row())
            path, T, tail = U4.psd_from_packets(x, n, r["tag"]); r["tail_number_fraction"] = tail
            params.append(U4.packet_params(path, T, min(r["f_nc"], 1.0))); rows.append(r)
            print(f"{r['tag']}: K={r['K_eff']:.1f} m_G={r['m_G_eV']:.2e} T_osc={r['T_osc_MeV']:.3f} MeV f_nc={r['f_nc']:.4f}", flush=True)
    for r, P in zip(rows, U4.pmap(params, workers, [r["tag"] for r in rows], "row")):
        if P is None: r.update(class_error=True, LYA_primary=False, LYA_dA_strict=False); continue
        r.update(U5.gates(P / P0, P0, fam, strict, wide))
        print(f"{r['tag']:<22} dA={r['deltaA']:.3f} fit m={r['shape']['m_keV']:.2f} f={r['shape']['f']:.3f} primary={r['LYA_primary']}", flush=True)
    out["rows"] = rows
    st = {}
    for b in BRANCHES:
        ok = [r["LYA_primary"] for r in rows if r["branch"] == b]
        st[b] = "PASA" if all(ok) else "FALLA" if not any(ok) else "PARCIAL"
    eta1_fail = lambda b: not any(r["LYA_primary"] for r in rows if r["branch"] == b and r["eta"] == 1.0)
    if st["B2"] == "PASA" and eta1_fail("B4a") and eta1_fail("B1"): v = "SOLO_XI_CONSTANTE_VIVE"
    elif all(s == "PASA" for s in st.values()): v = "TODAS_VIVEN"
    elif st["B2"] != "PASA": v = "XI_CONSTANTE_TAMBIEN_EN_PELIGRO"
    else: v = "DEPENDE_DE_ETA"
    out.update(status=st, verdict_U17b=v, runtime_h=(time.time() - t0) / 3600)
    (OUT / "U17B_RESULT.json").write_text(json.dumps(out, indent=1, default=float, allow_nan=False))
    print("U17b:", v, st, flush=True)


if __name__ == "__main__":
    main(sys.argv)
