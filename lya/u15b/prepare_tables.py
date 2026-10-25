"""U15b (docs/preregistros/ESPEC_U15b.md): transfer-function tables T^2(k; s) = P_model / P_LCDM at z = 0 for the
Lya_abgd_whale likelihood.

  rows        B2, q = 1, R_Gamma = 1.56 (the U12 rows), x0 in {10, 30}, K_sw in {0, 1}; omega_ncdm = s f_em Omega_DM,
              omega_cdm = Omega_DM - omega_ncdm, s in {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3} with s f_em <= 1
  validation  thermal mixtures m = 1, 2 keV with fraction f in {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6} (s = f)

    python -B lya/u15b/prepare_tables.py --workers 40       # WSL, CLASS with the quadrature patch
    python -B lya/u15b/prepare_tables.py --quick             # plumbing check: one validation table, 3 points, k <= 30 h/Mpc
    python -B lya/u15b/prepare_tables.py --u15c --workers 8  # ESPEC_U15c: the two living-branch cases (own manifest_u15c.json)
    python -B lya/u15b/prepare_tables.py --u15d --workers 4  # ESPEC_U15d: the four edge cases of U17f (own manifest_u15d.json)

Output: lya/u15b/tables/<name>.npz with k_h [h/Mpc], s [-], T2[s, k]; and tables/manifest.json."""
import itertools, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "lya")); sys.path.insert(0, str(ROOT / "chain"))
TAB = HERE / "tables"
S_ROWS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0)
F_VAL = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)
KH = np.geomspace(1e-3, 300.0, 600)


def pk0(params, kmax):
    from classy import Class
    p = dict(params); p.update({"z_pk": 0.0, "P_k_max_h/Mpc": kmax, "output": "mPk"})
    c = Class(); c.set(p); c.compute(); h = c.h()
    k = KH[KH <= kmax]
    P = np.array([c.pk_lin(x * h, 0.0) for x in k]) * h ** 3
    c.struct_cleanup(); c.empty(); return P


def job(args):
    name, i, params, kmax = args
    t = time.time()
    try:
        return name, i, pk0(params, kmax), time.time() - t
    except Exception as e:  # noqa: BLE001
        return name, i, None, f"{type(e).__name__}: {e}"


def row_cases():
    import u4a_class as U4, u5_class as U5, s6_measured as S6
    from s4_dark import lambda_wdm_reference
    from ucommon import benchmark
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    dneff = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text())["DeltaNeff_late_homogeneous_max"]; lam, _ = lambda_wdm_reference()
    out = {}
    for x0, K in itertools.product((10, 30), (0, 1)):
        eps = U5.eps_B2(1.56); r = S6.node(F, Trh, dneff, lam, eps, "powerlaw_q1", x0, K)
        tag = f"U15b_B2_R1.56_q1_x{x0}_K{K}"
        pk = S6.source_eps(F, r["m_G_eV"], Trh, x0, eps, spectrum="powerlaw_q1")
        path, T, _ = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], tag)
        s = [x for x in S_ROWS if x * r["f_em"] <= 1.0]
        out[tag] = dict(kind="row", x0=x0, K_sw=K, f_em=r["f_em"], m_G_eV=r["m_G_eV"], s=s,
                        params=[dict(U4.COSMO, omega_cdm=U4.OMEGA_DM) if x == 0 else U4.packet_params(path, T, r["f_em"] * x) for x in s])
    return out


def val_cases(quick=False):
    import u4a_class as U4
    out = {}
    for m in ((1.0,) if quick else (1.0, 2.0)):
        f = (0.0, 0.1, 0.3) if quick else F_VAL
        out[f"U15b_val_thermal_{m:g}keV"] = dict(kind="validation", m_keV=m, s=list(f),
                                                params=[dict(U4.COSMO, omega_cdm=U4.OMEGA_DM) if x == 0 else U4.thermal_params(m, x) for x in f])
    return out


def u15c_cases():
    """ESPEC_U15c, enmienda 3: the living branch at the edges of the U21/U20 band (B2, q = 0.87, x0 = 10, R_Gamma = 1.56, U17d collapse
    population c = 3, sigma_ln = 0.5, R_rho = 1.56, eta = 1): F = 0.5 F0 with r = 1 (primary edge) and F = 0.35 F0 with r = 3 (hot collapse edge)"""
    import u4a_class as U4, u17b_class as U17B
    from ucommon import benchmark
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; out = {}; U17B.SIG_LN = 0.5
    for fr, rd in ((0.5, 1.0), (0.35, 3.0)):
        r, x, n = U17B.row(F0 * fr, Trh, "B2", 3.0, 1.0, R_rho=1.56, r_decay=rd); tag = f"U15c_B2_F{fr:g}_r{rd:g}"
        path, T, _ = U4.psd_from_packets(x, n, tag); f = r["f_nc"]; smax = min(1.25, 1.0 / f)
        s = sorted({0.0, 0.25, 0.5, 0.75, 0.9, 1.0, round(smax, 6)})
        out[tag] = dict(kind="u15c", F_over_F0=fr, r_decay=rd, c=3.0, eta=1.0, R_rho=1.56, f_nc=f, m_G_eV=r["m_G_eV"], K_eff=r["K_eff"], s=s, s_max=smax,
                        params=[dict(U4.COSMO, omega_cdm=U4.OMEGA_DM) if x == 0 else U4.packet_params(path, T, min(f * x, 1.0)) for x in s])
    return out


U15D = (("U15d_F1_Y1.20", 1.0, 1.20), ("U15d_F1_Y1.70", 1.0, 1.70), ("U15d_F0.85_Y1.98", 0.85, 1.98), ("U15d_F0.7_Y1.98", 0.7, 1.98))


def u15d_cases():
    """ESPEC_U15d: the upper edge of U17f with the U19b-measured population (U17e measured_shape 'low', B2, q = 0.87, x0 = 10,
    U17 scaling with Y; sum plateau as U17e/U17f); s as in U15c"""
    import u4a_class as U4, u17e_class as E
    from ucommon import benchmark
    Bm = benchmark(); F0, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]; shape = E.measured_shape("low"); out = {}
    for tag, fr, Y in U15D:
        r, x, n = E.row(F0 * fr, Trh, "B2", shape, tag, Y)
        path, T, _ = U4.psd_from_packets(x, n, tag); f = r["f_nc"]; smax = min(1.25, 1.0 / f)
        s = sorted({0.0, 0.25, 0.5, 0.75, 0.9, 1.0, round(smax, 6)})
        out[tag] = dict(kind="u15d", F_over_F0=fr, Y=Y, X_univ=r["K"] * r["xbar"], f_nc=f, m_G_eV=r["m_G_eV"], K=r["K"], xbar=r["xbar"],
                        T_osc_MeV=r["T_osc_MeV"], s=s, s_max=smax, plateau=os.environ.get("WHALEFALL_PLATEAU"),
                        params=[dict(U4.COSMO, omega_cdm=U4.OMEGA_DM) if x == 0 else U4.packet_params(path, T, min(f * x, 1.0)) for x in s])
    return out


def main(argv):
    quick = "--quick" in argv; workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 2
    kmax = 30.0 if quick else 300.0
    global TAB
    if quick: TAB = HERE / "tables_quick"                       # never mixed with the production tables
    TAB.mkdir(parents=True, exist_ok=True)
    u15c = "--u15c" in argv                                     # ESPEC_U15c: only its two cases, own manifest
    u15d = "--u15d" in argv                                     # ESPEC_U15d: its four cases, own manifest
    cases = val_cases(quick) if quick else (u15c_cases() if u15c else u15d_cases() if u15d else {**val_cases(), **row_cases()})
    jobs = [(n, i, p, kmax) for n, c in cases.items() for i, p in enumerate(c["params"])]
    print(f"{len(jobs)} CLASS runs, {workers} workers", flush=True)
    res = {n: [None] * len(c["params"]) for n, c in cases.items()}
    with ProcessPoolExecutor(workers) as ex:
        for k, (n, i, P, dt) in enumerate(ex.map(job, jobs), 1):
            res[n][i] = P; print(f"  [{k}/{len(jobs)}] {n} s={cases[n]['s'][i]:g} " + (f"ok ({dt / 60:.1f} min)" if P is not None else f"FAILED {dt}"), flush=True)
    manifest = {}
    for n, c in cases.items():
        if any(P is None for P in res[n]): print("skip", n, "(failed runs)"); continue
        P0 = res[n][c["s"].index(0.0)]; T2 = np.array([P / P0 for P in res[n]])
        np.savez(TAB / f"{n}.npz", k_h=KH[KH <= kmax], s=np.array(c["s"]), T2=T2)
        manifest[n] = {k: v for k, v in c.items() if k != "params"}
        manifest[n]["T2_at_s1_k10"] = float(np.interp(math.log(10.0), np.log(KH[KH <= kmax]), T2[c["s"].index(1.0)])) if 1.0 in c["s"] else None
    (TAB / ("manifest_u15c.json" if u15c else "manifest_u15d.json" if u15d else "manifest.json")).write_text(json.dumps(manifest, indent=1, default=float))
    print("tables:", list(manifest), flush=True)


if __name__ == "__main__":
    if "--u15d" in sys.argv: os.environ["WHALEFALL_PLATEAU"] = "sum"          # ESPEC_U15d: the plateau of U17e/U17f, before any chain/ import
    main(sys.argv)
