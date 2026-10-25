"""U4a (docs/preregistros/ESPEC_U4.md): Lyman-alpha test of the S6 closure rows with CLASS (classy).

    python lya/u4a_class.py                 # full run (WSL, 32 cores: ~2 h of rows + C1, ~6-12 h on two cores)
    python lya/u4a_class.py --workers 4     # number of parallel CLASS runs (default: all cores)
    python lya/u4a_class.py --quick         # small grids + 10 central rows, no C1/C2 (pipeline check, ~1 h with 32 cores)

Steps: C0 (psd-file reading) -> thermal references + mixed thermal family (CLASS) -> rows (packets -> psd file -> CLASS)
-> delta A and shape-matching gates -> C1 (precision) and C2 (vs lya/boltz.py) -> verdict.  Every CLASS P(k) is cached in
lya/results/u4a/cache/ (the job can be restarted).  Output: lya/results/u4a/U4A_RESULT.json"""
import hashlib, json, math, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
OUT = HERE / "results" / "u4a"; CACHE = OUT / "cache"; PSD = OUT / "psd"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))

H = 0.674; OMEGA_DM = 0.120; Z = 3.0; KMAX = 1000.0
COSMO = {"h": H, "omega_b": 0.0224, "n_s": 0.9649, "A_s": 2.1e-9, "N_ur": 3.044, "output": "mPk",
         "P_k_max_h/Mpc": KMAX, "z_pk": Z}
NCDM_PREC = {"ncdm_fluid_approximation": 3, "l_max_ncdm": 30, "tol_ncdm_synchronous": 1e-4, "tol_ncdm_bg": 1e-6}
NCDM_PREC_HI = {"ncdm_fluid_approximation": 3, "l_max_ncdm": 50, "tol_ncdm_synchronous": 1e-5, "tol_ncdm_bg": 1e-7}
# ncdm_fluid_approximation = 3 is CLASS's 'ncdmfa_none' (enum mb, hu, CLASS, none): full hierarchy, no fluid switch.  It is
# required here: the fluid approximation (2) takes the sound speed from <v^2> of the whole species, which for these heavy
# tails (dN/dln q ~ 1/q) grows with the largest sampled q, so delta A does not converge (0.22 -> 0.29 -> 0.36 for 20/40/80 nodes).
# Row quadrature (ESPEC_U4 amendment 1).  CLASS's adaptive sampling (strategy 0) gives 125-230 momentum nodes for these broad
# distributions (a steep edge at low q and a dN/dln q ~ 1/q tail over three decades); with the full hierarchy to
# k = 1000 h/Mpc that is days per row.  CLASS's fixed rules (Laguerre, trapezoid in q or in 1/(1+q)) cannot place nodes
# uniformly in ln q and converge erratically (+-40 % in delta A at 20-40 nodes).  Rows therefore use strategy 4, a trapezoid
# uniform in ln q over the range of the psd file, added by env/class_qm_trapz_log.patch (18 lines, quadrature only): 45 nodes
# are converged to < 0.1 % in delta A in the lya/boltz.py tests (45 vs 90 vs 600 nodes).  C0b checks the patch on a
# Fermi-Dirac file; C1 doubles the nodes.
ROW_NQ = 45; ROW_NQ_HI = 90; QSCALE = 1.0
ROW_QUAD = {"ncdm_quadrature_strategy": 4, "ncdm_N_momentum_bins": ROW_NQ}
ROW_QUAD_HI = {"ncdm_quadrature_strategy": 4, "ncdm_N_momentum_bins": ROW_NQ_HI}
KB_TCMB = 8.617333262e-5 * 2.7255               # eV
M_SURR = 1000.0                                  # surrogate mass [eV] (scale invariance: only q/m and Omega matter)
KH = np.geomspace(0.05, KMAX, 120)               # h/Mpc
FAM_M = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0]          # keV
FAM_F = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.75, 1.0]
GG_TABLE = {1.0: 0.16, 2.0: 0.35, 3.0: 0.50, 4.0: 0.67}           # Garcia-Gallego et al. 2025, tabulated (wide)


# ------------------------------------------------------------------ CLASS wrappers
def _key(params):
    """cache key: parameters + classy version (so a cache made by another CLASS build, or by the test mock, is never reused)"""
    import classy
    tag = dict(params, _classy=getattr(classy, "__version__", "?"), _classy_file=getattr(classy, "__file__", "?"))
    if params.get("ncdm_psd_filenames"):          # the file's content, not only its name, enters the key
        tag["_psd_sha"] = hashlib.sha256(Path(params["ncdm_psd_filenames"]).read_bytes()).hexdigest()
    return hashlib.sha256(json.dumps(tag, sort_keys=True, default=str).encode()).hexdigest()[:20]


def class_pk(params):
    """P(k) [Mpc^3] at z = 3 on KH (h/Mpc); cached by parameter hash."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{_key(params)}.npz"
    if f.exists():
        try:
            return np.load(f)["P"]
        except Exception:  # noqa: BLE001  (corrupt cache from a killed run: recompute)
            f.unlink()
    from classy import Class
    c = Class(); c.set(params); t = time.time()
    try:
        c.compute()
    except Exception as e:  # noqa: BLE001
        (OUT / "class_errors.log").open("a").write(f"{f.stem} {type(e).__name__}: {e}\n{json.dumps(params, default=str)}\n")
        c.struct_cleanup(); c.empty(); return None
    P = np.array([c.pk(k * H, Z) for k in KH])
    meta = dict(seconds=time.time() - t)
    try:
        bg = c.get_background(); meta["Omega_ncdm_0"] = float(bg["(.)rho_ncdm[0]"][-1] / bg["(.)rho_crit"][-1]) if "(.)rho_ncdm[0]" in bg else None
    except Exception:  # noqa: BLE001
        pass
    c.struct_cleanup(); c.empty()
    tmp = f.with_name(f.stem + f".tmp{os.getpid()}.npz")
    np.savez(tmp, P=P, params=json.dumps(params, default=str), meta=json.dumps(meta)); os.replace(tmp, f)   # atomic
    return P


def p_cdm():
    return class_pk(dict(COSMO, omega_cdm=OMEGA_DM))


def thermal_params(m_keV, f_w, prec=NCDM_PREC):
    """thermal fermion WDM (Viel et al. convention Omega h^2 = (T_x/T_nu)^3 m / 94 eV), fraction f_w of the DM"""
    om = OMEGA_DM * f_w; Tnu_over_Tcmb = (4 / 11) ** (1 / 3)
    T = (om * 94.0 / (m_keV * 1e3)) ** (1 / 3) * Tnu_over_Tcmb
    p = dict(COSMO, omega_cdm=OMEGA_DM - om, N_ncdm=1, m_ncdm=m_keV * 1e3, T_ncdm=T, omega_ncdm=om, **prec)
    return p


def psd_from_packets(x, n, tag, npts=400, bw=0.15, qtail=0.999, refine=1, qscale=QSCALE):
    """KDE in ln(q/m) of number-weighted packets -> file 'q f0(q)' with q in units of T_ncdm = qscale * (number-weighted
    geometric median); returns (path, T_ncdm, tail)."""
    import boltz as B
    sp, tail = B.from_packets(x, n, 1.0, nq=npts * refine, bw=bw, qtail=qtail, lmax=4)
    xs, G = sp.x, sp.W / (np.log(sp.x[1]) - np.log(sp.x[0]))          # G = dN/dln x
    xmed = float(np.exp(np.sum(np.log(xs) * G) / np.sum(G))) * qscale
    q = xs / xmed; f0 = G / q ** 3; f0 = f0 / f0.max()
    PSD.mkdir(parents=True, exist_ok=True); path = PSD / f"{tag}.dat"
    np.savetxt(path, np.column_stack([q, f0]), fmt="%.10e")
    T = xmed * M_SURR / KB_TCMB
    return str(path), T, tail


def packet_params(path, T, f_em, prec=NCDM_PREC, quad=None):
    om = OMEGA_DM * f_em
    return dict(COSMO, omega_cdm=OMEGA_DM - om, N_ncdm=1, use_ncdm_psd_files=1, ncdm_psd_filenames=path,
                m_ncdm=M_SURR, T_ncdm=T, omega_ncdm=om, **prec, **(ROW_QUAD if quad is None else quad))


# ------------------------------------------------------------------ estimators
def deltaA(T2, P0):
    import boltz as B
    return B.deltaA(KH, T2, P0)


def f_max_strict(m):
    if m >= 5.3: return 1.0
    if m <= 4.0: return min(1.0, 0.14 * m ** 1.1)
    a = 0.14 * 4.0 ** 1.1
    return float(np.exp(np.interp(math.log(m), [math.log(4.0), math.log(5.3)], [math.log(a), 0.0])))


def f_max_wide(m):
    if m >= 5.3: return 1.0
    ms = sorted(GG_TABLE); lm = [math.log(v) for v in ms] + [math.log(5.3)]; lf = [math.log(GG_TABLE[v]) for v in ms] + [0.0]
    if m < 1.0: return 0.16 * m ** 1.1
    return float(np.exp(np.interp(math.log(m), lm, lf)))


def fit_shape(T2, fam):
    """least-squares in ln T^2 over 0.5-20 h/Mpc against the mixed-thermal family (bilinear in (ln m, f) + Nelder-Mead)."""
    from scipy.interpolate import RegularGridInterpolator
    from scipy.optimize import minimize
    sel = (KH >= 0.5) & (KH <= 20.0); y = np.log(np.maximum(T2[sel], 1e-30))
    lm = np.log(np.array(fam["m"])); ff = np.array(fam["f"]); G = np.log(np.maximum(np.array(fam["T2"])[:, :, sel], 1e-30))
    itp = RegularGridInterpolator((lm, ff), G, bounds_error=False, fill_value=None)
    cost = lambda p: float(np.sqrt(np.mean((itp([[p[0], p[1]]])[0] - y) ** 2)))
    i, j = np.unravel_index(np.argmin([[np.sqrt(np.mean((G[a, b] - y) ** 2)) for b in range(len(ff))] for a in range(len(lm))]), G.shape[:2])
    res = minimize(lambda p: cost([np.clip(p[0], lm[0], lm[-1]), np.clip(p[1], ff[0], ff[-1])]), [lm[i], ff[j]], method="Nelder-Mead",
                   options=dict(xatol=1e-4, fatol=1e-6, maxiter=400))
    m = float(np.exp(np.clip(res.x[0], lm[0], lm[-1]))); f = float(np.clip(res.x[1], ff[0], ff[-1]))
    return dict(m_keV=m, f=f, rms=float(res.fun), grid_best=dict(m_keV=float(np.exp(lm[i])), f=float(ff[j])))


# ------------------------------------------------------------------ jobs (parallel)
def _job(args):
    kind, payload = args
    os.environ["OMP_NUM_THREADS"] = "1"   # one CLASS thread per worker (workers = cores)
    if kind == "params":
        return class_pk(payload)
    raise ValueError(kind)


def pmap(params_list, workers, labels=None, what="CLASS"):
    """parallel CLASS runs; prints one line per finished run (cached runs return at once)."""
    n = len(params_list); labels = labels or [f"#{i}" for i in range(n)]; t0 = time.time(); out = [None] * n
    def report(i, done):
        print(f"   [{done}/{n}] {what} {labels[i]} {'ERROR' if out[i] is None else 'ok'}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
    if workers <= 1 or n == 1:
        for i, p in enumerate(params_list):
            out[i] = class_pk(p); report(i, i + 1)
        return out
    from multiprocessing import get_context
    with get_context("spawn").Pool(min(workers, n)) as pool:
        for done, (i, P) in enumerate(pool.imap_unordered(_job_i, list(enumerate(params_list))), 1):
            out[i] = P; report(i, done)
    return out


def _job_i(args):
    i, p = args
    return i, _job(("params", p))


def default_workers():
    """min(cores, available RAM / 1.2 GB per CLASS run, 32)"""
    avail = None
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"): avail = int(line.split()[1]) / 1e6
    except OSError:
        pass
    return max(1, min(os.cpu_count() or 2, int((avail or 8) / 1.2), 32))


# ------------------------------------------------------------------ main
def main(argv):
    workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else default_workers()
    quick = "--quick" in argv
    fam_m, fam_f = ([1.0, 3.0, 6.0], [0.1, 0.5, 1.0]) if quick else (FAM_M, FAM_F)
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    import classy
    out = dict(spec="docs/preregistros/ESPEC_U4.md (U4a)", classy_version=getattr(classy, "__version__", "?"), quick=quick, z=Z, kh=KH.tolist())
    print("classy", out["classy_version"], "workers", workers, flush=True)
    P0 = p_cdm()

    # C0: psd-file reading (Fermi-Dirac table vs analytic thermal ncdm, 3.5 keV)
    q = np.geomspace(1e-3, 40, 800); PSD.mkdir(parents=True, exist_ok=True); fd = PSD / "C0_fermi_dirac.dat"
    np.savetxt(fd, np.column_stack([q, 1 / (np.exp(q) + 1)]), fmt="%.10e")
    pa = thermal_params(3.5, 1.0); pf = dict(pa, use_ncdm_psd_files=1, ncdm_psd_filenames=str(fd))
    Pa, Pf = pmap([pa, pf], min(workers, 2), ["C0 analytic", "C0 file"]); sel = (KH >= 0.5) & (KH <= 50)
    if Pa is None or Pf is None:
        raise RuntimeError("CLASS failed in C0: see lya/results/u4a/class_errors.log")
    out["C0"] = dict(max_rel=float(np.max(np.abs(Pf[sel] / Pa[sel] - 1)))); out["C0"]["passed"] = out["C0"]["max_rel"] < 5e-3
    print("C0", out["C0"], flush=True)

    # references and mixed thermal family
    refs = {"thermal_5.3keV": thermal_params(5.3, 1.0)}
    refs.update({f"mixed_{m:g}keV_f{0.14 * m ** 1.1:.3f}": thermal_params(m, 0.14 * m ** 1.1) for m in (1.0, 2.0, 3.0, 4.0)})
    refs.update({f"table_{m:g}keV_f{f:.2f}": thermal_params(m, f) for m, f in GG_TABLE.items()})
    fam_list = [thermal_params(m, f) for m in fam_m for f in fam_f]
    # C0b (amendment 1): the row quadrature (patched strategy 4, 20 nodes) on the C0 Fermi-Dirac file must reproduce the analytic
    # run: |d deltaA| / deltaA < 0.5 % and T^2 within 0.5 % on 0.5-20 h/Mpc (the shape-gate range).  An unpatched CLASS fails.
    pq = dict(pf, ncdm_quadrature_strategy=4, ncdm_N_momentum_bins=20)
    Ps = pmap(list(refs.values()) + fam_list + [pq], workers, list(refs) + [f"family {m:g} keV f={f:g}" for m in fam_m for f in fam_f] + ["C0b"])
    Pq = Ps.pop(); s20 = (KH >= 0.5) & (KH <= 20)
    if Pq is None:
        out["C0b"] = dict(passed=False, error=True)
    else:
        dA_a, dA_q = deltaA(Pa / P0, P0), deltaA(Pq / P0, P0)
        out["C0b"] = dict(nodes=20, rel_deltaA=abs(dA_q / dA_a - 1), max_rel_T2_0p5_20=float(np.max(np.abs(Pq[s20] / Pa[s20] - 1))),
                          max_rel_T2_0p5_50=float(np.max(np.abs(Pq[sel] / Pa[sel] - 1))))
        out["C0b"]["passed"] = bool(out["C0b"]["rel_deltaA"] < 5e-3 and out["C0b"]["max_rel_T2_0p5_20"] < 5e-3)
    print("C0b", out["C0b"], flush=True)
    if not out["C0b"]["passed"] and out["classy_version"] != "MOCK":     # (the plumbing mock is not accurate)
        raise RuntimeError("C0b failed: CLASS without the quadrature patch? see env/INSTALAR_WSL.md section 3")
    if any(P is None for P in Ps):
        raise RuntimeError("CLASS failed on a thermal reference or on the mixed family: see lya/results/u4a/class_errors.log")
    ref_dA = {k: deltaA(P / P0, P0) for k, P in zip(refs, Ps[:len(refs)])}
    fam = dict(m=fam_m, f=fam_f, T2=np.array([P / P0 for P in Ps[len(refs):]]).reshape(len(fam_m), len(fam_f), -1).tolist())
    strict = min(v for k, v in ref_dA.items() if not k.startswith("table")); wide = max(v for k, v in ref_dA.items() if not k.startswith("table"))
    out["deltaA_table_max"] = max(v for k, v in ref_dA.items() if k.startswith("table"))   # informative (tabulated Garcia-Gallego fractions)
    out.update(references_deltaA=ref_dA, deltaA_strict=strict, deltaA_wide=wide)
    out["viel_check"] = {f"{m:g}keV": dict(deltaA_CLASS=deltaA(P / P0, P0)) for m, P in
                         zip((2.0, 3.5, 5.3), pmap([thermal_params(m, 1.0) for m in (2.0, 3.5, 5.3)], workers, ["Viel 2 keV", "Viel 3.5 keV", "Viel 5.3 keV"]))}
    print("refs", ref_dA, "strict", strict, "wide", wide, flush=True)

    # rows
    import s6_measured as S6
    from run_u3a import eps_fn_for
    from ucommon import benchmark, RESULTS
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    s6 = json.loads((RESULTS / "s6_measured_closure.json").read_text())
    U2 = json.loads((ROOT / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
    rows = [r for r in s6["rows"] if r.get("P1") and r["K_sw"] in (0, 1) and all(r.get(k) for k in ("P2", "P3", "P4"))]
    if quick:
        rows = [r for r in rows if r["K_sw"] == 0 and r["R_label"] == "central"]
    meta, params = [], []
    for r in rows:
        tag = f"{r['branch'][:2]}_{r['R_label']}_{r['spectrum'][:6]}_x{r['x0']}_K{r['K_sw']}"
        if r["f_em"] < 1e-4:
            meta.append((r, tag, None, 0.0)); continue
        pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
        path, T, tail = psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], tag)
        meta.append((r, tag, len(params), tail)); params.append(packet_params(path, T, r["f_em"]))
    print(f"{len(params)} CLASS row runs ({ROW_NQ} momentum nodes in ln q, full hierarchy; ~20-40 min each, one core each)", flush=True)
    Prows = pmap(params, workers, [m[1] for m in meta if m[2] is not None], "row")
    res = []
    for r, tag, idx, tail in meta:
        if idx is not None and Prows[idx] is None:      # CLASS failed for this row: recorded, never counted as allowed
            res.append(dict({k: r[k] for k in ("branch", "R_label", "R", "spectrum", "x0", "K_sw", "m_G_eV", "f_em")}, tag=tag, class_error=True,
                            LYA_dA_strict=False, LYA_dA_wide=False, LYA_shape_strict=False, LYA_shape_wide=False, LYA_primary=False, deltaA=float("nan")))
            print(f"{tag:<34} CLASS ERROR (see lya/results/u4a/class_errors.log)", flush=True); continue
        T2 = np.ones_like(KH) if idx is None else Prows[idx] / P0
        dA = deltaA(T2, P0); shp = fit_shape(T2, fam) if idx is not None else dict(m_keV=float("inf"), f=0.0, rms=0.0)
        rr = {k: r[k] for k in ("branch", "R_label", "R", "spectrum", "x0", "K_sw", "m_G_eV", "f_em", "DeltaNeff_CMB")}
        rr.update(tag=tag, tail_number_fraction=tail, deltaA=dA, shape=shp,
                  LYA_dA_strict=bool(dA <= strict), LYA_dA_wide=bool(dA <= wide),
                  LYA_shape_strict=bool(shp["f"] <= f_max_strict(shp["m_keV"])), LYA_shape_wide=bool(shp["f"] <= f_max_wide(shp["m_keV"])),
                  shape_representable=bool(shp["rms"] <= 0.1))
        rr["LYA_primary"] = rr["LYA_shape_strict"] if rr["shape_representable"] else rr["LYA_dA_strict"]
        res.append(rr)
        print(f"{tag:<34} f_em={r['f_em']:.3f} dA={dA:.3f} fit m={shp['m_keV']:.2f} keV f={shp['f']:.2f} rms={shp['rms']:.3f} "
              f"shape_strict={rr['LYA_shape_strict']} dA_strict={rr['LYA_dA_strict']}", flush=True)
    out["rows"] = res

    # C1 precision and C2 boltz comparison on the decisive rows (K_sw = 0: lowest delta A and closest-to-threshold allowed)
    dec = [x for x in res if x["K_sw"] == 0 and x["f_em"] >= 1e-3 and not x.get("class_error")]
    if quick:
        dec = []          # C1 (several hours) and C2 only in the full run
    dec = sorted(dec, key=lambda x: x["deltaA"])[:1] + sorted([x for x in dec if x["deltaA"] <= strict], key=lambda x: -x["deltaA"])[:1]
    c1, c2 = [], []
    import boltz as B
    hi = []
    for x in dec:
        r = next(m for m in meta if m[1] == x["tag"])[0]
        pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
        path, T, _ = psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], x["tag"] + "_hi", refine=2)
        hi.append(packet_params(path, T, r["f_em"], prec=NCDM_PREC_HI, quad=ROW_QUAD_HI))
    if hi:
        print(f"C1: {len(hi)} high-precision CLASS runs ({ROW_NQ_HI} nodes, l_max_ncdm 50, tolerances x0.1; the slowest step, ~6-12 h)", flush=True)
    Phis = pmap(hi, workers, [x["tag"] + " (C1)" for x in dec], "C1") if hi else []
    for x, Phi in zip(dec, Phis):
        r, _, idx, _ = next(m for m in meta if m[1] == x["tag"])
        pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
        if Phi is None:
            c1.append(dict(tag=x["tag"], deltaA=x["deltaA"], class_error=True, rel=float("inf")))
        else:
            dAh = deltaA(Phi / P0, P0); c1.append(dict(tag=x["tag"], deltaA=x["deltaA"], deltaA_hi=dAh, rel=abs(x["deltaA"] / dAh - 1)))
        sp, _ = B.from_packets(pk["q"] / pk["m_GeV"], pk["nc"], OMEGA_DM * r["f_em"] / B.h ** 2, nq=40, lmax=10)
        kb = np.geomspace(0.05, KMAX, 46); d0, _ = B.solve(kb * B.h); d1, _ = B.solve(kb * B.h, sp, phase_max=0.5)
        T2b = np.interp(np.log(KH), np.log(kb), (d1 / d0) ** 2); T2c = Prows[idx] / P0; s = (KH >= 0.5) & (KH <= 20)
        c2.append(dict(tag=x["tag"], deltaA_CLASS=x["deltaA"], deltaA_boltz=B.deltaA(kb, (d1 / d0) ** 2, B.p_cdm_shape(kb * B.h, d0)),
                       max_rel_T2_0p5_20=float(np.max(np.abs(T2b[s] / T2c[s] - 1)))))
    out["C1"] = dict(cases=c1, passed=bool(c1 and all(c["rel"] < 0.02 for c in c1))); out["C2"] = dict(cases=c2)
    print("C1", out["C1"]); print("C2", out["C2"], flush=True)

    # verdict (ESPEC_U2 §6 with P5 -> LYA_primary); network gate from S6 (or U4b if it passed)
    names = sorted({x["branch"] for x in res})
    def verdict(key):
        ok = lambda b, rl: any(x[key] for x in res if x["branch"] == b and x["R_label"] == rl and x["K_sw"] == 0)
        plus = {b: ok(b, "+1s") for b in names}; minus = {b: ok(b, "-1s") for b in names}
        v = ("CIERRE_FISICO" if all(plus.values()) else "FALSADA_EN_FORMA_ACTUAL" if not any(minus.values())
             else "CIERRE_SEGUN_RAMA" if any(plus.values()) else "INCONCLUSO")
        return dict(verdict=v, plus1s=plus, minus1s=minus, central={b: ok(b, "central") for b in names})
    out["verdict_lya"] = {k: verdict(k) for k in ("LYA_primary", "LYA_shape_strict", "LYA_shape_wide", "LYA_dA_strict", "LYA_dA_wide")}
    gate = bool((s6.get("U1_control_gate") or {}).get("passed"))
    u4b = ROOT / "network" / "results" / "u4" / "U4B_RESULT.json"
    if u4b.exists():
        blocks = json.loads(u4b.read_text()).get("blocks", [])
        if blocks and blocks[-1].get("status") == "GATES_PASSED": gate = True
    out["network_gate_passed"] = gate
    out["validated"] = bool(out["C0"]["passed"] and out["C0b"]["passed"] and out["C1"]["passed"])
    out["verdict_preregistered"] = (out["verdict_lya"]["LYA_primary"]["verdict"] if gate and out["validated"] else
                                    "INCONCLUSO (" + ("red: control U(1) fuera de rango" if not gate else "validación C0/C0b/C1 fallida") + ")")
    if quick:
        out["verdict_preregistered"] = "QUICK (comprobación de la tubería; sin veredicto)"
    out["n_class_errors"] = sum(bool(x.get("class_error")) for x in res)
    out["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U4A_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print("VERDICT (prereg):", out["verdict_preregistered"])
    if quick:   # only central rows exist: the +-1 sigma logic does not apply
        print("QUICK, filas centrales permitidas (LYA primary):", out["verdict_lya"]["LYA_primary"]["central"])
    else:
        print("LYA primary:", out["verdict_lya"]["LYA_primary"])
    return out


if __name__ == "__main__":
    main(sys.argv)
