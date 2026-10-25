"""U15 (docs/preregistros/ESPEC_U15.md): the warm tail of the whale against the public Lyman-ABGD likelihood (arXiv:2206.08188).

The likelihood code (montepython_public, montepython/likelihoods/Lya_abgd, commit 5b365b0 of 2024-06-20) is used unmodified,
through a minimal stand-in for the MontePython classes.  The only deviation (declared in the ESPEC): the LCDM-equivalent model is
built here (ncdm moved to CDM, same cosmology) and handed to the likelihood, without its eta rescaling for extra relativistic
species (our DeltaN_eff <= 0.0075).  Astro nuisances are profiled (bounded minimisation, 8 starts).

    MP_DIR=/path/to/montepython_public python -B lya/u15_lya_abgd.py --workers 2           # validation + rows
    MP_DIR=...  python -B lya/u15_lya_abgd.py --validate-only

Output: lya/results/u15/U15_RESULT.json"""
import itertools, json, math, os, re, sys, time, types
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "chain"))
OUT = HERE / "results" / "u15"
MP_DIR = Path(os.environ.get("MP_DIR", "/opt/src/montepython_public"))
ASTRO = {"T0a": (0.74, 0.0, 2.0), "T0s": (-4.38, -5.0, 5.0), "gamma_a": (1.45, 1.0, 1.7), "gamma_s": (-1.93, -5.0, 5.0),
         "Fz1": (0.35, 0.22, 0.54), "Fz2": (0.26, 0.15, 0.43), "Fz3": (0.18, 0.07, 0.31), "Fz4": (0.07, 0.024, 0.2), "F_UV": (0.0, 0.0, 1.0)}
LYA_EXTRA = {"output": "mPk", "z_max_pk": 3.0, "P_k_max_h/Mpc": 300.0}


# ---------------- minimal MontePython stand-in ----------------
def install_stubs():
    class ConfigurationError(Exception): pass
    class LikelihoodError(Exception): pass
    io_mp = types.ModuleType("io_mp"); io_mp.ConfigurationError = ConfigurationError; io_mp.LikelihoodError = LikelihoodError
    sys.modules["io_mp"] = io_mp

    class Likelihood:
        def __init__(self, path, data, command_line):
            self.name = self.__class__.__name__
            self.data_directory = ""
            for line in Path(path).read_text().splitlines():
                if line.strip().startswith("#") or (self.name + ".") not in line: continue
                m = re.match(r"\s*%s\.(\w+)\s*=\s*(.*)" % self.name, line)
                if m: setattr(self, m.group(1), eval(m.group(2), {"data": data, "np": np}))
        def need_cosmo_arguments(self, data, d): pass
    mp = types.ModuleType("montepython"); lc = types.ModuleType("montepython.likelihood_class"); lc.Likelihood = Likelihood
    sys.modules["montepython"] = mp; sys.modules["montepython.likelihood_class"] = lc
    init = types.ModuleType("initialise"); init.recover_cosmological_module = lambda data: data._lcdm_equiv
    sys.modules["initialise"] = init
    sys.path.insert(0, str(MP_DIR / "montepython" / "likelihoods" / "Lya_abgd"))


class FakeEquiv:
    """Pre-computed LCDM-equivalent CLASS object; ignores the parameters the likelihood passes (see module docstring)."""
    def __init__(self, cl): self.cl, self.state = cl, False
    def struct_cleanup(self): pass
    def empty(self): pass
    def set(self, p): pass
    def compute(self): pass
    def h(self): return self.cl.h()
    def pk_lin(self, k, z): return self.cl.pk_lin(k, z)


class Data:
    def __init__(self, folder):
        self.path = {"data": str(MP_DIR / "data"), "MontePython": str(MP_DIR / "montepython")}
        self.mcmc_parameters = {n: {"current": v[0], "scale": 1.0} for n, v in ASTRO.items()}
        self.mcmc_parameters.update({"sigma8": {"current": 0.8, "scale": 1.0}, "z_reio": {"current": 8.0, "scale": 1.0}})
        self.derived_lkl = {}; self.cosmo_arguments = {"N_ncdm": 1, "omega_ncdm": 0.0}   # only read for the 'something removed' check
        self.boundary_loglike = -1e30; self.log_flag = True; self.folder = folder
    def get_mcmc_parameters(self, roles):
        r = roles[0]
        if r == "varying": return list(ASTRO)
        if r == "derived": return ["sigma8", "z_reio"]
        if r == "derived_lkl": return ["lya_neff", "area_criterion", "self_distance_weighted", "weight_largest", "weightsum_2", "weightsum_4"]
        return []


_LKL = {}
def likelihood():
    if "L" not in _LKL:
        install_stubs()
        import importlib.util
        spec = importlib.util.spec_from_file_location("Lya_abgd", MP_DIR / "montepython" / "likelihoods" / "Lya_abgd" / "__init__.py")
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        folder = str(OUT / "mp_run"); os.makedirs(folder, exist_ok=True)
        data = Data(folder); cl = types.SimpleNamespace(folder=folder)
        L = mod.Lya_abgd(str(MP_DIR / "montepython" / "likelihoods" / "Lya_abgd" / "Lya_abgd.data"), data, cl)
        L.verbose = 0; _LKL.update(L=L, data=data)
    return _LKL["L"], _LKL["data"]


def compute(params):
    from classy import Class
    c = Class(); c.set(params); c.compute(); return c


def equiv_params(p):
    q = {k: v for k, v in p.items() if not (k.startswith("ncdm") or k in ("N_ncdm", "use_ncdm_psd_files", "m_ncdm", "T_ncdm", "omega_ncdm",
                                                                             "l_max_ncdm", "tol_ncdm_synchronous", "tol_ncdm_bg",
                                                                             "ncdm_fluid_approximation", "Quadrature strategy"))}
    q["omega_cdm"] = p.get("omega_cdm", 0.0) + p.get("omega_ncdm", 0.0)
    return q


def profile(model_params, lcdm=False, starts=8, seed=15):
    """chi2_min over the astro nuisances for one model (CLASS run once)."""
    from scipy.optimize import minimize
    L, data = likelihood()
    p = dict(model_params, **LYA_EXTRA)
    cm = compute(p); ce = cm if lcdm else compute(dict(equiv_params(p), **LYA_EXTRA))
    data._lcdm_equiv = FakeEquiv(ce); L.isLCDM = lcdm
    if lcdm: data.cosmo_arguments = {}
    names = list(ASTRO); lo = np.array([ASTRO[n][1] for n in names]); hi = np.array([ASTRO[n][2] for n in names])
    zmh = np.array([4.2, 4.6, 5.0, 5.4]); iz = slice(6, 10)                 # MIKE/HIRES z within redshift_list (ESPEC_U15, amendment 1)
    def in_thermal_grid(a):
        T0 = a["T0a"] * ((1 + zmh) / (1 + L.zp)) ** a["T0s"] * 1e4; g = a["gamma_a"] * ((1 + zmh) / (1 + L.zp)) ** a["gamma_s"]
        return bool(np.all(T0 >= L.t0_min[iz]) and np.all(T0 <= L.t0_max[iz]) and np.all(g >= L.slope_min[iz]) and np.all(g <= L.slope_max[iz]))
    def chi2(x):
        if not in_thermal_grid(dict(zip(names, x))): return 1e30
        for n, v in zip(names, x): data.mcmc_parameters[n]["current"] = float(v)
        ll = L.loglkl(cm, data)
        return -2 * ll if ll > -1e29 else 1e30
    rng = np.random.default_rng(seed); best = None
    x0s = [np.array([ASTRO[n][0] for n in names])]
    while len(x0s) < starts:
        x = lo + (hi - lo) * rng.random(len(names))
        if in_thermal_grid(dict(zip(names, x))): x0s.append(x)
    first = chi2(x0s[0])
    if first >= 1e29:
        err = open(os.path.join(data.folder, L.bin_file_name)).read().strip().splitlines()[-1].split()[-1] if os.path.exists(os.path.join(data.folder, L.bin_file_name)) else "error"
        return dict(chi2=None, error=str(err), derived=dict(data.derived_lkl))
    for x0 in x0s:
        r = minimize(chi2, x0, method="Powell", bounds=list(zip(lo, hi)), options=dict(xtol=1e-3, ftol=1e-4, maxfev=4000))
        if best is None or r.fun < best.fun: best = r
    chi2(best.x)
    return dict(chi2=float(best.fun), astro=dict(zip(names, map(float, best.x))), derived={k: float(v) for k, v in data.derived_lkl.items()},
                sigma8=float(cm.sigma8()))


def validation(U4):
    out = {}
    base = profile(dict(U4.COSMO, omega_cdm=U4.OMEGA_DM), lcdm=True); out["LCDM"] = base
    cases = {"WDM 2 keV": (2.0, 1.0), "WDM 3.5 keV": (3.5, 1.0), "WDM 5.3 keV": (5.3, 1.0), "WDM 10 keV": (10.0, 1.0),
             "mix 1 keV f=0.5": (1.0, 0.5), "mix 2 keV f=0.2": (2.0, 0.2)}
    for lab, (m, f) in cases.items():
        r = profile(U4.thermal_params(m, f)); r["dchi2"] = (r["chi2"] - base["chi2"]) if r["chi2"] is not None else None; out[lab] = r
        print(f"V {lab:<18} dchi2 = {r['dchi2']}", flush=True)
    d = lambda k: out[k]["dchi2"]
    ok = all(d(k) is not None for k in cases) and d("WDM 2 keV") > 4 and d("WDM 10 keV") < 4 and d("mix 1 keV f=0.5") > 4 and d("mix 2 keV f=0.2") < 4 \
        and d("WDM 2 keV") > d("WDM 3.5 keV") > d("WDM 5.3 keV") > d("WDM 10 keV")
    out["passed"] = bool(ok); return out


def row_specs():
    import u5_class as U5, u12_class as U12
    rows = [dict(set="B2_R1.56_q1", branch="B2", R=1.56, x0=x0, K_sw=K, spectrum="powerlaw_q1", eps=U5.eps_B2(1.56)) for x0, K in itertools.product((10, 30), (0, 1))]
    u7a = json.loads((HERE / "results" / "u7" / "U7C_RESULT.json").read_text())["U7a"]
    for fam, R, x0, K in itertools.product(("B2k", "B2z"), (1.39, 1.69, 2.09), (10, 30), (0, 1)):
        kap = u7a[str(x0)]["kappa"]
        eps = (lambda ell, R=R, k=kap: k * R * 1.20 / (0.24 * np.asarray(ell))) if fam == "B2k" else (lambda ell, R=R: R * 1.49 / (0.24 * np.asarray(ell)))
        rows.append(dict(set="U7c", family=fam, branch="B2", R=R, x0=x0, K_sw=K, spectrum="energy_logflat", eps=eps))
    return rows


def run_row(sp):
    import u4a_class as U4, s6_measured as S6
    from s4_dark import lambda_wdm_reference
    from ucommon import benchmark
    Bm = benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    dneff_late = json.loads((__import__("ucommon").RESULTS / "s3_sectors.json").read_text())["DeltaNeff_late_homogeneous_max"]
    lam_wdm, _ = lambda_wdm_reference()
    r = S6.node(F, Trh, dneff_late, lam_wdm, sp["eps"], sp["spectrum"], sp["x0"], sp["K_sw"])
    tag = f"U15_{sp['set']}_{sp.get('family', 'B2')}_R{sp['R']:g}_x{sp['x0']}_K{sp['K_sw']}"
    out = {k: v for k, v in sp.items() if k != "eps"}; out.update(tag=tag, f_em=r.get("f_em"), m_G_eV=r.get("m_G_eV"))
    if not (r.get("P1") and all(r.get(k) for k in ("P2", "P3", "P4"))): out.update(status="S6_EXCLUIDA"); return out
    pk = S6.source_eps(F, r["m_G_eV"], Trh, sp["x0"], sp["eps"], spectrum=sp["spectrum"])
    path, T, _ = U4.psd_from_packets(pk["q"] / pk["m_GeV"], pk["nc"], tag)
    out["lkl"] = profile(U4.packet_params(path, T, r["f_em"]))
    return out


def main(argv):
    import u4a_class as U4
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = dict(spec="docs/preregistros/ESPEC_U15.md", montepython_commit="5b365b093982099b926a0c0c169697c9dc6cbbc4")
    res["validation"] = validation(U4); print("validation passed:", res["validation"]["passed"], flush=True)
    (OUT / "U15_RESULT.json").write_text(json.dumps(res, indent=1, default=float))
    if "--validate-only" in argv or not res["validation"]["passed"]:
        if not res["validation"]["passed"]: res["verdict_U15"] = "INCONCLUSO (arnés no validado)"
        (OUT / "U15_RESULT.json").write_text(json.dumps(res, indent=1, default=float)); return res
    base = res["validation"]["LCDM"]["chi2"]; workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 1
    specs = row_specs()
    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(workers) as ex: rows = list(ex.map(run_row, specs))
    else:
        rows = [run_row(s) for s in specs]
    for r in rows:
        l = r.get("lkl")
        if l is None: continue
        if l["chi2"] is None: r["status"] = "FUERA_DE_VALIDEZ"; continue
        r["dchi2"] = l["chi2"] - base; r["status"] = "PERMITIDA" if r["dchi2"] <= 4 else "EXCLUIDA"
        r["self_distance_flag"] = bool(l["derived"].get("self_distance_weighted", 0) > 1)
        print(f"{r['tag']:<40} f_em={r['f_em']:.3f} dchi2={r['dchi2']:7.2f} {r['status']}", flush=True)
    res["rows"] = rows
    b2 = [r for r in rows if r["set"] == "B2_R1.56_q1"]; valid = [r for r in b2 if r.get("status") in ("PERMITIDA", "EXCLUIDA")]
    res["verdict_U15"] = ("INCONCLUSO" if len(valid) < 3 else "COLA_TIBIA_PERMITIDA" if all(r["status"] == "PERMITIDA" for r in valid)
                          else "COLA_TIBIA_EXCLUIDA" if all(r["status"] == "EXCLUIDA" for r in valid) else "COLA_TIBIA_PARCIAL")
    u7 = [r for r in rows if r["set"] == "U7c"]
    res["U7c_table"] = {f"{f}|R{R:g}": f"{sum(r.get('status') == 'PERMITIDA' for r in u7 if r['family'] == f and r['R'] == R)}/4"
                        for f in ("B2k", "B2z") for R in (1.39, 1.69, 2.09)}
    res["runtime_h"] = (time.time() - t0) / 3600
    (OUT / "U15_RESULT.json").write_text(json.dumps(res, indent=1, default=float))
    print("U15:", res["verdict_U15"], res["U7c_table"], flush=True); return res


if __name__ == "__main__":
    main(sys.argv)
