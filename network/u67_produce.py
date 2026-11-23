"""U67 (docs/preregistros/ESPEC_U67.md): the RP3 network born from the xi_Q R ordering of Q in a matter era, then radiation,
against the radiation-only production protocol, at 192^3 with 3 seeds.

    python -B network/u67_produce.py --checks      # G0 (bit identity with main's kern.c) and G1 (background, X0, Hessian)
    nohup python -B network/u67_produce.py > network/results/u67/produce.out 2>&1 &     # the 12 networks (skips done ones)
    python -B network/u67_analyze.py               # gates G2-G3 and the verdict"""
import ctypes, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import lattice as Lt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "u67"
N, TAU_END, DTAU, NSNAP = 192, 96.0, 0.25, 30
SEEDS = (101, 102, 103)
XI_Q, ELL0, ELL_D, AMP, LC_SEED = 0.25, 0.5, 3.8, 0.02, 1.0
M = math.sqrt(Lt.M2[1])
PROTOCOLS = {"RD_SPARSE": None, "RD_DENSE": None, "MD2": 2.0, "MD3": 3.0}


# ---------- the symmetric point of U_code and its Hessian ----------
def Vq(q):
    """U_code / Z_Q on the 10-comp Frobenius vector (as kern.c: 0.2 (TrQ - 1)^2 + 0.12 Tr[(Q^2 - Q)^2])"""
    Q = np.zeros((4, 4)); pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for i in range(4): Q[i, i] = q[i]
    for k, (i, j) in enumerate(pairs): Q[i, j] = Q[j, i] = q[4 + k] / math.sqrt(2)
    Mx = Q @ Q - Q
    return 0.25 * 0.8 * (np.trace(Q) - 1) ** 2 + 0.25 * 0.48 * np.trace(Mx @ Mx)


def x0_symmetric():
    dU = lambda x: 1.6 * (4 * x - 1) + 0.96 * (x * x - x) * (2 * x - 1)        # d/dx of U_code(x 1)
    x0 = brentq(dU, 0.1, 0.4, xtol=1e-16, rtol=1e-15)
    return x0, dU(x0)


def hessian(x0, h=1e-4):
    q0 = np.array([x0] * 4 + [0.0] * 6); H = np.zeros((10, 10))
    for i in range(10):
        for j in range(10):
            e = lambda di, dj: Vq(q0 + di * np.eye(10)[i] + dj * np.eye(10)[j])
            H[i, j] = (e(h, h) - e(h, -h) - e(-h, h) + e(-h, -h)) / (4 * h * h)
    return H


# ---------- background ----------
def md_params(ell_tr):
    """a_tr, tau_tr, ell_end with tau_end = N/2 and a(tau_end) = 1/m (ESPEC_U67 §2)"""
    c = math.exp(ell_tr / 2)
    x = (-c + math.sqrt(c * c + 4 * TAU_END)) / 2          # x = exp(ell_end / 2): x^2 + c x - tau_end = 0
    ell_end = 2 * math.log(x)
    return math.exp((ell_tr - ell_end) / 2) / M, 2 * math.exp((ell_tr + ell_end) / 2), ell_end


def a_of(t, atr, ttr): x = t / ttr; return atr * x * x if t <= ttr else atr * (2 * x - 1)
def H_of(t, atr, ttr):
    a = a_of(t, atr, ttr)
    return 2 / (a * t) if t <= ttr else 2 * atr / (ttr * a * a)
def ell_of(t, atr, ttr): return math.log(M / H_of(t, atr, ttr))
def tau_of(ell, atr, ttr): return brentq(lambda t: ell_of(t, atr, ttr) - ell, 1e-3, 10 * TAU_END, xtol=1e-12)


def lib():
    L = Lt.lib()
    L.set_bg_mdrd.argtypes = [ctypes.c_double, ctypes.c_double]; L.set_bg_mdrd.restype = None
    L.set_xiq.argtypes = [ctypes.c_double, ctypes.c_double]; L.set_xiq.restype = None
    L.a_mdrd.argtypes = [ctypes.c_double]; L.a_mdrd.restype = ctypes.c_double
    return L


def measure_md(f, p, tau, atr, ttr):
    L = lib(); a = a_of(tau, atr, ttr); H = H_of(tau, atr, ttr)
    n = np.empty((4, N ** 3), np.float32); gap = np.empty(N ** 3, np.float32)
    L.principal(N, f, n, gap.ctypes.data)
    seed = np.zeros(N ** 3, np.uint8); cnt = L.plaquettes(1, N, n, seed)
    E = np.zeros(3); L.energy(1, N, f, p, a, E)
    KG = np.empty(N ** 3, np.float32); KM = np.empty(N ** 3, np.float32); L.kinetic_split(1, N, f, p, n, a, KG, KM)
    return dict(tau=tau, a=a, H=H, ell=math.log(M / H), n_plaq=int(cnt), xi=(2 / 3) * cnt / (4 * (a * H) ** 2 * N ** 3),
                gap=float(gap.mean(dtype=np.float64)), E_kin=E[0], E_grad=E[1], E_pot=E[2],
                KG_all=float(KG.mean(dtype=np.float64)), KM_all=float(KM.mean(dtype=np.float64)), md=bool(tau < ttr))


def run_md(ell_tr, seed, x0, log=print, start="symmetric", ell0=None, amp=None, lc_seed=None):
    """start = "symmetric" (U67, and MDO2 of U67b with amp/lc_seed) or "vacuum" (U67b MDA: synthetic RP3 vacuum at ell0 = 1.3,
    correlation length from the producer's rule l_c = (tau0/18.65) sqrt(K_xi/0.5))"""
    if os.environ.get("WHALEFALL_BACKEND") == "jax":      # jax_kernel has neither switch: it would run a = s tau, no xi_Q R
        raise RuntimeError("U67 matter-era runs need the C kernel (unset WHALEFALL_BACKEND)")
    L = lib(); atr, ttr, ell_end = md_params(ell_tr)
    ell0 = ELL0 if ell0 is None else ell0; amp = AMP if amp is None else amp; lc_seed = LC_SEED if lc_seed is None else lc_seed
    tau0 = round(tau_of(ell0, atr, ttr) / DTAU) * DTAU
    if start == "vacuum":
        lc_seed = (tau0 / 18.65) * math.sqrt(Lt.KXI[1] / 0.5)
        f, p = Lt.initial(1, N, lc_seed, seed); amp = None
    else:
        rng = np.random.default_rng(seed)
        g = Lt.gaussian_field(N, 10, lc_seed, rng)
        g /= g.reshape(10, -1).std(axis=1)[:, None, None, None]
        f = (amp * g).reshape(10, -1).astype(np.float32); f[:4] += np.float32(x0); f = np.ascontiguousarray(f)
        p = np.zeros_like(f)
    scratch = np.empty_like(f)
    L.set_bg_mdrd(atr, ttr); L.set_xiq(XI_Q, x0); L.set_damp(tau_of(ELL_D, atr, ttr), M)
    taus = np.unique(np.round(np.concatenate([np.geomspace(tau0 * 1.05, TAU_END, NSNAP - 6),
                                              np.linspace(tau_of(3.9, atr, ttr), TAU_END, 6)]) / DTAU) * DTAU)
    rows, tau, t0 = [measure_md(f, p, tau0, atr, ttr)], tau0, time.time()
    try:
        for ts in taus:
            nst = int(round((ts - tau) / DTAU))
            if nst <= 0: continue
            L.evolve(1, N, f, p, scratch, tau, DTAU, 1.0, nst); tau += nst * DTAU
            r = measure_md(f, p, tau, atr, ttr); rows.append(r)
            log(f"[MD ell_tr={ell_tr} seed={seed}] tau={tau:.2f} ell={r['ell']:.2f} xi={r['xi']:.3f} gap={r['gap']:.3f} "
                f"plaq={r['n_plaq']} ({time.time() - t0:.0f}s)")
    finally:
        L.set_bg_mdrd(0.0, 0.0); L.set_xiq(0.0, 0.0); L.set_damp(0.0, 0.0)
    return dict(model="RP3", protocol=f"MD{ell_tr:g}", start=start, N=N, seed=seed, tau0=tau0, tau_end=TAU_END, dtau=DTAU, ell0=ell0,
                ell_tr=ell_tr, ell_end=ell_end, a_tr=atr, tau_tr=ttr, xi_q=XI_Q, x0=x0, amp=amp, lc=lc_seed, ell_d=ELL_D,
                threads=L.nthreads(), runtime_s=time.time() - t0, rows=rows)


def checks():
    """G0: bit identity with main's kern.c (both switches off, with and without drain); G1: background and X0"""
    out = {}
    old = OUT / "_kern_main.c"; so_old = OUT / "_kern_main.so"
    old.write_text(subprocess.check_output(["git", "show", "origin/main:network/kern.c"], cwd=HERE.parent).decode())
    subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", str(old), "-o", str(so_old), "-lm"])
    Lold = ctypes.CDLL(str(so_old)); Lnew = lib()
    fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS"); I, D = ctypes.c_int, ctypes.c_double
    for Lx in (Lold,):
        Lx.evolve.argtypes = [I, I, fp, fp, fp, D, D, D, I]; Lx.set_damp.argtypes = [D, D]
    n = 64; s = 1 / (M * 32.0)
    for drain in (False, True):
        f0, _ = Lt.initial(1, n, 3.0, 7); res = []
        for Lx in (Lold, Lnew):
            f = f0.copy(); p = np.zeros_like(f); sc = np.empty_like(f)
            Lx.set_damp(20.0 if drain else 0.0, M if drain else 0.0)
            Lx.evolve(1, n, f, p, sc, 10.0, DTAU, s, 40); Lx.set_damp(0.0, 0.0); res.append((f, p))
        out[f"G0_identical_drain{int(drain)}"] = bool(np.array_equal(res[0][0], res[1][0]) and np.array_equal(res[0][1], res[1][1]))
    x0, dU = x0_symmetric(); H = hessian(x0); w, vec = np.linalg.eigh(H)
    tr = np.array([0.5] * 4 + [0.0] * 6)
    trace_overlap = np.abs(vec.T @ tr)
    i_tr = int(np.argmax(trace_overlap))
    m2_tr = float(w[i_tr]); m2_tl = [float(x) for k, x in enumerate(w) if k != i_tr]
    cont = []
    for ell_tr in (2.0, 3.0):
        atr, ttr, _ = md_params(ell_tr)
        Lnew.set_bg_mdrd(atr, ttr)
        # the two branches evaluated at tau_tr itself (first version of this check compared a(tau_tr -+ e), which differ by
        # 2 e a' for any smooth a: a bug of the check, fixed before production and declared in the report)
        a_md, a_rd = atr * (ttr / ttr) ** 2, atr * (2 * ttr / ttr - 1)
        da_md, da_rd = 2 * atr * ttr / ttr ** 2, 2 * atr / ttr
        cont.append(dict(ell_tr=ell_tr, a_jump=abs(a_md - a_rd) / atr, da_jump=abs(da_md - da_rd) / abs(da_md),
                         a_c_at_tr=abs(Lnew.a_mdrd(ttr) - atr) / atr,
                         a_end_m=a_of(TAU_END, atr, ttr) * M, ell_end=md_params(ell_tr)[2],
                         a_python_vs_c=abs(Lnew.a_mdrd(0.7 * ttr) - a_of(0.7 * ttr, atr, ttr)) / atr))
        Lnew.set_bg_mdrd(0.0, 0.0)
    ell_crit = math.log(M / math.sqrt(-max(m2_tl) / (3 * XI_Q)))
    out.update(x0=x0, dU_at_x0=dU, m2_trace=m2_tr, m2_traceless=m2_tl, ell_crit=ell_crit, background=cont,
               G1=bool(abs(dU) < 1e-12 and m2_tr > 0 and all(x < 0 for x in m2_tl)
                       and all(c["a_jump"] < 1e-12 and c["da_jump"] < 1e-12 and c["a_c_at_tr"] < 1e-12 and abs(c["a_end_m"] - 1) < 1e-12 for c in cont)))
    out["G0"] = bool(out["G0_identical_drain0"] and out["G0_identical_drain1"])
    for x in (old, so_old): x.unlink()
    (OUT / "U67_CHECKS.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


def main_u67b():
    """ESPEC_U67b: MDA2, MDA3 (synthetic vacuum start in the matter era) and MDO2 (symmetric phase, seed 0.3 / 6 cells)"""
    out = HERE / "results" / "u67b"; (out / "runs").mkdir(parents=True, exist_ok=True)
    x0 = x0_symmetric()[0]
    log = open(out / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    lg(f"START {time.strftime('%Y-%m-%d %H:%M:%S')} threads={lib().nthreads()}")
    plan = {"MDA2": dict(ell_tr=2.0, start="vacuum", ell0=1.3), "MDA3": dict(ell_tr=3.0, start="vacuum", ell0=1.3),
            "MDO2": dict(ell_tr=2.0, start="symmetric", ell0=ELL0, amp=0.3, lc_seed=6.0)}
    for sd in SEEDS:
        for prot, kw in plan.items():
            path = out / "runs" / f"RP3_N{N}_{prot}_s{sd}.json"
            if path.exists(): lg(f"skip {path.name}"); continue
            kw = dict(kw); ell_tr = kw.pop("ell_tr")
            r = run_md(ell_tr, sd, x0, log=lambda x: None, **kw); r.update(protocol=prot, spec="docs/preregistros/ESPEC_U67b.md")
            tmp = path.with_suffix(".tmp"); tmp.write_text(json.dumps(r, default=float)); tmp.rename(path)
            lg(f"DONE {path.name} {r['runtime_s']:.0f}s xi_end={r['rows'][-1]['xi']:.3f} gap_end={r['rows'][-1]['gap']:.3f} plaq_end={r['rows'][-1]['n_plaq']}")
    lg(f"END {time.strftime('%Y-%m-%d %H:%M:%S')}")


def main(argv):
    OUT.mkdir(parents=True, exist_ok=True); (OUT / "runs").mkdir(exist_ok=True)
    if "--checks" in argv: checks(); return
    if "--u67b" in argv: main_u67b(); return
    x0 = x0_symmetric()[0]
    log = open(OUT / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    lg(f"START {time.strftime('%Y-%m-%d %H:%M:%S')} threads={lib().nthreads()}")
    for sd in SEEDS:
        for prot, ell_tr in PROTOCOLS.items():
            path = OUT / "runs" / f"RP3_N{N}_{prot}_s{sd}.json"
            if path.exists(): lg(f"skip {path.name}"); continue
            if ell_tr is None:
                r = Lt.run(1, N, TAU_END, sd, sfac=1.0, xi_init=0.5 if prot == "RD_SPARSE" else 2.0, dtau=DTAU, nsnap=NSNAP,
                           spectra_every=1000, log=lambda x: None, damp_ell=ELL_D)
                r.update(protocol=prot)
            else:
                r = run_md(ell_tr, sd, x0, log=lambda x: None)
            r["spec"] = "docs/preregistros/ESPEC_U67.md"
            tmp = path.with_suffix(".tmp"); tmp.write_text(json.dumps(r, default=float)); tmp.rename(path)
            lg(f"DONE {path.name} {r['runtime_s']:.0f}s xi_end={r['rows'][-1]['xi']:.3f} plaq_end={r['rows'][-1]['n_plaq']}")
    lg(f"END {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main(sys.argv)
