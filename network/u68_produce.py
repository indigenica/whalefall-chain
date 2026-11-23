"""U68 (docs/preregistros/ESPEC_U68.md): the RP3 string-wall collapse with a general quadratic breaking -kappa Tr(E Q), i.e. the
kern.c tilt -EPS_G sum_a w_a Q_aa (set_tilt_weights), masses m_c^2 = m_G^2 (w_0 - w_c); U19b sparse protocol at N = 192, l_G = 2.4,
seeds 91-93, protocols O3 (1,1,1), GEN (1,2,4) and O2 (1,1,4, informative).  The counter is U19b's, with one mass per mode.
    python -B network/u68_produce.py checks        G0-G2                      -> network/results/u68/U68_CHECKS.json
    python -B network/u68_produce.py produce       the 9 networks (skip-if-present) -> network/results/u68/runs/<TAG>.json"""
import ctypes, json, math, subprocess, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
import u19_produce as U19
import u19b_produce as B

HERE = Path(__file__).resolve().parent; OUT = HERE / "results" / "u68"
ELL_G, N_PROD, SEEDS = 2.4, 192, (91, 92, 93)
PROT = {"O3": (1.0, 0.0, 0.0, 0.0), "GEN": (2.0, 1.0, 0.0, -2.0), "O2": (1.75, 0.75, 0.75, -2.25)}
MR, DTAU, F = U19.MR, U19.DTAU, U19.F
PAIRS = U19.PAIRS


def L():
    lib = U19.L()
    if not hasattr(lib, "_u68"):
        D = ctypes.c_double
        lib.set_tilt_weights.argtypes = [D, D, D, D]; lib.set_tilt_weights.restype = None; lib._u68 = True
    return lib


def mass2(w): return [w[0] - w[c] for c in (1, 2, 3)]                       # m_c^2 / m_G^2


# ------------------------------------------------------------------ potential, vacuum and effective masses
def vacuum_modes(w, x):
    """B.vacuum with the weights: vacuum q* of V_sym - x^2 eps sum_a w_a q_a; m_eff,c^2 = eigenvalue of Hess(V_sym) whose eigenvector
    lies on the off-diagonal component Q_0c, divided by x^2; returns m_eff,c / (m_G sqrt(w0 - wc)) and zchi = Lambda - mu"""
    from scipy.optimize import minimize
    MG = B.set_mass(ELL_G); eps = MG ** 2; x = max(x, 1.0)
    q0 = np.zeros(10); q0[0] = 1.0
    W = lambda q: B.vsym(1, q) - x * x * eps * sum(w[a] * q[a] for a in range(4))
    q = minimize(W, q0, method="BFGS", options=dict(gtol=1e-13)).x
    h = 1e-4; Hs = np.zeros((10, 10))
    for i in range(10):
        for j in range(i, 10):
            e_i = np.eye(10)[i] * h; e_j = np.eye(10)[j] * h
            Hs[i, j] = Hs[j, i] = (B.vsym(1, q + e_i + e_j) - B.vsym(1, q + e_i - e_j) - B.vsym(1, q - e_i + e_j)
                                   + B.vsym(1, q - e_i - e_j)) / (4 * h * h)
    lam, vec = np.linalg.eigh(Hs); m2 = mass2(w); ratio = []
    for c in (1, 2, 3):
        k = int(np.argmax(np.abs(vec[3 + c, :])))                             # eigenvector on q[4 + (c - 1)] = sqrt2 Q_0c
        ratio.append(math.sqrt(max(lam[k], 0.0)) / x / (MG * math.sqrt(m2[c - 1])))
    d = q[:4]                                                               # ESPEC_U68b enmienda 1: z_c = Lambda - mu_c per mode
    return ratio, [float(d[0] - d[c]) for c in (1, 2, 3)]


def meff_table(w):
    # ESPEC_U68b: up to x = 2.6 (the WKB check lands at x ~ 2.1-2.3), same 0.05 step, so the rows in [1, 2] are those of U68
    xs = np.round(np.arange(1.0, 2.6 + 1e-9, 0.05), 10); v = [vacuum_modes(w, x) for x in xs]
    return dict(x=xs.tolist(), ratio=[r for r, _ in v], zchi=[z for _, z in v])


def meff(table, tau, tau_A):
    xx = max(tau / tau_A, 1.0)
    r = [float(np.interp(xx, table["x"], [row[c] for row in table["ratio"]])) for c in range(3)]
    z = [float(np.interp(xx, table["x"], [row[c] for row in table["zchi"]])) for c in range(3)]
    return r, z


def site_stats(N, f, w):
    """U19.site_stats with the weighted tilt: V_eps - V_min = eps (w0 - sum_a w_a Q_aa) (>= 0 on the vacuum manifold)"""
    q0 = f[0].astype(np.float64); cth = 2 * q0 - 1
    S = sum(w[a] * f[a].astype(np.float64) for a in range(4))
    return dict(wall_frac=float(np.mean(cth < -0.5)), far_frac=float(np.mean(cth < 0.9)), V_eps=float(np.mean(U19.MG ** 2 * (w[0] - S))))


# ------------------------------------------------------------------ the counter, one mass per mode
def linear_number_modes(N, f, p, a, kbins, masses, hubble=None, zchi=1.0):
    """U19.linear_number (RP3 branch) with the mass of mode c (chi_c = 2F n_c, c = 1, 2, 3) in w_k and in the WKB term;
    returns per mode (total, dict(n, kn))"""
    V = N ** 3; ia3 = 1.0 / a ** 3
    n = np.empty((4, V), np.float32); L().principal(N, np.ascontiguousarray(f), n, None)
    sg = np.where(n[0] < 0, -1.0, 1.0).astype(np.float32); n *= sg[None]; del sg
    k = np.fft.fftfreq(N) * 2 * np.pi; kz = np.fft.rfftfreq(N) * 2 * np.pi
    kl2 = (4 * np.sin(k / 2) ** 2)[:, None, None] + (4 * np.sin(k / 2) ** 2)[None, :, None] + (4 * np.sin(kz / 2) ** 2)[None, None, :]
    wz = np.ones(kl2.shape); wz[:, :, 1:(-1 if N % 2 == 0 else None)] = 2
    kk = np.sqrt(kl2).ravel(); idx = np.digitize(kk, kbins) - 1; nb = len(kbins) - 1; ok = (idx >= 0) & (idx < nb)
    out = []
    for c in (1, 2, 3):
        qd = p[c] * n[c]
        for kq, (i, j) in enumerate(PAIRS):
            if i == c: qd = qd + p[4 + kq] * (n[j] / math.sqrt(2))
            elif j == c: qd = qd + p[4 + kq] * (n[i] / math.sqrt(2))
        chi = (2 * F) * n[c]; chid = (2 * F * ia3) * qd
        m = masses[c - 1]; w = np.sqrt(kl2 / a ** 2 + m ** 2)
        zc = zchi[c - 1] if isinstance(zchi, (list, tuple)) else zchi
        # z_c only on the field: chid = 2F (Qdot n)_c already carries Lambda - mu_c (as U19b: zchi = Lambda - mu, zchid = 1);
        # the first production scaled chid too (bug found by the Codex review, ESPEC_U68b enmienda 2)
        X = np.fft.rfftn(chi.reshape(N, N, N)) * zc; Xd = np.fft.rfftn(chid.reshape(N, N, N))
        if hubble is not None: Xd = Xd + hubble * (1 + m ** 2 / (2 * w ** 2)) * X
        nk = (wz * (np.abs(Xd) ** 2 + w ** 2 * np.abs(X) ** 2) / (2 * w)).ravel()
        hn = np.bincount(idx[ok], nk[ok], nb); hk = np.bincount(idx[ok], (kk * nk)[ok], nb)
        out.append((a ** 3 * float(nk.sum()) / V ** 2, dict(n=(a ** 3 * hn / V ** 2).tolist(), kn=(a ** 3 * hk / V ** 2).tolist())))
    return out


# ------------------------------------------------------------------ checks
def homogeneous_mode(N, c, phi):
    n = np.zeros((4, N, N, N), np.float32); n[0] = math.cos(phi); n[c] = math.sin(phi)
    return np.ascontiguousarray(Lt.q_from_n(n).reshape(10, -1), np.float32), np.zeros((10, N ** 3), np.float32)


def vacuum_mode(w, c, x, phi, N=8):
    """homogeneous field at the vacuum q*(x) of the evolved theory, rotated by phi in the (0, c) plane"""
    from scipy.optimize import minimize
    MG = B.set_mass(ELL_G); x = max(x, 1.0); q0 = np.zeros(10); q0[0] = 1.0
    q = minimize(lambda q: B.vsym(1, q) - x * x * MG ** 2 * sum(w[a] * q[a] for a in range(4)), q0, method="BFGS", options=dict(gtol=1e-14)).x
    R = np.eye(4); R[0, 0] = R[c, c] = math.cos(phi); R[0, c] = -math.sin(phi); R[c, 0] = math.sin(phi)
    Qr = R @ B._Q(q) @ R.T; qq = np.zeros(10); qq[:4] = np.diag(Qr)
    for k, (i, j) in enumerate(PAIRS): qq[4 + k] = Qr[i, j] * math.sqrt(2)
    f = np.ascontiguousarray(np.repeat(qq[:, None], N ** 3, axis=1).astype(np.float32))
    return f, np.zeros_like(f)


def g1_minkowski(w, c):
    """V1 of U19 for mode c: small homogeneous oscillation in Minkowski; omega^2 / (m_G^2 (w0 - wc)) - 1"""
    MG = B.set_mass(ELL_G); Lb = L(); Lb.set_tilt(MG ** 2); Lb.set_tilt_weights(*w)
    N = 8; dt = 0.05; f, p = homogeneous_mode(N, c, 0.01); sc = np.empty_like(f)
    m = MG * math.sqrt(mass2(w)[c - 1]); T = 2 * math.pi / m; steps = 20; ns = int(4 * T / dt / steps)
    ys, ts, t = [], [], 0.0
    for _ in range(ns):
        Lb.evolve(1, N, f, p, sc, 1.0, dt, 0.0, steps); t += steps * dt; ys.append(float(f[3 + c].mean())); ts.append(t)
    ys, ts = np.array(ys), np.array(ts); zc = ts[1:][np.diff(np.sign(ys)) != 0]; om = math.pi / np.mean(np.diff(zc))
    Lb.set_tilt(0.0); Lb.set_tilt_weights(1.0, 0.0, 0.0, 0.0)
    return float(om ** 2 / m ** 2 - 1)


def g1_wkb(w, c, table, x_c):
    """V3 of U19b for mode c: WKB mass of the homogeneous oscillation around tau = x_c tau_A against the Hessian table"""
    MG = B.set_mass(ELL_G); g = B.geometry(ELL_G, N_PROD / 2); s = g["s"]; tA = g["tau_A"]
    Lb = L(); Lb.set_tilt(MG ** 2); Lb.set_tilt_weights(*w); Lb.set_core_late(tA)
    m_nom = MG * math.sqrt(mass2(w)[c - 1]); dt = 0.01
    t = round((x_c * tA - 2 * math.pi / (s * x_c * tA * m_nom)) / dt) * dt; ys, ts = [], []
    f, p = vacuum_mode(w, c, t / tA, 0.01); sc = np.empty_like(f)                # ESPEC_U68b enmienda 1: start at the displaced vacuum
    while True:
        Lb.evolve(1, 8, f, p, sc, t, dt, s, 1); t += dt; ys.append(float(f[3 + c].mean())); ts.append(t)
        y = np.array(ys); zc = np.flatnonzero(np.diff(np.sign(y)) != 0)
        if len(zc) >= 6: break
    T = np.array(ts); y = np.array(ys)
    tz = np.array([T[i] - y[i] * (T[i + 1] - T[i]) / (y[i + 1] - y[i]) for i in zc])[1:]
    m = float(np.mean(2 * math.pi / (s * np.diff(tz ** 2)))); xm = float(np.mean(tz)) / tA
    Lb.set_tilt(0.0); Lb.set_tilt_weights(1.0, 0.0, 0.0, 0.0); Lb.set_core_late(0.0)
    r, _ = meff(table, xm * tA, tA)
    return dict(m_over_nominal=m / m_nom, hessian=r[c - 1], rel=m / m_nom / r[c - 1] - 1, x_mid=xm)


def g2_counter(w, c, N=96):
    """V4 of U19b for mode c: two standing waves of known k; number-weighted <k> of mode c against the expectation"""
    MG = B.set_mass(ELL_G); g = B.geometry(ELL_G, N_PROD / 2); a = g["s"] * g["tau_B"]; m = MG * math.sqrt(mass2(w)[c - 1])
    x = np.arange(N); X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k1, k2 = 2 * math.pi * 3 / N, 2 * math.pi * 11 / N; A1, A2 = 0.02, 0.01
    th = A1 * np.cos(k1 * X) + A2 * np.cos(k2 * Y + k2 * Z)
    n = np.zeros((4, N, N, N), np.float32); n[0] = np.cos(th / 2); n[c] = np.sin(th / 2)
    f = np.ascontiguousarray(Lt.q_from_n(n).reshape(10, -1), np.float32); p = np.zeros_like(f)
    ms = [MG * math.sqrt(v) for v in mass2(w)]
    res = linear_number_modes(N, f, p, a, B.kbins(N), ms)
    tot, h = res[c - 1]; kmean = float(np.sum(h["kn"]) / np.sum(h["n"]))
    kl = lambda k: 2 * math.sin(k / 2); kl1, kl2 = kl(k1), math.sqrt(2) * kl(k2)
    w1, w2 = math.sqrt(kl1 ** 2 / a ** 2 + m ** 2), math.sqrt(kl2 ** 2 / a ** 2 + m ** 2)
    n1, n2 = w1 * A1 ** 2, w2 * A2 ** 2; expect = (kl1 * n1 + kl2 * n2) / (n1 + n2)
    others = sum(r[0] for i, r in enumerate(res) if i != c - 1) / tot
    return dict(k_mean=kmean, expected=expect, rel=kmean / expect - 1, leak_other_modes=others, PASS=bool(abs(kmean / expect - 1) <= 0.02))


def g0_identity():
    """bit identity with origin/main's kern.c, default weights: RP3 64^3 with tilt, 40 steps, with and without the late core;
    and U(1) with tilt"""
    old = OUT / "_kern_main.c"; so_old = OUT / "_kern_main.so"
    old.write_text(subprocess.check_output(["git", "show", "origin/main:network/kern.c"], cwd=HERE.parent).decode())
    subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", str(old), "-o", str(so_old), "-lm"])
    Lold = ctypes.CDLL(str(so_old)); Lnew = L()
    fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS"); I, D = ctypes.c_int, ctypes.c_double
    Lold.evolve.argtypes = [I, I, fp, fp, fp, D, D, D, I]; Lold.set_tilt.argtypes = [D]; Lold.set_core_late.argtypes = [D]
    MG = B.set_mass(ELL_G); out = {}
    for model in (1, 0):
        for late in (False, True):
            n = 64; s = 1 / (MR * 32.0); f0, _ = Lt.initial(model, n, 3.0, 7); res = []
            for Lx in (Lold, Lnew):
                f = f0.copy(); p = np.zeros_like(f); sc = np.empty_like(f)
                Lx.set_tilt(U19.eps_of(model)); Lx.set_core_late(12.0 if late else 0.0)
                Lx.evolve(model, n, f, p, sc, 10.0, DTAU, s, 40); Lx.set_tilt(0.0); Lx.set_core_late(0.0); res.append((f, p))
            out[f"G0_{'RP3' if model else 'U1'}_late{int(late)}"] = bool(np.array_equal(res[0][0], res[1][0]) and np.array_equal(res[0][1], res[1][1]))
    for x in (old, so_old): x.unlink()
    return out


def checks():
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = dict(spec="docs/preregistros/ESPEC_U68b.md", G0=g0_identity(), G1={}, G2={}, meff={})
    for nm, w in PROT.items():
        tab = meff_table(w); res["meff"][nm] = tab
        for c in (1, 2, 3):
            mk = g1_minkowski(w, c)
            hx1 = tab["ratio"][0][c - 1]                                         # ESPEC_U68b: Hessian mass at x = 1 against Minkowski
            mk_h = (1 + mk) / hx1 ** 2 - 1
            wk = [g1_wkb(w, c, tab, xc) for xc in (2.0, B.END_FACTOR[ELL_G] * 0.97)]
            res["G1"][f"{nm}_c{c}"] = dict(minkowski_rel_nominal=mk, hessian_ratio_x1=hx1, minkowski_rel_hessian=mk_h, wkb=wk,
                                           PASS=bool(abs(mk_h) < 0.01 and all(abs(x["rel"]) < 0.003 for x in wk)))
            res["G2"][f"{nm}_c{c}"] = g2_counter(w, c)
            print(nm, c, res["G1"][f"{nm}_c{c}"], res["G2"][f"{nm}_c{c}"], flush=True)
    res["PASS"] = dict(G0=all(res["G0"].values()), G1=all(v["PASS"] for v in res["G1"].values()),
                       G2=all(v["PASS"] for v in res["G2"].values()))
    res["runtime_s"] = time.time() - t0
    (OUT / "U68_CHECKS.json").write_text(json.dumps(res, indent=1, default=float))
    print("CHECKS", res["PASS"], flush=True)


# ------------------------------------------------------------------ production
def network(nm, seed, log, table):
    w = PROT[nm]; MG = B.set_mass(ELL_G); N = N_PROD; g = B.geometry(ELL_G, N / 2); s, tA, tB, l1 = g["s"], g["tau_A"], g["tau_B"], g["ell_1"]
    ms = [MG * math.sqrt(v) for v in mass2(w)]
    ell_i, xi_i = 1.0, 0.5; drain = min(3.8, l1 - 0.3)
    tau0 = round(math.sqrt(math.exp(ell_i) * tA) / DTAU) * DTAU
    lc = (tau0 / 18.65) * math.sqrt(Lt.KXI[1] / xi_i)
    f, p = Lt.initial(1, N, lc, seed); sc = np.empty_like(f)
    Lb = L(); Lb.set_tilt(MG ** 2); Lb.set_tilt_weights(*w); Lb.set_damp(math.sqrt(math.exp(drain) / (MR * s)), MR); Lb.set_core_late(tA)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, tB, B.NSNAP) / DTAU) * DTAU)
    kb = B.kbins(N); rows, tau, t0 = [], tau0, time.time()
    rows.append(dict(tau=tau0, ell=math.log(MR * s * tau0 * tau0), a=s * tau0, initial=True, **site_stats(N, f, w)))
    for ts in taus:
        nst = int(round((ts - tau) / DTAU))
        if nst > 0: Lb.evolve(1, N, f, p, sc, tau, DTAU, s, nst); tau += nst * DTAU
        ell = math.log(MR * s * tau * tau)
        r = Lt.measure(1, N, f, p, tau, s, spectra=False, tau_f=tA); r.update(site_stats(N, f, w))
        if ell >= l1 - 0.2:
            rat, z = meff(table, tau, tA); me = [ms[c] * rat[c] for c in range(3)]
            per = linear_number_modes(N, f, p, s * tau, kb, me, hubble=1 / (s * tau ** 2), zchi=z)
            r["m_eff"] = me; r["N_lin_modes"] = [x[0] for x in per]; r["N_lin"] = float(sum(x[0] for x in per))
            if ts == taus[-1]: r["lin_spectrum_modes"] = [x[1] for x in per]
        rows.append(r)
        log(f"[U68 {nm} s{seed}] ell={ell:.2f} xi={r['xi']:.4f} wall={r['wall_frac']:.4f} far={r['far_frac']:.4f} "
            f"N_lin={r.get('N_lin', float('nan')):.4e} ({time.time() - t0:.0f}s)")
    Lb.set_tilt(0.0); Lb.set_tilt_weights(1.0, 0.0, 0.0, 0.0); Lb.set_damp(0.0, 0.0); Lb.set_core_late(0.0)
    return dict(protocol=nm, weights=w, mass2_over_mG2=mass2(w), masses=ms, N=N, tau_A=tA, tau_B=tB, seed=seed, ell_G=ELL_G, s=s, tau0=tau0,
                lc=lc, m_G=MG, ell_1=l1, a_1=g["a_1"], drain_ell=drain, k_edges=kb.tolist(), threads=Lt.lib().nthreads(),
                runtime_s=time.time() - t0, rows=rows, spec="docs/preregistros/ESPEC_U68.md")


def produce():
    R = OUT / "runs"; R.mkdir(parents=True, exist_ok=True); logf = open(OUT / "produce.log", "a")
    def lg(msg): print(msg, flush=True); logf.write(msg + "\n"); logf.flush()
    chk = OUT / "U68_CHECKS.json"
    if not chk.exists() or not all(json.loads(chk.read_text())["PASS"].values()): lg("STOP: checks missing or failed (ESPEC_U68: INCONCLUSO)"); sys.exit(3)
    tables = json.loads(chk.read_text())["meff"]
    for seed in SEEDS:                                                        # paired: the three protocols of a seed together
        for nm in PROT:
            tag = f"RP3_N{N_PROD}_{nm}_s{seed}"
            if (R / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
            res = network(nm, seed, lg, tables[nm]); (R / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    import os
    if os.environ.get("WHALEFALL_BACKEND") == "jax": raise SystemExit("U68 needs kern.c (set_tilt_weights is not in the JAX kernel)")
    {"checks": checks, "produce": produce}[sys.argv[1]]()
