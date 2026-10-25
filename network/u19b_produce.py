"""U19b (docs/preregistros/ESPEC_U19b.md): eta and epoch of the RP3 string-wall collapse, with a late frozen-core phase
(kern.c set_core_late: for tau > tau_A the symmetric potential keeps a_p = s tau_A, the tilt stays physical) that takes the
box to m_G/H ~ 94, and a k-binned linear particle counter.
    python -B network/u19b_produce.py validate            V1-V5 (cloud)            -> network/results/u19b/U19B_VALIDATION.json
    python -B network/u19b_produce.py pilot               N = 192 RP3 (informative) -> network/results/u19b/pilot/
    python -B network/u19b_produce.py produce [N]         the 12 networks, priority order, skip-if-present (WSL; N = 512)
                                                                                    -> network/results/u19b/runs/<TAG>.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
import u19_produce as U19

HERE = Path(__file__).resolve().parent; OUT = HERE / "results" / "u19b"
MR, V0, F, DTAU = U19.MR, U19.V0, U19.F, U19.DTAU
ELL_GS = (3.2, 2.4); N_PROD = 512; NSNAP = 48; SEEDS = (91, 92)
END_FACTOR = {3.2: 3.0, 2.4: 2.0}                         # tau_B / tau_A: both reach m_G/H ~ 94 at N = 512
XI1_DESIGN = 0.65                                         # RP3 sparse xi near ell_1 (U19), only for the V5 correlation length


def kbins(N):
    """comoving k_lat bin edges: [0, 2pi/N) then log bins up to the lattice maximum 2 sqrt(3)"""
    return np.concatenate([[0.0], np.geomspace(2 * math.pi / N * 0.999, 2 * math.sqrt(3) * 1.0001, 72)])


def set_mass(ell_g):
    U19.MG = MR * math.exp(-ell_g); return U19.MG


def lib():
    return U19.L()                                        # set_core_late is declared in lattice.lib()


def geometry(ell_g, tau_A=N_PROD / 2):
    """s, tau_A, tau_B, tau_1 (H = m_G/3), a_1 for the production timeline (m_r a dx = 1 at tau_A)"""
    s = 1.0 / (MR * tau_A); ell_1 = ell_g + math.log(3.0)
    tau_1 = math.sqrt(math.exp(ell_1) / (MR * s)); tau_B = END_FACTOR[ell_g] * tau_A
    return dict(s=s, tau_A=tau_A, tau_B=tau_B, tau_1=tau_1, a_1=s * tau_1, ell_1=ell_1)


# ------------------------------------------------------------------ validation
def v1(ell_g):
    set_mass(ell_g); return U19.v1_mass()


def control(model, ell_g, release, ell_rel, table, N=128, theta0=0.3, noise=0.01, seed=5):
    """V2: homogeneous control on the production timeline (late frozen core from tau_A), released early (H = 30 m_G) or frozen
    until ell_rel (dense protocol); R = m_G N_lin(tau_B) / (<V_eps - V_min>_i a_1^3) vs the homogeneous ODE"""
    MG = set_mass(ell_g); g = geometry(ell_g); s = g["s"]
    U19.L().set_tilt(U19.eps_of(model)); lib().set_core_late(g["tau_A"])
    tau_i = math.sqrt(1 / (s * 30 * MG)) if release == "early" else math.sqrt(math.exp(ell_rel) / (MR * s))
    tau_i = round(tau_i / DTAU) * DTAU; tau_f = round(g["tau_B"] / DTAU) * DTAU
    rng = np.random.default_rng(seed); f, p = U19.homogeneous_field(model, N, theta0, noise, rng); sc = np.empty_like(f)
    V_i = U19.site_stats(model, N, f)["V_eps"]
    tau_mid = round(g["tau_A"] / DTAU) * DTAU
    U19.L().evolve(model, N, f, p, sc, tau_i, DTAU, s, int(round((tau_mid - tau_i) / DTAU)))
    N_A = U19.linear_number(model, N, f, p, s * tau_mid, **counter_kw(table, tau_mid, g["tau_A"], MG, s))
    U19.L().evolve(model, N, f, p, sc, tau_mid, DTAU, s, int(round((tau_f - tau_mid) / DTAU)))
    Nl = U19.linear_number(model, N, f, p, s * tau_f, **counter_kw(table, tau_f, g["tau_A"], MG, s))
    No, Vo = U19.ode_number(theta0, tau_i, tau_f, s)
    U19.L().set_tilt(0.0); lib().set_core_late(0.0)
    a1 = g["a_1"]
    return dict(R=MG * Nl / (V_i * a1 ** 3), R_ode=MG * No / (Vo * a1 ** 3), N_lin_A=N_A, N_lin_B=Nl, N_ode=No,
                late_number_drift=Nl / N_A - 1, tau_i=tau_i, tau_A=tau_mid, tau_f=tau_f)


ISQ2 = 1 / math.sqrt(2.0); OFF = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def vsym(model, q):
    """symmetric potential of kern.c (rp3_grad / u1_grad without the tilt), double precision"""
    if model == 0: return 0.06 * (q[0] ** 2 + q[1] ** 2 - 0.5) ** 2                    # LAM_U/4 (phi^2 - V2)^2
    Q = np.diag(q[:4]).astype(float)
    for k, (i, j) in enumerate(OFF): Q[i, j] = Q[j, i] = q[4 + k] * ISQ2
    M = Q @ Q - Q
    return 0.25 * 0.8 * (np.trace(Q) - 1) ** 2 + 0.25 * 0.48 * np.trace(M @ M)       # LAM_Z, LAMP_Z


def _Q(q):
    Q = np.diag(q[:4]).astype(float)
    for k, (i, j) in enumerate(OFF): Q[i, j] = Q[j, i] = q[4 + k] * ISQ2
    return Q


def vacuum(model, ell_g, x):
    """enmiendas 1-2: the evolved theory when a / a_p = x (late frozen core).  The conformal mass matrix is a_p^2 Hess(V_sym) at
    the vacuum q* of V_sym - x^2 eps q_0 (the tilt is linear); m_eff^2 = smallest eigenvalue / x^2.  Canonical Goldstone
    normalisation at q*: U(1) chi = r_0 theta (zchi = zchid = r_0/v); RP3 Q* = Lambda n n^T + mu (1 - n n^T),
    chi_a = 2f (Lambda - mu) n_a (zchi = Lambda - mu; the counter's chidot = 2f (Qdot n)_a already carries it)."""
    from scipy.optimize import minimize
    MG = set_mass(ell_g); eps = U19.eps_of(model); x = max(x, 1.0)
    q0 = np.zeros(10 if model else 2); q0[0] = 1.0 if model else V0
    W = lambda q: vsym(model, q) - x * x * eps * q[0]
    q = minimize(W, q0, method="BFGS", options=dict(gtol=1e-13)).x
    n = len(q); h = 1e-4; Hs = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            e_i = np.eye(n)[i] * h; e_j = np.eye(n)[j] * h
            Hs[i, j] = Hs[j, i] = (vsym(model, q + e_i + e_j) - vsym(model, q + e_i - e_j) - vsym(model, q - e_i + e_j)
                                   + vsym(model, q - e_i - e_j)) / (4 * h * h)
    ratio = math.sqrt(max(np.linalg.eigvalsh(Hs)[0], 0.0)) / x / MG
    if model == 0:
        z = math.hypot(q[0], q[1]) / V0; return ratio, z, z
    lam = np.linalg.eigvalsh(_Q(q))
    return ratio, float(lam[-1] - lam[:-1].mean()), 1.0


def meff_ratio(model, ell_g, x):
    return vacuum(model, ell_g, x)[0]


def meff_table(model, ell_g):
    xs = np.linspace(1.0, END_FACTOR[ell_g], 21); v = [vacuum(model, ell_g, x) for x in xs]
    return dict(x=xs.tolist(), ratio=[r for r, _, _ in v], zchi=[z for _, z, _ in v], zchid=[z for _, _, z in v])


def meff(table, tau, tau_A, key="ratio"):
    return float(np.interp(max(tau / tau_A, 1.0), table["x"], table[key]))


def counter_kw(table, tau, tau_A, MG, s):
    """keyword arguments of U19.linear_number for the U19b counter at conformal time tau"""
    return dict(m_eff=MG * meff(table, tau, tau_A), hubble=1 / (s * tau ** 2), zchi=meff(table, tau, tau_A, "zchi"),
                zchid=meff(table, tau, tau_A, "zchid"))


def wkb_mass(model, ell_g, x_c, N=8, halfper=4):
    """V3: Goldstone mass measured from the homogeneous oscillation around tau = x_c tau_A (late frozen core active), with the
    unbiased WKB estimator: between successive zero crossings (s m_eff / 2)(tau_k+1^2 - tau_k^2) = pi"""
    MG = set_mass(ell_g); g = geometry(ell_g); s = g["s"]
    U19.L().set_tilt(U19.eps_of(model)); lib().set_core_late(g["tau_A"])
    rng = np.random.default_rng(0); f, p = U19.homogeneous_field(model, N, 0.02, 0.0, rng); sc = np.empty_like(f)
    dt = 0.01; t = round((x_c * g["tau_A"] - 2 * math.pi / (s * x_c * g["tau_A"] * MG)) / dt) * dt
    ys, ts = [], []
    while True:
        U19.L().evolve(model, N, f, p, sc, t, dt, s, 1); t += dt
        ys.append(float(f[4].mean() if model == 1 else f[1].mean())); ts.append(t)
        y = np.array(ys); zc = np.flatnonzero(np.diff(np.sign(y)) != 0)
        if len(zc) >= halfper + 2: break
    T = np.array(ts); y = np.array(ys)
    tz = np.array([T[i] - y[i] * (T[i + 1] - T[i]) / (y[i + 1] - y[i]) for i in zc])[1:]        # drop the first (start-up)
    m = 2 * math.pi / (s * np.diff(tz ** 2)); tm = float(np.mean(tz))
    U19.L().set_tilt(0.0); lib().set_core_late(0.0)
    return dict(m_over_mG=float(np.mean(m) / MG), x_mid=tm / g["tau_A"])


def v4_spectrum(model, ell_g, N=96):
    """V4: two standing waves of known k in a defect-free field: number-weighted <k_lat> from the binned counter (<= 2 %)"""
    MG = set_mass(ell_g); g = geometry(ell_g); s = g["s"]; a = s * g["tau_B"]
    x = np.arange(N); X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k1, k2 = 2 * math.pi * 3 / N, 2 * math.pi * 11 / N; A1, A2 = 0.02, 0.01
    th = A1 * np.cos(k1 * X) + A2 * np.cos(k2 * Y + k2 * Z)
    if model == 1:
        n = np.stack([np.cos(th / 2), np.sin(th / 2), np.zeros_like(th), np.zeros_like(th)]).astype(np.float32)
        f = np.ascontiguousarray(Lt.q_from_n(n).reshape(10, -1), np.float32)
    else:
        f = np.ascontiguousarray(np.stack([V0 * np.cos(th), V0 * np.sin(th)]).reshape(2, -1), np.float32)
    p = np.zeros_like(f)
    tot, h = U19.linear_number(model, N, f, p, a, kbins(N))
    kmean = float(np.sum(h["kn"]) / np.sum(h["n"]))
    kl = lambda k: 2 * math.sin(k / 2)
    kl1, kl2 = kl(k1), math.sqrt(2) * kl(k2)
    w1, w2 = math.sqrt(kl1 ** 2 / a ** 2 + MG ** 2), math.sqrt(kl2 ** 2 / a ** 2 + MG ** 2)
    n1, n2 = w1 * A1 ** 2, w2 * A2 ** 2                                         # standing wave: number ~ w A^2
    expect = (kl1 * n1 + kl2 * n2) / (n1 + n2)
    return dict(k_mean=kmean, expected=expect, rel=kmean / expect - 1, binned_total_rel=float(np.sum(h["n"]) / tot - 1),
                PASS=bool(abs(kmean / expect - 1) <= 0.02))


def v5_kappa(model, ell_g, table, N=256, tau_A=16.0, theta_rms=0.3, seed=7):
    """V5: defect-free random small-angle field with the sparse network correlation length at ell_1 (tau_1/sqrt(xi_1)),
    released early; kappa = sum_{k<k_cut} k n_k / (m_G a_1 N) at tau_B (geometry rescaled to tau_A = 16 so that the box holds
    6-9 correlation lengths; without strings the physics depends only on m_G a l_c and tau/tau_A)"""
    MG = set_mass(ell_g); s = 1.0 / (MR * tau_A); ell_1 = ell_g + math.log(3.0)
    tau_1 = math.sqrt(math.exp(ell_1) / (MR * s)); a1 = s * tau_1; tau_B = END_FACTOR[ell_g] * tau_A
    lc = tau_1 / math.sqrt(XI1_DESIGN)
    rng = np.random.default_rng(seed)
    k = np.fft.fftfreq(N) * 2 * math.pi; K2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2
    fld = np.real(np.fft.ifftn(np.fft.fftn(rng.standard_normal((N, N, N))) * np.exp(-K2 * lc ** 2 / 2)))
    th = theta_rms * fld / np.std(fld)
    if model == 1:
        n = np.stack([np.cos(th / 2), np.sin(th / 2), np.zeros_like(th), np.zeros_like(th)]).astype(np.float32)
        f = np.ascontiguousarray(Lt.q_from_n(n).reshape(10, -1), np.float32)
    else:
        f = np.ascontiguousarray(np.stack([V0 * np.cos(th), V0 * np.sin(th)]).reshape(2, -1), np.float32)
    p = np.zeros_like(f); sc = np.empty_like(f)
    U19.L().set_tilt(U19.eps_of(model)); lib().set_core_late(tau_A)
    tau_i = round(math.sqrt(1 / (s * 30 * MG)) / DTAU) * DTAU; tau_f = round(tau_B / DTAU) * DTAU
    U19.L().evolve(model, N, f, p, sc, tau_i, DTAU, s, int(round((tau_f - tau_i) / DTAU)))
    kb = kbins(N); tot, h = U19.linear_number(model, N, f, p, s * tau_f, kb, **counter_kw(table, tau_f, tau_A, MG, s))
    U19.L().set_tilt(0.0); lib().set_core_late(0.0)
    kc = 0.5 * (kb[1:] + kb[:-1]); soft = kc < 10 * MG * a1
    kappa = float(np.sum(np.array(h["kn"])[soft]) / (MG * a1 * tot))
    return dict(kappa=kappa, l_c_cells=lc, m_G_a1_lc=MG * a1 * lc, N_lin=tot, far_frac_end=U19.site_stats(model, N, f)["far_frac"])


def validate():
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = dict(spec="docs/preregistros/ESPEC_U19b.md", V1={}, V2={}, V3={}, V4={}, V5={}, anharmonic={}, meff={})
    for lg in ELL_GS:
        g = geometry(lg); dense_ell = 3.5 if lg == 3.2 else 3.0
        res["V1"][f"{lg:g}"] = v1(lg); print("V1", lg, res["V1"][f"{lg:g}"], flush=True)
        for model in (1, 0):
            nm = "RP3" if model else "U1"
            tab = meff_table(model, lg); res["meff"][f"{nm}_{lg:g}"] = tab
            rB = tab["ratio"][-1]; chk = []                                        # enmienda 1: V3 = Hessian m_eff vs WKB oscillation
            for xc in (2.0, END_FACTOR[lg] * 0.97):
                w = wkb_mass(model, lg, xc); w["hessian"] = meff_ratio(model, lg, w["x_mid"]); w["rel"] = w["m_over_mG"] / w["hessian"] - 1
                chk.append(w)
            res["V3"][f"{nm}_{lg:g}"] = dict(meff_over_mG_at_tau_B=rB, wkb=chk,
                                             PASS=bool(abs(rB - 1) <= 0.03 and all(abs(w["rel"]) <= 0.003 for w in chk)))
            print("V3 m_eff", nm, lg, tab, flush=True)
            for rel in ("early", "dense"):
                c = control(model, lg, rel, dense_ell, tab); c["rel_vs_ode"] = c["R"] / c["R_ode"] - 1
                c["PASS"] = bool(abs(c["rel_vs_ode"]) <= 0.05 and abs(c["late_number_drift"]) <= 0.01)   # V2 and V3 (number)
                res["V2"][f"{nm}_{lg:g}_{rel}"] = c
                print("V2/V3", nm, lg, rel, {k: round(c[k], 5) for k in ("R", "R_ode", "rel_vs_ode", "late_number_drift")}, c["PASS"], flush=True)
                big = control(model, lg, rel, dense_ell, tab, N=32, theta0=2.5, noise=0.0); small = control(model, lg, rel, dense_ell, tab, N=32, theta0=0.3, noise=0.0)
                res["anharmonic"][f"{nm}_{lg:g}_{rel}"] = dict(ratio=big["R"] / small["R"], R_2p5=big["R"], R_0p3=small["R"])
            e, d = res["anharmonic"][f"{nm}_{lg:g}_early"]["ratio"], res["anharmonic"][f"{nm}_{lg:g}_dense"]["ratio"]
            res["anharmonic"][f"{nm}_{lg:g}_delta"] = d - e
            res["V4"][f"{nm}_{lg:g}"] = v4_spectrum(model, lg); print("V4", nm, lg, res["V4"][f"{nm}_{lg:g}"], flush=True)
            res["V5"][f"{nm}_{lg:g}"] = v5_kappa(model, lg, tab); print("V5", nm, lg, res["V5"][f"{nm}_{lg:g}"], flush=True)
    res["PASS"] = bool(all(v["PASS"] for grp in ("V2", "V3", "V4") for v in res[grp].values())
                       and all(v["PASS"] for lv in res["V1"].values() for v in lv.values()))
    res["runtime_s"] = time.time() - t0
    (OUT / "U19B_VALIDATION.json").write_text(json.dumps(res, indent=1, default=float, allow_nan=False))
    print("VALIDATION", "PASS" if res["PASS"] else "FAIL", flush=True)


# ------------------------------------------------------------------ production
def network(model, seed, protocol, ell_g, N, log, table):
    MG = set_mass(ell_g); g = geometry(ell_g, N / 2); s, tA, tB, l1 = g["s"], g["tau_A"], g["tau_B"], g["ell_1"]
    if protocol == "sparse": ell_i, xi_i = 1.0, 0.5
    else: ell_i, xi_i = (3.5 if ell_g == 3.2 else 3.0), 6.0
    drain = min(3.8, l1 - 0.3)
    tau0 = round(math.sqrt(math.exp(ell_i) * tA) / DTAU) * DTAU
    lc = (tau0 / 18.65) * math.sqrt(Lt.KXI[model] / xi_i)
    f, p = Lt.initial(model, N, lc, seed); sc = np.empty_like(f)
    U19.L().set_tilt(U19.eps_of(model)); U19.L().set_damp(math.sqrt(math.exp(drain) / (MR * s)), MR); lib().set_core_late(tA)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, tB, NSNAP) / DTAU) * DTAU)
    kb = kbins(N); rows, tau, t0 = [], tau0, time.time()
    rows.append(dict(tau=tau0, ell=math.log(MR * s * tau0 * tau0), a=s * tau0, initial=True, **U19.site_stats(model, N, f)))
    for ts in taus:
        nst = int(round((ts - tau) / DTAU))
        if nst > 0: U19.L().evolve(model, N, f, p, sc, tau, DTAU, s, nst); tau += nst * DTAU
        ell = math.log(MR * s * tau * tau)
        r = Lt.measure(model, N, f, p, tau, s, spectra=(ell >= l1 - 0.5), tau_f=tA); r.update(U19.site_stats(model, N, f))
        if ell >= l1 - 0.2:
            kw = counter_kw(table, tau, tA, MG, s); r["m_eff"], r["zchi"] = kw["m_eff"], kw["zchi"]
            r["N_lin"], r["lin_spectrum"] = U19.linear_number(model, N, f, p, s * tau, kb, **kw)
        if "spectrum" in r: r["N_kin"] = U19.kinetic_number(r, s * tau)
        rows.append(r)
        log(f"[U19b {'RP3' if model else 'U1'} lG{ell_g:g} {protocol} s{seed}] ell={ell:.2f} xi={r['xi']:.4f} wall={r['wall_frac']:.4f} "
            f"far={r['far_frac']:.4f} N_lin={r.get('N_lin', float('nan')):.4e} ({time.time() - t0:.0f}s)")
    U19.L().set_tilt(0.0); U19.L().set_damp(0.0, 0.0); lib().set_core_late(0.0)
    return dict(model="RP3" if model else "U1", N=N, tau_A=tA, tau_B=tB, seed=seed, protocol=protocol, ell_G=ell_g, s=s, tau0=tau0, lc=lc,
                xi_init=xi_i, ell_i=ell_i, m_G=MG, eps=U19.eps_of(model), ell_1=l1, a_1=g["a_1"], drain_ell=drain, k_edges=kb.tolist(),
                f2=Lt.F2, m2=Lt.M2[model], threads=Lt.lib().nthreads(), runtime_s=time.time() - t0, rows=rows, spec="docs/preregistros/ESPEC_U19b.md")


def plan():
    """priority order of ESPEC_U19b: RP3 l_G = 3.2, RP3 l_G = 2.4, U(1) l_G = 3.2"""
    return [(1, 3.2, pr, sd) for pr in ("sparse", "dense") for sd in SEEDS] + [(1, 2.4, pr, sd) for pr in ("sparse", "dense") for sd in SEEDS] \
         + [(0, 3.2, pr, sd) for pr in ("sparse", "dense") for sd in SEEDS]


def produce(N=N_PROD, sub="runs", only=None):
    R = OUT / sub; R.mkdir(parents=True, exist_ok=True); log = open(OUT / f"{sub}.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    val = OUT / "U19B_VALIDATION.json"
    if not val.exists() or not json.loads(val.read_text()).get("PASS"): lg("STOP: validation missing or failed (ESPEC_U19b: INVALIDO)"); sys.exit(3)
    tables = json.loads(val.read_text())["meff"]
    for model, ell_g, protocol, seed in (only or plan()):
        tag = f"{'RP3' if model else 'U1'}_N{N}_lG{ell_g:g}_{protocol}_s{seed}"
        if (R / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
        res = network(model, seed, protocol, ell_g, N, lg, tables[f"{'RP3' if model else 'U1'}_{ell_g:g}"])
        (R / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "validate": validate()
    elif cmd == "pilot": produce(N=192, sub="pilot", only=[(1, 2.4, "sparse", 91)])
    elif cmd == "produce": produce(N=int(sys.argv[2]) if len(sys.argv) > 2 else N_PROD)
