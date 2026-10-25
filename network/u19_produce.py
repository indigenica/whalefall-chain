"""U19 (docs/preregistros/ESPEC_U19.md, amendments 1-2): c/eta of the string-wall collapse measured on the lattice.
Explicit breaking (kern.c set_tilt): RP3 V -= eps Q_11 (m_G^2 = eps), U(1) V -= eps phi_1 (m_G^2 = eps / v); m_G = m_r e^-2.8.
    python -B network/u19_produce.py validate     V1 (mass, Minkowski) + V2 (particle counter vs homogeneous ODE, two releases)
                                                  + large-angle controls  ->  network/results/u19/U19_VALIDATION.json
    python -B network/u19_produce.py produce [U1|RP3]   8 networks (RP3, U(1); sparse and dense protocols; seeds 81, 82)
                                                  ->  network/results/u19/runs/<TAG>.json"""
import ctypes, json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results" / "u19"
MR = math.sqrt(0.24); V0 = math.sqrt(0.5); F = math.sqrt(0.5)           # m_r, v = f
MG = MR * math.exp(-2.8); ELL_G = 2.8; ELL_1 = 2.8 + math.log(3.0)
DTAU = 0.25; N_PROD, TAU_END = 320, 160; NSNAP = 40; DRAIN_ELL = 3.8
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def L():
    lib = Lt.lib()
    if not hasattr(lib, "_tilt"):
        lib.set_tilt.argtypes = [ctypes.c_double]; lib.set_tilt.restype = None; lib._tilt = True
    return lib


def eps_of(model): return MG ** 2 if model == 1 else MG ** 2 * V0


# ------------------------------------------------------------------ observables
def homogeneous_field(model, N, theta0, noise, rng):
    """field at angle theta0 from the vacuum (RP3: n = (cos theta0/2, sin theta0/2, 0, 0)) plus smooth noise of rms 'noise'"""
    k = np.fft.fftfreq(N) * 2 * np.pi; K2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2
    sm = lambda: np.real(np.fft.ifftn(np.fft.fftn(rng.standard_normal((N, N, N))) * np.exp(-K2 * 4.0)))
    if model == 1:
        ph = theta0 / 2
        n = np.stack([np.full((N, N, N), math.cos(ph)), np.full((N, N, N), math.sin(ph)), np.zeros((N, N, N)), np.zeros((N, N, N))])
        if noise: n += np.stack([noise * sm() / max(np.std(sm()), 1e-30) for _ in range(4)]) * 0.5
        n /= np.linalg.norm(n, axis=0); f = Lt.q_from_n(n.astype(np.float32))
    else:
        th = theta0 + (noise * sm() / max(np.std(sm()), 1e-30) if noise else 0.0)
        f = np.stack([V0 * np.cos(th) * np.ones((N, N, N)), V0 * np.sin(th) * np.ones((N, N, N))]).astype(np.float32)
    return np.ascontiguousarray(f.reshape(Lt.NCOMP[model], -1), np.float32), np.zeros((Lt.NCOMP[model], N ** 3), np.float32)


def site_stats(model, N, f):
    """wall volume fraction (cos theta < -1/2), far-from-vacuum fraction (cos theta < 0.9), <V_eps - V_min>"""
    if model == 1:
        q0 = f[0].astype(np.float64); cth = 2 * q0 - 1                        # Q_11 = cos^2(phi) = (1 + cos theta)/2
        Veps = eps_of(1) * (1 - q0)
    else:
        r = np.sqrt(f[0].astype(np.float64) ** 2 + f[1].astype(np.float64) ** 2) + 1e-30; cth = f[0] / r
        Veps = eps_of(0) * (V0 - f[0].astype(np.float64))
    return dict(wall_frac=float(np.mean(cth < -0.5)), far_frac=float(np.mean(cth < 0.9)), V_eps=float(np.mean(Veps)))


def linear_number(model, N, f, p, a, kbins=None, m_eff=None, hubble=None, zchi=1.0, zchid=1.0):
    """comoving number of free pseudo-Goldstones, field linearised around the vacuum (amendment 1): a^3 (1/V^2) sum_k
    (|chidot_k|^2 + w_k^2 |chi_k|^2) / (2 w_k), w_k^2 = k_lat^2/a^2 + m_G^2.
    kbins (U19b): comoving k_lat bin edges; then returns (total, dict(n=..., kn=...)) with the per-bin number and k*number.
    m_eff (U19b): Goldstone mass of the evolved theory (late frozen core), used in w_k instead of m_G.
    hubble (U19b): if given, chidot_k -> chidot_k + H (1 + m^2/(2 w_k^2)) chi_k: the first-order WKB number of u = a chi
    (u'' + w_tau^2 u = 0 in radiation; u' + (w_tau'/2w_tau) u = -i w_tau u for the adiabatic solution), whose error is
    O((H/m)^2); without it the counter carries an O(H/m) phase-dependent error.
    zchi, zchid (U19b): canonical normalisation of the Goldstone field and velocity when the vacuum is displaced (late frozen
    core): U(1) chi = r_0 theta (both scaled by r_0/v); RP3 Q = Lambda n n^T + mu (1 - n n^T), chi scaled by Lambda - mu."""
    V = N ** 3; ia3 = 1.0 / a ** 3                                      # physical d/dt = pi / a^3
    if model == 1:
        n = np.empty((4, V), np.float32); L().principal(N, np.ascontiguousarray(f), n, None)
        sg = np.where(n[0] < 0, -1.0, 1.0).astype(np.float32); n *= sg[None]; del sg
        chi, chid = [], []
        for c in (1, 2, 3):                                              # (Qdot n)_c = Qdot_cc n_c + sum_j Qdot_cj n_j (off-diagonal / sqrt2)
            qd = p[c] * n[c]
            for k, (i, j) in enumerate(PAIRS):
                if i == c: qd = qd + p[4 + k] * (n[j] / math.sqrt(2))
                elif j == c: qd = qd + p[4 + k] * (n[i] / math.sqrt(2))
            chi.append((2 * F) * n[c]); chid.append((2 * F * ia3) * qd)
        del n
    else:
        r2 = f[0].astype(np.float64) ** 2 + f[1].astype(np.float64) ** 2; th = np.arctan2(f[1], f[0])
        chi = [V0 * th]; chid = [V0 * ia3 * (f[0].astype(np.float64) * p[1] - f[1].astype(np.float64) * p[0]) / r2]
    k = np.fft.fftfreq(N) * 2 * np.pi; kz = np.fft.rfftfreq(N) * 2 * np.pi
    kl2 = (4 * np.sin(k / 2) ** 2)[:, None, None] + (4 * np.sin(k / 2) ** 2)[None, :, None] + (4 * np.sin(kz / 2) ** 2)[None, None, :]
    w = np.sqrt(kl2 / a ** 2 + (MG if m_eff is None else m_eff) ** 2)
    wz = np.ones(w.shape); wz[:, :, 1:(-1 if N % 2 == 0 else None)] = 2
    tot = 0.0
    if kbins is not None:
        kk = np.sqrt(kl2).ravel(); idx = np.digitize(kk, kbins) - 1; nb = len(kbins) - 1
        hn = np.zeros(nb); hk = np.zeros(nb); ok = (idx >= 0) & (idx < nb)
    for c, cd in zip(chi, chid):
        X = np.fft.rfftn(c.reshape(N, N, N)); Xd = np.fft.rfftn(cd.reshape(N, N, N))
        if zchi != 1.0: X = X * zchi
        if zchid != 1.0: Xd = Xd * zchid
        if hubble is not None: Xd = Xd + hubble * (1 + (MG if m_eff is None else m_eff) ** 2 / (2 * w ** 2)) * X
        nk = wz * (np.abs(Xd) ** 2 + w ** 2 * np.abs(X) ** 2) / (2 * w)
        tot += float(np.sum(nk))
        if kbins is not None:
            v = nk.ravel()
            hn += np.bincount(idx[ok], v[ok], nb); hk += np.bincount(idx[ok], (kk * v)[ok], nb)
    if kbins is None: return a ** 3 * tot / V ** 2
    return a ** 3 * tot / V ** 2, dict(n=(a ** 3 * hn / V ** 2).tolist(), kn=(a ** 3 * hk / V ** 2).tolist())


def kinetic_number(row, a):
    """comoving number from the kinetic spectrum of measure(): a^3 sum 2 dE / w(k_centre)"""
    sp = row["spectrum"]; e = np.array(sp["k_edges"]); kc = 0.5 * (e[1:] + e[:-1])
    w = np.sqrt(kc ** 2 / a ** 2 + MG ** 2)
    return float(a ** 3 * np.sum(2 * np.array(sp["dE"]) / w))


# ------------------------------------------------------------------ validation
def v1_mass():
    out = {}
    for model in (1, 0):
        L().set_tilt(eps_of(model)); N = 8; dt = 0.05; rng = np.random.default_rng(0)
        f, p = homogeneous_field(model, N, 0.02, 0.0, rng); sc = np.empty_like(f)
        T = 2 * math.pi / MG; nper = 4; steps_per = 20; nsamp = int(nper * T / dt / steps_per)
        ys, ts, t = [], [], 0.0
        for _ in range(nsamp):
            L().evolve(model, N, f, p, sc, 1.0, dt, 0.0, steps_per); t += steps_per * dt
            ys.append(float(f[4].mean() if model == 1 else f[1].mean())); ts.append(t)
        ys = np.array(ys); ts = np.array(ts); zc = ts[1:][np.diff(np.sign(ys)) != 0]
        w = math.pi / np.mean(np.diff(zc)); L().set_tilt(0.0)
        out["RP3" if model else "U1"] = dict(omega2=w ** 2, m_G2=MG ** 2, rel=w ** 2 / MG ** 2 - 1, PASS=bool(abs(w ** 2 / MG ** 2 - 1) <= 0.01))
    return out


def ode_number(theta0, tau_i, tau_f, s):
    """homogeneous theta'' + (2/tau) theta' + a^2 m^2 sin theta = 0 (a = s tau), released at tau_i; comoving number at tau_f"""
    sol = solve_ivp(lambda t, y: [y[1], -2 / t * y[1] - (s * t) ** 2 * MG ** 2 * math.sin(y[0])], (tau_i, tau_f), [theta0, 0.0],
                    rtol=1e-10, atol=1e-12, dense_output=True)
    th, thp = sol.y[0, -1], sol.y[1, -1]; a = s * tau_f
    rho = 0.5 * F ** 2 * (thp / a) ** 2 + MG ** 2 * F ** 2 * (1 - math.cos(th))
    return a ** 3 * rho / MG, MG ** 2 * F ** 2 * (1 - math.cos(theta0))


def control(model, N, theta0, noise, release, seed=5):
    """homogeneous control in radiation; release 'early' (H = 30 m_G) or 'dense' (frozen until ell_i = 3.5).
    Returns R = m_G N_lin(final) / (<V_eps - V_min>_initial a_1^3) and the ODE prediction of the same ratio."""
    L().set_tilt(eps_of(model)); s = 1.0 / (MR * 64)
    tau_i = math.sqrt(1 / (s * 30 * MG)) if release == "early" else math.sqrt(math.exp(3.5) / (MR * s))
    tau_i = round(tau_i / DTAU) * DTAU; tau_f = round(math.sqrt(20 / (MG * s)) / DTAU) * DTAU
    a1 = s * math.sqrt(3 / (MG * s))
    rng = np.random.default_rng(seed); f, p = homogeneous_field(model, N, theta0, noise, rng); sc = np.empty_like(f)
    V_i = site_stats(model, N, f)["V_eps"]
    L().evolve(model, N, f, p, sc, tau_i, DTAU, s, int(round((tau_f - tau_i) / DTAU)))
    Nl = linear_number(model, N, f, p, s * tau_f); L().set_tilt(0.0)
    No, Vo = ode_number(theta0, tau_i, tau_f, s)
    return dict(R=MG * Nl / (V_i * a1 ** 3), R_ode=MG * No / (Vo * a1 ** 3), N_lin=Nl, N_ode=No, V_i=V_i, V_ode=Vo,
                tau_i=tau_i, tau_f=tau_f, far_frac_final=site_stats(model, N, f)["far_frac"])


def validate():
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = dict(spec="docs/preregistros/ESPEC_U19.md", m_G=MG, V1=v1_mass(), V2={}, anharmonic={})
    print("V1", res["V1"], flush=True)
    for model in (1, 0):
        nm = "RP3" if model else "U1"
        for rel in ("early", "dense"):
            c = control(model, 128, 0.3, 0.01, rel); c["rel_vs_ode"] = c["R"] / c["R_ode"] - 1; c["PASS"] = bool(abs(c["rel_vs_ode"]) <= 0.05)
            res["V2"][f"{nm}_{rel}"] = c; print("V2", nm, rel, {k: c[k] for k in ("R", "R_ode", "rel_vs_ode", "PASS")}, flush=True)
            big = control(model, 32, 2.5, 0.0, rel); small = control(model, 32, 0.3, 0.0, rel)
            res["anharmonic"][f"{nm}_{rel}"] = dict(R_2p5=big["R"], R_0p3=small["R"], ratio=big["R"] / small["R"], R_ode_2p5=big["R_ode"], R_ode_0p3=small["R_ode"])
            print("anh", nm, rel, res["anharmonic"][f"{nm}_{rel}"], flush=True)
    for nm in ("RP3", "U1"):
        e, d = res["anharmonic"][f"{nm}_early"]["ratio"], res["anharmonic"][f"{nm}_dense"]["ratio"]
        res["anharmonic"][f"{nm}_delta"] = d - e
    res["PASS"] = bool(all(v["PASS"] for v in res["V1"].values()) and all(v["PASS"] for v in res["V2"].values()))
    res["runtime_s"] = time.time() - t0
    (OUT / "U19_VALIDATION.json").write_text(json.dumps(res, indent=1, default=float, allow_nan=False))
    print("VALIDATION", "PASS" if res["PASS"] else "FAIL", flush=True)
    return res


# ------------------------------------------------------------------ production
def network(model, seed, protocol, log):
    N, te = N_PROD, TAU_END; m = MR; s = 1.0 / (m * te)
    ell_i, xi_i = (1.0, 0.5) if protocol == "sparse" else (3.5, 6.0)
    tau0 = round(math.sqrt(math.exp(ell_i) * te) / DTAU) * DTAU
    lc = (tau0 / 18.65) * math.sqrt(Lt.KXI[model] / xi_i)
    f, p = Lt.initial(model, N, lc, seed); sc = np.empty_like(f)
    L().set_tilt(eps_of(model)); tau_d = math.sqrt(math.exp(DRAIN_ELL) / (m * s)); L().set_damp(tau_d, m)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, te, NSNAP) / DTAU) * DTAU)
    rows, tau, t0 = [], tau0, time.time()
    st0 = site_stats(model, N, f); rows.append(dict(tau=tau0, ell=math.log(m * s * tau0 * tau0), a=s * tau0, initial=True, **st0))
    for ts in taus:
        nst = int(round((ts - tau) / DTAU))
        if nst > 0: L().evolve(model, N, f, p, sc, tau, DTAU, s, nst); tau += nst * DTAU
        ell = math.log(m * s * tau * tau)
        r = Lt.measure(model, N, f, p, tau, s, spectra=(ell >= 3.5)); r.update(site_stats(model, N, f))
        if ell >= ELL_1 - 0.2: r["N_lin"] = linear_number(model, N, f, p, s * tau)
        if "spectrum" in r: r["N_kin"] = kinetic_number(r, s * tau)
        rows.append(r)
        log(f"[U19 {'RP3' if model else 'U1'} {protocol} s{seed}] ell={ell:.2f} xi={r['xi']:.3f} wall={r['wall_frac']:.4f} "
            f"far={r['far_frac']:.4f} N_lin={r.get('N_lin', float('nan')):.4e} ({time.time() - t0:.0f}s)")
    L().set_tilt(0.0); L().set_damp(0.0, 0.0)
    return dict(model="RP3" if model else "U1", N=N, tau_end=te, seed=seed, protocol=protocol, s=s, tau0=tau0, lc=lc, xi_init=xi_i, ell_i=ell_i,
                m_G=MG, eps=eps_of(model), ell_1=ELL_1, drain_ell=DRAIN_ELL, f2=Lt.F2, m2=Lt.M2[model], runtime_s=time.time() - t0, rows=rows,
                spec="docs/preregistros/ESPEC_U19.md")


def produce(models=(1, 0)):
    R = OUT / "runs"; R.mkdir(parents=True, exist_ok=True); log = open(OUT / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    val = json.loads((OUT / "U19_VALIDATION.json").read_text())
    if not val.get("PASS"): lg("STOP: validation failed (ESPEC_U19: INVALIDO)"); sys.exit(3)
    for model in models:
        for protocol in ("sparse", "dense"):
            for seed in (81, 82):
                tag = f"{'RP3' if model else 'U1'}_N{N_PROD}_{protocol}_s{seed}"
                if (R / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
                res = network(model, seed, protocol, lg)
                (R / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    if sys.argv[1] == "validate": validate()
    else: produce({"U1": (0,), "RP3": (1,)}.get(sys.argv[2], (1, 0)) if len(sys.argv) > 2 else (1, 0))   # RAM: RP3 runs one at a time
