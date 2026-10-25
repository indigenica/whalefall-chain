"""U19c (docs/preregistros/ESPEC_U19c.md): the U19b collapse networks in more masses, with the l_G = 3.2 protocol carried
relative to l_1: N = 2 tau_A with l(tau_A) - l_1 = 1.247, i.e. l_G = 3.2 + ln(N/512); tau_B = 3 tau_A; sparse xi_i = 0.5 at
l_i = 1, dense xi_i = 6 at l_i = l_1 - 0.8; drain to l_1 - 0.5; 48 snapshots; seeds 91, 92; late frozen core from tau_A.
Everything else (validation V1-V5, counter, measurements) is u19b_produce.py, with each mass on its own geometry.
The high mass (640^3, or 576^3 by enmienda 1; 768^3 by enmienda 2) uses evolve_lowmem (no force buffer, bit-identical to evolve:
V0'), in kern.c or, with WHALEFALL_BACKEND=jax, in network/jax_kernel.py on a GPU.
The four RP3 networks of U19b with l_G = 3.2 at 512^3 are the l_G = 3.2 point (already this protocol); not produced again.

    python -B network/u19c_produce.py validate            V1-V5 for every new mass (cloud) -> network/results/u19c/U19C_VALIDATION.json
    python -B network/u19c_produce.py produce N           the 4 RP3 networks of that mass, skip-if-present
                                                          (cloud: N = 256, 384; WSL: N = 640, or 576 if WSL has < 56 GB;
                                                          enmienda 2: 576 and 768 on a Colab A100 with WHALEFALL_BACKEND=jax)
    python -B network/u19c_produce.py control             enmienda 2 (b): the 256^3 control network on this backend (A100)
    -> network/results/u19c/runs/<TAG>.json and network/results/u19c/runs.log (or under U19C_OUT)"""
import json, math, os, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt                                                          # noqa: E402
import u19_produce as U19                                                     # noqa: E402
import u19b_produce as B                                                      # noqa: E402

HERE = Path(__file__).resolve().parent; OUT = Path(os.environ.get("U19C_OUT", HERE / "results" / "u19c"))   # U19C_OUT: e.g. Drive on Colab
MR, DTAU = U19.MR, U19.DTAU
MASSES = {N: round(3.2 + math.log(N / 512), 3) for N in (256, 384, 576, 640, 768)}  # 2.507, 2.912, 3.318, 3.423, 3.605 (enmienda 2)
LOWMEM_N = 576                                                                     # evolve_lowmem from this N up (enmienda 1)
SEEDS, NSNAP = B.SEEDS, B.NSNAP

for _N, _lg in MASSES.items():                    # each new mass on its own geometry: tau_A = N/2, tau_B = 3 tau_A
    B.END_FACTOR[_lg] = 3.0
_TAU_A = {lg: N / 2 for N, lg in MASSES.items()}; _geometry = B.geometry
B.geometry = lambda ell_g, tau_A=None: _geometry(ell_g, (_TAU_A.get(ell_g, B.N_PROD / 2) if tau_A is None else tau_A))


def protocol(ell_g):
    """ESPEC_U19c: the l_G = 3.2 protocol relative to l_1"""
    l1 = ell_g + math.log(3.0)
    return dict(ell_1=l1, sparse=(1.0, 0.5), dense=(l1 - 0.8, 6.0), drain_ell=l1 - 0.5)


def validate(masses=tuple(MASSES)):
    """V1-V5 of u19b_produce.validate, for each new mass on its own geometry; dense control released at l_1 - 0.8"""
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time(); path = OUT / "U19C_VALIDATION.json"
    res = json.loads(path.read_text()) if path.exists() else dict(spec="docs/preregistros/ESPEC_U19c.md", masses={})
    for N in masses:
        lg = MASSES[N]; key = f"{lg:g}"
        if res["masses"].get(key, {}).get("done"): print("validated", N, lg); continue
        r = dict(N=N, ell_G=lg, V1=B.v1(lg), V2={}, V3={}, V4={}, V5={}, anharmonic={}, meff={}); dense_ell = protocol(lg)["dense"][0]
        print("V1", lg, r["V1"], flush=True)
        for model in (1, 0):
            nm = "RP3" if model else "U1"; tab = B.meff_table(model, lg); r["meff"][nm] = tab
            rB = tab["ratio"][-1]; chk = []
            for xc in (2.0, B.END_FACTOR[lg] * 0.97):
                w = B.wkb_mass(model, lg, xc); w["hessian"] = B.meff_ratio(model, lg, w["x_mid"]); w["rel"] = w["m_over_mG"] / w["hessian"] - 1
                chk.append(w)
            r["V3"][nm] = dict(meff_over_mG_at_tau_B=rB, wkb=chk, PASS=bool(abs(rB - 1) <= 0.03 and all(abs(w["rel"]) <= 0.003 for w in chk)))
            for rel in ("early", "dense"):
                c = B.control(model, lg, rel, dense_ell, tab); c["rel_vs_ode"] = c["R"] / c["R_ode"] - 1
                c["PASS"] = bool(abs(c["rel_vs_ode"]) <= 0.05 and abs(c["late_number_drift"]) <= 0.01)
                r["V2"][f"{nm}_{rel}"] = c
                print("V2/V3", nm, lg, rel, {k: round(c[k], 5) for k in ("R", "R_ode", "rel_vs_ode", "late_number_drift")}, c["PASS"], flush=True)
                big = B.control(model, lg, rel, dense_ell, tab, N=32, theta0=2.5, noise=0.0)
                small = B.control(model, lg, rel, dense_ell, tab, N=32, theta0=0.3, noise=0.0)
                r["anharmonic"][f"{nm}_{rel}"] = dict(ratio=big["R"] / small["R"], R_2p5=big["R"], R_0p3=small["R"])
            r["anharmonic"][f"{nm}_delta"] = r["anharmonic"][f"{nm}_dense"]["ratio"] - r["anharmonic"][f"{nm}_early"]["ratio"]
            r["V4"][nm] = B.v4_spectrum(model, lg); r["V5"][nm] = B.v5_kappa(model, lg, tab)
            print("V4/V5", nm, lg, r["V4"][nm], r["V5"][nm], flush=True)
        r["PASS"] = bool(all(v["PASS"] for grp in ("V2", "V3", "V4") for v in r[grp].values()) and all(v["PASS"] for v in r["V1"].values()))
        r["done"] = True; res["masses"][key] = r
        path.write_text(json.dumps(res, indent=1, default=float, allow_nan=False))
        print("VALIDATION", N, lg, "PASS" if r["PASS"] else "FAIL", flush=True)
    res["runtime_s"] = res.get("runtime_s", 0) + time.time() - t0
    path.write_text(json.dumps(res, indent=1, default=float, allow_nan=False))


def network(model, seed, prot, N, log, table):
    """u19b_produce.network with the relative protocol; evolve_lowmem (no force buffer) for N >= LOWMEM_N"""
    lg = MASSES[N]; MG = B.set_mass(lg); g = B.geometry(lg, N / 2); s, tA, tB, l1 = g["s"], g["tau_A"], g["tau_B"], g["ell_1"]
    P = protocol(lg); ell_i, xi_i = P[prot]; drain = P["drain_ell"]; low = N >= LOWMEM_N
    tau0 = round(math.sqrt(math.exp(ell_i) * tA) / DTAU) * DTAU
    lc = (tau0 / 18.65) * math.sqrt(Lt.KXI[model] / xi_i)
    f, p = Lt.initial(model, N, lc, seed); sc = None if low else np.empty_like(f)
    L = U19.L(); L.set_tilt(U19.eps_of(model)); L.set_damp(math.sqrt(math.exp(drain) / (MR * s)), MR); L.set_core_late(tA)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, tB, NSNAP) / DTAU) * DTAU)
    kb = B.kbins(N); rows, tau, t0 = [], tau0, time.time()
    rows.append(dict(tau=tau0, ell=math.log(MR * s * tau0 * tau0), a=s * tau0, initial=True, **U19.site_stats(model, N, f)))
    for ts in taus:
        nst = int(round((ts - tau) / DTAU))
        if nst > 0:
            if low: L.evolve_lowmem(model, N, f, p, tau, DTAU, s, nst)
            else: L.evolve(model, N, f, p, sc, tau, DTAU, s, nst)
            tau += nst * DTAU
        ell = math.log(MR * s * tau * tau)
        r = Lt.measure(model, N, f, p, tau, s, spectra=(ell >= l1 - 0.5), tau_f=tA); r.update(U19.site_stats(model, N, f))
        if ell >= l1 - 0.2:
            kw = B.counter_kw(table, tau, tA, MG, s); r["m_eff"], r["zchi"] = kw["m_eff"], kw["zchi"]
            r["N_lin"], r["lin_spectrum"] = U19.linear_number(model, N, f, p, s * tau, kb, **kw)
        if "spectrum" in r: r["N_kin"] = U19.kinetic_number(r, s * tau)
        rows.append(r)
        log(f"[U19c {'RP3' if model else 'U1'} N{N} lG{lg:g} {prot} s{seed}] ell={ell:.2f} xi={r['xi']:.4f} wall={r['wall_frac']:.4f} "
            f"far={r['far_frac']:.4f} N_lin={r.get('N_lin', float('nan')):.4e} ({time.time() - t0:.0f}s)")
    L.set_tilt(0.0); L.set_damp(0.0, 0.0); L.set_core_late(0.0)
    return dict(model="RP3" if model else "U1", N=N, tau_A=tA, tau_B=tB, seed=seed, protocol=prot, ell_G=lg, s=s, tau0=tau0, lc=lc,
                xi_init=xi_i, ell_i=ell_i, m_G=MG, eps=U19.eps_of(model), ell_1=l1, a_1=g["a_1"], drain_ell=drain, k_edges=kb.tolist(),
                f2=Lt.F2, m2=Lt.M2[model], threads=Lt.lib().nthreads(), lowmem=low, backend=backend(), runtime_s=time.time() - t0, rows=rows,
                spec="docs/preregistros/ESPEC_U19c.md")


def backend():
    """kern.c, or the JAX kernel (WHALEFALL_BACKEND=jax, enmienda 2) with its device"""
    if os.environ.get("WHALEFALL_BACKEND") != "jax": return "kern.c"
    import jax
    return f"jax {jax.__version__} on {jax.devices()[0].device_kind}"


def produce(N):
    R = OUT / "runs"; R.mkdir(parents=True, exist_ok=True); log = open(OUT / "runs.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    key = f"{MASSES[N]:g}"; val = OUT / "U19C_VALIDATION.json"
    v = json.loads(val.read_text())["masses"].get(key) if val.exists() else None
    if not v or not v.get("PASS"): lg(f"STOP: validation of l_G = {key} missing or failed (ESPEC_U19c: that mass is INVALIDA)"); sys.exit(3)
    lg(f"START N={N} l_G={key} {time.strftime('%Y-%m-%d %H:%M:%S')} threads={Lt.lib().nthreads()} lowmem={N >= LOWMEM_N} backend={backend()}")
    for prot in ("sparse", "dense"):
        for seed in SEEDS:
            tag = f"RP3_N{N}_lG{key}_{prot}_s{seed}"
            if (R / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
            res = network(1, seed, prot, N, lg, v["meff"]["RP3"])
            (R / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


def control():
    """enmienda 2 (b): the control network RP3_N256_lG2.507_sparse_s91 again, on this backend -> OUT/control/<TAG>.json"""
    R = OUT / "control"; R.mkdir(parents=True, exist_ok=True); tag = "RP3_N256_lG2.507_sparse_s91"
    v = json.loads((OUT / "U19C_VALIDATION.json").read_text())["masses"]["2.507"]
    res = network(1, 91, "sparse", 256, lambda m: print(m, flush=True), v["meff"]["RP3"])
    (R / f"{tag}.json").write_text(json.dumps(res, default=float)); print(f"DONE control {tag} {res['backend']} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "validate": validate(tuple(int(x) for x in sys.argv[2:]) or tuple(MASSES))
    elif cmd == "produce": produce(int(sys.argv[2]))
    elif cmd == "control": control()
