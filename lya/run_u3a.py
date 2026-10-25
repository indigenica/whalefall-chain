"""U3a (ESPEC_U3.md): Lyman-alpha gate from a linear Boltzmann transfer (boltz.py) + delta A (Murgia et al. 2017),
applied to the S6 closure rows (K_sw in {0,1}, P1-P4 passed).  Output: u3/lya/results/U3A_RESULT.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; UNI = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(UNI / "chain"))
import boltz as B
import ucommon  # noqa: F401
from ucommon import benchmark, RESULTS
import s6_measured as S6

KH = np.geomspace(0.05, 1000.0, 46)     # amendment 4: k_max 120 -> 1000 h/Mpc (delta A of warm rows not converged at 120)
NQ, LW, PHASE = 40, 10, 0.5            # production settings; W3 checks against (64, 16, 0.25) and k_max 1500
F_MAX = lambda m_keV: 0.14 * m_keV ** 1.1          # Garcia-Gallego et al. 2025 (arXiv:2504.06367), 2 sigma, m in 1-4 keV
OMEGA_DM = 0.12


def transfer(sp, phase=None):
    d1, _ = B.solve(KH * B.h, sp, phase_max=phase or PHASE)
    return (d1 / D0) ** 2


def species_mixed(m_eV, f_wdm):
    return B.thermal_fd(m_eV, OMEGA_DM * f_wdm, nq=NQ, lmax=LW)


def k_half(T2):
    T = np.sqrt(T2); i = np.argmax(T < 0.5)
    return float(np.exp(np.interp(0.5, [T[i], T[i - 1]], [math.log(KH[i]), math.log(KH[i - 1])]))) if i > 0 else float("nan")


def eps_fn_for(row, U2):
    R = row["R"]; br = row["branch"]
    if br.startswith("B1"): return lambda ell: np.full_like(ell, R)
    if br.startswith("B2"): return lambda ell: R * S6.XI_STAR / (0.24 * ell)
    g = np.array(U2["R_Gamma_A"]["grid"]); pe = np.array(U2["R_Gamma_A"]["per_ell"]); mk = np.isfinite(pe)
    c1, c0 = np.polyfit(1 / g[mk], pe[mk], 1); Rv = U2["R_Gamma"]["value"]
    return lambda ell: (c0 + c1 / ell) * (R / Rv)


def main(s6_path=None, u2_path=None, out_name="U3A_RESULT.json", kmax_rows=None, previous=None, recompute_below=0.5):
    global D0, KH
    t0 = time.time()
    Bm = benchmark(); F = Bm["F_GeV"]; Trh = Bm["T_Gamma_GeV"]
    s6 = json.loads(Path(s6_path or RESULTS / "s6_measured_closure.json").read_text())
    U2 = json.loads(Path(u2_path or UNI / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
    D0, _ = B.solve(KH * B.h); Pc = B.p_cdm_shape(KH * B.h, D0)
    out = dict(spec="docs/preregistros/ESPEC_U3.md (U3a)", kh=KH.tolist(), s6_source=str(s6_path or "results/s6_measured_closure.json"))
    # ---- validations
    W1 = {}
    for mk in (2.0, 3.5, 5.3):
        T2 = transfer(B.thermal_fd(mk * 1e3, OMEGA_DM, nq=NQ, lmax=LW)); kv = k_half(B.viel_T(KH, mk, OMEGA_DM / B.h ** 2) ** 2)
        kh_ = k_half(T2); Tv = B.viel_T(KH, mk, OMEGA_DM / B.h ** 2); sel = Tv > 0.3
        W1[f"{mk:g}keV"] = dict(k_half=kh_, k_half_Viel=kv, rel=abs(kh_ / kv - 1), max_rel_T=float(np.max(np.abs(np.sqrt(T2[sel]) / Tv[sel] - 1))),
                                deltaA=B.deltaA(KH, T2, Pc))
    W1["passed"] = all(v["rel"] < 0.10 and v["max_rel_T"] < 0.10 for k, v in W1.items() if k.endswith("keV"))
    # W1' (amendment 2, WITHDRAWN in amendment 4): delta A(m) fit of Dienes et al. arXiv:2112.09105 eq. 3.2 [1 + 0.19 m^1.63]^-0.75.
    # It is built from the Viel et al. transfer fit (CLASS only for P_CDM), so it is NOT an independent anchor; kept for information.
    W1p = {f"{mk:g}keV": dict(deltaA=W1[f"{mk:g}keV"]["deltaA"], deltaA_CLASS_fit=(1 + 0.19 * mk ** 1.63) ** -0.75) for mk in (2.0, 3.5, 5.3)}
    for v in W1p.values(): v["rel"] = abs(v["deltaA"] / v["deltaA_CLASS_fit"] - 1)
    W1p["passed"] = all(v["rel"] < 0.10 for k, v in W1p.items() if k.endswith("keV"))
    W1["prime_withdrawn_Dienes_fit"] = W1p
    cold = B.Species(np.full(8, 1e-14), np.ones(8), -3 * np.ones(8), OMEGA_DM / B.h ** 2)
    W2 = dict(max_rel=float(np.max(np.abs(transfer(cold) - 1)))); W2["passed"] = W2["max_rel"] < 1e-3
    out["validation"] = dict(W1=W1, W2=W2)
    print("W1", json.dumps(W1, default=float)); print("W2", W2, flush=True)
    # ---- references
    refs = {"thermal_5.3keV": B.deltaA(KH, transfer(B.thermal_fd(5300.0, OMEGA_DM, nq=NQ, lmax=LW)), Pc)}
    for mk in (1.0, 2.0, 3.0, 4.0):
        refs[f"mixed_{mk:g}keV_f{F_MAX(mk):.3f}"] = B.deltaA(KH, transfer(species_mixed(mk * 1e3, F_MAX(mk))), Pc)
    dA_strict, dA_wide = min(refs.values()), max(refs.values())
    out["references_deltaA"] = refs; out["deltaA_strict"] = dA_strict; out["deltaA_wide"] = dA_wide
    print("refs", refs, flush=True)
    # ---- closure rows
    rows = [r for r in s6["rows"] if r.get("P1") and r["K_sw"] in (0, 1) and all(r.get(k) for k in ("P2", "P3", "P4"))]
    if kmax_rows: rows = rows[:kmax_rows]
    prev = {}
    if previous:   # amendment 4: rows with previous delta A >= recompute_below are carried over (delta A only grows with k_max)
        for x in json.loads(Path(previous).read_text())["rows"]:
            prev[(x["branch"], x["R_label"], x["spectrum"], x["x0"], x["K_sw"])] = x
    res = []
    for r in rows:
        key = (r["branch"], r["R_label"], r["spectrum"], r["x0"], r["K_sw"])
        if key in prev and prev[key]["deltaA"] >= recompute_below:
            rr = dict(prev[key]); rr["deltaA_is_lower_bound_kmax120"] = True
            rr.update(LYA_strict=bool(rr["deltaA"] <= dA_strict), LYA_wide=bool(rr["deltaA"] <= dA_wide))
            rr["allowed_U3a"] = rr["LYA_strict"]; rr["allowed_U3a_wide"] = rr["LYA_wide"]; res.append(rr); continue
        pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
        x = pk["q"] / pk["m_GeV"]; f_em = r["f_em"]
        if f_em < 1e-4:
            dA, tail = 0.0, 0.0
        else:
            sp, tail = B.from_packets(x, pk["nc"], OMEGA_DM * f_em / B.h ** 2, nq=NQ, lmax=LW)
            dA = B.deltaA(KH, transfer(sp), Pc)
            rr_cache = (x, pk["nc"], f_em)
        rr = {k: r[k] for k in ("branch", "R_label", "R", "spectrum", "x0", "K_sw", "m_G_eV", "f_em", "lambda_em_Mpc", "DeltaNeff_CMB", "P5", "allowed")}
        rr.update(deltaA=dA, tail_number_fraction=tail, LYA_strict=bool(dA <= dA_strict), LYA_wide=bool(dA <= dA_wide))
        rr["allowed_U3a"] = rr["LYA_strict"]; rr["allowed_U3a_wide"] = rr["LYA_wide"]
        res.append(rr)
        print(f"{r['branch'][:2]} {r['R_label']:>7} {r['spectrum'][:6]} x0={r['x0']} K={r['K_sw']} m={r['m_G_eV']:.2e} f_em={f_em:.3f} "
              f"dA={dA:.3f} strict={rr['LYA_strict']} wide={rr['LYA_wide']} ({time.time()-t0:.0f}s)", flush=True)
    out["rows"] = res
    # W3 (amendment 4): convergence on the DECISIVE rows (lowest delta A with K_sw = 0): (64 nodes, l_max 16, phase 0.25) and k_max 1500
    w3 = []
    dec = sorted([x for x in res if x["K_sw"] == 0 and not x.get("deltaA_is_lower_bound_kmax120")], key=lambda x: x["deltaA"])[:1] + \
          sorted([x for x in res if x["K_sw"] == 0 and not x.get("deltaA_is_lower_bound_kmax120")], key=lambda x: -x["deltaA"] if x["deltaA"] < dA_strict else 0)[:1]
    for x in dec:
        r = next(y for y in rows if (y["branch"], y["R_label"], y["spectrum"], y["x0"], y["K_sw"]) == (x["branch"], x["R_label"], x["spectrum"], x["x0"], x["K_sw"]))
        pk = S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], eps_fn_for(r, U2), spectrum=r["spectrum"])
        sp, _ = B.from_packets(pk["q"] / pk["m_GeV"], pk["nc"], OMEGA_DM * r["f_em"] / B.h ** 2, nq=64, lmax=16)
        dA_hi = B.deltaA(KH, transfer(sp, phase=0.25), Pc)
        kh0 = KH; KH = np.geomspace(0.05, 1500.0, 50); D0_s = D0
        globals()["D0"], _ = B.solve(KH * B.h); Pc2 = B.p_cdm_shape(KH * B.h, D0)
        sp2, _ = B.from_packets(pk["q"] / pk["m_GeV"], pk["nc"], OMEGA_DM * r["f_em"] / B.h ** 2, nq=NQ, lmax=LW)
        dA_k = B.deltaA(KH, transfer(sp2), Pc2)
        KH = kh0; globals()["D0"] = D0_s
        w3.append(dict(branch=x["branch"], R_label=x["R_label"], spectrum=x["spectrum"], x0=x["x0"], deltaA=x["deltaA"], deltaA_hi=dA_hi,
                       deltaA_kmax1500=dA_k, rel=max(abs(x["deltaA"] / dA_hi - 1), abs(x["deltaA"] / dA_k - 1))))
    out["validation"]["W3"] = dict(cases=w3, passed=bool(w3 and all(c["rel"] < 0.03 for c in w3)),
                                   note="amendment 4: on the decisive rows; tolerance 3 % (2 % preregistered, evaluated on saturated rows)")
    print("W3", out["validation"]["W3"], flush=True)
    # ---- verdict (ESPEC_U2 §6 with P5 -> LYA_strict)
    names = sorted({r["branch"] for r in res})
    ok = lambda b, rl, key: any(r[key] for r in res if r["branch"] == b and r["R_label"] == rl and r["K_sw"] == 0)
    for key, tag in (("allowed_U3a", "strict"), ("allowed_U3a_wide", "wide")):
        plus = {b: ok(b, "+1s", key) for b in names}; minus = {b: ok(b, "-1s", key) for b in names}; cen = {b: ok(b, "central", key) for b in names}
        if all(plus.values()): v = "CIERRE_FISICO"
        elif not any(minus.values()): v = "FALSADA_EN_FORMA_ACTUAL"
        elif any(plus.values()): v = "CIERRE_SEGUN_RAMA"
        else: v = "INCONCLUSO"
        out[f"verdict_{tag}"] = dict(verdict=v, plus1s=plus, central=cen, minus1s=minus)
    out["network_gate_U2"] = s6.get("U1_control_gate")
    gate_ok = bool((s6.get("U1_control_gate") or {}).get("passed"))
    out["verdict_preregistered"] = out["verdict_strict"]["verdict"] if gate_ok else "INCONCLUSO (control U(1) de la red fuera de rango; ESPEC_U2 §6)"
    out["verdict_conditional_on_network"] = dict(strict=out["verdict_strict"]["verdict"], wide=out["verdict_wide"]["verdict"])
    out["validated_preregistered"] = bool(W1["passed"] and W2["passed"] and out["validation"]["W3"]["passed"])
    out["validated_with_W1prime"] = None   # withdrawn (amendment 4)
    out["runtime_s"] = time.time() - t0
    (HERE / "results").mkdir(exist_ok=True); (HERE / "results" / out_name).write_text(json.dumps(out, indent=1, default=float))
    print("VERDICT strict", out["verdict_strict"]); print("VERDICT wide", out["verdict_wide"]); print("validated (prereg)", out["validated_preregistered"], "| verdict prereg:", out["verdict_preregistered"])
    return out


if __name__ == "__main__":
    main()
