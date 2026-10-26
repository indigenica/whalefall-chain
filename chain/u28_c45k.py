"""U28 = C45k (docs/preregistros/ESPEC_U28.md, amendments 1-6): multi-bounce coherence of the C45 parents with an
inhomogeneous chi, in the dust-comoving NPQT solver chain/u28_comoving.py; eternal and direct children, pole plateau.
    python -B chain/u28_c45k.py gates        ->  results/u28_runs/gates.json   (V0-V4)
    python -B chain/u28_c45k.py produce      ->  results/u28_runs/<child>_<family>_<LH>_<eps>.json (skip-if-present)
    python -B chain/u28_c45k.py analyze      ->  results/u28_c45k.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from ucommon import ROOT, benchmark, load, save
sys.path.insert(0, str(ROOT / "components" / "c45b_hom"))
import u28_comoving as C
import u27_direct_children as U27

import os
RUNS = Path(os.environ.get("U28_RUNS", ROOT / "results" / "u28_runs"))   # U28_RUNS: temporary folder while producing
GRID = (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
DOWN = (3e-6, 1e-6, 3e-7, 1e-7)
UP = (3e-2, 1e-1)
ZETA0 = 4.6e-5
X_HI = 0.95
N_PROD = 256                                   # fixed by V3 (amendment 6)
N_SU = 24                                      # separate-universe sample cells in the region (amendment 6)


def children():
    """the U27 pole grid, its classes, the dust tolerance of C45j per child, and the deterministic U28 selection (sec. 4)"""
    pole, et = U27.pole_grid()
    wP, CP = U27.classes(pole, et, "eternal")
    keys = {}
    for c in wP:
        a, b = U27.parent(c["lam"], c["sep"])
        for LH in (0.3, 1.0):
            keys.setdefault((round(a, 6), round(b, 12), LH), [a, b, LH])
    ks = sorted(keys); vals = U27.eps_table([keys[k] for k in ks]); E = {k: vals[tuple(keys[k])] for k in ks}
    out = {}
    for cls in ("E", "D"):
        rows = []
        for c in CP[cls]:
            a, b = U27.parent(c["lam"], c["sep"])
            rec = dict(c, AH=a, alphaH2=b, cls=cls, band=bool(10 <= c["chi_abs_max"] <= 50))
            for LH in (0.3, 1.0):
                rec[f"eps_est_LH{LH:g}"] = E[(round(a, 6), round(b, 12), LH)]
                rec[f"zeta0_LH{LH:g}"] = rec[f"eps_est_LH{LH:g}"] / max(c["n_b"], 1) / ((2 / 3) * math.pi ** 2 * (a / LH) ** 2)
            rows.append(rec)
        out[cls] = rows
    return out, len(wP), et


def pick(rows, n):
    rows = sorted(rows, key=lambda r: r["zeta0_LH1"])
    if len(rows) <= n:
        return rows
    return [rows[int(round(i))] for i in np.linspace(0, len(rows) - 1, n)]


def selection(ch):
    Ev = [r for r in ch["E"] if r["zeta0_LH1"] >= ZETA0]; En = [r for r in ch["E"] if r["zeta0_LH1"] < ZETA0]
    Dv = [r for r in ch["D"] if r["zeta0_LH1"] >= ZETA0]; Dn = [r for r in ch["D"] if r["zeta0_LH1"] < ZETA0]
    Dvb = [r for r in Dv if r["band"]]; Dvo = [r for r in Dv if not r["band"]]
    return {"E-v": pick(Ev, 10), "E-n": pick(En, 5), "D-v": pick(Dvb, 4) + pick(Dvo, 6), "D-n": pick(Dn, 5)}


def tag(r):
    return f"sep{r['sep']:g}_lam{r['lam']:.6g}_sig{r['sigma']:.6g}"


def state_of(r):
    from frw_c import c45_state
    B = benchmark()
    return c45_state(B["A_C45"], r["lam"], r["sigma"], 1, r["sep"])


def one_run(r, family, LH, eps, N=N_PROD, classify=True):
    """one perturbed parent: bounce-time spread over x <= X_HI and, if that passes, the separate-universe classification"""
    st = state_of(r); t0 = time.time()
    try:
        M, S, info = C.parent_data(st, N, family, eps, LH / r["AH"])
    except Exception as e:                                         # no admissible slice (e.g. 1 - K x^2 < 0)
        return dict(valid=False, reason=f"data: {e}")
    if info["rho_d_min"] < 0 or info["eL_min"] <= 0:
        return dict(valid=False, reason="negative dust or e^Lam on the slice", rho_d_min=info["rho_d_min"])
    out = C.evolve(M, S, mode="bounce", t_max=4 * (info["t_half"] * 2), track_rho=classify)
    sel = M.x <= X_HI
    res = dict(valid=out["status"] == "COMPLETE" and out["cmax"] <= 1e-2, status=out["status"], H_max=out["cmax"],
               steps=out["steps"], runtime_s=None, sqrt_al=math.sqrt(st["al"]), rho_d_min=info["rho_d_min"])
    if not res["valid"]:
        res["runtime_s"] = time.time() - t0; return res
    tB = out["tB"][sel]
    res["dtB"] = float(np.max(tB) - np.min(tB)); res["dtB_over_sqrtal"] = res["dtB"] / res["sqrt_al"]
    res["pass_dt"] = bool(res["dtB"] < 0.3 * res["sqrt_al"])
    if classify and res["pass_dt"]:
        res.update(classify_cells(M, out["S"], st, sel, r["cls"], out["rho_max"]))
        res["coherent"] = bool(res["f_loc"] >= 0.8)
    else:
        res["coherent"] = False
    res["runtime_s"] = time.time() - t0
    return res


def classify_cells(M, S, st, sel, cls, rho_max):
    """separate universes from the end state (every cell has bounced): f_loc = volume fraction of the region landing in the
    homogeneous child's class, fertile and sub-Planck (frw_core definitions; E_max of the continuation, see amendment 6)"""
    from ucommon import benchmark
    B = benchmark(); c2 = st["V0"] / B["V0_BENCH_C45"]
    pole_et = U27_CHI_ET
    r, U, G, phi, P, Y = S
    rp, eL, php, q, Pi, rho_d, rho, pr = M.derived(S)
    idx = np.nonzero(sel)[0]
    vol = (eL * r * r * M.h)[idx]
    pick = idx[np.unique(np.linspace(0, len(idx) - 1, N_SU).round().astype(int))]
    good = 0.0; tot = 0.0; labels = []
    # each sampled cell represents the cells up to the midpoint with its neighbours
    edges = np.concatenate(([idx[0]], ((pick[1:] + pick[:-1]) // 2), [idx[-1] + 1]))
    for k, i in enumerate(pick):
        w = float(np.sum((eL * r * r * M.h)[edges[k]:edges[k + 1]]))
        H = U[i] / r[i]; psi = rho[i] / 3 / math.sqrt(1 + (st["al"] * rho[i] / 3) ** 2)
        inv_a2 = psi - H * H
        if inv_a2 <= 0:
            labels.append("open"); tot += w; continue
        o = C.hom_continue(1 / math.sqrt(inv_a2), H, phi[i], Pi[i], rho_d[i], st)
        fert = o["outcome"] == "launched60" and o["chi_abs_max"] > 5
        sub = max(o["E_max"], float(rho_max[i])) / c2 < 1
        lab = ("E" if abs(o["chi_end"]) > pole_et else "D") if (fert and sub) else "none"
        labels.append(lab); tot += w; good += w if lab == cls else 0.0
    return dict(f_loc=good / tot, su_labels=labels)


U27_CHI_ET = None


V2_CODE = r"""
import sys, json, math, numpy as np
sys.path.insert(0, sys.argv[1])
from common.potential import H_INF
from ltb.growing_qt import build_qt
from ltb.growing import turnaround_radius
from ltb.shells import Shell
import u28_comoving as C
J = json.load(open(sys.argv[2]))["children"]
picks = sorted({(round(c["AH"], 6), c["alphaH2"]) for c in J if c["zeta0_max_LH1"] >= 4.6e-5})
picks = [picks[0], picks[len(picks) // 2], picks[-1]]
out = []
for AH, aH2 in picks:
    al = aH2 / H_INF ** 2; EPS = 1e-4
    x = np.linspace(0, 0.95, 20001); K2 = x ** 2 * (1 + EPS * np.sinc(x / (1 / AH)))
    sh, tw, xs, xm, mu, AHx = build_qt(x, K2, AH, al, n_sh=600)
    Kf = lambda xx: 1 + EPS * np.sinc(xx / (1 / AH))
    M, S, info = C.time_symmetric_data(256, al, 8 * math.pi, 0.0, 2 * mu, Kf, 0.0, lambda rho: (0 * rho, rho), "small", with_field=False)
    o = C.evolve(M, S, mode="bounce", t_max=1e12)
    sel = M.x <= 0.95
    ref = np.array([2 * Shell(mu * xx ** 3, turnaround_radius(mu * xx ** 3, Kf(np.array([xx]))[0] * xx * xx, al), al).tau_B for xx in M.x[sel]])
    out.append(dict(AH=AH, alphaH2=aH2, status=o["status"], max_rel=float(np.max(np.abs(o["tB"][sel] - ref) / ref)), H=o["cmax"]))
print(json.dumps(out))
"""


def gate_V1(picks):
    from frw_c import lib
    res = []
    for r in picks:
        st = state_of(r); fw = C.hom_to_bounce(st, +1)
        o = np.zeros(12); lib().set_deform(4, C.C_POLE, 0.0)
        try:
            lib().run(st["a0"], st["al"], st["V0"], st["chi0"], st["chid0"], 60 * 4.0 * st["a0"], 60.0, 1e-11, 400_000_000, 0.57, o, 0.0)
        finally:
            lib().set_deform(0, 1e30, 0.0)
        M, S, info = C.parent_data(st, N_PROD, None, 0.0, 1.0)
        out = C.evolve(M, S, mode="bounce", t_max=8 * info["t_half"])
        T = info["t_half"] + fw["t"]; a_end = float(np.mean(out["S"][0] / M.x))
        res.append(dict(child=tag(r), status=out["status"], tB_rel=float(abs(np.nanmean(out["tB"]) - T) / T),
                        a_rel_frw=abs(a_end - o[6]) / o[6], a_rel_ode=abs(a_end - fw["a"]) / fw["a"], H=out["cmax"]))
    return res, all(x["status"] == "COMPLETE" and x["tB_rel"] <= 1e-3 and x["a_rel_frw"] <= 1e-3 for x in res)


def gate_V4():
    res = []
    for N in (64, 128, 256):
        M, S, info = C.time_symmetric_data(N, 0.0, 8 * math.pi, 0.0, 1.0, lambda x: np.ones_like(x), 0.0,
                                           lambda rho: (0 * rho, rho), "large", with_field=False)
        out = C.evolve(M, S, mode="collapse", r_stop=1e-3 * S[0][0], t_max=10)
        a = float(np.mean(out["S"][0] / np.sin(M.chi))); eta = math.acos(2 * a - 1); te = (eta + math.sin(eta)) / 2
        res.append(dict(N=N, rel=abs(out["t"] - te) / te, H=out["cmax"]))
    return res, all(x["rel"] <= 1e-3 for x in res)


def gates():
    global U27_CHI_ET
    import subprocess
    ch, nfert, et = children(); U27_CHI_ET = et
    sel = selection(ch)
    out = dict(n_fertile=nfert, chi_et=et, selection={k: [tag(r) for r in v] for k, v in sel.items()})
    out["V4"], out["V4_pass"] = gate_V4(); print("V4", out["V4_pass"], flush=True)
    r = subprocess.run([sys.executable, "-B", "-c", V2_CODE, str(ROOT / "components" / "c45j_ltb"), str(ROOT / "data" / "c45b" / "c45j_results.json")],
                       cwd=ROOT / "components" / "c45j_ltb", capture_output=True, text=True, env={"PYTHONPATH": str(ROOT / "chain"), "PATH": "/usr/bin:/bin"})
    out["V2"] = json.loads(r.stdout.strip().splitlines()[-1]) if r.returncode == 0 else r.stderr[-2000:]
    out["V2_pass"] = r.returncode == 0 and all(x["status"] == "COMPLETE" and x["max_rel"] <= 1e-3 for x in out["V2"]); print("V2", out["V2_pass"], flush=True)
    picks = [sel["E-v"][0], sel["E-v"][-1], sel["D-v"][0], sel["D-v"][-1], sel["D-n"][2]]
    out["V1"], out["V1_pass"] = gate_V1(picks); print("V1", out["V1_pass"], flush=True)
    v3 = []
    for r in (sel["E-v"][0], sel["E-v"][-1], sel["D-v"][0], sel["D-v"][-1]):
        a = one_run(r, "a", 1.0, 1e-3, N=N_PROD); b = one_run(r, "a", 1.0, 1e-3, N=2 * N_PROD)
        rel = abs(a.get("dtB", np.nan) - b.get("dtB", np.nan)) / b.get("dtB", np.nan)
        v3.append(dict(child=tag(r), N=a, N2=b, rel_change=rel, same_label=a.get("coherent") == b.get("coherent"),
                       H_ratio=(a.get("H_max", np.nan) / b.get("H_max", np.nan))))
        print("V3", tag(r), rel, v3[-1]["same_label"], a.get("runtime_s"), b.get("runtime_s"), flush=True)
    out["V3"] = v3
    out["V3_pass"] = bool(all(v["rel_change"] < 0.1 and v["same_label"] for v in v3))
    out["V0_H_ratios"] = [v["H_ratio"] for v in v3]
    out["V0_pass"] = bool(all(v["H_ratio"] >= 3 for v in v3))          # order >= ~1.6 at least; ~4 expected for order 2
    RUNS.mkdir(parents=True, exist_ok=True)
    (RUNS / "gates.json").write_text(json.dumps(out, indent=1, default=float))
    print({k: out[k] for k in out if k.endswith("_pass")})


def groups_LH(group):
    return (1.0, 0.3) if group.endswith("-v") else (1.0,)


def run_file(r, family, LH, eps):
    return RUNS / f"{tag(r)}_{family}_LH{LH:g}_eps{eps:g}.json"


def get_run(r, family, LH, eps):
    f = run_file(r, family, LH, eps)
    if f.exists():
        return json.loads(f.read_text())
    o = one_run(r, family, LH, eps); o.update(child=tag(r), family=family, LH=LH, eps=eps, cls=r["cls"])
    f.write_text(json.dumps(o, default=float)); return o


def produce():
    global U27_CHI_ET
    ch, nfert, et = children(); U27_CHI_ET = et
    sel = selection(ch); RUNS.mkdir(parents=True, exist_ok=True)
    (RUNS / "selection.json").write_text(json.dumps({k: [dict(tag=tag(r), **{q: r[q] for q in ("sep", "lam", "sigma", "AH", "alphaH2", "n_b", "cls", "band", "chi_land", "chi_abs_max", "eps_est_LH1", "eps_est_LH0.3", "zeta0_LH1", "zeta0_LH0.3")}) for r in v] for k, v in sel.items()}, indent=1, default=float))
    for group in ("E-v", "D-v", "E-n", "D-n"):
        for r in sel[group]:
            for family in ("a", "b"):
                for LH in groups_LH(group):
                    runs = [get_run(r, family, LH, e) for e in GRID]
                    valid = [o for o in runs if o.get("valid")]
                    if valid and all(o["coherent"] for o in valid):
                        runs += [get_run(r, family, LH, e) for e in UP]
                    elif not any(o.get("coherent") for o in valid):
                        runs += [get_run(r, family, LH, e) for e in DOWN]
                    print(group, tag(r), family, LH, [(o["eps"], o.get("valid"), o.get("coherent")) for o in runs], flush=True)


def thresholds(runs, ok=lambda o: o["coherent"]):
    """amendment 6, point 3 (and amendment 1 censoring): (eps_lo, eps_hi) from the valid runs of one child/family/LH;
    `ok` is the coherence test (the preregistered one by default; the informative readings pass others)"""
    v = sorted([o for o in runs if o.get("valid")], key=lambda o: o["eps"])
    if not v:
        return None
    amps = sorted(set(GRID) | {o["eps"] for o in runs})
    coh = [o["eps"] for o in v if ok(o)]
    if not coh:
        return dict(lo=0.0, hi=min(o["eps"] for o in v), censored="below", monotone=True)
    lo = max(coh)
    inc_above = [o["eps"] for o in v if not ok(o) and o["eps"] > lo]
    if inc_above:
        hi = min(inc_above)
    else:
        above = [a for a in amps if a > lo]
        hi = above[0] if (above and any(not o.get("valid") for o in runs if o["eps"] > lo)) else math.inf
    mono = all(ok(o) for o in v if o["eps"] <= lo)
    return dict(lo=lo, hi=hi, censored="above" if hi == math.inf else None, monotone=mono)


def analyze():
    sel = json.loads((RUNS / "selection.json").read_text())
    gates_ = json.loads((RUNS / "gates.json").read_text())
    runs = [json.loads(f.read_text()) for f in RUNS.glob("sep*.json")]
    by = {}
    for o in runs:
        by.setdefault((o["child"], o["LH"]), {}).setdefault(o["family"], []).append(o)
    per = {}
    for group, rows in sel.items():
        for r in rows:
            for LH in groups_LH(group):
                fam = {f: thresholds(by.get((r["tag"], LH), {}).get(f, [])) for f in ("a", "b")}
                ok = {f: t for f, t in fam.items() if t is not None}
                if not ok:
                    per[(r["tag"], LH)] = dict(group=group, LH=LH, valid=False); continue
                f_min = min(ok, key=lambda f: ok[f]["lo"])
                t = ok[f_min]; e = r[f"eps_est_LH{LH:g}"]
                nruns = sum(len(v) for v in by.get((r["tag"], LH), {}).values())
                ninv = sum(1 for v in by.get((r["tag"], LH), {}).values() for o in v if not o.get("valid"))
                per[(r["tag"], LH)] = dict(group=group, cls=r["cls"], LH=LH, family=f_min, eps_lo=t["lo"], eps_hi=t["hi"], censored=t["censored"],
                                           rho_lo=t["lo"] / (0.3 * e), rho_hi=t["hi"] / (0.3 * e), thresholds=fam, n_runs=nruns, n_invalid=ninv,
                                           zeta0_dust=r[f"zeta0_LH{LH:g}"], band=r["band"])
    def cls_stats(group, LH=1.0):
        P = [v for (tg, L), v in per.items() if v["group"] == group and L == LH and v.get("valid", True) is not False and "rho_lo" in v]
        n_runs = sum(v["n_runs"] for v in P); n_inv = sum(v["n_invalid"] for v in P)
        if not P:
            return dict(n=0)
        mlo, mhi = float(np.median([v["rho_lo"] for v in P])), float(np.median([v["rho_hi"] for v in P]))
        verdict = "MEJOR_QUE_POLVO" if mlo >= 3 else ("PEOR" if mhi <= 1 / 3 else "IGUAL")
        cens = [v["censored"] for v in P]
        note = None
        for side in ("above", "below"):
            if sum(c == side for c in cens) > len(P) / 2:
                note = f"more than half censored {side}"
        inconcl = n_runs > 0 and n_inv / n_runs > 0.2
        return dict(n=len(P), median_rho_lo=mlo, median_rho_hi=mhi, verdict=("INCONCLUSA" if inconcl else verdict), censoring_note=note,
                    invalid_fraction=n_inv / n_runs if n_runs else None)
    stats = {g: cls_stats(g) for g in ("E-v", "D-v")}
    stats_LH03 = {g: cls_stats(g, 0.3) for g in ("E-v", "D-v")}
    # corrected viable fractions over all children of each class (U27 sets), L.H = 1
    ch, nfert, et = children()
    out_f = {}
    for cls, g in (("E", "E-v"), ("D", "D-v")):
        s_ = stats[g]
        if not s_.get("n"):
            continue
        own = {tg: v for (tg, L), v in per.items() if L == 1.0 and v.get("cls") == cls and "rho_lo" in v}
        nlo = nhi = 0; nlo_b = nhi_b = 0
        for r in ch[cls]:
            t_ = tag(r); z = r["zeta0_LH1"]
            rl, rh = (own[t_]["rho_lo"], own[t_]["rho_hi"]) if t_ in own else (s_["median_rho_lo"], s_["median_rho_hi"])
            vlo, vhi = z * rl >= ZETA0, z * rh >= ZETA0
            nlo += vlo; nhi += vhi
            if 10 <= r["chi_abs_max"] <= 50:
                nlo_b += vlo; nhi_b += vhi
        flo, fhi = nlo / nfert, nhi / nfert
        lab = lambda f: "VIABLE" if f >= 0.05 else ("MARGINAL" if f >= 0.01 else "NO_COHERENTE")
        out_f[cls] = dict(n_viable_lo=nlo, n_viable_hi=nhi, f_lo=flo, f_hi=fhi, verdict_D11=lab(flo), verdict_upper=lab(fhi),
                          band_viable_lo=nlo_b, band_f_lo=nlo_b / nfert)
    controls = [dict(child=tg, **{k: v[k] for k in ("group", "eps_lo", "eps_hi", "rho_lo")}) for (tg, L), v in per.items()
                if v["group"] in ("E-n", "D-n") and "rho_lo" in v and v["rho_lo"] >= 3 and v["zeta0_dust"] * v["rho_lo"] >= ZETA0]
    info = informative(sel, by)
    res = dict(spec="docs/preregistros/ESPEC_U28.md", gates={k: gates_[k] for k in gates_ if k.endswith("_pass")}, class_LH1=stats, class_LH03=stats_LH03,
               viable=out_f, controls_rescued=controls, per_child={f"{tg}|LH{L:g}": v for (tg, L), v in per.items()},
               n_runs=len(runs), n_fertile=nfert, informative_post_hoc=info)
    save("u28_c45k.json", res)
    print(json.dumps({k: res[k] for k in ("gates", "class_LH1", "class_LH03", "viable", "controls_rescued", "n_runs")}, indent=1, default=float))
    print(json.dumps({k: info[k] for k in ("class_dt_only", "class_fertile_any", "f_k_quartiles")}, indent=1, default=float))


def fertile_any(o):
    """informative: passes the bounce-time test and >= 80 % of the sampled cells land fertile and sub-Planck in either class
    (unweighted count of su_labels; the cell weights are not stored)"""
    if not o.get("pass_dt"):
        return False
    lab = o.get("su_labels") or []
    return bool(lab) and sum(l in ("E", "D") for l in lab) / len(lab) >= 0.8


def kinetic_fraction(r):
    """field kinetic fraction of the homogeneous child at its bounce: (chid^2/2) / (rho_dust + chid^2/2 + V)"""
    st = state_of(r); b = C.hom_to_bounce(st, +1)
    kin = 0.5 * b["chid"] ** 2
    return kin / (3 * b["Aeff"] / b["a"] ** 3 + kin + float(C.V_pole(b["chi"], st["V0"])))


def informative(sel, by):
    """readings decided after seeing the first runs (not preregistered, see INFORME_U28 sec. 5): the same thresholds with
    (i) the bounce-time test alone, comparable with C45j, and (ii) landing fertile in either class; and f_k per child"""
    out = dict(note="post hoc, decided after seeing the first runs; the verdict is the preregistered one")
    for name, ok in (("dt_only", lambda o: bool(o.get("pass_dt"))), ("fertile_any", fertile_any)):
        cl = {}
        for group, rows in sel.items():
            for LH in groups_LH(group):
                vals = []
                for r in rows:
                    fam = {f: thresholds(by.get((r["tag"], LH), {}).get(f, []), ok) for f in ("a", "b")}
                    fam = {f: t for f, t in fam.items() if t is not None}
                    if not fam:
                        continue
                    t = min(fam.values(), key=lambda t: t["lo"]); e = r[f"eps_est_LH{LH:g}"]
                    vals.append(dict(child=r["tag"], rho_lo=t["lo"] / (0.3 * e), rho_hi=t["hi"] / (0.3 * e), censored=t["censored"]))
                if vals:
                    cl[f"{group}|LH{LH:g}"] = dict(n=len(vals), median_rho_lo=float(np.median([v["rho_lo"] for v in vals])),
                                                  median_rho_hi=float(np.median([v["rho_hi"] for v in vals])), per_child=vals)
        out[f"class_{name}"] = cl
    fk = {r["tag"]: dict(group=g, f_k=kinetic_fraction(r)) for g, rows in sel.items() for r in rows}
    out["f_k"] = fk
    out["f_k_quartiles"] = {g: [float(q) for q in np.percentile([v["f_k"] for v in fk.values() if v["group"] == g], (0, 25, 50, 75, 100))]
                            for g in sel}
    return out


if __name__ == "__main__":
    {"gates": gates, "produce": produce, "analyze": analyze}[sys.argv[1]]()
