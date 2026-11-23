"""U68b analysis (docs/preregistros/ESPEC_U68.md and ESPEC_U68b.md): validity (U19b's four conditions), the mean velocity of the
soft quanta v = sum_c sum_{k < 10 m_c a_1} (k / m_c a_1) n_c / sum_c sum_{k < 10 m_c a_1} n_c, the paired D = v_GEN / v_O3 - 1 and the
verdict IGUAL / DISTINTO / INTERMEDIO / INCONCLUSO; informative f_L, D_rho and per-mode readings, against the (IND) expectations.
    python -B network/u68_analyze.py   ->  network/results/u68/U68_RESULT.json"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; OUT = HERE / "results" / "u68"
SEEDS = (91, 92, 93)
IND = dict(f_light=dict(O3=1.0, GEN=1 / (1 + 2 ** -0.25 + 4 ** -0.25), O2=2 / (2 + 4 ** -0.25)),
           rho_over_O3=dict(O3=1.0, GEN=(1 + 2 ** 0.25 + 2 ** 0.5) / 3, O2=(2 + 2 ** 0.5) / 3))


def per_run(r):
    meas = [x for x in r["rows"] if not x.get("initial")]
    em = np.array([x["ell"] for x in meas]); xi = np.array([x["xi"] for x in meas]); wall = np.array([x["wall_frac"] for x in meas])
    a = np.array([x["a"] for x in meas]); l1, a1 = r["ell_1"], r["a_1"]; ms = r["masses"]; mG = r["m_G"]
    xi1 = float(np.interp(l1, em, xi)); last = meas[-1]
    lin = [x for x in meas if "N_lin" in x]; NN = np.array([x["N_lin"] for x in lin]); Nend = last["N_lin"]
    Nhalf = float(np.interp(em[-1] - 0.5, [x["ell"] for x in lin], NN))
    checks = dict(strings=bool(xi[-1] < 0.01 * xi1), walls=bool(wall[-1] < 0.05 * wall.max()), counter=bool(last["far_frac"] < 0.01),
                  plateau=bool(abs(Nend / Nhalf - 1) < 0.03))
    e = np.array(r["k_edges"]); kc = 0.5 * (e[1:] + e[:-1])
    num = den = 0.0; modes = []
    for c, sp in enumerate(last["lin_spectrum_modes"]):
        n = np.array(sp["n"]); kn = np.array(sp["kn"]); soft = kc < 10 * ms[c] * a1
        Ns, Ps = float(n[soft].sum()), float(kn[soft].sum())
        num += Ps / (ms[c] * a1); den += Ns
        modes.append(dict(m_over_mG=ms[c] / mG, N_total=last["N_lin_modes"][c], N_soft=Ns, v_soft=Ps / (ms[c] * a1) / Ns if Ns > 0 else None,
                          f_hard_number=1 - Ns / max(float(n.sum()), 1e-300)))
    v = num / den
    rho = float(sum(ms[c] * last["N_lin_modes"][c] for c in range(3)) / mG)
    light = [c for c in range(3) if abs(ms[c] / mG - 1) < 1e-9]
    f_light = sum(modes[c]["N_soft"] for c in light) / den
    imax = int(np.argmax(wall)); after = np.flatnonzero((np.arange(len(wall)) > imax) & (wall < 0.1 * wall[imax]))
    r_W = float(a[after[0]] / a1) if len(after) else float("nan")
    return dict(tag=f"{r['protocol']}_s{r['seed']}", protocol=r["protocol"], seed=r["seed"], xi_1=xi1, checks=checks, valid=all(checks.values()),
                v_soft=v, rho_over_mG=rho, f_light=f_light, f_c1=modes[0]["N_soft"] / den, modes=modes, r_W=r_W, N_end=Nend,
                plateau_rel=Nend / Nhalf - 1, xi_end=float(xi[-1]), wall_end_over_max=float(wall[-1] / max(wall.max(), 1e-30)),
                far_frac_end=last["far_frac"], runtime_s=r["runtime_s"])


def main():
    chk = json.loads((OUT / "U68_CHECKS.json").read_text())
    runs = {}
    for f in sorted((OUT / "runs").glob("RP3_N192_*_s*.json")):
        r = json.loads(f.read_text()); runs[(r["protocol"], r["seed"])] = per_run(r)
    complete = all((p, s) in runs for p in ("O3", "GEN", "O2") for s in SEEDS)
    G3 = complete and all(x["valid"] for x in runs.values())
    gates = dict(G0=chk["PASS"]["G0"], G1=chk["PASS"]["G1"], G2=chk["PASS"]["G2"], G3=bool(G3))

    def paired(p, key):
        v = [runs[(p, s)][key] / runs[("O3", s)][key] - 1 for s in SEEDS if (p, s) in runs and ("O3", s) in runs]
        return (float(np.mean(v)), float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float("inf"), v)
    D, sig, Ds = paired("GEN", "v_soft")
    if not all(gates.values()): verdict = "INCONCLUSO"
    elif abs(D) <= 0.15: verdict = "IGUAL"
    elif abs(D) > 0.30 and abs(D) > 3 * sig: verdict = "DISTINTO"
    else: verdict = "INTERMEDIO"
    info = {}
    for p in ("GEN", "O2", "O3"):
        dv = paired(p, "v_soft"); dr = paired(p, "rho_over_mG")
        fl = [runs[(p, s)]["f_light"] for s in SEEDS if (p, s) in runs]
        info[p] = dict(D_v=dv[:2], D_v_seeds=dv[2], rho_over_O3=(1 + dr[0], dr[1]), rho_IND=IND["rho_over_O3"][p],
                       f_light=(float(np.mean(fl)), float(np.std(fl, ddof=1) / math.sqrt(len(fl))) if len(fl) > 1 else None), f_light_IND=IND["f_light"][p],
                       v_soft_mean=float(np.mean([runs[(p, s)]["v_soft"] for s in SEEDS if (p, s) in runs])),
                       r_W=float(np.nanmean([runs[(p, s)]["r_W"] for s in SEEDS if (p, s) in runs])),
                       per_mode_v=[float(np.mean([runs[(p, s)]["modes"][c]["v_soft"] for s in SEEDS if (p, s) in runs])) for c in range(3)],
                       per_mode_N_over_O3=[float(np.mean([runs[(p, s)]["modes"][c]["N_total"] / (runs[("O3", s)]["N_end"] / 3)
                                                          for s in SEEDS if (p, s) in runs])) for c in range(3)])
    res = dict(spec=["docs/preregistros/ESPEC_U68.md", "docs/preregistros/ESPEC_U68b.md"], gates=gates, D=D, sigma=sig, D_seeds=Ds,
               verdict=verdict, informative=info, runs={f"{p}_s{s}": v for (p, s), v in runs.items()})
    (OUT / "U68_RESULT.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(dict(gates=gates, D=D, sigma=sig, D_seeds=Ds, verdict=verdict, informative=info), indent=1, default=float))


if __name__ == "__main__":
    main()
