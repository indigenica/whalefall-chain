"""C45b-dn bank: collapse of chi configurations in spherical NPQT gravity with a
double-null code that crosses the horizon.  Question: after the horizon forms,
does the core bounce and LAUNCH chi onto the plateau (V-dominated expansion)?

python -m dn_bank.run_c45b_dn [--quick]
"""
import sys, time, json, warnings
import numpy as np
from common.provenance import new_run_dir, write_json, write_manifest, source_hashes, environment
from common.potential import V, dV, V0, H_INF
from dn.solver import Grid, initial_ray, evolve
from dn.analysis import centre_history, summary
from dn.convergence import trusted_mask
from dn.profiles import gaussian, chistar

warnings.filterwarnings("ignore")
CHI_END = 0.6555 / np.sqrt(8 * np.pi)          # end of slow roll, G = 1 units (0.1308)

PROTOCOL = {
    "status": "written AFTER exploratory probes (declared): gaussian A=0.25 w=3 alpha=0.1 (N=800/1600/3200) "
              "and chi-star s=1.5, 2.0, 2.5, 3.0 alpha=0.1 (N=1600; s=2.0 also N=800/3200). Probe outcome seen: "
              "horizon, QT-saturated centre that bounces and is recompressed, kinetic-dominated. "
              "The criteria below were fixed before running this bank; a --quick smoke run (3 configs, alpha = 0.1) "
              "was used only to fix the centre-density estimator and the trust tolerance (declared).",
    "units": "G = 1, m_chi = 1; potential common/potential.py (V0 = 8.695e-3, H_inf = 0.270)",
    "code": "dn/solver.py: diamond scheme, regular centre, psi from generalized mass, causal masking "
            "(failed nodes -> NaN, only their causal future is lost)",
    "families": {
        "star": "phi(u=0, r) = sqrt2 * s * phi_star(r), C59a ground state (phi0 = 0.12), tapered tail; "
                "s in {1.75, 2.0, 2.5, 3.0}",
        "gauss": "A exp(-(r/w)^2), (A, w) in {(0.25, 3), (0.15, 5), (0.1, 8)}",
        "tophat": "A / (1 + exp((r - R)/0.5)), A = 0.3, R in {3, 4, 6} (closest analogue of the homogeneous patch)",
    },
    "alpha": [0.1, 0.03],
    "grid": "L = 80 (initial cone r <= 40), N in {800, 1600, 3200}; --fine (WSL job): {1600, 3200, 6400}, "
            "trust and events then use the two finest",
    "trust": "centre trusted while |phi_c(N=1600) - phi_c(3200)| < 0.02 and |delta rho_c| < 0.05 max(rho_c) "
             "(rho_c intrinsic: phi_v = (1/2) d phi_c/du on the centre line) "
             "(same u); events are evaluated on the N = 3200 centre history up to that time only. "
             "2D trusted mask (dn/convergence.py) reported for information",
    "events": {
        "H": "trapped node (r_u < 0, r_v <= 0) anywhere",
        "B": "centre QT bounce: comoving H_c changes sign - to + while alpha psi_c >= 0.9",
        "L": "launch: an interval of centre proper time >= 1/H_inf during which kinetic fraction < 0.5, "
             "rho_c < 10 V0, H_c > 0 and |phi_c| > CHI_END hold continuously (with or without B; "
             "L_after_B recorded separately)",
        "R": "recompression after B: rho_c exceeds 10x its post-bounce minimum",
        "MI": "inner-horizon diagnostic: max interior (r < 2) generalized mass / M_ADM, over the 2D trusted "
              "mask (N = 1600 vs 3200) and, for information, over all finite nodes",
    },
    "verdict": {
        "PASS": "L in at least one run with H, for each alpha, in the star family",
        "FAIL": "no L in any run with H, and every run with H and B shows R inside the trusted domain",
        "INCONCLUSIVE": "otherwise",
        "note": "the specification's PASS (lambda in the C54 window, e-folds within 10 % of C45) can only be "
                "evaluated if L occurs; then the trusted centre state is continued with the homogeneous model",
    },
}


def fermi(A, R, d=0.5):
    F = lambda r: A / (1 + np.exp((np.asarray(r) - R) / d))
    dF = lambda r: -F(r) * (1 - F(r) / A) / d
    return F, dF


def configs(quick):
    out = []
    for al in ([0.1] if quick else [0.1, 0.03]):
        for s in ([2.0] if quick else [1.75, 2.0, 2.5, 3.0]):
            out.append(dict(alpha=al, family="star", s=s))
        for A, w in ([(0.25, 3)] if quick else [(0.25, 3), (0.15, 5), (0.1, 8)]):
            out.append(dict(alpha=al, family="gauss", A=A, w=w))
        for R in ([4] if quick else [3, 4, 6]):
            out.append(dict(alpha=al, family="tophat", A=0.3, R=R))
    return out


def profile(c):
    if c["family"] == "star":
        F, dF, _ = chistar(c["s"], c["alpha"]); return F, dF
    if c["family"] == "gauss":
        return gaussian(c["A"], c["w"])
    return fermi(c["A"], c["R"])


RES = (800, 1600, 3200)


def one(c, N, L=80.0):
    F, dF = profile(c)
    g = Grid(N, L / N)
    initial_ray(g, F, dF, V, dV, c["alpha"], 1.0)
    Madm = float(g.M[0, -1])
    t0 = time.time()
    st = evolve(g, V, dV, c["alpha"], 1.0, stop_on_fail=False)
    ch = centre_history(g, V, c["alpha"])
    sm = summary(g, V, c["alpha"], M_adm=Madm)
    sm.update(status=st, M_adm=Madm, wall=time.time() - t0)
    return g, ch, sm


def centre_trust(chA, chB):
    """chA at N, chB at 2N: last u where they agree (interpolating B on A's u)."""
    uA, uB = chA["u"], chB["u"]
    n = len(uA)
    last = 0.0
    for q in range(n):
        if uA[q] > uB[-1]:
            break
        pb = np.interp(uA[q], uB, chB["phi"]); rb = np.interp(uA[q], uB, chB["rho"])
        if not (abs(chA["phi"][q] - pb) < 0.02 and abs(chA["rho"][q] - rb) < 0.05 * max(chA["rho"][q], rb)):
            break
        last = uA[q]
    return float(last)


def events(ch, u_trust, alpha):
    m = ch["u"] <= u_trust
    u, rho, H, ap, kin, phi, tau = (ch[k][m] for k in ("u", "rho", "H", "apsi", "kin_frac", "phi", "tau"))
    ev = dict(B=False, L=False, R=False)
    if len(u) < 3:
        return ev
    Hs = np.where(np.isfinite(H), H, 0.0)
    for q in range(1, len(u)):
        if Hs[q - 1] < 0 <= Hs[q] and ap[q] >= 0.9:
            ev.update(B=True, u_B=float(u[q]), tau_B=float(tau[q]), rho_B=float(rho[q]), phi_B=float(phi[q]))
            break
    good = (kin < 0.5) & (rho < 10 * V0) & (Hs > 0) & (np.abs(phi) > CHI_END)
    run_start = None
    for q in range(len(u)):
        if good[q]:
            if run_start is None:
                run_start = q
            if tau[q] - tau[run_start] >= 1 / H_INF:
                ev.update(L=True, u_L=float(u[run_start]))
                ev["L_after_B"] = bool(ev["B"] and u[run_start] >= ev["u_B"])
                break
        else:
            run_start = None
    if ev["B"]:
        post = u >= ev["u_B"]
        rp, kp, Hp, pp, up = rho[post], kin[post], Hs[post], phi[post], u[post]
        runmin = np.minimum.accumulate(rp)
        Rm = rp > 10 * runmin
        if Rm.any():
            k = int(np.argmax(Rm))
            ev.update(R=True, u_R=float(up[k]), rho_min_before_R=float(runmin[k]), rho_R=float(rp[k]))
        ev["phi_abs_max_post"] = float(np.max(np.abs(pp)))
        ev["kin_min_post"] = float(np.min(kp))
        ev["rho_min_post"] = float(np.min(rp))
    ev["apsi_max"] = float(np.max(ap))
    return ev


def main(quick=False, fine=False):
    global RES
    if fine:
        RES = (1600, 3200, 6400)
    d = new_run_dir("c45b_dn" + ("_quick" if quick else "") + ("_fine" if fine else ""))
    write_json(d / "PROTOCOL.json", dict(PROTOCOL, sources=source_hashes(), env=environment()))
    rows = []
    for c in configs(quick):
        res = {}
        for N in RES:
            g, ch, sm = one(c, N)
            res[N] = (g, ch, sm)
            np.savez_compressed(d / f"centre_{c['family']}_{c['alpha']}_{c.get('s', c.get('A'))}_{c.get('w', c.get('R', ''))}_N{N}.npz",
                                **ch)
        ut = centre_trust(res[RES[1]][1], res[RES[2]][1])
        ut_lo = centre_trust(res[RES[0]][1], res[RES[1]][1])
        tm, _ = trusted_mask(dict(R=res[RES[1]][0].r, S=res[RES[1]][0].s, P=res[RES[1]][0].p, M=res[RES[1]][0].M),
                             dict(R=res[RES[2]][0].r, S=res[RES[2]][0].s, P=res[RES[2]][0].p, M=res[RES[2]][0].M))
        frac_trusted = float(tm.sum() / np.isfinite(res[RES[1]][0].r).sum())
        Mt = np.where(tm & (res[RES[1]][0].r < 2), res[RES[1]][0].M, np.nan)
        MI_trusted = float(np.nanmax(Mt) / res[RES[1]][2]["M_adm"]) if np.isfinite(Mt).any() else None
        ev = events(res[RES[2]][1], ut, c["alpha"])
        sm = res[RES[2]][2]
        row = dict(config=c, M_adm=sm["M_adm"], H=bool(sm["trapped"]),
                   first_trapped=[sm.get("first_trapped_u"), sm.get("first_trapped_v"), sm.get("first_trapped_r")],
                   first_trapped_N1600=[res[RES[1]][2].get("first_trapped_u"), res[RES[1]][2].get("first_trapped_r")],
                   u_trust_centre=ut, u_trust_centre_800_1600=ut_lo, frac_trusted_2d=frac_trusted,
                   MI_ratio={N: res[N][2].get("mass_inflation_ratio") for N in RES[1:]},
                   MI_ratio_trusted=MI_trusted,
                   sigma_min={N: res[N][2]["sigma_min"] for N in RES[1:]},
                   ok_until_uv={N: res[N][2]["ok_until_uv"] for N in RES},
                   events=ev, wall={N: res[N][2]["wall"] for N in RES})
        rows.append(row)
        print(json.dumps(dict(cfg=c, M=round(sm["M_adm"], 3), H=row["H"], ut=ut, ev={k: v for k, v in ev.items()
                              if k in ("B", "L", "R", "u_B", "u_R", "phi_abs_max_post", "kin_min_post", "apsi_max")},
                              MI_trusted=MI_trusted), default=float), flush=True)
        write_json(d / "results_partial.json", dict(rows=rows))
        del res
    verdict = "INCONCLUSIVE"
    star = [r for r in rows if r["config"]["family"] == "star" and r["H"]]
    alphas = sorted(set(r["config"]["alpha"] for r in rows))
    if all(any(r["events"]["L"] for r in star if r["config"]["alpha"] == a) for a in alphas):
        verdict = "PASS"
    else:
        withH = [r for r in rows if r["H"]]
        noL = not any(r["events"]["L"] for r in withH)
        allR = all(r["events"]["R"] for r in withH if r["events"]["B"])
        anyB = any(r["events"]["B"] for r in withH)
        if noL and allR and anyB:
            verdict = "FAIL"
    write_json(d / "results.json", dict(rows=rows, verdict=verdict, CHI_END=CHI_END, V0=V0, H_INF=H_INF))
    write_manifest(d)
    print("VERDICT", verdict)
    return d


if __name__ == "__main__":
    main(quick="--quick" in sys.argv, fine="--fine" in sys.argv)
