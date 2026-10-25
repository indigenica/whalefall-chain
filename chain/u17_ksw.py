"""U17 (docs/preregistros/ESPEC_U17.md): the string-wall collapse contribution K_sw made consistent with the energy left in the
network at H = m_G/3,  K_eff = n_collapse / n_mis = (8 pi / 27) xi(ell) ell / (c/eta),  with ell = ln(3 m_r / m_G).
m_G is solved self-consistently: Omega_mis(m_G) (1 + K_eff(m_G)) + Omega_em(m_G) = 0.12 (emission of U12 unchanged).
Collapse quanta: c is the mean ENERGY per quantum / m_G (ESPEC_U17), so with eta = 1 their momentum is sqrt(c^2 - 1) m_G at
H = m_G/3; free-streaming vs the 5.3 keV WDM.  (Before 14 Oct 2026 the free-streaming used momentum c m_G, i.e. energy and
momentum were mixed: codex audit, audit/codex_20261014; the number per energy was always c m_G.)
R_rho (default 1) multiplies the network energy: RP3/U(1) energy-density ratio at the collapse (ESPEC_U17d).
    python -B chain/u17_ksw.py      ->  results/u17_ksw.json"""
import itertools, json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from ucommon import ROOT, RESULTS, benchmark
import c35b_pseudo_goldstone_dm_rp3 as C35
from s4_dark import abundance, at_epoch, lambda_fs, lambda_wdm_reference, scale
from s6_measured import source_eps
sys.path.insert(0, str(ROOT / "lya"))

R_GAMMA = 1.56; BRANCHES = ("B2", "B4a", "B4b", "B1"); QS = (0.87, 1.0); X0S = (10, 30); C_ETA = (3.0, 10.0, 30.0)


def xi_branch(branch, ell):
    if branch == "B2": return 1.20
    x = max(0.616 + 0.133 * (ell - 5.65), 0.3)
    if branch == "B4a": return x
    if branch == "B4b": return x * 1.20 / 0.616
    if branch == "B1": return 0.24 * ell
    raise ValueError(branch)


def eps_fn(branch):
    if branch == "B2": return lambda ell: R_GAMMA * 1.20 / (0.24 * np.asarray(ell))
    if branch == "B1": return lambda ell: np.full_like(np.asarray(ell, dtype=float), R_GAMMA)
    k = 1.0 if branch == "B4a" else 1.20 / 0.616
    return lambda ell: R_GAMMA * k * np.maximum(0.616 + 0.133 * (np.asarray(ell) - 5.65), 0.3) / (0.24 * np.asarray(ell))


def k_eff(branch, F, m_eV, c_eta, R_rho=1.0):
    """c_eta = (mean energy per collapse quantum / m_G) / eta"""
    if c_eta is None: return None
    ell = math.log(3 * (F / 57.735) / (m_eV * 1e-9))
    return 8 * math.pi / 27 * R_rho * xi_branch(branch, ell) * ell / c_eta


def solve(F, Trh, branch, q, x0, c_eta=None, K_fixed=None, R_rho=1.0):
    spec = f"powerlaw_q{q:g}"; ef = eps_fn(branch)
    K = (lambda m: K_fixed) if c_eta is None else (lambda m: k_eff(branch, F, m, c_eta, R_rho))
    def budget(lm):
        m = 10 ** lm
        return C35.omega(m, F, K(m))[0] + abundance(source_eps(F, m, Trh, x0, ef, spectrum=spec))
    m = 10 ** brentq(lambda lm: math.log(budget(lm) / 0.12), -27, -6, xtol=1e-7)
    pk = source_eps(F, m, Trh, x0, ef, spectrum=spec); em = abundance(pk)
    Kv = K(m); _, Tosc = C35.omega(m, F, Kv)
    out = dict(branch=branch, q=q, x0=x0, c_over_eta=c_eta, K=Kv, m_G_eV=m, T_osc_MeV=Tosc * 1e3, f_em=em / 0.12,
               M_minicluster_Msun=C35.minicluster_Msun(m), DeltaNeff_CMB_emission=at_epoch(pk, 0.26e-9)["DeltaNeff_radiation_moment"],
               ell_osc=math.log(3 * (F / 57.735) / (m * 1e-9)))
    if c_eta is not None:
        a = float(scale(Tosc)); p = math.sqrt(max(c_eta ** 2 - 1, 0.0))          # eta = 1: energy c m -> momentum sqrt(c^2-1) m
        lam = float(lambda_fs(np.array([p * m * a]), m, np.array([a]))[0])
        out.update(lambda_collapse_Mpc=lam)
    return out


def main():
    B = benchmark(); F, Trh = B["F_GeV"], B["T_Gamma_GeV"]; lam_wdm, _ = lambda_wdm_reference()
    rows = []
    for br, q, x0 in itertools.product(BRANCHES, QS, X0S):
        for K in (0.0, 1.0):
            rows.append(solve(F, Trh, br, q, x0, K_fixed=K))
        for ce in C_ETA:
            r = solve(F, Trh, br, q, x0, c_eta=ce); r["lambda_collapse_over_WDM"] = r["lambda_collapse_Mpc"] / lam_wdm; rows.append(r)
        print(br, q, x0, [(r["c_over_eta"], round(r["K"], 2), f"{r['m_G_eV']:.2e}", round(r["T_osc_MeV"], 3), f"{r['f_em']:.2e}",
                           round(r.get("lambda_collapse_over_WDM", float('nan')), 4)) for r in rows[-5:]], flush=True)
    sel = [r for r in rows if r["c_over_eta"] is not None]
    central_over = any(r["K"] > 1 for r in sel if r["c_over_eta"] == 10.0)
    cold = all(r["lambda_collapse_over_WDM"] < 1 for r in sel)
    q087_cold = all(r["f_em"] < 0.01 for r in sel if r["q"] == 0.87)
    if not central_over: verdict = "K_SW_DENTRO_DEL_RANGO"
    elif cold and q087_cold: verdict = "K_SW_FUERA_DEL_RANGO_FRIA"
    elif not cold: verdict = "K_SW_FUERA_DEL_RANGO_TIBIA"
    else: verdict = "K_SW_FUERA_DEL_RANGO_FRIA (condición ii no cumplida: ver filas q = 0.87)"
    rng = lambda key, ce: [min(r[key] for r in sel if r["c_over_eta"] == ce), max(r[key] for r in sel if r["c_over_eta"] == ce)]
    summary = dict(verdict_U17=verdict, condition_i_collapse_cold=cold, condition_ii_q087_fem_below_1pc=q087_cold,
                   lambda_WDM_5p3keV_Mpc=lam_wdm,
                   ranges={str(ce): dict(K=rng("K", ce), m_G_eV=rng("m_G_eV", ce), T_osc_MeV=rng("T_osc_MeV", ce),
                                         M_minicluster_Msun=rng("M_minicluster_Msun", ce),
                                         lambda_collapse_over_WDM=rng("lambda_collapse_over_WDM", ce)) for ce in C_ETA},
                   ranges_Ksw_0_1=dict(m_G_eV=[min(r["m_G_eV"] for r in rows if r["c_over_eta"] is None), max(r["m_G_eV"] for r in rows if r["c_over_eta"] is None)],
                                       T_osc_MeV=[min(r["T_osc_MeV"] for r in rows if r["c_over_eta"] is None), max(r["T_osc_MeV"] for r in rows if r["c_over_eta"] is None)]),
                   B2_q1_tail={f"x{r['x0']}_c{r['c_over_eta']}_K{r['K']:.3g}": dict(f_em=r["f_em"], m_G_eV=r["m_G_eV"])
                               for r in rows if r["branch"] == "B2" and r["q"] == 1.0})
    out = dict(spec="docs/preregistros/ESPEC_U17.md", plateau=B["plateau"], F_GeV=F, T_start_GeV=Trh, R_Gamma=R_GAMMA, summary=summary, rows=rows)
    (RESULTS / "u17_ksw.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
