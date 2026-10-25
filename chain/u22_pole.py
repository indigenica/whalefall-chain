"""U22 (docs/preregistros/ESPEC_U22.md): the plateau as pole inflation of order 4, V = V0 [1 - C/sqrt(chi^2 + C^2)] (exact V0 (1 - C/chi)
tail), through the same S2 pipeline (exact endpoint, continuous reheating, pivot matching; cadena_cierre/thermal.normalize_continuous with
the potential, m_chi and predict swapped) and the S1 birth runs (C port, frw_core.c mode 4).  C = 0.3516 (shape fit of the sum's tail).
    python -B chain/u22_pole.py      ->  results/u22_pole.json"""
import json, math, sys, time
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from ucommon import ROOT, benchmark, load, save
import common, scalar_end, thermal, c46_chi_reheating_relic as C46
sys.path.insert(0, str(ROOT / "components" / "c45b_hom"))

from pole_potential import C_FIT, PolePotential  # noqa: E402  (the class moved there for U26, unchanged)
ZETA0 = 4.6e-5


def predict_pole(C):
    def predict(q, s, f, L, N):                             # same algorithm as c44i.predict, pole potential
        P = PolePotential(1., C); v = lambda x: float(P.V(x)); dv = lambda x: float(P.dV(x)); d2 = lambda x: float(P.values(x)[2])
        eps = lambda x: 0.5 * (dv(x) / v(x)) ** 2
        xend = brentq(lambda x: eps(x) - 1, 1e-3, 50)
        def Nof(x):
            pts = [xend * (x / xend) ** (k / 40) for k in range(41)]
            return sum(quad(lambda y: v(y) / dv(y), pts[k], pts[k + 1], limit=200)[0] for k in range(40))
        hi = xend * 2
        while Nof(hi) < N: hi *= 2
        xs = brentq(lambda x: Nof(x) - N, xend, hi, xtol=1e-10)
        return dict(n_s=1 - 6 * eps(xs) + 2 * d2(xs) / v(xs), r=16 * eps(xs), chi_star=xs, chi_end=xend, V_star_over_V0=v(xs))
    return predict


def s2_pole(C, rtol=2e-10):
    old = (scalar_end.Potential, thermal.predict, C46.chi_numbers)
    try:
        scalar_end.Potential = lambda V0=1., L=10: PolePotential(V0, C)
        thermal.predict = predict_pole(C); C46.chi_numbers = lambda: 1.0 / C ** 2
        return thermal.normalize_continuous(rtol)
    finally:
        scalar_end.Potential, thermal.predict, C46.chi_numbers = old


def chi_eternal_pole(V0, C):
    P = PolePotential(V0, C)
    return brentq(lambda x: 3 * math.log(float(P.V(x))) - math.log(12 * math.pi ** 2) - 2 * math.log(float(P.dV(x))), 5.0, 1e7, xtol=1e-6)


def births(C, chi_et):
    from frw_c import lib, c45_state
    B = benchmark(); H = load("c45b/c45h_results.json"); J = load("c45b/c45j_results.json")
    rows = [r for r in H["rows"] if r["frac_rad"] == 0.0]
    viable = [c for c in J["children"] if c["zeta0_max_LH1"] >= ZETA0]
    ctrl = [r for r in rows if r["fertile"]][::13]
    P = PolePotential(1., C); xend = brentq(lambda x: 0.5 * (float(P.dV(x)) / float(P.V(x))) ** 2 - 1, 1e-3, 50)
    lib().set_deform(4, C, 0.0)
    out = []
    try:
        for kind, r, c in [("viable", next(r for r in rows if (r["sep"], r["lam"], r["sigma"]) == (c["sep"], c["lam"], c["sigma"])), c) for c in viable] \
                + [("control", r, None) for r in ctrl]:
            st = c45_state(B["A_C45"], r["lam"], r["sigma"], 1, r["sep"]); o = np.zeros(12)
            code = lib().run(st["a0"], st["al"], st["V0"], st["chi0"], st["chid0"], 60 * 4.0 * st["a0"], 60.0, 1e-11, 400_000_000, xend, o, 0.0)
            c2 = st["V0"] / B["V0_BENCH_C45"]; land = abs(o[2])
            out.append(dict(kind=kind, sep=r["sep"], lam=r["lam"], sigma=r["sigma"], code=int(code), n_b=int(o[0]), n_b_sum=r["n_bounces"],
                            chi_land=land, chi_land_sum=r["chi_land"], E_phys_max=o[10] / c2,
                            fertile=bool(code == 1 and o[1] > 5 and o[4] < 1e-3), eternal=bool(land > chi_et and o[10] / c2 < 1)))
    finally:
        lib().set_deform(0, 1e30, 0.0)
    return out, xend


def main():
    t0 = time.time(); B = benchmark()
    th = load("cadena_cierre/thermal.json")["continuous"] if (ROOT / "data" / "cadena_cierre" / "thermal.json").exists() else None
    out = dict(spec="docs/preregistros/ESPEC_U22.md", C_fit=C_FIT, reference_sum=dict(N_star=th["N_star"], n_s=th["n_s"], r=th["r"],
               T_rh_GeV=th["T_rh_GeV"], H_star_GeV=th["H_star_GeV"], F_GeV=th["F_GeV"], m_chi_GeV=th["m_chi_GeV"]), variants={})
    for lab, C in (("primary", C_FIT), ("C_x0.8", 0.8 * C_FIT), ("C_x1.25", 1.25 * C_FIT)):
        nr = s2_pole(C); chi_et = chi_eternal_pole(nr["V0"], C)
        v = dict(C=C, N_star=nr["N_star"], n_s=nr["n_s"], r=nr["r"], T_rh_GeV=nr["T_rh_GeV"], H_star_GeV=nr["H_star_GeV"], F_GeV=nr["F_GeV"],
                 m_chi_GeV=nr["m_chi_GeV"], V0=nr["V0"], chi_et=chi_et, dn_s=nr["n_s"] - th["n_s"], dr_rel=nr["r"] / th["r"] - 1)
        if lab == "primary":
            b, xend = births(C, chi_et); vb = [x for x in b if x["kind"] == "viable"]
            v.update(chi_end_pole=xend, births=b, n_viable=len(vb), n_viable_eternal=sum(x["eternal"] and x["fertile"] for x in vb),
                     viable_chi_land=[round(x["chi_land"], 1) for x in vb], max_rel_change_chi_land=max(abs(x["chi_land"] / x["chi_land_sum"] - 1) for x in b),
                     n_b_changed=sum(x["n_b"] != x["n_b_sum"] for x in b))
        out["variants"][lab] = v
        print(lab, {k: (round(x, 6) if isinstance(x, float) else x) for k, x in v.items() if k not in ("births", "viable_chi_land")}, flush=True)
    p = out["variants"]["primary"]
    cmb = abs(p["dn_s"]) <= 0.002 and abs(p["dr_rel"]) <= 0.10
    out["verdict_U22"] = ("POLO_CAMBIA_CMB" if not cmb else "POLO_EQUIVALENTE" if p["n_viable_eternal"] >= 15 else "POLO_CAMBIA_NACIMIENTO")
    out["runtime_s"] = time.time() - t0
    save("u22_pole.json", out); print("U22:", out["verdict_U22"], flush=True)


if __name__ == "__main__":
    main()
