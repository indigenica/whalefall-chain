"""Exploratory orders of magnitude for black-hole interiors as seed factories (docs/informes/NOTA_SEMILLAS_AGUJEROS_NEGROS.md).

No verdict and no preregistration: necessary conditions and dimensional estimates only.
1. What a whalefall parent is: the turnaround state of the 30 U28 children (chi0, kinetic/dust fraction, A.H, density).
2. Turnaround of a closed dust + Lambda region (sympy): recollapse iff Lambda <= 4/(9 A^2), i.e. rho_dust >= 2 rho_Lambda.
3. Late-universe candidates: residual inflaton energy after reheating, and shear growth in a global recollapse.
4. Dimensional estimate of the kinetic fraction left by quantum fluctuations of chi at the Hubble scale of the plateau,
   f_q ~ (H^4 / 8 pi^2) / rho_dust, against the kinetic fraction of each fertile parent.

    python3 -B chain/explora_semillas_bh.py  ->  results/explora_semillas_bh.json"""
import json, math
import numpy as np
import sympy as sp
from ucommon import ROOT, benchmark
import u28_c45k as K
import u28_comoving as C


def parents():
    sel = json.loads((ROOT / "results" / "u28_runs" / "selection.json").read_text())
    B = benchmark(); rows = []
    for g, rs in sel.items():
        for r in rs:
            st = K.state_of(r); rho_d = 3 * C.hom_aeff(st) / st["a0"] ** 3
            fk = 0.5 * st["chid0"] ** 2 / rho_d                          # a ratio: invariant under the scale symmetry
            # the C45 states are rescaled (energies x c2, c2 = V0 / V0_BENCH_C45, as in u28_c45k.classify_cells); the quantum
            # estimate is not scale invariant, so it is evaluated with the physical densities (Codex, PR 47)
            c2 = st["V0"] / B["V0_BENCH_C45"]; rho_phys = rho_d / c2; H2 = B["V0_BENCH_C45"] / 3
            fq = (H2 ** 2 / (8 * math.pi ** 2)) / rho_phys
            rows.append(dict(child=r["tag"], group=g, chi0=st["chi0"], f_kin=fk, AH=r["AH"], rho_d_over_V0=rho_d / st["V0"],
                             c2=c2, f_quantum=fq, ratio_q_over_kin=fq / fk))
    return rows


def turnaround():
    a, A, L = sp.symbols("a A Lambda", positive=True)
    f = A / a + L * a ** 2 / 3 - 1                                   # adot^2 for k = +1, 8 pi G = 1, A = rho_d a^3 / 3
    astar = sp.solve(sp.diff(f, a), a)[0]
    Lmax = sp.solve(sp.Eq(sp.simplify(f.subs(a, astar)), 0), L)[0]
    return dict(Lambda_max=str(Lmax), rho_dust_over_rho_Lambda_at_critical_turnaround=str(sp.simplify(3 * A / astar ** 3 / L)))


def late_universe():
    T_rh, M_pl, t0 = 2.3e9, 2.435e18, 4.35e17 / 6.582e-25            # GeV, GeV, GeV^-1 (13.8 Gyr)
    Gam = 3 * T_rh ** 2 / M_pl                                       # ~ H(T_rh), the inflaton decay rate at reheating
    R0, l_qt = 4.4e26, 1e3 * 1.616e-35                               # m: observable radius today and l_QT ~ 10^3 l_P
    return dict(Gamma_chi_GeV=Gam, Gamma_t0=Gam * t0, contraction_factor=R0 / l_qt, log10_shear_growth=3 * math.log10(R0 / l_qt))


def main():
    P = parents(); fk = np.array([p["f_kin"] for p in P]); fq = np.array([p["f_quantum"] for p in P])
    out = dict(note="docs/informes/NOTA_SEMILLAS_AGUJEROS_NEGROS.md", parents=P,
               summary=dict(chi0_max=max(abs(p["chi0"]) for p in P), f_kin_range=[float(fk.min()), float(fk.max())],
                            AH_range=[min(p["AH"] for p in P), max(p["AH"] for p in P)],
                            rho_d_over_V0_range=[min(p["rho_d_over_V0"] for p in P), max(p["rho_d_over_V0"] for p in P)],
                            f_quantum_range=[float(fq.min()), float(fq.max())],
                            log10_q_over_kin_range=[float(np.log10(fq / fk).min()), float(np.log10(fq / fk).max())]),
               turnaround=turnaround(), late_universe=late_universe())
    (ROOT / "results" / "explora_semillas_bh.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("summary", "turnaround", "late_universe")}, indent=1))


if __name__ == "__main__":
    main()
