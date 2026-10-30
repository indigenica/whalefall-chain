"""The result numbers quoted in main-v3.tex, recomputed from the versioned result files and looked up in the text.

    python3 -B paper/zenodo_whale/paper_numbers_v3.py      # exit 1 if anything does not match

Each row computes a value from the stored results, formats it as the text does, and checks that the LaTeX fragment
is in main-v3.tex (whitespace-insensitive). It covers the dark-sector results of version 2 (pocket inflation, the
collapse, the band and the full likelihood) and the coherence of the birth added in version 3 (U27, U28); the network numbers are those of the companion paper, checked by
paper/prd_network/paper_numbers.py.  Unlike that script it does not search for unchecked numbers."""
import json, math, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEX = re.sub(r"\s+", " ", (HERE / "main-v3.tex").read_text())


def J(*p):
    return json.loads(ROOT.joinpath(*p).read_text())


S1, S2 = J("results", "s1_birth.json"), J("results", "s2_pocket.json")
U19B, U19D = J("network", "results", "u19b", "U19B_RESULT.json"), J("network", "results", "u19d", "U19D_RESULT.json")
U17F, U17G = J("lya", "results", "u17f", "U17F_RESULT.json"), J("lya", "results", "u17g", "U17G_RESULT.json")
U15B, U15C, U15D = (J("lya", "results", s, f) for s, f in (("u15b", "U15B_RESULT.json"), ("u15c", "U15C_RESULT.json"), ("u15d", "U15D_RESULT.json")))
U22B = J("results", "u22b_births.json")

ob = S2["observables"]; pv = S2["pole_vs_sum"]
d = U19D["power_law"]; mods = U19D["models"]; ms = sorted(U19D["masses"], key=int)


def row(res, fam, F):
    return next(r for r in res["rows"] if r["kind"] == fam and abs(r["F_over_F0"] - F) < 1e-9)


def e(x, d=1):
    m, k = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(k)}}}"


def band_row(res, fam):
    b = res["bands"][fam]; up = row(res, fam, b["F_max_over_F0"]); lo = row(res, fam, b["F_min_contiguous_over_F0"])
    return f"{up['deltaA']:.3f}", e(up["m_G_eV"]), e(lo["m_G_eV"])


CHECKS = [
    ("n_s", f"n_s = {ob['n_s']:.4f}"),
    ("r", f"r = {ob['r']:.4f}"),
    ("Planck pull", f"$+{ob['pull_Planck18']:.1f}\\sigma$ from Planck 2018"),
    ("ACT pull", f"${ob['pull_ACT_DR6']:.1f}\\sigma$ from ACT DR6"),
    ("N_*", f"N_* = {S2['pivot_rerun']['N_star']['rerun']:.1f}"),
    ("T_rh", f"T_{{\\rm rh}}\\approx{S2['pivot_rerun']['T_rh_GeV']['rerun'] / 1e9:.1f}\\times10^{{9}}"),
    ("Delta n_s pole", f"\\Delta n_s={pv['dn_s']:.1e}".replace("e-04", "\\times10^{-4}")),
    ("Delta r pole", f"\\Delta r/r={pv['dr_rel'] * 100:.1f}\\%"),
    ("T_rh change", f"reheating temperature (${pv['rel_change']['T_rh_GeV'] * 100:+.0f}\\%$)"),
    ("m_chi change", f"inflaton mass (${pv['rel_change']['m_chi_GeV'] * 100:+.0f}\\%$)"),
    ("f_ES pole/sum", f"${S1['f_ES_pole']:.3f}$ with the pole against ${S1['f_ES_sum']:.3f}$"),
    ("eternal children pole", f"({S1['viable_eternal_pole']} of {S1['n_viable_C45j']},"),
    ("chi_et pole", f"\\chi_{{\\rm et}}\\simeq {S1['chi_et_pole']:.0f}\\,\\Mpl"),
    ("chi_et sum", f"${S1['chi_et_sum']:.0f}\\,\\Mpl$ with the hierarchy"),
    ("flatness ratio 1e3", f"by a factor $\\approx{1 / S2['eternal']['ratio_over_bound_at']['1e3']:.0f}$ at $\\chi=10^3\\Mpl$"),
    ("f_viable", f"${S1['f_viable_C45j'] * 100:.1f}\\%$ of fertile children bounce coherently"),
    ("U19b Y", f"$Y={U19B['blocks']['RP3_3.2']['Y']:.2f}\\pm{U19B['blocks']['RP3_3.2']['sigma']:.2f}$"),
    ("U19b U(1) Y", f"$Y={U19B['blocks']['U1_3.2']['Y']:.2f}\\pm{U19B['blocks']['U1_3.2']['sigma']:.2f}$"),
    ("X0 per mass", "X_0 = " + ",\\ ".join(f"{U19D['masses'][k]['affine']['X0']:.2f}\\pm{U19D['masses'][k]['affine']['X0_err']:.2f}" for k in ms)),
    ("ell_1", "($\\ell_1=" + "$, $".join(f"{U19D['masses'][k]['ell_1']:.2f}" for k in ms) + "$;"),
    ("p", f"p={d['p']:.2f}\\pm{d['sigma_p']:.2f}"),
    ("X_univ", f"\\Xu \\simeq {d['X_univ']:.0f}\\qquad ({d['X_lo']:.0f}\\text{{--}}{d['X_hi']:.0f}"),
    ("string term", f"least favourable string term (${d['strings_term']:.0f}$)"),
    ("models", f"$\\Xu\\simeq{mods['X0_proporcional_a_l1']['X_univ']:.0f}$ ($X_0\\propto\\ell_1$) and ${mods['X0_lineal_en_l1']['X_univ']:.0f}$ (linear"),
    ("constant model", f"a constant $X_0$ ($\\Xu\\simeq{mods['X0_constante']['X_univ']:.0f}$) is disfavoured by $\\Delta\\chi^2\\approx{mods['X0_constante']['chi2'] - mods['X0_proporcional_a_l1']['chi2']:.0f}$"),
    ("validation", f"${U15B['validation']['U15b_val_thermal_1keV']['f_95']:.2f}$ (1\\,keV) and ${U15B['validation']['U15b_val_thermal_2keV']['f_95']:.2f}$ (2\\,keV)"),
]
for k in ("U15c_B2_F0.5_r1", "U15c_B2_F0.35_r3"):
    CHECKS.append((f"{k} p1", f"{U15C['cases'][k]['p1_over_pmax']:.2f}"))
for k, t in U15D["tables"].items():
    c = U15D["cases"][k]; b = c["p1_over_pmax_bins"]; jk = c["p1_over_pmax_jackknife"]
    F = "$F_0$" if t["F_over_F0"] == 1 else f"${t['F_over_F0']:g}F_0$"
    lab = "excluded" if c["status"] == "EXCLUIDA" else "allowed"
    CHECKS.append((f"{k} row", f"{F} & {t['X_univ']:.0f} & {t['T2_at_s1_k10']:.3f} & ${c['p1_over_pmax']:.3f}$ & ${b['20']:.3f}$ / ${b['45']:.3f}$ & "
                   f"${min(jk):.3f}$--${max(jk):.3f}$ & {lab}"))
for fam, X in (("Y1.36", "217"), ("Y1.70", "272"), ("Y1.98", "315")):
    dA, mu, ml = band_row(U17F, fam); CHECKS.append((f"band {fam}", f"$\\simeq{X}$ & ${U17F['bands'][fam]['F_max_over_F0']:g}F_0$ & ${dA}$"))
    CHECKS.append((f"band {fam} masses", f"${mu}$ & ${U17F['bands'][fam]['F_min_contiguous_over_F0']:g}F_0$ & ${ml}$ & ${U17F['bands_with_SMBH'][fam]['F_max_over_F0']:g}F_0$"))
for fam, X in (("Y2.50", "396"), ("Y3.00", "477"), ("Y3.43", "541"), ("Y4.25", "674"), ("Y5.50", "859"), ("Y8.44", "1333")):
    dA, mu, ml = band_row(U17G, fam); Fm = U17G["bands"][fam]["F_max_over_F0"]
    Ftxt = f"\\mathbf{{{Fm:g}F_0}}" if fam == "Y3.43" else f"{Fm:g}F_0"
    CHECKS.append((f"band {fam}", f"$\\simeq{X}$ & ${Ftxt}$ & ${dA}$"))
    CHECKS.append((f"band {fam} masses", f"${mu}$ & ${U17G['bands'][fam]['F_min_contiguous_over_F0']:g}F_0$ & ${ml}$ & ${U17G['bands_with_SMBH'][fam]['F_max_over_F0']:g}F_0$"))
c3 = row(U17G, "Y3.43", 0.5); c3l = row(U17G, "Y3.43", 0.0625)
# reply to the review of 24 October: pole coordinate, l of the collapse, miniclusters (literature scaling of C35b), Lambda_H, Z_N order
Cp, ce = 0.3516, S1["chi_et_pole"]
CHECKS.append(("pole phi range", f"$\\phi\\simeq{Cp / 1e3 * 1e4:.1f}$--${Cp / ce * 1e4:.1f}\\times10^{{-4}}$"))
CHECKS.append(("pole c2 slope", f"$2c_2C^2/\\chi^3\\simeq{2 * Cp ** 2 / ce ** 3 * 1e9:.1f}\\times10^{{-9}}\\,c_2$"))
CHECKS.append(("ell at F0 and lower edge", f"from $\\simeq{row(U17G, 'Y3.43', 1.0)['ell_osc']:.0f}$ at $F_0$ to $\\simeq{c3l['ell_osc']:.0f}$ at the lower edge"))
CHECKS.append(("exponent effect", f"at most ${math.ceil((1 - (c3l['ell_osc'] / row(U17G, 'Y3.43', 1.0)['ell_osc']) ** (d['p'] - 1)) * 100):.0f}\\%$ lower"))
Mmc = lambda m: 2.3e-7 * (m / 1e-10) ** -1.5
CHECKS.append(("minicluster central", f"(${Mmc(c3['m_G_eV']) / 10 ** math.floor(math.log10(Mmc(c3['m_G_eV']))):.0f}\\times10^{math.floor(math.log10(Mmc(c3['m_G_eV'])))}$ to ${e(Mmc(c3l['m_G_eV']), 0)}\\,M_\\odot$ at the central $\\Xu$)"))
lam = lambda r: (r["m_G_eV"] * r["f_a_GeV"] * 1e9) ** 0.5 / 1e3
CHECKS.append(("Lambda_H", f"\\approx{lam(c3):.0f}$--${lam(c3l):.0f}$\\,keV"))
nZ = lambda r: math.log((r["m_G_eV"] * r["f_a_GeV"] * 1e9) ** 2 / 2.435e27 ** 4) / math.log(r["F_GeV"] / 2.435e18)
H45 = J("data", "c45b", "c45h_results.json")
fw = [r for r in H45["rows"] if r["frac_rad"] == 0.0 and r["fertile"] and r["sep"] in (1e-4, 1e-5, 1e-6)]
# chi_land in the archived grid is the field after ~60 e-folds; the post-bounce landing is chi_abs_max (Codex, PR 34)
fr = lambda hi: sum(1 for r in fw if r["chi_abs_max"] < hi) / len(fw)
CHECKS.append(("direct births", f"${S1['f_ES_window_common_V0'] * 100:.0f}\\%$ of fertile births land above $\\chi_{{\\rm et}}$ below the Planck density, while for ${fr(50) * 100:.0f}\\%$ the field reaches at most $50\\,\\Mpl$ after the bounce and for ${fr(10) * 100:.0f}\\%$ at most $10\\,\\Mpl$"))
from scipy.integrate import quad
from scipy.optimize import brentq
_V = lambda x: 1 - Cp / math.sqrt(x * x + Cp * Cp); _dV = lambda x: Cp * x / (x * x + Cp * Cp) ** 1.5
_xe = brentq(lambda x: 0.5 * (_dV(x) / _V(x)) ** 2 - 1, 1e-3, 50); _N = lambda x0: quad(lambda y: _V(y) / _dV(y), _xe, x0)[0]
assert min(r["chi_abs_max"] for r in fw) > 5
CHECKS.append(("observable bounce", f"a landing below $\\approx{brentq(lambda x: _N(x) - 60, 3, 10):.1f}\\,\\Mpl$, with at most $\\approx60$ e-folds"))
CHECKS.append(("fertility e-folds", f"exceed $5\\,\\Mpl$ after the bounce, which already gives $\\approx{_N(5.0):.0f}$ e-folds"))
al = S1["viable_ranges"]["alpha_phys"]; psi = S2["pivot_rerun"]["V0"]["rerun"] / 3
CHECKS.append(("NPQT during inflation", f"$\\alpha\\psi=\\alpha V_0/3\\lesssim{al[1] * psi * 1e4:.0f}\\times10^{{-4}}$ for the viable children, so the background is that of general relativity to $\\alpha^2\\psi^2\\lesssim{(al[1] * psi) ** 2 * 1e7:.0f}\\times10^{{-7}}$"))
CHECKS.append(("R_rho string term", f"(${d['strings_term']:.0f}$ of the central ${d['X_univ']:.0f}$)"))
# informative reading of the collapse (review of 25 October): c_p and r_W in the sparse networks of the four lattice sizes
sys.path.insert(0, str(ROOT / "network"))
import numpy as np
import u19d_analyze as U19DA
_sp = {N: [r for r in U19DA.runs(N, U19D["seeds"]) if r["protocol"] == "sparse"] for N in (256, 384, 512, 576)}
_cp = [float(np.mean([r["c_p"] for r in v])) for v in _sp.values()]; _rw = [float(np.mean([r["r_W"] for r in v])) for v in _sp.values()]
_l1 = [U19D["masses"][str(N)]["ell_1"] for N in _sp]
CHECKS.append(("c_p sparse", f"$c_p\\simeq{min(_cp):.2f}$--${max(_cp):.2f}$"))
CHECKS.append(("r_W slope", f"$r_W\\propto\\ell_1^{{{np.polyfit(np.log(_l1), np.log(_rw), 1)[0]:.2f}}}$"))
CHECKS.append(("energy budget r", f"if $r(90)\\gtrsim{d['X_univ'] / 162:.1f}$"))
CHECKS.append(("Z_N order", f"$n\\approx{nZ(c3l):.0f}$--${nZ(c3):.0f}$"))
CHECKS.append(("central m_G range", f"m_G\\approx{e(c3['m_G_eV'])}$--${e(c3l['m_G_eV'], 0)}"))
CHECKS.append(("central T_osc", f"T\\approx{c3['T_osc_MeV']:.1f}$--${c3l['T_osc_MeV']:.0f}$\\,MeV"))

# version 3: coherence of the birth (U27 dust model, U28 inhomogeneous chi)
U27, U28 = J("results", "u27_direct_children.json"), J("results", "u28_c45k.json")
pc = lambda x: f"{100 * x:.1f}\\%"
E1, D1, B1 = (U27[k]["LH1"] for k in ("R1_pole_E", "R2_pole_D", "R3_pole_D10_50"))
vE, vD = U28["viable"]["E"], U28["viable"]["D"]; cE, cD = U28["class_LH1"]["E-v"], U28["class_LH1"]["D-v"]
CHECKS += [
    ("U27 fertile", f"of the ${U27['n_window_fertile_pole']}$ fertile births of the grid"),
    ("U27 eternal", f"${pc(E1['f_v'])}$ of the ${U27['n_window_fertile_pole']}$ fertile births of the grid are eternal children"),
    ("U27 direct", f"${pc(D1['f_v'])}$ are direct children that do"),
    ("U27 band", f"${pc(B1['f_v'])}$ are direct children that also land"),
    ("U27 LH0.3", f"leaves ${pc(U27['R2_pole_D']['LH0.3']['f_v'])}$ direct and ${pc(U27['R1_pole_E']['LH0.3']['f_v'])}$ eternal children"),
    ("U28 table direct", f"Direct & ${pc(D1['f_v'])}$ & ${100 * vD['f_lo']:.1f}$--${pc(vD['f_hi'])}$ & viable"),
    ("U28 table band", f"& ${pc(B1['f_v'])}$ & $\\ge{pc(vD['band_f_lo'])}$"),
    ("U28 table eternal", f"Eternal & ${pc(E1['f_v'])}$ & ${100 * vE['f_lo']:.1f}$--${pc(vE['f_hi'])}$ & marginal"),
    ("U28 abstract", f"$1.8$--$2.9\\%$ for eternal children" if (round(100 * vE["f_lo"], 1), round(100 * vE["f_hi"], 1)) == (1.8, 2.9) else "MISMATCH"),
    ("U28 abstract direct", f"${100 * vD['f_lo']:.0f}$--${100 * vD['f_hi']:.0f}\\%$ of fertile births are direct"),
    ("U28 median D", f"lies between ${cD['median_rho_lo']:.2f}$ and ${cD['median_rho_hi']:.2f}$, compatible with one"),
    ("U28 median E", f"lies between ${cE['median_rho_lo']:.2f}$ and ${cE['median_rho_hi']:.2f}$: lower"),
    ("U28 verdicts", "\\textsc{igual} in both classes; direct branch \\textsc{viable}" if (cE["verdict"], cD["verdict"], vD["verdict_D11"], vE["verdict_D11"]) == ("IGUAL", "IGUAL", "VIABLE", "MARGINAL") else "MISMATCH"),
    ("U28 runs", f"in ${U28['n_runs']}$ runs"),
    ("U28 measure", f"(${100 * vD['f_lo']:.1f}$--${pc(vD['f_hi'])}$ of fertile births against ${100 * vE['f_lo']:.1f}$--${pc(vE['f_hi'])}$"),
    ("U28 field range", f"$10$--$50\\,\\Mpl$ for at least ${pc(vD['band_f_lo'])}$ of fertile births, still"),
]
_g = J("results", "u28_runs", "gates.json")
CHECKS += [
    ("U28 V1", f"homogeneous bounce to ${e(max(x['tB_rel'] for x in _g['V1']))}$ in the bounce time"),
    ("U28 V2", f"shell by shell to ${e(max(x['max_rel'] for x in _g['V2']))}$"),
    ("U28 V3", f"by at most ${100 * max(x['rel_change'] for x in _g['V3']):.1f}\\%$ when the resolution"),
    ("U28 V4", f"collapse to ${e(max(x['rel'] for x in _g['V4']), 0)}$"),
]

# version 3.1: open assumptions "Seed" (docs/informes/NOTA_SEMILLAS_AGUJEROS_NEGROS.md) and "Late-time acceleration".
# The check of the quantum shortfall (10^2.4--10^11) was dropped on 29 October 2026 with the sentence that quoted it:
# it held only for the grid of parents (correction of 29 October in that note)
SB = J("results", "explora_semillas_bh.json"); sm = SB["summary"]
CHECKS += [
    ("seed AH", f"($AH\\simeq{sm['AH_range'][0]:.1f}$--${sm['AH_range'][1]:.0f}$)"),
    ("seed chi0", "with $\\chi$ at rest at its minimum" if sm["chi0_max"] == 0 else "MISMATCH"),
    ("seed kinetic", f"kinetic energy of ${e(sm['f_kin_range'][0], 0)}$--${e(sm['f_kin_range'][1], 0)}$ of the dust"),
]
_H0 = 67.4e3 / 3.0857e22 * 6.582e-16; _rhoL = 0.685 * 3 * _H0 ** 2 * (2.435e27) ** 2       # eV, eV^4 (8 pi G = 1/M_Pl^2)
CHECKS.append(("Lambda shortfall", f"about {['zero','one','two','three','four','five','six','seven','eight','nine'][round(math.log10(_rhoL / (_H0 ** 2 * (2.7e23) ** 2)))]} orders of magnitude below the dark energy"))

# version 3.1: reasons of the lost cells in U28 (chain/u28_motivos.py, informative)
MO = J("results", "u28_motivos.json"); ms = MO["summary"]
_pc = lambda g, k: 100 * ms[g]["reasons"].get(k, 0) / sum(ms[g]["reasons"].values())
_runs = [json.loads(f.read_text()) for f in ROOT.joinpath("results", "u28_runs").glob("sep*.json")]
_sel = json.loads(ROOT.joinpath("results", "u28_runs", "selection.json").read_text()); _grp = {x["tag"]: g for g, v in _sel.items() for x in v}
_inc = {g: sum(1 for o in _runs if _grp[o["child"]] == g and o["LH"] == 1.0 and o.get("valid") and not o["coherent"]) for g in ("E-v", "D-v")}
CHECKS += [
    ("U28 reasons runs", f"Repeating the ${MO['n_runs']}$ runs that pass the bounce-time test" if MO["labels_mismatches"] == 0 else "MISMATCH"),
    ("U28 reasons flips", f"(${_pc('E-v', 'D'):.0f}\\%$ of the cells in both classes)" if round(_pc("E-v", "D")) == round(_pc("D-v", "E")) else "MISMATCH"),
    ("U28 reasons Planck", f"(${_pc('E-v', 'super_planck'):.0f}\\%$ for eternal and ${_pc('D-v', 'super_planck'):.0f}\\%$ for direct children)"),
    ("U28 reasons frequency", f"(${ms['E-v']['n_runs']}$ of their ${_inc['E-v']}$ incoherent runs, against ${ms['D-v']['n_runs']}$ of ${_inc['D-v']}$)"),
]

fails = 0
for name, frag in CHECKS:
    if frag == "SKIP":
        print(f"SKIP {name}"); continue
    ok = re.sub(r"\s+", " ", frag) in TEX
    fails += not ok
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"   expected: {frag}"))
print(f"{len(CHECKS)} checks, {fails} failures")
sys.exit(1 if fails else 0)
