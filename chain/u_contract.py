"""Unified closure contract: one row per link of the chain, merging the statuses of c45b (ballena 5), cadena_cierre,
cierre_accion, cierre_red, cierre_cosmologia and the U1 stages run here.  Output: results/CONTRATO_CIERRE.json."""
import json
from ucommon import ROOT, RESULTS, PLATEAU, load, save

S = lambda n: json.loads((RESULTS / n).read_text())


def main():
    s1, s2, s3, s4, s5 = (S(f) for f in ("s1_birth.json", "s2_pocket.json", "s3_sectors.json", "s4_dark_checks.json", "s5_closure.json"))
    s6 = S("s6_measured_closure.json") if (RESULTS / "s6_measured_closure.json").exists() else None
    _u3a = ROOT / "lya" / "results" / "U3A_RESULT.json"
    u3a = json.loads(_u3a.read_text()) if _u3a.exists() else None
    if u3a and "verdict_preregistered" not in u3a: u3a = None
    _u4a = ROOT / "lya" / "results" / "u4a" / "U4A_RESULT.json"; _u4b = ROOT / "network" / "results" / "u4" / "U4B_RESULT.json"
    u4a = json.loads(_u4a.read_text()) if _u4a.exists() else None
    if u4a and u4a.get("quick"): u4a = None
    u4b = json.loads(_u4b.read_text()) if _u4b.exists() else None
    u4b_rp3 = [b for b in (u4b or {}).get("blocks", []) if b.get("status") == "GATES_PASSED" and "R_Gamma_summary" in b]
    u4b_rp3 = max(u4b_rp3, key=lambda b: b["N"]) if u4b_rp3 else None
    _u5 = ROOT / "lya" / "results" / "u5" / "U5_RESULT.json"
    u5 = json.loads(_u5.read_text()) if _u5.exists() else None
    if u5 and u5.get("quick"): u5 = None
    _u6 = ROOT / "network" / "results" / "u4" / "U6_RESULT.json"; _u7 = ROOT / "lya" / "results" / "u7" / "U7_RESULT.json"
    u6 = json.loads(_u6.read_text()) if _u6.exists() else None
    u7 = json.loads(_u7.read_text()) if _u7.exists() else None
    _u7c = ROOT / "lya" / "results" / "u7" / "U7C_RESULT.json"
    u7c = json.loads(_u7c.read_text()) if _u7c.exists() else None
    _u10b = ROOT / "network" / "results" / "u4" / "U10B_RESULT.json"
    u10b = json.loads(_u10b.read_text()) if _u10b.exists() else None
    _u11 = ROOT / "network" / "results" / "u4" / "U11_RESULT.json"
    u11 = json.loads(_u11.read_text()) if _u11.exists() else None
    _u10 = ROOT / "network" / "results" / "u4" / "U10A_RESULT.json"
    u10 = json.loads(_u10.read_text()) if _u10.exists() else None
    _u9 = ROOT / "network" / "results" / "u4" / "U9_RESULT.json"
    u9 = json.loads(_u9.read_text()) if _u9.exists() else None
    _j = lambda *parts: (lambda f: json.loads(f.read_text()) if f.exists() else None)(ROOT.joinpath(*parts))
    u13 = _j("action", "results_u13", "U13_RESULT.json"); u14 = _j("lya", "results", "u14", "U14_RESULT.json")
    u15 = _j("lya", "results", "u15", "U15_RESULT.json"); u15b = _j("lya", "results", "u15b", "U15B_RESULT.json"); u15c = _j("lya", "results", "u15c", "U15C_RESULT.json")
    u15c_final = bool(u15c and str(u15c.get("verdict_U15c", "")).startswith("RAMA_VIVA"))      # PENDIENTE / INCONCLUSO: no ratios
    u15c_allowed = bool(u15c_final and u15c["verdict_U15c"] == "RAMA_VIVA_PERMITIDA_LYA_COMPLETO" and all("p1_over_pmax" in v for v in u15c["cases"].values())); u12 = _j("lya", "results", "u12", "U12_RESULT.json")
    aud = _j("audit", "u4_u11", "a1_network.json"); u17 = (lambda f: json.loads(f.read_text()) if f.exists() else None)(RESULTS / "u17_ksw.json"); u16 = _j("network", "results", "u16", "U16_RESULT.json")
    u17b = _j("lya", "results", "u17b", "U17B_RESULT.json"); u18 = _j("network", "results", "u18", "U18_RESULT.json"); u19 = _j("network", "results", "u19", "U19_RESULT.json")
    u17c = _j("lya", "results", "u17c", "U17C_RESULT.json"); u17d = _j("lya", "results", "u17d", "U17D_RESULT.json")
    u21 = _j("lya", "results", "u21", "U21_RESULT.json"); u22 = _j("results", "u22_pole.json"); u22b = _j("results", "u22b_births.json"); u20 = _j("results", "u20_iso.json"); u19b = _j("network", "results", "u19b", "U19B_RESULT.json"); u17e = _j("lya", "results", "u17e", "U17E_RESULT.json"); u23 = _j("results", "u23_superradiance.json"); u20b = _j("results", "u20b_iso.json"); u26 = _j("results", "u26_pole_default.json"); u17f = _j("lya", "results", "u17f", "U17F_RESULT.json"); u19c = _j("network", "results", "u19c", "U19C_RESULT.json"); u19d = _j("network", "results", "u19d", "U19D_RESULT.json"); u17g = _j("lya", "results", "u17g", "U17G_RESULT.json")
    u15d = _j("lya", "results", "u15d", "U15D_RESULT.json"); u15d_final = bool(u15d and u15d.get("converged") and str(u15d.get("verdict_U15d", "")).startswith("F0_"))      # PENDIENTE: no ratios
    u24 = _j("network", "results", "u24", "U24_RESULT.json"); u25 = _j("network", "results", "u25", "U25_RESULT.json")
    u27 = _j("results", "u27_direct_children.json"); u28 = _j("results", "u28_c45k.json")
    cc = {c["connection"]: c for c in load("cadena_cierre/contract.json")["connections"]}
    J = load("c45b/c45j_results.json")
    L = [
        dict(link="1 semilla cerrada en contracción", status="ASSUMED", evidence="docs/evidencia/nacimiento_c45/NOTA_C45F_SEMILLAS.md; data/cadena_cierre/contract.json",
             open="origen y probabilidad de la semilla (D1); tamaño A·H ≈ 0.5–1.4 requerido por C45j"),
        dict(link="2 rebote(s) NPQT coherente(s) en el padre que se contrae",
             status=("PASS_TOY · D11: RAMA DIRECTA VIABLE, ETERNA MARGINAL (U28, χ inhomogéneo)" if u28 and u28["viable"]["D"]["verdict_D11"] == "VIABLE"
                     and u28["viable"]["E"]["verdict_D11"] == "MARGINAL" else "PASS_TOY / MARGINAL (D11)"),
             evidence=(f"C45d–C45j; U1-S1 reproduce los {s1['n_viable_C45j']} hijos viables bit a bit; f_viable = {J['frac_viable_LH1_over_fertile']:.3f}"
                       if PLATEAU == "sum" else
                       "U26-S1 (meseta de polo): S1 omitido sin gcc (el puerto C de C45 lo necesita)" if "f_ES_pole" not in s1 else
                       f"U26-S1 (meseta de polo): f_ES = {s1['f_ES_pole']:.3f} frente a {s1['f_ES_sum']:.3f} con la suma (R = {s1['R_f_ES_pole_over_sum']:.2f}) "
                       f"en las {s1['grid_pole']['n']} filas de C45h; {s1['viable_eternal_pole']} de los {s1['n_viable_C45j']} hijos viables de C45j son eternos con el polo; "
                       f"el puerto C (suma, modo 0) reproduce bit a bit; el filtro ζ₀ de C45j (f_viable = {J['frac_viable_LH1_over_fertile']:.3f}) es de la suma")
             + (f"; U27 (polvo, meseta de polo): {u27['R2_pole_D']['LH1']['f_v']:.1%} directos y {u27['R1_pole_E']['LH1']['f_v']:.1%} eternos viables de "
                f"{u27['n_window_fertile_pole']} fértiles ({u27['verdict_U27']})" if u27 else "")
             + (f"; U28 (χ inhomogéneo, un rebote): {u28['class_LH1']['D-v']['verdict']} en D y {u28['class_LH1']['E-v']['verdict']} en E; "
                f"directos {u28['viable']['D']['f_lo']:.1%}–{u28['viable']['D']['f_hi']:.1%} ({u28['viable']['D']['verdict_D11']}), "
                f"eternos {u28['viable']['E']['f_lo']:.1%}–{u28['viable']['E']['f_hi']:.1%} ({u28['viable']['E']['verdict_D11']})" if u28 else ""),
             open=("U28: un solo rebote extrapolado a varios (ε/n_rebotes), simetría esférica y una escala (L·H = 1; con L·H = 0.3, INCONCLUSA "
                   "por corridas inválidas); qué rama nos alberga es la medida (D9)" if u28 else "coherencia con χ inhomogéneo (C45k montado, sin ejecutar)")),
        dict(link="3 bomba χ + disparo por la garganta → aterrizaje χ_land", status="PASS_NUMERICS (homogéneo)",
             evidence="C45g/C45h (R_throat ≈ 1); cadena_cierre 'garganta fuerte' = " + cc["garganta fuerte -> coherencia del disparo"]["state"],
             open="lanzamiento global inhomogéneo"),
        dict(link="4 fase eterna (P_ζ > 1) en la misma meseta", status="PASS_EFFECTIVE_TOY · D3 FRÁGIL · UV (b)",
             evidence=f"U1-S2{' (polo)' if PLATEAU == 'pole' else ''}: χ_et = {s2['eternal']['chi_Pzeta_1']:.1f} con el V0 común; cota de planitud {s2['eternal']['flatness_bound_per_MPl']:.2e}/M_Pl; C45i FRÁGIL{' (calculado con la suma)' if PLATEAU == 'pole' else ''}; docs/literatura/NOTA_D3_UV.md",
             open="protección UV de la meseta (incompatible con las conjeturas de la ciénaga); medida (D9); islas (D12)"),
        dict(link="5 bolsillo: inflación observable, fin, recalentamiento σ, pivote", status="PASS_NUMERICS",
             evidence=f"U1-S2 rehace el pivote (máx. rel. {s2['pivot_max_rel']:.1e}): n_s = {s2['observables']['n_s']:.4f}, r = {s2['observables']['r']:.4f}"
                      + (f"; meseta de polo (U26): Δn_s = {s2['pole_vs_sum']['dn_s']:+.1e} y Δr/r = {s2['pole_vs_sum']['dr_rel']:+.1%} frente a la suma" if PLATEAU == "pole" else ""),
             open="precalentamiento no lineal; espectros a primer orden"),
        dict(link="6 χ → curvatura → ordenación Q en RP³", status=cc["χ -> curvatura -> ordenación Q"]["state"],
             evidence="data/cadena_cierre/contract.json (q_bridge); action/ (acción común, residuos simbólicos nulos)", open="Q espectador; normalización heredada"),
        dict(link="7 red de cuerdas ℤ₂ (cicatrices)", status="LOCAL_3D_RESOLVED · RED_REPRESENTATIVA_ABIERTA",
             evidence=f"docs/evidencia/cierre_*; data/cierre_cosmologia/closure_status.json: {s3['network']['local_certified_perimeters']} perímetros certificados en el octree 3D",
             open="red representativa con núcleos resueltos y AMR bidireccional"),
        dict(link="8 red → fuente de Goldstones libres S(t,p)",
             status=((f"MEDIDO_EN_EL_ATRACTOR (U4b, N = {u4b_rp3['N']}): R_Γ = {u4b_rp3['R_Gamma_summary']['value']:.2f} ± {u4b_rp3['R_Gamma_summary']['total']:.2f}"
                      f" (U2: 1.69 ± 0.29); R_ξ = {u4b_rp3['R_xi']['value']:.2f}"
                      + (f"; U6 ({u6['n_rp3']} RP³ + {u6['n_u1']} U(1)): R_Γ = {u6['R_Gamma']['value']:.3f} ± {u6['R_Gamma']['boot_std']:.3f}, "
                         f"R_ξ = {u6['R_xi']['value']:.2f}, q_RP3(ℓ=5) = {u6['q_RP3']['q0_at_ell5']:.2f} ± {u6['q_RP3']['q0_std']:.2f}" if u6 else "")
                      + (f"; U9 (redes densas, ξ_i = 2): R_Γ = {u9['R_Gamma_dense']['value']:.2f}, R_ξ = {u9['R_xi_dense']['value']:.2f}, "
                         f"q_RP3 = {u9['q_RP3_dense']['q0_at_ell5']:.2f}: R_ξ robusto, R_Γ y q dependen de la densidad inicial" if u9 else "")
                      + (f"; U10a (radiación masiva inicial drenada, {u10['N']}³): R_Γ = {u10['R_Gamma']['drained_dense']['value']:.2f} (densas) / "
                         f"{u10['R_Gamma']['drained_sparse']['value']:.2f} (poco densas), brecha {u10['gap_G']:.3f}; q ≈ 0.6–0.7; {u10['verdict_A']}; "
                         f"la conversión del baño era el 25–60 % de la emisión medida" if u10 and "verdict_A" in u10 else "")
                      + (f"; U11 (U6 drenado, {u11['n']['rp3']} RP³ + {u11['n']['u1']} U(1), 576³): R_Γ limpia = {u11['R_Gamma']['value']:.2f} ± {u11['R_Gamma']['boot_std']:.2f} "
                         f"({u11['verdict_R']}), R_ξ = {u11['R_xi']['value']:.2f}, q_RP3(ℓ=5) = {u11['q_RP3']['q0_at_ell5']:.2f} ± {u11['q_RP3']['q0_std']:.2f} ({u11['verdict_q']})"
                         if u11 else "")
                      + (f"; U10b (desde arriba): atractor acotado en ℓ = {u10b['window'][1]:.2f}: RP³ ∈ [{u10b['RP3']['xi_low_end']:.2f}, {u10b['RP3']['xi_ab_end']:.2f}] "
                         f"({u10b['RP3']['verdict']}), U(1) ∈ [{u10b['U1']['xi_low_end']:.2f}, {u10b['U1']['xi_ab_end']:.2f}] ({u10b['U1']['verdict']}); "
                         f"R_ξ observado {min(u10b['R_xi_bracket_end']['low'], u10b['R_xi_bracket_end']['high']):.2f} (desde arriba) y {max(u10b['R_xi_bracket_end']['low'], u10b['R_xi_bracket_end']['high']):.2f} (desde abajo), no una cota" if u10b else "")) if u4b_rp3 else
                     (f"MEDIDO_EN_COCIENTE (U2): R_Γ = {s6['R_Gamma']['value']:.2f} ± {s6['R_Gamma']['total']:.2f}; control U(1) fuera de rango"
                      if s6 else s3["network"]["source_contract_status"])),
             evidence="docs/informes/INFORME_U4.md, network/results/u4/U4B_RESULT.json; antes U2/U3 (network/results/u2, u3); data/cierre_red/source_contract.json",
             open=(("q_RP3 < 1 en ℓ ≈ 5 por muy poco (q + 2σ = 0.99) y bajando con ℓ; el estimador de q es el menos robusto (0.59 a 512³ con 2 semillas). "
                    + ("Auditoría U4–U11 (audit/u4_u11): todo se reproduce, pero q depende del convenio de la ventana en k (0.87–0.91; con media geométrica "
                       "o inicio del par, Q_COMPATIBLE_CON_1): Q_MENOR_QUE_1 no es robusto. " if aud else "") +
                    "Extrapolación al ℓ físico ≈ 25–87 (U5): con q < 1 la materia oscura sería fría y sin firma tibia") if u11 else
                   ("R_Γ limpia ≈ 1.6 ± 0.2 sin el baño masivo inicial (U10a); q depende del tratamiento del baño (0.6–1.0) y el estimador es frágil a ℓ ≈ 5; "
                    "falta repetir U6 drenado a 576³. Extrapolación al ℓ físico ≈ 25–87 (U5)") if u10 else
                   ("sistemática de protocolo en la emisión: R_Γ = 1.39 (ξ_i = 0.5) frente a 1.69 (ξ_i = 2), q = 1.03 frente a 0.74; la relajación de la radiación masiva "
                    "de la condición inicial no termina en la ventana (U9, informe). Extrapolación al ℓ físico ≈ 25–87 (U5)") if u9 else
                   "extrapolación al ℓ físico ≈ 50–85: ξ crece en la red (dξ/dℓ ≈ 0.06–0.18) y la pendiente espectral de RP³ es q ≈ 1.1–1.2 en ℓ ≈ 5.4 (U5)"
                   if u4b_rp3 else "redes en el atractor de escalado: con N = 256 no se llega ni con preparación de núcleo (U3b); hace falta N ≳ 1024. Extrapolación al ℓ físico (GHV frente a Hindmarsh)")),
        dict(link="9 especies, desintegraciones, fluido tardío", status="PASS (vacío, nivel árbol)",
             evidence=f"U1-S3: T_tardía/T_rh − 1 = {s3['late_vs_S2_rel']:.1e}; ΔN_eff homogéneo {s3['DeltaNeff_late_homogeneous_max']:.1e}",
             open="emisión por defectos, tasas térmicas"),
        dict(link="10 materia oscura pseudo-Goldstone (desalineamiento + emisión)", status=f"U1: {s5['verdict']}",
             evidence=f"U1-S4 (libro verificado, {s4['status']}) + U1-S5 ({s5['n_conservative']}/{s5['n_nodes']} nodos permitidos)",
             open=("m_G es ruptura explícita añadida (no está en la acción común); " +
                   (f"U13: {u13['verdict_U13']} (acción exactamente O(4), sin anomalía; ℤ_N con n ≈ 25 daría m_G ≈ 10⁻¹⁵ eV); " if u13 else "") +
                   ("Boltzmann completo con CLASS en U4a–U17e; isocurvatura en U20 y U20b (no limita)" if u20b else "Boltzmann completo e isocurvatura ausentes"))),
        dict(link="10b cierre con ε medida (S6, puerta P5 = proxy de free-streaming, superada por 10c/10d)", status=(f"S6: {s6['verdict']} (R_Γ de {s6['R_Gamma'].get('source', 'U2')})" if s6 else "pendiente"),
             evidence="results/s6_measured_closure.json; docs/informes/INFORME_U2.md §4",
             open="con las puertas de U1, B1/B2/B3 quedan excluidas incluso a −1σ; con la cota Lyα mixta (exploratoria) sobreviven B2 y B3 con el espectro log-plano"),
        dict(link="10c Lyman-α con transferencia de Boltzmann (U3a)",
             status=((f"U3a: {u3a['verdict_preregistered']}; condicional a la red: {u3a['verdict_conditional_on_network']['strict']}"
                      f" (estricto) / {u3a['verdict_conditional_on_network']['wide']} (amplio); validación preregistrada: {'SÍ' if u3a['validated_preregistered'] else 'NO'}")
                     if u3a else "pendiente"),
             evidence="lya/results/U3A_RESULT.json; docs/informes/INFORME_U3.md; audit/u3/AUDIT_U3.md",
             open="superada por 10d (CLASS)"),
        dict(link="10d Lyman-α con CLASS (U4a)",
             status=(f"U4a: {u4a['verdict_preregistered']}; validada (C0, C0b, C1): {'SÍ' if u4a['validated'] else 'NO'}; "
                     f"permitidas {sum(r['LYA_primary'] for r in u4a['rows'])}/{len(u4a['rows'])} filas (B2 y B3 con espectro log-plano)") if u4a else "pendiente",
             evidence="lya/results/u4a/U4A_RESULT.json; docs/informes/INFORME_U4.md; docs/preregistros/ESPEC_U4.md (enmienda 1)",
             open="filas calculadas con R_Γ de U2 (1.69 ± 0.29); con el del atractor (|ΔR_Γ| = 0.30, en el límite de ESPEC_U4) solo baja la emisión. Simulaciones de Lyα para la distribución; ramas con ξ creciente y q > 1 (U5)"),
        dict(link="10e extrapolación al ℓ físico: ξ creciente y pendiente espectral q (U5)",
             status=(f"U5a: {u5['verdict_U5a']}; U5f: {u5['verdict_U5f']} (q_crit = {u5['q_crit_B2']:.3f}); U5c isocurvatura: "
                     f"{'PASA' if u5['U5c']['passed'] else 'FALLA'}; dξ/dℓ (U(1), 768³) = {u5['U5d']['mean_slope']:.3f} ± {u5['U5d']['sem']:.3f}") if u5 else "pendiente",
             evidence="lya/results/u5/U5_RESULT.json; docs/informes/INFORME_U5.md; docs/preregistros/ESPEC_U5.md" + ("; lya/results/u12/U12_RESULT.json; docs/informes/INFORME_U12.md" if u12 else ""),
             open=((f"U12 ({u12['verdict_U12']}): con q ≤ 0.87 la rama fría sobrevive aunque ξ crezca (B4a, B4b y B1/GHV; q_crit ≈ "
                    f"{min(v for v in u12['q_crit'].values() if v):.2f}–{max(v for v in u12['q_crit'].values() if v):.2f}), con ΔN_eff hasta ≈ 0.07; con q = 1 solo sobrevive ξ constante (q_crit = 1.025, U5f). "
                    "Decide q a ℓ ≈ 70; la medida a ℓ = 5 (0.87–0.94 según el convenio, auditoría) está en la frontera" +
                    (f". U17 ({u17['summary']['verdict_U17']}): U12 supone K_sw ≤ 1, pero la red guarda en H = m_G/3 entre 100 y 2000 veces la energía del "
                     "desalineamiento (K_eff ≈ 3–770 según c/η ∈ [3, 30]); con ξ creciente, m_G cae a 10⁻²⁰–5×10⁻¹⁸ eV y los cuantos del colapso son más "
                     "tibios que la WDM de 5.3 keV, así que la rama fría con ξ creciente no está asegurada; con ξ constante sigue fría, con m_G ≈ 5×10⁻¹⁸–5×10⁻¹⁶ eV" if u17 else "")) if u12 else
                   "el sector oscuro vive solo si ξ deja de crecer y q ≤ 1.02 en ℓ ≈ 70 (disputa GHV/Buschmann frente a Hindmarsh); nuestras redes limpias (ℓ ≤ 6) "
                   "son compatibles con las dos condiciones (U10b: atractor acotado; U11: q < 1), sin poder extrapolar a ℓ ≈ 70" if (u10b and u11) else
                   "el sector oscuro vive solo si ξ deja de crecer y q ≤ 1.02 en ℓ ≈ 70 (disputa GHV/Buschmann frente a Hindmarsh); nuestras redes (ℓ ≤ 6) apuntan al lado excluido; U6")),
        dict(link="10f emisión calibrada con arXiv:2512.13653 (U7)",
             status=(f"U7: {u7['verdict_U7']}; κ(x₀=10) = {u7['U7a']['10']['kappa']:.2f}, κ(x₀=30) = {u7['U7a']['30']['kappa']:.2f}; "
                     f"{sum(bool(r.get('LYA_primary')) for r in u7['rows'])}/{len(u7['rows'])} filas permitidas"
                     + (f"; U7c (R_Γ = 1.39/1.69/2.09): {u7c['verdict_U7c']}, {sum(bool(r.get('LYA_primary')) for r in u7c['rows'])}/{len(u7c['rows'])} "
                        f"(solo δA: {sum(bool(r.get('LYA_dA_strict')) for r in u7c['rows'])}/{len(u7c['rows'])})" if u7c else "")) if u7 else "pendiente",
             evidence="lya/results/u7/U7_RESULT.json; docs/informes/INFORME_U6_U7.md; docs/preregistros/ESPEC_U7.md",
             open=("con R_Γ = 2.09 y x₀ = 30 falla (U7c): el margen depende de la R_Γ sin relajación inicial (U10a) y de x₀; la puerta de forma está extrapolada por debajo de 1 keV"
                   + (f"; U15 (Lyman-ABGD): {u15.get('verdict_U15')}" if u15 else "") + "; la puerta δA se reproduce en la auditoría U4–U11"
                   if u7c else "ν_ax es un número de U(1); la traslación a RP³ usa R_Γ (cociente de energías). Filas calculadas con R_Γ ≤ 1.53; U9 da 1.69 con redes densas (U7c)"
                   if u9 else "ν_ax es un número de U(1); la traslación a RP³ usa R_Γ (cociente de energías). Crecimiento de ξ_RP3 en la ventana (U9)")),
        dict(link="11 observables: n_s, r, ΔN_eff, ondas, minicúmulos", status="PREDICCIONES CONDICIONALES",
             evidence="docs/informes/INFORME_U1_CIERRE.md §3; docs/PREDICCIONES.md" + ("; docs/informes/INFORME_U14.md" if u14 else ""),
             open="SNR de LISA; verosimilitud CMB" + (f"; minicúmulos (U14): {u14['verdict_U14']}, f_MC ≲ {u14['min_f_MC_max']['estricta (EROS-2, 0.08)']:.2f} con EROS-2 estricta" if u14 else "")),
        dict(link="12 islas → agujeros negros → descendencia (cierre del ciclo)", status=cc["región eterna -> islas/agujeros negros -> descendencia"]["state"],
             evidence="D12 (argumento, docs/evidencia/nacimiento_c45/INFORME_C45IJ.md §4); data/cadena_cierre/contract.json", open="unión, radio crítico, función de masas"),
    ]
    out = dict(kind="unified_closure_contract", plateau=PLATEAU, cosmology_closed=False, conditional_closure=s5["verdict"],
               stage_status=dict(S1=s1["status"], S2=s2["pivot_status"], S3=s3["status"], S4=s4["status"], S5=s5["verdict"], S6=(s6["verdict"] if s6 else None), U3a=(u3a["verdict_preregistered"] if u3a else None),
                                 U4a=(u4a["verdict_preregistered"] if u4a else None), U4b=("GATES_PASSED, R_Γ en el atractor" if u4b_rp3 else None),
                                 U5a=(u5["verdict_U5a"] if u5 else None), U5f=(u5["verdict_U5f"] if u5 else None),
                                 U6=(u6["q_RP3_trend"] if u6 else None), U7=(u7["verdict_U7"] if u7 else None),
                                 U7c=(u7c["verdict_U7c"] if u7c else None),
                                 U10a=(u10.get("verdict_A") if u10 else None),
                                 U11=(f"{u11['verdict_R']}, {u11['verdict_q']}" + (" (q dependiente del convenio, auditoría)" if aud else "") if u11 else None),
                                 U12=(u12.get("verdict_U12") if u12 else None), U13=(u13["verdict_U13"] if u13 else None),
                                 U14=(u14["verdict_U14"] if u14 else None), U15=(u15.get("verdict_U15") if u15 else None),
                                 U15b=(u15b.get("verdict_U15b") if u15b else None),
                                 U15c=((u15c["verdict_U15c"] + " (Planck + BAO + Lyman-ABGD, MCMC; bordes F = 0.5 F₀ (r = 1) y 0.35 F₀ (r = 3) con la población de U17d: "
                                        + ", ".join(f"p(s=1)/p_max = {v['p1_over_pmax']:.2f}" for v in u15c["cases"].values())
                                        + "; posterior de s plana: la supresión en k = 10 h/Mpc es < 0.15 %, indistinguible de CDM; "
                                        + ("el borde 0.7–1 F₀ lo prueba U15d)" if u15d_final else "el borde actual (U17f, 0.7–1 F₀) no se ha probado)"))
                                       if u15c_allowed else (u15c.get("verdict_U15c") if u15c else None)),
                                 U15d=((u15d["verdict_U15d"] + " (Planck + BAO + Lyman-ABGD, MCMC, 4 cadenas por caso; borde superior con la población medida de U19b: "
                                        + ", ".join(f"{t['F_over_F0']:g} F₀ con X ≈ {t['X_univ']:.0f}: p(s=1)/p_max = {u15d['cases'][k]['p1_over_pmax']:.3f} ({u15d['cases'][k]['status']})"
                                                    for k, t in u15d["tables"].items())
                                        + "; umbral e⁻² = 0.135; la puerta de cribado dejaba pasar F₀ con X ≈ 193 solo por forma (δA 0.378) y Lyα completo lo excluye; "
                                        "informativo: Lyα completo corta entre δA = 0.32 y 0.38)")
                                       if u15d_final else (u15d.get("verdict_U15d") if u15d else None)),
                                 U16=((f"resolución: R_ξ {u16['verdict_resolution']['R_xi']}, R_Γ {u16['verdict_resolution']['R_Gamma']}; "
                                       f"volumen: R_ξ {u16['verdict_volume']['R_xi']}, R_Γ {u16['verdict_volume']['R_Gamma']}") if u16 and "verdict_volume" in u16 else None),
                                 U17=(u17["summary"]["verdict_U17"] if u17 else None),
                                 U17b=(u17b.get("verdict_U17b") if u17b else None),
                                 U17c=((u17c["verdict_U17c"] + f" (η_crit B4a = {u17c['eta_crit']['B4a|c=3|sigma=0.5']:.2g}, B1 ≈ {u17c['eta_crit']['B1|c=3|sigma=0.5']:.2g}; "
                                        "Lyα depende de ξℓη, no de c: la rama con ξ creciente exige que ≥ 90 % de la energía del colapso vaya a cuantos con p ≳ 10³ m_G)") if u17c else None),
                                 U17d=((u17d["verdict_U17d"] + " (balance corregido y R_ρ = 1.56: B2 con η = 1 pasa solo por forma, δA = 0.32 > 0.234; "
                                        "falla con r = a_decay/a_osc ≥ 2; η ≤ 0.3 pasa en todas las variantes)") if u17d else None),
                                 U21=((u21["verdict_U21"] + f" (rama viva fría con F ≤ {u21['bands']['P']['F_max_over_F0']:g} F₀; rama cómoda con F ≤ "
                                        f"{u21['bands']['A']['F_max_over_F0']:g} F₀; isocurvatura de la literatura; superradiancia con autointeracción por calcular; "
                                        "con F < F₀ se pierde la señal de cuerdas en LISA)") if u21 else None),
                                 U22=((u22["verdict_U22"] + " (inflación de polo: CMB equivalente, Δn_s = −1.2e-4; cambian los hijos concretos, no la estadística: ver U22b)") if u22 else None),
                                 U20=((u20["verdict_U20"] + f" (estimación propia de ruido blanco, C = 1: F ≤ {u20['variant_P']['F_max_over_F0_C1']:g} F₀; "
                                       f"m_min ≈ {u20['variant_P']['m_min_eV']['1']:.1e} eV, como la literatura)") if u20 else None),
                                 U22b=((u22b["verdict_U22b"] + f" (f_ES polo/suma = {u22b['R']:.2f}: el parche D2 se puede sustituir por inflación de polo de orden 4)") if u22b else None),
                                 U24=((u24["verdict_U24"] + f" (drenaje a 384³, 3 semillas: max |δR_Γ| = {u24['max_abs_delta_R_Gamma']:.1%}; R_ξ cambia "
                                       f"{u24['delta_vs_D0']['Dstrong']['R_xi']['delta']:+.1%} con el doble de drenaje)") if u24 and not u24.get("provisional") else None),
                                 U25=((u25["verdict_U25"] + f" (dτ/2 a 384³, 3 semillas: ΔR_Γ = {u25['delta']['R_Gamma']['delta']:+.2%}, ΔR_ξ = {u25['delta']['R_xi']['delta']:+.2%}; "
                                       "sesgo del paso de producción < 1 %)") if u25 and not u25.get("provisional") else None),
                                 U26=((u26["verdict_U26"] + f" (la meseta de la cadena es la inflación de polo de orden 4; la suma de 10 niveles queda como opción, WHALEFALL_PLATEAU=sum; "
                                       f"Δm_G ≤ {u26['rows']['C5_dark']['m_G_max_abs_rel']:.1%}, los veredictos de CLASS se arrastran)") if u26 else None),
                                 U18=((u18.get("verdict_U18") + " (con instantáneas densas los estimadores local e integral de R_Γ coinciden a < 1 %; el efecto de resolución de U16 era de muestreo)") if u18 and u18.get("verdict_U18") else None),
                                 U19=(u19.get("verdict_U19") if u19 else None),
                                 U17e=((u17e["verdict_U17e"] + f" (población del colapso medida en U19b: Lyα solo deja hasta {u17e['bands']['Y']['F_max_lya_only']:g} F₀; el corte en "
                                        f"{u17e['bands']['Y']['F_max_over_F0']:g} F₀ es del criterio de isocurvatura de U20; filas X (lectura (b)) pasan toda la rejilla, F₀ incluido)") if u17e else None),
                                 U17f=((u17f["verdict_U17f"] + " (borde superior de la banda en todo el intervalo del colapso de U19b: 0.7 F₀ con X ≈ 320, 0.85 F₀ con X ≈ 220–275, "
                                        "F₀ con X ≲ 194, este solo por forma (δA 0.35–0.38, sobre la referencia ancha); con la cota secundaria de supermasivos, 0.5 F₀ (X ≳ 275) a 0.7 F₀ (X ≲ 220); "
                                        "borde inferior 0.09 F₀ sin cambios)") if u17f else None),
                                 U23=((u23["verdict_U23"] + " (agujeros estelares excluyen 0.0625 y 0.045 F₀; 0.03 F₀ en la región menos robusta; supermasivos (Ünal et al., secundaria) "
                                       "excluyen m_G < 1.7e-17 eV, F ≥ 0.7 F₀)") if u23 and u23.get("verdict_U23") else None),
                                 U20b=((u20b["verdict_U20b"] + " (con la cota observacional f_iso < 0.0084 de Lyα la isocurvatura no recorta; el criterio de U20 era ~500 veces más estricto)") if u20b else None),
                                 U19b=((u19b["verdict_U19b"] + f" (RP³ ℓ_G = 3.2, 4 redes: Y = {u19b['blocks']['RP3_3.2']['Y']:.2f} ± {u19b['blocks']['RP3_3.2']['sigma']:.2f}, "
                                        "Y − σ = 1.15 frente a 1.14, en la frontera; colapso tardío r_W ≈ 1.8, c_p ≈ 4.3; " + (f"control U(1) de la misma red (4 redes): Y = {u19b['blocks']['U1_3.2']['Y']:.2f} ± {u19b['blocks']['U1_3.2']['sigma']:.2f}, c_p ≈ 5.8: "
                                              "RP³ rinde como U(1) por energía de cuerdas (Y_RP³/Y_U(1) ≈ 1.07) y no es más caliente; " if (u19b.get("blocks") or {}).get("U1_3.2") else "") + "el momento blando X es casi el mismo "
                                        "con los dos protocolos aunque ξ₁ varíe ×3; el veredicto es frágil: con el desalineamiento aún congelado o restando el sesgo de anarmonicidad sería INTERMEDIO; "
                                        "con ℓ_G = 2.4 (4 redes más), enmienda 4: lectura (c) mixta: X depende poco de ξ y más de ℓ (X₀ crece con ℓ₁ a ≈ 2σ; dos masas no separan las extrapolaciones); "
                                        "X del universo ≈ 170–320 si crece con ℓ (U17: 275), 31 si X₀ fuera constante; γ − 1 = −1.2 ± 2.4)")
                                       if u19b and u19b.get("verdict_U19b") and not u19b["verdict_U19b"].startswith("INDETERMINADO") else None),
                                 U19c=((u19c["verdict_U19c"] + f" (X₀(ℓ₁) con el protocolo de ℓ_G = 3.2 en {len(u19c['masses_in_fit'])} masas válidas "
                                        f"(N = {', '.join(str(n) for n in u19c['masses_in_fit'])}" + ("; la de 576³ en una A100 con JAX" if 576 in u19c["masses_in_fit"] else "") + "): X₀ = "
                                        + ", ".join(f"{u19c['masses'][str(n)]['affine']['X0']:.1f}" for n in u19c["masses_in_fit"])
                                        + f" en ℓ₁ = " + ", ".join(f"{u19c['masses'][str(n)]['ell_1']:.1f}" for n in u19c["masses_in_fit"])
                                        + f"; potencia p = {u19c['power_law']['p']:.2f} ± {u19c['power_law']['sigma_p']:.2f}"
                                        + (f" con χ²/ν = {u19c['power_law']['chi2'] / u19c['power_law']['nu']:.1f}" if u19c["power_law"]["nu"] > 0 else " (sin grados de libertad)")
                                        + "; X del universo según el modelo, informativo: " + ", ".join(f"{k} {v['X_univ']:.0f}" for k, v in u19c["models"].items())
                                        + "; 640³ y 768³ no validan (V2, RP³ denso > 5 %))") if u19c and u19c.get("power_law") else None),
                                 U19d=((u19d["verdict_U19d"] + f" (U19c con 12 redes por masa, semillas {u19d['seeds'][0]}–{u19d['seeds'][-1]}: X₀ = "
                                        + ", ".join(f"{u19d['masses'][k]['affine']['X0']:.1f}" for k in sorted(u19d["masses"], key=int))
                                        + f" en ℓ₁ = 3.6–4.4; p = {u19d['power_law']['p']:.2f} ± {u19d['power_law']['sigma_p']:.2f}; X del universo según el modelo, informativo: "
                                        + ", ".join(f"{k} {v['X_univ']:.0f}" for k, v in u19d["models"].items())
                                        + f"; los modelos que ajustan (χ²/ν ≤ 2) dan X ≥ {min(v['X_univ'] for v in u19d['models'].values() if v['nu'] > 0 and v['chi2'] / v['nu'] <= 2):.0f}"
                                        + "; X > 320 queda fuera de las filas de U17f)")
                                       if u19d and u19d.get("power_law") else None),
                                 U17g=((u17g["verdict_U17g"] + " (borde superior sobre el intervalo de X de U19d, con la tubería de U17f: "
                                        + ", ".join(f"X ≈ {162 * float(f[1:]):.0f}: {b['F_max_over_F0']:g} F₀" for f, b in u17g["bands"].items())
                                        + f"; en el X central de U19d, {u17g['F_max_at_U19d_central']:g} F₀"
                                        + ("; lo pone Lyα en todas" if all(b["F_max_lya_only"] == b["F_max_over_F0"] for b in u17g["bands"].values()) else "")
                                        + "; borde inferior "
                                        + f"{min(b['F_min_contiguous_over_F0'] for b in u17g['bands'].values()):g}–{max(b['F_min_contiguous_over_F0'] for b in u17g['bands'].values()):g} F₀ (superradiancia estelar))")
                                       if u17g and u17g.get("bands") and u17g.get("pipeline_ok") else (u17g.get("verdict_U17g") if u17g else None)),
                                 U10b=(f"RP3 {u10b['RP3']['verdict']}, U1 {u10b['U1']['verdict']}" if u10b else None),
                                 U9=(u9["classification"] + " (dentro de la ventana ℓ ≤ 5.66)" if u9 else None)),
               decisive_open_link=("10e: banda combinada preregistrada (U17e + U20b + U23 + U17f): rama viva con F ≈ 0.09–0.7 F₀ en todo el intervalo del colapso "
                                   "(X del universo ≈ 170–320); hasta 0.85 F₀ con X ≲ 275 y hasta F₀ con X ≲ 194, este solo por forma; la extrapolación ∝ ℓ (X ≈ 195) cae en la frontera. "
                                   "Con la cota secundaria de supermasivos, 0.09–0.5 F₀ (0.7 F₀ con X ≲ 220). Materia oscura fría dominada por el colapso; LISA solo si F₀ y sin la cota de supermasivos. "
                                   + (("La verosimilitud completa de Lyα (U15c) permite los bordes 0.35–0.5 F₀, indistinguibles de CDM; "
                                       + (f"con la población medida, U15d: {u15d['verdict_U15d']} (F₀ excluido también con X ≈ 193; con X ≈ 320, 0.85 F₀ excluido y 0.7 F₀ permitido); "
                                          "las filas que pasan la puerta de cribado solo por forma (δA ≳ 0.34) quedan dudosas, no el X central (0.5 F₀, δA 0.29). "
                                          if u15d_final else "el borde 0.7–1 F₀ con la población medida está sin probar (una U15d). ")) if u15c_allowed else "")
                                   + ("La extrapolación del colapso a ℓ ≈ 90 (U19d, cuatro masas con 12 redes cada una): " + u19d["verdict_U19d"]
                                      + f"; con X₀ ∝ ℓ₁, lineal o en potencia, X ≈ {u19d['models']['X0_proporcional_a_l1']['X_univ']:.0f}, {u19d['models']['X0_lineal_en_l1']['X_univ']:.0f} "
                                      + f"y {u19d['models']['potencia_libre']['X_univ']:.0f}, todos > 194: F₀ queda fuera; el borde superior es 0.85 F₀ (X ≲ 275), 0.7 F₀ (X ≈ 320)"
                                      + ((f"; U17g: {u17g['verdict_U17g']}, 0.6 F₀ con X ≈ 400–490, {u17g['F_max_at_U19d_central']:g} F₀ en el X central (556) y 0.35 F₀ con X ≈ 900–1370. "
                                          "Abiertos: C de la isocurvatura") if u17g and u17g.get("bands") and u17g.get("pipeline_ok") else
                                         "o menor si X ≫ 320 (sin filas de U17f; una U17g). Abiertos: C de la isocurvatura")
                                      if u19d and u19d.get("models") else
                                      "Decide la extrapolación a ℓ ≈ 90, que U19c (cuatro masas con el mismo protocolo) deja " + u19c["verdict_U19c"].split(" (")[0]
                                      + ": la palanca ℓ₁ = 3.6–4.4 y la dispersión entre masas no fijan el exponente (con X₀ ∝ ℓ₁, X ≈ "
                                      + f"{u19c['models']['X0_proporcional_a_l1']['X_univ']:.0f}, en la frontera de F₀). Abiertos: C de la isocurvatura"
                                      if u19c and u19c.get("models") else
                                      "Decide la extrapolación a ℓ ≈ 90: una tercera masa de U19b con el mismo protocolo (ℓ_G ≈ 3.6 a 768³; docs/NOTA_TERCERA_MASA_U19.md). Abiertos: C de la isocurvatura")
                                   + ("" if u15c_final else ", U15c")
                                   + ("; la meseta es la inflación de polo (U26)" if PLATEAU == "pole" else "")
                                   if (u17f and u17e and u23 and u20b) else
                                   "10e: banda combinada preregistrada (U17e + U20b + U23): rama viva con F ≈ 0.09–0.7 F₀ (m_G ≈ 1.4e-17–8e-14 eV); con la cota secundaria y poco robusta de supermasivos, 0.09–0.5 F₀ (f_a ≈ 0.5–2.7e14 GeV, m_G ≈ 6e-17–8e-14 eV); "
                                   "materia oscura fría dominada por el colapso; sin LISA. U19b con ℓ_G = 2.4: lectura (c) mixta, la (b) no se confirma; la extrapolación del colapso a ℓ ≈ 90 "
                                   "sigue abierta (X del universo ≈ 170–320; hace falta una tercera masa y una U17f). Abiertos: C de la isocurvatura" + ("" if u15c_final else ", U15c")
                                   + ("" if (u24 and not u24.get("provisional") and u25 and not u25.get("provisional")) else ", U24/U25")
                                   + ("; la meseta es la inflación de polo (U26)" if PLATEAU == "pole" else "")
                                   if (u17e and u23 and u20b) else
                                   "10e: U19b mide un colapso de RP³ caliente (Y = 1.7 ± 0.55, COLAPSO_BLANDO en la frontera): la banda conservadora de la rama viva es "
                                   "F ≲ 0.35 F₀ (entre 0.35 y 0.5 F₀ según la época medida; falta pasar el espectro medido por CLASS, U17e). Decide la extrapolación a ℓ ≈ 90: "
                                   "si el rendimiento blando no sigue a las cuerdas (U19b §4, prueba con ℓ_G = 2.4), la rama viva sería fría incluso en F₀. "
                                   "Superradiancia con autointeracción por calcular (U23); con F < F₀ se pierde LISA; la meseta puede ser inflación de polo (U22/U22b)"
                                   if (u19b and u19b.get("verdict_U19b") and not u19b["verdict_U19b"].startswith("INDETERMINADO") and u21 and u20) else
                                   "10e: elegir F dentro de la banda de U21/U20 (F ≲ 0.5 F₀ para la rama viva, ≲ 0.35 F₀ con colapso caliente, ≲ 0.18 F₀ para la cómoda) "
                                   "y comprobar ahí la superradiancia con autointeracción (U23); el colapso de RP³ (U19b) fija la variante. Con F < F₀ se pierde LISA. "
                                   "Eslabón 1–2: la meseta puede ser inflación de polo (U22/U22b: mismo CMB y misma estadística de nacimiento)" if (u17d and u21 and u22) else
                                   "10e: que ξ se estabilice a ℓ ≈ 70 y que el colapso cuerdas–paredes de RP³ sea compatible con Lyα e isocurvatura "
                                   "(U17d: con el balance corregido y R_ρ = R_Γ, ξ constante con colapso tipo U(1) pasa solo por forma y sin margen, y falla si el colapso "
                                   "termina en H ≲ m_G/12; con η ≲ 0.3 pasa siempre; ξ creciente exige η ≲ 0.1). Falta medir en la red de RP³ la fracción η de la energía "
                                   "del colapso en cuantos de p ~ m_G y la época del colapso; isocurvatura por calcular (U20)" if (u17 and u17b and u17c and u17d) else
                                   "10e: que ξ se estabilice a ℓ ≈ 70 (U17b/U17c: con el colapso cuerdas–paredes contado vive ξ constante, como materia oscura casi tibia "
                                   "equivalente a WDM de ~7–8 keV con independencia de c; ξ creciente solo vive si el colapso de RP³ manda ≥ 90 % de su energía a cuantos con p ≳ 10³ m_G, "
                                   "lo que U(1) no hace) y la fracción de la energía del colapso de RP³ en k ≲ 10 m_G (U19: no medible en la nube); fuera del alcance de la red" if (u17 and u17b and u17c) else
                                   "10e: que ξ se estabilice a ℓ ≈ 70 (U17b: con el colapso cuerdas–paredes contado, solo vive ξ constante, como materia oscura casi tibia "
                                   "equivalente a WDM de ~8 keV si c/η ≈ 3) y c/η del colapso de RP³ (U19: no medible en la nube); fuera del alcance de la red" if (u17 and u17b) else
                                   "10e: q a ℓ ≈ 70 y la energía de la red en H ≈ m_G (U17: con K_sw coherente con ξℓ, la rama fría con ξ creciente no está asegurada; "
                                   "falta medir c/η del colapso cuerdas–paredes); fuera del alcance de la red" if (u12 and u17) else
                                   "10e: q a ℓ ≈ 70 (U12: con q ≲ 0.91 el sector oscuro sobrevive aunque ξ crezca; con q ≈ 1 necesita ξ constante); fuera del alcance de la red" if u12 else
                                   "10e: ξ constante y q ≤ 1.02 en ℓ ≈ 70 son condiciones necesarias del sector oscuro (U5); fuera del alcance de la red" if u5 else
                                   "8 → 10: extrapolación de ξ(ℓ) y de la pendiente espectral q al ℓ físico (U5)" if u4b_rp3 else
                                   "8 red → fuente (ε): medido en cociente, falta el atractor de escalado"), measured_closure=(s6["verdict"] if s6 else None),
               note_historical_labels=("S6 = FALSADA_EN_FORMA_ACTUAL es la puerta proxy antigua (antes de CLASS); el estado actual de la inferencia, "
                                       "con su grado de certeza, está en docs/ESTADO_ACTUAL.md"), requirement_on_source=s5["requirement_on_source"], predictions_K_sw_0_1=s5["predictions_K_sw_0_1"], links=L)
    if PLATEAU != "pole": out["stage_status"].pop("U26", None)       # U26 is the pole's verdict; the sum contract stays as before U26
    save("CONTRATO_CIERRE.json", out)
    print("contract written:", out["stage_status"])


if __name__ == "__main__":
    main()
