# U1 · Cierre cosmológico unificado · Preregistro

*6 de octubre de 2026, 06:55 UTC · claude_dual/unification. Escrito antes de ejecutar código. Las cotas observacionales salen de la literatura citada, no de nuestros resultados.*

## 0. Qué se quiere saber

Las campañas de cierre (4–5 de octubre) dejaron la cadena conectada en casi todos los eslabones, pero con un eslabón sin calibrar: **cuánto emite la red RP³ de cuerdas ℤ₂ en Goldstones libres** (`source_contract.json`: NOT_CALIBRATED). La ballena 5 añadió el nacimiento por delante.

U1 no calibra la fuente, porque eso exige una red 3D representativa con núcleos resueltos. Lo que hace es **convertir la incógnita en un requisito**. Para cada eficiencia ε de la fuente, referida a la normalización asintótica GHN de cuerdas U(1) que ya usa `cadena_cierre/dark.py`, calcula:
- la masa m_G que da Ω_DM h² = 0.12;
- todos los observables tardíos.

El cierre cosmológico se declara **condicional**: existe una región de ε compatible con todas las cotas, y la cadena hace predicciones dentro de ella.

## 1. Etapas (un único conjunto de parámetros)

| Etapa | Qué recalcula | Puerta numérica |
|---|---|---|
| S1 nacimiento | Rehace con `frw_core.c` (copia) los hijos eternos sub-Planck de C45j que toleran ζ0 = 4.6e-5. Necesita gcc; si no lo hay, la etapa queda SKIPPED, no se aprueba. | χ_land y n_rebotes idénticos a C45h (diferencia relativa < 1e-8, mismo código determinista) |
| S2 bolsillo | Pivote continuo de `cadena_cierre/thermal.py` (copia) con rtol 2e-8: N\*, n_s, r, V0, T_rh, F | Concordancia con `thermal.json` < 1e-5 relativo |
| S2b eterno / D3 | Umbral P_ζ = 1 con el V0 común; condición de planitud \|V'/V\| < √(V0/3)/(2π) en la meseta sin deformar | Informativa: compara χ_et(V0 común) con χ_et(V0_BENCH) de C45 |
| S3 sectores | Lee especies (`species.json`), anchuras (`decays.json`) y fluido tardío (`late_fluid.json`); comprueba que T_final del fluido coincide con T_rh de S2 | < 1e-3 relativo |
| S4 materia oscura | Libro de `dark.py`, adaptado y vectorizado, con dos añadidos: (i) longitud comóvil de flujo libre λ_fs de la componente emitida; (ii) barrido extendido | Reproduce las raíces de `dark.json` < 1e-4 en log m |
| S5 cierre | Región permitida y predicciones | §2 |

## 2. Barrido y puertas de S5 (preregistradas)

**Malla:**
- espectro ∈ {ir_delta, energy_logflat};
- x0 ∈ {10, 30, 100};
- K_sw ∈ {0, 1, 10} (horquilla de muros de C35b sobre el desalineamiento);
- ε ∈ {0} ∪ 10^{−4 … 0} en 17 pasos (0.25 dex).

Para cada nodo se busca la raíz Ω_mis(m) + Ω_em(m, ε) = 0.12.

**Puertas por nodo:**
- **P1 Ω_DM:** la raíz existe en m ∈ [1e-27, 1e-10] eV.
- **P2 ΔN_eff (CMB):** momento de radiación 3P/ρ_ν1 a T = 0.26 eV, más el ΔN_eff del fluido tardío. Tiene que ser < 0.29 (Planck 2018 TT,TE,EE+lowE+lensing+BAO, N_eff = 2.99 ± 0.17 frente a 3.044, al 95 %). Es un diagnóstico, no una verosimilitud.
- **P3 ΔN_eff (BBN):** el mismo momento a T = 1 MeV tiene que ser < 0.38 (Fields et al. 2020, N_eff = 2.88 ± 0.27 al 95 %).
- **P4 Lyα, materia oscura ondulatoria:** m_G ≥ 2×10⁻²⁰ eV (Rogers & Peiris 2021, 95 %). Toda la materia oscura es esta partícula, así que la cota se aplica entera.
- **P5 flujo libre de la componente emitida:**
  - Pasa si f_em = Ω_em/0.12 ≤ 0.1.
  - Si no, pasa si λ_fs(em), promediada en masa, es menor que la λ_fs de un fermión térmico de 5.3 keV con Ω h² = 0.12 (Iršič et al. 2017, al 2σ). Las dos se calculan con el mismo integrador: ⟨p⟩ = 3.15 T_x y (T_x/T_ν)³ = 0.12·94 eV/m.
- **P6 superradiancia:** ventanas de Ünal et al. 2021 (2009.07206, 68 %): [2.5e-21, 1.2e-20], [5.5e-20, 1.3e-16] y [3.8e-14, 3.4e-11] eV.
  - Esas cotas solo valen si f_a ≳ 1.5×10¹⁵ GeV; por debajo, la autointeracción (bosenova) las apaga.
  - Nuestra f_a = F/2 ≈ 5.35×10¹⁴ GeV queda por debajo. P6 se informa como **TENSIÓN**, no como exclusión.
  - Se dan dos regiones: la conservadora (P1–P5) y la estricta (P1–P6).

**Veredicto U1:**
- **CIERRE CONDICIONAL** si la región conservadora no es vacía para al menos una combinación (espectro, x0, K_sw) con K_sw = 0. Se informa el intervalo de ε permitido y las predicciones (m_G, ΔN_eff, fracción emitida, λ_fs).
- **CIERRE ESTRECHO** si además la región estricta no es vacía.
- **SIN CIERRE** si la región conservadora es vacía en toda la malla.

**Predicciones que se informan siempre**, para cualquier ε:
- n_s y r frente a Planck 2018 y ACT DR6;
- f_a frente al umbral de LISA de GHN (5e14 GeV, nota D8);
- masa de los minicúmulos (escala de Fairbairn, como C35b);
- detectabilidad de ΔN_eff en CMB-S4 (σ ≈ 0.03: se marca si ΔN_eff ≥ 0.06).

## 3. Qué no hace U1 (declarado)

- No calibra ε. **La fuente RP³ sigue sin calibrar**; U1 dice qué valor tendría que salir.
- No cambia la acción común. La masa de los pseudo-Goldstone, m_G, es un **término de ruptura explícita añadido** (como en C35/C35b). En la acción de `cierre_accion`, los 3 Goldstone son exactamente sin masa (`species.json`).
- No resuelve D3, D11 ni la semilla. Se arrastran como hipótesis con su estado.
- El momento 3P no es N_eff en sentido de verosimilitud. Aquí no hay Boltzmann completo ni isocurvatura.
- Las ondas gravitacionales se dan solo como comparación con el umbral; no se calcula la SNR.

## Erratas posteriores (auditoría independiente, `audit_u1/AUDIT_U1.md`; añadidas después de ejecutar, sin cambiar ninguna puerta)

1. **P6, cita.** Las ventanas citadas vienen de **M. J. Stott (2020), arXiv:2009.07206**, al 68 %; las cifras coinciden exactamente. No son de «Ünal et al. 2021», que es arXiv:2012.12790. El umbral de validez tampoco es único:
   - f_a ≳ 1.2×10¹⁵ GeV para binarias de rayos X;
   - f_a ≳ 2.5×10¹⁴ GeV para fusiones de ondas gravitacionales;
   - f_a ≳ 2.0×10¹⁶ GeV para agujeros negros supermasivos.

   Con f_a = 5.35×10¹⁴ GeV, la ventana SMBH queda aún más claramente apagada. Las masas permitidas con K_sw ≤ 1 no tocan la ventana estelar.
2. **P2 y P3.** Las barras ±0.17 y ±0.27 son al 68 %. Los umbrales 0.29 y 0.38 corresponden a 2σ, así que los números se mantienen.
3. **P5.** La franquicia f_em ≤ 0.1 es la que fija el borde superior de ε en 6 de las 12 configuraciones con K_sw ∈ {0, 1}. Sin ella, ε_max(K_sw = 0) baja de 5.6e-2 a 1.8e-2.
   - Prueba de cola frente a las cotas de materia oscura mixta fría + tibia: pasa en las 110 filas permitidas, con un uso máximo del 11 % del límite.
   - Con K_sw = 10 (componente casi caliente), la franquicia es probablemente demasiado permisiva. No afecta a la cabecera.
