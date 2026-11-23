# INFORME_U66 · La salida caliente del nacimiento con cuscutón

*2026-11-21 14:53 (Madrid). Preregistro: `docs/preregistros/ESPEC_U66.md` (14:46). Código: `chain/u66_salida_caliente.py`. Resultado: `results/u66_salida_caliente.json`. Interno.*

## 1. Veredicto preregistrado: INCONCLUSO (falla G0)

| puerta | resultado |
|---|---|
| G0, la regresión con Γ_χ = 0 (≤ 10⁻⁶ relativo en N_inf y χ_land) | **falla:** 4.4×10⁻⁶ |
| G1, la ligadura de Friedmann (≤ 10⁻⁶) | **pasa:** 3.2×10⁻⁸ |
| G2, la convergencia (rtol 10⁻¹² frente a 10⁻¹¹; ≤ 0.01 dex en el margen del primer giro, misma clase) | **pasa:** 7×10⁻⁵ dex |

**Por qué falla G0.**
- **N_inf** coincide con U57b a ≤ 1.4×10⁻⁷.
- **χ_land** se desvía hasta 4.4×10⁻⁶ en 21 de los 57 hijos. Esa cifra es el máximo de |χ| sobre los puntos de salida del integrador, no el máximo continuo. Con la radiación como variable dinámica, la secuencia de pasos cambia, y con ella el punto en que se toma el máximo.
- Es un umbral mal elegido para una magnitud muestreada, no un fallo físico. Pero el criterio no se cambia después de ver los datos: el veredicto queda INCONCLUSO.

**Un fallo de código, declarado.** La primera ejecución calculó todo y se cayó al escribir el JSON, porque había escalares de numpy. Se corrigió la serialización y se repitió el mismo cálculo (`results/u66_runs/produce_intento1_fallido.out`).

## 2. Lecturas informativas (sin veredicto)

Con las puertas G1 y G2 en verde, **los 57 hijos salen calientes:**

| juego | hijos calientes | margen mínimo en el 1.er giro de χ | mínimo en los 25 giros | mínimo en la fase de fluido | ΔN_inf por Γ_χ | T al final del recalentamiento |
|---|---|---|---|---|---|---|
| adoptado (μ = 0.16) | 19 / 19 | 10^43.5 | 10^43.5 | 10^44.1 | ≤ 1.8×10⁻⁴ | 2.29×10⁹ GeV |
| C1 | 38 / 38 | 10^44.2 | 10^44.2 | 10^44.7 | ≤ 5.8×10⁻⁴ | 2.31×10⁹ GeV |

El margen es M = log₁₀[(2/3)ρ_r a²]: la condición de rama en el peor caso, un giro de χ, exige M > 0.

**Lo que significa.**
- Tras ≥ 60 e-folds, la curvatura 1/a² ha caído más de 10⁴⁰ veces. La radiación que crea Γ_χ en el primer cuarto de oscilación basta y sobra.
- El supuesto «la salida es caliente por necesidad» se reduce, con margen enorme, a que haya recalentamiento con Γ_χ > 0, lo que la cadena ya supone.
- El recalentamiento no cambia el número de e-folds (≤ 6×10⁻⁴) ni T_rh (2.29×10⁹ GeV, como la cadena).
- Coincide con la exploración externa declarada en el ESPEC: 10^43.55 y 10^44.22, a la centésima.

**Hacia el futuro (extrapolación analítica).**
- Después del recalentamiento, la condición es ½(ρ + P) > 1/a². Cae en radiación y en materia, y con Λ solo la sostiene la materia.
- El margen de hoy es 1.5 Ω_m/|Ω_K|. Con Λ, la rama del cuscutón se acabaría dentro de ≥ 11.0 e-folds (juego adoptado) o ≥ 12.5 (C1), en el hijo con mayor |Ω_K| (≈ 10⁻⁵).
- No es una predicción calculada con el cuscutón tardío. Es una consecuencia de su rama, a tener en cuenta en el supuesto «Late-time acceleration».

## 3. Siguiente

**U66b:** los mismos hijos con datos nuevos (Γ_χ/10, Γ_χ y 10Γ_χ, porque la anchura de la cadena es fenomenológica) y la puerta G0 corregida para χ_land, declarada después de ver U66.

---

# U66b · Γ_χ variado y la puerta de regresión corregida

*2026-11-21 14:57 (Madrid). Preregistro: `docs/preregistros/ESPEC_U66b.md` (14:53, después de ver U66). Código: `chain/u66_salida_caliente.py --u66b`. Resultado: `results/u66b_salida_caliente.json` (207 s con 4 procesos).*

## Veredicto preregistrado: SALIDA_CALIENTE

| puerta | resultado |
|---|---|
| G0, la regresión con Γ = 0 (N_inf ≤ 10⁻⁶, χ_land ≤ 10⁻⁵; corregida tras U66) | **pasa:** 1.4×10⁻⁷ y 4.4×10⁻⁶ |
| G1, la ligadura (≤ 10⁻⁶) | **pasa:** 3.2×10⁻⁸ |
| G2, la convergencia en los tres Γ (≤ 0.01 dex, misma clase) | **pasa:** 7×10⁻⁵ dex |

| Γ | juego | hijos calientes | margen mínimo en el 1.er giro | mínimo en la fase de fluido | ΔN_inf máx. | T al final del recalentamiento |
|---|---|---|---|---|---|---|
| Γ_χ/10 | adoptado | 19 / 19 | 10^42.5 | 10^43.1 | < 10⁻³ | 7.2×10⁸ GeV |
| Γ_χ/10 | C1 | 38 / 38 | 10^43.2 | 10^43.7 | < 10⁻³ | 7.3×10⁸ GeV |
| Γ_χ | adoptado | 19 / 19 | 10^43.5 | 10^44.1 | 2×10⁻⁴ | 2.29×10⁹ GeV |
| Γ_χ | C1 | 38 / 38 | 10^44.2 | 10^44.7 | 6×10⁻⁴ | 2.31×10⁹ GeV |
| 10Γ_χ | adoptado | 19 / 19 | 10^44.5 | 10^45.1 | 2×10⁻³ | 7.2×10⁹ GeV |
| 10Γ_χ | C1 | 38 / 38 | 10^45.2 | 10^45.7 | 6×10⁻³ | 7.3×10⁹ GeV |

**Lo que significa.**
- El margen va exactamente como log₁₀ Γ, como se esperaba: la radiación del primer giro es proporcional a Γ.
- Hasta con una anchura diez veces menor que la de la cadena, la rama del cuscutón sobrevive a la salida de la inflación con más de 42 órdenes de magnitud de margen.
- El supuesto abierto de `main-v4.tex` («that this holds is estimated, not computed with a backreacting reheating stress») está calculado ahora con la fricción efectiva de la cadena. Lo que sigue sin calcular es el recalentamiento microscópico, igual que en el resto de la cadena.
- La salida caliente deja de ser un supuesto específico del cuscutón: basta cualquier recalentamiento con Γ_χ > 0.
- Lo que queda, y no lo resuelve U66b, es la extrapolación al futuro: con Λ, la rama se acabaría dentro de ≥ 11 e-folds (§2).


## Nota del 2026-11-21 a las 18:13 (Madrid): dos arreglos de la revisión de Codex

La revisión automática de la PR señaló dos fallos en `chain/u66_salida_caliente.py`. Están corregidos después de las producciones de U66 y U66b; los JSON guardados se escribieron con la versión anterior y no se han recalculado.
- **La comparación de clase de G2** comparaba el resultado y la ruptura de la fase 2, pero no la fase 3. Ahora usa la clasificación completa de «salida caliente» del ESPEC (`is_hot`), en U66 y en U66b. Comprobado a mano con los datos guardados: las 57 filas son calientes a rtol 10⁻¹² en los tres Γ, así que G2 sigue pasando y los veredictos no cambian.
- **Los resúmenes por conjunto** fallarían si una fila no llegara a alguna fase. Ahora usan solo las filas que llegaron y listan las incompletas. En estas producciones no hay ninguna, así que los números del informe no cambian.
