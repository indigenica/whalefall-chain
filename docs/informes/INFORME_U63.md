# INFORME U63 · El orden cúbico del cuscutón cerca de φ̇ → 0

**Preregistro:** `docs/preregistros/ESPEC_U63.md` (17 nov, 16:57), con la enmienda 1 (antes del código).
**Código:** `chain/u63_cubico_cuscuton.py`. **Diagnóstico posterior:** `results/u63_runs/diagnostico_G1.{py,out}`.
**Registro:** `results/u63_runs/produce.out` (el primer intento, que se cayó por una tolerancia de `brentq`, está en `produce_intento1_fallido.out`).

Interno: no se exporta.

---

## 2026-11-17 17:59 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO (falla G1, y la producción no terminó)

| puerta | resultado |
|---|---|
| G0, el L₂ de este código frente al de U59 | **pasa:** 2.6×10⁻¹⁵ |
| G1, la integral de L₂ sobre S³ frente a la forma de U59 | **falla:** 1.4×10³⁸ |
| G2 y G3 | no se llegaron a calcular |

**La producción se paró** porque, en alguna fila, el vacío de entrada en A1 (0.1 e-fold antes del giro de la madre) no tiene energía congelada E′ positiva ni acercándose 3 e-folds hacia el rebote. No hay filas de producción y no hay medidas.

### 2. Por qué falla G1 (diagnóstico posterior, informativo)

- **El L₂ de este código es el de U59 punto a punto,** también junto a los polos de S³, a ≤ 10⁻¹⁵. No es un error de la derivación.
- **La integral de L₂ coincide con la forma de U59 a ≤ 4×10⁻¹⁵** en la contracción a 2–4 e-folds del rebote, en la ventana y en la expansión hasta ρ_r/ρ = 10⁻⁵, ya en inflación con φ̇₀/H ≈ 10⁻⁸.
- **Solo se separa en A1, junto al giro de la madre,** con n ≥ 10: la integral directa da 10⁷³–10⁷⁶ y la de U59, 10³⁶–10³⁸.
  - Allí φ̇₀ ≈ 3×10⁻¹⁴ y φ̇₀/H ≈ 2×10⁻⁷. La madre está tan diluida (W ~ 1/a²) que φ está pegado a −φ∞ y apenas se mueve.
  - Los términos del cuscutón proporcionales a potencias de 1/φ̇₀ son entonces enormes, y se cancelan entre sí sobre la solución (las ligaduras). En doble precisión, las dos formas dan solo el ruido de esa cancelación.
- **Es el régimen φ̇ → 0 que buscaba la prueba,** pero aparece donde no lo esperaba: en el giro de la madre, no al principio de la inflación. Al principio de la inflación φ̇₀/H es igual de pequeño, pero las cantidades siguen bien condicionadas.
- **Que E′ no sea positiva** junto al giro, en algunas filas, apunta a lo mismo: en el gauge de radiación uniforme, el estado de vacío junto al giro no está bien definido.

### 3. Lo que esto enseña (lectura)

- **La pregunta del acoplamiento fuerte cerca de φ̇ → 0 tiene dos sitios,** el giro de la madre y la inflación.
  - En la inflación, las cantidades del gauge de U59 están bien condicionadas incluso con φ̇₀/H ~ 10⁻⁸. Es lo que se espera si los términos 1/φ̇₀ de este gauge son artefactos que se cancelan sobre la solución.
  - En el giro, la cancelación es tan grande (de 10⁷⁶ a lo que sea el valor físico) que el gauge con δφ ≠ 0 no sirve para medir nada.
- **La herramienta adecuada es el gauge unitario del cuscutón** (δφ = 0). Allí los términos del cuscutón son ∝ μ²φ̇₀·δN^k y se anulan cuando φ̇₀ → 0, como en Dehghani, Geshnizjani y Quintin. El acoplamiento fuerte, si lo hay, aparece como una transformación al gauge unitario no lineal: ξ⁰ = f/φ̇₀, con H·ξ⁰ o χ̇·ξ⁰ de orden uno.
- **Ese cociente f/φ̇₀ se calcula sin cancelaciones,** porque es un solo cociente y no una suma de términos enormes.

### 4. Siguiente, como propuesta (U63b, con datos nuevos)

- **Las medidas,** invariantes de gauge a primer orden y sin sumas de términos grandes:
  - ζ_u = ζ − H·f/φ̇₀ y x_u = x − χ̇·f/φ̇₀, la curvatura y el inflatón sobre las hojas de φ constante (bandas de vacío; > 1 significa no lineal);
  - la inclinación k·f/(a·φ̇₀) de esas hojas.
- **Los estados:** el vacío en A2 (2 e-folds antes de la ventana) como estado principal, y A1 solo donde E′ > 0.
- **El cúbico completo** se evaluaría en el gauge unitario, pidiendo la derivación de la acción cuadrática en ese gauge, que no está hecha. Primero bastan las medidas anteriores.
- **Lo esperado:** perturbativo en la inflación y en la ventana. En el giro de la madre, f/φ̇₀ puede ser grande, lo que diría que la teoría efectiva del cuscutón no describe bien la madre cerca de su máximo, con φ pegado a −φ∞ y φ̇₀ → 0.

---

## U63b · 2026-11-17 23:41 (Madrid): el acoplamiento fuerte en variables del gauge unitario

**Preregistro:** `docs/preregistros/ESPEC_U63b.md` (17 nov, 18:14), con la enmienda 1: las filas principales sin vacío definido son NO_DEFINIDAS, y se quita una puerta que el código tenía y el ESPEC no. La enmienda se escribió antes de mirar ninguna medida.
**Código:** `chain/u63b_gauge_unitario.py`. **Resultados:** `results/u63b_gauge_unitario.json` y `results/u63b_runs/`: las 144 filas una a una en `rows/`, más `produce.out`.
**Incidencias, declaradas:**
- el primer intento se perdió por un reinicio del contenedor (`produce_intento1_perdido.out`);
- desde entonces cada fila se guarda al terminar;
- la producción se paró una vez más por un `ValueError` que el código no atrapaba (la enmienda 1).

### 1. Veredicto preregistrado: INCONCLUSO (G1, G2 y G3 fallan, solo en n = 3×10⁴)

| puerta | resultado |
|---|---|
| G0, la normalización del vacío en A2 | **pasa:** 1.1×10⁻¹⁵ |
| G1, el wronskiano en S5 (≤ 10⁻⁶) | **falla:** 4.8×10⁻⁴ en n = 3×10⁴; para n ≤ 10⁴, ≤ 4.5×10⁻⁷ |
| G2, la tolerancia | **falla** en n = 3×10⁴ (cambio relativo 1.0); en n = 100 y 3×10³, ≤ 3×10⁻¹⁰ |
| G3, el muestreo | **falla** en n = 3×10⁴ (27 %); en n = 100 y 3×10³, 0 |

**Lo que pasa en n = 3×10⁴:** es el mismo suelo de precisión que en U62b con n = 10⁵. Son decenas de miles de oscilaciones en un camino de unos 8 e-folds, y en 5 de los 8 hijos las medidas salen de 10⁵⁷–10⁸¹. Esa fila no sirve.

**Si se descarta n = 3×10⁴** (lectura informativa, posterior al veredicto), las tres puertas pasan, y el veredicto con la regla escrita sería **ACOPLAMIENTO_FUERTE**. Lo provoca la condición infrarroja:
- en dos hijos (A·H = 100, f = 0.184 y A·H = 1 000, f = 0.898), el n de ℓ = 30 (5 247 y 12 972) está por encima de la escala fuerte (3×10³ y 10⁴);
- en todos, n_sc > n_mem.

**Las ocho filas NO_DEFINIDAS** del camino principal son n = 2 en siete hijos y n = 10 en uno: modos de escala del horizonte cuyo vacío no está definido en A2.

### 2. Lo que se ve (n ≤ 10⁴)

| hijo (A·H, f) | n_mem | n_sc | n_sc/n_mem | k_sc en el rebote (M_Pl) | n de ℓ = 30 |
|---|---|---|---|---|---|
| 30, 0.096 | 415 | 3×10³ | 7.2 | 0.070 | 1 408 |
| 48, 0.104 | 516 | 3×10³ | 5.8 | 0.056 | — |
| 100, 0.184 | 738 | 3×10³ | 4.1 | 0.039 | **5 247** |
| 100, 0.965 | 818 | 10⁴ | 12.2 | 0.12 | — |
| 300, 0.246 | 1 243 | 10⁴ | 8.0 | 0.077 | — |
| 1 000, 0.542 | 2 190 | 10⁴ | 4.6 | 0.044 | 1 851 |
| 1 000, 0.835 | 2 114 | 10⁴ | 4.7 | 0.046 | 2 588 |
| 1 000, 0.898 | 2 106 | 10⁴ | 4.7 | 0.046 | **12 972** |

Las referencias son la H máxima de la ventana, H_b = 2.05×10⁻³; k_mem ≈ 0.0096; y la escala de DGQ, (H_b²)^{1/3} = 0.016.

- **La curvatura sobre las hojas de φ constante (M1) es pequeña siempre:** ≤ 5×10⁻⁴ en todos los n y fases.
- **La medida que limita es M3,** la no linealidad del paso al gauge unitario, y crece como ~n^{1.5}. M2, la inclinación, llega a 1 solo hacia n ~ 10⁴–3×10⁴.
- **M3 pasa de 1 al entrar en el rebote,** en la contracción, a ln a ≈ ln a_b + 0.1, justo antes de la ventana, con φ̇₀/H ≈ 0.2–0.3. No pasa en los sitios donde φ̇₀ → 0:
  - **en el giro de la madre** (el camino informativo, con φ̇₀/H ≈ 2×10⁻⁵ en la contracción), las medidas son ≤ 0.3 para n ≤ 3×10³;
  - **al principio de la inflación** (φ̇₀/H hasta 10⁻¹⁰), el máximo de M3 en la expansión es 0.15–0.7 para n ≤ 10⁴.
- **La causa:** allí φ̇₀ cambia muy deprisa (φ̈₀/φ̇₀ ≫ H, en la escala del rebote). El desplazamiento al gauge unitario de un modo de n alto se vuelve comparable a ese tiempo.

### 3. Lo que significa (lectura)

- **El cuscutón no está fuertemente acoplado donde φ̇₀ → 0** (el giro de la madre, la inflación). Lo que U63 vio allí eran cancelaciones del gauge de radiación uniforme, no física.
- **La escala fuerte del cuscutón está en el rebote, en k_sc ≈ 0.04–0.12 M_Pl** (físico, en el rebote):
  - es 20–60 veces la escala del fondo (H_b ≈ 2×10⁻³);
  - es 2.5–7 veces la estimación de Dehghani, Geshnizjani y Quintin (≈ 0.016);
  - es **4–12 veces la escala de memoria** (k_mem ≈ 0.01).
- **Los modos que llevan la memoria del rebote** (n ≲ 2–3·n_mem, U62b) **están dentro de la teoría efectiva.** Los de n ≳ 4–12·n_mem cruzan el rebote fuera de ella. Su paso adiabático (U62b, N ≲ 10⁻⁴) es entonces el supuesto habitual de que la completitud UV los deja en el vacío, no una predicción del cuscutón.
- **En dos de los cinco hijos visibles,** los modos de ℓ ≈ 20–30 tienen n > n_sc. Su espectro en esos ℓ depende de ese supuesto UV. En los otros tres, todo ℓ ≤ 30 está dentro de la teoría efectiva.
- **El veredicto escrito sería ACOPLAMIENTO_FUERTE** por esa regla infrarroja, que se eligió antes de saber dónde caía n_sc. La lectura física es más matizada: es un corte de la teoría efectiva bien por encima de las escalas del rebote, como el de cualquier rebote no singular de la literatura.

### 4. Para la v4

- **Se puede decir:** «la teoría efectiva del cuscutón es perturbativa durante el nacimiento hasta k ≈ 0.04–0.1 M_Pl en el rebote, 4–12 veces la escala de los modos que conservan su memoria. Con φ̇₀ → 0 (el giro de la madre y la inflación) no hay acoplamiento fuerte».
- **Lo que queda abierto:**
  - el cúbico completo en el gauge unitario, que daría la escala con más precisión que M3;
  - los n ~ 3×10⁴ (la precisión numérica);
  - la completitud UV por encima de k_sc.

---

## Nota · 2026-11-18 10:23 (Madrid): U64 corrige la escala fuerte de U63b

U64 evaluó el orden cúbico completo en el gauge unitario, sobre la misma solución.
- **Lo que mostró:** el lapso en ese gauge, α_u = α − ḟ/φ̇₀ + fφ̈₀/φ̇₀², es pequeño (≲ 10⁻³ en banda), porque M3 se cancela con los demás términos.
- **Lo que sale del cúbico:** ρ₃ ≤ 0.53 para n ≤ 10⁴ en los 8 hijos, y llegaría a 1 hacia k ~ 0.3 M_Pl, dominado por la radiación.
- **La escala de U63b** (k_sc ≈ 0.04–0.12 M_Pl) era, por tanto, una cota pesimista y no la escala fuerte. Ver `INFORME_U64.md`.
