# INFORME U62 · ¿Excita el rebote con cuscutón los modos observables?

**Preregistro:** `docs/preregistros/ESPEC_U62.md` (16 nov 2026, Madrid), con la enmienda 1 (cómo se cuentan los hijos, qué es N_χ y las coordenadas), escrita antes del código.
**Código:** `chain/u62_bogoliubov_rebote.py`, que importa la maquinaria de U59 sin cambios. El diagnóstico posterior al veredicto está en `chain/u62_diagnostico.py`.
**Resultados:** `results/u62_bogoliubov_rebote.json`, `results/u62_runs/produce.out` y `results/u62_runs/diagnostico.{json,out}`.
**Coste:** la producción tardó 26 s.

Interno: no se exporta.

---

## 2026-11-16 09:06 (Madrid)

### 1. Veredicto preregistrado: REBOTE_EXCITA

| puerta | resultado |
|---|---|
| G0, la normalización de los modos de entrada (≤ 10⁻¹⁰) | **pasa:** 6×10⁻¹³ |
| G1, la identidad de Bogoliubov (≤ 10⁻⁶) | **pasa:** 1.8×10⁻¹⁰ |
| G2, la convergencia en rtol (10⁻⁹ frente a 10⁻¹⁰) | **pasa:** como mucho el 7×10⁻⁵ de la tolerancia, en n = 10², 10⁴ y 10⁵ |

El montaje, en breve:
- los modos salen del vacío instantáneo en φ = −0.2;
- cruzan el rebote con la acción de U59;
- en φ = +0.2 se descomponen en los modos instantáneos de salida.

N es el número de cuantos creados en cada modo de salida.

**N del modo de mayor frecuencia (el que el ESPEC llamaba «de χ») en los n_obs de los cinco hijos visibles:**

| hijo (A·H, f) | Ω_K,0 | ℓ = 2 | ℓ = 10 | ℓ = 30 |
|---|---|---|---|---|
| 30, 0.096 | −4.7×10⁻⁵ | 0.55 (n = 94) | 0.067 (469) | 1.49 (1 408) |
| 100, 0.184 | −3.4×10⁻⁶ | 0.39 (350) | 0.79 (1 749) | 1.96 (5 247) |
| 1 000, 0.542 | −2.7×10⁻⁵ | 0.62 (123) | 0.53 (617) | 0.12 (1 851) |
| 1 000, 0.835 | −1.4×10⁻⁵ | 0.60 (173) | 0.43 (863) | 0.086 (2 588) |
| 1 000, 0.898 | −5.5×10⁻⁷ | 0.42 (865) | 0.61 (4 324) | 2.34 (12 972) |

- **Los cinco hijos cuentan como excitados** (N ≥ 0.1 en algún ℓ), y tres de ellos en los tres ℓ.
- **Los factores de potencia de ese modo** van de 0.09–11 (en el peor caso, N = 2.3) a 0.8–1.5 (N = 0.07).
- **En todos los n del barrido** (10² a 10⁵), N oscila entre 0.005 y 2.3, sin decaer con n. Es la ganancia de energía de U59b vista como creación de partículas.

**El otro modo** (el de menor frecuencia) se excita poco:
- N ≤ 0.08 en todos los n_obs;
- y cae a 10⁻³–10⁻⁵ para n ≳ 2×10³.
- **Con la misma regla, sería INTERMEDIO:** ningún N ≥ 0.1, pero no todos ≤ 0.01.

### 2. Lo que el diagnóstico posterior cambia en la lectura (informativo, hecho después de ver el veredicto)

**a) El modo de mayor frecuencia no es el de χ.** El ESPEC suponía que, en los bordes de la ventana, las frecuencias serían ≈ k (χ) y ≈ k/√3 (la radiación), y llamó «de χ» al de mayor frecuencia. No es así:

| | modo rápido | modo lento |
|---|---|---|
| ω/k para n ≳ 10⁴ | ≈ 1.88–1.95 | ≈ 0.99–1.00 |
| fracción en δχ, ponderada con la energía cinética | 0.42–0.48 | 0.95–1.00 |

- Para n bajos, las frecuencias las domina un término de masa o de curvatura: ω/k llega a 80 y 11 en n = 123.
- **El modo que se parece a χ es el lento,** que es casi δχ puro, con velocidad 1. El rápido es una mezcla a medias de δχ y de ζ, con una velocidad de propagación ≈ 1.9 cerca del rebote. Esa velocidad superlumínica también habrá que entenderla.
- **La identificación del ESPEC era, por tanto, errónea en su premisa.** El veredicto se queda como está escrito (no se cambia el criterio después de ver los datos), pero su lectura física es otra: **el rebote excita a orden uno el modo mixto rápido, y apenas el modo de χ.**

**b) En los ℓ bajos, el vacío instantáneo no está bien definido en los bordes de la ventana.**
- El parámetro de adiabaticidad, |ω̇|/ω², llega a 0.1–1.2 en el modo lento cuando k/aH ≲ 5 al final de la ventana. Eso pasa en casi todos los ℓ = 2 y ℓ = 10.
- Allí los N dependen de la prescripción del estado, como ya ocurría en U56.
- **Las filas limpias** (adiabaticidad ≤ 0.1 en los dos bordes y en los dos modos) son cinco:

| hijo | ℓ | N rápido | N lento |
|---|---|---|---|
| 30, 0.096 | 30 | 1.49 | 0.0037 |
| 100, 0.184 | 10 | 0.79 | 0.0004 |
| 100, 0.184 | 30 | 1.96 | 0.0014 |
| 1 000, 0.898 | 10 | 0.61 | 0.0003 |
| 1 000, 0.898 | 30 | 2.34 | 0.0024 |

En todas, **el modo rápido está excitado a orden uno y el de χ, a ≤ 0.4 %.**

**c) El supuesto de §3 (que N se conserva después de la ventana) no se puede comprobar con esta maquinaria.**
- La energía congelada E del gauge de radiación uniforme deja de ser definida positiva ya a Δln a ≈ 0.016–0.025 después de φ = +0.2, en todos los hijos y en todos los n, incluso con k/aH ≈ 1 000.
- Es el mismo problema que en U56: con H ≠ 0, la energía congelada en este gauge no es una buena energía.

### 3. Lo que queda establecido, y lo que no

**Establecido:**
- el rebote con cuscutón (juego adoptado) **crea partículas a orden uno** (N ≈ 0.4–2.3) en un modo escalar, en todos los armónicos calculados hasta 10⁵;
- **el modo casi puro de δχ cruza casi en el vacío** donde el cálculo está limpio (N ≤ 0.4 %).

**No establecido:**
- **qué acaba en P_ζ.** Si el modo rápido, mitad δχ, se convierte después en el modo de sonido de la radiación y se diluye, el espectro del inflatón vendría casi solo del modo lento, y los rasgos serían ≲ ±13 % (N ≤ 0.004 en las filas limpias). Si, en cambio, su mitad δχ pasa al inflatón, habría rasgos de orden uno. La respuesta pide seguir los modos desde el rebote hasta que salen del horizonte en la inflación, en un gauge que no degenere (U62b).

### 4. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0–G2 pasan»:** acertado.
- **«REBOTE_EXCITA, con N ~ 1–6»:** acertado en el veredicto, con N de 0.4 a 2.3, pero en el modo equivocado. El modo de χ casi no se excita.
- **«Consecuencia: rasgos de orden uno, excluidos»:** **no se puede afirmar todavía,** por el §2.

### 5. Siguiente, como propuesta

- **U62b:**
  - la evolución de los dos modos desde la contracción, con el vacío puesto donde el modo sea adiabático, hasta la salida del horizonte en la inflación;
  - en el gauge plano o con variables invariantes de gauge;
  - comparando P_ζ con el de un control sin rebote, para medir los rasgos frente a ℓ.
  - Antes, identificar qué es el modo rápido de velocidad ≈ 1.9: radiación modificada por el cuscutón, o un modo de la ligadura.
- **En la v4,** lo que se puede decir ya: el rebote con cuscutón crea partículas a orden uno en un modo mixto, y casi no en el modo del inflatón. Si eso llega a P_ζ es una pregunta abierta.

---

## Nota de corrección · 2026-11-16 13:40 (Madrid): el «modo rápido» es un artefacto del marco

Al preparar U62b, y tras leer C3c del explorador (`whalefall-birth-explorer`, PR #11), comprobé que **el hamiltoniano congelado que usé no es único.**
- La acción reducida tiene un término q̇ᵀB_s q, con B_s la parte simétrica de B.
- Si se integra por partes (una transformación canónica dependiente del tiempo, p′ = p − a³B_s q), la dinámica es la misma, pero el vacío instantáneo cambia.

**En el marco sin B_s, con φ = +0.2 y n ≥ 10⁴,** las frecuencias congeladas son 0.578k y 1.000k: la radiación y χ, como pide la física local. El operador de C3c da lo mismo. En el marco de U62 eran 1.90k y 1.00k.

**Consecuencias:**
- **El «modo rápido» de §2 a), con velocidad ≈ 1.9, no es físico.**
- **Los N de §1 se midieron contra un vacío que no corresponde a los modos locales,** y no son un número de partículas físico. Esto vale también para el modo lento.
- **El veredicto REBOTE_EXCITA se queda como está escrito,** pero no tiene lectura física.
- **La afirmación del §3, «el rebote crea partículas a orden uno en un modo escalar», queda retirada.** La sustituye U62b, que pone los estados en el marco sin B_s y lejos de la ventana.
- **La indefinición de E justo después de la ventana (§2 c) también es del marco:** en el marco sin B_s, E′ sigue positiva durante varios e-folds para n ≳ 10³.

Los scripts de esta comprobación están en `results/u62b_runs/previo/`.

---

## U62b · 2026-11-16 14:05 (Madrid): los modos físicos, en el marco que recupera la física local

**Preregistro:** `docs/preregistros/ESPEC_U62b.md` (16 nov, 13:40). Las comprobaciones de frecuencias previas al ESPEC están declaradas allí.
**Código:** `chain/u62b_bogoliubov_fisico.py`. **Resultados:** `results/u62b_bogoliubov_fisico.json` y `results/u62b_runs/produce.out`. Tardó 19 min.

### 1. Veredicto preregistrado: INCONCLUSO (falla G2 en un punto)

| puerta | resultado |
|---|---|
| G0, la normalización | **pasa:** 9×10⁻¹⁶ |
| G1, la identidad de Bogoliubov | **pasa:** 2.5×10⁻⁷ |
| G2, la convergencia en rtol | **falla** en n = 10⁵ del primer hijo (A·H = 30). En n = 10³ y 10⁴ pasa, con 10⁻⁶ y 10⁻⁷ de la tolerancia |
| G3, las etiquetas (χ ≈ k, radiación ≈ k/√3) | **pasa** en todas las filas definidas |
| G4, la ligadura en los tramos nuevos (≤ 10⁻⁶) | **pasa:** 5.4×10⁻⁷ |

**Sobre G2:**
- **El punto que falla** es la única fila de n alto que ya salía mal por otro lado: N_rad ≈ 2×10⁶ con los extremos en Δln a = 2, frente a 10⁻⁸ con Δln a = 1.5. Por eso estaba como NO_DEFINIDA.
- **Es un fallo numérico,** no físico: unas 10⁴ oscilaciones en el tramo más largo, con un error que en el marco del vacío se convierte en «partículas».
- No interviene en ninguna fila de n_obs. Pero la puerta estaba escrita sobre ese n, y no se cambia.

**Con las puertas en verde, el veredicto habría sido INTERMEDIO:** de las 7 filas definidas de n_obs, ninguna tiene N_χ ≥ 0.1, pero tres pasan de 0.01.

| hijo (A·H, f) | ℓ = 2 | ℓ = 10 | ℓ = 30 |
|---|---|---|---|
| 30, 0.096 | no definida (n = 94) | **0.034** (469) | 0.0012 (1 408) |
| 100, 0.184 | no definida (350) | 0.0046 (1 749) | 1.2×10⁻⁴ (5 247) |
| 1 000, 0.542 | no definida (123) | no definida (617) | **0.052** (1 851) |
| 1 000, 0.835 | no definida (173) | no definida (863) | inestable (2 588) |
| 1 000, 0.898 | no definida (865) | **0.0115** (4 324) | 1.2×10⁻⁵ (12 972) |

Las filas no definidas lo son por falta de adiabaticidad en los extremos (|ω̇|/ω² > 0.1), y una por inestabilidad. Son 8 de 15, entre ellas todas las de ℓ = 2.

### 2. Lo que se ve (informativo)

**a) El marco importa, y el de U62 era el equivocado.**
- En el marco sin B_s, los extremos dan exactamente ω_χ/k = 1.000 y ω_rad/k = 0.577 en todas las filas definidas.
- Con n alto, la mezcla física desaparece: N_χ y N_rad ≲ 10⁻⁵ para n ≳ 10⁴. En U62, en el mismo n, N era de 0.4 a 2.

**b) La memoria del rebote sigue la escala n_mem de U61b, con una curva casi universal.** N_χ en las filas definidas, frente a n/n_mem:

| n/n_mem | N_χ (cinco hijos) | rasgo de potencia ≈ ±2√(N(1+N)) |
|---|---|---|
| ≈ 0.7–1.1 | 0.03–0.10 | ±35–65 % |
| ≈ 2.0–2.4 | 0.005–0.012 | ±14–22 % |
| ≈ 3.4 | 0.0012 | ±7 % |
| ≈ 4.6–4.8 | 3×10⁻⁵–3×10⁻⁴ | ±1–3.5 % |
| ≳ 12 | ≲ 10⁻⁵ | ≲ ±0.6 % |

- **El rebote con cuscutón borra su memoria** por encima de unas pocas n_mem, como con NPQT (U56: borrada para ν > 2.5·n_mem).
- **El criterio de U61/U61b (n_mem) era el correcto,** y la lectura de U61 §3 («la mezcla no decae con n, y n_mem subestima la huella») venía del artefacto del marco. Queda corregida aquí.
- **Para n ≲ n_mem los modos no son adiabáticos** en ningún punto a menos de Δln a = 2 del rebote. Allí el número de partículas no está definido (N de 0.3 a 16 con esta prescripción, y cambia con Δ). Solo un cálculo del espectro hasta la salida del horizonte en la inflación diría qué pasa en P_ζ.

**c) La ganancia de energía de U59b (≈ 14 sin decaer con n) también es del marco.** Se medía con la energía del marco de U62. El veredicto ESCALARES_SANOS de U59b (sin fantasmas, sin crecimiento con n) no cambia. La frase «una mezcla de orden uno que no decae con n» sí, y se corrige con una nota en `INFORME_U59.md`.

### 3. Lo que esto significa para los hijos visibles (lectura)

**En los cinco hijos visibles de U61b, la huella del rebote en el espectro primordial está en ℓ ≲ ℓ_mem:**
- **en ℓ ≈ ℓ_mem**, rasgos de decenas de por ciento: ±46 % en ℓ = 30 del hijo de ℓ_mem ≈ 36, ±37 % en ℓ = 10 del hijo de ℓ_mem ≈ 9;
- **a unas pocas ℓ_mem**, nada medible;
- **por debajo de ℓ_mem** (los ℓ = 2 de casi todos), modos que nunca son adiabáticos cerca del rebote, cuyo resultado depende del estado de la madre.

Es la misma estructura que con NPQT, con una escala de memoria cien veces mayor.
- **Lo que decide** si estos hijos son compatibles con el CMB es la potencia en ℓ ≲ ℓ_mem: ¿rasgos de ±40 %, o más, en los multipolos más bajos, donde la varianza cósmica ya es del 20–60 %?
- **La pregunta pide P_ζ propagado hasta la inflación (U62c)** y una proyección a C_ℓ.

### 4. Frente a lo que se esperaba (§6 del ESPEC)

- **«G0–G4 pasan»:** falló en G2, por un punto numérico.
- **«Entre 3 y 8 filas definidas»:** acertado, con 7.
- **«REBOTE_ADIABATICO»:** no. Habría sido INTERMEDIO: las filas de n ≈ n_mem tienen N_χ de 0.01 a 0.05.

### 5. Siguiente, como propuesta

- **U62c, el espectro.**
  - Los modos de ℓ = 2–40 de los hijos visibles, desde la contracción, con el vacío puesto donde el modo sea adiabático en la madre (o con dos prescripciones), hasta la salida del horizonte en la inflación.
  - El cambio al gauge de inflatón uniforme (C3d da la transformación, regular mientras χ̇ ≠ 0), y P_ζ(n)/P_ζ,control frente a ℓ.
  - G2 con un integrador distinto en los n más altos, o con un n máximo de 3×10⁴.
  - Cabe en la nube.
- **Corregir U62b como U62b′** (otro integrador o el n máximo de G2 en 3×10⁴) solo tendría sentido como repetición con datos nuevos (otros hijos o más n). No aporta información nueva frente a U62c.

**Añadido (informativo, después del veredicto), sobre el punto de G2:**
- Repetí la fila (A·H = 30, n = 10⁵) con rtol 10⁻¹¹ y 10⁻¹² (`results/u62b_runs/fila_n1e5.{py,json,out}`), y **empeora:** N_rad ≈ 10¹¹–10¹², y la identidad de Bogoliubov se rompe a 4–7×10⁻⁴.
- **No es la tolerancia del integrador,** sino un suelo de precisión del método: unas 10⁴ oscilaciones en tramos de 2 e-folds, con matrices de n = 10⁵ que llevan el ruido de cuadratura de U59 (≈ 10⁻¹⁰ relativo).
- En los otros cuatro hijos, n = 10⁵ sale bien (N ≲ 10⁻⁶), porque allí k/aH es menor.
- Para U62c: tope en n ≈ 3×10⁴, o tramos más cortos para n altos.

---

## U62c · 2026-11-16 18:15 (Madrid): el espectro primordial, propagado hasta la inflación

**Preregistro:** `docs/preregistros/ESPEC_U62c.md` (16 nov, 17:35).
**Código:** `chain/u62c_espectro_rebote.py`. **Resultados:** `results/u62c_espectro_rebote.json` y `results/u62c_runs/`. Tardó 23 min.

**Una corrección de código, declarada.** En el primer intento, el quinto hijo (A·H = 1 000, f = 0.898) hizo fallar la búsqueda de la referencia: vive en la hoja χ < 0, y la referencia buscaba un χ_× negativo en un fondo con χ > 0. Se corrigió usando |χ_×| (el potencial es par) y se repitió todo. El primer intento está en `results/u62c_runs/intento1_fallido/`, con una nota.

**Una prueba de funcionamiento, declarada:** una fila (n = 200 del primer hijo), mirando solo los wronskianos y el tiempo (`results/u62c_runs/previo/prueba_humo.py`).

### 1. Veredicto preregistrado: INCONCLUSO (falla G5)

| puerta | resultado |
|---|---|
| G0, la normalización | **pasa:** 9×10⁻¹⁶ |
| G1, el wronskiano al final (i/2 a ≤ 10⁻³) | **pasa:** 4×10⁻⁸. La acción ADM de la etapa 1 y la de Boruah, Kim y Geshnizjani de la etapa 2 tienen la misma normalización |
| G2, la tolerancia | **pasa:** ≤ 6×10⁻⁸ |
| G3, el control n_c = 10·n_mem (\|F_tot − 1\| ≤ 0.03) | **pasa:** ≤ 0.6 % en los cinco hijos |
| G4, el cambio de etapa (10⁻⁴ frente a 10⁻⁵) | **pasa:** ≤ 4×10⁻⁴ |
| G5, la ligadura del fondo (≤ 3×10⁻⁵) | **falla:** hasta 1.3 |

**Por qué falla G5:**
- El tramo de contracción, integrado hacia atrás hasta 0.5 e-folds del giro de la madre, pierde la ligadura de Friedmann a medida que se acerca al giro. Allí H ~ 10⁻⁸ y φ está a 10⁻⁸ de −φ∞.

| hijo | residuo hasta Δln a = 2 (A2) | hasta 3 | hasta 4 | hasta A1 |
|---|---|---|---|---|
| 30, 0.096 | 4×10⁻⁸ | 2×10⁻⁶ | 7×10⁻⁵ | 7×10⁻⁵ |
| 100, 0.184 | 8×10⁻⁸ | 3×10⁻⁶ | 2×10⁻⁴ | 1.5×10⁻³ |
| 1 000, 0.542 | 2×10⁻⁷ | 9×10⁻⁶ | 5×10⁻⁴ | 0.46 |
| 1 000, 0.835 | 4×10⁻⁷ | 2×10⁻⁵ | 1×10⁻³ | 0.93 |
| 1 000, 0.898 | 5×10⁻⁷ | 2×10⁻⁵ | 1×10⁻³ | 1.1 |

(`results/u62c_runs/ligadura_contraccion.out`, calculado tras el veredicto)

- **Mi comprobación previa solo miró el tramo de expansión.** No vi esto.
- **Consecuencia:** el estado A1 de los tres hijos de A·H = 1 000 está puesto sobre un fondo que no cumple la ligadura y no es fiable. El estado A2 sí lo es en todos.

**Con las puertas en verde, el veredicto habría sido RASGOS_FUERTES,** y lo sería también solo con A2: D_A2 = 1.97, 0.66, 16.6, 14.5 y 1.00. Todos son ≥ 0.5.

### 2. Lo que sale (A2, el estado fiable)

**F_tot = P_ζ/P_atractor,** con el hijo en el vacío a 2 e-folds del rebote en la contracción:

| hijo (A·H, f; ℓ_mem) | ℓ = 2–5 | ℓ ≈ 6–15 | ℓ ≈ 15–30 | ℓ = 40 |
|---|---|---|---|---|
| 30, 0.096 (8.9) | 0.67–2.97 | 0.78–1.38 | ±10–20 % | 1.00 |
| 100, 0.184 (4.2) | 0.77–1.66 | 0.85–1.20 | ≤ 1 % | 1.00 |
| 1 000, 0.542 (36) | **6.5–17.6** | 0.7–1.9 | 1.5–2.0 | 0.80 |
| 1 000, 0.835 (25) | **3.8–15.5** | 0.77–2.7 | 0.78–1.04 | 1.00 |
| 1 000, 0.898 (4.9; hoja χ < 0) | 1.4–2.0 | ±20 % | ≤ 2 % | 1.00 |

- **El hijo que empieza tras el rebote** (B) se queda en F_tot = 0.83–1.26, casi siempre a ≤ 5 %. La era de radiación previa a la inflación casi no deja huella. **Los rasgos vienen de la contracción y del rebote.**
- **Los rasgos son sobre todo un exceso de potencia en los ℓ más bajos,** de ×1.7–3 con A·H = 30–100 y de ×4–18 con A·H = 1 000, y oscilan con n. Encaja con la amplificación de los modos que salen del radio de Hubble en la contracción, el efecto conocido de los universos con rebote.
- **Se apagan hacia ≈ 2–4 ℓ_mem.** Es la escala de memoria de U61b y U62b.
- **Dependencia del estado:**
  - con A·H = 30 y 100, A1 y A2 coinciden a ~1 % (mediana de |P_A1/P_A2 − 1|: 1.4 % y 0.9 %), y el resultado no depende de dónde se ponga el vacío en la madre;
  - con A·H = 1 000 difieren mucho, pero allí A1 no es fiable (§1).

### 3. Lo que significa (lectura, sin proyectar a C_ℓ)

- **Los hijos visibles del cuscutón tienen un espectro primordial con excesos de orden uno o mayores en ℓ ≲ 5–30.** El cielo observado tiene, si acaso, poca potencia en esos ℓ: el cuadrupolo bajo.
  - **Con A·H = 1 000** (×4–18 en ℓ = 2–5), es casi seguro que quedan excluidos.
  - **Con A·H = 30–100** (×1.7–3 en algunos ℓ ≤ 8, con una varianza cósmica del 30–60 % por multipolo), habría tensión, pero hace falta proyectar a C_ℓ para decirlo.
- **El hijo de f = 0.898 está en la hoja χ < 0,** donde el portal σχ|H|² hace taquiónico al Higgs: según el supuesto de medida del preprint, no alojaría observadores. De los cinco visibles de U61b, cuatro están en la hoja χ > 0.
- **Si esto se confirma, la lectura de U61b cambia de signo.** Con el cuscutón, un hijo que nos aloje necesita el rebote escondido (ℓ_mem < 2, es decir, más e-folds). La curvatura sigue siendo indetectable, y la huella visible pasa a ser una exclusión, no una predicción.

### 4. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0–G5 pasan; la menos segura, G3»:** G3 pasó con holgura (≤ 0.6 %). Falló G5, que no había mirado bien.
- **«RASGOS_FUERTES o MODERADOS»:** RASGOS_FUERTES con las puertas en verde, y también solo con A2.
- **«Dependencia del estado importante»:** no con A·H = 30–100. Con A·H = 1 000 está por decidir.

### 5. Siguiente, como propuesta

- **U62d**, con datos nuevos:
  - el tramo de contracción integrado hacia delante desde el giro (como hace `u61b.host`), con la proyección de la ligadura;
  - A1 en el giro;
  - y una proyección a C_ℓ con Sachs–Wolfe en el S³ cerrado, para compararla con la potencia observada en ℓ = 2–30 de Planck.

  Diría si los hijos visibles del cuscutón quedan excluidos. Cabe en la nube.

---

## U62d · 2026-11-16 21:05 (Madrid): el espectro desde el giro de la madre, proyectado a C_ℓ

**Preregistro:** `docs/preregistros/ESPEC_U62d.md` (16 nov, 18:29).
**Código:** `chain/u62d_cl_rebote.py`, que usa la propagación de U62c sin cambios. **Resultados:** `results/u62d_cl_rebote.json` y `results/u62d_runs/`. Tardó 88 min.

**Una corrección de código, declarada.** En el primer intento, el control de los hijos de A·H = 1 000 dio F = 3 356 y 10 722.
- Las diferencias finitas de U62b (Ḃ_s y la adiabaticidad) tienen pasos relativos a |t|, y aquí el origen de tiempos estaba en el giro de la madre: cerca del rebote, |t| ≈ 4×10⁶, y el paso valía ≈ 40 unidades, del orden de la ventana entera.
- Se movió el origen al rebote, como en U62b y U62c, y se repitió todo. El primer intento está en `results/u62d_runs/intento1_fallido/`, con una nota.
- U62b y U62c no estaban afectadas: su origen ya estaba en el rebote.

**Una prueba de funcionamiento, declarada:** dos filas del primer hijo, mirando solo los wronskianos, dónde sale el modo y el tiempo (`results/u62d_runs/previo/prueba_humo.py`).

### 1. Veredicto preregistrado: INCONCLUSO (falla G6)

| puerta | resultado |
|---|---|
| G0, la normalización | **pasa:** 8×10⁻¹⁶ |
| G1, el wronskiano | **pasa:** 4.6×10⁻⁷ |
| G2, la tolerancia | **pasa:** ≤ 8×10⁻⁸ |
| G3, el control n_c = 10·n_mem | **pasa:** ≤ 0.4 % |
| G4, el cambio de etapa | **pasa:** ≤ 3×10⁻⁴ |
| G5, la ligadura del fondo | **pasa:** 8.8×10⁻⁸ (en U62c llegaba a 1.3) |
| G6, el muestreo (media rejilla) | **falla** con A1 en los dos hijos de A·H = 1 000: Δχ² pasa de 46.6 a 32.0 y de 31.7 a 55.0. Con A2 pasa en todos, y con A1 en los otros dos |

**Con las puertas en verde, el veredicto habría sido MIXTO.**

**Los cocientes R_ℓ = C_ℓ/C_ℓ(atractor)** (Sachs–Wolfe) **y la distancia a ΛCDM en ℓ = 2–30** (Δχ² de Wishart, f_sky = 0.8, con dilución w = 0.6 y sin ella, w = 1):

| hijo (A·H, f; ℓ_mem) | estado | R_2 | R_3 | R_5 | R_10 | R_20 | Δχ² (w = 0.6) | Δχ² (w = 1) |
|---|---|---|---|---|---|---|---|---|
| 30, 0.096 (8.9) | A1 | 1.82 | 1.44 | 1.11 | 1.07 | 1.03 | 0.7 | 1.7 |
| | A2 | 1.87 | 1.45 | 1.10 | 1.08 | 1.03 | 0.8 | 1.8 |
| 100, 0.184 (4.2) | A1 | 1.25 | 1.02 | 0.98 | 0.98 | 1.00 | 0.05 | 0.14 |
| | A2 | 1.22 | 1.02 | 0.98 | 0.99 | 1.00 | 0.05 | 0.12 |
| 1 000, 0.542 (36) | A1 | 18.2 | 7.0 | 3.2 | 1.5 | 1.6 | 46.6* | 90* |
| | A2 | 10.5 | 5.7 | 2.9 | 1.5 | 1.2 | 21.4 | 39 |
| 1 000, 0.835 (25) | A1 | 10.6 | 4.5 | 2.5 | 2.0 | 1.4 | 31.7* | 63* |
| | A2 | 6.0 | 3.3 | 2.0 | 1.04 | 1.1 | 7.9 | 15 |

\* Falla el muestreo (G6).

- **A·H = 30 y 100: compatibles.**
  - Δχ² ≤ 1.8 incluso sin dilución.
  - No dependen del estado: la mediana de |P_A1/P_A2 − 1| es del 2 % y el 1 %.
  - **La proyección suaviza los rasgos oscilantes de P_ζ** (×0.7–3 en n): en C_ℓ queda un exceso del 22–87 % en el cuadrupolo, del 2–45 % en el octupolo y ≤ 10 % a partir de ℓ = 5. La varianza cósmica allí es del 63 % y el 53 %.
- **A·H = 1 000: muy desfavorecidos, pero no excluidos con la regla preregistrada.**
  - C_2 sale ×6–18 y C_ℓ, ×1.5–3 hasta ℓ ≈ 10–20.
  - Δχ² va de 8 a 90 según el estado y la dilución. La regla exigía ≥ 25 con A1 y con A2, y A2 da 21 y 8.
  - Dependen mucho del estado (mediana del 44–47 %). En estos hijos, A1 tuvo que acercarse hasta 1.5 e-folds al rebote para que E′ fuera positiva.
- **La fracción del núcleo de C_2 que queda fuera de la rejilla** (n < n_lo) es ≤ 2×10⁻³.

### 2. Lo que significa (lectura)

- **La lectura de U62c, «los hijos visibles quedan excluidos», era demasiado fuerte.**
  - **Con A·H = 30–100** la huella del rebote es un exceso modesto en el cuadrupolo (+20–90 %), dentro de la varianza cósmica. **Es una predicción compatible**, pero apenas comprobable: una sola realización del cielo no distingue un +50 % en C_2.
  - **Con A·H = 1 000** (madres grandes), el exceso es enorme y está muy desfavorecido.
- **El cuadrupolo observado está por debajo de ΛCDM, no por encima.** Con el dato real, cualquier exceso en C_2 estaría algo más desfavorecido que aquí, donde el dato es el propio ΛCDM. No se ha calculado.
- **De los cinco hijos visibles de U61b:**
  - uno está en la hoja χ < 0, sin observadores;
  - dos (A·H = 1 000) están muy desfavorecidos;
  - dos (A·H = 30 y 100) son compatibles, con un cuadrupolo algo alto.
- **El balance para el cuscutón** es el mismo que para NPQT en cuanto a la curvatura: indetectable, y todos cerrados frente a DESI. La huella del rebote sobrevive solo como un exceso pequeño en ℓ = 2–3 en los hijos de madres de A·H ≲ 100.

### 3. Siguiente, como propuesta

- **Si se quiere un veredicto limpio en A·H = 1 000:** U62e, con la rejilla de A1 densificada en n ≲ n_mem y datos nuevos (otros f).
- **Para una comparación real:** el C_ℓ completo con CLASS (geometría cerrada, espectro externo, ISW) y la verosimilitud de los ℓ bajos de Planck. Es más trabajo, y en la nube cabe en lotes pequeños.

---

## U62e · 2026-11-16 22:25 (Madrid): el C_ℓ completo con CLASS, frente a Planck 2018

**Preregistro:** `docs/preregistros/ESPEC_U62e.md`.
**Código:** `chain/u62e_cl_class.py`. **Resultados:** `results/u62e_cl_class.json` y `results/u62e_runs/produce.out`. Tardó 8 s.
**El dato:** el espectro TT de Planck 2018 (Commander en ℓ = 2–29), congelado en `data/planck2018/` con su procedencia.

### 1. Veredicto preregistrado: MIXTO

| puerta | resultado |
|---|---|
| G0, la tabla externa frente a la ley de potencias analítica | **pasa:** 8×10⁻¹³ |
| G1, la precisión (q_linstep 0.01 frente a 0.02) | **pasa:** Δχ² cambia ≤ 0.2 |
| G2, la cobertura (un estado definido por hijo) | **pasa** |

**La distancia a Planck,** Δχ²_Planck = −2 ln L(hijo) + 2 ln L(ΛCDM), con Wishart en ℓ = 2–29 y f_sky = 0.8:

| hijo (A·H, f; ℓ_mem) | estado | ¿definido? | C_2/ref | C_3 | C_5 | C_10 | Δχ²_Planck | sus ℓ = 2–5 | Δχ²_ΛCDM |
|---|---|---|---|---|---|---|---|---|---|
| 30, 0.096 (8.9) | A1 | no (muestreo: 1.8 → 3.9) | 1.56 | 1.34 | 1.11 | 1.07 | 1.8 | 1.5 | 1.1 |
| | A2 | sí | 1.60 | 1.35 | 1.10 | 1.08 | **2.1** | 1.6 | 1.2 |
| 100, 0.184 (4.2) | A1 | sí | 1.19 | 1.03 | 0.98 | 0.99 | **0.7** | 0.7 | 0.1 |
| | A2 | sí | 1.16 | 1.02 | 0.98 | 0.99 | **0.7** | 0.6 | 0.1 |
| 1 000, 0.542 (36) | A1 | no | 13.4 | 5.9 | 2.9 | 1.4 | 43 | 19 | 67 |
| | A2 | sí | 7.6 | 4.7 | 2.7 | 1.5 | **22.6** | 15.4 | 32 |
| 1 000, 0.835 (25) | A1 | no | 8.0 | 3.9 | 2.3 | 1.8 | 35 | 13 | 51 |
| | A2 | sí | 4.6 | 2.8 | 1.9 | 1.0 | **10.1** | 8.5 | 11.5 |

- **A·H = 30 y 100: compatibles** con Planck (Δχ² = 2.1 y 0.7).
- **A·H = 1 000: desfavorecidos** (Δχ² = 22.6 y 10.1, sin parámetros libres). No llegan al umbral de exclusión preregistrado (≥ 25).
- **Ninguno es preferido.** El cuadrupolo de Planck es 0.22 veces el de ΛCDM, y estos hijos lo predicen más alto, no más bajo.

### 2. Frente a U62d (Sachs–Wolfe contra ΛCDM)

- **El ISW y el Doppler diluyen el exceso,** como se esperaba: C_2 ×1.60 frente a ×1.87 en A·H = 30, y ×7.6 frente a ×10.5 en A·H = 1 000.
- **La dilución w = 0.6 que suponía U62d resultó algo pesimista para el exceso:** el C_ℓ completo da Δχ²_ΛCDM = 1.2, 0.1, 32 y 11.5, frente a 0.8, 0.05, 21 y 7.9.
- **Con el dato real** (el cuadrupolo bajo de Planck) los Δχ² cambian poco, porque los ℓ = 2–5 pesan poco en la suma (2ℓ+1).

### 3. Lo que significa

- **La huella visible del rebote con cuscutón** es un exceso en el cuadrupolo y el octupolo.
  - **Con madres pequeñas** (A·H ≲ 100) es modesto (+16–60 % en C_2) y compatible con Planck.
  - **Con madres grandes** (A·H ≈ 1 000) es enorme (×5–8) y está desfavorecido a nivel de Δχ² ≈ 10–23.
- **El cielo no pide esta huella.** El cuadrupolo bajo observado va en la dirección contraria.
- **El cuscutón no ofrece una explicación de las anomalías de los ℓ bajos,** como tampoco NPQT. Lo que deja es una cota: los hijos que nos alojen tienen el rebote escondido, o vienen de madres pequeñas.
- **Para la v4,** el resultado honesto es que el nacimiento con cuscutón:
  - es compatible con el CMB (U58b, U62e);
  - no tiene fantasmas (U59b);
  - borra la memoria del rebote por encima de unas pocas n_mem (U62b);
  - y no predice una señal observable distinta de la curvatura positiva, que DESI DR2 ya pone a 2σ.

### 4. Lo que no cubre

- La verosimilitud oficial de Commander (clik).
- La polarización.
- El ajuste conjunto de parámetros.
- La geometría cerrada en CLASS: no fue viable con |Ω_K| ~ 10⁻⁵.
- El hijo de la hoja χ < 0.
