# INFORME_U68 · El colapso de RP³ con un triplete no degenerado (U68, U68b)

*2026-11-21 23:06 (Madrid). Interno.*
- **Preregistros:**
  - `docs/preregistros/ESPEC_U68.md` (21:43);
  - `ESPEC_U68b.md` (21:47), con su enmienda 1, anterior a cualquier red.
- **Origen:** `docs/informes/NOTA_UV_MG.md` (PR #102), §8.
- **Código:** `network/kern.c` (`set_tilt_weights`), `network/u68_produce.py` y `network/u68_analyze.py`.
- **Resultados:** `network/results/u68/`, con 9 redes RP³ de 192³ en la nube a 4 hilos (7–10 min cada una).

> **Nota del 2026-11-21 a las 23:10 (Madrid).** Los números de este informe (§0–§4) salen de una primera producción con un error del contador: la velocidad se multiplicaba dos veces por la normalización z_c. Lo encontró la revisión de Codex. Esa producción está en `network/results/u68/_superseded/produccion1_z_doble/`. Las 9 redes se están repitiendo con el contador corregido (`ESPEC_U68b.md`, enmienda 2), y el veredicto de U68b será el de la repetición. El texto de abajo se conserva como registro.

## 0. En una frase

**IGUAL.** Con masas (1 : √2 : 2), la velocidad media de los cuantos blandos queda un 13 ± 1 % por debajo de la del triplete degenerado, dentro del ±15 % preregistrado.
- **La materia oscura que deja el colapso no se calienta** si la ruptura no es O(3).
- **Lo que sí cambia es el reparto:** el colapso se canaliza hacia el modo ligero, más de lo que darían especies independientes.

## 1. U68: INCONCLUSO (falla G1, antes de calcular ninguna red)

| puerta | resultado |
|---|---|
| G0, la identidad del núcleo con los pesos por defecto | **pasa:** campos idénticos bit a bit a los de `main`, en RP³ y U(1), con y sin núcleo congelado tardío |
| G1, las masas | **falla**, por dos motivos |
| G2, el contador por modo | **pasa:** ⟨k⟩ a 10⁻⁵, sin fugas a otros modos |

Los dos motivos de G1:
- **Un fallo de código de la comprobación.** La medida de WKB caía en x = τ/τ_A ≈ 2.1–2.3, fuera de la tabla de masas efectivas, que acababa en x = 2.0.
- **Un criterio demasiado estricto.** La fórmula de primer orden m_c² = m_G²(w₀ − w_c) se desvía −2.8 % para los modos de m² = 4m_G², por el desplazamiento del vacío radial.

## 2. U68b: IGUAL

**Lo que cambió** (`ESPEC_U68b.md`, declarado antes de cualquier red):
- G1 se compara con la masa exacta de la hessiana, con la tabla extendida hasta x = 2.6.
- La prueba de WKB arranca en el vacío desplazado (enmienda 1). La primera ejecución de U68b, que arrancaba en Q = nnᵀ, dio +1.1–1.5 % en los modos pesados y está guardada.
- La normalización canónica va por modo, z_c = Λ − μ_c.

**Las puertas:**

| puerta | resultado |
|---|---|
| G0 | pasa (la de U68) |
| G1 | pasa: Minkowski frente a la hessiana ≤ 0.4 %, WKB ≤ 0.02 % |
| G2 | pasa |
| G3, la validez de las 9 redes | pasa: cuerdas < 1 %, paredes < 5 %, lejos del vacío < 1 %, meseta del número a ≤ 2.7 % |

**El veredicto.** v̄ es la velocidad media de los cuantos blandos (k < 10 m_c a₁), en unidades de a₁.

| semilla | v̄, O3 | v̄, GEN | D = v̄_GEN/v̄_O3 − 1 |
|---|---|---|---|
| 91 | 2.88 | 2.43 | −0.158 |
| 92 | 2.78 | 2.47 | −0.113 |
| 93 | 3.07 | 2.68 | −0.127 |
| **media** | | | **−0.133 ± 0.013** |

Como |D| ≤ 0.15, **IGUAL**. El margen no es grande: D está a 1.3σ del umbral.

## 3. Lo que se ve: el colapso va al modo ligero

| | O3 (1, 1, 1) | GEN (1, 2, 4) | (IND) | O2 (1, 1, 4) | (IND) |
|---|---|---|---|---|---|
| fracción de cuantos blandos en los modos de masa m_G | 1 | **0.58 ± 0.04** | 0.39 | **0.88 ± 0.01** | 0.74 |
| densidad Σ m_c N_c frente a O3 | 1 | **1.07 ± 0.02** | 1.20 | **1.018 ± 0.004** | 1.14 |
| v̄ de cada modo (ligero → pesado) | 2.8 / 3.1 / 2.8 | 2.9 / 2.1 / 1.5 | — | 2.7 / 3.0 / 1.5 | — |
| número de cada modo frente al de un modo de O3 | 1.09 / 0.99 / 0.92 | 1.51 / 0.71 / 0.34 | 1 / 0.84 / 0.71 | 1.28 / 1.12 / 0.32 | 1 / 1 / 0.71 |
| r_W (paredes al 10 %) | 1.56 | 1.42 | — | 1.48 | — |

**Lectura.**
- **(IND), especies independientes, no describe el colapso.** El modo ligero se queda 1.5 veces el número de un modo de O3, y el más pesado, menos de la mitad de lo que daría (IND).
- **(CANAL), la hipótesis declarada en el ESPEC, sí lo describe:** las cuerdas giran su plano de giro hacia la dirección ligera, de pared más barata, y el colapso va sobre todo a ese modo.
- **El modo ligero conserva la velocidad de O3** (2.9 frente a 2.8–3.1). Los modos pesados salen más lentos (1.5–2.1): reciben menos energía del colapso y una parte mayor de sus cuantos es de desalineamiento.
- **Por eso el conjunto sale algo más frío (D < 0)** y con menos densidad de la que darían especies independientes.
- **El número de cada modo de O3** varía en ±8 % entre los modos, que es el ruido estadístico de esta lectura.

## 4. Qué cambia para la banda y para la ruta 6 de la nota ultravioleta

**Para la banda,** en orden de magnitud y sin CLASS:
- **La temperatura:** la materia oscura sale un 13 % más fría, en el sentido favorable para Lyα.
- **La abundancia:** un 7 % mayor a F fijo. Como ρ ∝ m^{1/2}F², la m_G de la banda bajaría un ≈ 12 %.

Los dos efectos son pequeños frente al factor 2.3 que movería la banda.

**Para la ruta 6 de `NOTA_UV_MG.md`** (la calidad imperfecta, que pierde la O(3) exacta):
- **En el retículo, perder la degeneración del triplete no estropea el colapso,** en una ruptura cuadrática general con un cociente de masas de 2.
- **Lo que deja es una predicción distinta:** la materia oscura iría sobre todo en el Goldstone más ligero, con los otros dos más escasos y más fríos.
- Eso cambiaría las nubes de superradiancia. Ya no habría tres especies degeneradas y una carga O(3) conservada, sino especies de masas distintas.

## 5. Lo que no prueba

- **El alcance del retículo:** 192³, ℓ₁ ≈ 3.5, y el cociente de masas solo llega a 2. Con jerarquías mayores o a ℓ ≈ 90, la canalización podría ser más fuerte. Su efecto sobre v̄ también podría cambiar.
- **Solo rupturas cuadráticas.** El caso con simetría residual discreta y potencial de grado alto, con varios vacíos, no se prueba (ESPEC_U68, «Una precisión sobre la nota»).
- **No pasa nada por CLASS.** La lectura para la banda es de orden de magnitud.

## 6. El resultado vigente: la producción repetida con el contador corregido (2026-11-22 a las 00:23 (Madrid))

Las 9 redes se repitieron con el contador corregido (`ESPEC_U68b.md`, enmienda 2: z_c solo multiplica el campo, no la velocidad), con las mismas semillas, protocolos y medidas. **Este es el veredicto de U68b.** Los §0–§4 de arriba son de la producción superada y se conservan como registro.

**Las puertas:** G0–G2 pasan (no dependían del error) y G3 pasa en las 9 redes, con la meseta del número a ≤ 1.6 %.

| semilla | v̄, O3 | v̄, GEN | D |
|---|---|---|---|
| 91 | 2.88 | 2.44 | −0.154 |
| 92 | 2.78 | 2.48 | −0.109 |
| 93 | 3.07 | 2.70 | −0.121 |
| **media** | | | **−0.128 ± 0.013** |

**IGUAL** (|D| ≤ 0.15). D está a 1.7σ del umbral.

**Las lecturas informativas, corregidas:**

| | O3 (1, 1, 1) | GEN (1, 2, 4) | (IND) | O2 (1, 1, 4) | (IND) |
|---|---|---|---|---|---|
| fracción de cuantos blandos en los modos de masa m_G | 1 | **0.59 ± 0.04** | 0.39 | **0.88 ± 0.01** | 0.74 |
| densidad Σ m_c N_c frente a O3 | 1 | **1.04 ± 0.02** | 1.20 | **1.006 ± 0.004** | 1.14 |
| v̄ de cada modo (ligero → pesado) | 2.8 / 3.1 / 2.8 | 2.9 / 2.2 / 1.5 | — | 2.7 / 3.0 / 1.5 | — |
| número de cada modo frente al de un modo de O3 | 1.09 / 0.99 / 0.92 | 1.52 / 0.69 / 0.32 | 1 / 0.84 / 0.71 | 1.29 / 1.13 / 0.30 | 1 / 1 / 0.71 |

**Lo que cambia respecto de la producción superada:** casi nada. D pasa de −0.133 a −0.128, la fracción ligera de 0.58 a 0.59 y la densidad de 1.07 a 1.04.

**Las conclusiones de §3 y §4 se mantienen,** con estos números:
- el colapso se canaliza al modo ligero;
- la materia oscura blanda sale un 13 % más fría;
- la abundancia sube un 4 % a F fijo, lo que bajaría la m_G de la banda un ≈ 8 % (antes, ≈ 12 %).
