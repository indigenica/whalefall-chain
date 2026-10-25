# U19b · El colapso de RP³ medido en la red: más caliente que U(1), pero el rendimiento blando no sigue a las cuerdas

*17 de octubre de 2026 (Madrid). Nota del 19 de octubre: el control U(1) (§8) muestra que, frente al U(1) de la misma red, el colapso de RP³ no es «más caliente»: rinde lo mismo (Y ≈ 1.07 veces el de U(1)) y es algo más blando por cuanto. El título se conserva como registro.*
- *Preregistro:* `docs/preregistros/ESPEC_U19b.md`, enmiendas 1–3 anteriores a la producción y enmienda 4 (informativa) posterior a estas cuatro redes; enmienda 5 (informativa) antes de las redes densas de ℓ_G = 2.4. Desde el 18 de octubre: ocho redes, cuatro con ℓ_G = 3.2 y cuatro con ℓ_G = 2.4 (§7).
- *Código:* `network/u19b_produce.py`, `network/u19b_analyze.py`.
- *Resultado:* `network/results/u19b/U19B_RESULT.json`; redes en `network/results/u19b/runs/`.
- *Ejecución:* cuatro redes RP³ de 512³ en WSL, de 9.4 a 10.7 h cada una (16 hilos las poco densas y 34 las densas, sin diferencia de tiempo: el kernel está limitado por la memoria).

## 0. En una frase

**Veredicto preregistrado: COLAPSO_BLANDO, por un margen mínimo.** Y = 1.70 ± 0.55, y la condición era Y − 1σ > 1.14: sale 1.150.
- Con las hipótesis de U17d en F = F₀, la rama viva falla la puerta primaria de Lyα. Eso ya lo sabíamos desde U21: el benchmark F₀ quedaba fuera por la isocurvatura.
- **El error lo domina la diferencia entre protocolos, no la estadística:** las redes poco densas dan Y ≈ 2.6 y las densas Y ≈ 0.76.
- **Esa diferencia sale casi entera de la normalización:** el momento blando del colapso, X, es casi el mismo en las cuatro redes, aunque la densidad de cuerdas varía un factor 3. El rendimiento blando no escala con la energía de las cuerdas, como supone U17.

## 1. Las cuatro redes (RP³, ℓ_G = 3.2, ℓ₁ = 4.30)

| Red | ξ₁ | X_lat | **Y_lat** | K_lat | c_p | r_W | r_N | η_lat | f_E,duro | ¿válida? |
|---|---|---|---|---|---|---|---|---|---|---|
| poco densa, s91 | 0.73 | 8.20 | **2.81** | 1.85 | 3.84 | 1.82 | 2.59 | 2.76 | 0.44 | sí (4/4) |
| poco densa, s92 | 0.78 | 7.66 | **2.45** | 1.61 | 4.47 | 1.70 | 2.42 | 2.58 | 0.43 | sí (4/4) |
| densa, s91 | 2.48 | 7.25 | **0.73** | 1.63 | 4.67 | 1.80 | 2.24 | 1.15 | 0.53 | sí (4/4) |
| densa, s92 | 2.36 | 7.50 | **0.79** | 1.73 | 4.30 | 1.80 | 2.34 | 1.18 | 0.52 | sí (4/4) |

Las cuatro condiciones de validez se cumplen en todas:
- al final quedan menos del 1 % de las cuerdas;
- las paredes caen por debajo del 5 % de su máximo;
- hay 0 % de sitios lejos del vacío;
- el contador queda plano a < 3 %.

## 2. El veredicto

| | Valor |
|---|---|
| Y (media de las 4) | **1.70** |
| σ total | 0.55 = estadística 0.545 ⊕ κ 0.047 ⊕ anarmonicidad 0.008 |
| σ/Y | 0.32 (≤ 0.5: el resultado es determinado) |
| Y − σ | **1.150** > 1.14 → **COLAPSO_BLANDO** |
| Rama cómoda (B4a): Y + σ ≤ 0.10 | 2.24: **no** |
| Sesgo anarmónico (cota superior, informativo) | +0.11 |

**Lectura honesta:**
- **El veredicto pasa su umbral por 0.01.** Con la lógica preregistrada es COLAPSO_BLANDO, pero está en la frontera con COLAPSO_INTERMEDIO.
- **La σ es casi entera la diferencia sistemática entre protocolos.** Con 2 semillas por protocolo, las semillas coinciden bien: 2.81/2.45 y 0.73/0.79.
- **El preregistro promedia los dos protocolos** porque acotan el atractor por arriba y por abajo. Se mantiene.

## 3. Época y momento del colapso

**El colapso termina tarde y es caliente,** como anunciaba el piloto:
- **r_W ≈ 1.8** (las paredes caen al 10 % en a ≈ 1.8·a₁) y **r_N ≈ 2.4** (el número se estabiliza en a ≈ 2.4·a₁). Son iguales en los dos protocolos.
- **c_p ≈ 4.3**: el momento medio por cuanto, en unidades de m_G a la época r_W. En U(1), Hiramatsu et al. dan ≈ 3.1. *(Nota, 19 de octubre: el U(1) de nuestra propia red da c_p = 5.8, así que RP³ es algo más blando que U(1) en la misma red; ver §8.)*
- **La fracción de energía en cuantos duros** (k ≥ 10·m_G·a₁) es 0.43–0.53. Es una cota inferior: la red no tiene la jerarquía física.
- **η_lat ≈ 1.2–2.8 > 1.** Como dijo la enmienda 3, a ℓ₁ ≈ 4 la energía de las paredes es comparable a la de las cuerdas, así que η_lat no es una fracción.

## 4. El hallazgo: el rendimiento blando no sigue a las cuerdas

| | Poco densas | Densas | Cociente densa/poco densa |
|---|---|---|---|
| ξ₁ | 0.75 | 2.42 | 3.2 |
| ξℓ, promediado sobre el colapso | 3.8 | 9.9 | 2.6 |
| **X_lat** (momento blando por cuanto de desalineamiento) | **7.9** | **7.4** | **0.93** |
| K_lat (número del colapso) | 1.73 | 1.68 | 0.97 |
| Momento en cuantos duros, por m_G·a₁·N_mis | 14.8 | 24.3 | 1.65 |
| Y_lat | 2.63 | 0.76 | 0.29 |

- **Las redes densas tienen hasta tres veces más cuerdas durante el colapso, pero producen el mismo número y el mismo momento de cuantos blandos.** La energía de más va a cuantos duros, que no cuentan para Lyα.
- **La diferencia de Y entre protocolos sale casi entera de dividir por ξ₁.** No es un cambio físico en el colapso.
- **Una lectura natural:** los cuantos blandos los fija el colapso de las paredes, una por dominio de Hubble más o menos, y no la longitud de las cuerdas. Las cuerdas radian sobre todo cuantos duros, que es justo el «colapso duro» que U17c y U17d necesitaban.
- **No es todavía una conclusión:** se ha visto con una sola masa (ℓ_G = 3.2). Por eso la enmienda 4 fija, **antes de ver las redes de ℓ_G = 2.4**, cómo se leerá la prueba con la otra masa:
  - si X no cambia mientras ξ₁ℓ₁ sí, el rendimiento blando no depende de las cuerdas;
  - si X sigue a ξ₁ℓ₁, se mantiene el escalado de U17.

**Por qué importa.** U17d y U21 traducen Y al universo con X = (8π/27)·R_ρ·ξℓ·Y y ξℓ ≈ 110 (B2, ℓ ≈ 90).
- Si X escala con ξℓ, X_B2 ≈ 162·1.7 ≈ 275, y la rama viva necesita F < F₀ (lo que ya dice U21).
- Si el rendimiento blando lo fijan las paredes, el X del universo no crecería con ℓ. Quedaría cerca del X_lat ≈ 8, muy por debajo del umbral estricto (≈ 95), **incluso en F₀**. Sería muy favorable para el modelo, pero **no se puede afirmar sin la prueba de ℓ_G = 2.4.**
- *(Nota, 18 de octubre: las redes de ℓ_G = 2.4 no lo confirman. X₀ crece con ℓ₁, aunque solo a ≈ 2σ, y el X del universo sería ≈ 170–320; ver §7 y su corrección. El párrafo anterior, sobre X ≈ 8 «incluso en F₀», queda superado.)*

## 5. Qué cambia en la banda de F

Sin calcular nada nuevo, por acotación entre las variantes que ya se pasaron por CLASS (U21, U20):
- **El Y medido (1.70, con 1σ superior 2.24)** cae entre la variante P de U21 (Y efectivo 0.94) y la H (2.82). Por Lyα, la rama viva pasa con F hasta entre 0.5·F₀ (H) y 0.7·F₀ (P).
- **La época medida (r ≈ 1.8–2.4)** cae entre P (r = 1) y H (r = 3) de U20. Por isocurvatura, la banda llega a entre 0.35·F₀ (H) y 0.5·F₀ (P).
- **Banda conservadora: F ≲ 0.35·F₀** (≈ 3.7×10¹⁴ GeV), la de H. El límite exacto, entre 0.35 y 0.5·F₀, necesita pasar el espectro medido por CLASS (la U17e prevista en el preregistro).
- **La rama cómoda (B4a):** U21 le daba F ≤ 0.18·F₀ con Y = 0.94. Con el Y medido, su banda se estrecha, sin calcular.

## 6. Lo que U19b no resuelve todavía

- **La tendencia con ℓ (γ)** necesita las 4 redes de ℓ_G = 2.4; la primera está corriendo. El veredicto usa γ = 1, el lado conservador. Con esas mismas redes se lee el diagnóstico de la enmienda 4.
- **La calibración con U(1)** frente a Hiramatsu et al. (4 redes, informativas) va después.
- **La extrapolación de ℓ₁ ≈ 4 a ℓ ≈ 90** sigue siendo el salto grande. Si el rendimiento blando no sigue a las cuerdas, este salto se vuelve mucho más favorable.
- **La isocurvatura con la normalización C medida** y **el espectro medido pasado por CLASS** (U17e) siguen pendientes.

## 7. Las cuatro redes de ℓ_G = 2.4 (18 de octubre): el rendimiento blando sigue a ℓ, no a ξ

*Redes de WSL, RP³ 512³, ℓ₁ = 3.50; las cuatro válidas. Análisis con `network/u19b_analyze.py`, sin cambios. Las redes U(1) de control siguen corriendo.*

| ℓ_G = 2.4 | Poco densas (s91, s92) | Densas (s91, s92) | Cociente densa/poco densa |
|---|---|---|---|
| ξ₁ | 0.56, 0.57 | 3.16, 3.13 | 5.5 (en ξ₁ℓ₁) |
| **X_lat** | **6.54, 6.74** | **7.68, 7.55** | **1.15** |
| Y_lat | 3.58, 3.62 | 0.75, 0.74 | 0.21 |

**El veredicto no cambia: COLAPSO_BLANDO.** Se da con ℓ_G = 3.2. Con ℓ_G = 2.4, Y = 2.17 ± 0.83.
- La tendencia con ℓ es γ − 1 = −1.2 ± 2.4: indeterminada. El signo iría a favor del modelo, porque Y baja al subir ℓ, pero el error es el doble del valor.
- **El veredicto sigue usando γ = 1,** el lado conservador.

**Diagnóstico de la enmienda 4: lectura (c), mixta,** con los umbrales fijados antes de ver estas redes.

| | ρ_X = X(2.4)/X(3.2) | ρ_ξℓ | (a): \|ρ_X/ρ_ξℓ − 1\| ≤ 0.15 | (b): \|ρ_X − 1\| ≤ 0.15 y \|ρ_ξℓ − 1\| > 0.3 |
|---|---|---|---|---|
| Poco densas | 0.84 | 0.61 | no (ρ_X/ρ_ξℓ = 1.37) | no (0.16) |
| Densas | 1.03 | 1.06 | sí | no (ρ_ξℓ no cambia) |

**Descomposición afín de la enmienda 5 (informativa):**

| ℓ_G | ℓ₁ | X₀ (parte que no depende de las cuerdas) | y_s (por unidad de energía de cuerdas) |
|---|---|---|---|
| 3.2 | 4.30 | 8.20 ± 0.30 | −0.085 ± 0.042 |
| 2.4 | 3.50 | 6.43 ± 0.10 | +0.116 ± 0.013 |

**Lo que dicen las ocho redes juntas:**
1. **El momento blando casi no depende de cuántas cuerdas hay.**
   - X va de 6.5 a 8.2 en las ocho redes, mientras ξ₁ℓ₁ cambia un factor 5.5.
   - Con ℓ_G = 2.4, cinco veces más cuerdas dan un 15 % más de X.
   - El término de cuerdas, y_s, es pequeño, y su signo cambia entre las dos masas.
2. **Pero sí depende de ℓ, y descarta la lectura más favorable.**
   - X₀ sube de 6.43 a 8.20 (+28 ± 5 %) cuando ℓ₁ pasa de 3.50 a 4.30 (+23 %). Es decir, X₀ ∝ ℓ₁^(1.2 ± 0.2).
   - **Con eso, «X₀ constante» queda excluida a unas 6σ,** aunque estadísticas de solo dos semillas por protocolo.
   - Una lectura natural: la tensión de una cuerda global crece como ln(m_r/H) ≈ ℓ. El colapso de las paredes libera una energía por dominio proporcional a esa tensión, con independencia de cuántas cuerdas haya.
3. **El X del universo** (B2, ℓ ≈ 90) con las tres extrapolaciones fijadas en la enmienda 5, cada una más el término de cuerdas conservador (162·(y_s + 2σ) = 23):

| Extrapolación de X₀ | X_univ | ¿Umbral estricto de δA (X ≤ 95)? | ¿Puerta primaria (X ≤ 185)? | ¿La apoyan los datos? |
|---|---|---|---|---|
| Constante | 31 | sí | sí | **no** (excluida a ~6σ) |
| ∝ ℓ | 195 | no | no, por poco | **sí** (exponente 1.2 ± 0.2) |
| Lineal en ℓ₁ | 221 | no | no | compatible |

Sin el término de cuerdas, la extrapolación ∝ ℓ da X ≈ 172, dentro de la puerta primaria. El escalado de U17 (filas Y de U17e) daba 275.

**Qué significa para el modelo:**
- **La vuelta de F₀ y de LISA (lectura (b)) no se confirma.** El rendimiento blando no se queda en el valor de la red; crece con ℓ.
- **La extrapolación es más favorable que la de U17,** en un factor 1.4–1.6 en X, porque no crece con ξ.
- **La banda de F probablemente sube algo por encima de 0.5·F₀,** pero F₀ queda en la frontera de la puerta primaria de Lyα, no dentro con margen.
- Cuantificarlo necesita pasar X_univ ≈ 172–195 por CLASS: **una U17f, que hay que preregistrar antes** (enmienda 5).

**Lo que queda:**
- las redes U(1) de control, que están corriendo;
- U17f, si el autor lo decide;
- la palanca en ℓ es pequeña (de 3.5 a 4.3, frente a 90). Un tercer ℓ_G, por ejemplo 4.0, a 512³ en WSL, mediría el exponente de X₀ con mucha más palanca.

### Corrección de §7 · 18 de octubre de 2026, 20:53 (Madrid), tras una auditoría independiente

Una auditoría reimplementó el análisis desde los ficheros de las redes, sin el código de `u19b_analyze.py`. **Reproduce todos los números de §7**, pero muestra que la lectura de §7 se pasaba de firme en varios puntos. Lo corregido:

1. **«X₀ constante excluida a ~6σ» es falso como afirmación estadística.**
   - Cada X₀ sale de 4 redes con 2 grados de libertad. Con una t de Student, t ≈ 5.5–6 da p ≈ 0.03, **unos 2σ**.
   - Además, el error no incluye las sistemáticas:
     - la incertidumbre de κ;
     - la resta de los cuantos anteriores a ℓ₁;
     - cuándo se toma el «desalineamiento congelado». Esta elección sola mueve X₀(3.2) de 8.20 a 7.96, y el exponente de 1.18 a 1.07.
   - Que y_s cambie de signo entre las dos masas (−0.085 ± 0.042 frente a +0.116 ± 0.013) indica que esos errores son optimistas, o que el modelo afín no es el mismo en las dos.
   - **Lectura correcta:** los datos prefieren que X₀ crezca con ℓ, a ≈ 2σ. No excluyen que sea constante.
2. **Dos masas no distinguen «∝ ℓ» de «lineal en ℓ₁», y la tabla las ordenaba al revés.**
   - La lineal pasa exactamente por los dos puntos. ∝ ℓ queda a 0.9σ.
   - La ley de potencias ajustada (exponente 1.18 ± 0.2) da X₀(90) ≈ 300, y ≈ 320 con el término de cuerdas: **peor que el escalado de U17 (275).** Su banda de 1σ va de ≈ 165 a ≈ 540.
   - **Lectura correcta:** con las extrapolaciones que crecen con ℓ, el X del universo cae entre ≈ 170 y ≈ 320. Eso va de un factor ≈ 1.6 mejor a ≈ 1.2 peor que U17. Solo «X₀ constante» (31), que los datos no favorecen, sería mucho mejor.
   - **Se retira la frase «la banda de F probablemente sube algo por encima de 0.5·F₀».** La banda sigue siendo la de U17e, U20b y U23, F ≈ 0.09–0.5·F₀. Una U17f tendría que cubrir X entre ≈ 170 y ≈ 320.
3. **La dependencia de ξ era menor de lo que decía §7.**
   - «Cinco veces más cuerdas» usa ξ₁ℓ₁. Con el ξℓ promediado sobre el colapso, que es lo que define la enmienda 4 para D1, el cociente es 2.9–3.7 en ℓ_G = 2.4.
   - y_s(2.4) = +0.116 ± 0.013 no es cero: las cuerdas sí aportan algo de rendimiento blando.
   - El título de §7, «sigue a ℓ, no a ξ», debería decir **«depende poco de ξ y más de ℓ».**
   - (El código calcula D1 con ξ₁ℓ₁ y no con ξ̄ℓ. D1 es informativo y no entra en la lectura. La ventana «ℓ₁ + ln r_W» del ESPEC debería ser ℓ₁ + 2·ln r_W, porque ℓ = 2 ln a + const.)
4. **«La tensión de la cuerda» no es la única explicación de que X₀ crezca con ℓ₁.** Todo lo siguiente cambia a la vez que ℓ₁ entre las dos masas, y no se puede separar con solo dos:
   - **la época del colapso:** X lleva un factor r por construcción, y r_N es 2.40 frente a 2.02. Con t_col ∝ μ/σ ∝ ℓ, esto da por sí solo X ∝ ℓ;
   - **el tamaño de la caja:** ≈ 1.5 volúmenes de Hubble al final con ℓ_G = 3.2;
   - **el núcleo congelado,** activo solo con ℓ_G = 3.2;
   - **la jerarquía m_r/m_G,** que cambia el corte de k_cut;
   - **la resta de los cuantos anteriores a ℓ₁:** quita el 44 % frente al 29 %;
   - **los protocolos,** que no son idénticos entre masas: el ℓ_i de las densas es 3.0 frente a 3.5, y el fin del drenaje de las poco densas también cambia.
5. **El signo de γ − 1 sale solo de las redes poco densas** (−1.5 frente a +0.1 en las densas). Es el artefacto de normalizar por 1/ξ₁ que identifica §4, así que no apunta a favor del modelo.
6. **El término de cuerdas de 23** (162·(y_s + 2σ)) sale solo de ℓ_G = 2.4. Con ℓ_G = 3.2 sería 0.
7. **El veredicto de §2 es frágil, y no solo «por un margen mínimo».** Por la letra del preregistro sigue siendo **COLAPSO_BLANDO** (Y − σ = 1.150 frente a 1.14). Pero dos elecciones razonables lo pasan a **COLAPSO_INTERMEDIO**:
   - tomar el desalineamiento aún congelado: Y − σ = 1.139;
   - restar el sesgo de anarmonicidad que el propio análisis informa (+0.11): Y − σ ≈ 1.04.

**Para medir de verdad la dependencia con ℓ** hace falta una tercera masa, con más palanca (por ejemplo ℓ_G = 4.0), y protocolos idénticos entre masas. Su lógica se preregistraría antes.

## 8. El control U(1) (19 de octubre): el colapso de RP³ es como el de U(1) en nuestra red

*Cuatro redes U(1) de 512³ con ℓ_G = 3.2 y los mismos protocolos y semillas que las primarias (WSL); las cuatro válidas (4/4). Con ellas, U19b termina: 12 de 12 redes. Es el control «U(1) frente a Hiramatsu et al., informativo» del preregistro: **no cambia el veredicto.***

| ℓ_G = 3.2, media de 2 semillas | Poco densas: RP³ | Poco densas: U(1) | RP³/U(1) | Densas: RP³ | Densas: U(1) | RP³/U(1) |
|---|---|---|---|---|---|---|
| ξ₁ | 0.76 | 0.50 | 1.52 | 2.42 | 1.71 | 1.41 |
| X_lat | 7.93 | 4.85 | 1.63 | 7.38 | 4.94 | 1.49 |
| **Y_lat** | **2.63** | **2.44** | **1.08** | **0.76** | **0.72** | **1.06** |
| c_p | 4.15 | 5.18 | 0.80 | 4.49 | 6.45 | 0.70 |
| r_W | 1.76 | 1.59 | 1.11 | 1.81 | 1.65 | 1.09 |

**Con las cuatro redes de cada modelo:** Y(U(1)) = **1.58 ± 0.50** y Y(RP³) = 1.70 ± 0.55. El cociente es 1.07 ± 0.49 con los errores completos, aunque en cada protocolo, con las mismas semillas, sale 1.06 y 1.08.

**Lo que dice:**
1. **En nuestra red, el colapso de RP³ rinde lo mismo que el de U(1) por unidad de energía de cuerdas.**
   - RP³ da más momento blando X (×1.5–1.6) porque tiene más cuerdas (ξ₁ ×1.4–1.5, el R_ξ de siempre), y el cociente de Y es ≈ 1.
   - Además, RP³ es algo **más blando** por cuanto: c_p es un 20–30 % menor y r_W un 10 % más tarde.
2. **El U(1) de la red no es el U(1) supuesto en U17d.** Aquel colapso «tipo U(1)» (η = 1, c = 3, r = 1) daba Y = 0.94, y nuestra red da 1.58 (×1.7). Su c_p = 5.8 duplica el k̄/m ≈ 3.1 de Hiramatsu et al., con dos diferencias:
   - nuestra red tiene ℓ₁ ≈ 4.3 y m_r/m_G ≤ 25;
   - nuestra masa es constante, y la de ellos es la del axión QCD, que depende de la temperatura.

**Qué cambia en la lectura (sin cambiar el veredicto preregistrado, COLAPSO_BLANDO):**
- **Se corrige lo de «más caliente que U(1)»** (título, §3 y el contrato). Era una comparación con el U(1) analítico (Y = 0.94) y con Hiramatsu, no con el U(1) de nuestra propia red, que es lo que se compara igual con igual. Frente a ese U(1), **RP³ no es más caliente: Y es el mismo y c_p, algo menor.**
- **Dos lecturas posibles del valor absoluto:**
  - **(i) La red es fiable:** a ℓ ≈ 4 el colapso es igual de blando en U(1) que en RP³, y el Y ≈ 0.94 «tipo U(1)» de U17d no lo reproduce ningún modelo. Valen las filas Y de U17e.
  - **(ii) La red infla Y por igual en los dos modelos** (ℓ pequeño, poca jerarquía): lo robusto es el cociente RP³/U(1) ≈ 1.07. Calibrando con un U(1) de Y ≈ 0.94, RP³ quedaría en Y ≈ 1.0 ± 0.5, la variante P de U21. Esa variante pasa Lyα hasta 0.7·F₀ (por forma).
- **Ninguna de las dos lecturas estaba preregistrada como veredicto.** Elegir entre ellas, o pasar la lectura (ii) por CLASS, exige una **U17f preregistrada antes.** Mientras tanto, la mejor estimación sigue siendo la de U17e/U20b/U23: la banda combinada preregistrada, 0.09–0.7·F₀.
- **Para el artículo de la red,** es un resultado limpio y útil: con la misma red y el mismo contador, el colapso cuerdas–paredes de RP³ y el de U(1) dan el mismo rendimiento blando por unidad de energía de cuerdas. RP³ produce más cuantos solo porque tiene más cuerdas.
