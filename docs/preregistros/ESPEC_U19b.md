# U19b · η y la época del colapso cuerdas–paredes de RP³, medidos en la red · Preregistro

*15 de octubre de 2026, ≈09:10 (Madrid). Escrito antes de tocar el kernel y antes de producir ninguna red de U19b (la hora que cuenta es la de su commit). La validación y un piloto pequeño se hacen en la nube; la producción, en WSL.*

## Motivo

U17d (`docs/informes/INFORME_U17d.md`) dejó el sector oscuro pendiente de dos números del colapso de RP³ que nadie ha medido:
- **η:** la fracción de la energía de la red en H ≈ m_G/3 que acaba en cuantos de momento ~m_G;
- **la época:** r = a_decay/a_osc.

Lyα solo depende de la combinación **Y ≡ η·(c/c_E)·r**, con la que la variable de U17c/U17d es X = (8π/27)·R_ρ·ξℓ·Y. Para la rama viva (B2, hipótesis primarias de U17d) se cumple X = 162·Y, y:
- **pasa también el umbral estricto de δA** (X ≲ 95) si Y ≲ 0.59;
- **pasa solo la puerta primaria** (X ≲ 185) si 0.59 < Y ≲ 1.14;
- **falla** si Y > 1.14.

Un colapso tipo U(1) (η = 1, c = 3, r = 1) da Y = 0.94.

U19 no pudo medirlo: al final de la caja (m_G/H ≈ 10) quedaban cuerdas y un 12–31 % de sitios lejos del vacío, así que el contador de partículas no era válido.

## Qué cambia respecto a U19

**1. Fase tardía con el núcleo congelado** (cambio del kernel, apagado por defecto y bit a bit igual cuando lo está). Hasta τ_A, el momento en que m_r·a·dx = 1, la red es la de siempre: cuerdas físicas, y ℓ_A = ln(N/2). Desde τ_A:
- **el núcleo se congela:** la parte simétrica del potencial usa a_p = s·τ_A, con anchura comóvil constante, así que las cuerdas que queden siguen resueltas;
- **la inclinación sigue siendo física:** el término −ε se multiplica por a², no por a_p², y m_G·a sigue creciendo.

Así la caja llega a **m_G/H ≈ 95** sin perder resolución. El coste es que, en la fase tardía, la tensión de las cuerdas que queden deja de crecer (ℓ_eff = ℓ_A). Se declara como sistemática: con ℓ_G = 2.4 el colapso termina casi entero antes de τ_A, y la comparación con ℓ_G = 3.2 la acota.

**2. Contador espectral.** El número lineal de U19 (enmienda 1) se reparte por bins de k comóvil, y se guarda dN/d ln k en cada instantánea con ℓ ≥ ℓ₁ − 0.2.

**3. Observable de Lyα.** Por red, al final:
- X_lat = Σ_{k < k_cut} k·n_k / (m_G·a₁·N_mis) − κ, con k_cut = 10·m_G·a₁;
- **Y_lat = X_lat / [(8π/27)·ξ₁ℓ₁]**, la estimación en la red de Y = η·(c/c_E)·r.

Los ingredientes:
- **N_mis** = R_ctrl·⟨V_ε − V_min⟩·a₁³/m_G es el número del desalineamiento, con el control homogéneo de su arranque, como U19.
- **κ** es el momento que aporta el propio desalineamiento de un campo aleatorio, medido con un control sin defectos (V5). Su incertidumbre es el 100 % de κ.
- **Los cuantos con k > k_cut** se informan aparte: su fracción de energía f_E,duro y su fracción de número.

## Parámetros y producción (`network/u19b_produce.py`)

- **Masas:** ℓ_G ∈ {**3.2** (primaria), 2.4}. Es decir, m_G = m_r·e^(−ℓ_G), con m_r/m_G = 24.5 y 11.0, y ℓ₁ = ℓ_G + ln 3 = 4.30 y 3.50.
- **Caja:** N = 512, τ_A = 256, dτ = 0.25.
  - Fin: τ_B = 3τ_A con ℓ_G = 3.2 y 2τ_A con ℓ_G = 2.4. Las dos llegan a m_G/H ≈ 94.
  - Al final, m_r·a_p/(m_G·a) es 8.2 y 5.5.
- **Protocolos de U19:** poco denso (ξ_i = 0.5 en ℓ_i = 1, drenaje hasta ℓ_d = min(3.8, ℓ₁ − 0.3)) y denso (ξ_i = 6 en ℓ_i = 3.5; con ℓ_G = 2.4, en ℓ_i = 3.0, porque ℓ₁ = 3.5).
- **Semillas** 91 y 92.
- **Instantáneas:** 48, geométricas en τ desde τ₀·1.08 hasta τ_B. El espectro de `measure` se toma desde ℓ = ℓ₁ − 0.5 y el contador lineal desde ℓ₁ − 0.2.
- **Redes, por orden de prioridad:**
  1. RP³, ℓ_G = 3.2: 2 protocolos × 2 semillas (**primaria**);
  2. RP³, ℓ_G = 2.4: ídem;
  3. U(1), ℓ_G = 3.2: ídem (calibración con Hiramatsu et al., informativa).

  **12 redes.** Coste estimado en WSL: 2–5 h cada RP³ y ~1 h cada U(1), según los núcleos que dejen libres las cadenas de U15b/U15c.

## Validación (en la nube, antes de cualquier red de producción; `u19b_produce.py validate`)

- **V0:** con la congelación tardía apagada, `make test` pasa, incluida la reproducción bit a bit de una red de U2.
- **V1 (masa)** para cada ℓ_G: la de U19 (≤ 1 %).
- **V2 (contador frente a la EDO homogénea)** para cada ℓ_G y cada arranque: la de U19 (≤ 5 %), evaluada al final de la fase tardía.
- **V3 (fase tardía):** un campo sin defectos (θ₀ = 0.3, con ruido) evolucionado de 0.8τ_A a τ_B con la congelación activada.
  - El número lineal debe conservarse a ≤ 1 % entre τ_A y τ_B, porque es un invariante adiabático.
  - m_G efectivo debe seguir siendo físico: la frecuencia del modo homogéneo en τ_B, a ≤ 1 %.
- **V4 (contador espectral):** un campo sin defectos con dos paquetes de ondas de k conocido. El ⟨k⟩ por número debe recuperarse a ≤ 2 %.
- **V5 (κ):** un campo U(1) y otro RP³ con ángulos aleatorios pequeños (θ_rms = 0.3, sin cuerdas), con la longitud de correlación de la red poco densa en ℓ₁. κ = Σk·n_k/(m_G·a₁·N) al final. Se informa, y se resta.

Si falla alguna de V0–V4, U19b es **INVALIDO** y no se producen redes.

**Piloto** (en la nube, informativo, no entra en el veredicto): una RP³ de N = 192, con ℓ_G = 2.4 y protocolo poco denso, hasta τ_B = 2τ_A. Comprueba el coste, que las cuerdas desaparezcan y que el contador se estabilice. Si el piloto muestra que las condiciones de validez no se pueden cumplir, se dirá en una enmienda **antes** de producir.

## Condiciones de validez de cada red (se comprueban antes de clasificar)

1. **Cuerdas:** ξ(τ_B) < 1 % de ξ₁.
2. **Paredes:** fracción de pared en τ_B < 5 % de su máximo.
3. **Contador:** sitios lejos del vacío (cos θ < 0.9) < 1 % en τ_B.
4. **Meseta:** |N_lin(τ_B)/N_lin(ℓ_B − 0.5) − 1| < 3 %.

## Época y η (se informan, sin veredicto propio)

- **Época:**
  - r_W = a/a₁ cuando la fracción de pared cae por debajo del 10 % de su máximo;
  - r_N = a/a₁ cuando N_lin llega al 90 % de su valor final.
- **η_lat** = (27/8π)·K_lat·c_E/(ξ₁ℓ₁), con:
  - K_lat = [N_lin(τ_B) − N_kin(ℓ₁)]/N_mis − 1, que es el R/R_V2 − 1 de U19;
  - c_E = ⟨√(1 + (k/(m_G·a₁·r_W))²)⟩, por número, con los cuantos de k < k_cut.
- **c_p** = ⟨k⟩/(m_G·a₁·r_W), por número.
- **Tendencia con ℓ:** γ − 1 = ln[Y(ℓ_G = 3.2)/Y(ℓ_G = 2.4)] / ln(4.30/3.50), con su error.
  - γ = 1 significa que η no depende de ℓ, como supone U17.
  - γ < 1 significa que η baja con ℓ, y extrapolar a ℓ ≈ 90 favorecería al modelo.

## Veredicto (RP³, ℓ_G = 3.2; Y = media de las 4 redes ± error)

El error suma en cuadratura la dispersión entre redes/√4, κ y la anarmonicidad δ de U19. La extrapolación a ℓ ≈ 90 se hace con γ = 1 (η constante), **el lado conservador:** si η baja con ℓ, el resultado solo mejora.

- **COLAPSO_DURO:** Y + 1σ ≤ 0.59. La rama viva pasa incluso el umbral estricto de δA.
- **COLAPSO_INTERMEDIO:** ni DURO ni BLANDO. La rama viva pasa solo la puerta primaria.
- **COLAPSO_BLANDO:** Y − 1σ > 1.14. La rama viva falla la puerta primaria con las hipótesis de U17d.
- **INDETERMINADO:** falla alguna condición de validez en alguna red primaria, o el error relativo de Y supera el 50 %.

**Se informa también:**
- **la rama cómoda (B4a):** con X = 1800·Y, haría falta Y + 1σ ≤ 0.10;
- Y con ℓ_G = 2.4;
- η_lat, c_p, r_W y r_N;
- f_E,duro;
- U(1) frente a Hiramatsu et al. (k̄/m ≈ 3.1);
- el espectro final dN/d ln(k/m_G·a₁)/N_mis. Será la entrada de una U17e que pase la población medida por CLASS.

## Lo que U19b no resuelve

- **ℓ₁ ≈ 3.5–4.3, no ≈ 90,** y m_r/m_G ≤ 25 en la red, frente a ~e⁹⁰ en el universo. La parte «dura» del colapso físico, con p ≳ 10³ m_G, no cabe en la red: f_E,duro es una cota inferior. La extrapolación a ℓ ≈ 90 se apoya en γ, medido con una palanca pequeña (ℓ₁ de 3.5 a 4.3).
- **R_ρ** se usa para traducir Y en X de B2. La red mide RP³ directamente, pero la ξ de B2 a ℓ ≈ 90 sigue siendo la hipótesis de U17d.
- **Masa constante,** como en el modelo. No se calibra con el caso de masa dependiente de la temperatura del axión QCD.
- **Caja de 2τ_A:** en la fase tardía la caja es menor que el horizonte. Para los modos con k ≳ m_G·a no importa, pero puede afectar a los últimos defectos.

## Enmienda 1 · 15 de octubre de 2026 (Madrid; la hora es la de su commit): masa efectiva y contador adiabático, antes de la validación formal y de cualquier red

Al preparar la validación se hicieron pruebas de diagnóstico con campos homogéneos (N = 8 y N = 128, sin defectos). No se produjo ninguna red. Se declaran dos cosas que obligan a cambiar el método antes de validar.

**1. El núcleo congelado baja un poco la masa de los Goldstones.** Con a_p < a, la inclinación desplaza el vacío fuera de la variedad simétrica, porque la jerarquía m_r·a_p/(m_G·a) se reduce. Eso cambia la masa del Goldstone.
- **Cálculo exacto** (`meff_ratio`, con el hessiano del potencial de `kern.c` en el vacío desplazado, en doble precisión): m_eff/m_G = 0.9928 en τ_B con ℓ_G = 3.2, y 0.985 con ℓ_G = 2.4. Es igual en RP³ y en U(1). Sin congelación (a_p = a), 0.9992 y 0.996.
- **La V3 original estaba mal planteada:** medía la frecuencia con el intervalo medio entre ceros sobre 5 periodos, con una frecuencia que crece un ~30 % en la ventana, y ese estimador tiene un sesgo del ~0.8 %. El estimador WKB, (s·m/2)(τ_{k+1}² − τ_k²) = π entre ceros sucesivos, no tiene ese sesgo. Coincide con el hessiano a ≤ 0.17 % en las ocho comprobaciones (dos modelos, dos ℓ_G y dos ventanas).
- **Cambio:** el contador usa m_eff(τ) del hessiano, y no m_G.
- **Nueva V3:** el hessiano y el WKB coinciden a ≤ 0.3 % en x = 2 y en x ≈ τ_B/τ_A, y |m_eff/m_G − 1| ≤ 3 % en τ_B.
- La emisión del colapso ocurre sobre todo con x ≲ 2, donde m_eff difiere de m_G menos de un 0.5 %, así que el efecto en c_p y c_E es despreciable.

**2. El contador necesita el número adiabático de primer orden.** El número lineal de U19, (|χ̇|² + ω²|χ|²)/2ω, tiene un error que oscila con la fase, de orden H/2m. En el control homogéneo, el número entre τ_A (m_G/H = 10) y τ_B variaba un −12 %. Con la variable conforme u = aχ bajaba al −4.4 %, y con la corrección WKB de primer orden, **χ̇_k → χ̇_k + H·(1 + m²/2ω_k²)·χ_k** (`linear_number(..., hubble=H)`), queda en **+0.2 %**, con R/R_EDO − 1 = +0.8 %. El error restante es O((H/m)²).

**V2 y V3** se evalúan con este contador. Ni el criterio ni los umbrales del veredicto ni las redes cambian.

## Enmienda 2 · 15 de octubre de 2026 (Madrid; la hora es la de su commit): la validación con el método de la enmienda 1 falló; normalización canónica del contador

**La validación formal con el método de la enmienda 1 FALLÓ.**
- V1, V4 y V5 pasan.
- En V2/V3, el número de los controles homogéneos varía entre τ_A y τ_B un −1.8 % a −4.6 % en 6 de los 8 casos (umbral: 1 %).
- En U(1) con ℓ_G = 2.4, R/R_EDO − 1 = −5.2 % y −5.0 %.

Está archivada con su tabla en `network/results/u19b/_superseded/validacion_enmienda1/`. Por la regla del preregistro, U19b quedó INVALIDO en esa forma, y **no se produjo ninguna red.**

**Causa.** Con el vacío desplazado, el campo de Goldstone canónico no es el que usaba el contador:
- **U(1):** χ = r₀·θ, con r₀ = |φ_vacío| > v. El contador usaba v·θ, con lo que el número salía bajo por un factor (v/r₀)², que llega al −6 % con ℓ_G = 2.4.
- **RP³:** Q = Λ·nnᵀ + μ·(1 − nnᵀ), así que χ_a = 2f·(Λ − μ)·n_a. La velocidad del contador, 2f·(Q̇n)_a, ya llevaba el factor; el campo no.

**Cambio** (`linear_number(..., zchi, zchid)`). r₀ y Λ − μ salen de la misma minimización del potencial que da m_eff (`vacuum`). Valen 1.015 con ℓ_G = 3.2 y 1.031 con ℓ_G = 2.4, en τ_B.

**Prueba antes de revalidar,** con tres controles con arranque temprano:

| Control | R/R_EDO − 1 | Deriva |
|---|---|---|
| RP³, ℓ_G = 3.2 | +0.84 % | +0.10 % |
| U(1), ℓ_G = 2.4 | +0.84 % | +0.01 % |
| RP³, ℓ_G = 2.4 | +0.70 % | −0.02 % |

Se repite la validación completa, V1–V5, con los umbrales de siempre. Ni el veredicto ni las redes cambian.

## Enmienda 3 · 15 de octubre de 2026 (Madrid; la hora es la de su commit): dos definiciones del análisis, cambiadas tras el piloto y antes de producir

**El piloto** (RP³, N = 192, ℓ_G = 2.4, poco densa; 14 min en la nube; `network/results/u19b/pilot/`, `U19B_PILOT.json`) funciona y cumple las cuatro condiciones de validez:
- las cuerdas desaparecen en ℓ ≈ 5.5;
- las paredes también;
- al final, el 0.06 % de los sitios está lejos del vacío;
- el contador queda plano a < 1 %.

**El análisis tenía dos fallos de definición.** Al verlos se cambian, y se declara que se cambiaron **después de ver el piloto**; sus números, antes y después, quedan aquí:

1. **X contaba también los cuantos radiados antes de ℓ₁.** Son los del régimen de escalado, que luego se vuelven lentos. K_lat ya los descontaba (N_fin − N_kin(ℓ₁)), y en U17 esa emisión es una componente aparte (f_em). Ahora X resta su momento Σk·n_kin, tomado del espectro cinético de `measure` interpolado en ℓ₁, con el mismo corte k < k_cut, y c_p y c_E también los descuentan. En el piloto, Y_lat pasa de 4.44 a **3.23**.
2. **r_N saltaba antes de tiempo.** El primer instante en que N_lin ≥ 0.9·N_fin caía donde el contador aún no es válido. Ahora es el primer a/a₁ desde el que N_lin queda a menos del 10 % de su valor final. En el piloto, r_N pasa de 0.91 a 2.05; r_W = 1.51.

**Lo demás del piloto,** a título informativo: K_lat = 0.89, c_p = 4.8, η_lat = 3.1 y f_E,duro = 0.23.
- η_lat > 1: el número del colapso supera el que permite la energía de las cuerdas, (8π/27)ξ₁ℓ₁. Con ℓ₁ = 3.5, la energía de las paredes es comparable a la de las cuerdas: η_lat mide el colapso respecto a esa normalización, no una fracción.
- **El piloto avisa de que Y_lat puede salir muy por encima del caso U(1)** (Y ≈ 0.94). Si es así, el veredicto con γ = 1 (η constante, la extrapolación conservadora) sería COLAPSO_BLANDO, y la tendencia con ℓ (γ) decidirá cuánto vale esa extrapolación.

Los umbrales, el veredicto y las redes no cambian.

## Enmienda 4 · 17 de octubre de 2026, ≈17:50 (Madrid; la hora es la de su commit): un diagnóstico informativo para ℓ_G = 2.4 y U(1), escrito tras ver las cuatro redes primarias y antes de ver ninguna otra

**Qué se ha visto ya.** Las cuatro redes primarias (RP³, ℓ_G = 3.2) están terminadas y analizadas. Su veredicto se da con la lógica de arriba, sin cambios, en `docs/informes/INFORME_U19b.md`. La red RP³ con ℓ_G = 2.4, poco densa, semilla 91, está corriendo y no se ha mirado.

**Lo que motiva el diagnóstico.** En las cuatro redes, X_lat (el momento blando del colapso por cuanto de desalineamiento) es casi el mismo con los dos protocolos:
- poco densas, 8.20 y 7.66; densas, 7.25 y 7.50;
- en cambio, ξ₁ difiere un factor 3.2 (0.73–0.78 frente a 2.36–2.48), y ξℓ, promediado sobre el colapso (de ℓ₁ a ℓ₁ + ln r_W), un factor 2.6.

Y_lat = X_lat/[(8π/27)·ξ₁ℓ₁] supone que el rendimiento blando escala con la energía de las cuerdas, como en U17. Las redes de ℓ_G = 3.2 no lo muestran. La energía de más que llevan las redes densas va sobre todo a cuantos duros (k ≥ 10·m_G·a₁: momento total 24 frente a 15, en unidades de m_G·a₁·N_mis).

**Diagnóstico (informativo; no cambia el veredicto ni sus umbrales).** Con las redes de ℓ_G = 2.4 (y U(1) con ℓ_G = 3.2 como control), por protocolo y con la media de semillas:
- **D1 (protocolo):** el cociente X(densa)/X(poco densa) frente al cociente de ξ̄ℓ sobre el colapso.
- **D2 (masa):** ρ_X = X(2.4)/X(3.2), frente a ρ_ξℓ = (ξ₁ℓ₁)(2.4)/(ξ₁ℓ₁)(3.2), por protocolo.

**Lecturas, fijadas ahora:**
- **(a) Escalado de U17 (Y universal):** |ρ_X/ρ_ξℓ − 1| ≤ 0.15 en los dos protocolos.
- **(b) Rendimiento blando independiente de las cuerdas:** |ρ_X − 1| ≤ 0.15 en los dos protocolos, con |ρ_ξℓ − 1| > 0.3.
- **(c) Mixta:** cualquier otro caso.

**Qué se hará con cada lectura.**
- Con (a), la extrapolación de U19b (Y constante, γ medido) queda como está.
- Con (b), el veredicto de U19b sigue siendo el que es, pero se preregistrará una U17e con dos filas por cada F: la conservadora, con el Y medido, y otra con el X medido directamente, sin escalar con ξℓ. Su lógica se fija antes de calcular.
- La U17e con el espectro medido estaba ya prevista en este preregistro.

**Aclaración a la enmienda 4 · 17 de octubre de 2026, ≈19:06 (Madrid), antes de ver ninguna red de ℓ_G = 2.4.** Si se cumplen a la vez (a) y (b), la lectura es (c), mixta. ρ_ξℓ usa ξ₁ℓ₁ de cada red, como dice el texto. El diagnóstico lo calcula `network/u19b_analyze.py` (`diagnostic_enmienda4`), y da su lectura solo con las 4 redes de ℓ_G = 2.4.

## Enmienda 5 · 18 de octubre de 2026, ≈01:09 (Madrid): descomposición afín del rendimiento blando, informativa, antes de ver las redes densas de ℓ_G = 2.4

**Origen.** La revisión de Gemini (`audit/gemini_20261018/`) propuso leer X como X = X₀ + y_s·Z, con Z = (8π/27)·ξ₁ℓ₁:
- X₀ sería la parte que no depende de las cuerdas, el colapso de las paredes;
- y_s, el rendimiento por unidad de energía de cuerdas.

Es la versión cuantitativa de las lecturas (a) y (b) de la enmienda 4. Los números de la revisión estaban mal: suponía que las redes densas tienen más X, y tienen menos.

**Lo ya visto, y se declara.**
- **Con las 4 redes de ℓ_G = 3.2:** X₀ = 8.20 ± 0.30 e y_s = −0.085 ± 0.042, por mínimos cuadrados con 2 grados de libertad. Es decir, con más cuerdas no hay más cuantos blandos.
- **Con ℓ_G = 2.4** solo existe la red poco densa de la semilla 91 (X = 6.54). Las densas no se han visto.

**Qué se calcula** (`u19b_analyze.py`, `diagnostic_enmienda4.affine`): el mismo ajuste en cada ℓ_G con las dos semillas de cada protocolo, con sus errores.

**El X del universo** (B2, ℓ ≈ 90, donde Z vale 162·Y en unidades de X) se informa con tres extrapolaciones de X₀, fijadas ahora:
1. **X₀ constante:** el término de paredes no depende de ℓ. Es la expectativa física, A·σ/t frente a ρ_mis.
2. **X₀ lineal en ℓ₁,** con la pendiente entre ℓ_G = 2.4 y 3.2. Es la pesimista.
3. **X₀ ∝ ℓ:** la ordenada también escala como las cuerdas. Es el límite que equivale al escalado de U17.

En las tres se suma 162·max(y_s + 2σ, 0) para la parte de cuerdas. Cada X_univ se compara con los umbrales de U19b: X ≤ 95 pasa el umbral estricto de δA, y X ≤ 185 la puerta primaria.

**Lo que no cambia:** el veredicto de U19b, sus umbrales y la lectura de la enmienda 4. Si se usa una de estas extrapolaciones para unas filas nuevas de CLASS (una U17f), se preregistrará antes cuál.
