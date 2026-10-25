# U19 · c/η del colapso cuerdas–paredes, medido en la red · Preregistro

*14 de octubre de 2026, 01:10 (Madrid). Escrito antes de modificar el kernel y antes de producir ninguna red de U19. Se ejecuta en la nube.*

## Motivo

U17 mostró que la contribución del colapso de la red en H ≈ m_G vale K_eff = (8π/27)·ξℓ·η/c, y que el veredicto depende de c/η. c/η es la energía de la red en H = m_G/3 dividida por m_G veces el número de Goldstones que produce el colapso. Para U(1) con N_DW = 1, arXiv:1202.5851 da c ≈ 3.3 con η = 1 por hipótesis. Para RP³ no hay nada. U19 lo mide en la red, sin suponer η.

## Modelo

**Ruptura explícita añadida al kernel** (`network/kern.c`, `set_tilt`; 0 por defecto, y entonces la evolución es idéntica bit a bit a la de siempre):
- RP³: V → V − ε Q₁₁. En el vacío Q = nnᵀ, V_ε = −ε cos²φ (φ, el ángulo de n con e₁), y los tres Goldstones tienen m_G² = ε.
- U(1): V → V − ε′ φ₁, es decir −ε′ v cos θ, con m_G² = ε′/v (v² = 1/2).
- En los dos casos N_DW = 1: cada cuerda acaba en una pared.

**Masa:** m_G = m_r e^(−2.8). H = m_G/3 en ℓ₁ = 2.8 + ln 3 = 3.90. Al final, m_G/H ≈ 9.8 (N = 320, ℓ_fin = ln 160 = 5.08) y la pared mide ~16 celdas.

## Validación (antes de la producción)

- **V0:** con ε = 0, `make test` pasa (incluida la reproducción bit a bit de una red de U2).
- **V1 (masa):** en Minkowski (s = 0), una rotación homogénea pequeña desde el vacío oscila con ω² = m_G² a ≤ 1 %, en RP³ y en U(1).
- **V2 (contador de partículas):** red sin defectos de N = 128, con un ángulo casi homogéneo de θ₀ = 0.3 rad más un ruido pequeño, en radiación desde H = 30 m_G. El número comóvil final de Goldstones masivos, N = a³ Σ_k 2 dE_k/ω_k con ω_k = √(k²/a² + m_G²), debe coincidir a ≤ 5 % con el de la EDO homogénea θ'' + (2/τ)θ' + a²m_G² sin θ = 0 integrada desde el mismo θ₀.

Si V1 o V2 fallan, U19 se informa como **INVALIDO** y no se producen redes.

## Producción (`network/u19_produce.py`)

N = 320, τ_fin = 160, m_r·a_fin·dx = 1, dτ = 0.25 y 40 instantáneas con espectro desde ℓ = 3.5. Dos protocolos, para tener dos densidades de red en ℓ₁:
- **poco densa:** el de U11 (ξ_i = 0.5 en ℓ_i = 1, drenaje hasta ℓ_d = 3.8);
- **densa:** el de U10b (ξ_i = 6 en ℓ_i = 3.5, drenaje hasta 3.8).

Semillas 81 y 82 de cada protocolo, en RP³ y en U(1): **8 redes**, ~4 h con 2 hilos.

## Medidas (`network/u19_analyze.py`)

Por red:
- ξ₁ = ξ(ℓ₁), por interpolación;
- la fracción de volumen de pared W(ℓ): sitios con cos θ < −1/2 (RP³: Q₁₁ < 1/4);
- el número comóvil de Goldstones N(ℓ), con el espectro de `measure` (máscara de 2.5 anchuras de núcleo);
- la energía potencial de inclinación antes de que el campo se mueva, ⟨V_ε − V_min⟩, medida en la primera instantánea con H > m_G (campo congelado).

Con eso se define:
- **R ≡ m_G·[N(ℓ_fin) − N(ℓ₁)] / (⟨V_ε − V_min⟩ a₁³):** el número producido después de ℓ₁, en unidades del desalineamiento con la normalización de la cadena;
- **K_lat ≡ R/R_V2 − 1:** respecto al número de un campo sin defectos (V2), que absorbe el factor de arranque y la anarmonicidad débil.

**Ajuste**, por modelo: K_lat = A + B·(ξ₁ℓ₁) con las 4 redes (2 protocolos × 2 semillas). Entonces **c/η = (8π/27)/B**. El error sale de un remuestreo sobre las redes y de la propagación.

Se informa también la fracción de pared en ℓ_fin respecto a su máximo: si es > 5 %, el colapso no terminó y se dice.

## Veredicto (para RP³; U(1) se compara con arXiv:1202.5851, a título informativo)

- **C_ETA_BAJO:** c/η + 1σ ≤ 6. Es el lado tibio de U17: el veredicto K_SW_FUERA_DEL_RANGO_TIBIA se refuerza.
- **C_ETA_MEDIO:** c/η entre 6 y 20, con 1σ dentro.
- **C_ETA_ALTO:** c/η − 1σ ≥ 20. El colapso pesa menos de lo que U17 suponía en su caso central.
- **INDETERMINADO:** B < 2σ_B (sin señal de la red), el colapso no terminó en alguna red, o el error relativo de c/η supera el 50 %.

## Lo que U19 no resuelve

- **ℓ₁ ≈ 3.9, no 85.** Se supone que c/η no depende de ℓ: el número del colapso crece con la energía de la red, ∝ ξℓ, y c/η es la constante de proporcionalidad. No se comprueba.
- **Una sola masa,** m_G/m_r = e^(−2.8).
- **Solo cuatro redes por modelo** y dos densidades: el ajuste lineal es la hipótesis más simple.

## Enmienda 1 · 14 de octubre de 2026, 01:40 (Madrid): contador de partículas al final, antes de ejecutar nada

Al escribir el código, antes de cualquier ejecución (incluidas V1 y V2), se vio un problema del contador previsto para N(ℓ_fin). El espectro de `measure` solo tiene la energía cinética de la velocidad tangente, y 2·dE/ω cuenta bien las ondas relativistas. Tras el colapso, en cambio, los Goldstones masivos no relativistas oscilan en fase (sobre todo el modo k ≈ 0 del desalineamiento), y su energía pasa de cinética a potencial y vuelta. Una sola instantánea daría un número que oscila con la fase.

**Cambio.** N(ℓ_fin) se mide con el número de un campo libre linealizado en torno al vacío n = ±e₁ (RP³) o φ = v (U(1)), que no depende de la fase:
- RP³: χ_a = 2f·n_a (a = 2, 3, 4), con n orientado para que n₁ > 0, y χ̇_a = 2f (Q̇ n)_a;
- U(1): χ = v·θ;
- con ω_k² = k_lat²/a² + m_G², donde k_lat² = Σ_i 4 sin²(k_i/2) es la dispersión del laplaciano de la red,
- n = (1/V²) Σ_k (|χ̇_k|² + ω_k² |χ_k|²) / (2ω_k), y N = a³ n.

Las normalizaciones son las del kernel: ½χ̇² = |ṅ|², ½m_G²χ² es la inclinación cerca del vacío, y ½(∇χ)² el gradiente.
- **N(ℓ₁) se sigue midiendo con el espectro cinético**, como estaba previsto: en ℓ₁ el campo no está cerca del vacío y lo que se cuenta son cuantos relativistas ya emitidos.
- **V2 usa el mismo contador lineal al final.**
- **Comprobación añadida:** en ℓ_fin se informan los dos contadores. La fracción de sitios lejos del vacío (RP³: Q₁₁ < 0.9; U(1): cos θ < 0.9) debe ser < 1 %, o el contador lineal no es válido y se dice.

No cambian ni el veredicto, ni el ajuste, ni las redes.


## Corrección de horas · 14 de octubre de 2026, 00:44 (Madrid)

Las horas de la cabecera (01:10) y de la enmienda 1 (01:40) están mal. Las reales son las de los commits: preregistro en `4f1d9d3` a las 00:38:34 y enmienda 1 (con el cambio del kernel) en `db326f6` a las 00:41:00. No se ha ejecutado nada de U19 antes de la enmienda 1. El orden no cambia.

## Enmienda 2 · escrita antes de ejecutar nada (la hora es la de su commit): un control por arranque

En el protocolo denso la red arranca en ℓ_i = 3.5, cuando ya H = m_G/2. El campo, congelado por la condición inicial (π = 0), empieza a oscilar en cuanto se suelta, antes que en el protocolo poco denso, donde lo hace hacia 3H ≈ m_G. El rendimiento del desalineamiento por unidad de energía inicial depende de ese arranque, y la diferencia se colaría en la pendiente B.

**Cambios:**
1. **V2 se hace con los dos arranques:**
   - (a) desde H = 30 m_G, como estaba;
   - (b) congelado hasta ℓ_i = 3.5 y soltado ahí, como el protocolo denso.

   Cada protocolo usa su control: K_lat = R/R_V2(arranque) − 1. El criterio de V2 (≤ 5 % frente a la EDO) se aplica a los dos.
2. **⟨V_ε − V_min⟩** se mide en la última instantánea con ℓ ≤ 2.5 (H ≥ 1.35 m_G, campo aún congelado) en el protocolo poco denso, y en la condición inicial en el denso.
3. **Sistemática de anarmonicidad.** Los controles tienen θ₀ = 0.3 (régimen armónico), mientras que en las redes los ángulos son aleatorios y grandes. Se añade un control homogéneo de ángulo grande, θ₀ = 2.5 rad (N = 32), con los dos arranques. Si R(2.5)/R(0.3) cambia entre arranques en una cantidad δ, δ se suma en cuadratura al error de K_lat de las redes densas.

No cambian ni el ajuste, ni los umbrales, ni las redes.

## Enmienda 3 · 14 de octubre de 2026 (Madrid; la hora es la de su commit): una sola RP³ de confirmación

*Escrita después de ver las cuatro redes U(1) (informativas) y un minuto después de lanzar la red RP³ que describe. Se declara.*

**Lo que dieron las U(1)** (`network/results/u19/U19_RESULT.json`, provisional). Al final de la caja (m_G/H ≈ 10) las paredes ya no están (≤ 4 % de su máximo), pero entre el 12 y el 29 % de los sitios siguen lejos del vacío (cos θ < 0.9). El desalineamiento oscila todavía con amplitud de ~1 rad. Por tanto, la comprobación de la enmienda 1 (< 1 %) falla, y el contador lineal no es válido. K_lat sale entre 3 y 8, con una dispersión entre semillas del orden de su valor.

**Cambio de la producción.** Las cuatro redes RP³ costarían ~4 h y casi seguro fallarían la misma comprobación, porque el desalineamiento evoluciona igual en los dos modelos. Se produce una sola RP³ (poco densa, semilla 81) para confirmarlo:
- si su fracción lejos del vacío al final es > 1 %, la producción se detiene y **U19 se informa como INDETERMINADO (contador no válido)**, sin ajuste;
- si fuera < 1 %, se producen las otras tres, como estaba previsto.

No cambian ni los umbrales ni el veredicto.

## Nota 1 · 14 de octubre de 2026, ≈21:10 (Madrid): vetos antes de clasificar el ajuste (auditoría de Codex)

`u19_analyze.verdict` podía clasificar un ajuste aunque fallaran dos cosas que el análisis guardaba aparte:
- la validación;
- la validez del contador lineal, `far_frac_end_ok`.

Codex lo mostró con un ejemplo sintético: un 50 % de sitios lejos del vacío daba C_ETA_BAJO. Ahora, antes de mirar el ajuste:
- con la validación fallida, el veredicto es INDETERMINADO (validación fallida);
- con far_frac_end ≥ 0.01 en alguna ejecución, INDETERMINADO (contador no válido);
- si falta N_kin al final o K_lat no es finito, INDETERMINADO (ejecución incompleta).

**El veredicto publicado no cambia:** INDETERMINADO (contador no válido). El submodelo U(1) pasa de «INDETERMINADO (B < 2 sigma)» a «INDETERMINADO (contador no válido)», que es la razón previa: sus cuatro ejecuciones tienen far_frac_end > 0.01.
