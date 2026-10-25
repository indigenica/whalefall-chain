# U21 · ¿Qué F permite la cadena? Banda de la escala de Q frente a Lyα, isocurvatura y las demás cotas · Preregistro

*15 de octubre de 2026, ≈19:35 (Madrid). Escrito antes de calcular ninguna fila de U21 (la hora que cuenta es la de su commit). Se ejecuta en la nube, con la caché de CLASS de U12/U17d.*

## Motivo

**F no está fijado por la teoría.** Según C49 (`C49_RESULTS.json`, `NOTA_C46_C47_C49_BALLENA4.md`, que el equipo aportó el 15 de octubre):
- la restauración de Q durante la inflación solo exige μ/H_* ≤ √(12ξ − 9/4) = 0.87 con ξ = 0.25, es decir, F ≤ 2.7×10¹⁵ GeV;
- el límite de cuerdas en el CMB da F ≤ 1.26×10¹⁵ GeV;
- el borde inferior de C34 desaparece;
- **μ/H_* = 0.30, y con él F = 1.07×10¹⁵ GeV, es el benchmark declarado de C50,** pegado al techo.

**Por qué importa.** Un diagnóstico sin veredicto (`audit/referee_20261015/RESPUESTA.md` §11) muestra que m_G ∝ F⁻⁴. Con un factor 2–3 menos en F, la población del colapso se vuelve fría. Elegir F a mano, después de ver Lyα, sería ajustar a posteriori. **U21 calcula la banda de F que permiten todos los eslabones, sin elegir.**

## Filas (`lya/u21_fband.py`)

**Rejilla:** F/F₀ ∈ {1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03}, con F₀ = 1.07×10¹⁵ GeV. Equivale a μ/H_* ∈ [0.009, 0.30]. Todas las filas usan la población del colapso de U17d: c = 3, σ_ln = 0.5, el número normalizado por la energía media, q = 0.87, x₀ = 10 y R_Γ = 1.56.

**Variantes:**
- **P (primaria):** B2, η = 1, R_ρ = 1.56, r = 1. Son las hipótesis primarias de U17d.
- **H (colapso caliente):** B2, η = 1, R_ρ = 1.56, r = 3. Imita el Y ≈ 3 del piloto de U19b.
- **D (colapso duro):** B2, η = 0.3, R_ρ = 1.56, r = 1.
- **A (rama cómoda):** B4a, η = 1, R_ρ = 1.56, r = 1.

**Regla de frío trivial, fijada aquí.** Si λ_proxy/λ_WDM(5.3 keV) < 0.05, la fila pasa Lyα sin CLASS. λ_proxy es la longitud de libre recorrido de U17 con el momento medio, multiplicada por r.
- **Control:** la primera fila de la variante P que cumpla la regla se calcula también con CLASS, y debe dar δA < 0.05.
- Si el control falla, todas las filas que la cumplían se calculan con CLASS.

## Puertas por fila

1. **Lyα (propia):** la puerta primaria de U12 (forma estricta; δA estricto si la forma no es representable), sin cambios.
2. **Isocurvatura (de la literatura, no propia):** m_G ≥ 1×10⁻¹⁷ eV.
   - Es la cota de Chathirathas y Schwetz (arXiv:2511.15790) para ALPs postinflacionarios con la radiación de cuerdas como canal dominante. Solo se ha leído el resumen.
   - Se informa también la cota laxa, m > 1.73×10⁻¹⁸ eV (arXiv:2603.04401).
   - Esta puerta **se sustituirá por U20** (cálculo propio) cuando exista.
   - No se aplica la cota de modos tensoriales del mismo artículo (~5×10⁻⁷ eV): supone restauración térmica (f < T_max), y en la ballena la restauración es por curvatura (ξQ²R, C49).
3. **Radiación oscura:** ΔN_eff(emisión) ≤ 0.3.
4. **Cuerdas en el CMB:** F ≤ 1.26×10¹⁵ GeV. Se cumple en toda la rejilla.

**Se informa, sin puerta:**
- **Superradiancia,** con las bandas sin autointeracción de la literatura: agujeros negros estelares 6×10⁻¹³–10⁻¹¹ eV (Cardoso et al. 2018); supermasivos, indirectas, ~10⁻¹⁹–10⁻¹³ eV. Con F ~ 10¹⁴–10¹⁵ GeV, la autointeracción (λ ~ m²/F²) las debilita (Baryakhtar et al. 2021), y no hay una cota cuantitativa verificada para esta F. **Queda como cuestión abierta, que puede ser decisiva.**
- **LISA:** f_a = F/2 frente al umbral de `s5_closure` (5×10¹⁴ GeV).
- **Otros:** Λ_H = √(m_G·f_a), la escala del sector oculto que daría m_G; m_G; T_osc; las escalas de minicúmulos; f_nc.

## Veredicto (variante P)

La banda es el conjunto de valores de la rejilla que pasan las puertas 1–4.
- **BENCHMARK_PERMITIDO:** F₀ pasa.
- **BANDA_BAJO_EL_BENCHMARK:** F₀ no pasa, pero algún F < F₀ sí. Se dan F_max y F_min de la rejilla.
- **SIN_BANDA:** ningún F pasa.

H, D y A se informan con su propia banda.

## Lo que se espera, escrito antes de calcular

- **En F₀, la variante P** pasa Lyα solo por forma (U17d), pero m_G = 2.6×10⁻¹⁸ eV falla la isocurvatura estricta. Con F = 0.5·F₀, m_G ≈ 5.5×10⁻¹⁷ eV y λ/λ_WDM ≈ 0.22, y debería pasar.
- **Se espera BANDA_BAJO_EL_BENCHMARK** con F_max ≈ 0.5–0.7·F₀.
- **Coste anunciado:** con F < F₀, f_a queda por debajo del umbral de LISA en `s5_closure`, y **la señal de ondas gravitacionales de cuerdas deja de ser una predicción detectable.**

## Lo que U21 no resuelve

- **La isocurvatura** es la de la literatura (U20, pendiente).
- **La superradiancia con autointeracción** queda abierta.
- **Ordenamiento tardío de Q.** Con μ/H_* pequeño, Q se ordena más tarde tras la inflación (C49: H_orden ∝ μ). Las condiciones iniciales de la red cambian, pero el atractor de escalado no debería.
- **El colapso sigue con las hipótesis de U17d** (η, c, r, R_ρ): U19b las medirá.
