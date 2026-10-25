# U17b · Con el colapso de la red dentro, solo vive ξ constante, y como materia oscura casi tibia

*14 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17b.md` (commit `c7ff41b`).
- *Código:* `lya/u17b_class.py`. *Resultado:* `lya/results/u17b/U17B_RESULT.json`.
- CLASS v3.4.1 con el parche (C0b superada), las puertas de U12 sin cambios y las ~150 ejecuciones de referencia y familia recalculadas en la nube (~5.5 h con 2 procesos).

## 0. En una frase

**Veredicto preregistrado: SOLO_XI_CONSTANTE_VIVE.**
- **B2 (ξ constante) pasa con las tres variantes del colapso.**
- **B1 (ξ ∝ ℓ, GHV) falla con las tres.**
- **B4a (crecimiento medido) solo pasa si el colapso deja cuantos fríos** (c = 3 con η = 0.1).

Con c/η ≈ 3, el valor de la literatura para U(1), la rama que vive tiene casi toda la materia oscura en la población del colapso. Esa población se comporta en Lyα como una WDM de ~8 keV, justo en el borde de las cotas. Es una predicción nueva y comprobable.

## 1. Referencias

δA de las referencias: WDM de 5.3 keV, 0.336; mezclas de 1–4 keV, 0.234–0.311. Umbral estricto 0.2341, el mismo valor que en U12.

## 2. Filas (q = 0.87, x₀ = 10, R_Γ = 1.56)

| Fila | K_eff | m_G (eV) | T_osc | f no fría | δA | forma (m, f) | puerta primaria | solo δA |
|---|---|---|---|---|---|---|---|---|
| B2, c = 3, η = 1 | 34 | 5.4×10⁻¹⁸ | 82 keV | 0.97 | 0.247 | 7.8 keV, 0.97 | **pasa** | falla |
| B2, c = 30, η = 1 | 3.3 | 5.2×10⁻¹⁶ | 0.64 MeV | 0.77 | 0.211 | 3.7 keV, 0.26 | **pasa** | pasa |
| B2, c = 3, η = 0.1 | 3.3 | 5.2×10⁻¹⁶ | 0.64 MeV | 0.77 | 0.025 | 10 keV, 0.13 | **pasa** | pasa |
| B4a, c = 3, η = 1 | 383 | 5.0×10⁻²⁰ | 8 keV | 0.997 | 0.819 | 1.3 keV, 1.0 | falla | falla |
| B4a, c = 30, η = 1 | 35 | 5.3×10⁻¹⁸ | 82 keV | 0.97 | 0.812 | 1.3 keV, 0.98 | falla | falla |
| B4a, c = 3, η = 0.1 | 35 | 5.3×10⁻¹⁸ | 82 keV | 0.97 | 0.249 | 7.5 keV, 0.93 | **pasa** | falla |
| B1, c = 3, η = 1 | 716 | 1.4×10⁻²⁰ | 4 keV | 0.999 | 0.907 | (rms 0.57, no representable) | falla | falla |
| B1, c = 30, η = 1 | 65 | 1.6×10⁻¹⁸ | 45 keV | 0.985 | 0.901 | (rms 0.19, no representable) | falla | falla |
| B1, c = 3, η = 0.1 | 65 | 1.6×10⁻¹⁸ | 45 keV | 0.985 | 0.378 | 4.0 keV, 0.80 | falla | falla |

## 3. Lectura

- **El indicador de U17 acertó.** Las filas con λ_colapso/λ_WDM(5.3 keV) ≈ 5–13 (B4a y B1 con η = 1) dan δA = 0.81–0.91, como una WDM pura de ~1 keV, excluida de sobra. Las de B2 (λ/λ_WDM ≈ 0.6–0.8) pasan.
- **ξ constante vive, pero cambia de naturaleza.** Con c/η ≈ 3:
  - la materia oscura es en un 97 % la población del colapso (m_G ≈ 5×10⁻¹⁸ eV, T_osc ≈ 80 keV);
  - para Lyα equivale a una WDM de 7.8 keV;
  - pasa la puerta de forma (la cota térmica es 5.3 keV), pero su δA (0.247) supera en un 6 % el umbral estricto de las mezclas;
  - **es una materia oscura casi tibia, justo en el borde: la próxima generación de cotas de Lyα (m_WDM ≳ 6–10 keV) la decidiría.**
  - Con c/η = 30 vuelve a ser fría con una cola tibia del 26 % a 3.7 keV.
- **ξ creciente muere salvo que el colapso sea frío.** Solo con η = 0.1 (el 90 % de la energía de la red a cuantos duros, no a cuantos de p ~ m_G) pasa B4a, también por forma y en el borde. B1 no pasa ni así.
- **Lo que esto hace con U12:** «la rama fría sobrevive aunque ξ crezca» no se sostiene con el colapso contado: con η = 1, B4a y B1 fallan. **El sector oscuro vuelve a depender de que ξ se estabilice**, ahora casi con independencia de q, que pesa poco frente al colapso.

## 4. Lo que U17b no resuelve

- **c y η de RP³ no están medidos.** U19 no pudo medirlos con la caja de la nube (INDETERMINADO).
- **El espectro del colapso** se supone log-normal de anchura 0.5. Un espectro más ancho hacia el infrarrojo enfriaría la población; uno más duro la calentaría.
- **Solo q = 0.87, x₀ = 10 y R_Γ = 1.56.** Con el colapso dominando, q y R_Γ importan poco (f_em ≤ 10⁻⁴ en todas las filas).
- **m_G ≈ 10⁻²⁰ eV** (B4a y B1 con c = 3, η = 1) roza además la cota de la materia oscura ultraligera (2×10⁻²⁰ eV), que no se aplicó como puerta.
- **Minicúmulos:** 19–2×10⁴ M_⊙ en las filas que viven, muy lejos de los 0.2–10 M_⊙ de U14. U14 habría que rehacerlo para la rama viva.

## Nota 1 · 14 de octubre de 2026 (Madrid): lo que U17c añade a esta lectura

U17c (`docs/informes/INFORME_U17c.md` §4) muestra que, mientras K_eff ≫ 1, la longitud de libre recorrido de la población del colapso va como c·K_eff ∝ ξℓη: **c se cancela.** No cambia ningún número de este informe, pero sí dos frases:
- **«Con c/η ≈ 3, la rama viva equivale a una WDM de ~8 keV».** Vale para cualquier c con η = 1 mientras c/η ≲ 10: el control B2 con c = 1 de U17c da δA = 0.246, frente a 0.247 con c = 3. Con c/η = 30, K_eff ≈ 3, la población del colapso deja de dominar y la fila sale más fría.
- **«B4a solo pasa si el colapso deja cuantos fríos».** Lo que cuenta no es que los cuantos sean fríos (c pequeño), sino η: qué fracción de la energía va a cuantos de p ~ m_G. La rama necesita η ≤ 0.1, con el resto en cuantos con p ≳ 10³ m_G.

## Nota 2 · 14 de octubre de 2026 (Madrid): dos correcciones de la auditoría de Codex (las filas se rehacen en U17d)

`audit/codex_20261014/AUDITORIA.md` y `RESPUESTA.md`:

1. **Energía y momento mezclados.** Las filas de este informe contaron los cuantos con la energía c·m_G, pero les dieron **momento medio** c·m_G. Con c = 3 la energía de la población supera la asignada en un 6.6 %, y con c = 30, en un 0.07 %. Corregido en `row()`. Los resultados de arriba son los de la construcción antigua y se regeneran con `legacy_energy=True`. Con el número bien normalizado, B2 con c = 3 y η = 1 da m_G = 6.1×10⁻¹⁸ eV (en vez de 5.4) y la población es algo más fría.
2. **La energía de la red de RP³.** El inventario usa ξ de U(1), mientras que la emisión de RP³ se escala con R_Γ = 1.56. Lo coherente es escalar también la energía de la red, R_ρ ≈ R_Γ (ESPEC_U17d). Con eso, B2 con η = 1 queda más caliente: X ≈ 152 frente al umbral ≈ 95 de U17c. Se espera que falle Lyα.

**Lo que esto hace con la frase «solo vive ξ constante, como una WDM de ~8 keV»:** queda pendiente de U17d. Además:
- las masas de minicúmulos (19–2×10⁴ M_⊙) son **escalas candidatas**, no una población demostrada, porque la longitud de libre recorrido de la población del colapso (~0.05 Mpc) supera su radio comóvil (~0.005 Mpc);
- la comprobación de isocurvatura de U5c no vale para estas masas y queda pendiente.
