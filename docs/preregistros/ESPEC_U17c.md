# U17c · ¿Cuánto de «duro» tiene que ser el colapso para que viva la rama con ξ creciente? · Preregistro

*Escrito antes de calcular ninguna fila de U17c (la hora es la de su commit). Se ejecuta en WSL, en paralelo con U15b.*

## Motivo

U17b (`docs/informes/INFORME_U17b.md`): con el colapso cuerdas–paredes contado, B4a (crecimiento medido) solo pasa Lyα con η = 0.1, y B1 (GHV) no pasa ni así (δA = 0.38). Aquí η es la fracción de la energía de la red en H = m_G/3 que acaba en Goldstones de momento ~c·m_G; el resto se va a cuantos duros, que son radiación. La pregunta es qué η crítico necesita cada rama que crece. De eso depende si la «rama cómoda» (q < 0.91 con ξ creciente) sigue siendo razonable.

## Filas (`lya/u17c_class.py`; misma construcción que U17b)

- Ramas **B4a** y **B1**, con q = 0.87, x₀ = 10 y R_Γ = 1.56.
- c ∈ {1, 3} (momento medio del colapso / m_G en H = m_G/3).
- η ∈ {1, 0.3, 0.1, 0.03, 0.01}.
- Anchura de la log-normal σ_ln ∈ {0.5, 1.0}.
- Control: B2 con c = 1 y η = 1 (σ_ln = 0.5).
- Total: 2 × 2 × 5 × 2 + 1 = **41 filas**.

m_G autoconsistente con K_eff = (8π/27)ξℓ·η/c, la emisión de U12 y las puertas de U12 sin cambios (`lya/u5_class.gates`; umbral estricto δA = 0.2341; puerta primaria de forma si es representable). Las referencias y la familia salen de la caché de CLASS de U12 en WSL.

## Veredicto (B4a con c = 3 y σ_ln = 0.5; el resto se informa)

η_crit(B4a) es el mayor η de la rejilla con la fila permitida por la puerta primaria, interpolando δA en ln η donde haya cambio.
- **RAMA_COMODA_RAZONABLE:** η_crit ≥ 0.3. Basta con que una parte apreciable de la energía de la red se vaya a cuantos duros.
- **RAMA_COMODA_FORZADA:** 0.03 ≤ η_crit < 0.3. Exige que el colapso sea mayoritariamente duro.
- **RAMA_COMODA_MUERTA:** ninguna fila con η ≥ 0.03 pasa.

**Se informa también**, sin veredicto:
- η_crit de B1 y de B4a con c = 1, y con σ_ln = 1;
- m_G, T_osc, minicúmulos y δA de cada fila;
- si alguna fila con η pequeño queda con m_G > 10⁻¹⁵ eV y T_osc > 1 MeV (el régimen antiguo, con minicúmulos de U14).

## Lo que U17c no resuelve

- **η y c de RP³ siguen sin medirse** (U19 no pudo).
- **La energía que no va a la población del colapso** se supone radiación oscura sin efecto en Lyα. No se calcula su ΔN_eff; con η pequeño puede ser apreciable, y se dirá.
