# ESPEC_U67b · La memoria de la era de materia en una red de RP³ ya ordenada

*Escrito el 2026-11-21 a las 16:14 (Madrid). U67 sigue en producción (4 de 12 redes hechas). Este ESPEC se escribe después de ver sus primeras redes y una prueba de viabilidad, y antes de calcular U67b.*

## Por qué

**Lo que se ha visto en U67 (MD2, semilla 101).**
- Desde la fase simétrica O(4), con una semilla de amplitud 0.02 y correlación de 1 celda, Q no se ordena hasta ℓ ≈ 4.1: la brecha media sigue en 0.004–0.03 hasta ℓ ≈ 3.9.
- El modo taquiónico del punto simétrico es débil, m² = −0.019 frente a m_r² = 0.24. Solo crecen longitudes de onda de ≳ 100 celdas.
- U67 no podrá decir nada de la memoria: su puerta G2 va a fallar.

**La prueba de viabilidad (96³, MD con ℓ_tr = 2, una semilla).** Se hizo después de ver eso, en el directorio temporal y sin veredicto:

| semilla | ℓ al que la brecha llega a 0.5 |
|---|---|
| amplitud 0.02, correlación 1 celda | nunca, dentro de la ventana |
| 0.1, 1 celda | nunca |
| 0.3, 1 celda | 3.82 |
| 0.1, 6 celdas | 3.37 |
| 0.3, 6 celdas | 3.00 |

Ordenarse desde la fase simétrica deja, como mucho, ~1 unidad de ℓ de radiación en estas cajas.

**Lo que hace U67b.** La parte del supuesto abierto que sí se puede probar en el retículo es la otra: que una red ya formada, que pasa parte de su vida en la era de materia (la del inflatón) antes de llegar a la radiación, olvida esa era. Eso hace U67b. El ordenamiento desde la fase simétrica queda como lectura informativa.

## Los protocolos (192³, dτ = 0.25, τ_end = 96, a(τ_end) = 1/m_r; semillas 101, 102 y 103)

- **RD_SPARSE y RD_DENSE:** se reutilizan las redes de U67, el mismo protocolo de producción. Si U67 no las termina, se producen aquí.
- **MDA2 y MDA3 (las del veredicto):** el fondo de U67, de materia hasta ℓ_tr = 2.0 o 3.0 y de radiación después.
  - Arranque sintético en el vacío de RP³ (`lattice.initial`), en ℓ_i = 1.3, después de ℓ_crit = 1.11, con π = 0.
  - Correlación l_c = (τ₀/18.65) √(K_ξ/0.5), la regla del productor.
  - Con ξ_Q R encendido (ξ_Q = 1/4, X₀ del punto simétrico).
  - Drenaje de la velocidad masiva hasta ℓ_d = 3.8 a la tasa m_r.
- **MDO2 (informativo):** como MD2 de U67, pero con una semilla de amplitud 0.3 y correlación de 6 celdas.

## Las puertas

- **G0 y G1:** las de U67 (`U67_CHECKS.json`), que ya pasaron.
- **G2:** en las 6 redes MDA, la brecha media es ≥ 0.5 en todas las instantáneas con ℓ ≥ 2: la red sigue en el vacío.
- **G3:** ≥ 200 plaquetas atravesadas en la última instantánea de las redes MDA y RD.

## El veredicto

Se define D_k = ξ_MDAk(ℓ_c)/ξ_RD_SPARSE(ℓ_c) − 1, con ℓ_c = 4.10, la misma convención de ξ y el mismo bootstrap que U67.

- **OLVIDA.** Pasan G0–G3, |D_2| ≤ 0.15 y |D_3| ≤ 0.15.
- **RECUERDA.** Pasan G0–G3, y para algún k, |D_k| > 0.30 y |D_k| > 3σ_k.
- **INTERMEDIO.** Pasan G0–G3 y no se cumple ninguno de los dos.
- **INCONCLUSO.** Falla alguna puerta.

**Lo esperado:** OLVIDA en MDA2. En MDA3 es posible INTERMEDIO, porque solo pasa 1.1 unidades de ℓ en radiación.

## Salidas

- `network/u67_produce.py --u67b` escribe las redes en `network/results/u67b/runs/`.
- `network/u67_analyze.py --u67b` escribe `network/results/u67b/U67B_RESULT.json`.
- El informe va en `docs/informes/INFORME_U67.md`.
