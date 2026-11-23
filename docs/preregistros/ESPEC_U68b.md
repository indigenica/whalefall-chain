# ESPEC_U68b · U68 con la puerta de las masas corregida

*Escrito el 2026-11-21 a las 21:47 (Madrid). Es posterior a las comprobaciones G0–G2 de U68 y anterior a calcular cualquier red.*

## Por qué

Las comprobaciones de U68 (`network/results/u68/U68_CHECKS_U68.json`) dieron:
- **G0 pasa:** campos idénticos bit a bit al núcleo de `main`, en RP³ y U(1), con y sin núcleo congelado tardío.
- **G2 pasa:** el contador por modo da la ⟨k⟩ esperada a 10⁻⁵, sin fugas a los otros modos.
- **G1 falla, por dos motivos:**
  1. **Un fallo de la comprobación.** La medida de WKB cae en x = τ/τ_A ≈ 2.1–2.3, pero la tabla de masas efectivas solo llegaba a x = 2.0, y la interpolación se quedaba fija en su último valor. A N = 192, τ_A es menor que en U19b (N = 512) y el periodo ocupa más fracción de τ_A. Incluso en O3, que es el caso de U19b, salía −0.4 %.
  2. **Un criterio demasiado estricto.** La fórmula m_c² = m_G²(w₀ − w_c) es de primer orden. Para los modos de m²/m_G² = 4 (c = 3 en GEN y en O2), la oscilación en Minkowski se desvía −2.8 %. La inclinación grande desplaza el vacío radial, y la corrección es O(κ), como en `NOTA_ESPURION_E.md`.

  El contador no usa esa fórmula, sino la masa de la hessiana. Lo que hace falta comprobar es que esa masa sea la de la dinámica.

**Por tanto, U68 tal como se preregistró queda INCONCLUSO por G1.** U68b repite U68 con este único cambio.

## El cambio

- **G1 (U68b):**
  - En Minkowski, para cada modo y protocolo, la frecuencia homogénea coincide con la masa de la hessiana en el vacío de x = 1 (la primera fila de la tabla), con un error < 1 %.
  - En el universo con núcleo congelado tardío, la masa de WKB coincide con la de la hessiana a < 0.3 %, con la tabla extendida hasta x = 2.6. El paso es el mismo, 0.05, así que los valores entre x = 1 y 2, que son los que usa la producción, no cambian.
- **No cambia nada más:** los protocolos, las semillas, las medidas, G0, G2, G3 y toda la lógica del veredicto (D, los umbrales 0.15 y 0.30, y 3σ) son los de ESPEC_U68.
- Las comprobaciones de U68b se guardan en `U68_CHECKS.json`. Las de U68 quedan en `U68_CHECKS_U68.json`.

## Enmienda 1 (2026-11-21 a las 21:50 (Madrid), antes de calcular cualquier red)

Las comprobaciones de U68b (primera ejecución) dieron G1 bien en 7 de los 9 modos. En los dos de m²/m_G² = 4 (c = 3 en GEN y en O2), la masa de WKB se separaba de la de la hessiana en +1.1–1.5 %.

**La causa es la comprobación, no la producción.** La prueba de WKB arrancaba en Q = nnᵀ y no en el vacío desplazado por la inclinación. Los modos radiales oscilaban y arrastraban la frecuencia. Arrancando en el vacío de la teoría evolucionada (girado un ángulo pequeño en el plano (0, c)), la desviación baja a −0.015 % en GEN, c = 3 (probado en el directorio temporal, sin veredicto).

En la producción, el campo llega a ese vacío mucho antes del colapso, con el drenaje.

**Dos cambios, antes de cualquier red:**
1. **La prueba de WKB de G1** arranca en el vacío desplazado q*(x) (el mínimo de V_sym − x²ε Σ w_a q_a), girado 0.01 rad en el plano (0, c). El umbral sigue siendo 0.3 %.
2. **La normalización canónica, por modo.** Con masas distintas, los autovalores del vacío son Λ (en e₁) y μ_c (en e_c), distintos para cada c. U19b usaba z = Λ − μ, con un solo μ. Aquí se usa z_c = Λ − μ_c para el modo c.
   - En GEN, a x = 2, es z = (1.023, 1.058, 1.112), frente a la media 1.064. Con la media, el número de cada modo estaría mal en hasta un ±8 %.
   - En O3 los dos coinciden.

La lógica del veredicto no cambia.

## Enmienda 2 (2026-11-21 a las 23:10 (Madrid), después de la primera producción)

**Un error de código en el contador, encontrado por la revisión de Codex en la PR #103.**
- `linear_number_modes` multiplicaba por z_c = Λ − μ_c tanto el campo χ_c = 2F n_c como su velocidad χ̇_c = 2F (Q̇n)_c/a³.
- La velocidad ya lleva ese factor, porque Q̇n = (Λ − μ_c) ṅ en el vacío desplazado. U19b lo distinguía: `zchi = Λ − μ` y `zchid = 1`.
- Así, la parte cinética del número de cada modo quedaba multiplicada por z_c² de más: hasta un ≈ 24 % en el modo pesado de GEN al final.

**Qué se hace.**
- Se corrige el código para que z_c solo multiplique el campo.
- Se repiten las 9 redes, con las mismas semillas, protocolos y medidas.
- La primera producción y su análisis (IGUAL, D = −0.133 ± 0.013) se mueven a `network/results/u68/_superseded/produccion1_z_doble/`, con una nota. **El veredicto de U68b es el de la producción repetida.**

**Lo que no cambia:** G0–G2, que no usan esta normalización (G2 usa z = 1 y π = 0), las puertas y la lógica del veredicto.
