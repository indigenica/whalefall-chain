# INFORME U60 · La anisotropía con cuscutón: un Bianchi IX con Einstein + cuscutón

**Preregistro:** `docs/preregistros/ESPEC_U60.md` (17 nov 2026, 14:24, Madrid).
**Código:** `chain/u60_bianchi_cuscuton.py`, que usa `explora_cuscuton` sin cambios (Cusc.all, el potencial y H_meseta).
**Resultados:** `results/u60_runs/runs.json` (1 924 corridas), `results/u60_runs/produce.out` y `results/u60_bianchi_cuscuton.json`. Tardó 21 min.

Interno: no se exporta.

---

## 2026-11-17 14:52 (Madrid)

### 1. Veredicto preregistrado: ANISOTROPIA_INOCUA

| puerta | resultado |
|---|---|
| B0, el límite isótropo frente al `host` de U61b | **pasa:** χ_land a ≤ 3×10⁻⁸ y N_acc a ≤ 5.5×10⁻⁶ |
| B1, el modo de forma lineal | **pasa:** la deformación de Bianchi IX sigue a δ̈ + 3Hδ̇ + 8δ/a² = 0 a ≤ 2.4×10⁻⁷. Crece ×7.4–8.4 desde el giro hasta el rebote |
| B2, la ligadura (≤ 10⁻⁷) | **pasa:** solo el 0.4 % de las corridas con ε ≤ 10⁻² son inválidas |
| B3, la tolerancia | **pasa:** la misma clasificación, con \|ΔN\| ≤ 2×10⁻⁸ |

**Todos los hijos aguantan más que ζ₀** (f_rob = 1.0, 52 de 52). La tolerancia ε\* va de 10⁻⁴ a 3×10⁻³ (mediana 10⁻³) y baja con el tamaño de la madre:

| A·H | hijos | ε\* mínima | mediana | máxima |
|---|---|---|---|---|
| 30 | 7 | 10⁻³ | 3×10⁻³ | 3×10⁻³ |
| 48 | 11 | 10⁻³ | 3×10⁻³ | 3×10⁻³ |
| 100 | 13 | 3×10⁻⁴ | 10⁻³ | 3×10⁻³ |
| 300 | 14 | 3×10⁻⁴ | 3×10⁻⁴ | 10⁻³ |
| 1 000 | 7 | 10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |

### 2. Qué limita ε\*: no la supervivencia, sino los e-folds

- **Por encima de ε\*, casi todas las corridas (147 de 156) siguen dando un hijo** isótropo al final (Ω_σ ~ 10⁻¹⁰⁰). Lo que pasa es que su número de e-folds cambia en más de uno.
  - ΔN escala como ε², con un coeficiente de mediana ≈ −4.5×10⁵ y de hasta 10⁷ con las madres grandes.
  - **El mecanismo:** la cizalla en el rebote (σ²/ρ_r ≈ 1.5ε²) cambia lo que el rebote bombea en χ, y con ello dónde aterriza.
  - Un hijo con ΔN = −1.4 en ε = 3×10⁻³ sigue siendo un hijo con su inflación. Lo que cambia es su curvatura y su ℓ_mem, no su existencia.
- **La amplitud hasta la que el hijo nace e isotropiza** en las tres direcciones (sin la condición de ΔN, solo informativa) va de 3×10⁻³ (A·H = 1 000) a 3×10⁻² (A·H = 30–100). Es la misma escala que U37 con NPQT (3×10⁻³–3×10⁻²).
- **Los fallos de verdad** (RAMA_ROTA, ANISOTROPO, SIN_ATERRIZAJE y RECOLAPSO) aparecen solo con ε ≳ 10⁻²: 15 corridas en total con ε ≤ 10⁻², y la mayoría desde ε ≈ 3×10⁻².
- **La cizalla ayuda a la rama del cuscutón:** entra en W con signo positivo. Las roturas de rama vienen de la curvatura anisótropa (R3) cuando la deformación es grande.

### 3. Frente a lo que se esperaba (§8 del ESPEC)

- **«B0–B3 pasan»:** acertado.
- **«ε\* ~ 10⁻³–10⁻², como en U37»:** acertado a medias. Con la regla de |ΔN| ≤ 1, ε\* es 10⁻⁴–3×10⁻³, unas diez veces menor. Para que el hijo nazca, es la misma escala que en U37.
- **«Empeora con el tamaño de la madre»:** acertado, como decía la servilleta de NOTA_CUSCUTON_V4 §3.5.

### 4. Lo que no se hizo, declarado

- **El evento COLAPSO** (a_i < 10⁻³·a_b,iso) del §2 del ESPEC no se implementó. Solo se usó ANISOTROPO (máx q_i − mín q_i > 10), que corta antes cualquier colapso en una dirección.
- **La fila informativa con cizalla inicial con velocidad** (§7 del ESPEC) no se ejecutó.

### 5. Lo que significa para la completitud 4D del cuscutón

- **El problema clásico de los rebotes, la cizalla en la contracción, no rompe el nacimiento con cuscutón** para anisotropías de forma de hasta 10⁻²–3×10⁻² en la madre. La madre es cerrada, y la curvatura de Bianchi IX convierte la anisotropía en una oscilación de forma que crece solo ≈ ×8 hasta el rebote, no como a⁻⁶.
- **El precio es que el número de e-folds del hijo depende de la anisotropía de la madre** (ΔN ∝ ε²). Eso mueve la fracción de hijos con el rebote a la vista y la curvatura de cada uno, pero no la existencia de hijos ni el atractor de inflación.
- **Con esto, de la lista de lo que faltaba para la completitud 4D** (PLAN_V4_CUSCUTON) queda resuelta la anisotropía homogénea. Siguen abiertos:
  - la acción cúbica y el acoplamiento fuerte cerca de φ̇ → 0;
  - los vectores y los tensores;
  - los modos n = 0 y 1;
  - el pasado de la madre;
  - el recalentamiento.
