# U59b · Los escalares a través del rebote con cuscutón, repetido con datos nuevos y reglas mejor elegidas · Preregistro

*2026-11-15 23:54 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U59b y lo que queda de la prueba C»). Escrito después de ver U59 (INCONCLUSO) y antes de cualquier cálculo de U59b. Interno: no se exporta.*

## 1. Qué cambia respecto a U59, declarado (las reglas se eligen después de ver U59)

1. **Datos nuevos: otra madre.** El juego adoptado (μ = 0.16, m = 0.02, φ∞ = 0.25) con A·H = 300 y el hijo de f más bajo de U57b (f = 0.604, 92 e-folds, a_b ≈ 1.3×10⁵). U59 usó A·H = 100.
2. **G1 desde el suelo de redondeo:** con el doble de puntos, las integrales cambian ≤ **10⁻⁸** relativo en n = 37 y en n = 3×10⁵. En U59 el redondeo medido fue 2×10⁻¹⁰.
3. **La ganancia es la de energía, E(n) = G(n)²,** como en C3e.
4. **El crecimiento se juzga en un barrido denso,** no en dos puntos:
   - n log-uniforme con 8 valores por década, de 10² a 3×10⁵ (28 valores);
   - M_d es el máximo de E en la década [10^d, 10^(d+1)), con d = 2, 3, 4, y en [10⁵, 3×10⁵] para d = 5.

## 2. El resto, como en U59

- La acción (`chain/u59_escalares_rebote.py`, importado sin cambios).
- El gauge de radiación uniforme y la ventana φ ∈ [−0.2, 0.2].
- La positividad en 81 instantes × 9 armónicos (los de U59).
- G0, la reproducción de C2 a ≤ 10⁻⁸.
- G2 (convergencia en rtol, ≤ 10⁻⁴, en n = 10², 10⁴ y 3×10⁵), G3 (simplecticidad, ≤ 10⁻⁶) y G4 (energía definida positiva en los extremos, en todos los n).

## 3. Veredicto

- **INCONCLUSO:** falla una puerta.
- **ESCALARES_SANOS:**
  - los dos autovalores de K_red son > 0 y S es invertible en todo el muestreo;
  - y la energía no crece más de un factor 2 por década: M₄ ≤ 2·M₃ y M₅ ≤ 2·M₄.
- **INESTABILIDAD:**
  - algún autovalor de K_red < −10⁻⁹ relativo;
  - o M₃ < M₄ < M₅ con M₅ > 10·M₃ (crecimiento sostenido de un orden de magnitud en dos décadas).
- **INCONCLUSO** también si no se cumple ninguna de las dos.

**Se informan:** E(n) en todos los n, el k_fís en el rebote y, como contraste, el barrido denso de U59 (A·H = 100).

## 4. Partes informativas, sin veredicto: lo que queda de la prueba C

- **C.3, n = 0 y n = 1.** Con la misma acción, en el rebote de C1:
  - el rango y los autovalores de la matriz auxiliar;
  - los autovalores cinéticos tras eliminar las auxiliares (en n = 0 sin β, que no existe: Y′ = 0);
  - comparados con lo publicado por el explorador en C4 (`revenge/C4/RESULTS_C4.json`).

  Y lo mismo en el rebote del juego adoptado.
- **C.5, el comienzo de la inflación,** con la acción de Boruah, Kim y Geshnizjani (U58), en el fondo plano del juego adoptado durante los primeros 60 e-folds de la inflación:
  - el signo de z² y de c_s² para k/(aH) ∈ [10⁻³, 10³] (25 valores) en 200 instantes;
  - el valor mínimo de c_s².

  Es la versión de C.5 con la acción ya reducida. El gauge de radiación uniforme degenera ahí (U56), así que no se usa.
- **C.4, los vectores y los tensores:** no se calculan. La radiación es irrotacional por construcción, y los tensores son los de Einstein (el explorador, C4, de forma estructural; no auditado).

## 5. Lo que se espera

- **G0–G4 pasan.**
- **ESCALARES_SANOS,** con confianza media: el barrido denso de U59 mostraba una ganancia acotada.
- **C.3:** reproduce la tabla de C4.
- **C.5:** z² > 0 y c_s² > 0 en todo el rango.

## 6. Salidas

- **Código:** `chain/u59b_escalares_rebote.py`.
- **Resultados:** `results/u59b_escalares_rebote.json`.
- **Informe:** una sección de `INFORME_U59.md`.
