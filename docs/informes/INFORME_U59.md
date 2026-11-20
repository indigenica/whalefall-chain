# INFORME U59 · Los escalares a través del rebote con cuscutón: auditoría independiente

**Preregistro:** `docs/preregistros/ESPEC_U59.md` (15 nov 2026, 22:23, Madrid), con la enmienda 1 (cuadratura directa en lugar del ajuste en k²) y la enmienda 2 (G1 ya medida antes de la producción, y fallida, declarada).
**Código:** `chain/u59_escalares_rebote.py`. Deriva la acción cuadrática ADM con sympy desde la acción covariante, con un armónico explícito de S³. No usa el código del explorador.
**Resultados:** `results/u59_escalares_rebote.json`, más `results/u59_runs/` (la caché simbólica, `produce.out` y el barrido informativo `dense_n_adoptado.json`).
**Coste:** la producción tardó 78 s.

Interno: no se exporta.

---

## 2026-11-15 22:52 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO

| puerta | resultado |
|---|---|
| G0, la reproducción de C2 (rebote de C1) | **pasa:** autovalores de K_red a ≤ 3.3×10⁻¹¹ de los publicados (n = 2, 10 y 1 000) |
| G1, la convergencia de la cuadratura | **falla**, ya declarado en la enmienda 2: 7.7×10⁻¹⁵ en n = 37, pero 2.0×10⁻¹⁰ en n = 3×10⁵ frente a ≤ 10⁻¹⁰. Es el suelo de redondeo |
| G2, la convergencia de la transferencia (rtol 10⁻⁹ frente a 10⁻¹⁰) | **pasa:** ≤ 1.9×10⁻⁹ |
| G3, la simplecticidad | **pasa:** ≤ 10⁻²⁵ |
| G4, la energía definida positiva en los extremos | **pasa:** en todos los n |

**Aunque G1 hubiera pasado, el veredicto sería también INCONCLUSO.** Con el juego adoptado:
- G(3×10⁵)/G(10⁴) = 2.68, por encima de 2, así que no se cumple la regla de ESCALARES_SANOS;
- G no crece de forma monótona en los tres últimos n (1.97, 3.75 y 3.47), así que tampoco se cumple la de INESTABILIDAD.

### 2. Lo que queda establecido por la auditoría

**La acción cuadrática del explorador queda confirmada de forma independiente.** Derivada aquí desde la acción covariante, sin su código:
- **los autovalores cinéticos de C2** en el rebote de C1 coinciden a 10⁻¹¹;
- **el coeficiente de ζḟ** es −3μ²a³ (a 10⁻¹³), y ḟ no aparece en ningún otro término, como dice C2;
- **la ganancia de energía de C3e**, en φ ∈ [−0.2, 0.2] con el juego de C1, se reproduce en sus cuatro cifras. Con el cuadrado de la ganancia en amplitud que da este código:

| n | ganancia de energía (U59, G²) | C3e (publicado) |
|---|---|---|
| 2 | 12.303 | 12.3027 |
| 1 000 | 2.555 | 2.5548 |
| 10 000 | 1.711 | 1.7114 |
| 30 000 | 1.603 | 1.6026 |

**Nota sobre la definición:** el ESPEC define G(n) como σ_max de E_f^(1/2)·T·E_i^(−1/2), que es una ganancia **de amplitud**. La de energía, la que publica C3e, es G². Las reglas de §5 se aplicaron a G tal como estaba escrito.

### 3. El juego adoptado (μ = 0.16, m = 0.02, φ∞ = 0.25), informativo

- **La madre:** A·H = 100 y f = 0.294 (el hijo de f más bajo con A·H = 100 en U57b). Rebota en a_b = 80 397; la ventana de φ = −0.2 a +0.2 dura 46 unidades de tiempo de Planck.
- **Sin fantasmas:** en 81 instantes × 9 armónicos, los dos autovalores de K_red (gauge de radiación uniforme) son positivos. El cociente mínimo es 2.3×10⁻⁶, y la matriz auxiliar es siempre invertible (|det| relativo ≥ 10⁻⁸).
- **Las ganancias de amplitud:**

| n | k_fís en el rebote (M_Pl) | G (amplitud) | G² (energía) |
|---|---|---|---|
| 2 | 3.5×10⁻⁵ | 2.49 | 6.2 |
| 100 | 1.3×10⁻³ | 2.39 | 5.7 |
| 1 000 | 0.012 | 1.66 | 2.7 |
| 10⁴ | 0.12 | 1.29 | 1.7 |
| 3×10⁴ | 0.37 | 1.97 | 3.9 |
| 10⁵ | 1.24 | 3.75 | 14.1 |
| 3×10⁵ | 3.73 | 3.47 | 12.0 |

- **El barrido denso** (posterior, sin veredicto; 16 valores de n entre 10⁴ y 3×10⁵, en `dense_n_adoptado.json`): **G oscila entre 1.3 y 3.75 sin tendencia a crecer.**
  - El máximo, 3.75 (energía 14), aparece ya en n ≈ 1.3×10⁴, con k_fís ≈ 0.16 M_Pl, y se repite en varios n.
  - En n = 3.2×10⁵ G vale 1.35.
  - **La regla de los dos puntos** (3×10⁵ frente a 10⁴) **cayó en un valle de la oscilación en 10⁴.** No era un buen estimador de tendencia.
- **Lectura, sin probar:** con el juego adoptado, el rebote en φ es unas 40 veces más brusco que con C1 (D(0) = 1.8×10⁻⁴). Los modos hasta k_fís ~ 0.1–4 M_Pl lo cruzan sin adiabaticidad, con una mezcla finita (energía ≤ 14), que depende de la fase y no crece con n. Es la firma de un transitorio, no de una inestabilidad de gradiente. Pero el máximo cae ya por debajo de la escala de Planck.

### 4. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0 la menos segura»: falló la previsión, en el buen sentido.** G0 pasa a 10⁻¹¹, y la derivación independiente coincide con la del explorador.
- **«G0–G4 pasan»: fallado en G1,** por un umbral por debajo del redondeo.
- **«ESCALARES_SANOS»: no se alcanza**, por la regla de los dos puntos. El barrido denso apunta a una mezcla acotada.

### 5. Siguiente, como propuesta

- **U59b, con datos nuevos** (otra madre: A·H = 300 o la de f más alto) y tres cambios fijados de antemano:
  1. G1 desde el suelo de redondeo, por ejemplo ≤ 10⁻⁹, o sumas compensadas;
  2. la regla de crecimiento sobre un barrido denso: la pendiente de log G frente a log n en n ∈ [10⁴, 3×10⁵], o el máximo de G en décadas sucesivas, en lugar de dos puntos;
  3. la ganancia de energía G², como en C3e.
- **Para la v4,** lo que ya se puede decir sin veredicto: la acción cuadrática del rebote con cuscutón está derivada de dos maneras independientes que coinciden, y no tiene fantasmas en el rango calculado.

---

## U59b · 2026-11-15 23:57 (Madrid): la repetición con otra madre y reglas mejor elegidas

**Preregistro:** `docs/preregistros/ESPEC_U59b.md`. Las reglas se eligieron después de ver U59, y así se declara.
**Código:** `chain/u59b_escalares_rebote.py`, que importa U59 y U58 sin cambios.
**Resultados:** `results/u59b_escalares_rebote.json` y `results/u59b_runs/produce.out`.

### 1. Veredicto preregistrado: ESCALARES_SANOS

| puerta | resultado |
|---|---|
| G0, la reproducción de C2 | **pasa:** ≤ 3.3×10⁻¹¹ |
| G1, la cuadratura (≤ 10⁻⁸) | **pasa:** 7.7×10⁻¹⁵ con n = 37 y 2.0×10⁻¹⁰ con n = 3×10⁵ |
| G2, la convergencia en rtol (n = 10², 10⁴ y 3×10⁵) | **pasa:** 1.4×10⁻¹¹, 6.8×10⁻¹⁰ y 3.2×10⁻⁹ |
| G3, la simplecticidad | **pasa** |
| G4, la energía definida positiva en los extremos | **pasa** |

**Una corrección de código, declarada:** n = 10⁴ no estaba en la rejilla densa (los valores redondeados de la rejilla no lo contienen), así que el script evaluó G2 solo en dos de los tres n. El tercero se calculó después con el mismo código y está añadido al JSON con una nota.

**La madre:** el juego adoptado, A·H = 300, f = 0.604, 92 e-folds y a_b = 1.31×10⁵.

**Sin fantasmas:** el cociente mínimo de los autovalores de K_red es 6.3×10⁻⁶ > 0, y la matriz auxiliar es siempre invertible.

**La ganancia de energía por décadas** (28 valores de n):

| década de n | máximo de E |
|---|---|
| 10²–10³ | 5.8 |
| 10³–10⁴ | 13.7 |
| 10⁴–10⁵ | 13.3 |
| 10⁵–3×10⁵ | 13.2 |

No crece: M₄/M₃ = 0.97 y M₅/M₄ = 1.00. Dentro de cada década, E oscila entre ≈ 1 y ≈ 13.7, según la fase con la que cada modo cruza el rebote. Es el mismo comportamiento que el barrido informativo de U59 con A·H = 100 (máximo 14.1).

### 2. Lo que queda de la prueba C (informativo)

**C.3, los modos n = 0 y n = 1,** con la acción de esta auditoría, en el rebote de C1 y frente a C4 del explorador:
- **los autovalores cinéticos** de n = 0, 1 y 2 coinciden a ≤ 10⁻¹⁵:
  - n = 0: 0, 0.7266 y 8.58;
  - n = 1: 0, 0 y 0.7266, con la segunda dirección nula que C4 atribuye a que en S³ la parte sin traza de D_iD_jY se anula;
  - n = 2: 0, 0.7266 y 2.6095;
- **los de la matriz auxiliar** coinciden al reescalar β → a²β. Aquí N_i = D_i(βY); en C4, N^i = D^iβ. Con ese cambio de convención coinciden a 10⁻¹⁵.
- **Con el juego adoptado** la estructura es la misma: en n = 1 hay dos direcciones nulas y una positiva (0.583).
- **La lectura física de esas direcciones nulas** (gauge o ligadura global de S³) sigue abierta, como en C4. Haría falta un análisis de Dirac adaptado a n = 0 y n = 1.

**C.5, el comienzo de la inflación,** con la acción de Boruah, Kim y Geshnizjani en el fondo plano del juego adoptado (los primeros 60 e-folds, 200 instantes × 25 escalas con k/aH ∈ [10⁻³, 10³]):
- **z² > 0 en todo el rango:** z²/(a²α) ≥ 0.99999988;
- **c_s² ≥ 0.9948** en todo el rango.

No hay inestabilidad de gradiente ni fantasma en la inflación.

**C.4, los vectores y los tensores:** no se calculan; están en el §4 del ESPEC.

### 3. Lo que esto deja para la v4

**La parte escalar lineal del rebote con cuscutón está calculada y es sana:**
- con el juego adoptado: sin fantasmas, con una ganancia de energía acotada (≤ 14) que no crece con el armónico hasta n = 3×10⁵, y regular al empezar la inflación;
- con una acción cuadrática derivada de dos maneras independientes que coinciden.

**Lo que no se puede decir todavía:** nada de los modos n = 0 y 1 (su interpretación), de los tensores y vectores fuera de lo estructural, de la acción cúbica y el acoplamiento fuerte, ni de la escala de corte de la EFT. El máximo de ganancia aparece ya a k_fís ≈ 0.05 M_Pl.

---

## Nota · 2026-11-16 14:03 (Madrid): la ganancia de energía que no decae es un efecto del marco (U62 y U62b)

**La «energía» E** de U59, U59b y C3e del explorador es el hamiltoniano congelado de la acción reducida con el término q̇ᵀB_s q tal cual (B_s, la parte simétrica de B).
- Integrando ese término por partes (una transformación canónica dependiente del tiempo), la dinámica es la misma, pero el hamiltoniano congelado cambia.
- En ese otro marco, lejos del rebote, las frecuencias son las locales: k para χ y k/√3 para la radiación. En el de E salían ≈ 1.9k y 1.0k.
- **Medida en el marco físico** (U62b, con los estados a Δln a = 2 del rebote), la mezcla cae rápidamente con n: N ≲ 10⁻⁵ para n ≳ 10⁴. En U62, con E, N era de orden uno.

**Por tanto:**
- **El veredicto ESCALARES_SANOS de U59b no cambia:** sin fantasmas y sin crecimiento con n.
- **La frase «una mezcla finita (energía ≤ 14) que no crece con n» se mantiene como afirmación sobre la norma E,** que es correcta pero no física.
- **«No decae con n» no describe la física.** La física es que la memoria del rebote se borra por encima de unas pocas n_mem (`INFORME_U62.md`, sección U62b).

---

## Nota · 2026-11-17 14:29 (Madrid): el control de G2 en n = 10⁴, ya en el código

La revisión de Codex en la PR #89 señaló que `u59b_escalares_rebote.py` filtraba los n de control de G2 que no están en la rejilla densa. Por eso, al repetir el script se perdía el control de n = 10⁴, que en el JSON había añadido a mano (§1 de U59b).
- **El código ya lo calcula:** para un n de control fuera de la rejilla, hace una referencia propia con rtol 10⁻¹⁰.
- **Se repitió U59b entera.** El resultado es idéntico:
  - G2 = 1.4×10⁻¹¹, 6.8×10⁻¹⁰ y 3.2×10⁻⁹ (los mismos valores);
  - los máximos por década no cambian;
  - ESCALARES_SANOS.
- **El JSON sale ahora del código,** sin añadidos a mano.
