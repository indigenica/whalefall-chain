# INFORME U64 · El orden cúbico completo en el gauge unitario del cuscutón

**Preregistro:** `docs/preregistros/ESPEC_U64.md` (18 nov, 09:44).
**Código:** `chain/u64_cubico_unitario.py`. **Resultados:** `results/u64_cubico_unitario.json` y `results/u64_runs/`: las filas una a una y `produce.out`.
**Una corrección de código, declarada:** el primer intento se cayó en G1 al pedir n = 10 del hijo A·H = 1 000, f = 0.542, cuyo vacío en A2 no está definido, como ya se sabía por U63b. G1 salta ahora esas filas. Solo se había calculado G0 (`produce_intento1.out`).

Interno: no se exporta.

---

## 2026-11-18 10:23 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO (fallan G1 y G3, por poco)

| puerta | resultado |
|---|---|
| G0, el L₂ frente a U59 | **pasa:** 1.2×10⁻¹⁵ |
| G1, las ligaduras lineales (α, β) en el gauge unitario (≤ 10⁻⁶) | **falla:** 5.6×10⁻⁶. El paso de gauge es correcto, signos incluidos, a 6×10⁻⁶ de la suma de \|contribuciones\| |
| G2, la tolerancia (≤ 10⁻³) | **pasa:** ≤ 7×10⁻⁹ |
| G3, 128 frente a 64 puntos (≤ 10⁻³) y el doble de instantes (≤ 5 %) | **falla en los puntos:** 2.7×10⁻³. Pasa en los instantes: ≤ 5×10⁻⁴ |

**Con las puertas en verde, el veredicto habría sido LIMITE.**
- En ningún hijo llega ρ₃ a 1 con n ≤ 10⁴.
- Pero en el hijo A·H = 100, f = 0.965 (no visible), ρ₃ = 0.21 en un n ≤ n(ℓ = 30), por encima del 0.1 que pedía PERTURBATIVO.

### 2. Lo que sale

**El máximo de ρ₃ = n^{3/2}·⟨|L₃|⟩/⟨|L₂|⟩** (lagrangiano ADM completo, gauge unitario), sobre todo el camino:

| hijo (A·H, f) | n = 100 | 10³ | 3×10³ | 10⁴ | n_sc (U63b) |
|---|---|---|---|---|---|
| 30, 0.096 | 6×10⁻⁴ | 8×10⁻³ | 0.049 | **0.53** | 3×10³ |
| 48, 0.104 | 4×10⁻⁴ | 5×10⁻³ | 0.032 | 0.33 | 3×10³ |
| 100, 0.184 | 1×10⁻³ | 9×10⁻³ | 0.017 | 0.17 | 3×10³ |
| 100, 0.965 | 3×10⁻⁴ | 8×10⁻³ | 0.020 | 0.21 | 10⁴ |
| 300, 0.246 | 3×10⁻⁴ | 2×10⁻³ | 0.008 | 0.058 | 10⁴ |
| 1 000, 0.542 | 1×10⁻⁴ | 1×10⁻³ | 0.004 | 0.021 | 10⁴ |
| 1 000, 0.835 | 1×10⁻⁴ | 1×10⁻³ | 0.008 | 0.022 | 10⁴ |
| 1 000, 0.898 | 1×10⁻⁴ | 1×10⁻³ | 0.005 | 0.023 | 10⁴ |

- **ρ₃ no llega a 1 en ningún hijo con n ≤ 10⁴.** Crece como ~n² entre 3×10³ y 10⁴. Extrapolando, llegaría a 1 hacia n ≈ 1.4×10⁴ (A·H = 30) y 7×10⁴ (A·H = 1 000), es decir, **k ≈ 0.3 M_Pl en el rebote** en todos.
- **El lapso en el gauge unitario es pequeño:** la banda de |α_u| ≤ 9×10⁻⁴ en todos los n y fases. Esto corrige la lectura de U63b (§3).
- **El máximo está en la ventana del rebote** para n ≥ 10³, y en la expansión o la contracción para n bajos.
- **Qué sector domina L₃:**
  - para n ≥ 10³, la radiación (85–100 %): el fluido P(X) = X², cuyo potencial de velocidades en el gauge unitario es w_u = −q·ξ⁰;
  - **el cuscutón aporta ≲ 0.1–0.34** en n ≲ 10³ y nada en n alto;
  - para n bajos, la gravedad y χ se reparten.
- **Donde φ̇₀ → 0** (el giro de la madre y la inflación), nada: todos los máximos con n ≲ 300 están por debajo de 10⁻³.

### 3. Lo que esto corrige de U63b (lectura)

- **M3 de U63b** (|ξ⁰|·|φ̈₀/φ̇₀|, que pasaba de 1 con n ≈ 3×10³–10⁴) **es solo uno de los términos del lapso en el gauge unitario,** α_u = α − ḟ/φ̇₀ + fφ̈₀/φ̇₀².
  - Sobre la solución, esos términos se cancelan casi del todo: |α_u| ≲ 10⁻³ cuando M3 ~ 1. La cancelación es de unos 3 órdenes, sin problema en doble precisión, y G1 confirma que los campos transformados cumplen las ligaduras.
  - **La escala fuerte que daba U63b (k ≈ 0.04–0.12 M_Pl) era una sobrestimación de la no linealidad.**
- **Con el cúbico completo, la teoría es perturbativa al menos hasta n = 10⁴,** k ≈ 0.05–0.23 M_Pl en el rebote. La extrapolación pone el límite hacia k ~ 0.3 M_Pl.
  - **Lo fija el fluido de radiación, no el cuscutón.** La radiación descrita como superfluido P(X) es a su vez una teoría efectiva, y un plasma real no tiene ese límite en la misma forma.
  - Es la misma escala que la de la gravedad cuántica (k_fís ~ M_Pl), muy por encima de las del rebote: H_b = 2×10⁻³, k_mem ≈ 10⁻², DGQ ≈ 0.016.
- **Hay que anotar en U63b** que su escala fuerte era una cota pesimista. La lectura buena es la de U64.

### 4. Para la v4

- **Lo que se puede decir:** «el orden cúbico completo (gravedad, inflatón, radiación y cuscutón), evaluado en el gauge unitario del cuscutón sobre los modos de vacío, se mantiene por debajo del cuadrático para todos los modos calculados (k ≲ 0.05–0.2 M_Pl en el rebote). La extrapolación pone el acoplamiento fuerte hacia k ~ 0.3 M_Pl, dominado por el fluido de radiación. No hay acoplamiento fuerte donde φ̇ → 0.»
- **Lo que queda abierto:**
  - las tríadas de modos distintos (el C9 del explorador, con la exportación de U63b);
  - n = 3×10⁴, por la precisión;
  - los bucles;
  - el gauge plano de la radiación como comprobación cruzada.
