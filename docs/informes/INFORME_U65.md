# INFORME U65 · Los tensores y los vectores a través del rebote con cuscutón

**Preregistro:** `docs/preregistros/ESPEC_U65.md` (18 nov, 13:07).
**Código:** `chain/u65_tensores_vectores.py`. **Resultados:** `results/u65_tensores_vectores.json` y `results/u65_runs/`, con las filas una a una y `produce.out`.
**Prueba de funcionamiento, declarada:** `results/u65_runs/previo/prueba_humo.py`, que solo mira la normalización, el wronskiano y el tiempo.

Interno: no se exporta.

---

## 2026-11-18 13:18 (Madrid)

### 1. Veredicto preregistrado: TENSORES_CON_RASGOS

| puerta | resultado |
|---|---|
| G0, la normalización del vacío | **pasa:** 3×10⁻¹⁶ |
| G1, el wronskiano al final (≤ 10⁻³) | **pasa:** 4×10⁻⁸ |
| G2, la tolerancia | **pasa:** ≤ 1.2×10⁻⁸ |
| G3, el control n = 10·n_mem (≤ 3 %) | **pasa:** ≤ 0.6 % |
| G4, el cambio de etapa | **pasa:** ≤ 6×10⁻⁴ |

**F_t = P_t/P_t,atractor** en los cinco hijos visibles:

| hijo (A·H, f; ℓ_mem) | máx F_t en ℓ ≤ 5 | (escalar, U62c) | F_t en ℓ ≈ 30–40 | D_A1 / D_A2 | hijo que empieza tras el rebote (B) |
|---|---|---|---|---|---|
| 30, 0.096 (8.9) | 3.1 | (3.0) | 0.97–0.99 | 1.8 / 2.1 | 0.87–1.05 |
| 100, 0.184 (4.2) | 2.4 | (1.7) | 0.99–1.00 | 1.3 / 1.4 | 0.97–1.05 |
| 1 000, 0.542 (36) | 36 | (18) | 0.53–0.79 | 24 / 35 | 0.93–1.27 |
| 1 000, 0.835 (25) | 27 | (15) | 0.66–1.14 | 20 / 26 | 0.83–1.17 |
| 1 000, 0.898 (4.9) | 2.6 | (2.0) | 1.00–1.02 | 1.6 / 1.6 | 0.99–1.06 |

Los cinco tienen mín(D_A1, D_A2) ≥ 0.5, de ahí el veredicto.

**El número de cuantos tensoriales entre A2 y B,** N_t:
- A·H = 30: 0.95 (n = 100), 0.15 (300), 9×10⁻³ (10³), 1.3×10⁻⁴ (3×10³), 2×10⁻⁷ (10⁴);
- A·H = 1 000: 33 (n = 100), 2.2–2.5 (300), 0.33–0.41 (10³), 0.04 (3×10³), (0.4–3)×10⁻⁴ (10⁴).

### 2. Lo que significa

- **Los tensores cruzan el rebote sanos,** porque su acción es la de la relatividad general. No hay patología, solo la misma física de memoria que los escalares:
  - amplificación en la contracción de los modos de n ≲ n_mem;
  - rasgos en ℓ ≲ ℓ_mem;
  - nada a partir de unas pocas ℓ_mem;
  - nada en el hijo que empieza tras el rebote.
- **Los rasgos tensoriales son del orden de los escalares o algo mayores:** ×1–2 el cociente F_t/F_s en ℓ ≤ 5. Así, r_eff(ℓ) en los ℓ más bajos cambiaría en un factor de orden uno a dos.
- **No se observa:** con r = 0.005, el modo B en ℓ ≲ 5 está muy por debajo de la varianza cósmica y de los primeros planos. Es una consecuencia del modelo, no una prueba disponible.
- **Con A·H ≳ 100 no hay nada nuevo,** porque esos hijos ya están desfavorecidos por los escalares (U62e).

### 3. Los vectores (informativo)

- **En la acción de la ballena no hay grados de libertad vectoriales.** Toda la materia es escalar, y la radiación es el superfluido irrotacional de U59. La ligadura de momento vectorial anula el potencial vectorial invariante de la métrica.
- **Para un plasma real con vorticidad,** el potencial vectorial crecería en la contracción como (a_giro/a_rebote)²: entre **8×10³ y 3.7×10⁵** (mediana 3×10⁴) en los 52 anfitriones de U61b.
  - Para que Φ_V < 0.1 en el rebote, la vorticidad en el giro de la madre tiene que ser **≲ 3×10⁻⁷–10⁻⁵**, por debajo de ζ₀ = 4.6×10⁻⁵ en un factor de 4 a 170.
  - Es una condición sobre el estado inicial de la madre. Las perturbaciones escalares lineales no generan vorticidad, así que se añade al supuesto abierto de la semilla.
  - No es un fallo, pero hay que decirla.

### 4. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0–G4 pasan»:** acertado.
- **«Rasgos parecidos a los escalares; CON_RASGOS o MODERADOS»:** acertado, con CON_RASGOS.
- **«Ninguna patología»:** acertado.
- **«Vectores: (a_t/a_b)² ~ 10⁴–4×10⁵»:** acertado, con 8×10³–3.7×10⁵.

### 5. Para la v4

- **Lo que se puede decir:** «los tensores cruzan el rebote con su acción de relatividad general. En los hijos visibles, P_t tiene en ℓ ≲ ℓ_mem rasgos de orden uno, del tamaño de los escalares, inobservables con r = 0.005».
- **En los supuestos abiertos:** «en la acción no hay vectores; un plasma con vorticidad pediría una vorticidad inicial ≲ 10⁻⁶–10⁻⁵ en el giro de la madre».
