# U65 · Los tensores y los vectores a través del rebote con cuscutón · Preregistro

*2026-11-18 13:07 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con los vectores y tensores a través del rebote»). Es la parte C.4 del plan de la v4: lo que U59b y C4 del explorador dejaron solo en lo estructural. Escrito antes del código. Interno: no se exporta.*

## 1. Las preguntas

1. **Tensores:**
   - ¿Pasan las ondas gravitatorias el rebote con cuscutón sin problemas?
   - ¿Qué huella dejan en el espectro tensorial P_t de los hijos visibles, en ℓ = 2–40, frente al atractor?
2. **Vectores:** ¿hay grados de libertad vectoriales en esta acción? Y si la radiación tuviera vorticidad, ¿cuánto crecería la perturbación vectorial de la métrica en la contracción?

## 2. Lo que se sabe, y se declara

**Los tensores.**
- El cuscutón, χ y la radiación (el superfluido P(X) = X² de U59) son escalares, sin tensión anisótropa.
- Las ondas gravitatorias obedecen la ecuación de la relatividad general sobre el fondo: ḧ + 3Hḣ + n(n+2)h/a² = 0, con los armónicos tensoriales transversos y sin traza de S³ (n ≥ 2, k² + 2K = n(n+2)).
- La acción por polarización es (1/8)a³(ḣ² − n(n+2)h²/a²), positiva: no hay fantasmas ni inestabilidad de gradiente, también en el rebote.
- El cuscutón entra solo a través de a(t) y H(t). En la inflación, la acción tensorial de Boruah, Kim y Geshnizjani es la de la relatividad general (U58).

**Los vectores.**
- **En esta acción no hay grados de libertad vectoriales.** Todos los campos de materia son escalares y la radiación es irrotacional. La ligadura de momento vectorial fuerza que la perturbación vectorial invariante de la métrica sea nula para n(n+2) ≠ 2.
- **Para un plasma real con vorticidad,** la conservación del momento angular da (ρ + p)a⁴v_V = cte (v_V constante para la radiación), y la ligadura da un potencial vectorial de la métrica Φ_V ∝ 1/a², que crece en la contracción.

## 3. El montaje

**Tensores (con veredicto):**
- **Los hijos:** los cinco visibles de U61b, con el fondo de U62d (hacia delante desde el giro).
- **Los estados:** el vacío instantáneo del par canónico (h, π = a³ḣ/4), que no tiene la ambigüedad de marco de los escalares (ω = k/a), en tres sitios:
  - **A1:** 0.1 e-fold antes del giro;
  - **A2:** Δln a = 2 antes de la ventana;
  - **B:** Δln a = 2 después.
- **La propagación:**
  - sobre el fondo del hijo hasta ρ_r/ρ = 10⁻⁴;
  - después, sobre el fondo plano de U58 que arranca en ese punto (la etapa 2 de U62c), con la ecuación tensorial de U58 (z² = a², c = 1), hasta k/aH = 10⁻³.
- **El espectro:** P_t(n) = n³|h|²/π², la misma normalización que U58.mode con tensor=True.
- **La referencia:** el U58 tensorial con Bunch–Davies en el atractor, en el |χ_×| donde sale el modo.
- **Los cocientes:** F_t(n) = P_t/P_t,att, en 16 n de ℓ = 2–40 por hijo (como en U62c), más un control en n_c = 10·n_mem.
- **El número de cuantos (informativo):** N_t entre A2 y B, en n ∈ {100, 300, 10³, 3×10³, 10⁴}.

**Vectores (informativo, sin veredicto):**
- el factor de crecimiento (a_giro/a_rebote)² del potencial vectorial de un plasma con vorticidad, para los 52 anfitriones de U61b;
- la vorticidad inicial máxima, en el giro, que deja Φ_V(rebote) < 0.1.

## 4. Puertas (tensores; si falla alguna, U65 es INCONCLUSO)

- **G0, la normalización del vacío:** Ω(h, h*) = i a ≤ 10⁻¹⁰.
- **G1, el wronskiano al final:** a ≤ 10⁻³ relativo.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en el primer hijo y en n(ℓ = 2), n(ℓ = 10) y n(ℓ = 40). P_t,A2 cambia ≤ 10⁻³.
- **G3, el control:** |F_t(n_c) − 1| ≤ 0.03 con A2.
- **G4, el cambio de etapa:** con 10⁻⁵ en lugar de 10⁻⁴, P_t cambia ≤ 1 % en los tres n de G2.

## 5. Veredicto (tensores)

D_X = máx |F_t − 1| en ℓ = 2–40, para X = A1 y A2.
- **INCONCLUSO:** falla una puerta.
- **TENSORES_LISOS:** D_A1 y D_A2 ≤ 0.1 en los cinco hijos.
- **TENSORES_CON_RASGOS:** mín(D_A1, D_A2) ≥ 0.5 en al menos 3 hijos.
- **TENSORES_MODERADOS:** cualquier otro caso.

**Se informan siempre:**
- F_t(n) para A1, A2 y B, y el cociente de los rasgos tensoriales frente a los escalares de U62c/U62d en los mismos n;
- N_t;
- la tabla vectorial.

## 6. Lo que no cubre

- **La proyección a BB**, la polarización del CMB y r_eff en ℓ bajos.
- **Los tensores de n = 2** (el modo más bajo) fuera de la rejilla de ℓ.
- **La vorticidad de un plasma real,** solo como estimación de servilleta.
- **El acoplamiento tensor-escalar no lineal.**

## 7. Lo que se espera

- **G0–G4 pasan.**
- **Rasgos tensoriales en ℓ ≲ ℓ_mem** parecidos a los escalares (la misma amplificación en la contracción), y **TENSORES_CON_RASGOS o MODERADOS.**
- **Ninguna patología:** la acción tensorial es la de la relatividad general.
- **Vectores:** un crecimiento de (a_t/a_b)² ~ 10⁴–4×10⁵, lo que pediría una vorticidad inicial ≲ 10⁻⁶–10⁻⁵ en el giro.

## 8. Salidas

- **Código:** `chain/u65_tensores_vectores.py`.
- **Resultados:** `results/u65_tensores_vectores.json` y `results/u65_runs/`, con filas guardadas una a una.
- **Informe:** `docs/informes/INFORME_U65.md`.
- Interno.
