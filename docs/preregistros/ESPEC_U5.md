# U5 · ¿Sobrevive el sector oscuro si ξ sigue creciendo? · Preregistro

*8 de octubre de 2026, ≈16:30 (Madrid). Escrito antes de calcular nada de U5. Se conocen U4a (CIERRE_SEGÚN_RAMA) y U4b fase 1 (GATES_PASSED). La fase 2 de U4b (R_Γ en el atractor) está en marcha.*

## Motivo

U4b fase 1 mide, en el atractor, una pendiente dξ_U1/dℓ = 0.062, 0.181, 0.093 y 0.085 en [4.95, 5.65] (media 0.105), con ξ_U1(5.65) = 0.62. Las ramas de S6 suponen dos extremos en ℓ ≈ 50–85:
- **B2:** ξ constante = 1.2 (Hindmarsh et al.);
- **B1:** ξ = 0.24 ℓ (GHV, Buschmann et al.), excluida por Lyα en U3a y U4a.

Hay dos lecturas posibles de nuestra pendiente:
- **En valor absoluto,** 0.105 queda entre las dos.
- **En valor relativo,** d ln ξ/dℓ ≈ 0.105/0.62 ≈ 0.17, casi el 1/ℓ = 0.177 de un crecimiento ξ ∝ ℓ, pero con la mitad de normalización que GHV (ξ ≈ 0.11 ℓ).

Nuestro ξ a ℓ ≈ 5–6 parece, además, del orden de la mitad de los valores publicados. U5 decide si el sector oscuro sobrevive a esa extrapolación.

## U5a · rama B4, pendiente medida

**Definición.** Para todo ℓ de la emisión,
- ξ_B4(ℓ) = max(ξ_ref + s (ℓ − ℓ_ref), 0.3), con ξ_ref = 0.616 y ℓ_ref = 5.65 (U4B_RESULT);
- ε(ℓ) = R_Γ · ξ_B4(ℓ)/(0.24 ℓ). Es la misma convención de S6: emisión ∝ ξ y normalizada a la de GHV.

**Pendientes:** s ∈ {0.062, 0.105, 0.181}, es decir, el mínimo, la media y el máximo de G3′. Quedan fijadas aquí. Si U5d añade semillas, se informa aparte, pero estas filas no cambian.

**Dos normalizaciones:**
- **B4a (primaria):** la de nuestra red, como arriba. Es coherente con R_Γ, medido en la misma red y con el mismo recuento de cuerdas.
- **B4b (sistemática):** ξ_B4b = (1.2/ξ_ref) · ξ_B4a. Tiene el mismo crecimiento relativo, pero normalizado a ξ(ℓ_ref) = 1.2, el valor de B2. Cubre la posibilidad de que nuestro recuento subestime ξ por un factor ~2 frente a la literatura.

**R_Γ:** el de U4b fase 2 si pasa sus puertas (valor ± total); si no, el de U2 (1.69 ± 0.29). U5a se ejecuta después de la fase 2.

**Filas:**
- espectro log-plano (el IR ya está excluido con menos emisión);
- x₀ ∈ {10, 30}, K_sw ∈ {0, 1}, R ∈ {−1σ, central, +1σ};
- tres pendientes y dos normalizaciones: 72 filas.

Cada fila pasa por S6 (raíz de abundancia, P1–P4) y luego por la puerta Lyα de U4a, sin cambios: CLASS con el parche, 45 nodos, umbrales de U4a (estricto 0.234, amplio 0.336) y forma estricta como puerta primaria.

**Veredicto U5a**, con B4a y s = 0.105:
- **SOBREVIVE_A_LA_PENDIENTE** si alguna fila con R a +1σ está permitida;
- **DEPENDE_DEL_ESCALADO** si ninguna fila con R a −1σ está permitida. En ese caso, el sector oscuro solo vive si ξ deja de crecer entre ℓ ≈ 6 y 70, que es la disputa de la literatura;
- **MARGINAL** en cualquier otro caso.

Se informa también con s = 0.062 y 0.181 y con B4b. Se interpola δA en s para dar la pendiente crítica s_crit de la fila central, x₀ = 10, K_sw = 0.

## U5b · malla de forma ampliada (informativa)

- La familia de mezclas térmicas de U4a se amplía con m ∈ {0.2, 0.3, 0.4} keV y f ∈ {0.01, 0.02, 0.03}, en producto con la malla de U4a.
- Con ella se reajustan las 62 filas de U4a y las de U5a.
- Se informa cuántas filas cambian de veredicto de forma y dónde cae el mejor ajuste.
- El límite f ≤ 0.14 m^1.1 por debajo de 1 keV se marca como extrapolado.
- No cambia el veredicto preregistrado de U4a.

## U5c · isocurvatura de ruido blanco (cota)

- Estimación según Gorghetto, Trifinopoulos y Valogiannis (2025). El campo del Goldstone es aleatorio a escalas mayores que el horizonte cuando la masa entra, en H(T_*) = m_G, con k_* = a_* H_*.
- Por debajo de k_*, Δ²_iso(k) ≈ f_DM² (k/k_*)³ (cota de orden uno).
- Se compara con la potencia de materia adiabática en z = 3 y k = 20 h/Mpc, ambas con su crecimiento lineal (la isocurvatura, congelada hasta la igualdad).
- **Pasa** si el cociente es < 1 % para todas las filas permitidas de U4a y U5a.

## U5d · más semillas U(1) a 768³ (informativa)

- ξ_i = 0.5, semillas 63 y 64, tras la fase 2, con el código acelerado. Sirven para estrechar la pendiente.
- Se informa la pendiente media ± error con las cuatro semillas de ξ_i = 0.5 que haya.

## U5e · normalización de ξ frente a la literatura (informativa)

- Se compara nuestro ξ_U1(ℓ) en ℓ = 4–6 con los valores publicados: GHV 2021, Buschmann et al. 2022, Hindmarsh et al. 2021 y Kim y Saikawa 2024.
- Si nuestro ξ es sistemáticamente menor, la sistemática B4b pasa a ser la referencia en la discusión, no en el veredicto.

---

## Enmienda 1 · 8 de octubre de 2026, ≈21:05 (Madrid) · R_Γ del atractor y pendiente espectral

*Escrita tras la fase 2 de U4b y antes de calcular nada de U5.*

**1. R_Γ.** La fase 2 (RP³ y U(1) a 576³) da **R_Γ = 1.40 ± 0.13** en [4.36, 5.61]: ventana de U4b, bootstrap de semillas y sistemáticas de U2. Lo usan U5a y U5f.
- Además da R_ξ = 1.65 ± 0.10: la red RP³ tiene un 65 % más de longitud de cuerda que la U(1) en la misma ventana, y emite un 15 % menos por unidad de longitud.
- Frente a U2 (1.69 ± 0.29), |ΔR_Γ| = 0.296, por debajo del umbral de 0.3 de ESPEC_U4, así que U4a no se repite. Con menos emisión, las filas permitidas solo pueden aumentar.

**2. Corrección del análisis de U4b** (`network/u4b_analyze.py`, `inherit_gates`).
- ESPEC_U4 condiciona la fase 2 a las puertas de la fase 1 (N = 768). El código evaluaba cada N por separado y exigía, a 576³, unas redes con ξ_i = 0.25 y 2.0 que la fase 2 no produce. Por eso el bloque de 576 salía GATES_FAILED.
- Ahora ese bloque hereda las puertas del de 768, como dice el preregistro.
- Para información, la memoria de la condición inicial de las redes de 768 en la ventana de 576 es 19 % (ℓ = 4.36), 5 % (5.36) y 7.5 % (5.61).

**3. Nuevo U5f · pendiente espectral.**
- La pendiente del espectro de emisión de RP³ en ℓ ≈ 5.4 es q = 1.10 y 1.22 (las dos semillas); en ℓ ≈ 4.5–5.0 es 0.75–0.97.
- El espectro log-plano superviviente es q = 1. Con q > 1, el peso por intervalo logarítmico es ∝ x^(1−q), y sobre las ~60 e-foldings de x en ℓ ≈ 70 eso desplaza la energía hacia el infrarrojo. Un q apenas mayor que 1 cambia el resultado.
- **Filas:** espectro de potencia dE/dln x ∝ x^(1−q) en [x₀, m_r/H], el mismo soporte que el log-plano, con q ∈ {0.9, 1.1, 1.2}. Ramas B2 y B4a (s = 0.105); R ∈ {−1σ, central, +1σ}; x₀ ∈ {10, 30}; K_sw ∈ {0, 1}. Son 72 filas.
- **Referencia:** B2 con log-plano (q = 1) y el R_Γ nuevo, 12 filas. Son también las filas de U4a recalculadas con el R_Γ del atractor.
- **Puertas y umbrales:** los de U4a.
- **Veredicto U5f**, sobre B2:
  - **SOBREVIVE_A_q** si alguna fila con q = 1.2 y R central está permitida;
  - **DEPENDE_DE_q** si ninguna fila con q = 1.1 y R a −1σ está permitida (solo vive si q ≤ 1 en el ℓ físico);
  - **MARGINAL** en otro caso.
- Se da q_crit, interpolando δA en q, para la fila central con x₀ = 10 y K_sw = 0.

**4. U5d** se ejecuta con `u4b_produce.py 3` (U(1) 768³, ξ_i = 0.5, semillas 63 y 64).
