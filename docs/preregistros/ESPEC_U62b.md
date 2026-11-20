# U62b · Los modos físicos a través del rebote con cuscutón: Bogoliubov en el marco que recupera la física local · Preregistro

*2026-11-16 13:40 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U62b»). Escrito antes del código de U62b. Interno: no se exporta.*

## 1. Por qué una U62b distinta de la que proponía INFORME_U62

**U62 midió el número de partículas contra un vacío instantáneo mal elegido.**
- La acción reducida en el gauge de radiación uniforme, a³(½q̇ᵀKq̇ + q̇ᵀBq + ½qᵀCq), tiene una B con parte simétrica B_s.
- El hamiltoniano congelado depende de si q̇ᵀB_s q se deja así (U62) o se integra por partes: ½(a³ qᵀB_s q)˙ − ½qᵀ(a³B_s)˙q.
- Los dos marcos están relacionados por una transformación canónica dependiente del tiempo, p′ = p − a³B_s q. Tienen la misma dinámica, pero vacíos instantáneos distintos.

**Las comprobaciones previas a este ESPEC** (declaradas: frecuencias congeladas y adiabaticidad, sin calcular ningún Bogoliubov; scripts en `results/u62b_runs/previo/`, que se ejecutan con `PYTHONPATH=chain`):
- **en el marco de U62,** con φ = +0.2 y n ≥ 10⁴, las frecuencias son ≈ 1.90k y 1.00k;
- **en el marco sin B_s,** son 0.578k y 1.000k: la radiación (1/√3) y χ (1), como pide la física local;
- **el operador de segundo orden del explorador** (C3c, `whalefall-birth-explorer` PR #11) da lo mismo: 0.578 / 1.000.
- **El «modo rápido» de U62 es, por tanto, un artefacto del marco.** C3c ya advertía que las frecuencias congeladas de este gauge no son velocidades físicas cerca del rebote, y que los límites H → 0 y n → ∞ no son uniformes.
- **En el marco sin B_s, la energía congelada E′ sigue positiva y es adiabática** (|ω̇|/ω² ≲ 0.1) durante Δln a ≈ 1–4 después de la ventana, para n ≳ 10³ en los tres hijos mirados. Para n ≈ 100, nunca: allí k/aH ≲ 8 siempre.

Por eso U62b pone los estados en ese marco, y lejos del rebote, en la contracción y en la expansión, donde los modos ya son los locales.

## 2. El montaje

- **Hijos:** los cinco visibles de U61b, como en U62.
- **Fondo:**
  - la ventana φ ∈ [−0.2, 0.2] de U59.Trajectory;
  - un tramo de contracción integrado hacia atrás en el modo −1 de `explora_cuscuton`, desde φ = −0.2 hasta Δln a = +2 (la contracción dura 4.5–6.4 e-folds);
  - un tramo de expansión en el modo +1, desde φ = +0.2 hasta Δln a = +2.
  - Sin proyección de la ligadura en tramos tan cortos; su residuo se informa.
- **La evolución:** la misma de U59, en (q, p) con las matrices K, B y C, sobre todo el intervalo.
- **Los estados:**
  - en cada extremo, la transformación canónica p′ = p − a³B_s q, el reescalado (a^{3/2}q, a^{−3/2}p′) y los modos normales de E′ = energy_matrix(K, B_a, C − Ḃ_s − 3HB_s);
  - Ḃ_s por diferencias centradas con paso relativo 10⁻⁵ en t;
  - vacío instantáneo a la entrada y descomposición a la salida, con las fórmulas de U62.
- **Las etiquetas:** el modo de salida de frecuencia más próxima a k es **χ**, y el más próximo a k/√3, la **radiación**.
- **Los modos:** n ∈ {100, 200, 500, 10³, 2×10³, 5×10³, 10⁴, 3×10⁴, 10⁵}, más los n_obs de ℓ = 2, 10 y 30 de cada hijo, como en U62.
- **Filas definidas.** Una fila (hijo, n) está **definida** si se cumplen las tres condiciones siguientes. Las demás son NO_DEFINIDAS: se informan y no cuentan.
  1. E′ > 0 en los dos extremos;
  2. adiabaticidad |ω̇|/ω² ≤ 0.1 en los dos extremos y en los dos modos, con ω̇ por diferencias con paso relativo 10⁻³;
  3. con los extremos en Δln a = 1.5 en lugar de 2, N_χ cambia ≤ máx(0.005, 0.25·N_χ).

## 3. Puertas (si falla alguna, U62b es INCONCLUSO)

- **G0, la normalización:** como en U62, ≤ 10⁻¹⁰.
- **G1, la identidad de Bogoliubov:** ≤ 10⁻⁶.
- **G2, la convergencia** (rtol 10⁻⁹ frente a 10⁻¹⁰, en n = 10³, 10⁴ y 10⁵ del primer hijo): N_χ y N_rad cambian ≤ 10⁻⁴ + 10⁻⁴·N.
- **G3, las etiquetas:** en las filas definidas, ω_χ/k ∈ [0.95, 1.05] y ω_rad/k ∈ [0.95, 1.05]/√3 en los dos extremos.
- **G4, el fondo:** el residuo relativo de Friedmann en los tramos nuevos es ≤ 10⁻⁶.

## 4. Veredicto

Sobre N_χ en las filas definidas de los n_obs (5 hijos × 3 ℓ):
- **INCONCLUSO:** falla una puerta, o hay menos de 3 filas definidas.
- **REBOTE_ADIABATICO:** N_χ ≤ 0.01 en todas las filas definidas (rasgos de potencia ≤ ±20 %).
- **REBOTE_EXCITA:** N_χ ≥ 0.1 en al menos un tercio de las filas definidas.
- **INTERMEDIO:** cualquier otro caso.

**Se informan siempre:**
- N_χ y N_rad frente a n en las filas definidas y no definidas del barrido;
- R_min y R_max para χ;
- la comparación con U62 (el marco de U62).

## 5. Lo que no cubre

- **El espectro P_ζ.** Para los n que nunca son adiabáticos (ℓ = 2 en casi todos los hijos), el número de partículas no está definido, y solo un cálculo del espectro hasta la salida del horizonte en la inflación (U62c) diría algo.
- La conversión radiación → ζ.
- La acción cúbica y el corte de la EFT.
- Otros estados de la madre: su pasado está abierto.

## 6. Lo que se espera

- **G0–G4 pasan.**
- **Entre 3 y 8 filas definidas,** las de ℓ = 30 y algunas de ℓ = 10.
- **REBOTE_ADIABATICO,** con confianza media: en U62, el modo casi puro de δχ tenía N ≤ 0.004 en las filas limpias, y el artefacto del marco no debería afectarle mucho.

## 7. Salidas

- **Código:** `chain/u62b_bogoliubov_fisico.py`.
- **Resultados:** `results/u62b_bogoliubov_fisico.json` y `results/u62b_runs/`.
- **Informe:** sección U62b de `docs/informes/INFORME_U62.md`.
- Interno.
