# U59 · Los escalares a través del rebote con cuscutón: auditoría independiente · Preregistro

*2026-11-15 22:23 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U59»). Es la prueba C de `paper/zenodo_whale/PLAN_V4_CUSCUTON.md`, en su versión central (C.1, C.2 y la positividad). Escrito antes del código de U59. Interno: no se exporta.*

## 1. La pregunta

Con Einstein + cuscutón (el juego adoptado: μ = 0.16, m = 0.02, φ∞ = 0.25) + χ + radiación, en S³:
- ¿es sana la acción cuadrática de los escalares a través del rebote? Es decir: energía cinética positiva tras eliminar las ligaduras y matriz auxiliar invertible;
- ¿crece la transferencia de energía con el número de armónico, que sería la firma de una inestabilidad de gradiente?

## 2. Independencia y lo que ya se sabe (declarado)

- **La acción cuadrática se deriva aquí sin el código del explorador** (`revenge/C2`, `C3d`, `C3e`), del que solo se leen la convención de variables y los **números publicados**:
  - h_ij = a²e^(2ζ)γ_ij, N = 1 + α, N^i = D^iβ, (δχ, δθ, δφ) = (x, w, f);
  - la radiación como P(X_θ) = X_θ²;
  - el término −3μ²a³ζḟ integrado por partes para que f sea auxiliar;
  - la normalización ∫Y²√γ = 1 y L₂ por volumen a³.
- **El método:** el lagrangiano ADM se escribe con un armónico explícito de S³, Y(r) = sen((n + 1)r)/sen r (el de ℓ = 0 alrededor de un punto), y se desarrolla a segundo orden con sympy, punto a punto. Se integra sobre S³ con cuadratura de Gauss y se identifica la dependencia en k² = n(n + 2) como un polinomio, comprobado en un n que no entra en el ajuste. Como el fondo es homogéneo, el resultado vale para todos los armónicos con el mismo n.
- **Lo publicado por el explorador para el juego de C1:**
  - C2: los autovalores cinéticos en el gauge de radiación uniforme (bloque (ζ, x)) en el rebote: 0.71759660 y 2.60952408 (n = 2), 3.47931201 (n = 10) y 11.56368148 (n = 1 000);
  - C3e: las ganancias de energía en φ ∈ [−0.2, 0.2], sin crecimiento con n hasta 3×10⁴.
- **Yo ya vi** el diagnóstico de `chain/explora_cuscuton_modos.py` (con las matrices del explorador, juego de C1), sin crecimiento con n hasta 3×10⁵. **Con el juego adoptado no se ha calculado nada.**

## 3. El montaje

- **El fondo:** el de `chain/explora_cuscuton.py` (sin cambios), en el reloj φ dentro de la rama central |φ| ≤ φ∞/2.
  - Para G0, el juego de C1 con los datos de C1 (`c1_row`).
  - Para la producción, el juego adoptado, con la madre A·H = 100 y el primer f del barrido de U57b que da HIJO (el de f más bajo).
- **Las variables:** q = (ζ, x, w) y las auxiliares u = (α, β, f). Las matrices K₀, B₀, C₀, T, R y S se reducen eliminando u. El gauge temporal se fija con w = 0 (radiación uniforme), y queda q = (ζ, x).
- **La positividad:** en 81 instantes con φ ∈ [−0.2, 0.2] y n ∈ {2, 3, 10, 100, 10³, 10⁴, 3×10⁴, 10⁵, 3×10⁵}:
  - los dos autovalores de K_red tienen que ser > 0;
  - la matriz auxiliar tiene que ser invertible (|det| relativo > 10⁻¹²).
- **La transferencia:**
  - se integran las ecuaciones de Euler–Lagrange del sistema reducido, d/dt[a³(Kq̇ + Bq)] = a³(Bᵀq̇ + Cq), de φ = −0.2 a +0.2, con sus cuatro soluciones independientes;
  - la ganancia de energía es G(n) = σ_max(E_f^(1/2)·T·E_i^(−1/2)), con E la forma cuadrática del hamiltoniano reducido en los extremos (exigida definida positiva). Es la medida de C3e, reimplementada.

## 4. Puertas (si falla alguna, U59 es INCONCLUSO)

- **G0, la reproducción de C2:** con el juego y el fondo de C1, en el rebote, los autovalores de K_red coinciden con los publicados (§2) a ≤ 10⁻⁸ relativo en n = 2, 10 y 1 000.
- **G1, el ajuste en k²:** el polinomio identificado reproduce la cuadratura en un n de control (n = 37) a ≤ 10⁻¹⁰ relativo, en todas las matrices.
- **G2, la convergencia:** con rtol 10⁻⁹ y 10⁻¹⁰, G(n) cambia a ≤ 10⁻⁴ relativo (n = 2, 10³, 10⁵).
- **G3, la simplecticidad:** |TᵀΩ_fT − Ω_i| relativo ≤ 10⁻⁶ en todas las transferencias, con Ω la forma simpléctica a³(…).
- **G4, la energía:** E_i y E_f son definidas positivas en todos los n.

## 5. Veredicto (juego adoptado)

- **ESCALARES_SANOS:**
  - los dos autovalores de K_red son > 0 y S es invertible en todo el muestreo;
  - y G(3×10⁵)/G(10⁴) ≤ 2, la regla del plan.
- **INESTABILIDAD:**
  - algún autovalor de K_red < −10⁻⁹ relativo (un fantasma);
  - o G(3×10⁵)/G(10⁴) > 2 con G creciendo en los tres últimos n (gradiente).
- **INCONCLUSO:** falla una puerta, o no se cumple ninguna de las dos anteriores.

**Se informan siempre:**
- G(n) para todos los n, también con el juego de C1, frente a C3e;
- el mínimo de los autovalores de K_red;
- k_fís en el rebote para cada n. Los n con k_fís ≳ 1 M_Pl son prueba de crecimiento, no física.

## 6. Lo que no cubre

- n = 0 y n = 1 (C.3).
- Los vectores y los tensores (C.4; el explorador, C4, lo hizo de forma estructural).
- Las primeras decenas de e-folds de inflación (C.5).
- La acción cúbica y el acoplamiento fuerte.
- La continuación de la rama φ̇ → 0.

## 7. Lo que se espera

- **G0–G4 pasan.** G0 es la menos segura: una derivación independiente puede discrepar en una convención, y entonces hay que encontrar cuál.
- **ESCALARES_SANOS,** con confianza media. Con el juego adoptado, D(0) es unas 40 veces menor que con el de C1, así que el rebote es más brusco en φ, pero no tiene por qué ser inestable.

## 8. Coste y salidas

- **Coste:** la derivación simbólica, unos minutos; las transferencias de n alto, unos minutos cada una.
- **Código:** `chain/u59_escalares_rebote.py`.
- **Resultados:** `results/u59_escalares_rebote.json`.
- **Informe:** `docs/informes/INFORME_U59.md`.
- Interno.

---

## Enmienda 1 · 2026-11-15 22:35 (Madrid): cuadratura directa en lugar del ajuste en k², antes de la producción

**Qué se ha calculado ya, en la depuración:**
1. La derivación simbólica: 52 términos en L₂ y 39 integrandos radiales.
2. El ajuste por polinomios en k², que **falla:** la comprobación en n = 37 da un error relativo de 14. Los integrandos sueltos (por ejemplo, Y′²cos²r o Y′Y″sen 2r) no son polinómicos en k² por separado; solo lo son sus sumas.
3. **Una primera comparación con C2 en el rebote de C1,** con cuadratura directa. Los autovalores de K_red salen 0.7175966 y 2.60952408 (n = 2), 3.47931201 (n = 10) y 11.56368148 (n = 1 000), iguales a los publicados en las cifras impresas. El coeficiente de ζḟ es −3μ²a³ a 10⁻¹³, y ḟ no aparece en ningún otro término. Esta comparación es la de G0, hecha antes de su versión formal (≤ 10⁻⁸), que se repite en la producción.

**Los cambios:**
- **Las integrales sobre S³ se calculan por cuadratura de Gauss–Legendre compuesta en cada n,** con paneles de 16 puntos y unos 8(n + 2) puntos en total (mínimo 400). No se ajusta ningún polinomio.
- **G1 pasa a ser la convergencia de la cuadratura:** con el doble de puntos, las integrales cambian ≤ 10⁻¹⁰ relativo en n = 37 y en n = 3×10⁵.

Lo demás no cambia.

---

## Enmienda 2 · 2026-11-15 22:48 (Madrid): G1 ya se ha medido, y falla; se declara antes de la producción

**Lo que se midió:** al cronometrar la cuadratura antes de la producción se evaluó la comprobación de G1:
- **n = 37:** 7.7×10⁻¹⁵; pasa.
- **n = 3×10⁵:** con 8(n + 2) puntos frente a 16(n + 2), **2.0×10⁻¹⁰**, frente a ≤ 10⁻¹⁰. **Falla.**

Con 16 frente a 32 da 1.0×10⁻¹⁰. El error baja como el número de puntos, no exponencialmente: es redondeo acumulado en sumas de ~2×10⁶ términos con cancelaciones, no falta de resolución. **El umbral de la enmienda 1 estaba por debajo del suelo de redondeo para n = 3×10⁵.**

**Decisión:**
- **No se cambia el umbral.** Por la regla, U59 queda **INCONCLUSO** por G1, pase lo que pase con el resto.
- **La producción se hace igual,** con la cuadratura de la enmienda 1, y todo lo demás se informa como **informativo**: G0, la positividad, las transferencias y G2–G4.
- **Una U59b,** con un G1 fijado desde el suelo de redondeo y datos nuevos (otra madre), daría el veredicto.
