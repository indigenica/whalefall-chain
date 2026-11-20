# U64 · El orden cúbico completo en el gauge unitario del cuscutón · Preregistro

*2026-11-18 09:44 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con el cúbico completo en el gauge unitario»). Es la continuación de U63 y U63b. Escrito antes del código. Interno: no se exporta.*

## 1. La pregunta y la idea

**Lo que dejaron U63 y U63b:**
- U63 no pudo medir el cociente cúbico/cuadrático en el gauge de radiación uniforme, por las cancelaciones de los términos 1/φ̇₀.
- U63b estimó la escala fuerte con la no linealidad del paso al gauge unitario (M3): k_sc ≈ 0.04–0.12 M_Pl al entrar en el rebote.

**Por qué el gauge unitario** (δφ = 0):
- allí 2X = φ̇₀²/N² y la acción del cuscutón es −√h·(μ²φ̇₀ + N·U(φ₀)), sin ningún 1/φ̇₀;
- el orden cúbico completo (gravedad, χ, radiación y cuscutón) se evalúa sin cancelaciones grandes.

**La pregunta:** ¿dónde llega a 1 el cociente cúbico/cuadrático del lagrangiano ADM completo, y coincide con la escala de U63b?

## 2. El montaje

- **Los hijos, el fondo y los modos:** los 8 hijos de U63b y el fondo de U62d.
  - n ∈ {10, 30, 100, 300, 10³, 3×10³, 10⁴}. Sin n = 2, cuyo vacío no está definido en A2 en 7 hijos (U63b). Sin n = 3×10⁴, por el suelo de precisión de U63b.
  - El vacío en A2 y el camino hasta ρ_r/ρ = 10⁻⁵, con `u63b_gauge_unitario.integrate` y `u63_cubico_cuscuton.fields`, que dan la solución lineal en el gauge de U59 (w = 0): (ζ, x, α, β, f) y sus derivadas.
- **El paso al gauge unitario,** lineal, con ξ⁰ = f/φ̇₀ y ξ̇⁰ = ḟ/φ̇₀ − fφ̈₀/φ̇₀²:
  - ζ_u = ζ − Hξ⁰;
  - x_u = x − χ̇ξ⁰;
  - w_u = −qξ⁰ (q = θ̇₀);
  - α_u = α − ξ̇⁰;
  - β_u = β + ξ⁰;
  - f_u = ḟ_u = 0.
  - Las derivadas temporales vienen de las anteriores, con Ḣ, χ̈ = −3Hχ̇ − V′, q̇ = −Hq y φ̈₀ (por diferencias, con el paso de U63).
- **Las densidades:**
  - L₂ y L₃ del lagrangiano ADM completo de U59 (el Ltot de `u63_cubico_cuscuton.build_L23`), evaluados en esos campos;
  - ρ₃ = n^{3/2}·⟨|L₃|⟩/⟨|L₂|⟩ sobre una longitud de onda local en r ≈ π/2, con 64 puntos, como en la enmienda 1 de U63;
  - las dos cuadraturas y los dos modos de entrada se evalúan por separado, y se toma el máximo.
- **El muestreo** en el tiempo es el de U63b: 120, 120 y 160 instantes.
- **La escala fuerte:** n_sc es el menor n en el que máx_t ρ₃ ≥ 1, y k_sc = n_sc/a_b.

## 3. Puertas (si falla alguna, U64 es INCONCLUSO)

- **G0, el L₂ de este código frente al de U59** (como G0 de U63): ≤ 10⁻¹².
- **G1, las ligaduras lineales en el gauge unitario:** sobre la solución transformada, las ecuaciones de Euler–Lagrange de α y β de la acción cuadrática integrada (la forma de U59, Σ items) se anulan.
  - El residuo, relativo a la suma de |contribuciones|, debe ser ≤ 10⁻⁶.
  - Se comprueba en 3 hijos × n = 10, 300 y 10⁴ × 5 instantes.
  - Valida el paso de gauge, signos incluidos.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en el primer hijo y n = 100, 3×10³ y 10⁴. El máximo de ρ₃ cambia ≤ 10⁻³ relativo.
- **G3, la cuadratura y el muestreo:**
  - con 128 puntos en r, el máximo de ρ₃ cambia ≤ 10⁻³;
  - con el doble de instantes, cambia ≤ 5 %;
  - en las mismas filas que G2.

## 4. Veredicto

Las reglas son las de U63b, sobre n ≤ 10⁴.
- **INCONCLUSO:** falla una puerta.
- **ACOPLAMIENTO_FUERTE:** en algún hijo, n_sc ≤ n_mem, o máx ρ₃ ≥ 1 en algún n ≤ n(ℓ = 30).
- **PERTURBATIVO:** en todos los hijos, n_sc ≥ 10·n_mem (o no se alcanza hasta 10⁴), y ρ₃ ≤ 0.1 en todos los n ≤ n(ℓ = 30).
- **LIMITE:** cualquier otro caso.

**Se informan siempre:**
- ρ₃ frente a n y frente a la fase, con φ̇₀/H y |α_u|;
- k_sc frente al de U63b, H_b, k_mem y la escala de DGQ;
- qué parte del lagrangiano domina L₃ (gravedad, χ, radiación o cuscutón), cuando se pueda separar evaluando cada sector por su cuenta.

## 5. Lo que no cubre

- **El acoplamiento entre modos distintos** (las tríadas de S³; el explorador prepara C9 con la exportación de U63b) y el biespectro.
- **Los bucles.**
- **n = 3×10⁴.**
- **La inflación profunda.**
- **Los vectores y los tensores.**

## 6. Lo que se espera

- **G0–G3 pasan.**
- **ρ₃ ~ M3/2 de U63b,** porque el término que limita es α_u, como allí. Por tanto, n_sc parecido: 3×10³–10⁴.
- **Veredicto: el mismo que tendría U63b sin n = 3×10⁴,** es decir, ACOPLAMIENTO_FUERTE por la regla infrarroja en 2 hijos, o LIMITE. Confianza media.

## 7. Salidas

- **Código:** `chain/u64_cubico_unitario.py`.
- **Resultados:** `results/u64_cubico_unitario.json` y `results/u64_runs/`, con filas guardadas una a una.
- **Informe:** `docs/informes/INFORME_U64.md`.
- Interno.
