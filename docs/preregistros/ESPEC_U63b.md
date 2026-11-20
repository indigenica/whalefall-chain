# U63b · El acoplamiento fuerte del cuscutón en variables del gauge unitario · Preregistro

*2026-11-17 18:14 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U63b»). Es la repetición de U63 con datos nuevos y otras medidas, elegidas después de ver su fallo, y así se declara. Escrito antes del código. Interno: no se exporta.*

## 1. Qué cambia frente a U63, y por qué

**Lo que se vio en U63:**
- U63 falló en G1 porque, junto al giro de la madre (φ̇₀/H ≈ 2×10⁻⁷), los términos 1/φ̇₀ del gauge de radiación uniforme se cancelan en ~70 órdenes sobre la solución.
- La doble precisión no puede medir allí ρ₃ = ⟨|L₃|⟩/⟨|L₂|⟩.
- En la ventana y en la inflación (φ̇₀/H ~ 10⁻⁸), todo estaba bien condicionado.

**Las medidas nuevas** no suman términos grandes. Usan el desplazamiento temporal que lleva al gauge unitario del cuscutón (δφ = 0), ξ⁰ = f/φ̇₀, donde f es el δφ del gauge de U59.
- En el gauge unitario, el cuscutón entra con μ²φ̇₀·δN^k y no se refuerza cuando φ̇₀ → 0.
- Lo que se hace no lineal es la propia transformación, si ξ⁰ es grande.

**ρ₃ no se calcula.** El cúbico completo pediría la acción en el gauge unitario, que no está derivada.

## 2. El montaje

- **Hijos, modos y fondo:** los 8 hijos de U63, con n ∈ {2, 10, 30, 100, 300, 10³, 3×10³, 10⁴, 3×10⁴}, y el fondo de U62d (hacia delante desde el giro).
- **Los caminos:**
  - **el principal** empieza con el vacío del marco sin B_s en A2 (Δln a = 2 antes de la ventana) y llega hasta ρ_r/ρ = 10⁻⁵, ya en inflación;
  - **el informativo** empieza en A1 (0.1 e-fold antes del giro) y se calcula solo si el vacío está definido allí (E′ > 0, con los desplazamientos de U62c). Si no lo está, la fila se marca como NO_DEFINIDA.
- **Los campos:** q, q̇ y u = (α, β, f) como en U63.
- **Las medidas,** en cada instante y con amplitud de banda de vacío (×n^{3/2}/(π√2), como en la enmienda 1 de U63):
  - **M1** = |ζ − H·f/φ̇₀|: la curvatura sobre las hojas de φ constante;
  - **M2** = k_fís·|f|/φ̇₀: la inclinación de esas hojas (≥ 1, ∂φ deja de ser temporal);
  - **M3** = |f/φ̇₀|·|φ̈₀/φ̇₀|: el desplazamiento temporal frente al ritmo de cambio de φ̇₀ (≥ 1, la transformación al gauge unitario es no lineal).
- **El muestreo:** 120 instantes en la contracción, uniformes en ln a desde el inicio del camino; 120 en la ventana; 160 en la expansión, uniformes en ln a.
- **La escala fuerte:** n_sc es el menor n en el que máx_t máx(M1, M2, M3) ≥ 1, en el camino principal, con k_sc = n_sc/a_b.

## 3. Puertas (si falla alguna, U63b es INCONCLUSO)

- **G0, la normalización del vacío en A2:** ≤ 10⁻¹⁰.
- **G1, el wronskiano a lo largo del camino:** Ω(f_j, f_j*) se conserva a ≤ 10⁻⁶ en S5, en todas las filas del camino principal.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en el primer hijo y n = 100, 3×10³ y 3×10⁴. El máximo de máx(M1, M2, M3) cambia ≤ 10⁻³ relativo.
- **G3, el muestreo:** con el doble de instantes, en las mismas filas, el máximo cambia ≤ 5 %.

## 4. Veredicto

Sobre el camino principal.
- **INCONCLUSO:** falla una puerta.
- **ACOPLAMIENTO_FUERTE:** en algún hijo, n_sc ≤ n_mem, o máx(M1, M2, M3) ≥ 1 en algún n ≤ n(ℓ = 30).
- **PERTURBATIVO:** en todos los hijos, n_sc ≥ 10·n_mem (o no se alcanza hasta 3×10⁴), y en todos los n ≤ n(ℓ = 30), máx(M1, M2, M3) ≤ 0.1.
- **LIMITE:** cualquier otro caso.

**Se informan siempre:**
- cada medida frente a n y frente a la fase (contracción, ventana, inflación), con φ̇₀/H;
- k_sc frente a H_b, k_mem y (H_b²)^{1/3};
- el camino desde el giro: dónde y para qué n se pasa de 1, y en qué filas el vacío no está definido.

## 5. Lo que no cubre

- **El cúbico completo** en el gauge unitario.
- **El acoplamiento de modos distintos** y el biespectro.
- **Los bucles.**
- **La inflación profunda** más allá de ρ_r/ρ = 10⁻⁵.
- **Los vectores y los tensores.**

## 6. Lo que se espera

- **G0–G3 pasan.**
- **En el camino principal, PERTURBATIVO:** M1 ~ ζ ~ 10⁻⁵ en el infrarrojo, y M2 crece con n hasta n_sc ≫ n_mem.
- **En el giro de la madre,** M2 y M3 grandes, incluso con n bajos: el cuscutón no describe bien la madre cerca de su máximo.
- **Confianza media.**

## 7. Salidas

- **Código:** `chain/u63b_gauge_unitario.py`.
- **Resultados:** `results/u63b_gauge_unitario.json` y `results/u63b_runs/`.
- **Informe:** sección U63b de `docs/informes/INFORME_U63.md`.
- Interno.

---

## Enmienda 1 · 2026-11-17 23:15 (Madrid), tras ver qué filas no se podían calcular y antes de mirar ninguna medida

La producción se paró porque, en 8 filas del camino principal (n = 2 en siete hijos y n = 10 en el hijo A·H = 1 000, f = 0.542), la energía congelada E′ no es positiva en A2 ni acercándose al rebote, en pasos de 0.1 e-fold, hasta la ventana. El código solo atrapaba RuntimeError, y aquí saltó un ValueError.

1. **Esas filas se marcan NO_DEFINIDAS,** como las del camino informativo. Se excluyen de n_sc y del máximo infrarrojo, y se informan.
2. **El código tenía una puerta extra** que no está en este ESPEC («todas las filas principales definidas»). Se quita: valen G0–G3 tal como están escritas, calculadas sobre las filas definidas.
3. **Todavía no se ha mirado ningún valor de M1, M2 ni M3.** Solo se sabe qué filas fallan.
