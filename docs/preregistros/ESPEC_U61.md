# U61 · La curvatura y la huella del rebote con los hijos del cuscutón · Preregistro

*2026-11-15 23:58 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, siguiendo el «lo que puedas avanzar» del autor. Es la prueba E de `paper/zenodo_whale/PLAN_V4_CUSCUTON.md`, la versión con el fondo (como U55b). La numeración salta U60, la prueba D (anisotropía), que queda para después. Escrito antes del código de U61. Interno: no se exporta.*

## 1. La pregunta

Los hijos del cuscutón aterrizan bajo (4–8 M_Pl), y varios tienen menos de 70 e-folds. Con el juego adoptado:
- ¿qué fracción de los nacimientos que pueden alojarnos tiene una curvatura permitida por Planck y su rebote a la vista en ℓ ≥ 2?
- ¿qué Ω_K,0 predicen?

Es U55b con los hijos del cuscutón en lugar de los de NPQT.

## 2. Lo que ya se sabe, y se declara

- **U55b (NPQT):** HUELLA_RARA, f_vis = 3.8×10⁻⁴; y ningún nacimiento visible con la curvatura de DESI DR2.
- **U57 y U57b (este juego):** 52 filas con al menos 51.6 e-folds según `explora_cuscuton`, todas con A·H ≥ 30. La nota estimó |Ω_K,0| ≤ 2×10⁻⁵ en su muestra, con una aproximación (N_HORIZON = 57).
- **No se ha calculado** n_mem para el cuscutón, ni Ω_K,0 con el pivote de la cadena.

## 3. El montaje

**Los anfitriones:** las filas de U57 y U57b con el juego adoptado y N_inf ≥ N_* = 51.59 según `explora_cuscuton` (52 filas).

**La medida:** la rejilla conjunta de f de U57 y U57b (129 valores log-uniformes por A·H), con el mismo peso por A·H. f_vis es el promedio, sobre los A·H con anfitriones, de la fracción de anfitriones visibles y permitidos.

**El fondo** se integra aquí con `explora_cuscuton.deriv` (sin cambios), por tramos como su `run`:
1. el modo −1 hasta φ = −φ∞/2;
2. el modo 0 hasta φ = +φ∞/2, que cruza el rebote en φ = 0;
3. el modo +1 hasta el final de la inflación, cuando |χ| = χ_end = 0.336 tras el aterrizaje (el final de la cadena, como en U55b).

Se usa DOP853 con rtol = 10⁻¹², y se guarda la solución densa.

**Lo que se calcula por anfitrión:**
- **el aterrizaje:** el primer χ̇ = 0 tras el rebote con |χ| > 1;
- **N_tot,** de ese aterrizaje al final de la inflación;
- **a_*H_*** cuando faltan N_* e-folds para el final, y con él **Ω_K,0 = −(222.53/(a_*H_*))²**;
- **n_mem**, el máximo de aH y de a√(|R|/6) entre el rebote y el aterrizaje, con R = 6(Ḣ + 2H² + 1/a²);
- **ℓ_mem** = √(n_mem² − 1)·sen(D√|Ω_K,0|), con D = 3.1158 (U55).

**Visible** si ℓ_mem ≥ 2. **Permitido** si Ω_K,0 ≥ −0.0031 (Planck 2018 + BAO, a 2σ). Se informa también la lectura con DESI DR2 + CMB (Ω_K ≥ +10⁻⁴).

## 4. Puertas (si falla alguna, U61 es INCONCLUSO)

- **G0, el fondo frente a explora_cuscuton:**
  - χ_land a ≤ 10⁻⁶ relativo;
  - los e-folds entre el comienzo y el final de la fase acelerada (ä = 0 en los dos) a ≤ 0.05 de N_inf;
  - en todos los anfitriones.
- **G1, la ligadura:** residuo de Friedmann ≤ 10⁻⁸ relativo, hasta el final de la inflación, en todos los anfitriones.
- **G2, la tolerancia:** con rtol = 10⁻¹¹, Ω_K,0 cambia ≤ 1 % y n_mem ≤ 10⁻⁴ en los 5 anfitriones de menos e-folds.

## 5. Veredicto

- **INCONCLUSO:** falla una puerta, o hay menos de 20 anfitriones.
- **HUELLA_OCULTA:** f_vis = 0.
- **HUELLA_RARA:** 0 < f_vis < 10⁻².
- **HUELLA_POSIBLE:** f_vis ≥ 10⁻².

**Se informan siempre:**
- Ω_K,0, n_mem, ℓ_mem y N_tot por anfitrión;
- la fracción excluida por la curvatura (f_excl);
- la fracción con una curvatura medible, |Ω_K,0| ≥ 10⁻³ y permitida (f_curv);
- la lectura con DESI DR2.

## 6. Lo que se espera

- **G0–G2 pasan.**
- **n_mem** de orden 1–30, como con NPQT.
- **f_excl y f_curv mayores que con NPQT**, por los aterrizajes bajos, pero con solo 52 anfitriones y la medida discreta en f. **HUELLA_RARA o POSIBLE,** con confianza baja.

## 7. Salidas

- **Código:** `chain/u61_cuscuton_huella.py`.
- **Resultados:** `results/u61_cuscuton_huella.json`.
- **Informe:** `docs/informes/INFORME_U61.md`.
