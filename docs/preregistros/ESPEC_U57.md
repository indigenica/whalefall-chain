# U57 · El nacimiento con cuscutón, fondo homogéneo · Preregistro

*15 de noviembre de 2026, 19:46 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («si lo ves asumible dale caña»). Es la prueba A de `paper/zenodo_whale/PLAN_V4_CUSCUTON.md`. Escrito antes del código de U57. La hora que cuenta es la de su commit. Interno: no se exporta.*

## 1. La pregunta

Con la acción de REVENGE C1 (Einstein + cuscutón con el potencial de Boruah et al. + χ en la meseta de polo de la cadena + radiación, en un universo cerrado), ¿qué madres cerradas dan hijos (≥ 60 e-folds de inflación) sin salirse de la rama φ̇ > 0?

Se prueba con **dos juegos de parámetros del cuscutón**:
- **el de C1:** μ = 0.3, m = 0.05, φ∞ = 0.25 (ε_G = 0.049);
- **uno con ε_G pequeño:** μ = 0.16, m = 0.02, φ∞ = 0.25 (ε_G = 0.0039). Es el que probó la nota; el autor no ha elegido otro.

## 2. Lo que ya se sabe, y se declara (no es una prueba ciega)

- **`docs/informes/NOTA_CUSCUTON_V4.md`** y `results/explora_cuscuton.json` ya contienen esta misma rejilla con el juego de C1 (8 valores de A·H × 65 de f). Hay hijos con A·H ≥ 30 y ninguno con A·H ≤ 20.
- **Con el segundo juego** la nota solo cubrió A·H ∈ {48, 100, 300}, y en los tres hay hijos. **Lo nuevo de U57** es ese juego con A·H ∈ {3, 10, 20, 30, 1 000}, más las puertas y la regla de veredicto.
- Por lo anterior, **CUSCUTON_FERTIL es casi seguro** salvo que fallen las puertas. Lo que U57 añade es convertir la nota en un resultado preregistrado y contrastado: la regresión frente a C1, la ligadura, la tolerancia y D > 0.

## 3. El montaje

- **El integrador:** `chain/explora_cuscuton.py`, importado sin cambios (su `run`). Las ecuaciones son las de su cabecera.
- **La madre:** un universo cerrado de radiación en su giro (H = 0, φ = −φ∞), de tamaño A·H en radios de Hubble de la meseta, con χ en el mínimo y una fracción f = K_χ/ρ_r.
- **Las filas (1 040):**
  - A·H ∈ {3, 10, 20, 30, 48, 100, 300, 1 000};
  - f = 10^(k/16), con k = −64…0 (de 10⁻⁴ a 1);
  - los dos juegos de parámetros.
- **Las clases** de `explora_cuscuton`: HIJO (≥ 60 e-folds), INFLACION_CORTA, SIN_INFLACION, RAMA_ROTA, NUMERICO, ERROR.
- **La tolerancia:** rtol = 10⁻¹¹, y además todas las filas a rtol = 10⁻¹² para G2.

## 4. Puertas (si falla alguna, U57 es INCONCLUSO)

- **G0, la regresión frente a C1:** con los datos iniciales de C1 (`c1_row`), el rebote reproduce a, χ y χ̇ de `data/explorador/cuscuton_c1/result.json` (copiado sin cambios del explorador, commit `76d6713`) a ≤ 10⁻⁸ relativo.
- **G1, la ligadura:**
  - una fila es **válida** si su residuo de Friedmann máximo (`constraint_max`) es ≤ 10⁻⁶;
  - todas las filas HIJO tienen que ser válidas;
  - al menos el 90 % de las filas de cada juego tienen que serlo.

  Las filas inválidas no cuentan en el veredicto y se informan.
- **G2, la tolerancia:**
  - todas las filas HIJO a rtol = 10⁻¹¹ siguen siendo HIJO a 10⁻¹², con N a ≤ 0.1 e-folds;
  - al menos el 95 % de las filas válidas tienen la misma clase con las dos tolerancias.
- **G3:** D > 0 en todas las filas (ninguna con `D_negative`).

## 5. Veredicto

Para cada juego y cada A·H, sobre las filas válidas a rtol = 10⁻¹¹: hay **hijos en un intervalo de f** si hay al menos 2 hijos y log₁₀(f_max/f_min) ≥ 0.1 entre ellos.
- **INCONCLUSO:** falla una puerta.
- **CUSCUTON_FERTIL:** con los dos juegos hay al menos un A·H con hijos en un intervalo de f.
- **SOLO_C1:** solo con el juego de C1.
- **SOLO_EPS_PEQUENO:** solo con el segundo juego.
- **ESTERIL:** con ninguno.

**Se informan siempre**, por juego y A·H:
- las clases;
- el intervalo de f, el aterrizaje y los e-folds de los hijos;
- el A·H mínimo con hijos;
- la densidad del rebote;
- el residuo de Friedmann;
- −U/(3H²) al final de la inflación, frente a ε_G.

## 6. Lo que no cubre

- Las perturbaciones (U59).
- La cota observacional sobre ε_G (U58).
- La anisotropía (U60).
- La curvatura y la huella (U61).
- La continuación de la rama φ̇ → 0.

## 7. Lo que se espera

- **G0–G3 pasan.**
- **CUSCUTON_FERTIL**, con alta confianza, por §2.
- **Con el segundo juego, sin hijos en A·H ≤ 20:** la condición de rama de §3.3 de la nota no depende de ε_G. Confianza media.

## 8. Coste y salidas

- **Coste:** unos 6–8 min con 4 procesos.
- **Código:** `chain/u57_cuscuton_fondo.py`.
- **Resultados:** `results/u57_runs/` y `results/u57_cuscuton_fondo.json`.
- **Informe:** `docs/informes/INFORME_U57.md`.
- Interno.
