# U15 · La cola tibia frente a la verosimilitud Lyman-ABGD (arXiv:2206.08188) · Preregistro

*11 de octubre de 2026, 22:23 (Madrid). Tarea 2 de `docs/PLAN_NUBE_2026-10-11.md`. El plan la llamaba ESPEC_U12, pero ese número ya lo usa la tarea 6.*

**Lo que se sabe antes de calcular.**
- Del artículo (Hooper et al. 2022) se ha leído su resumen de las mezclas fría + tibia: m_WDM ≳ 7.2 keV (f_WDM − 0.1). Con esa fórmula, fracciones por debajo de ~0.1 no estarían restringidas a ninguna masa. Esa fórmula es un ajuste a sus cadenas; aquí se usa su verosimilitud, no la fórmula.
- No se ha ejecutado ningún código de U15 ni se ha evaluado la verosimilitud sobre ningún modelo.

## Motivo

La versión tibia del sector oscuro (q ≈ 1) deja una cola de un 1–4 % de materia oscura con temperatura efectiva de 0.2–0.35 keV. Con la cota extrapolada de Garcia-Gallego et al. 2025 (f < 0.14 m^1.1), eso es marginal (plan, §«Comprobación de reliquias calientes»). La verosimilitud Lyman-ABGD está pensada justo para supresiones generales de P(k), mezclas incluidas, y es pública: `montepython_public/montepython/likelihoods/Lya_abgd`, con los datos en `data/Lya_grid` (datos MIKE/HIRES, z = 4.2–5.4). Permite evaluar nuestras distribuciones tal cual, sin pasar por la familia térmica.

## Cálculo (`lya/u15_lya_abgd.py`)

- **Arnés:** el código de la verosimilitud se usa sin modificar, con un sustituto mínimo de las clases de MontePython (`Likelihood`, `data`, `io_mp`). Opciones del fichero `.data` por defecto (`smoothing_method = "None"`, `use_least_square = "None"`).
- **Equivalente ΛCDM:** la misma cosmología, con la especie ncdm pasada a CDM (ω_cdm + ω_ncdm). Se omite el reescalado η de la verosimilitud para grados de libertad relativistas extra, porque nuestro ΔN_eff ≤ 0.0075 da |η − 1| < 10⁻³. Es la única desviación del código original, y se dice.
- **Cosmología:** la de U4/U5 (`lya/u4a_class.py`, COSMO), con la salida a z = 0 y P_k_max = 300 h/Mpc. La verosimilitud exige σ8 y z_reio derivados; se toman de CLASS.
- **Perfil sobre los 9 parámetros astrofísicos** (T0a, T0s, γ_a, γ_s, Fz1–4, F_UV), con los rangos del `.param` de ejemplo de MontePython (`input/base2018TTTEEE_lya_abgd.param`): minimización acotada desde 8 puntos de partida, quedándose con el mejor.
- **Δχ² = χ²_min(modelo) − χ²_min(ΛCDM),** con el mismo arnés y la misma cosmología.

## V · Validación del arnés (antes de leer ninguna fila)

Con modelos térmicos (`thermal_params` de U4):
- WDM pura de 2 keV: **Δχ² > 4**;
- WDM pura de 10 keV: **Δχ² < 4**;
- mezcla de m = 1 keV y f = 0.5: **Δχ² > 4**;
- mezcla de m = 2 keV y f = 0.2: **Δχ² < 4**;
- Δχ² decreciente en m para la WDM pura en {2, 3.5, 5.3, 10} keV.

Si falla cualquiera de estos criterios, U15 se informa como **INCONCLUSO (arnés no validado)** y no se lee ninguna fila.

## Filas

- **B2, q = 1, R_Γ = 1.56** (la R_Γ limpia de U11): x₀ ∈ {10, 30}, K_sw ∈ {0, 1}. Son 4 filas, las mismas especificaciones que las filas B2 con q = 1 de U12.
- **Las 24 filas de U7c:** B2κ y B2ζ, R_Γ ∈ {1.39, 1.69, 2.09}, x₀ ∈ {10, 30}, K_sw ∈ {0, 1}, como en `lya/u7_class.py --u7c`.

Cada fila se clasifica como:
- **PERMITIDA:** Δχ² ≤ 4;
- **EXCLUIDA:** Δχ² > 4;
- **FUERA_DE_VALIDEZ:** si la verosimilitud devuelve un error de comprobación (`Error_cosmo` o `Error_equiv`), o si el ajuste de pesos no es fiable según sus derivadas (`self_distance_weighted` > 1, informativo).

Se informa también Δχ² frente a 6.18 (dos parámetros a 2σ), sin veredicto.

## Veredicto (sobre las 4 filas B2 con R_Γ = 1.56)

- **COLA_TIBIA_PERMITIDA:** las 4 filas válidas, PERMITIDAS.
- **COLA_TIBIA_EXCLUIDA:** las 4 filas válidas, EXCLUIDAS.
- **COLA_TIBIA_PARCIAL:** cualquier otra combinación con al menos 3 filas válidas.
- **INCONCLUSO:** menos de 3 filas válidas.

Las 24 filas de U7c se informan como tabla, comparadas con la puerta de forma de U7c (18/24) y con δA (22/24), sin veredicto propio.

## Lo que U15 no resuelve

- **La cosmología está fija.** La verosimilitud está pensada para usarse junto con Planck en una cadena; el perfil solo varía los parámetros astrofísicos.
- **El rango de la rejilla.** Las simulaciones cubren α ≤ 0.17 y δ ≤ 1. Una cola muy caliente da una supresión casi plana (δ), que puede caer en una zona poco poblada de la rejilla. Las comprobaciones propias de la verosimilitud lo detectan solo en parte.

## Enmienda 1 · 11 de octubre de 2026, 22:49 (Madrid)

**Qué se vio antes de esta enmienda.** La validación V falló: todos los Δχ² perfilados salieron cerca de 0. WDM pura de 2 keV: −0.57; de 3.5 keV: −2.09; de 5.3 keV: +0.37; de 10 keV: −0.10. Mezcla de 1 keV con f = 0.5: −2.68. La mezcla de 2 keV con f = 0.2 no llegó a calcularse. **No se ha evaluado ninguna fila del modelo.**

**Diagnóstico.** El perfil sobre los parámetros astrofísicos usaba solo los rangos del `.param` de ejemplo. Con ellos, el mínimo cae en historias térmicas fuera de la rejilla de simulaciones: en el ajuste ΛCDM, T0s = −5 (en el borde) y γ(z = 5.4) ≈ 0.58, por debajo del mínimo de la rejilla (0.776). Ahí el kriging extrapola, y la temperatura absorbe la supresión. Hooper et al. (2022, §2) dicen que rechazan los puntos cuyas T0(z) y γ(z) caen fuera de la rejilla térmica. Esa comprobación no está en `loglkl` del código público.

**Cambio.** El perfil impone, en los cuatro z de MIKE/HIRES (4.2, 4.6, 5.0, 5.4), las cotas de la rejilla que el propio código define: T0(z) ∈ [t0_min, t0_max] y γ(z) ∈ [slope_min, slope_max]. Los puntos fuera tienen χ² = ∞. Los rangos de F(z) del `.param` ya coinciden con `F_prior_min/max` del código. Nada más cambia: ni las filas, ni los umbrales, ni los criterios de V.

**Regla.** Si V vuelve a fallar con esta enmienda, U15 se informa como INCONCLUSO (arnés no validado) y no se evalúa ninguna fila.
