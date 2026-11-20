# U57b · El nacimiento con cuscutón, repetido con datos nuevos y dos reglas corregidas · Preregistro

*2026-11-15 19:57 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, siguiendo el «dale caña» del autor. Escrito después de ver U57 (INCONCLUSO por G1 y G2) y antes de cualquier cálculo de U57b. Interno: no se exporta.*

## 1. Qué cambia respecto a U57, y por qué (declarado: las dos reglas se eligen después de ver U57)

1. **Datos nuevos.** f = 10^((k + ½)/16), con k = −64…−1: los 64 puntos medios, en escala logarítmica, de la rejilla de U57. Los A·H (8) y los dos juegos del cuscutón son los mismos. Son 1 024 filas.
2. **La validez de la ligadura, solo donde decide.** El veredicto positivo (CUSCUTON_FERTIL, SOLO_…) solo usa hijos. Así que:
   - **G1** exige que **todas las filas HIJO** tengan un residuo de Friedmann ≤ 10⁻⁶;
   - el porcentaje de filas válidas en el resto se informa, y solo cuenta para un veredicto negativo: ESTERIL exige además ≥ 90 % de filas válidas en el juego.

   En U57, las 274 filas inválidas eran todas SIN_INFLACION o INFLACION_CORTA, cuya clase no decide fertilidad.
3. **El borde de la rama, aparte.** Un hijo con el mínimo de W relativo tras el rebote < 10⁻² está «en el borde de la rama». En U57 cambiaron de clase los dos hijos con 0.0018 y 0.0022.
   - **G2** exige que los demás hijos, a rtol = 10⁻¹², sigan siendo HIJO con N a ≤ 0.1 e-folds.
   - Los del borde se informan, y **no cuentan** para los intervalos de f del veredicto.

## 2. Lo demás, como en U57

- **El integrador** (`explora_cuscuton.run`, sin cambios), las clases, rtol = 10⁻¹¹ (y 10⁻¹² para G2).
- **G0:** la regresión frente a C1, a ≤ 10⁻⁸.
- **G3:** D > 0 en todas las filas.
- **La regla de los intervalos:** al menos 2 hijos, sin los del borde, y log₁₀(f_max/f_min) ≥ 0.1.
- **Los veredictos** CUSCUTON_FERTIL / SOLO_C1 / SOLO_EPS_PEQUENO / ESTERIL / INCONCLUSO.

## 3. Lo que se espera

- **CUSCUTON_FERTIL,** con confianza alta (U57, informativo).
- **El A·H mínimo con hijos,** 30 con los dos juegos.
- **Unos pocos hijos en el borde de la rama.**

## 4. Coste y salidas

- **Coste:** unos 8 min.
- **Código:** `chain/u57b_cuscuton_fondo.py`.
- **Resultados:** `results/u57b_runs/` y `results/u57b_cuscuton_fondo.json`.
- **Informe:** una sección de `docs/informes/INFORME_U57.md`.
