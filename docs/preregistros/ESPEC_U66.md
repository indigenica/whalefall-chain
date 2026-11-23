# ESPEC_U66 · La salida caliente del nacimiento con cuscutón

*Escrito el 2026-11-21 a las 14:46 (Madrid), antes de escribir el código y de calcular nada en este repositorio.*

**Declaración previa.** El 20 de noviembre, en una exploración con otro asistente de IA, fuera de este repositorio y con código propio, se calculó algo parecido («Blindaje 1A»): los 19 hijos del juego adoptado y los 38 del juego C1 salen calientes, con un margen de ≈ 10^43.5 en el primer giro de χ. El autor me pasó el código y las salidas, y los he leído. Esta prueba no usa ese código. Lo reimplementa sobre el integrador del repositorio (`chain/explora_cuscuton.py`, la clase `Cusc` y el mismo plan de tramos de U57b). Las expectativas de abajo se escriben, por tanto, conociendo ese resultado.

## 1. La pregunta

En `main-v4.tex` el supuesto abierto «The cuscuton as an effective theory» dice:

> en el primer punto de giro de χ tras la inflación, donde χ̇ = 0, la condición de rama necesita radiación, ρ_r > 3/(2a²), que el recalentamiento de la cadena tiene que haber producido para entonces. Los hijos calculados se detienen al final de la inflación, y que esto se cumpla se ha estimado, no calculado.

U66 lo calcula. La pregunta es si, con el recalentamiento efectivo de la cadena (Γ_χ de `data/cadena_cierre/thermal_pole.json`), los hijos de U57b se mantienen en la rama del cuscutón (W > 0) desde el giro de la madre hasta el final del recalentamiento.

## 2. Las filas

**Las madres.** Los 57 hijos de U57b (`results/u57b_runs/rows_rtol11.json`, `outcome = HIJO`): 19 del juego adoptado (μ = 0.16, m = 0.02, φ∞ = 0.25) y 38 del juego C1. Se toman sus mismos parámetros (A·H, f) y su mismo dato inicial en el giro de la madre.

**El acoplo.** Γ_χ = 1.4811708028419583×10⁻¹⁶ M_Pl, el valor de la cadena (`continuous.Gamma_MPl`). Está encendido durante toda la historia, como en la cadena, porque la fricción efectiva no tiene interruptor:

  χ̈ + (3H + Γ_χ) χ̇ + V′(χ) = 0,  ρ̇_r = −4Hρ_r + Γ_χ χ̇².

La radiación pasa a ser una variable dinámica, y deja de ser C_rad/a⁴. Lo demás es el sistema de `explora_cuscuton.py`: la ligadura U′(φ) = 3μ²H, φ̇ = 3μ²W/D, con W = ⅔ρ_r + ½χ̇² − 1/a² y D = 3μ⁴/2 − U″. La ligadura de Friedmann se comprueba y no se impone, salvo la proyección de d al empezar los tramos tardíos, como en U57b.

**Las fases.**
1. **Del giro de la madre al final de la inflación,** con el mismo plan de tramos de U57b. El final es el segundo cruce de ä = 0.
2. **La fase exacta de oscilación.** Desde el final de la inflación se integra exactamente hasta el 25.º punto de giro de χ (χ̇ = 0) y se vigila W/W_escala. La rama se rompe si W/W_escala < 10⁻¹². En cada giro se registra el margen M = log₁₀[(2/3)ρ_r a²].
3. **La fase de fluido.** Desde el final de la fase 2, χ se trata como un fluido con w = 0 y se integran ρ̇_χ = −3Hρ_χ − Γ_χρ_χ, ρ̇_r = −4Hρ_r + Γ_χρ_χ y 3(H² + 1/a²)(1 + ε_G) = ρ_χ + ρ_r, hasta que ρ_χ/ρ_r < 10⁻¹⁰.
   - Como en un giro de χ el término cinético se anula, la condición que se exige en esta fase es la pésima: (2/3)ρ_r a² > 1, es decir, M > 0.

## 3. Las puertas

- **G0, la regresión.** Con Γ_χ = 0, el integrador reproduce, en los 57 hijos, el N_inf y el χ_land de U57b con un error relativo ≤ 10⁻⁶, y la misma clase.
- **G1, la ligadura.** El residuo relativo de Friedmann, con la normalización de `explora_cuscuton.py`, es ≤ 10⁻⁶ en las fases 1 y 2 de todos los hijos.
- **G2, la convergencia.** Con rtol = 10⁻¹² frente a 10⁻¹¹, en todos los hijos el margen del primer giro cambia ≤ 0.01 dex y la clase no cambia.

## 4. El veredicto

- **SALIDA_CALIENTE.** Pasan G0–G2, y en los 57 hijos la rama sigue en la fase 2 (W/W_escala ≥ 10⁻¹²) y M > 0 en toda la fase 3.
- **SALIDA_FRIA.** Pasan G0–G2, pero al menos un hijo rompe la rama en la fase 2 o tiene M ≤ 0 en la fase 3. Se informa de cuántos y dónde.
- **INCONCLUSO.** Falla alguna puerta.

**Lo esperado,** por la exploración declarada y por una cuenta de servilleta: SALIDA_CALIENTE, con M ≈ 43 en el primer giro.
- **La servilleta:** tras ≈ 60 e-folds, 1/a² ha caído un factor ≈ e⁻¹²⁰ respecto al giro de la madre. En cambio, ρ_r ~ (Γ_χ/H) ρ_χ ~ 10⁻²¹ en el primer giro.
- **Si sale así,** el supuesto queda reducido a la existencia de un recalentamiento con Γ_χ > 0, que la cadena ya supone.

## 5. Lecturas informativas, sin veredicto

- **El efecto de Γ_χ en el número de e-folds:** ΔN_inf frente a U57b, que se espera ≲ 10⁻³.
- **T_rh,** con la convención de la cadena (g* = 106.75), frente a la de la cadena (2.29×10⁹ GeV).
- **El margen después del recalentamiento.** Con la historia estándar (radiación hasta la igualdad, luego materia y Λ), la condición pasa a ser ½(ρ + P) > 1/a².
  - En radiación el margen cae como a⁻² y en materia como a⁻¹.
  - Con Λ, ρ + P cae como a⁻³ y la rama del cuscutón se acabaría en el futuro.
  - Se estima, por hijo, el margen hoy y cuántos e-folds de futuro quedan con Ω_m = 0.315 y su Ω_K. Es una extrapolación analítica: no calcula la cosmología tardía con el cuscutón.

## 6. Qué no hace

- **No calcula** el potencial periódico del cuscutón (los rebotes repetidos de `NOTA_CUSCUTON_V4.md` §3.7), ni el «Blindaje 1B» de la exploración.
- **No calcula** el recalentamiento microscópico: usa la fricción efectiva de la cadena, como el resto de la cadena.
- **No modifica** `main-v4.tex`. Si el veredicto es SALIDA_CALIENTE, el supuesto abierto puede reescribirse en una versión posterior; lo decide el autor.

## 7. Salidas

- `chain/u66_salida_caliente.py` escribe `results/u66_salida_caliente.json`, con las puertas, el veredicto, una fila por hijo y el resumen por juego.
- El informe va en `docs/informes/INFORME_U66.md`.

## Enmienda del 2026-11-21 a las 18:34 (Madrid): el nombre del asistente

A petición del autor, la declaración previa ya no nombra al asistente de la exploración externa, porque no se sabe con certeza qué modelo se usó. Donde lo nombraba entre paréntesis, dice ahora «otro asistente de IA». Es el único cambio en el texto preregistrado: no toca las filas, las puertas, los umbrales ni la lógica del veredicto. El texto original sigue en el historial de git.
