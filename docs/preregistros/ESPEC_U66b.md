# ESPEC_U66b · La salida caliente del cuscutón, con Γ_χ variado y la puerta de regresión corregida

*Escrito el 2026-11-21 a las 14:53 (Madrid), después de ver U66 (INCONCLUSO por G0) y antes de calcular U66b.*

## Qué cambia respecto a U66, y por qué

**La regla de G0 se cambia después de ver U66.** En U66, G0 exigía ≤ 10⁻⁶ relativo en N_inf y en χ_land.
- N_inf coincidió con U57b a ≤ 1.4×10⁻⁷.
- χ_land se desvió hasta 4.4×10⁻⁶. χ_land es el máximo de |χ| sobre los puntos de salida del integrador, una magnitud muestreada que depende de la secuencia de pasos.
- En U66b, G0 exige ≤ 10⁻⁶ en N_inf y ≤ 10⁻⁵ en χ_land.

**Datos nuevos.** La anchura Γ_χ de la cadena es fenomenológica, así que U66b calcula los 57 hijos de U57b con tres valores: Γ_χ/10, Γ_χ y 10Γ_χ. La fila con Γ_χ repite la producción de U66, que ya se ha visto.

El resto es como en U66: el sistema, las tres fases (25 giros exactos de χ y después el fluido hasta ρ_χ/ρ_r < 10⁻¹⁰), G1 (≤ 10⁻⁶) y G2 (rtol 10⁻¹², ≤ 0.01 dex y la misma clase, ahora en los tres valores de Γ).

## Veredicto

- **SALIDA_CALIENTE.** Pasan G0–G2, y para los tres valores de Γ los 57 hijos siguen en la rama en la fase 2 y tienen M > 0 en toda la fase 3.
- **SALIDA_FRIA.** Pasan las puertas, pero algún hijo rompe la rama con algún Γ. Se informa del Γ y del hijo.
- **INCONCLUSO.** Falla alguna puerta.

**Lo esperado:** SALIDA_CALIENTE. El margen del primer giro debería ir más o menos como log₁₀ Γ: unos 42.5 con Γ_χ/10 y unos 44.5 con 10Γ_χ, en el hijo más justo.

## Salidas

- `chain/u66_salida_caliente.py`, con la opción `--u66b`, escribe `results/u66b_salida_caliente.json`.
- El informe va en una sección U66b de `docs/informes/INFORME_U66.md`.
