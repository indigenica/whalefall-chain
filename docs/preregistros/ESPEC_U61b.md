# U61b · La curvatura y la huella con los hijos del cuscutón, con la proyección de la ligadura y datos nuevos · Preregistro

*2026-11-16 08:51 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («empezaría por U61b»). Escrito después de ver U61 (INCONCLUSO por G0 y G1) y antes de cualquier cálculo de U61b. Interno: no se exporta.*

## 1. Qué cambia respecto a U61, declarado (elegido después de ver U61)

1. **Datos nuevos.** f = 10^((k + ¼)/16) y 10^((k + ¾)/16), con k = −64…−1: 128 valores por A·H, los cuartos de la rejilla conjunta de U57 y U57b, con los 8 A·H y el juego adoptado. Cada fila se clasifica con `explora_cuscuton.job` (sin cambios), y los anfitriones son las filas con N_inf ≥ N_* y clase HIJO o INFLACION_CORTA.
2. **La proyección de la ligadura, como en `explora_cuscuton.run`.** El tramo del modo +1 se parte en los suyos:
   - **decel**, hasta que ä = 0 hacia abajo (o el giro φ → φ∞);
   - **onset**, hasta que ä = 0 hacia arriba, empezando con la proyección;
   - **end**, hasta que ä = 0 hacia abajo, empezando con la proyección;
   - **la cola**, hasta |χ| = χ_end, también empezando con la proyección.

   La proyección es su misma fórmula, reimplementada: d se reajusta para que se cumpla Friedmann.
3. **Las puertas, ajustadas a lo que hace el propio explora_cuscuton:**
   - **G0:**
     - los e-folds de la fase acelerada a ≤ 0.01 de su N_inf;
     - χ_land a ≤ 10⁻⁵ relativo (su χ_land es el máximo sobre los pasos aceptados, no el aterrizaje exacto);
   - **G1:** el residuo de Friedmann a ≤ 10⁻⁷, el orden que alcanza su propio integrador en los hijos (5.6×10⁻⁸).

## 2. Lo demás, como en U61

- La rama central del rebote.
- El aterrizaje, N_tot, el pivote, Ω_K,0, n_mem y ℓ_mem.
- Los criterios de visible y permitido (Planck), y la lectura con DESI.
- G2: rtol de 10⁻¹¹ frente a 10⁻¹², en los 5 anfitriones de menos e-folds.
- La medida: el promedio por A·H.
- Los veredictos HUELLA_OCULTA / RARA / POSIBLE / INCONCLUSO, con un mínimo de 20 anfitriones.

## 3. Lo que se espera (informado por U61)

- **G0–G2 pasan.**
- **n_mem de unos cientos a unos miles.**
- **f_vis de unos pocos por ciento, o 0:** la ventana visible es estrecha y la muestra pequeña. **HUELLA_POSIBLE o RARA,** con confianza baja.

## 4. Salidas

- **Código:** `chain/u61b_cuscuton_huella.py`, que importa `u61_cuscuton_huella` para las constantes.
- **Resultados:** `results/u61b_runs/` y `results/u61b_cuscuton_huella.json`.
- **Informe:** una sección de `INFORME_U61.md`.
