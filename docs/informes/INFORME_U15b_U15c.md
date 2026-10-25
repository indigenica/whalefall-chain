# U15b y U15c · La verosimilitud completa de Lyman-α (Planck + BAO + Lyman-ABGD): la rama viva, PERMITIDA, porque en esos bordes es indistinguible de CDM

*20 de octubre de 2026, 09:55 (Madrid).*
- *Preregistros:* `docs/preregistros/ESPEC_U15b.md` (notas 1–12) y `ESPEC_U15c.md` (enmiendas 1–4).
- *Código:* `lya/u15b/` (tablas, cadenas y análisis) y `lya/u15c/analyze_u15c.py`.
- *Resultados:*
  - `lya/results/u15b/U15B_RESULT.json` (análisis en WSL, 20 de octubre, 06:37);
  - `lya/results/u15c/U15C_RESULT.json` (06:48);
  - las tablas T²(k; s), en `lya/u15b/tables/`.
- *Ejecución:* MCMC con MontePython en WSL, 4 cadenas por caso, del 12 al 20 de octubre (U15b) y del 16 al 20 de octubre (U15c).
  - Las cadenas no están en el repositorio: el autor las conserva en WSL (`~/u15b/chains`).
  - El análisis lo ejecutó el autor con los analizadores de `main`, y aquí se suben sus JSON.

## 0. En una frase

**Veredictos preregistrados:**
- **U15b: DETENIDO** (las filas físicas, superadas por U17–U21; la validación PASA);
- **U15c: RAMA_VIVA_PERMITIDA_LYA_COMPLETO.**

**Los dos bordes de la banda que se probaron** (F = 0.5·F₀ con r = 1, y F = 0.35·F₀ con r = 3, con la población del colapso de U17d) **pasan la verosimilitud completa con un margen enorme:** p(s = 1)/p_max = 0.87 y 0.91, frente al umbral de 0.135.

**Pero la razón importa más que el número.** En esos casos el modelo suprime la potencia menos de un 0.15 % en k = 10 h/Mpc. Los datos solo detectan supresiones del 6–15 % ahí. Para la verosimilitud completa, esa materia oscura **es CDM**, y la posterior de s es plana.

## 1. La validación (U15b, mezclas térmicas)

| Mezcla | f₉₅ medido | Ventana preregistrada (Hooper et al. ×[0.6, 1.6]) | ¿Pasa? |
|---|---|---|---|
| 1 keV | 0.271 | [0.14, 0.38] (centro 0.24) | sí |
| 2 keV | 0.422 | [0.23, 0.61] (centro 0.38) | sí |

- Además, f₉₅(1 keV) < f₉₅(2 keV).
- La convergencia: R − 1 ≤ 0.021 (≤ 0.023 por mitades), con ESS_s ≈ 1350–1400.
- **La tubería reproduce, con un 10–12 % de holgura, el límite publicado de la mezcla térmica.**

**Las cuatro filas físicas de U15b** (la cola tibia del escenario K_sw ∈ {0, 1}) se detuvieron el 16 de octubre (Nota 7), sin leer su posterior, porque ese escenario dejó de ser el físico.
- Sus cadenas no convergieron (R − 1 = 0.05–0.17).
- El analizador no calcula nada de s con ellas, y **no se interpretan.**

## 2. U15c: los dos bordes

| Caso | F/F₀ | r | m_G (eV) | f_nc | δA (U21, cribado) | p(s = 1)/p_max (30 bins) | 20 / 45 bins | Jackknife (quitando una cadena) | Estado |
|---|---|---|---|---|---|---|---|---|---|
| `U15c_B2_F0.5_r1` | 0.5 | 1 | 5.5×10⁻¹⁷ | 0.980 | 0.094 | **0.874** | 0.888 / 0.889 | 0.853–0.893 | **PERMITIDA** |
| `U15c_B2_F0.35_r3` | 0.35 | 3 | 2.6×10⁻¹⁶ | 0.979 | 0.121 | **0.913** | 0.916 / 0.916 | 0.881–0.944 | **PERMITIDA** |

- **Convergencia:** R − 1 ≤ 0.017 (≤ 0.025 por mitades), ESS_s ≈ 2500–2900 y unos 140 000 puntos por caso.
- **La etiqueta no cambia con los bins ni al quitar una cadena:** no hay INDECISO.
- **La posterior de s es casi plana en todo el prior [0, 1.02]:**
  - la densidad varía entre 0.88 y 1.04 (en unidades de la media);
  - p(s = 0)/p_max = 0.94 y 0.92.

  Los datos no prefieren ni CDM pura (s = 0) ni el modelo (s = 1).

**Δχ²_min entre los puntos con s ≈ 1 y el mínimo global** (informativo): 5.9 en F = 0.5 y 1.9 en F = 0.35.
- Es el mínimo de unas 6 000 muestras frente al de 140 000, en 25 dimensiones, así que es un estadístico ruidoso.
- **El orden va al revés de lo físico:** el caso más frío (F = 0.5, δA = 0.094) da el Δχ² mayor. Eso apunta a ruido de muestreo, no a una señal.
- La posterior marginal en s = 1 solo baja un 13 %, lo que equivale a Δχ² ≈ 0.3.

## 3. Por qué la posterior es plana

T²(k) = P_modelo/P_ΛCDM a z = 0, de las tablas de las cadenas:

| k (h/Mpc) | 1 | 3 | 5 | 10 | 20 | 50 |
|---|---|---|---|---|---|---|
| U15c F = 0.5·F₀, r = 1, s = 1 | 1.000 | 1.000 | 1.000 | **0.999** | 0.996 | 0.968 |
| U15c F = 0.35·F₀, r = 3, s = 1 | 1.000 | 1.000 | 1.000 | **0.999** | 0.993 | 0.944 |
| 1 keV con f = 0.27 (su f₉₅) | 0.999 | 0.988 | 0.962 | **0.847** | 0.620 | 0.410 |
| 2 keV con f = 0.42 (su f₉₅) | 1.000 | 0.996 | 0.987 | **0.938** | 0.763 | 0.366 |

**La sensibilidad de los datos** (las mezclas en su límite del 95 %) es una supresión del 6–15 % en k = 10 h/Mpc. **La de la rama viva en esos bordes es del 0.07–0.14 %:** unas 50–100 veces por debajo.
- El colapso pone ≈ 98 % de la materia oscura en cuantos «no fríos». Pero nacen muy pronto, cuando T ≈ 0.4–0.8 MeV, y con p de unas pocas m_G, así que su escala de libre recorrido queda muy por debajo de lo que Lyman-α resuelve.

## 4. Qué significa

1. **La banda de F de U17e, U20b, U23 y U17f no se recorta por el lado de Lyman-α en F ≤ 0.5·F₀.** Lo confirma ahora la verosimilitud completa, marginalizando la cosmología y la astrofísica del IGM, no solo la puerta de cribado (δA).
2. **La puerta de cribado era conservadora en esos casos.** U21 les daba δA = 0.094 y 0.121, y ya pasaban. La verosimilitud completa no ve nada.
3. **Lo que U15c no prueba es el borde superior actual.** Los casos son de U17d (log-normal, c = 3, η = 1), definidos antes de que U19b midiera la población. La banda vigente llega a 0.7·F₀ y, con X ≲ 194, a F₀ (U17f), con δA de 0.28 a 0.38. Esos casos **no se han pasado por la verosimilitud completa.**
4. **La pregunta abierta de U17f (si vuelve F₀) podría decidirla Lyman-α completo**, además de la tercera masa de U19b. Ver §5.
5. **Lo que pidió la enmienda 2 de ESPEC_U15c:** U15c prueba el sector oscuro con la cosmología marginalizada, n_s incluido, no la ballena entera, que fija n_s ≈ 0.973. Como la posterior de s es plana, la conclusión no depende de n_s.

## 5. Una calibración informativa de la puerta de cribado (sin veredicto)

**¿A qué δA corresponden los límites del 95 % de la validación?** Si es bastante más que el umbral estricto de la puerta (0.234), las filas del borde actual (δA ≈ 0.28–0.38) podrían pasar la verosimilitud completa.

Calculado con CLASS en la nube (`U4.deltaA`, la misma métrica de U12–U17f), *después* de ver la validación:

| Mezcla térmica | δA |
|---|---|
| 1 keV, f = 0.14: el umbral estricto de la puerta de cribado (referencia «mixed 1 keV») | 0.234 |
| 2 keV, f = 0.30: referencia de la puerta | 0.290 |
| **1 keV, f = 0.27: su f₉₅ en la verosimilitud completa** | **0.429** |
| **2 keV, f = 0.42: su f₉₅ en la verosimilitud completa** | **0.389** |

- **El límite del 95 % de la verosimilitud completa cae en δA ≈ 0.39–0.43.** Es ≈ 1.7 veces el umbral estricto de la puerta, y está por encima también de la referencia ancha (0.336).
- **Con esa vara de medir, las filas del borde de U17f quedarían:**

  | | δA en F₀ | δA en 0.85·F₀ |
  |---|---|---|
  | Y = 1.05 / 1.20 (X ≈ 170 / 195) | 0.347 / 0.378: dentro | 0.262 / 0.295: dentro |
  | Y = 1.36 (X ≈ 220) | 0.409: en la frontera | 0.328: dentro |
  | Y = 1.70 / 1.98 (X ≈ 275 / 320) | 0.468 / 0.510: fuera | 0.383 / 0.421: dentro / en la frontera |

- **Ojo:** δA es un número solo. Las mezclas térmicas y la población del colapso no tienen la misma forma de T²(k), así que la comparación es orientativa.

**Esto no da ningún veredicto ni cambia ninguno.** Sugiere el siguiente paso:
- una **U15d preregistrada:** la verosimilitud completa en el borde actual, con la población medida de U19b (U17f: F = 0.7, 0.85 y 1·F₀ para Y ≈ 1.2–1.98);
- se ejecutaría en WSL, con la misma tubería y unos 4 días de cadenas.

Si F₀ pasara ahí con X ≈ 195–320, la vuelta del benchmark ya no dependería solo de la tercera masa.

## 6. Lo que U15b y U15c no resuelven

- **El borde superior actual de la banda,** con la población medida (§5).
- **T² se calcula en la cosmología fiducial.** La auditoría de U12/U15b midió un error de |ΔlnL| ≈ 0.4–0.7 a 4–5σ de Planck, muy por debajo de lo que haría falta para cambiar una posterior plana.
- **La cota de supermasivos, la isocurvatura (C) y la extrapolación del colapso,** que siguen como estaban.
