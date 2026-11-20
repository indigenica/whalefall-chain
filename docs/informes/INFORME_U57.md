# INFORME U57 · El nacimiento con cuscutón, fondo homogéneo

**Preregistro:** `docs/preregistros/ESPEC_U57.md` (15 nov 2026, 19:46, Madrid).
**Código:** `chain/u57_cuscuton_fondo.py`, que importa `chain/explora_cuscuton.py` sin cambios.
**Referencia de G0:** `data/explorador/cuscuton_c1/result.json` (REVENGE C1, commit `76d6713`, con su procedencia).
**Resultados:** `results/u57_runs/` (con `produce.out`) y `results/u57_cuscuton_fondo.json`.
**Todo en la nube:** 1 040 filas a rtol = 10⁻¹¹ en 213 s y las mismas a 10⁻¹² en 284 s.

Interno: no se exporta.

---

## 2026-11-15 19:55 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO (fallan G1 y G2)

| puerta | resultado |
|---|---|
| G0, la regresión frente a C1 | **pasa:** a, χ y χ̇ del rebote a 2×10⁻¹¹, 7×10⁻¹¹ y 5×10⁻¹¹ |
| G1, la ligadura (≥ 90 % de filas válidas por juego) | **falla:** 410 de 520 (79 %) con C1 y 356 de 520 (68 %) con ε_G pequeño. Las 69 filas HIJO son válidas (residuo ≤ 5.6×10⁻⁸) |
| G2, la tolerancia | **falla:** 2 de los 69 hijos pasan a RAMA_ROTA con rtol = 10⁻¹². Las clases coinciden en 762 de 766 filas válidas |
| G3, D > 0 | **pasa** |

**Por qué fallan:**
- **G1:** las 274 filas inválidas son todas SIN_INFLACION (167) o INFLACION_CORTA (107).
  - En esas madres el universo no infla, se diluye y el error absoluto de la ligadura, congelado en el rebote, crece en términos relativos. El residuo mediano es 9×10⁻⁵ y el máximo 7×10⁻².
  - La nota ya lo advertía («trampa ii» del plan). La regla de validez de ≤ 10⁻⁶ en todas las clases era demasiado exigente para filas cuya clase no depende de ese residuo. **Pero era la regla preregistrada.**
- **G2:** los dos hijos que cambian son los de margen más pequeño en la condición de rama: el mínimo de W relativo tras el rebote es 0.0018 y 0.0022. Uno es de C1 con A·H = 30, f = 0.075 y 85 e-folds; el otro, del juego con ε_G pequeño, con A·H = 30, f = 0.178 y 80 e-folds. Rozan W = 0, y la tolerancia decide de qué lado caen. **Es física en el borde, no ruido:** es el borde de la rama φ̇ > 0.

### 2. Lo que se ve, informativo

Solo filas válidas, a rtol = 10⁻¹¹. Un hijo es una madre que da al menos 60 e-folds:

| A·H | hijos (C1) | f de los hijos | aterrizaje | e-folds | hijos (ε_G pequeño) | f | aterrizaje | e-folds |
|---|---|---|---|---|---|---|---|---|
| 3, 10, 20 | 0 | — | — | — | 0 | — | — | — |
| 30 | 5 | 0.024–0.15 | 4.16–4.66 | 61–85 | 4 | 0.10–0.49 | 4.20–4.49 | 65–80 |
| 48 | 7 | 0.056–0.21 | 4.44–5.42 | 74–135 | 3 | 0.15–0.27 | 4.55–5.38 | 83–138 |
| 100 | 12 | 0.042–0.56 | 4.19–7.48 | 62–358 | 10 | 0.18–0.87 | 4.09–6.02 | 60–194 |
| 300 | 13 | 0.12–1 | 4.28–8.03 | 66–444 | 6 | 0.32–1 | 4.22–5.44 | 66–142 |
| 1 000 | 7 | 0.27–1 | 4.21–6.89 | 63–279 | 2 | 0.75–1 | 4.51–4.78 | 81–96 |

- **Con las puertas en verde, el veredicto habría sido CUSCUTON_FERTIL.** Con los dos juegos hay A·H con hijos en un intervalo de f de entre 0.12 y 1.12 décadas.
- **El A·H mínimo con hijos es 30 con los dos juegos,** como esperaba el ESPEC: la condición de rama no depende de ε_G. Por debajo de A·H ≈ 24–30 no hay hijos con ninguno de los dos.
- **La densidad del rebote la fija |U(0)| y es la misma en todas las filas:** 1.50×10⁻⁴ M_Pl⁴ con C1 y 2.48×10⁻⁵ con ε_G pequeño. Siempre por debajo de la de Planck.
- **Reproduce la nota** con el juego de C1, fila a fila. Con el segundo juego añade A·H = 30 y 1 000, que sí tienen hijos, y 3, 10 y 20, que no.

### 3. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0–G3 pasan»: fallado** en G1 (una regla demasiado exigente) y en G2 (el borde de la rama).
- **«CUSCUTON_FERTIL»:** sería el resultado con las puertas en verde.
- **«Sin hijos en A·H ≤ 20 con el segundo juego»: acertado.**

### 4. Qué haría falta para un veredicto (propuesta)

Una U57b con datos nuevos (otra rejilla de f, intercalada) y dos cambios fijados de antemano:
1. **La validez de la ligadura solo donde decide la clase.** Por ejemplo: residuo ≤ 10⁻⁶ hasta el final de la fase acelerada en las filas que inflan, y en las que no inflan, solo hasta el giro de χ. O con la proyección del script aplicada en todas las etapas.
2. **G2 sin las filas en el borde de la rama.** Se excluirían las que tengan el mínimo de W relativo < 10⁻², contándolas aparte como «borde de la rama».

### 5. Para la prueba B (U58): lo que dice la literatura sobre ε_G

Leído hoy en el artículo (Afshordi, Chung, Doran y Geshnizjani 2007, «Cuscuton cosmology: dark energy meets modified gravity», arXiv:astro-ph/0702002, ecs. 14 y 34):
- **El cuscutón cuadrático, V = V₀ + ½m²φ², es «early dark energy» con fracción constante:** Ω_Q = −ΔM_p²/M_p² = 3μ⁴/(2M_p²m²). Está acotado al 95 % de confianza:
  - **Ω_Q < 3.6 %** con WMAP3 sola;
  - **< 2.7 %** añadiendo SNe Ia y estructura a gran escala;
  - **< 1.6 %** añadiendo además el bosque Lyα de SDSS.
- **El potencial de Boruah,** cerca de φ∞, es una parábola **invertida:** U″(φ∞) = −4φ∞². Con m² → −4φ∞², la misma fórmula da |Ω_Q| = 3μ⁴/(8φ∞²) = ε_G, pero **con el signo opuesto:** Ω_Q = −ε_G < 0. La masa de Planck efectiva crece en lugar de bajar.
- **Consecuencia:** la cota publicada es para Ω_Q > 0 y no se aplica tal cual. En magnitud apunta a |ε_G| ≲ 2 %, lo que coincide con la estimación de la nota: el juego de C1 (4.9 %) estaría fuera y el pequeño (0.4 %) dentro. Pero hay que **rehacer el contraste con el signo negativo**: crecimiento de las perturbaciones aumentado en lugar de suprimido, ISW y BBN. Es la tarea de U58.
- **El mapeo con λ de Hořava–Lifshitz** (Afshordi 2009, arXiv:0907.5201, |λ − 1| < 0.014 al 95 %): el resumen no da la fórmula. Queda por leer en el cuerpo del artículo.

---

## U57b · 2026-11-15 20:05 (Madrid): la repetición con datos nuevos y las dos reglas corregidas

**Preregistro:** `docs/preregistros/ESPEC_U57b.md` (19:57). Las dos reglas se eligieron después de ver U57, y así se declara.
**Código:** `chain/u57b_cuscuton_fondo.py`. **Resultados:** `results/u57b_runs/` y `results/u57b_cuscuton_fondo.json`.
**Coste:** 1 024 filas en 217 s, y en 272 s a rtol = 10⁻¹².

### 1. Veredicto preregistrado: CUSCUTON_FERTIL

| puerta | resultado |
|---|---|
| G0, la regresión frente a C1 | **pasa:** 7×10⁻¹¹ |
| G1, los 57 hijos con residuo ≤ 10⁻⁶ | **pasa** |
| G2, los hijos fuera del borde de la rama, estables a rtol = 10⁻¹² | **pasa:** ninguno cambia; tampoco hay ninguno en el borde (W_rel < 10⁻²) |
| G3, D > 0 | **pasa** |

**Las filas válidas en el resto** (para un veredicto negativo, que no hace falta) son el 79 % con C1 y el 68 % con ε_G pequeño, como en U57.

**Los hijos de U57b**, sin los del borde:

| A·H | C1: hijos | aterrizaje | e-folds | ε_G pequeño: hijos | aterrizaje | e-folds |
|---|---|---|---|---|---|---|
| 3, 10, 20 | 0 | — | — | 0 | — | — |
| 30 | 2 (en 0.06 décadas de f) | 4.21–4.23 | 63–64 | 3 | 4.11–4.40 | 61–75 |
| 48 | 6 | 4.19–5.37 | 62–131 | 4 | 4.33–5.08 | 71–116 |
| 100 | 15 | 4.21–7.36 | 63–341 | 7 | 4.84–6.51 | 100–245 |
| 300 | 10 | 4.21–8.03 | 63–444 | 3 | 4.71–5.52 | 92–149 |
| 1 000 | 5 | 4.15–7.00 | 60–293 | 2 (en 0.06 décadas) | 4.22–4.37 | 66–74 |

### 2. Lo que queda establecido

- **Con el cuscutón, una madre cerrada de radiación con un condensado de χ da hijos de ≥ 60 e-folds,** con los dos juegos de parámetros, incluido uno con ε_G = 0.4 %.
- **Solo si es grande:** A·H ≥ 30. Por debajo, la rama φ̇ > 0 se rompe al entrar en la inflación. Es el régimen opuesto al del padre de radiación de NPQT (A·H ≤ 1).
- **Los aterrizajes van de 4.1 a 8 M_Pl.** No hay rama eterna. Varios hijos tienen 60–75 e-folds, justo donde la curvatura y la huella del rebote pueden no estar escondidas (prueba E del plan, U61).

### 3. Siguiente (el plan)

- **U58,** la cota sobre ε_G. Ver §5 de arriba: la literatura acota el signo opuesto, así que hay que rehacer el contraste.
- **U59,** la auditoría independiente de los escalares a través del rebote.

---

## Nota · 2026-11-15 20:06 (Madrid): la cota sobre ε_G, leída en Afshordi (2009); corrige §5

Leído en el cuerpo de Afshordi (2009, arXiv:0907.5201, ecs. 6, 13 y 14):
- **El equivalente de Hořava–Lifshitz:** el cuscutón cuadrático equivale a Hořava–Lifshitz de baja energía con μ² = −m² = (1 − λ)/(8πG_N) (ec. 6). El potencial es una parábola **invertida**, como la de Boruah cerca de φ∞.
- **El desajuste de masas de Planck:** M²_p,IR = M²_p,UV − 3μ⁴/(2m²) = [1 + (3/2)(1 − λ)]·M²_p,UV (ec. 13). Comparando con 3H²(1 + ε_G) = ρ: **ε_G = (3/2)(1 − λ)**, es decir, **λ − 1 = −(2/3)·ε_G.**
- **La cota es de dos lados:** λ − 1 = 0.003 ± 0.014 (95 %; ec. 14, de su referencia [6]). En ε_G da **−0.025 ≤ ε_G ≤ 0.017** al 95 %.

**Corrige §5 de esta nota:** el signo no impide aplicar la cota, porque la de λ es de dos lados.

**Consecuencias para los dos juegos de U57/U57b:**

| juego | ε_G | λ − 1 | frente a la cota |
|---|---|---|---|
| C1 | 0.049 | −0.032 | **excluido**, a unas 5σ con σ ≈ 0.007 |
| ε_G pequeño | 0.0039 | −0.0026 | **dentro** |

**Una advertencia del mismo artículo:** el resumen y el final dicen que con λ < 1, que es nuestro signo, las perturbaciones de longitud de onda corta pueden tener energía arbitrariamente negativa, una posible inestabilidad no perturbativa en el ultravioleta. Es una pregunta para U59 y para el supuesto «completación UV» de la v4. No la he analizado.

**Para U58** queda la parte cuantitativa: Δr/r y Δn_s con las perturbaciones con cuscutón durante el slow-roll, y la comprobación de que la cota de su referencia [6] se obtuvo sin imponer un signo. La lectura de hoy indica que sí (0.003 ± 0.014).
