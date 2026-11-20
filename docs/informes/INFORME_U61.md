# INFORME U61 · La curvatura y la huella del rebote con los hijos del cuscutón

**Preregistro:** `docs/preregistros/ESPEC_U61.md` (15 nov 2026, 23:57, Madrid).
**Código:** `chain/u61_cuscuton_huella.py`. **Resultados:** `results/u61_cuscuton_huella.json` y `results/u61_runs/produce.out`. Tardó 20 s.

Interno: no se exporta.

---

## 2026-11-16 00:00 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO (fallan G0 y G1)

| puerta | resultado |
|---|---|
| G0, el fondo frente a explora_cuscuton | **falla:** χ_land difiere hasta 5.5×10⁻⁶, frente a ≤ 10⁻⁶ (en 33 de 52 anfitriones). Los e-folds de la fase acelerada coinciden a ≤ 3.7×10⁻⁴ |
| G1, la ligadura (≤ 10⁻⁸) | **falla:** el residuo de Friedmann va de 5×10⁻⁷ a 5.8×10⁻⁵ en los 52 anfitriones |
| G2, la tolerancia | **pasa:** Ω_K,0 cambia ≤ 0.1 % y n_mem ≤ 2×10⁻⁸ |

**Por qué:** `explora_cuscuton.run` proyecta la ligadura (reajusta d para que se cumpla Friedmann) al empezar la fase acelerada y al terminarla, y mi integración por tramos no lo hace. Sin esa proyección, el error absoluto de la ligadura, congelado en el rebote, crece en términos relativos a medida que el universo se diluye. Es la «trampa ii» del plan, y la había olvidado al escribir este código. Las diferencias de 10⁻⁶–10⁻⁵ no cambian ningún orden de magnitud de lo que sigue, pero las puertas estaban escritas así.

### 2. Lo que se ve (informativo)

**La memoria del rebote llega mucho más lejos que con NPQT:** n_mem va de 381 a 2 199 (con NPQT era de 6 a 26).
- En el cuscutón, el rebote ocurre con un radio de curvatura grande (a_b ~ 10⁵ ℓ_P), y su escala física (la de |U(0)|, ≈ 5×10⁻³ M_Pl) cabe muchas veces en él.
- Por eso la huella llegaría a ℓ ≥ 2 con curvaturas minúsculas: basta |Ω_K,0| ≳ (2/(n_mem·D))², es decir **≈ 10⁻⁷**.

**Las fracciones** (52 anfitriones, promediadas por A·H):

| | Planck 2018 + BAO (Ω_K ≥ −0.0031) | DESI DR2 + CMB (Ω_K ≥ +10⁻⁴) |
|---|---|---|
| excluidos por la curvatura (f_excl) | 15 % | **100 %** (todos los hijos son cerrados) |
| permitidos y con la huella en ℓ ≥ 2 (f_vis) | **4.6 %** (3 de 52) | 0 |
| permitidos y con curvatura medible, \|Ω_K\| ≥ 10⁻³ (f_curv) | 0 | 0 |

- **Los tres visibles:**
  - A·H = 100, f = 0.154: Ω_K = −8.6×10⁻⁴, ℓ_mem ≈ 67;
  - A·H = 100, f = 0.178: Ω_K = −2.4×10⁻⁵, ℓ_mem ≈ 11;
  - A·H = 30, f = 0.340: Ω_K = −7.8×10⁻⁶, ℓ_mem ≈ 3.3.
- **Con las puertas en verde** habría sido **HUELLA_POSIBLE** (f_vis ≥ 10⁻²), frente a HUELLA_RARA con NPQT (U55b).
- **Con la cota de DESI DR2 a 2σ** no queda ningún hijo, porque todos son cerrados. Ese falsador es el mismo que con NPQT.
- **Siete anfitriones tienen |Ω_K,0| ≫ 10⁻³** (hasta −97: no son universos como el nuestro; están en f_excl). Para ellos la fórmula de ℓ_mem, con sen(D√|Ω_K|), pierde el sentido, pero no importa, porque están excluidos.

### 3. Lo que esto sugiere, junto con U59b (lectura, sin probar)

U59b mostró que, con el juego adoptado, **la mezcla de los escalares al cruzar el rebote no decae con el armónico** (energía de hasta ≈ 14 en todos los n hasta 3×10⁵, con k_fís de hasta unas M_Pl). Con NPQT, la memoria se borraba por encima de 2.5·n_mem (U56).

Si eso es así, el criterio de U61 (n_mem) **subestima** la huella:
- un hijo cerrado de cuscutón con pocos e-folds no solo tendría el rebote a la vista en ℓ bajos, sino también un espectro primordial con rasgos de orden uno allí donde los modos observables tenían k_fís ≲ M_Pl en el rebote;
- y hasta qué n llega eso depende de cuántos e-folds infló el hijo.

**Es una predicción potencialmente fuerte del nacimiento con cuscutón,** y también un riesgo, porque el espectro observado es casi perfectamente liso. Pide la versión con perturbaciones de esta prueba: la de U56, con las matrices de U59, que ahora sí están definidas por la acción.

### 4. Siguiente, como propuesta

- **U61b:** la misma prueba con la proyección de la ligadura de `explora_cuscuton` (o importando su `run` para el fondo y reintegrando solo el tramo entre el rebote y el aterrizaje para n_mem), y con datos nuevos (más f por A·H).
- **U62, el espectro primordial a través del rebote con cuscutón:** los modos de U59, con un estado inicial preregistrado en la contracción, propagados hasta la inflación con la acción de BKG. La pregunta es si la mezcla de orden uno de U59b deja rasgos en P_ζ para los hijos con N_tot ≈ 55–65, y para qué ℓ. Es la prueba que decidiría si el cuscutón con este potencial es compatible con un espectro liso.
- **U60, la anisotropía (Bianchi IX),** pendiente.

---

## U61b · 2026-11-16 08:55 (Madrid): con la proyección de la ligadura y datos nuevos

**Preregistro:** `docs/preregistros/ESPEC_U61b.md`. Los cambios de método y de umbral se eligieron después de ver U61, y así se declara.
**Código:** `chain/u61b_cuscuton_huella.py`. **Resultados:** `results/u61b_runs/rows.json` y `results/u61b_cuscuton_huella.json`. Tardó 3 min.

### Veredicto preregistrado: HUELLA_POSIBLE

| puerta | resultado |
|---|---|
| G0, frente a explora_cuscuton | **pasa:** χ_land a ≤ 6.2×10⁻⁶ (frente a ≤ 10⁻⁵); e-folds de la fase acelerada a ≤ 10⁻⁴ (frente a ≤ 0.01) |
| G1, la ligadura (≤ 10⁻⁷) | **pasa:** residuo máximo 4.4×10⁻⁹ |
| G2, la tolerancia | **pasa:** Ω_K,0 cambia ≤ 4×10⁻⁶ y n_mem ≤ 6×10⁻⁹ |

**Los datos:** 1 024 filas nuevas (128 valores de f por A·H), de las que 52 son anfitriones; 51 llegan al pivote en la integración propia.

| | Planck 2018 + BAO | DESI DR2 + CMB |
|---|---|---|
| excluidos por la curvatura | 12.7 % | 100 % |
| **permitidos y con el rebote en ℓ ≥ 2 (f_vis)** | **13.0 %** (5 de 51) | 0 |
| permitidos y con \|Ω_K\| ≥ 10⁻³ | 0 | 0 |

- **n_mem va de 400 a 2 205,** como en U61.
- **Los cinco visibles:**
  - A·H = 1 000, f = 0.54: Ω_K = −2.7×10⁻⁵, ℓ_mem ≈ 36;
  - A·H = 1 000, f = 0.84: Ω_K = −1.4×10⁻⁵, ℓ_mem ≈ 25;
  - A·H = 30, f = 0.097: Ω_K = −4.7×10⁻⁵, ℓ_mem ≈ 9;
  - A·H = 1 000, f = 0.90: Ω_K = −5.5×10⁻⁷, ℓ_mem ≈ 5;
  - A·H = 100, f = 0.18: Ω_K = −3.4×10⁻⁶, ℓ_mem ≈ 4.
- **Los visibles tienen todos 58–61 e-folds** de inflación tras aterrizar y una curvatura **indetectable** (|Ω_K| ≤ 5×10⁻⁵, por debajo de la sensibilidad de cualquier sondeo previsto): **el rebote podría verse en el CMB sin que la curvatura se vea.** Con NPQT era al revés: las dos cosas caían en el mismo sitio (U55).
- **El 13 % es un promedio por A·H** y lo domina A·H = 1 000 (3 de 7). Con la medida conjunta de todos los anfitriones sería 5/51 ≈ 10 %.
- **Con DESI DR2 a 2σ, todos excluidos,** porque todos son cerrados.

### Lo que esto significa (lectura)

- **Con el cuscutón, la huella del rebote es, para unos pocos de cada diez hijos que pueden alojarnos, una predicción y no algo escondido.** Su forma la decide U62: el espectro primordial con la mezcla de U59b.
- **El precio:** si U62 confirma rasgos de orden uno en el espectro, esa fracción de hijos quedaría excluida por la lisura del espectro observado.

---

## Nota · 2026-11-16 14:03 (Madrid): U62b corrige la lectura del §3

**El §3 de U61 sugería que el criterio n_mem subestima la huella,** porque la mezcla de U59b no decae con n. Esa mezcla resultó ser un efecto del marco del hamiltoniano congelado (`INFORME_U62.md`, nota de corrección y sección U62b).
- En el marco físico, la memoria del rebote se borra por encima de unas pocas n_mem: N_χ ≈ 0.03–0.1 en n ≈ n_mem, ≈ 0.01 en 2.4·n_mem y ≲ 10⁻⁴ en 5·n_mem.
- **El criterio de visibilidad de U61b (ℓ_mem ≥ 2) es, por tanto, el adecuado,** y lo que queda abierto es la forma de la huella en ℓ ≲ ℓ_mem (U62c).
