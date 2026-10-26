# ESPEC_U27 · La coherencia de los hijos directos (no eternos) con la meseta de polo

*25 de octubre de 2026, 18:51 (Madrid). Claude, en la nube, sobre `main` tras la publicación de la v2 del preprint. Preregistro escrito antes de cualquier código o cálculo de U27.*

## 1. Pregunta

La v2 sitúa nuestro bolsillo en la rama eterna: el campo aterriza por encima de χ_et ≈ 541 M_Pl. Eso tiene dos costes:
- exige una meseta plana hasta ~10³ M_Pl (§8, «Field range»);
- trae consigo el problema de la medida (§8, «Measure»).

En la misma rejilla de nacimientos, la mayoría de los hijos fértiles son **directos**: aterrizan por debajo de χ_et.
- Con ≳ 100 e-folds, su cosmología observable es la misma, porque los últimos ~52 e-folds recorren el mismo atractor.
- Con N ≈ χ³/(3C), un aterrizaje en 10–50 M_Pl da 950–1.2×10⁵ e-folds.

Lo único que la v2 da para su nacimiento es la tolerancia a los varios rebotes del padre (D11). Esa tolerancia, el 3.6 % de C45j, **solo se calculó para los hijos eternos**.

**U27 pregunta:** ¿qué fracción de los hijos fértiles son directos y toleran perturbaciones de amplitud inflacionaria a lo largo de los rebotes del padre, en el mismo modelo que dio el 3.6 %? ¿Cuántos de ellos aterrizan en 10–50 M_Pl?

## 2. Qué no cambia

- **El modelo de coherencia es el de C45j:** polvo LTB-QT y ε\*_n = ε_est(A·H, αH², L·H)/n_rebotes, con la amplitud equivalente ε_phys = (2/3)·ζ0·π²·(A·H/L·H)² de C45d.
  - Un hijo es **viable** si ζ0_max ≥ ζ0 = 4.6×10⁻⁵.
  - U27 no mejora ese modelo: la coherencia con χ inhomogéneo (C45k) sigue pendiente, y el aviso de C45j (rebotes dominados por la energía cinética de χ) vale igual para los hijos directos.
- **La rejilla** es la de C45h:
  - polvo (`frac_rad = 0`);
  - la ventana sep ∈ {1e-4, 1e-5, 1e-6};
  - los (λ, σ) archivados en `data/c45b/c45h_results.json`.
- **La fertilidad** es la de `frw_c.py`: lanzamiento a 60 e-folds, `chi_abs_max > 5` y ligadura < 10⁻³.
- **Sub-Planck** quiere decir E_phys_max < 1, como en S1 y U22b.
- **Las clases**, con χ_land = `chi_end` (el campo tras ≈ 60 e-folds) y χ_et de la meseta, como en S1 y U22b:
  - **E (eterno):** fértil, sub-Planck y χ_land > χ_et;
  - **D (directo):** fértil, sub-Planck y χ_land ≤ χ_et;
  - **D₁₀₋₅₀:** los D con aterrizaje `chi_abs_max` (el máximo de |χ| tras el rebote; ver `RESPUESTA_2.md`) entre 10 y 50 M_Pl, ambos incluidos.

## 3. Paso 0: el código de C45j y su validación (puerta)

El código de C45j (`ltb/multibounce_viability.py`, con sus dependencias) y sus entradas (`results/c45h_content_*`) están en el proyecto local del autor (`claude_dual/c45b`). En este repositorio solo está su resultado (`data/c45b/c45j_results.json`).

1. **El autor copia** a `components/c45j_ltb/` el script, lo que importa y los ficheros `c45h_content_dust*` que lea, sin modificarlos. Se añaden con su sha256 a `components/MANIFEST.json`.
2. **Validación V1 (puerta).** Con la meseta de suma y el código tal cual, se reproduce `data/c45b/c45j_results.json`:
   - los mismos 91 hijos;
   - ε_est y ζ0_max por hijo, con una diferencia relativa ≤ 10⁻⁶;
   - `frac_viable_LH1_over_fertile` = 0.035907 (exacto a 10⁻⁹).
3. **Si V1 falla**, U27 se para y se informa como INVÁLIDO. Reimplementar ε_est desde los informes de C45d exigiría una enmienda fechada antes de calcular.
4. **Para la meseta de polo,** la clasificación (fértil, χ_land, n_rebotes, E_phys_max) se recalcula con `chain/u22b_births.grid` en modo 4 (C = 0.3516).
   - **Validación V2 (puerta):** el resumen tiene que coincidir con el `pole` guardado en `results/u22b_births.json` (n_fertile = 886, n_window_fertile = 556, f_ES = 0.152878).
   - Si (A·H, αH²) dependen de la meseta en el código de C45j, se recalculan con el polo. Si no, se usan los del padre. El informe dirá cuál de los dos casos se dio.

## 4. Filas

Todas en la ventana, con polvo y ζ0 = 4.6×10⁻⁵:

| Fila | Meseta | Clase | L·H | Papel |
|---|---|---|---|---|
| R1 | polo | E | 1 | referencia: el 3.6 % de C45j, rehecho con el polo |
| R2 | polo | D | 1 | **veredicto principal** |
| R3 | polo | D₁₀₋₅₀ | 1 | **veredicto de la banda** |
| R4 | polo | E, D, D₁₀₋₅₀ | 0.3 y 3 | sensibilidad a la escala de la perturbación (informativa) |
| R5 | suma | E, D, D₁₀₋₅₀ | 0.3, 1 y 3 | continuidad con C45j y la v1 (informativa) |

Por fila se informa:
- el número de hijos y de viables;
- **f_v(X) = viables de la clase X / fértiles de la ventana**, con la misma normalización que el 3.6 %;
- la fracción condicionada (viables de X / hijos de X);
- el desglose por separación;
- la distribución de n_rebotes y de A·H;
- la mediana de ζ0_max.

## 5. Veredicto

**Veredicto principal (R2, polo, L·H = 1):**
- **RAMA_DIRECTA_VIABLE** si f_v(D) ≥ 5 %, el umbral que ESPEC_C45K fijó para que D11 pase a VIABLE.
- **RAMA_DIRECTA_MARGINAL_COMO_ETERNA** si f_v(E) ≤ f_v(D) < 5 %, con f_v(E) de R1.
- **RAMA_DIRECTA_MARGINAL_PEOR** si 1 % ≤ f_v(D) < f_v(E).
- **RAMA_DIRECTA_NO_COHERENTE** si f_v(D) < 1 %.

**Veredicto de la banda (R3, polo, L·H = 1):**
- **BANDA_10_50_POBLADA** si hay al menos 5 hijos viables en D₁₀₋₅₀ y f_v(D₁₀₋₅₀) ≥ 1 %.
- **BANDA_10_50_VACÍA** en cualquier otro caso.

**Qué cambiaría en el texto** (en una v3, decisión del autor):
- **Con VIABLE o COMO_ETERNA y la banda POBLADA:** la rama directa pasa a ser una alternativa con el mismo peso en D9 y en §8. El requisito de planitud baja de ~10³ M_Pl a 10–50 M_Pl.
  - Para nuestro bolsillo, las predicciones dejan de depender de los pesos dentro de la inflación eterna. La fracción de nacimientos de cada rama sigue dependiendo del prior sobre los padres, como ya dice §8.
  - Sigue abierto, y se dirá: 10–50 M_Pl sigue siendo muy super-Planckiano frente a la conjetura de distancia, y la coherencia con χ inhomogéneo sigue sin calcular.
- **Con MARGINAL_PEOR o NO_COHERENTE:** la rama eterna sigue siendo la preferida, y la razón queda calculada.

## 6. Expectativa escrita antes de calcular

La tolerancia de C45j es una propiedad del padre (A·H, αH², n_rebotes), no de dónde aterriza el hijo. Los hijos directos son unas cinco veces más numerosos que los eternos en la ventana. Por eso espero que f_v(D) ≥ f_v(E), probablemente RAMA_DIRECTA_VIABLE, y la banda 10–50 POBLADA.

**El riesgo:** que los hijos directos vengan de padres más grandes (A·H mayor) o necesiten más rebotes. Las dos cosas bajan ζ0_max deprisa, como (A·H)⁻³ y como 1/n_rebotes.

## 7. Coste y salidas

- **Coste:**
  - la rejilla de polo de U22b (1 100 integraciones homogéneas en C), del orden de minutos;
  - C45j, unos 2 minutos según su informe.
  
  Cabe en la nube. No necesita CLASS ni redes.
- **Código:** `chain/u27_direct_children.py`.
- **Resultado:** `results/u27_direct_children.json`.
- **Informe:** `docs/informes/INFORME_U27.md`.
- **Contrato:** U27 no cambia ningún veredicto de la cadena, así que no entra en `make chain`. Lo cita la entrada D9/D11 del contrato solo si el autor decide llevarlo a una v3.

## Enmienda 1 (25 de octubre de 2026, 18:53, Madrid): las puertas comparan con los valores archivados completos

*Antes de cualquier código o cálculo de U27. La señaló Codex en la PR 42.*

- **Puerta V1, con el valor redondeado:** §3 daba 0.035907 con una tolerancia de 10⁻⁹. El valor archivado es 0.03590664272890485, es decir, 20 viables de 557 fértiles. Una reproducción exacta no habría pasado la puerta.
  - V1 compara ahora `frac_viable_LH1_over_fertile` con el valor archivado completo, el que lee `data/c45b/c45j_results.json`, a 10⁻⁹. Además exige los mismos enteros: 20 viables y 557 fértiles.
- **Puerta V2, con el mismo problema:** f_ES = 0.152878 también estaba redondeado.
  - V2 compara ahora los enteros (n = 1 100, n_fertile = 886, n_window_fertile = 556, y por separación n_fertile y n_eternal) y f_ES con el valor completo de `results/u22b_births.json` (0.1528776978417266, es decir, 85/556), a 10⁻⁹.
- No cambia nada más: ni las filas, ni los umbrales, ni la lógica del veredicto.

## Enmienda 2 (25 de octubre de 2026, 19:09, Madrid): el código de C45j llegó, y lo que dice antes de ejecutarlo

*El autor envió `c45j_ltb.zip`, copiado sin cambios de `claude_dual/c45b` y con `SHA256SUMS.txt`, y los hashes coinciden. Se leyó el código antes de ejecutar nada. Todavía no se ha calculado nada de U27.*

- **Ficheros:**
  - `ltb/multibounce_viability.py`, con lo que importa: `common/provenance.py`, `common/potential.py` y `ltb/growing_qt.py`;
  - lo que arrastran esos módulos: `ltb/shells.py`, `ltb/growing.py`, `ltb/bank_pbh.py` y `ltb/coherence.py`;
  - `common/__init__.py`.
- **La entrada no es `results/c45h_content_dust*`,** como decía §3, sino `results/c45h_content_20261004T153219Z/` completa. El script lee el último `results/c45h_content_2*` y se queda con las filas de polvo. Su `PROTOCOL.json` confirma que C45j usó esa carpeta.
  - Sus 3 300 filas son idénticas, una a una, a las de `data/c45b/c45h_results.json`, que es la entrada de S1.
- **Dónde va y cómo se ejecuta:** todo va a `components/c45j_ltb/` con las rutas originales, y se ejecuta desde esa carpeta (`python -B -m ltb.multibounce_viability`), porque busca `results/` por ruta relativa. V1 se hace así, con el código tal cual.
- **Lo que dice el código, y que decide la rama de §3.4:**
  - A·H = 0.385·√λ y αH² = sep·4λ/27 solo dependen del padre (λ, sep).
  - ε_est(A·H, αH², L·H) solo usa `H_INF`, la de la meseta de jerarquía de `common/potential.py`, para pasar α a las unidades de LTB.
  - Por tanto, ni A·H, ni αH², ni ε_est dependen de la meseta ni de dónde aterriza el hijo.
  - Se toma la segunda rama de §3.4: con el polo se usan los mismos (A·H, αH²) y ε_est. Solo cambian, con la rejilla de polo, la clasificación (fértil, χ_land, E_phys_max) y n_rebotes.
- **Concreciones de §2 y §4,** que no cambian ningún umbral:
  - **La clase E con la suma (R5)** es la bandera `eternal_sub` de las filas archivadas, con χ_et = 417.28 de C45h, como en C45j. Se comprobó que es exactamente fértil, χ_land > χ_et y E_phys_max < 1, sin ninguna discrepancia en la ventana. La clase D es fértil, E_phys_max < 1 y χ_land ≤ 417.28.
  - **Con el polo (R1–R4),** χ_et = 540.95, el valor de U22b. La rejilla se recalcula con la misma integración que `chain/u22b_births.grid` (modo 4, C = 0.3516), guardando además `chi_abs_max` (`out[1]`) para la banda.
  - **V2** compara el resumen de esa rejilla, hecho con `u22b_births.summary`, con el guardado.
- **Coste previsto:** 22 valores de λ, 3 separaciones y 3 L·H dan como mucho 198 evaluaciones de ε_est, porque los resultados se reutilizan y no dependen de σ. Más la rejilla de polo, unos 50 s. Todo, en la nube.
