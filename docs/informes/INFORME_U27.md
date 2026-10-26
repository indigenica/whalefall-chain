# INFORME_U27 · La coherencia de los hijos directos con la meseta de polo

*25 de octubre de 2026, 19:15 (Madrid). Claude, en la nube.*
- **Preregistro:** `docs/preregistros/ESPEC_U27.md`, escrito a las 18:51 y antes de cualquier cálculo, con dos enmiendas:
  - la enmienda 1 corrige las puertas de validación;
  - la enmienda 2 se escribió al recibir el código de C45j, que se leyó antes de ejecutarlo.
- **Código:** `chain/u27_direct_children.py`, con el código de C45j copiado sin cambios en `components/c45j_ltb/`.
- **Resultado:** `results/u27_direct_children.json`.

## En una frase

**Veredicto preregistrado: RAMA_DIRECTA_VIABLE, y BANDA_10_50_POBLADA.** En el mismo modelo que dio el 3.6 % para los hijos eternos, el 14.9 % de los hijos fértiles son directos y aguantan perturbaciones de amplitud inflacionaria a lo largo de los rebotes del padre (L·H = 1). El 4.9 % de los fértiles son además directos con aterrizaje entre 10 y 50 M_Pl. Para nuestro bolsillo, la rama eterna deja de ser la única opción: con un hijo directo, la meseta solo tiene que ser plana hasta su aterrizaje, no hasta ~10³ M_Pl.

## 1. Qué se calculó

- **El modelo** es el de C45j, sin cambios. El padre se recolapsa en LTB-QT con polvo, y su tolerancia a la desincronización es ε\*_n = ε_est(A·H, αH², L·H)/n_rebotes. La amplitud equivalente de unas perturbaciones de curvatura de amplitud inflacionaria es ε_phys = (2/3)·ζ0·π²·(A·H/L·H)², con ζ0 = 4.6×10⁻⁵. Un hijo es **viable** si ζ0_max ≥ ζ0.
- **ε_est solo depende del padre (λ, sep) y de L·H,** no de la meseta ni de dónde aterriza el hijo (enmienda 2). Con el polo solo cambian la clasificación de cada nacimiento y su número de rebotes.
- **Las clases,** en la ventana sep ∈ {1e-4, 1e-5, 1e-6}, con polvo:
  - **E (eterno):** fértil, sub-Planck y χ_land > χ_et;
  - **D (directo):** fértil, sub-Planck y χ_land ≤ χ_et;
  - **D₁₀₋₅₀:** los D cuyo aterrizaje (máximo de |χ| tras el rebote) está entre 10 y 50 M_Pl.
- Con el polo, χ_et = 540.95 M_Pl.

## 2. Validaciones

| Puerta | Qué compara | Resultado |
|---|---|---|
| **V1** | C45j, ejecutado tal cual, frente a `data/c45b/c45j_results.json` | **PASA, idéntico**: los mismos 91 hijos, diferencia relativa 0 en ε_est y ζ0_max, 20 viables de 557 fértiles (0.03590664272890485) |
| **V2** | la rejilla de polo, frente al resumen guardado de U22b | **PASA**: 1 100 nacimientos, 886 fértiles, 556 en la ventana, los recuentos por separación y f_ES = 85/556 |

Además, R5 (con la suma, clase E) reproduce el 3.6 % de C45j: 20 viables de 557.

## 3. Resultados (meseta de polo; f_v = viables / 556 fértiles de la ventana)

| Fila | Clase | Hijos | L·H = 0.3 | **L·H = 1** | L·H = 3 |
|---|---|---|---|---|---|
| R1 | E | 85 | 6 (1.1 %) | 20 (**3.6 %**) | 36 (6.5 %) |
| **R2** | **D** | 305 | 12 (2.2 %) | 83 (**14.9 %**) | 164 (29.5 %) |
| **R3** | **D₁₀₋₅₀** | 93 | 4 (0.7 %) | 27 (**4.9 %**) | 54 (9.7 %) |

**Fracción condicionada,** es decir, cuántos de cada clase aguantan, con L·H = 1:
- eternos: 24 %;
- directos: 27 %;
- directos de la banda: 29 %.

Por hijo, la tolerancia es prácticamente la misma. Lo que cambia es el número: hay 3.6 veces más hijos directos sub-Planck que eternos.

**Desglose de los directos viables (L·H = 1) por separación:**
- 45 con 1e-4;
- 28 con 1e-5;
- 10 con 1e-6.

En los eternos eran 18, 2 y 0. Los directos aguantan también con las separaciones pequeñas.

**El riesgo previsto no se dio.** Los directos no vienen de padres más grandes, sino de más pequeños: la mediana de A·H es 3.6, frente a 5.3 en los eternos. Tampoco necesitan más rebotes: el 76 % lanza en el primero, frente al 55 % de los eternos.

**Hay 166 fértiles de la ventana que no son E ni D,** porque superan la densidad de Planck en el rebote. Cuentan en el denominador, igual que en C45j.

**Con la suma (R5, informativa),** las cifras son casi iguales:
- E: 20 de 557 (3.6 %);
- D: 82 (14.7 %);
- D₁₀₋₅₀: 31 (5.6 %).

## 4. Veredicto

- **Principal: RAMA_DIRECTA_VIABLE.** f_v(D) = 14.9 % ≥ 5 %, el umbral de ESPEC_C45K para que D11 pase a VIABLE.
- **Banda: BANDA_10_50_POBLADA.** Hay 27 viables, y f_v = 4.9 % ≥ 1 %.
- **Se cumplió la expectativa escrita antes de calcular** (ESPEC §6).

## 5. Lo que esto significa, y lo que no

**Lo que significa:**
- La rama directa ya no es una nota al margen. En el mismo modelo de coherencia, hay cuatro veces más hijos viables directos que eternos, y la banda de 10–50 M_Pl tiene por sí sola más viables que toda la rama eterna (27 frente a 20).
- **D3 (rango de campo):** si nuestro bolsillo es un hijo directo, la meseta solo tiene que ser plana hasta su aterrizaje, típicamente 10–50 M_Pl, en lugar de ~10³ M_Pl. Las predicciones del CMB no cambian, porque los últimos ~52 e-folds son los mismos.
- **D9 (medida):** para un hijo directo, las predicciones no dependen de pesos dentro de la inflación eterna. La fracción de nacimientos de cada rama sigue dependiendo del prior sobre los padres, como ya dice §8 de la v2.

**Lo que no significa:**
- **No resuelve el conflicto con el swampland.** 10–50 M_Pl sigue siendo muy super-Planckiano frente a la conjetura de distancia, que pide un recorrido del orden de 1 M_Pl. El conflicto baja en un factor 20–100, pero no desaparece.
- **No mejora el modelo de coherencia.** Es el mismo modelo de polvo de C45j, y su aviso vale igual. En los hijos viables, el rebote está dominado por la energía cinética de χ, y la coherencia con χ inhomogéneo (C45k) sigue sin calcular.
- **Depende de la escala de la perturbación:**
  - con L·H = 0.3, f_v(D) = 2.2 %, que habría sido MARGINAL;
  - con L·H = 3, sube al 29.5 %.

  El veredicto se fijó con L·H = 1, como el 3.6 % de C45j. Los eternos son igual de sensibles: con L·H = 0.3 bajan al 1.1 %.
- **La huella del rebote sigue fuera de alcance.** Todos estos hijos tienen ≳ 108 e-folds, así que el rebote queda fuera del universo observable, igual que en la v2.

## 6. Qué cambia en el repositorio

- **La cadena y el contrato no cambian:** U27 no está en `make chain`, como decía el ESPEC.
- **La entrada D11 del contrato** («PASS_TOY / MARGINAL») podría citar U27 para la rama directa. Igual que D3 y D9 en el texto, lo decide el autor de cara a una v3 del preprint.
- **El código de C45j** queda en `components/c45j_ltb/`, con su entrada original y la ejecución de V1 (`results/c45j_multibounce_20261025T181201Z/`), y con el sha256 de cada fichero en `components/MANIFEST.json`.

## Reproducir

```bash
python -B chain/u27_direct_children.py          # V1 (C45j tal cual, ~1.5 min) + la rejilla de polo (~1 min) + ε_est (~1 min)
cd components/c45j_ltb && python -B -m ltb.multibounce_viability    # solo C45j
```
