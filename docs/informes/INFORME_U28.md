# INFORME_U28 · C45k: la coherencia de los rebotes del padre con χ inhomogéneo

*26 de octubre de 2026, 02:51 (Madrid). Claude, en la nube.*
- **Preregistro:** `docs/preregistros/ESPEC_U28.md`, escrito el 25 de octubre a las 21:38, antes de cualquier código de U28, con siete enmiendas fechadas, todas anteriores a la producción:
  - la 1 cambia los umbrales a intervalos y fija la censura;
  - la 2 suspende el solver preregistrado (`c58c/regular_center`), que no puede cubrir el parche;
  - la 3, elegida por el autor (vía A), pasa a un solver comóvil con el polvo y fija la calibración frente al polvo;
  - la 4 hace adiabático el reparto de la familia (a);
  - la 5 pasa a la formulación de Misner–Sharp;
  - la 6 cierra el universo en la «panza», sin atenuador, y fija las reglas de umbral y de clasificación;
  - la 7 corrige un error de código en Klein–Gordon cerca del polo, antes de repetir las puertas.
- **Código:** `chain/u28_comoving.py` (el solver) y `chain/u28_c45k.py` (puertas, producción y análisis).
  - El código de C45b en doble nulo, de donde salen las ecuaciones, queda como copia sin cambios en `components/c45b_dn/`. U28 no lo importa.
  - `c58c/regular_center` no se copió: la enmienda 2 lo descartó antes de escribir código.
- **Resultados:**
  - `results/u28_c45k.json`;
  - las 780 corridas, una por fichero, en `results/u28_runs/`, con `gates.json`, `selection.json` y `produce.out`.

## En una frase

**Veredicto preregistrado: IGUAL en las dos clases** (a L·H = 1).
- **En los directos (D),** el umbral de coherencia con χ inhomogéneo coincide con el del polvo, dentro de la resolución de la rejilla de amplitudes (ρ entre 0.58 y 1.83).
- **En los eternos (E),** sale por debajo (ρ entre 0.16 y 0.49), sin llegar al umbral de PEOR.
- **El tiempo de rebote solo,** en una lectura informativa, coincide con el polvo en las dos clases (§4.1). Lo que rebaja a los eternos es dónde aterrizan las celdas perturbadas (§4.2).

Con la fracción viable corregida:
- **la rama directa sigue VIABLE:** entre el 12.4 % y el 18.0 % de los fértiles, frente al 14.9 % del polvo;
- **la rama eterna pasa a MARGINAL:** entre el 1.8 % y el 2.9 %, frente al 3.6 %;
- **la banda de 10–50 M_Pl** conserva al menos el 4.0 %, frente al 4.9 %.

Se cumple la expectativa escrita antes de calcular: IGUAL o MEJOR, y la rama directa VIABLE. Pero el campo no frena la desincronización, como suponía su argumento (§7): la iguala.

## 1. Qué se calculó

- **El padre** es el de C45 de cada hijo: con λ, sep y σ, la meseta de polo (C = 0.3516), polvo y χ.
  - Se evoluciona en la gauge síncrona comóvil con el polvo, ds² = −dt² + e^{2Λ}dχ² + r²dΩ², sobre medio universo cerrado, con paridad en el polo y espejo en el ecuador.
  - Las variables son r, U = ṙ, Γ = e^{−Λ}r', φ, P (el momento de χ, en volúmenes finitos) y la masa generalizada Y = r³h(ψ).
  - Las ecuaciones son las NPQT esféricas de C45b (Bueno, Cano, Hennigar y Murcia, ec. 33–34), en forma covariante 2D. N = 256 celdas, con RK4 y CFL de 0.25.
- **Datos iniciales:** la rodaja simétrica en el tiempo del rebote anterior (U = 0). Es la sincronización «rebote previo uniforme» de C45d y C45j.
- **Las dos familias de perturbación,** con el sinc de C45j y escala L·H:
  - **(a)** la curvatura, K = 1 + ε·sinc, con reparto adiabático entre campo y polvo;
  - **(b)** el momento del campo, δΠ/Π = ε·sinc, con el polvo sin perturbar.
- **Amplitudes:** ε ∈ {10⁻⁵ … 10⁻²}, ampliadas hacia abajo hasta 10⁻⁷ o hacia arriba hasta 10⁻¹ según la enmienda 1.
- **Una corrida es COHERENTE** si cumple las dos condiciones en la región x ≤ 0.95:
  - la dispersión de los tiempos de rebote es Δt_B < 0.3·√α;
  - al menos el 80 % del volumen (f_loc ≥ 0.8) aterriza en la clase del hijo homogéneo, fértil y por debajo de la densidad de Planck. Se mide con 24 universos separados.
- **ρ = ε_lo/(0.3·ε_est):** compara el umbral del campo con la tolerancia continua del polvo de C45j en el mismo hijo, llevada al criterio de U28. ρ = 1 quiere decir «igual que el polvo».
- **Las filas:**
  - E-v: 10 de los 20 eternos viables;
  - D-v: 10 de los 83 directos viables, 4 de ellos de la banda 10–50;
  - controles: E-n y D-n, 5 hijos cada uno.

  Los grupos -v se corren con L·H = 1 (veredicto) y con L·H = 0.3 (informativo).

## 2. Validaciones (puertas registradas, con el código corregido de la enmienda 7)

| Puerta | Qué compara | Resultado |
|---|---|---|
| **V0** | las ligaduras al doblar N (64 → 128 → 256) | **PASA**: el residuo de la hamiltoniana baja un factor 15.4–16.2 al doblar N (orden ≈ 4; se pedía ≈ 2) |
| **V1** | el límite homogéneo frente a `frw_core.c`, en 5 hijos | **PASA**: t_B con diferencia relativa ≤ 1.4×10⁻⁷ y a en el rebote ≤ 2.1×10⁻⁵ (se pedía 10⁻³); la misma rama |
| **V2** | el límite de solo polvo frente a LTB-QT de C45j (`components/c45j_ltb`), en 3 padres | **PASA**: t_B capa a capa ≤ 1.4×10⁻⁸ |
| **V3** | N = 256 frente a 512, en 4 corridas (2 E y 2 D, ε = 10⁻³) | **PASA**: Δt_B cambia como mucho un 4.6 %, y ninguna etiqueta cambia |
| **V4** | Oppenheimer–Snyder en el límite GR | **PASA**: ≤ 8.9×10⁻¹³ |

La salida de las puertas anteriores a la enmienda 7 queda en `results/_superseded/u28_gates_pre_amendment7/`, con su nota.

## 3. Resultados (meseta de polo; L·H = 1)

### 3.1 Veredicto por clase

| Clase | Hijos | Mediana de ρ_lo | Mediana de ρ_hi | Corridas inválidas | **Veredicto** |
|---|---|---|---|---|---|
| **E** (E-v) | 10 | 0.16 | 0.49 | 26/144 (18.1 %) | **IGUAL** |
| **D** (D-v) | 10 | 0.58 | 1.83 | 27/152 (17.8 %) | **IGUAL** |

- **En D, el intervalo [0.58, 1.83] contiene ρ = 1.** La rejilla va en pasos de un factor ≈ 3, así que no se puede distinguir del polvo.
- **En E, el intervalo [0.16, 0.49] queda por debajo de 1,** pero no llega a PEOR, que exige ρ_hi ≤ 1/3. La clase eterna sale peor que el polvo, sin cruzar el umbral preregistrado.
- **Ningún hijo de E-v ni de D-v tiene el umbral censurado:** todos tienen una amplitud coherente y una incoherente.
- **Las inválidas se quedan por debajo del 20 %** que haría la clase INCONCLUSA. 50 de las 53 son de ε ≥ 3×10⁻³: cruces de capas o pasos de tiempo que se hunden, siempre por encima del umbral de coherencia de su hijo.

### 3.2 Fracción viable corregida y D11 (sobre los 556 fértiles de U27)

| Rama | Polvo (U27) | Con χ: f_lo | Con χ: f_hi | **D11 (con f_lo)** | Con f_hi |
|---|---|---|---|---|---|
| **D (directa)** | 83 (14.9 %) | 69 (**12.4 %**) | 100 (18.0 %) | **VIABLE** | VIABLE |
| **D₁₀₋₅₀** | 27 (4.9 %) | 22 (**4.0 %**) | — | (sin umbral propio) | — |
| **E (eterna)** | 20 (3.6 %) | 10 (**1.8 %**) | 16 (2.9 %) | **MARGINAL** | MARGINAL |

- Los veredictos con f_lo y con f_hi coinciden en las dos ramas.
- **Por qué baja f_lo:** sale de la cota inferior ρ_lo, que la rejilla redondea hacia abajo. En la rama D, el valor del polvo (14.9 %) cae dentro de [f_lo, f_hi]. En la E, no: queda por encima de f_hi.

### 3.3 Por hijo (L·H = 1; la familia es la más desfavorable)

| Grupo | Hijo (sep, λ, σ) | A·H | χ_land | f_k | Familia | [ε_lo, ε_hi] | ρ (preregistrado) | ρ (solo Δt_B) |
|---|---|---|---|---|---|---|---|---|
| E-v | 1e-4, 18.39, 3.1e-4 | 1.65 | 2628 | 2×10⁻⁴ | a | [3e-5, 1e-4] | [0.058, 0.19] | [0.58, 1.9] |
| E-v | 1e-4, 12.43, 3.1e-5 | 1.36 | 594 | 2×10⁻⁵ | a | [1e-4, 3e-4] | [0.19, 0.56] | [0.56, 1.9] |
| E-v | 1e-4, 18.39, 0.3 | 1.65 | 1183 | 1.00 | a | [1e-4, 3e-4] | [0.19, 0.58] | [0.19, 0.58] |
| E-v | 1e-5, 5.68, 9.7e-3 | 0.92 | 636 | 0.99 | b | [3e-5, 1e-4] | [0.13, 0.43] | [4.3, 13] |
| E-v | 1e-4, 5.68, 3.1e-5 | 0.92 | 1457 | 3×10⁻⁵ | a | [3e-4, 1e-3] | [0.40, 1.3] | [0.40, 1.3] |
| E-v | 1e-4, 3.84, 9.8e-4 | 0.75 | 1790 | 3×10⁻⁴ | a | [3e-4, 1e-3] | [0.31, 1.0] | [0.31, 1.0] |
| E-v | 1e-4, 2.59, 0.095 | 0.62 | 1007 | 1.00 | a | [1e-3, 3e-3] | [0.76, 2.3] | [0.76, 2.3] |
| E-v | 1e-4, 1.75, 0.095 | 0.51 | 1063 | 1.00 | b | [1e-4, 3e-4] | [0.054, 0.16] | [0.54, 1.6] |
| E-v | 1e-4, 1.18, 1e-5 | 0.42 | 608 | 4×10⁻⁶ | a | [3e-5, 1e-4] | [0.011, 0.038] | [0.38, 1.1] |
| E-v | 1e-4, 1.18, 0.095 | 0.42 | 611 | 1.00 | b | [3e-4, 1e-3] | [0.11, 0.38] | [1.1, 3.8] |
| D-v | 1e-5, 12.43, 3.1e-4 (banda) | 1.36 | 16 | 0.05 | a | [1e-4, 3e-4] | [0.59, 1.8] | [0.59, 1.8] |
| D-v | 1e-4, 12.43, 9.7e-3 (banda) | 1.36 | 49 | 0.61 | a | [3e-4, 1e-3] | [0.56, 1.9] | [0.56, 1.9] |
| D-v | 1e-5, 2.59, 3.1e-4 (banda) | 0.62 | 27 | 0.13 | a | [3e-4, 1e-3] | [0.72, 2.4] | [0.72, 2.4] |
| D-v | 1e-4, 1.18, 3.1e-3 (banda) | 0.42 | 24 | 0.27 | a | [3e-3, 1e-2] | [1.1, 3.8] | [1.1, 3.8] |
| D-v | 1e-5, 12.43, 3.1e-3 | 1.36 | 154 | 0.84 | a | [1e-4, 3e-4] | [0.59, 1.8] | [0.59, 1.8] |
| D-v | 1e-4, 18.39, 0.030 | 1.65 | 54 | 0.67 | a | [3e-4, 1e-3] | [0.58, 1.9] | [0.58, 1.9] |
| D-v | 1e-4, 5.68, 1e-5 | 0.92 | 289 | 3×10⁻⁶ | a | [3e-4, 1e-3] | [0.40, 1.3] | [0.40, 1.3] |
| D-v | 1e-4, 5.68, 0.030 | 0.92 | 202 | 0.96 | b | [1e-3, 3e-3] | [1.3, 4.0] | [4.0, 13] |
| D-v | 1e-5, 1.75, 3.1e-3 | 0.51 | 336 | 0.96 | b | [1e-5, 3e-5] | [0.017, 0.051] | [0.51, 1.7] |
| D-v | 1e-4, 1.18, 0.030 | 0.42 | 231 | 0.97 | a | [1e-3, 3e-3] | [0.38, 1.1] | [1.1, 3.8] |

- **La familia más desfavorable es casi siempre la (a), la curvatura,** como en el modelo de polvo. La (b), el momento del campo, manda en 5 de los 20 hijos, todos con f_k ≥ 0.96.
- **El conjunto coherente no es monótono** en 5 de los 60 pares hijo–familia con L·H = 1, y en 3 más con L·H = 0.3. Se informa, como pide la enmienda 6, y ε_lo es la mayor amplitud coherente.

### 3.4 Controles y L·H = 0.3

- **Controles (E-n y D-n):** ninguno queda rescatado a amplitud inflacionaria (`controls_rescued` está vacío). Con el criterio de Δt_B solo, sus medianas de ρ son [0.57, 1.9], las mismas que en los grupos -v.
- **L·H = 0.3 (informativo):** las dos clases salen **INCONCLUSA** por la regla de §3, con el 39.6 % de inválidas en E y el 30.4 % en D. A esa escala hay inválidas también a amplitudes pequeñas (22 en E-v y 6 en D-v con ε < 3×10⁻³). Las medianas de ρ_lo (0.19 en E y 0.58 en D) coinciden con las de L·H = 1, pero no cuentan.

## 4. Lecturas informativas, decididas después de ver las primeras corridas

*No son preregistradas ni cambian el veredicto. Se decidieron al ver el primer hijo E-v, donde Δt_B igualaba al polvo y aun así la corrida era incoherente. Están en `informative_post_hoc` de `results/u28_c45k.json`.*

### 4.1 El campo iguala al polvo en el tiempo de rebote

Con el criterio de Δt_B solo, comparable directamente con C45j, la mediana de ρ es **[0.55, 1.75] en E-v** y **[0.59, 1.90] en D-v**. En los controles también es [0.57, 1.9]. Todos contienen ρ = 1.

Hijo a hijo, el campo cambia cosas, pero en los dos sentidos:
- **Donde domina el polvo** (f_k ≲ 10⁻³), Δt_B es el del polvo. En el primer hijo E-v, por ejemplo, Δt_B/√α = 0.0576 frente a 0.0575.
- **Donde domina el campo** (f_k ≳ 0.9), ρ va de [0.19, 0.58] a [4.3, 13]. En unos hijos la presión de los gradientes frena la desincronización, y en otros la acelera.

### 4.2 Lo que rebaja la clase E es el aterrizaje, no la desincronización

- **En E-v,** 30 de las 42 corridas válidas incoherentes pasan la prueba de Δt_B y fallan solo en f_loc. En D-v son 13 de 27.
- **Las celdas muestreadas de esas 30 corridas de E-v aterrizan así:** 335 como «ninguna» (no fértiles o por encima de Planck; el clasificador no los separa), 325 como D y solo 60 como E.
  - Un hijo eterno con una perturbación que el rebote tolera ya no aterriza, en buena parte de su volumen, más allá de χ_et = 540.95.
  - Pasa tanto en hijos con χ_land cerca de χ_et (608 y 611) como en hijos lejos de él (1063 y 2628). No es solo un cambio de lado en el borde.
- **Si se acepta aterrizar fértil en cualquiera de las dos clases,** la mediana de E-v solo sube a [0.21, 0.67]. En la mayoría de esas corridas, más del 20 % del volumen aterriza como «ninguna».
- **En la clase D,** el aterrizaje no cambia la mediana: la preregistrada y la de Δt_B solo coinciden. En las 13 corridas de D-v que fallan por f_loc, las celdas aterrizan como «ninguna» (157), como E (141) y solo 14 como D: los cambios de clase van en los dos sentidos.

### 4.3 La fracción cinética del campo en el rebote

ESPEC §1 daba, con el dato de C45j, «entre 46 y 4 200 veces el polvo» en la energía cinética de χ en el rebote. En la rodaja simétrica de U28, la fracción f_k = (χ̇²/2)/ρ varía mucho de un hijo a otro:

| Grupo | mín. | cuartil 1 | mediana | cuartil 3 | máx. |
|---|---|---|---|---|---|
| E-v | 4×10⁻⁶ | 8×10⁻⁵ | 0.49 | 0.998 | 0.999 |
| D-v | 3×10⁻⁶ | 0.17 | 0.64 | 0.93 | 0.97 |
| E-n | 2×10⁻⁸ | 3×10⁻⁶ | 0.97 | 0.996 | 0.999 |
| D-n | 3×10⁻⁸ | 0.04 | 0.74 | 0.89 | 0.91 |

La selección cubre, por tanto, tanto padres dominados por el polvo como padres dominados por el campo. El veredicto IGUAL no depende de que unos u otros falten.

## 5. Qué cambia

- **D11 (varios rebotes) deja de ser solo una estimación de polvo en la rama directa:** calculada con χ inhomogéneo, sigue VIABLE (≥ 12.4 % de los fértiles). La banda de 10–50 M_Pl conserva al menos el 4.0 %.
  - Con U27, eso sostiene la rama directa en una v3 del preprint.
  - D3 (el rango de campo de 10–50 M_Pl para un hijo directo) y D9 (la medida) siguen como los dejó U27.
- **La rama eterna queda MARGINAL** (1.8–2.9 %). Lo que la rebaja es sobre todo dónde aterrizan las celdas perturbadas (§4.2). Ya era la menos favorable en U27.
- **Lo que no se probó:**
  - una sola escala de veredicto (L·H = 1); con L·H = 0.3 las clases salen INCONCLUSA;
  - un solo rebote, extrapolado a varios como en C45j (ε/n_rebotes);
  - simetría esférica;
  - el clasificador de universos separados, que no distingue «no fértil» de «por encima de Planck» (§4.2).
- **La cadena y el contrato no cambian** (ESPEC §8). La línea D11 del contrato se actualizará si el autor lleva U28 a la v3, y llevarlo lo decide el autor.

---

*Nota (26 de octubre de 2026, 09:42, Madrid).* El autor llevó U27 y U28 a la v3 del preprint (`paper/zenodo_whale/main-v3.tex`). La línea D11 del contrato (eslabón 2) recoge ahora este veredicto, como preveía §8 del ESPEC. El texto de arriba no cambia.
