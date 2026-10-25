# U17g · El borde superior de la banda de F sobre el intervalo de X de U19d (244–1 367) · Preregistro

*22 de octubre de 2026, 21:33 (Madrid). Escrito antes de calcular ninguna fila de CLASS de U17g; la hora que cuenta es la de su commit. Se ejecuta en la nube con la tubería de U17f y su caché de CLASS.*

**Lo que ya se ha visto, y se declara:**
- **U17f completo** (`INFORME_U17f.md`): BORDE_0.7_CONFIRMADO; F₀ PERMITIDO CON X ≤ 194. Sus filas llegan hasta Y = 1.98 (X ≈ 320 en F₀).
- **U19d completo** (`INFORME_U19d.md`): X_UNIV_INDECISO, entre 244 y 1 367, con X_univ ≈ 556 y p = 1.35 ± 0.31. Los modelos que ajustan dan 218, 273 y 556.
- **El `--dry` de este script,** hecho hoy antes de escribir este documento. Solo calcula, sin CLASS:
  - las masas, T_osc y λ/λ_WDM;
  - la isocurvatura f_iso, la superradiancia y la ruta de cada fila.

  Por ejemplo, f_iso supera 0.0084 en Y = 4.25 con F ≥ 0.85, y en Y = 8.44 con F ≥ 0.5. **No se ha calculado ningún δA ni ninguna forma de Lyα** de las filas nuevas.

## Motivo

U19d dejó el X del universo entre 244 y 1 367. U17f solo dice dónde está el borde superior de la banda hasta X ≈ 320 (0.7·F₀). Por encima no hay filas, y el valor central de U19d (556) cae fuera.

**U17g calcula el borde con la misma tubería** para X de 405 a 1 367.

## Filas (`lya/u17g_class.py`)

**La población, el número y la rama** son los de U17f, sin cambios:
- la forma medida `measured_shape("low")`;
- B2, q = 0.87, x₀ = 10;
- K = (8π/27)·R_ρ·ξ_B2·ℓ·Y/x̄, con R_ρ = 1.56 y m_G autoconsistente;
- meseta de suma, como la enmienda 1 de U17f.

**Familias primarias,** con X ≈ 162·Y en F₀:

| Y | 2.5 | 3.0 | 3.43 | 4.25 | 5.5 | 8.44 |
|---|---|---|---|---|---|---|
| X en F₀ | ≈ 405 | ≈ 486 | ≈ 556 (X_univ de U19d) | ≈ 689 | ≈ 891 | ≈ 1 367 (X_hi de U19d) |

**Control de la tubería, sin veredicto:** la familia Y = 1.98 de U17f, en F ∈ {0.85, 0.7, 0.5}. Su δA tiene que coincidir con el de U17f a 10⁻³. Si no coincide, el veredicto es INVALIDO y se busca el error.

**Rejilla de F:** la de U17f más 0.6, para resolver el borde entre 0.7 y 0.5: F/F₀ ∈ {1, 0.85, 0.7, **0.6**, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03}.

**Qué pasa por CLASS:**
- **Todas las filas que no cumplen la regla de frío trivial de U21** (λ_proxy/λ_WDM(5.3 keV) < 0.05).
- **No se usa la dominancia de U17f:** U17e no tiene filas con Y > 2.24 con las que comparar.
- Las filas frías triviales pasan Lyα sin CLASS.
- En total, **52 filas primarias más las 3 del control.**

## Puertas por fila (las de U17f, sin cambios)

1. **Lyα,** con la puerta primaria de U12 (`U5.gates`).
2. **Isocurvatura,** f_iso < 0.0084 con C = 1 y r = 2.40.
3. **Radiación oscura:** ΔN_eff ≤ 0.3.
4. **Cuerdas en el CMB.**
5. **Superradiancia estelar** de U23: la fila no puede estar «excluida».

**Se informan también:**
- la cota secundaria de supermasivos;
- LISA;
- r = 1.78;
- el δA de cada fila.

## Veredicto

Para cada familia, F_max es el F más alto de la rejilla que pasa las cinco puertas. **F_max,bajo** es el menor de las seis familias primarias; en principio, el de Y = 8.44.

| Veredicto | Condición |
|---|---|
| **BORDE_0.7_SOSTENIDO** | F_max,bajo ≥ 0.7: el borde de U17f vale en todo el intervalo de U19d |
| **BORDE_0.6** | F_max,bajo = 0.6 |
| **BORDE_0.5** | F_max,bajo = 0.5 |
| **BORDE_BAJA_A_x** | F_max,bajo = x ≤ 0.35 |
| **SIN_BANDA** | alguna familia primaria no pasa en ningún F |
| **INVALIDO** | falla el control de la tubería |

**Se informan también, sin veredicto:**
- el F_max de cada familia;
- **el F_max en Y = 3.43,** el valor central de U19d;
- la banda completa de cada familia (F_min, contigüidad y la superradiancia), y la banda con supermasivos;
- la traducción a la banda combinada: el borde en función de X, de 194 a 1 367, uniendo U17f y U17g.

## Lo que se espera, escrito antes de calcular

- **Lyα empeora con Y:** en U17f, Y = 1.98 pasaba en 0.7·F₀ con δA ≈ 0.32.
  - En Y = 2.5–3.43, **se espera que 0.7·F₀ falle y que pase 0.6 o 0.5·F₀.**
  - En Y = 8.44, el borde podría bajar a 0.35·F₀ o menos.
- **La isocurvatura,** según el `--dry`, recorta F ≥ 0.5 en Y = 8.44 y F ≥ 0.85 en Y = 4.25. No llega por debajo de lo que ya recortaría Lyα.
- **Expectativa: BORDE_0.5 o BORDE_BAJA_A_0.35.** En el valor central de U19d, F_max ≈ 0.5–0.6·F₀.

## Coste

55 filas de CLASS de 20–40 min cada una, con 4 núcleos: ≈ 5–9 h en la nube. El productor salta las filas en caché y se puede relanzar.

## Lo que U17g no resuelve

- **Cuál es el X del universo:** U19d lo acota entre 244 y 1 367, y U17g da el borde en todo ese intervalo.
- **Lo mismo que U17f:** la normalización C de la isocurvatura, la forma de la población (la de ℓ₁ = 4.3) y la superradiancia de supermasivos.
- **Lyα completo,** que es U15d. Solo prueba X ≤ 320.

## Nota · 23 de octubre de 2026, 04:36 (Madrid): el resultado

- **Las 55 filas de CLASS terminaron sin fallos,** en 7.0 h en la nube.
- **El control de la tubería pasa:** δA idéntico al de U17f en las tres filas.
- **`u17g_class.py`: BORDE_BAJA_A_0.35.** F_max = 0.6, 0.6, 0.5, 0.5, 0.35 y 0.35·F₀ de Y = 2.5 a Y = 8.44, y 0.5·F₀ en el X central de U19d. El informe es `docs/informes/INFORME_U17g.md`.


## Nota · 23 de octubre de 2026, 04:50 (Madrid): una comprobación del parche de CLASS dentro del script, añadida tras el cálculo (Codex)

- **El problema:** las filas usan la cuadratura del parche `env/class_qm_trapz_log.patch` (strategy 4), y `u17g_class.py` no comprobaba él mismo que CLASS lo llevara.
  - Antes de lanzar, `env/check_env.py` lo verificó: «ENTORNO LISTO», con la comprobación del parche incluida.
  - Pero un CLASS sin parche no se queja, y la caché no distingue entre compilaciones.
- **El cambio:** `quadrature_patch()` repite la comprobación de `check_env.py` (el fondo ncdm con strategy 4 frente al muestreo por defecto, < 2×10⁻³). Corre antes de cualquier fila, nunca sale de la caché, y si falla el script se detiene.
- **Se volvió a ejecutar U17g desde la caché:** el parche pasa (2.2×10⁻⁵), y el resultado es idéntico salvo el tiempo y el campo nuevo.
- **El contrato** solo publica las bandas de U17g si pasa el control de la tubería (`pipeline_ok`).
- **No cambia:** ninguna fila, ninguna puerta ni el veredicto.
