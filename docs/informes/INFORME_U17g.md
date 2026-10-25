# U17g · El borde superior de la banda sobre el intervalo de U19d: BORDE_BAJA_A_0.35; 0.5·F₀ en el X central

*23 de octubre de 2026, 04:36 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17g.md`.
- *Código:* `lya/u17g_class.py` (la tubería de U17f, con las mismas puertas).
- *Resultado:* `lya/results/u17g/U17G_RESULT.json`; el registro, en `u17g_run.log`.
- *Ejecución:* 55 filas de CLASS en la nube, con 4 núcleos, en 7.0 h, el 22 y 23 de octubre. Meseta de suma, como U17f.

## 0. En una frase

**Veredicto preregistrado: BORDE_BAJA_A_0.35.** Con la población del colapso medida en U19b, la banda de F baja según el X del universo:
- 0.6·F₀ con X ≈ 400–490;
- **0.5·F₀ con X ≈ 556,** el valor central de U19d, y hasta X ≈ 690;
- 0.35·F₀ con X ≈ 890–1 370.

En todas las familias, quien pone el borde superior es Lyα. **El control de la tubería reproduce U17f exactamente.**

«Por forma» marca las filas que pasan Lyα con δA por encima de la referencia ancha (0.336), solo por la forma del espectro, como en U17f.

## 1. Controles

- **La tubería:** la familia Y = 1.98 de U17f da δA = 0.4209, 0.3228 y 0.1770 en F = 0.85, 0.7 y 0.5·F₀, idénticos a los de U17f.
- **Las referencias de la puerta de cribado** son las de U17f: estricta 0.234 y ancha 0.336.

## 2. Por familia

| Y | X en F₀ | F_max (cinco puertas) | F_max solo con Lyα | Borde inferior contiguo | F_max con supermasivos |
|---|---|---|---|---|---|
| 2.50 | ≈ 405 | **0.6** | 0.6 | 0.0625 | 0.5 |
| 3.00 | ≈ 486 | **0.6** | 0.6 | 0.0625 | 0.5 |
| 3.43 | ≈ 556 (central de U19d) | **0.5** | 0.5 | 0.0625 | 0.35 |
| 4.25 | ≈ 689 | **0.5** | 0.5 | 0.045 | 0.35 |
| 5.50 | ≈ 891 | **0.35** | 0.35 | 0.045 | 0.35 |
| 8.44 | ≈ 1 367 (X_hi de U19d) | **0.35** | 0.35 | 0.03 | 0.25 |

**δA en el borde** (primera fila que falla → última que pasa):

| Y | Falla | Pasa |
|---|---|---|
| 2.50 | 0.7·F₀, δA = 0.381 | 0.6·F₀, δA = 0.301 |
| 3.00 | 0.7·F₀, δA = 0.427 | 0.6·F₀, δA = 0.350 (por forma) |
| 3.43 | 0.6·F₀, δA = 0.382 | 0.5·F₀, δA = 0.288 |
| 4.25 | 0.6·F₀, δA = 0.437 | 0.5·F₀, δA = 0.346 (por forma) |
| 5.50 | 0.5·F₀, δA = 0.408 | 0.35·F₀, δA = 0.231 |
| 8.44 | 0.5·F₀, δA = 0.524 | 0.35·F₀, δA = 0.336 (por forma) |

**Las demás puertas:**
- **La isocurvatura** solo falla donde Lyα ya falla: F ≥ 0.85 con Y ≥ 3, y F ≥ 0.6 con Y = 8.44.
- **La radiación oscura y las cuerdas del CMB** no recortan nada.
- **La superradiancia estelar** pone el borde inferior, que baja con Y porque m_G baja: de 0.09·F₀ en U17f a 0.0625–0.03·F₀.

## 3. La banda en todo el intervalo de X (U17f + U17g)

| X del universo | Borde superior (cinco puertas) | Con supermasivos | Origen |
|---|---|---|---|
| ≲ 194 | F₀ (solo por forma) | 0.7·F₀ | U17f; excluido por U19d |
| ≈ 220–275 | 0.85·F₀ | 0.5–0.7·F₀ | U17f |
| ≈ 320 | 0.7·F₀ | 0.5·F₀ | U17f |
| ≈ 400–490 | 0.6·F₀ | 0.5·F₀ | U17g |
| ≈ 556–690 | **0.5·F₀** | 0.35·F₀ | U17g |
| ≈ 890–1 370 | 0.35·F₀ | 0.25–0.35·F₀ | U17g |

**Con el intervalo de U19d** (244–1 367, central 556), el borde superior queda entre 0.85·F₀ y 0.35·F₀. **En el valor central es 0.5·F₀.**

Según los modelos de extrapolación que ajustan los datos de U19d:

| Modelo | X | Borde superior |
|---|---|---|
| ∝ ℓ₁ | ≈ 218 | 0.85·F₀ |
| Lineal | ≈ 273 | 0.85·F₀ |
| Potencia libre | ≈ 556 | 0.5·F₀ |

## 4. Lo que se esperaba

El ESPEC esperaba BORDE_0.5 o BORDE_BAJA_A_0.35, con F_max ≈ 0.5–0.6·F₀ en el X central. Ha salido BORDE_BAJA_A_0.35, con 0.5·F₀ en el X central.

## 5. Qué significa

1. **La rama viva sigue permitida en todo el intervalo de U19d.** Ninguna familia queda sin banda.
2. **El borde superior depende de la extrapolación del colapso:** va de 0.85·F₀ (∝ ℓ₁ o lineal) a 0.5·F₀ (la potencia central) y a 0.35·F₀ en el extremo alto. **La mejor estimación preregistrada (la potencia de U19d) da 0.5·F₀.**
3. **La banda combinada, la estimación más prudente:**
   - **F ≈ 0.03–0.85·F₀** en todo el intervalo de X;
   - **≈ 0.06–0.5·F₀** con el X central;
   - con la cota secundaria de supermasivos, ≈ 0.06–0.35·F₀ en el X central.
4. **F₀ y LISA siguen fuera** (U19d). Ahora hay además un argumento de Lyα: con el X central, F₀ tiene δA = 0.665, muy por encima de cualquier referencia.
5. **U15d,** con Lyα completo en X ≤ 320, decidirá si la puerta de cribado es conservadora en el borde. Si lo es, como en U15c, los bordes de esta tabla podrían subir. Eso pediría una U15e con la verosimilitud completa en X ≈ 556.

## 6. Lo que U17g no resuelve

- **Cuál es el X del universo:** U19d lo deja en 244–1 367, y la forma de la extrapolación sigue abierta.
- **Lyα completo con X > 320:** U15d solo llega a 320.
- **Lo mismo que U17f:** la normalización C de la isocurvatura (C = 1), la forma de la población (ℓ₁ = 4.3) y la superradiancia de supermasivos.

## Nota · 24 de octubre de 2026, 09:24 (Madrid): lo que dijo U15d

- **El punto 5 de §5 ya tiene respuesta.** U15d (`INFORME_U15d.md`) muestra que la puerta de cribado **no es conservadora** en las filas que pasan solo por forma.
  - Lyα completo excluye F₀ con X ≈ 193 (δA = 0.378), que aquí y en U17f pasaba por forma.
  - Lectura informativa: el corte está entre δA = 0.32 y 0.38.
- **Las tres filas «por forma» de §2** (0.6·F₀ con Y = 3.0, 0.5·F₀ con Y = 4.25 y 0.35·F₀ con Y = 8.44, con δA de 0.336 a 0.350) **quedan dudosas.** Sus bordes podrían bajar un escalón.
- **No cambian:**
  - el veredicto preregistrado de U17g;
  - el borde en el X central (0.5·F₀, δA = 0.288);
  - las filas con δA ≤ 0.32.
- **Una U15e** no haría subir los bordes, como suponía §5. Lo que decidiría es si esas tres filas bajan.
