# U22b · Con la inflación de polo, el nacimiento eterno es igual de frecuente: el parche D2 se puede quitar

*16 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U22b.md` (commit `33f76d3`, anterior a la ejecución).
- *Código:* `chain/u22b_births.py`. *Resultado:* `results/u22b_births.json`. En la nube, en 50 s.

## 0. En una frase

**Veredicto preregistrado: NACIMIENTO_ROBUSTO,** con R = f_ES(polo)/f_ES(suma) = **0.94**.

En las 1100 filas de C45h (frac_rad = 0), la fracción de nacimientos fértiles que acaban eternos es:
- **con la suma de 10 niveles:** 0.163, con umbral eterno χ_et = 416 M_Pl;
- **con el polo:** 0.153, con χ_et = 541 M_Pl.

El modo 0 reproduce lo guardado: diferencia máxima 0.0e+00.

## 1. Por separación de escalas

| sep | fértiles (suma) | eternos (suma) | fértiles (polo) | eternos (polo) |
|---|---|---|---|---|
| 1e-07 | 205 | 6 | 203 | 5 |
| 1e-06 | 190 | 19 | 194 | 17 |
| 1e-05 | 187 | 34 | 188 | 33 |
| 0.0001 | 180 | 38 | 174 | 35 |
| 0.001 | 138 | 24 | 127 | 17 |

| | suma | polo |
|---|---|---|
| fértiles en total | 900 | 886 |
| mediana de χ_land entre los fértiles | 341 | 284 |

## 2. Qué significa

- **El nacimiento no depende estadísticamente de la forma del potencial.** Lo que cambia es qué combinación concreta de condiciones iniciales da un hijo eterno: el aterrizaje es muy sensible al potencial caso a caso, y por eso U22, con 20 hijos concretos, dio POLO_CAMBIA_NACIMIENTO. La fracción no cambia.
- **Con U22 (CMB equivalente) y U22b (nacimiento equivalente),** la meseta puede ser la de la **inflación de polo de orden 4**, un lagrangiano estándar con protección radiativa, en lugar de la suma postulada de 10 niveles. **El parche D2 se puede quitar.**
- **La expectativa preregistrada** (R ≈ 0.2, NACIMIENTO_REDUCIDO) **no se cumplió.** Venía de una cuenta errónea en U22 (Nota 1 de INFORME_U22).

## 3. Lo que U22b no resuelve

- **El filtro de perturbaciones de C45j** (ζ₀,max, la «fracción viable» de 3.6 %) no se puede rehacer sin el código de C45j. f_ES es la fracción antes de ese filtro.
- **La conjetura de la distancia de la Ciénaga** sigue sin evitarse.
- **Adoptar el polo como potencial de la cadena** (cambiar S1–S2 por defecto) es una decisión del autor: S2 cambiaría en la cuarta cifra de n_s.
