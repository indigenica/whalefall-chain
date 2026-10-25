# U19d · Más semillas por masa para el X₀(ℓ₁) de U19c · Preregistro

*21 de octubre de 2026, 10:41 (Madrid). Escrito antes de producir ninguna red de U19d; la hora que cuenta es la de su commit. Las masas pequeñas se producen en la nube con `kern.c` y las grandes, en una A100 de Colab con el núcleo JAX.*

**Lo que ya se ha visto, y se declara:**
- **U19c completo** (`docs/informes/INFORME_U19c.md`): las 16 redes, sus X_lat y el veredicto **INDETERMINADO**, con σ_p = 0.76 > 0.5.
  - X₀ = 6.13, 6.03, 8.20 y 7.68 en ℓ₁ = 3.61, 4.01, 4.30 y 4.42;
  - p = 1.49 y χ²/ν = 7.2.
- **Las lecturas informativas de U19c §5,** hechas después de ver los datos. La más importante es que la potencia ajustada al X_lat medio por masa da p = 0.42 ± 0.24, con χ²/ν = 1.1.
- **La simulación de la sección «Lo que se espera»,** hecha hoy con la dispersión entre redes de U19c antes de escribir este documento.
- **No se ha producido ninguna red con las semillas 93–96.**

## Motivo

U19c no fijó el exponente p de X₀ ∝ ℓ₁^p. El error de X₀ por masa salía de un ajuste de 4 redes con 2 grados de libertad, y quedó muy por debajo de la dispersión entre masas (χ²/ν = 7.2). Hay dos lecturas, y U19c no las distingue:
- **(i) La dispersión es ruido entre redes,** y los errores de 2 grados de libertad la subestimaban. En ese caso, más semillas por masa la reducen y fijan p.
- **(ii) X₀(ℓ₁) tiene una estructura que no es una potencia,** por ejemplo el cambio de signo de y_s entre 384³ y 512³. Con más semillas, χ²/ν seguiría alto y el veredicto, INDETERMINADO.

U19d triplica el número de redes por masa con el mismo protocolo, y deja que los datos elijan entre (i) y (ii) con el criterio de U19c, sin cambiarlo.

## Las redes

Cuatro semillas más, **93, 94, 95 y 96**, en los dos protocolos y las cuatro masas de U19c: **8 redes nuevas por masa, 32 en total.**

| ℓ_G | N | Protocolo | Dónde | Motor | Coste estimado |
|---|---|---|---|---|---|
| 2.507 | 256³ | el de U19c (`u19c_produce.network`) | nube | `kern.c` | 8 × 0.75 h |
| 2.912 | 384³ | el de U19c | nube | `kern.c` | 8 × 2.4 h |
| 3.2 | 512³ | **el de U19b** (`u19b_produce.network`), el de sus 4 redes primarias | A100 | JAX `evolve` | 8 × ≈ 1.6 h |
| 3.318 | 576³ | el de U19c | A100 | JAX `evolve_lowmem` | 8 × 2.35 h |

- **Las validaciones (V1–V5) de cada masa son las de U19c y U19b,** y ya pasaron. No se repiten.
- **512³ usa el productor de U19b** para que sus redes nuevas tengan exactamente el protocolo de las cuatro que ya existen.
- **En la A100 rigen las condiciones de la enmienda 2 de U19c:**
  - (a) en cada sesión, antes de producir, `test_jax_kernel.py` tiene que dar `ALL PASSED`, y se guarda la salida;
  - (b) la red de control de U19c (ξ idéntico a `kern.c`) sigue valiendo, porque el código no cambia.
- **Cada red guarda `backend`.** Se escribe en `network/results/u19d/runs/` (o en `U19D_OUT`, en Drive).

## Análisis (`network/u19d_analyze.py`)

**Las redes:** las de U19c (en 512³, las 4 primarias de U19b) más las de U19d, por masa.
- **Si no se completan las 32,** por créditos o por tiempo, el análisis usa en todas las masas solo las semillas completas en las cuatro (los dos protocolos). Así ninguna masa entra con más semillas que otra por haber acabado antes.

**Lo que se calcula, sin cambios respecto a U19c:**
- por red: `u19b_analyze.per_run`;
- por masa: la descomposición afín X = X₀ + y_s·Z (`u19c_analyze.affine`);
- el ajuste ponderado de ln X₀ frente a ln ℓ₁, con pivote, y σ_p escalado por √(χ²/ν) si χ²/ν > 1;
- X_univ = X₀(90) + 162·max(y_s + 2σ, 0), con el y_s más desfavorable;
- la comparación de modelos (constante, ∝ ℓ₁, lineal y potencia libre).

## Veredicto (el de U19c, sin cambios)

| Veredicto | Condición |
|---|---|
| **X_UNIV_BAJO** | X_hi ≤ 194 |
| **X_UNIV_MEDIO** | 194 < X_lo y X_hi ≤ 275 |
| **X_UNIV_ALTO** | X_lo > 275 |
| **X_UNIV_INDECISO** | el intervalo cruza 194 o 275 |
| **INDETERMINADO** | menos de 3 masas válidas, falta la de 576³, o σ_p > 0.5 |

**La traducción a la banda de F** es la de U17f, como en U19c.

**Se informan también, sin entrar en el veredicto:**
- el veredicto solo con las semillas nuevas (93–96), como prueba de estabilidad;
- el ajuste al X_lat medio por masa (la lectura de U19c §5), que **se eligió después de ver U19c** y no tiene estatus de veredicto;
- p quitando cada masa por turno.

## Lo que se espera, escrito antes de calcular

**Una simulación** con las Z de U19c, ruido gaussiano entre redes de σ = 0.36 (la dispersión de los residuos afines de U19c), X₀ = 7·(ℓ₁/4)^p, y_s = 0.05 y 2 000 repeticiones por caso, con el análisis de U19c. Frecuencia de cada veredicto:

| p real | Redes por protocolo y masa | BAJO | MEDIO | ALTO | INDECISO | INDETERMINADO | σ_p mediano |
|---|---|---|---|---|---|---|---|
| 0.4 | 2 (U19c) | 0.60 | 0 | 0.01 | 0.21 | 0.18 | 0.35 |
| 0.4 | **6 (U19d)** | **0.93** | 0 | 0 | 0.06 | 0.01 | 0.22 |
| 1.0 | 2 | 0.17 | 0 | 0.13 | 0.51 | 0.19 | 0.35 |
| 1.0 | **6** | 0.17 | 0 | 0.05 | **0.77** | 0.01 | 0.22 |
| 1.5 | 2 | 0.02 | 0 | 0.42 | 0.35 | 0.21 | 0.36 |
| 1.5 | **6** | 0 | 0 | **0.70** | 0.29 | 0.01 | 0.22 |

**Lo que dice la simulación:**
- **Si la lectura (i) es la buena, INDETERMINADO casi desaparece** (≈ 1 %). Queda:
  - **BAJO** si p ≲ 0.5: F₀ vuelve, solo por forma;
  - **ALTO** si p ≳ 1.4: el borde se queda en 0.7·F₀;
  - **INDECISO** si p ≈ 1, porque X ≈ 200 queda justo en la frontera de 194.
- **MEDIO es prácticamente inalcanzable** con esta palanca.
- **Si la lectura (ii) es la buena,** χ²/ν seguirá por encima de ≈ 3 y el veredicto seguirá siendo INDETERMINADO, ahora con muchas más redes. Eso también es un resultado: X₀(ℓ₁) no se puede extrapolar con una potencia desde este intervalo.
- **El resultado más probable** depende de p, que es justo lo que no se sabe. Con la lectura informativa de U19c (p ≈ 0.4), sería BAJO; con la extrapolación de U17 (p ≈ 1), INDECISO.

## Coste

- **En la nube:** ≈ 25 h (256³ y 384³), sin créditos. El productor salta las redes hechas.
- **En la A100 de 40 GB:** ≈ 32 h (512³ y 576³). Con el ritmo de U19c (≈ 6.4 unidades por hora), son ≈ 200 unidades.
  - Si los créditos no llegan, se aplica la regla de las semillas completas.
  - El orden de producción es por semilla (93 en las dos masas, después 94…), para que la regla pierda lo menos posible.

## Lo que U19d no resuelve

- **Lo mismo que U19c:** la extrapolación de ℓ₁ ≈ 4.4 a 90; R_ρ y la ξ de B2; la parte dura del colapso; la caja menor que el horizonte.
- **Una palanca mayor en ℓ₁:** por encima de ℓ_G ≈ 3.35 falla V2.

## Nota · 22 de octubre de 2026, 20:45 (Madrid): el resultado

- **Las 32 redes están completas,** y las cuatro semillas nuevas también en las cuatro masas: la regla de las semillas completas no quita ninguna.
  - Nube: cuatro reinicios del contenedor, anotados en `produce.out`; cada red cortada se repitió entera.
  - A100: dos sesiones en una A100 de 80 GB, con `ALL PASSED` en las dos (`network/results/u19d/a100/`).
- **`u19d_analyze.py`: X_UNIV_INDECISO (entre 244 y 1 367),** con p = 1.35 ± 0.31 y χ²/ν = 0.007. El informe es `docs/informes/INFORME_U19d.md`.

