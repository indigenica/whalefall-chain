# U19c · La dependencia con ℓ del rendimiento blando del colapso de RP³, con el mismo protocolo en cuatro masas · Preregistro

*20 de octubre de 2026, 10:51 (Madrid). Escrito antes de producir ninguna red de U19c y antes de tocar el kernel para la variante de memoria reducida; la hora que cuenta es la de su commit. Las masas pequeñas se producen en la nube y la grande, en WSL.*

**Lo que ya se ha visto, y se declara:**
- **U19b completo**, con sus 12 redes. La descomposición afín X = X₀ + y_s·Z da:
  - X₀ = 8.20 ± 0.30 con ℓ_G = 3.2;
  - X₀ = 6.43 ± 0.10 con ℓ_G = 2.4, pero con un protocolo distinto;
  - una ley de potencias con exponente 1.18 ± 0.2.
- **U17f completo:** el borde superior de la banda según el X del universo es F₀ con X ≲ 194, 0.85·F₀ con X ≲ 275 y 0.7·F₀ con X ≈ 320.
- **U15d está preregistrado** y sus tablas se están calculando. No se ha visto ninguna cadena.
- **`docs/NOTA_TERCERA_MASA_U19.md`:** la geometría y la memoria. Hoy se ha medido además que el pico de memoria de una red de U19b es ≈ 2.2 veces el estado (f, p y fuerza), con una red de 128³ y 12 instantáneas.
- **Ninguna red de las masas nuevas existe todavía.**

## Motivo

La corrección de §7 de U19b dejó el X del universo entre ≈ 170 y ≈ 320, porque dos masas no distinguen las extrapolaciones. La auditoría señaló además que, entre ℓ_G = 2.4 y 3.2, cambiaban a la vez varias cosas:
- la fase de núcleo congelado;
- la caja en volúmenes de Hubble;
- el comienzo de las redes densas;
- el final del drenaje.

U17f mostró que esa extrapolación decide el borde superior de la banda, y en particular si vuelve F₀ (X ≲ 194).

**U19c mide X₀(ℓ₁) en cuatro masas con el protocolo de ℓ_G = 3.2 trasladado a cada una,** de modo que solo cambien ℓ_G y, con él, ℓ₁.

## El protocolo, relativo a ℓ₁ (el de las redes primarias de U19b con ℓ_G = 3.2)

- **La masa fija la caja:** N = 2τ_A y ℓ(τ_A) − ℓ₁ = ln 256 − 3.2 − ln 3 = 1.247, igual en todas. Es decir, **ℓ_G = 3.2 + ln(N/512)**.
- **El final:** τ_B = 3τ_A, de modo que todas acaban en m_G/H ≈ 94. La fase de núcleo congelado empieza en τ_A, como en U19b.
- **Los protocolos:**
  - poco denso: ξ_i = 0.5 en ℓ_i = 1, el mismo comienzo absoluto que en U19b;
  - denso: ξ_i = 6 en ℓ_i = ℓ₁ − 0.8.
- **El drenaje,** hasta ℓ_d = ℓ₁ − 0.5. En ℓ_G = 3.2 coincide con el 3.8 de U19b.
- **El redondeo de U19b:** allí ℓ₁ = 4.2986, así que ℓ₁ − 0.8 y ℓ₁ − 0.5 dan 3.499 y 3.799, frente a los 3.5 y 3.8 que se usaron. La diferencia, 0.001, está por debajo del paso de la malla temporal y se ignora.
- **El resto, igual que en U19b:** 48 instantáneas; el espectro desde ℓ₁ − 0.5; el contador desde ℓ₁ − 0.2; las semillas 91 y 92; dτ = 0.25.
- **Lo único que cambia con la masa** es ℓ_G, que es la jerarquía m_r/m_G y fija ℓ₁. En el universo, ℓ ≈ 90 corresponde a esa misma jerarquía, así que la variable de la extrapolación es la buena.

**Las masas:**

| ℓ_G | ℓ₁ | m_r/m_G | N | Dónde | Redes |
|---|---|---|---|---|---|
| 2.507 | 3.605 | 12.3 | 256 | nube | 4 nuevas |
| 2.912 | 4.011 | 18.4 | 384 | nube | 4 nuevas |
| 3.2 | 4.299 | 24.5 | 512 | — | **las 4 primarias de U19b**, reutilizadas: ya tienen este protocolo |
| 3.423 | 4.522 | 30.7 | 640 | WSL, kernel de memoria reducida | 4 nuevas |

**La masa más alta practicable es 3.423 (640³).**
- 768³ (ℓ_G = 3.605) necesitaría ≈ 95 GB, o ≈ 77 GB con la variante de memoria reducida.
- 640³ necesita ≈ 45 GB con ella, y cabe en WSL.
- Si hubiera una máquina con ≥ 96 GB, 768³ se añadiría por enmienda **antes** de verla.

Las redes de ℓ_G = 2.4 de U19b se informan, pero no entran en el ajuste: su protocolo es distinto.

## La variante de memoria reducida del kernel (`kern.c`, `evolve_lowmem`)

- **Qué hace:** el mismo leapfrog KDK, sin el vector de fuerza intermedio. La fuerza se recalcula en cada medio paso, así que hay el doble de cálculo de fuerza y un tercio menos de memoria.
- **Se exige que sea idéntica bit a bit a `evolve`:** el estado (f, p) tras 40 pasos en 32³, para RP³ y U(1), con y sin inclinación, núcleo congelado y drenaje. Es la prueba **V0′**, y entra en `make test`.
- **Se usa solo en la masa de 640³.**

## Validación (en la nube, antes de cualquier red de producción de cada masa nueva)

V1–V5 de U19b, con sus mismos umbrales, para cada ℓ_G nuevo y con su geometría (τ_A = N/2):
- **V1:** la masa, ≤ 1 %.
- **V2:** el contador frente a la EDO homogénea, ≤ 5 %, con los arranques temprano y denso.
- **V3:** la fase tardía, con el número conservado a ≤ 1 % y la frecuencia física a ≤ 1 %.
- **V4:** el contador espectral, ≤ 2 %.
- **V5:** κ, que se informa y se resta.

Más **V0′**. Si falla alguna, esa masa es **INVALIDA** y no se produce.

## Condiciones de validez de cada red

Las de U19b, sin cambios:
- ξ(τ_B) < 1 % de ξ₁;
- la fracción de pared en τ_B, < 5 % de su máximo;
- los sitios lejos del vacío, < 1 %;
- la meseta de N_lin, < 3 %.

Si alguna red de una masa no las cumple, esa masa no entra en el ajuste, y se informa.

## Lo que se mide (`network/u19c_analyze.py`, con las funciones de `u19b_analyze.py`)

- **Por red:** X_lat, Y_lat, ξ₁ℓ₁, r_W, r_N y c_p, con las definiciones de U19b (enmiendas 1–3).
- **Por masa,** con 2 protocolos × 2 semillas: la descomposición afín de la enmienda 5 de U19b, X = X₀ + y_s·Z, con Z = (8π/27)·ξ₁ℓ₁, X₀ ± σ y y_s ± σ.

**El ajuste**, sobre las masas válidas de la tabla: ln X₀ = ln A + p·ln ℓ₁, por mínimos cuadrados ponderados con σ(X₀)/X₀. El error de p se multiplica por √(χ²/ν) si χ²/ν > 1.

**El X del universo** (B2, ℓ = 90, como en la enmienda 5 de U19b):
- X_univ = A·90^p + 162·max(y_s + 2σ, 0), con el y_s más desfavorable de las masas;
- se informan sus valores con p − σ_p, p y p + σ_p.

**Comparación de modelos, informativa:**
- X₀ constante;
- X₀ ∝ ℓ₁ (p = 1);
- X₀ lineal en ℓ₁ con ordenada;
- la potencia libre.

Para cada una se dan χ² y su X_univ.

## Veredicto

Con X_lo y X_hi, el X_univ con p − σ_p y p + σ_p, y los umbrales de U17f (194, 275 y 320):

| Veredicto | Condición | Qué significa (con U17f) |
|---|---|---|
| **X_UNIV_BAJO** | X_hi ≤ 194 | el borde llega a F₀, solo por forma en el cribado |
| **X_UNIV_MEDIO** | 194 < X_lo y X_hi ≤ 275 | el borde es 0.85·F₀ y F₀ queda fuera |
| **X_UNIV_ALTO** | X_lo > 275 | el borde es 0.7·F₀ |
| **X_UNIV_INDECISO** | el intervalo cruza 194 o 275 | se informa entre qué bordes queda |
| **INDETERMINADO** | menos de 3 masas válidas, falta la de 640³, o σ_p > 0.5 | |

**El veredicto es sobre X.** Su traducción al borde de F usa la puerta de cribado de U17f.
- Si U15d cambia el criterio de Lyα en el borde (por ejemplo, F0_PERMITIDO_LYA_COMPLETO), la traducción se rehace con U15d en el informe, sin cambiar el veredicto de U19c.
- La cota secundaria de supermasivos, que excluye F₀, no cambia.

## Lo que se espera, escrito antes de calcular

- Si la tendencia de U19b se mantiene (p ≈ 1.2 ± 0.2), con una palanca algo mayor (ℓ₁ de 3.6 a 4.5) y cuatro puntos, σ_p debería bajar a ≈ 0.15.
- X_univ ≈ 300, con X_lo ≈ 230 y X_hi ≈ 400.
- **Se espera X_UNIV_INDECISO** (entre 0.85 y 0.7·F₀) **o X_UNIV_ALTO.** Si el protocolo de U19b infló el exponente, p bajaría hacia 1 y X_univ hacia ≈ 200: X_UNIV_MEDIO o INDECISO.

## Coste

- **En la nube** (4 hilos), con el productor que salta las redes hechas:
  - 256³: ≈ 0.5 h por red;
  - 384³: ≈ 2 h por red;
  - en total, ≈ 10 h.
- **En WSL:** 640³ con `evolve_lowmem`, ≈ 15 h por red con la máquina dedicada (sin U15c ni U15d a la vez). Son ≈ 60 h para las 4.
- La GPU de 8 GB no da para 640³, y no se usa en U19c.

## Lo que U19c no resuelve

- **La extrapolación sigue siendo larga:** de ℓ₁ ≈ 4.5 a 90.
- **R_ρ y la ξ de B2 a ℓ ≈ 90,** como en U19b.
- **La parte dura del colapso** (p ≳ 10³ m_G) sigue sin caber en la red.
- **La caja menor que el horizonte en la fase tardía,** como en U19b, aunque ahora igual en todas las masas.

## Enmienda 1 · 20 de octubre de 2026, 10:55 (Madrid): la memoria de 640³ corregida, y una masa de reserva; antes de cualquier red

- **La medida.** Se repitió la medida de memoria separando las fases, con una red RP³ de 128³:
  - las condiciones iniciales tienen un pico de ≈ 0.21 GB;
  - la evolución con medidas, de ≈ 0.48 GB, de los que ≈ 0.08 GB son el vector de fuerza.
- **Escalado a 640³ con `evolve_lowmem`:** ≈ 50 GB, no los ≈ 45 GB del texto.
- **La regla, fijada ahora e independiente de los datos:** antes de lanzar la masa alta en WSL, se mira la memoria disponible para WSL (`free -g`).
  - Con ≥ 56 GB, se produce **640³ (ℓ_G = 3.423)**.
  - Si no, **576³ (ℓ_G = 3.2 + ln(576/512) = 3.318, ℓ₁ = 4.417)**, con ≈ 36 GB.

  En los dos casos el protocolo es el mismo, relativo a ℓ₁, y se dice cuál se usó.
- **La validación** (V1–V5) se hace en la nube para las dos, 3.423 y 3.318, antes de saber cuál se produce.
- **V0′ pasó:** `evolve_lowmem` es idéntico bit a bit a `evolve`, en una llamada y en dos, para RP³ y U(1), sin opciones, con inclinación y núcleo congelado, y con drenaje. Es `test_evolve_lowmem_bit_identical`, y entra en `make test`.

## Nota · 20 de octubre de 2026, 11:51 (Madrid): la regla de la enmienda 1, aplicada; antes de cualquier red

- **La memoria de WSL:** el autor informa `memory=52GB` (más 8 GB de swap) en `.wslconfig`, y `free -g` da **49 GB disponibles**.
- **La regla de la enmienda 1:** con menos de 56 GB, **la masa alta es 576³, con ℓ_G = 3.318 y ℓ₁ = 4.417.**
- **Por qué no se sube el límite de WSL para producir 640³:** el anfitrión tiene 64 GB, y Windows se quedaría con ≈ 6 GB durante unas 60 h, con riesgo de quedarse sin memoria a mitad de una red. La palanca solo ganaría 0.1 en ℓ₁.
- **La decisión no depende de ningún dato de U19c.**
- **Las cuatro masas del ajuste:** 2.507 (256³), 2.912 (384³), 3.2 (512³, de U19b) y 3.318 (576³).
- **Coste en WSL:** ≈ 10 h por red con `evolve_lowmem` y la máquina dedicada, ≈ 40 h para las cuatro.

## Nota · 20 de octubre de 2026, 13:22 (Madrid): la validación de las masas nuevas; antes de cualquier red

La validación se hizo en la nube, con `u19c_produce.py validate 256 384 576 640`. El resultado está en `network/results/u19c/U19C_VALIDATION.json` y el registro, en `validate.out`.

| Masa (ℓ_G) | N | V1 | V2 | V3 | V4 | Resultado |
|---|---|---|---|---|---|---|
| 2.507 | 256³ | pasa | pasa | pasa | pasa | **PASA** |
| 2.912 | 384³ | pasa | pasa | pasa | pasa | **PASA** |
| 3.318 | 576³ | pasa | pasa | pasa | pasa | **PASA** |
| 3.423 | 640³ | pasa | **falla** (RP³ denso: 5.4 % > 5 %) | pasa | pasa | **INVALIDA** |

- **640³ no se iba a producir** (la regla de la enmienda 1, aplicada en la nota de las 11:51). Su fallo no cambia nada, pero se informa. Con esta validación, 640³ sería INVALIDA aunque WSL tuviera 56 GB.
- **El margen de V2 en RP³ con el arranque denso se estrecha al subir la masa,** y se declara. El contador frente a la EDO, |R/R_EDO − 1|, vale:

  | ℓ_G | 2.4 (U19b) | 2.507 | 2.912 | 3.2 (U19b) | 3.318 | 3.423 |
  |---|---|---|---|---|---|---|
  | RP³ denso | 0.026 | 0.033 | 0.036 | 0.044 | 0.049 | 0.054 |

  - En el arranque temprano, y en U(1) con los dos arranques, se queda en ≤ 2 %.
  - **576³ pasa, pero al borde del umbral.** En el informe se dirá qué cambia en X₀ si la masa alta se pondera con ese sesgo del contador (informativo, sin cambiar el veredicto).
- **Las redes de 256³ y 384³** se lanzaron en la nube a las 13:20, con `produce 256` y después `produce 384`. Su registro está en `produce.out` y `runs.log`.

## Enmienda 2 · 20 de octubre de 2026, 14:28 (Madrid): la masa alta en una A100 de Colab con el núcleo JAX, y 768³ si valida; antes de cualquier red de 576³ o 768³

**Motivo.** El autor tiene acceso a una A100 de 80 GB en Colab, con 167 GB de RAM. Producir allí la masa alta deja libre la máquina WSL para las cadenas de U15d y permite, además, una masa más: el texto original ya preveía añadir 768³ por enmienda «si hubiera una máquina con ≥ 96 GB».

**Lo que ya se ha visto, y se declara:**
- **La prueba de acuerdo en la A100**, hecha por el autor con `network/test_jax_kernel.py` en `main`: todo pasa. La diferencia con `kern.c` es **exactamente cero** en f y p (40 pasos a 32³, con inclinación, núcleo tardío y drenaje) y en ξ (redes de 48³ de U(1) y RP³).
- **La velocidad por paso** de JAX en la A100, frente a `kern.c` con los 12 vCPU de Colab:

  | N | RP³ (JAX) | RP³ (kern.c) | U(1) (JAX) |
  |---|---|---|---|
  | 256³ | 0.021 s | 0.73 s | 0.002 s |
  | 384³ | 0.072 s | 2.6 s | 0.006 s |
  | 512³ | 0.244 s | 6.0 s | 0.015 s |
  | 576³ | 0.236 s | 5.0 s | 0.020 s |

  En RP³ es entre 21 y 36 veces más rápido.
- **768³ con `evolve` de JAX no cabe en la GPU.** La compilación pedía 67.7 GB: f, p y dos copias de la fuerza.
- **Ninguna red de U19c** en 576³ ni en 768³. Las de 256³ siguen en marcha en la nube, con `kern.c`. La primera, `RP3_N256_lG2.507_sparse_s91`, que es la referencia del control (b), ya terminó. **No se ha mirado su X_lat:** `u19c_control.py` solo se probó comparándola consigo misma, y sin imprimir valores.
- **La tendencia de V2** (nota de las 13:22): en RP³ denso, la discrepancia del contador sube con la masa hasta el 5.4 % en ℓ_G = 3.423, que falló. **Se espera que 768³ (ℓ_G = 3.605) falle V2** y quede INVALIDA. El umbral del 5 % no se cambia.

**Lo que se añade:**
1. **`evolve_lowmem` en JAX** (`network/jax_kernel.py`), con el mismo esquema que el de `kern.c` (sin el vector de fuerza, con dos evaluaciones de la fuerza por paso). En CPU es **idéntico bit a bit** a `kern.c`, solo y a través del conmutador de `lattice.lib()` que usa `u19c_produce.py`; lo comprueban las pruebas nuevas de `test_jax_kernel.py`. En la A100, el estado de 768³ baja a ≈ 35 GB.
2. **El motor:** con `WHALEFALL_BACKEND=jax`, `u19c_produce.py` evoluciona en la GPU y mide en la CPU, como siempre. Cada red guarda `backend`. `U19C_OUT` permite escribir en Google Drive.
3. **768³ (ℓ_G = 3.605, ℓ₁ = 4.704)**, con el mismo protocolo relativo. Se valida (V1–V5) en la nube antes de producirla. Si pasa, se produce en la A100.

**Las condiciones, fijadas ahora:**
- **(a) La prueba de acuerdo en la misma sesión.** Antes de cada producción en la A100, `test_jax_kernel.py`, ya con `evolve_lowmem`, tiene que dar `ALL PASSED` en esa sesión, y se guarda su salida.
- **(b) Una red de control.** La primera red en la A100 es `RP3_N256_lG2.507_sparse_s91`, que ya existe con `kern.c` en la nube. Pasa si:
  - el ξ de cada instantánea coincide con el de `kern.c` dentro del 2 %. En las instantáneas en que la red de `kern.c` ya no tiene cuerdas (ξ = 0), el ξ de la A100 tiene que ser ≤ 2 % del ξ máximo;
  - el X_lat de las dos (`u19b_analyze.per_run`) coincide dentro del 2 %.

  La comparación la hace `network/u19c_control.py`, que escribe `U19C_CONTROL.json`. Si no pasa, **nada de la A100 entra en U19c**, y la masa alta vuelve a ser 576³ en WSL, con `kern.c`. La red de control no sustituye a la de la nube: solo se usa para esta comparación.
- **(c) La masa alta:** es 768³ si valida y sus cuatro redes terminan. Si no, es 576³, producida en la A100 si (a) y (b) pasan, o en WSL si no. En el ajuste entran todas las masas válidas y completas: 4 o 5.
- **(d) Lo que no cambia:** las definiciones, los umbrales de validez por red, el ajuste, los umbrales del veredicto y el término de las cuerdas.

**Coste:** unas 58 h de A100 disponibles. La evolución de una red de 576³ cuesta ≈ 3 300 pasos × 0.24 s ≈ 13 min; la de una de 768³, ≈ 4 500 pasos × (≈ 1 s con `evolve_lowmem`) ≈ 1.3 h. El resto son las medidas en la CPU: `network/bench_u19c_snapshot.py N` da su coste por instantánea antes de lanzar nada.

## Nota · 20 de octubre de 2026, 15:54 (Madrid): 768³ no valida; la masa alta es 576³, en la A100

- **La validación de 768³ (ℓ_G = 3.605):** V1, V3 y V4 pasan, y **V2 falla** en RP³ con el arranque denso: |R/R_EDO − 1| = 0.065 > 0.05. Es lo que se esperaba en la enmienda 2 (la tendencia de la nota de las 13:22). **768³ es INVALIDA y no se produce.** El registro está en `validate_768.out` y el resultado, en `U19C_VALIDATION.json`.
- **Por la regla (c),** la masa alta es **576³**, producida en la A100 si pasan (a) y (b), o en WSL si no.
- **Un arreglo de código, antes de cualquier red con `evolve_lowmem` de JAX** (commit `e9e368b`). En la A100, el compilador seguía guardando copias completas de f o de p, y 768³ no cabía. Por eso la deriva y el drenaje de `evolve_lowmem` pasan a ser bucles por capas que escriben en el sitio, como la fuerza.
  - Sigue siendo idéntico bit a bit a `kern.c` en CPU (`test_jax_kernel.py`), y `make test` pasa.
  - Con 768³ fuera, solo importa para 576³, que ya cabía.

## Nota · 20 de octubre de 2026, 16:27 (Madrid): la red de control en la A100 pasa

- **La GPU es una A100-SXM4 de 40 GB**, no de 80, con 83 GB de RAM. 576³ cabe con `evolve_lowmem` (≈ 15 GB de estado). 768³ no cabría, pero ya era INVALIDA.
- **El control (b)**, ejecutado por el autor en Colab con el código `c6ac303` y JAX 0.11.1, es **PASA**:
  - el ξ de las 48 instantáneas es idéntico al de `kern.c` (diferencia máxima, 0);
  - X_lat difiere en 3.6×10⁻⁹ relativo.

  El resultado está en `network/results/u19c/U19C_CONTROL.json`, copiado de la salida de `u19c_control.py`.
- **Se declara:** esa salida imprime el X_lat de la red de referencia `RP3_N256_lG2.507_sparse_s91`, que vale 6.43. Es el único dato de una red de U19c que se ha visto antes del análisis, y no cambia nada de lo fijado.
- **La red de control tardó 581 s en la A100,** frente a 2 128 s con `kern.c` y 4 hilos en la nube. Escalando por volumen, cada red de 576³ costará ≈ 2 h.
- **Antes de producir 576³,** la condición (a) pide pasar la prueba de acuerdo en la misma sesión. La primera línea de la celda de producción lo hace.

## Nota · 21 de octubre de 2026, 09:58 (Madrid): el resultado

- **Las cuatro redes de 576³** se calcularon en la A100 (JAX 0.11.1, `evolve_lowmem`), en 2.2–2.5 h por red.
  - La prueba de acuerdo pasó en la misma sesión antes de producir (`network/results/u19c/a100/a100_acuerdo.log`).
  - `u19c_control.py` se volvió a ejecutar sobre la red de control subida: PASA.
- **`u19c_analyze.py`: INDETERMINADO (σ_p > 0.5).** p = 1.49 ± 0.76, con χ²/ν = 7.2. El informe es `docs/informes/INFORME_U19c.md`.

