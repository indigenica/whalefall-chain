# U19c · El rendimiento blando del colapso de RP³ en cuatro masas con el mismo protocolo: INDETERMINADO, porque X₀ no sigue una ley limpia en ℓ₁

*21 de octubre de 2026, 09:58 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U19c.md`, con las enmiendas 1 (memoria; 576³ como reserva) y 2 (la masa alta en una A100 de Colab con el núcleo JAX) y sus notas fechadas.
- *Código:* `network/u19c_produce.py`, `u19c_analyze.py` y `u19c_control.py`; el núcleo, `network/kern.c` (`evolve_lowmem`) y `network/jax_kernel.py` (`evolve_lowmem`).
- *Resultados:* `network/results/u19c/`:
  - `U19C_RESULT.json` (el análisis), `U19C_VALIDATION.json` y `U19C_CONTROL.json`;
  - las 12 redes nuevas, en `runs/`;
  - la red de control de la A100, en `control/`;
  - los registros de Colab, en `a100/`.
- *Ejecución:*
  - 256³ y 384³ en la nube, con `kern.c` y 4 hilos, el 20 de octubre;
  - 512³: las 4 redes primarias de U19b, en WSL;
  - 576³ en una A100-SXM4 de 40 GB de Colab, con JAX 0.11.1 y la medida en 12 vCPU, el 20 de octubre (las ejecutó el autor).

## 0. En una frase

**Veredicto preregistrado: INDETERMINADO (σ_p > 0.5).**

Con el mismo protocolo en cuatro masas (ℓ₁ = 3.6, 4.0, 4.3 y 4.4), el momento blando sin cuerdas sale X₀ = 6.1, 6.0, 8.2 y 7.7. Ese patrón no es una ley de potencias:
- el ajuste da p = 1.49 con χ²/ν = 7.2;
- con el error escalado por √(χ²/ν), σ_p = 0.76.

El X del universo queda **sin fijar**. Según el modelo, informativo, va de ≈ 55 a ≈ 770. Con X₀ ∝ ℓ₁ da ≈ 200, en la frontera de F₀, como ya decía U17f.

**La banda de F no cambia:** F ≈ 0.09–0.7·F₀, o 0.09–0.5·F₀ con la cota de supermasivos.

## 1. Qué se produjo

| ℓ_G | N | ℓ₁ | Validación (V1–V5) | Redes | Motor |
|---|---|---|---|---|---|
| 2.507 | 256³ | 3.606 | PASA | 4 | `kern.c`, nube (0.6–1.0 h por red) |
| 2.912 | 384³ | 4.011 | PASA | 4 | `kern.c`, nube (2.3–2.4 h) |
| 3.2 | 512³ | 4.299 | PASA (U19b) | 4, de U19b | `kern.c`, WSL (9.4–10.7 h) |
| 3.318 | 576³ | 4.417 | PASA | 4 | JAX `evolve_lowmem`, A100 (2.2–2.5 h) |
| 3.423 | 640³ | 4.522 | **INVALIDA** (V2) | — | — |
| 3.605 | 768³ | 4.704 | **INVALIDA** (V2) | — | — |

**V2 compara el contador con la EDO homogénea, y su margen se estrecha con la masa** (nota de las 13:22). En RP³ con el arranque denso, |R/R_EDO − 1| vale:

| ℓ_G | 2.4 | 2.507 | 2.912 | 3.2 | 3.318 | 3.423 | 3.605 |
|---|---|---|---|---|---|---|---|
| RP³ denso | 0.026 | 0.033 | 0.036 | 0.044 | 0.049 | **0.054** | **0.065** |

640³ y 768³ superan el umbral del 5 %, que no se cambió. **El protocolo deja de validar por encima de ℓ_G ≈ 3.35.** Esto limita la palanca en ℓ₁ (§5).

**La A100** (enmienda 2):
- **(a) La prueba de acuerdo** pasó en la misma sesión antes de producir (`a100/a100_acuerdo.log`). Las diferencias con `kern.c` son exactamente cero en la GPU, también con el drenaje de tamaño de producción y con `evolve_lowmem`.
- **(b) La red de control** (256³, poco densa, semilla 91) da un ξ idéntico al de `kern.c` en las 48 instantáneas, y X_lat a 4×10⁻⁹. `u19c_control.py` lo ha vuelto a comprobar con la red subida.
- **Dos arreglos del núcleo JAX, hechos antes de producir y documentados en el ESPEC:**
  - la deriva y el drenaje por capas, en el sitio;
  - los autovalores del drenaje en trozos de 2¹⁵ matrices, porque cuSolver rechazaba los ≈ 4×10⁶ de una capa.
- **La masa alta es 576³**, por la regla (c) de la enmienda 2: 768³ no validó.

Todas las redes cumplen los criterios de validez por red de U19b.

## 2. Por red

| Red | ℓ_G | ℓ₁ | Protocolo | Semilla | ξ₁ | X_lat | Y_lat | r_W | c_p |
|---|---|---|---|---|---|---|---|---|---|
| 256³ | 2.507 | 3.606 | poco denso | 91 | 0.57 | 6.43 | 3.34 | 1.60 | 3.75 |
| 256³ | 2.507 | 3.606 | poco denso | 92 | 0.63 | 6.62 | 3.11 | 1.60 | 3.37 |
| 256³ | 2.507 | 3.606 | denso | 91 | 1.81 | 7.04 | 1.16 | 1.65 | 3.60 |
| 256³ | 2.507 | 3.606 | denso | 92 | 2.01 | 7.57 | 1.12 | 1.65 | 3.56 |
| 384³ | 2.912 | 4.011 | poco denso | 91 | 0.67 | 6.72 | 2.69 | 1.67 | 4.21 |
| 384³ | 2.912 | 4.011 | poco denso | 92 | 0.64 | 6.36 | 2.65 | 1.67 | 3.33 |
| 384³ | 2.912 | 4.011 | denso | 91 | 2.16 | 7.52 | 0.93 | 1.80 | 3.90 |
| 384³ | 2.912 | 4.011 | denso | 92 | 2.20 | 7.91 | 0.96 | 1.80 | 3.59 |
| 512³ | 3.2 | 4.299 | poco denso | 91 | 0.73 | 8.20 | 2.81 | 1.82 | 3.84 |
| 512³ | 3.2 | 4.299 | poco denso | 92 | 0.78 | 7.66 | 2.45 | 1.70 | 4.47 |
| 512³ | 3.2 | 4.299 | denso | 91 | 2.48 | 7.25 | 0.73 | 1.80 | 4.67 |
| 512³ | 3.2 | 4.299 | denso | 92 | 2.36 | 7.50 | 0.79 | 1.80 | 4.30 |
| 576³ | 3.318 | 4.417 | poco denso | 91 | 0.75 | 6.93 | 2.24 | 1.79 | 4.17 |
| 576³ | 3.318 | 4.417 | poco denso | 92 | 0.80 | 8.03 | 2.43 | 1.79 | 3.81 |
| 576³ | 3.318 | 4.417 | denso | 91 | 2.55 | 6.83 | 0.65 | 1.89 | 3.96 |
| 576³ | 3.318 | 4.417 | denso | 92 | 2.50 | 7.17 | 0.70 | 1.97 | 3.84 |

**El X de la red apenas cambia** en todo el intervalo: 6.4–8.2. Entre las dos semillas de un mismo protocolo varía hasta un 15 %. Lo que sí cambia de forma ordenada con la masa es el instante del colapso, r_W, que sube de 1.60 a 1.97.

## 3. Por masa: la descomposición afín X = X₀ + y_s·Z

| ℓ₁ | X₀ ± σ | y_s ± σ | X_lat medio (± desv.) | Y (bloque de U19b) |
|---|---|---|---|---|
| 3.606 | 6.13 ± 0.24 | +0.185 ± 0.050 | 6.91 ± 0.50 | 2.18 ± 0.61 |
| 4.011 | 6.03 ± 0.26 | +0.209 ± 0.044 | 7.13 ± 0.71 | 1.81 ± 0.50 |
| 4.299 | 8.20 ± 0.30 | −0.085 ± 0.042 | 7.65 ± 0.40 | 1.70 ± 0.55 |
| 4.417 | 7.68 ± 0.62 | −0.065 ± 0.081 | 7.24 ± 0.55 | 1.51 ± 0.48 |

**La descomposición cambia de signo entre las dos masas pequeñas y las dos grandes:**
- en 256³ y 384³, las redes densas (con más cuerdas) dan más X que las poco densas, y y_s > 0;
- en 512³ y 576³ ocurre al revés, y y_s < 0.

Con solo dos grupos de Z por masa, X₀ y y_s están anticorrelados: el salto de X₀ entre ℓ₁ = 4.0 y 4.3 (de 6.0 a 8.2) es en buena parte ese cambio de signo. **Los errores de X₀** salen de un ajuste de 4 redes con 2 grados de libertad, y quedan muy por debajo de la dispersión entre masas.

## 4. El ajuste preregistrado y la comparación de modelos

**Ley de potencias** (pivote en ℓ₁ = 4.00):

| Magnitud | Valor |
|---|---|
| p | 1.49 |
| σ_p (escalado) | 0.76 |
| χ²/ν | 14.4/2 = 7.2 |
| término de las cuerdas, 162·max(y_s + 2σ, 0) | 47.9 (de 384³) |
| X_univ con p | 766 |
| X_lo (p − σ_p) | 115 |
| X_hi (p + σ_p) | 7 731 |

σ_p > 0.5, así que el veredicto es **INDETERMINADO**. Además, X_lo < 194 < 275 < X_hi: el intervalo cruza los dos umbrales de U17f.

**Los modelos, informativos:**

| Modelo | χ² / ν | X del universo |
|---|---|---|
| X₀ constante | 39.0 / 3 | 55 |
| X₀ ∝ ℓ₁ | 17.7 / 3 | 201 |
| X₀ lineal en ℓ₁ | 15.7 / 2 | 261 |
| Potencia libre | 14.4 / 2 | 766 |

- **Ningún modelo ajusta bien** (χ²/ν ≥ 5.9).
- **X₀ constante es el peor,** con Δχ² ≈ 21 frente a ∝ ℓ₁: como en U19b, los datos prefieren que X₀ crezca con ℓ₁, pero no dicen cómo.
- **∝ ℓ₁ y lineal** quedan a Δχ² = 3.3 y 1.3 de la potencia libre; ∝ ℓ₁ tiene un parámetro menos.

**Lo que se esperaba** (escrito en el ESPEC antes de calcular): σ_p ≈ 0.15 y X_UNIV_INDECISO o ALTO. Ha salido INDETERMINADO, porque la dispersión entre masas es ≈ 2.7 veces la que dan los errores por masa.

## 5. Lecturas informativas (sin veredicto)

1. **El sesgo del contador de V2,** prometido en la nota de las 13:22. Si se divide el X de cada red por su (1 + |R/R_EDO − 1|) de V2:
   - X₀ = 6.08, 6.04, 8.24 y 7.74;
   - p = 1.61 ± 0.76, con χ²/ν = 7.3.

   **No cambia nada:** el sesgo es del 1–5 % y casi igual en todas las masas.
2. **El X_lat medio por masa** no es la magnitud preregistrada, porque incluye el término de las cuerdas en ℓ₁ y no separa X₀. Aun así, varía de forma más suave:
   - una potencia da p = 0.42 ± 0.24, con χ²/ν = 1.1;
   - extrapolado sin más a ℓ = 90, daría X ≈ 27 (de 13 a 56).

   Se da para mostrar cuánto depende la extrapolación de cómo se separa la parte de las cuerdas, no como estimación.
3. **El Y de bloque** (definición de U19b) baja con la masa: 2.18, 1.81, 1.70 y 1.51. La ordenación es clara, pero cada valor tiene un error de ≈ 0.5.
4. **Por qué no basta con cuatro masas.** La palanca es ln(4.417/3.606) = 0.20. Con un error real de X₀ de ≈ 7–10 % por masa, σ_p ≳ 0.4 incluso con un buen ajuste. Ampliarla hacia arriba no es posible con este protocolo: por encima de ℓ_G ≈ 3.35 falla V2.

## 6. Qué significa

- **La pregunta que U19c debía cerrar,** cómo extrapolar el colapso a ℓ ≈ 90, **sigue abierta.** Ahora se sabe por qué: en ℓ₁ = 3.6–4.4 la dependencia de X₀ con ℓ₁ queda por debajo de la dispersión entre redes, y la separación entre X₀ y la parte de las cuerdas no es estable.
- **U17f no cambia.** El borde superior sigue en 0.7·F₀ en todo el intervalo; llega a 0.85·F₀ con X ≲ 275 y a F₀ con X ≲ 194. La lectura ∝ ℓ (X ≈ 200) cae, otra vez, en la frontera de F₀.
- **U15d** (Lyα completo en el borde, en WSL) se lee igual que antes: dirá si F₀ pasa con X ≈ 195 y con X ≈ 275.
- **La mejor estimación no cambia:** F ≈ 0.09–0.7·F₀ (0.09–0.5·F₀ con supermasivos).
- **El contrato añade la fila U19c,** y el eslabón abierto dice que U19c deja la extrapolación INDETERMINADA.

## 7. Lo que U19c no resuelve, y qué podría resolverlo

**Lo que queda igual que en U19b:**
- la extrapolación de ℓ₁ ≈ 4.4 a 90;
- R_ρ y la ξ de B2 a ℓ ≈ 90;
- la parte dura del colapso (p ≳ 10³ m_G);
- la caja menor que el horizonte en la fase tardía.

**Lo que haría falta para fijar p.** No es otra masa, que no validaría, sino **más semillas por masa**: el error de X₀ hoy está dominado por la dispersión entre redes y por la degeneración X₀–y_s. Con la A100, una red cuesta ≈ 10 min en 256³, ≈ 35 min en 384³, ≈ 1.5 h en 512³ (estimado) y ≈ 2.3 h en 576³. Ocho semillas más por protocolo y masa serían ≈ 75 h de A100, más de lo disponible (≈ 50 h con los créditos que quedan). Cuatro más por protocolo, ≈ 38 h.

Eso exigiría una **U19d preregistrada**, con:
- la regla de lectura: X₀ y y_s por ajuste conjunto, o X_lat medio con su error entre semillas;
- un umbral de σ_p fijado de antemano.

**Lo decide el autor.**
