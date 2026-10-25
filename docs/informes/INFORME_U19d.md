# U19d · Más semillas por masa: X_UNIV_INDECISO entre 244 y 1 367. F₀ queda por debajo del intervalo

*22 de octubre de 2026, 20:45 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U19d.md`.
- *Código:* `network/u19d_produce.py` y `u19d_analyze.py` (el análisis de U19c, sin cambios).
- *Resultados:* `network/results/u19d/`:
  - `U19D_RESULT.json`;
  - las 32 redes nuevas, en `runs/`;
  - los registros de Colab, en `a100/`.
- *Ejecución:*
  - **256³ y 384³:** nube, `kern.c` con 4 hilos, del 21 al 22 de octubre. Hubo cuatro reinicios del contenedor, anotados en `produce.out`; las redes cortadas se repitieron enteras.
  - **512³ y 576³:** A100-SXM4 de 80 GB en Colab, JAX 0.11.1, del 21 al 22 de octubre, en dos sesiones; las ejecutó el autor. La prueba de acuerdo dio `ALL PASSED` en las dos.

## 0. En una frase

**Veredicto preregistrado: X_UNIV_INDECISO (entre 244 y 1 367).**

**Con 12 redes por masa en lugar de 4, la dispersión entre masas de U19c resultó ser ruido:** X₀ crece ahora de forma ordenada con ℓ₁ (6.5, 7.5, 8.2, 8.6). La potencia da:
- p = 1.35 ± 0.31, con χ²/ν = 0.007;
- X del universo ≈ 560, con un intervalo de 244 a 1 367.

**Cómo leerlo:**
- **El extremo bajo (244) ya está por encima de 194,** el umbral de U17f por debajo del cual vuelve F₀. Todos los modelos que ajustan los datos dan X ≥ 214. **F₀ queda por debajo del intervalo.**
- **El extremo alto (1 367)** cruza 275 y sale muy por encima de 320, el X más alto que llegó a probar U17f.
- **El borde superior de la banda** queda entre 0.85·F₀ (si X ≲ 275) y 0.7·F₀, o por debajo de 0.7·F₀ si X ≫ 320. Eso último no se ha calculado.

## 1. Qué se produjo

| ℓ_G | N | ℓ₁ | Redes (U19c/U19b + U19d) | Motor de las nuevas | Tiempo por red |
|---|---|---|---|---|---|
| 2.507 | 256³ | 3.606 | 4 + 8 | `kern.c`, nube | 0.4–0.6 h |
| 2.912 | 384³ | 4.011 | 4 + 8 | `kern.c`, nube | 1.9–2.6 h |
| 3.2 | 512³ | 4.299 | 4 + 8 | JAX `evolve`, A100 (protocolo de U19b) | 1.3–1.6 h |
| 3.318 | 576³ | 4.417 | 4 + 8 | JAX `evolve_lowmem`, A100 | 1.9–2.4 h |

- **Las 48 redes cumplen la validez por red de U19b.**
- **Las cuatro semillas nuevas están completas en las cuatro masas,** así que la regla de las semillas completas no quita ninguna.
- **Tiempo de A100:** ≈ 29 h.

## 2. Por masa

| ℓ₁ | X₀ ± σ (U19d, 12 redes) | X₀ en U19c (4 redes) | y_s ± σ | X_lat medio ± desv. |
|---|---|---|---|---|
| 3.606 | **6.49 ± 0.32** | 6.13 ± 0.24 | +0.169 ± 0.066 | 7.21 ± 0.65 |
| 4.011 | **7.50 ± 0.43** | 6.03 ± 0.26 | +0.063 ± 0.072 | 7.84 ± 0.70 |
| 4.299 | **8.20 ± 0.35** | 8.20 ± 0.30 | −0.076 ± 0.049 | 7.72 ± 0.59 |
| 4.417 | **8.57 ± 0.42** | 7.68 ± 0.62 | −0.115 ± 0.054 | 7.78 ± 0.78 |

**Cambios respecto a U19c:**
- **El X₀ de 384³ sube de 6.0 a 7.5,** y el de 576³, de 7.7 a 8.6. El escalón de U19c entre ℓ₁ = 4.0 y 4.3 desaparece.
- **Ahora X₀ crece de forma monótona con ℓ₁.** Es la lectura (i) del ESPEC: la dispersión era ruido entre redes, y los errores de 2 grados de libertad la subestimaban.
- **y_s sigue cambiando de signo,** de positivo en las masas pequeñas a negativo en las grandes, pero ahora de forma gradual.

## 3. El ajuste preregistrado

| Magnitud | Valor |
|---|---|
| p | **1.35** |
| σ_p | **0.31** (sin escalar, porque χ²/ν < 1) |
| χ² / ν | 0.014 / 2 |
| término de las cuerdas, 162·max(y_s + 2σ, 0) | 48.7 (de 256³) |
| X_univ con p | 556 |
| X_lo (p − σ_p) | **244** |
| X_hi (p + σ_p) | **1 367** |

- σ_p ≤ 0.5, así que el resultado ya no es INDETERMINADO.
- X_lo > 194 y X_hi > 275, de modo que el intervalo cruza 275: **X_UNIV_INDECISO.**

**Un χ² tan bajo es poco probable:** P(χ² ≤ 0.014 | 2 g.l.) = 0.7 %. O el ajuste afín sobrestima un poco los errores de X₀, o es casualidad. En los dos casos σ_p es conservador: no se escala hacia abajo.

**Los modelos, informativos:**

| Modelo | χ² / ν | X del universo |
|---|---|---|
| X₀ constante | 20.4 / 3 | 56 |
| X₀ ∝ ℓ₁ | 1.35 / 3 | **218** (214–222 con el error del coeficiente) |
| X₀ lineal en ℓ₁ | 0.018 / 2 | **273 ± 48** |
| Potencia libre | 0.014 / 2 | **556** (244–1 367) |

**Qué dicen los modelos:**
- **X₀ constante queda descartado:** Δχ² ≈ 19 frente a ∝ ℓ₁.
- **Los otros tres ajustan bien** y dan X del universo de 218, 273 y 556. Los datos no los distinguen, porque la extrapolación multiplica ℓ por 22, de ≈ 4 a 90.
- **Todos los que ajustan dan X > 194.**

## 4. Estabilidad (sin veredicto)

- **Solo con las semillas nuevas (93–96):** X_UNIV_INDECISO, entre 178 y 1 483, con p = 1.27 ± 0.40. Es la misma categoría; la cota baja cae por debajo de 194 porque hay menos redes.
- **Quitando cada masa por turno,** p vale 1.36, 1.36, 1.37 y 1.33. Sin la de 256³, σ_p sube a 0.77, porque se pierde la mitad de la palanca.
- **El X_lat medio por masa** (informativo; la lectura que se eligió después de ver U19c) da p = 0.36 ± 0.17, con χ²/ν = 1.0. Extrapolado sin más a ℓ = 90, X ≈ 23.
  - Esta lectura mezcla la parte de las cuerdas, y y_s cambia con ℓ₁. Por eso no es la magnitud del veredicto.
  - La diferencia entre p ≈ 0.4 (X_lat) y p ≈ 1.35 (X₀) es justo la parte de las cuerdas: **la extrapolación depende de separarla bien.**

## 5. Lo que se esperaba

La simulación del ESPEC predecía, con p ≈ 1.5:
- ALTO, el 70 % de las veces;
- INDECISO, el 29 %;
- INDETERMINADO, el 1 %.

Ha salido INDECISO, con p = 1.35 y σ_p = 0.31. Coincide con la mediana simulada de σ_p, que era 0.22–0.36.

## 6. Qué significa

1. **El paso de U19c a U19d funcionó:** el INDETERMINADO de U19c era ruido. Con las mismas cuatro masas y tres veces más redes, X₀(ℓ₁) es suave y creciente.
2. **F₀, y con él el benchmark y LISA, queda fuera del intervalo preregistrado** (244–1 367 > 194). Es también lo que dan los tres modelos que ajustan, con 218 como valor más bajo. El margen es pequeño en el modelo ∝ ℓ₁ (218 frente a 194), y la cota secundaria de supermasivos ya excluía F₀ por su lado.
3. **El borde superior de la banda,** con la puerta de cribado de U17f:
   - 0.85·F₀ si X ≲ 275, como en los modelos ∝ ℓ₁ y lineal;
   - 0.7·F₀ si X ≈ 320;
   - **por debajo de 0.7·F₀ si X ≫ 320,** como en la potencia libre (556) o su extremo alto. U17f no tiene filas para X > 320 (Y > 1.98), así que ese tramo **no está calculado.**
4. **La lectura más prudente de la banda** pasa a ser F ≈ 0.09–0.85·F₀, con el borde superior pendiente de X. **F₀ queda descartado por la extrapolación del colapso,** además de por la cota de supermasivos.
5. **U15d** (Lyα completo en F₀ con X ≈ 195 y 275, y en 0.85 y 0.7·F₀ con X ≈ 320) se lee ahora así: con X ≥ 244, F₀ ya no está en juego por este lado. Lo que U15d decide es si el borde de 0.85 o 0.7·F₀ resiste con Lyα completo.

## 7. Lo que U19d no resuelve

- **La forma de la extrapolación.** El exponente queda en 1.35 ± 0.31, y los modelos ∝ ℓ₁ y lineal ajustan igual de bien. De ℓ₁ ≈ 4.4 a ℓ = 90, eso separa X ≈ 220 de X ≈ 560.
- **El borde con X > 320.** Haría falta una **U17g preregistrada:** las filas de CLASS de U17f con Y hasta ≈ 4 (X ≈ 650), en 0.5–0.85·F₀. Se pueden hacer en la nube en unas horas.
- **Lo que ya quedaba igual en U19c:**
  - R_ρ y la ξ de B2 a ℓ ≈ 90;
  - la parte dura del colapso;
  - la caja menor que el horizonte;
  - el límite de V2 por encima de ℓ_G ≈ 3.35.
