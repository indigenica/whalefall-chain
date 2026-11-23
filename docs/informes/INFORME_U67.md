# INFORME_U67 · La red de RP³ y la era de materia (U67, U67b)

*2026-11-21 17:21 (Madrid). Interno.*
- **Preregistros:** `docs/preregistros/ESPEC_U67.md` (15:52) y `ESPEC_U67b.md` (16:14).
- **Código:** `network/kern.c` (conmutadores `set_bg_mdrd` y `set_xiq`), `network/u67_produce.py` y `network/u67_analyze.py`.
- **Resultados:** `network/results/u67/` y `network/results/u67b/`, con 21 redes de 192³ en la nube a 4 hilos (≈ 4 min cada una).

## 1. U67: INCONCLUSO (falla G2)

| puerta | resultado |
|---|---|
| G0, la identidad del núcleo con los conmutadores apagados | **pasa:** campos idénticos bit a bit al núcleo de `main`, con y sin drenaje; `make test` pasa |
| G1, el fondo y el punto simétrico | **pasa** (ver abajo) |
| G2, Q ordenado (brecha ≥ 0.5) en la primera instantánea con ℓ ≥ 3.5 | **falla:** brecha de 0.004–0.005 en las 6 redes MD |
| G3, ≥ 200 plaquetas al final | pasa |

**Lo que mide G1:**
- X₀ = 0.235716, con dU/dx = −4×10⁻¹⁷;
- los 9 modos sin traza tienen m² = −0.0194 (taquiónicos) y el de traza, +1.58;
- ℓ_crit = 1.11;
- el fondo es continuo y C¹ en τ_tr.

**Un fallo de código, declarado.** La primera versión de la comprobación de continuidad comparaba a(τ_tr ∓ ε), que difieren en 2ε a′ aunque a sea suave. Se corrigió antes de la producción para comparar las dos ramas en τ_tr.

**Por qué falla G2.** En el punto simétrico el taquión es débil: |m²| = 0.019, frente a m_r² = 0.24. Solo crecen longitudes de onda de ≳ 100 celdas, del orden de la caja. Con una semilla de amplitud 0.02 y correlación de 1 celda:
- MD2 se ordena a ℓ ≈ 4.2;
- MD3 no llega a ordenarse.

Mientras la brecha es ≈ 0, las plaquetas «atravesadas» son ruido de fase, no cuerdas, y el ξ de esas redes no significa nada.

## 2. U67b: RECUERDA

**Diseño.** U67b se escribió después de ver U67 y una prueba de viabilidad a 96³, ambas declaradas.
- **MDA2 y MDA3 (deciden):** la red ya en el vacío de RP³, con el arranque sintético del productor en ℓ_i = 1.3. Evoluciona en la era de materia hasta ℓ_tr = 2 o 3, con ξ_Q R encendido, y luego en la de radiación.
- **Control:** RD_SPARSE de U67.

| puerta | resultado |
|---|---|
| G0, G1 | las de U67: pasan |
| G2, brecha ≥ 0.5 en toda instantánea con ℓ ≥ 2 (MDA) | **pasa:** mínimo 0.70 |
| G3 | pasa |

**El resultado:** D = ξ_MDA/ξ_RD_SPARSE − 1, con media de 3 semillas.

| ℓ | ξ, RD_SPARSE | D, MDA2 (ℓ_tr = 2) | D, MDA3 (ℓ_tr = 3) | D, RD_DENSE (informativo) |
|---|---|---|---|---|
| 2.5 | 0.94 | −0.79 | −0.83 | +12.5 |
| 3.0 | 0.54 | −0.52 | −0.66 | +7.3 |
| 3.5 | 0.61 | −0.43 | −0.58 | +0.16 |
| 3.8 | 0.68 | −0.38 | −0.54 | +0.10 |
| 4.1 | 0.75 | **−0.35 ± 0.04** | **−0.50 ± 0.02** | +0.04 ± 0.03 |

Con |D| > 0.30 y > 3σ en los dos protocolos, el veredicto preregistrado es **RECUERDA**.

**Lo que significa.**
- **La red que pasa parte de su vida en la era de materia** llega a la radiación más diluida que una nacida en radiación. Al final de la ventana conserva un déficit del 35–50 % en la densidad de cuerdas.
- **El déficit se cierra despacio,** con una escala de ≈ 3 unidades de ℓ en MDA2 y ≈ 4 en MDA3.
- **En cambio, la red de radiación pura** olvida su densidad inicial enseguida: de +1250 % a +4 % entre ℓ = 2.5 y 4.1.
- **Una lectura plausible:** en radiación, la red se acerca al atractor desde abajo más despacio que desde arriba, como ya se vio en el artículo de la red. La red que viene de la era de materia queda muy por debajo.
- **Con la escala medida,** 57–66 unidades de ℓ de radiación bastarían de sobra para borrar el déficit (e^{−57/4} ~ 10⁻⁶). Pero eso es una extrapolación fuera del retículo, no un cálculo.

**Qué cambia en el artículo.**
- **El supuesto abierto «Network before radiation domination» sigue siendo un supuesto.** U67b no confirma el olvido en las escalas accesibles: confirma que hay memoria y que decae.
- **La cota que da el artículo sigue en pie:** la banda solo se movería si la memoria cambiara la parte independiente de la cuerda del rendimiento en un factor > 2.3.
  - La densidad de cuerdas entra solo en el término de cuerdas de X_univ, 49 de 556.
  - Un déficit del 50 % movería X_univ ≈ 4 %.
- **Sería honesto añadir al supuesto, en una versión posterior,** que el retículo ve memoria de la era de materia (−35 a −50 % en ξ a ℓ ≈ 4), que se cierra con una escala de 3–4 unidades de ℓ. Lo decide el autor.

## 3. Lecturas informativas

**MDO2: el ordenamiento desde la fase simétrica,** con una semilla de amplitud 0.3 y correlación de 6 celdas.
- Q se ordena a ℓ ≈ 3.1 en las tres semillas.
- La red nace densa (D = +1.2 a ℓ = 3) y baja deprisa: +0.48 a ℓ = 3.5 y **+0.14 ± 0.03** a ℓ = 4.1.
- Desde arriba, la red olvida rápido, igual que RD_DENSE.
- La dirección y el tamaño de la memoria dependen, por tanto, de cómo nace la red: más densa (ordenamiento desde la fase simétrica) o más diluida (red ya formada que atraviesa la era de materia).

**La energía masiva frente a la de Goldstone al final:** 0.09–0.14 en las redes MD y 0.14–0.22 en las RD.

**La exploración externa declarada** (otro asistente de IA, 192³; su protocolo de «atractor MD→RD» es el equivalente de MDA) daba −76 % → −32 % entre ℓ = 2.5 y 4.3. Coincide con MDA2.

**Su protocolo de «ordenamiento»** daba −6 % a ℓ = 4.3. Aquí, MDO2 da +14 % a ℓ = 4.1 con otra semilla, y con la semilla pequeña de U67 la red no se ordena. Su resultado dependía de la semilla elegida.
