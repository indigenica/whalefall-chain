# ESPEC_U67 · La red de RP³ nacida en la era de materia: ¿olvida su nacimiento al llegar a la radiación?

*Escrito el 2026-11-21 a las 15:52 (Madrid), antes de modificar el núcleo y de calcular nada en este repositorio.*

**Declaración previa.** El 20 de noviembre, una exploración con otro asistente de IA, fuera de este repositorio, con un núcleo C propio y sin preregistro, hizo algo parecido a 192³ («Blindaje 2A»). Encontró que la red ordenada por ξ_Q R en un fondo de materia y luego de radiación queda a −5.7 % de la red de radiación pura de arranque disperso, a ℓ = 4.3. Lo he leído. U67 no usa ese código: extiende `network/kern.c` y el productor del repositorio. Las expectativas se escriben, por tanto, conociendo ese resultado.

## 1. La pregunta

En `main-v4.tex`, el supuesto abierto «Network before radiation domination» dice que Q se ordena mientras domina el inflatón oscilante (un fondo de tipo materia), y que nuestros retículos solo cubren la era de radiación. Que la red haya olvidado ese nacimiento al llegar al colapso se supone, no se calcula.

U67 lo prueba en el retículo:
- Q parte de la fase simétrica O(4) sostenida por el acoplo ξ_Q R.
- Se ordena en un fondo de materia (w = 0) y pasa después a radiación.
- La pregunta es si, al final de la ventana, su densidad de cuerdas coincide con la de las redes de radiación pura que usa el artículo de la red.

## 2. El modelo y el fondo

**El núcleo.** `network/kern.c` se extiende con dos conmutadores. Los dos están apagados por defecto, y entonces la evolución es idéntica bit a bit a la actual.
- **`set_bg_mdrd(a_tr, τ_tr)`:** un fondo de materia hasta τ_tr, a = a_tr (τ/τ_tr)², y de radiación después, a = a_tr (2τ/τ_tr − 1). Es continuo y C¹ en τ_tr.
- **`set_xiq(ξ, X₀)`:** suma a la fuerza el término ξ R (q − X₀ e_traza), la variación de ½ ξ R Tr[(Q − X₀𝟙)²] en la base de Frobenius de 10 componentes. Usa R = 3H² en la era de materia y R = 0 en la de radiación.
  - ξ = ξ_Q = 1/4, el valor de la cadena.
  - X₀ es el punto estacionario O(4)-simétrico de U_code, Q = X₀𝟙.

**Las ecuaciones de movimiento** son las del núcleo actual en tiempo conforme (π = a²q′, π′ = a²(∇²q − a²∂V)), con el a(τ) del fondo.

**Las escalas.** Las redes tienen N = 192, dτ = 0.25 y τ_end = N/2 = 96, como en producción, con a(τ_end) = 1/m_r. Así se cumple m_r a dx ≤ 1, y el horizonte no pasa de media caja.
- **La transición:** el instante τ_tr se fija con ℓ_tr = ln(m_r/H) al cambiar de era.
- **La densidad:** ξ ≡ L/(4 (aH)² V), es decir, con t ≡ 1/(2H) en cualquier fondo. Coincide con la convención de radiación del productor.

## 3. Los protocolos (3 semillas cada uno: 101, 102, 103)

- **RD_SPARSE:** el protocolo de producción de la red, a 192³. Radiación pura, arranque sintético en el vacío de RP³ con ξ_i = 0.5 en ℓ_i = 1 (`lattice.run`), y drenaje de la velocidad masiva hasta ℓ_d = 3.8 a la tasa m_r.
- **RD_DENSE:** lo mismo con ξ_i = 2. Es informativo.
- **MD2:** fondo de materia hasta ℓ_tr = 2.0 y de radiación después (ℓ_end = 4.29).
  - Arranque en ℓ₀ = 0.5 en la fase simétrica: Q = X₀𝟙 + δ, con δ un campo gaussiano de correlación 1 celda y desviación típica 0.02 en cada una de las 10 componentes, y π = 0.
  - Ese arranque está antes de ℓ_crit, donde ξ_Q R = |m²_taquión|, que se espera ≈ 1.1.
  - Con ξ_Q R encendido, mismo drenaje hasta ℓ_d = 3.8.
- **MD3:** lo mismo con ℓ_tr = 3.0 (ℓ_end = 4.11). Más tiempo en materia y menos en radiación.

La amplitud 0.02 de la semilla es una elección del retículo: la cuántica real es mucho menor y retrasaría el ordenamiento. Se declara como tal.

**Las medidas.** Unas 30 instantáneas por red, espaciadas geométricamente en τ y más densas al final. En cada una:
- ℓ, ξ y el número de plaquetas atravesadas;
- la brecha media λ₁ − λ₂ de los autovalores de Q;
- las energías.

## 4. Las puertas

- **G0, la identidad.** Con los dos conmutadores apagados, el núcleo nuevo da los mismos campos, bit a bit, que el de `main` en una red RP³ de 64³ con 40 pasos (con y sin drenaje). Y `make test` pasa.
- **G1, el fondo y el punto simétrico.**
  - a(τ) y a′(τ) son continuos en τ_tr (error relativo < 10⁻¹²).
  - X₀ cumple dU/dx = 0 a < 10⁻¹².
  - En X₀, los 9 modos sin traza tienen m² < 0 (taquiónicos) y el modo de traza, m² > 0.
- **G2, el ordenamiento.** En las 6 redes MD, la brecha media es ≥ 0.5 en la primera instantánea con ℓ ≥ 3.5: Q está en el vacío de RP³ antes del final.
- **G3, la estadística.** En las 12 redes hay ≥ 200 plaquetas atravesadas en la última instantánea.

## 5. El veredicto

Se define D_k = ξ_MDk(ℓ_c)/ξ_RD_SPARSE(ℓ_c) − 1, con las medias sobre semillas y ξ interpolado linealmente en ℓ entre instantáneas.
- **ℓ_c = 4.10** es el mayor ℓ que alcanzan todas las redes.
- **σ_k** sale de un bootstrap (10⁴ remuestreos) de las semillas, independiente en cada protocolo.

- **OLVIDA.** Pasan G0–G3, |D_2| ≤ 0.15 y |D_3| ≤ 0.15.
- **RECUERDA.** Pasan G0–G3, y para algún k |D_k| > 0.30 y |D_k| > 3σ_k.
- **INTERMEDIO.** Pasan G0–G3 y no se cumple ninguno de los dos.
- **INCONCLUSO.** Falla alguna puerta.

**Lo esperado:** OLVIDA en MD2, por la exploración declarada (−5.7 % a ℓ = 4.3 con ℓ_tr = 1.84), y quizá INTERMEDIO por MD3, que tiene solo 1.1 unidades de ℓ en radiación.

**Lo que no prueba.**
- El artículo supone ≈ 14 e-folds en la era del inflatón, ≈ 21 unidades de ℓ, y 57–66 en radiación antes del colapso. El retículo solo cubre ℓ ≲ 4.3.
- U67 mide si el mecanismo de olvido funciona en las escalas accesibles; no lo extrapola.
- Tampoco mide el rendimiento del colapso (eso sería otra prueba).

## 6. Lecturas informativas

- D(ℓ) a lo largo de la ventana, para MD2 y MD3.
- La comparación con RD_DENSE.
- El ℓ al que se completa el ordenamiento (brecha ≥ 0.5).
- La energía masiva frente a la de Goldstone al final.

## 7. Salidas

- `network/u67_produce.py` escribe las redes en `network/results/u67/runs/` y lleva `produce.log`.
- `network/u67_analyze.py` escribe `network/results/u67/U67_RESULT.json`.
- El informe va en `docs/informes/INFORME_U67.md`.
- Corre en la nube, con 4 hilos (≈ 1–2 h). Se relanza tras un reinicio, saltando las redes terminadas.

## Enmienda del 2026-11-21 a las 18:34 (Madrid): el nombre del asistente

A petición del autor, la declaración previa ya no nombra al asistente de la exploración externa, porque no se sabe con certeza qué modelo se usó. Donde lo nombraba entre paréntesis, dice ahora «otro asistente de IA». Es el único cambio en el texto preregistrado: no toca las filas, las puertas, los umbrales ni la lógica del veredicto. El texto original sigue en el historial de git.
