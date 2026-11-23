# ESPEC_U68 · El colapso cuerdas–paredes de RP³ con una ruptura general, −κ Tr(E·Q), y un triplete no degenerado

*Escrito el 2026-11-21 a las 21:43 (Madrid), antes de modificar el núcleo y de calcular nada.*

**Origen.** `docs/informes/NOTA_UV_MG.md` (PR #102), §8 y su corrección de las 20:05. Las rutas ultravioleta con protección frente a la gravedad cuántica no dan de forma genérica la ruptura −εQ₁₁ con O(3) exacta, sino una ruptura general. Su caso cuadrático es −κ Tr(E·Q) = −κ nᵀEn, con E = diag(e₁, e₂, e₃, e₄) y e₁ el mayor de los cuatro. Según `NOTA_ESPURION_E.md`, las tres masas son m_c² = κ(e₁ − e_c).

**Una precisión sobre la nota.** Dentro de la familia cuadrática, las masas solo son iguales si e₂ = e₃ = e₄, es decir, con O(3). El caso «degenerado en masa con una simetría residual discreta» de la nota necesita potenciales de grado alto, con varios vacíos (§3 de la nota), y no se prueba aquí.

## 1. La pregunta

¿El colapso con un triplete no degenerado deja la misma materia oscura blanda que el de O(3)?

**Dos hipótesis, declaradas antes:**
- **(IND) especies independientes.** Cada Goldstone se comporta como en O(3), con su propia masa. Entonces la velocidad media de los cuantos blandos, en unidades de su propia masa, no cambia.
- **(CANAL) canalización hacia el modo ligero.** El giro de π de una cuerda de RP³ puede estar en cualquier plano que contenga el vacío. Con masas distintas, las cuerdas pueden girar ese plano hacia la dirección ligera, de pared más barata, y el colapso iría sobre todo a ese modo.

Lo que importa para la banda es la velocidad de los cuantos blandos (Lyα) y la abundancia.

## 2. El modelo y el núcleo

**El núcleo.** `network/kern.c` gana `set_tilt_weights(w₀, w₁, w₂, w₃)`. La inclinación pasa de −ε_G q₀ a −ε_G Σ_a w_a q_a (a = 0–3, la diagonal de Q).
- Por defecto w = (1, 0, 0, 0), y los pesos nulos se saltan. Así, la evolución es idéntica bit a bit a la actual.
- Lo mismo vale para la corrección de U19b con a_p ≠ a.
- El modelo U(1) no se toca.

**Los pesos.** Se elige w = E − (tr E/4)·𝟙 + 𝟙/4, que tiene la misma parte de traza (Σw = 1) que el caso actual y no desplaza el modo de traza. Con ε_G = m_G², las masas son m_c² = m_G²(w₀ − w_c).

**Los protocolos (3 semillas: 91, 92 y 93).** Son los de U19b a N = 192, «sparse», con ℓ_G = 2.4: los mismos τ_A = N/2 y τ_B = 2τ_A, el núcleo congelado tardío, el drenaje y las instantáneas. Es la geometría del piloto de U19b a 192³, que fue válido (4/4).

| protocolo | E | w | m_c²/m_G² (c = 2, 3, 4) | papel |
|---|---|---|---|---|
| **O3** | diag(1, 0, 0, 0) | (1, 0, 0, 0) | (1, 1, 1) | referencia (el de U19b) |
| **GEN** | diag(4, 3, 2, 0) | (2, 1, 0, −2) | (1, 2, 4) | decide |
| **O2** | diag(4, 3, 3, 0) | (1.75, 0.75, 0.75, −2.25) | (1, 1, 4) | informativo |

- Las condiciones iniciales dependen solo de la semilla, así que cada semilla da tres redes **emparejadas**.
- El modo más ligero tiene siempre la masa m_G: la época del colapso, H ≈ m_G/3, es la misma.

**Las medidas.**
- Las de U19b en cada instantánea: ξ, la fracción de paredes, la fracción lejos del vacío y la energía.
- **Al final (τ_B):** el espectro lineal de cada modo c, contado con su propia masa efectiva. Es el contador de U19b, con la masa de cada modo en lugar de una sola.
- La masa efectiva de cada modo sale de la hessiana del potencial en el vacío del núcleo congelado tardío, como en la enmienda 1 de U19b.

## 3. Las puertas

- **G0, la identidad.** Con los pesos por defecto, el núcleo nuevo da los mismos campos, bit a bit, que el de `main`: una red RP³ de 64³ con inclinación, 40 pasos, con y sin núcleo congelado tardío. Y `make test` pasa.
- **G1, las masas.** Para cada modo y protocolo, la oscilación homogénea pequeña (como V1 de U19b) da m_c² = m_G²(w₀ − w_c) con un error < 1 %. La masa efectiva de la hessiana coincide con la oscilación de WKB a < 0.3 % en τ ≈ 2τ_A (como V3 de U19b).
- **G2, el contador por modo.** Con ondas estacionarias de k conocido en un campo sin defectos (como V4 de U19b), la ⟨k⟩ ponderada por número de cada modo coincide con la esperada a ≤ 2 %, con su propia masa.
- **G3, la validez** (las cuatro condiciones de U19b) en las 9 redes:
  - ξ final < 1 % de ξ₁;
  - paredes < 5 % de su máximo;
  - fracción lejos del vacío < 1 %;
  - número total plano a < 3 % en la última media unidad de ℓ.

## 4. El veredicto

**La cantidad.** En cada red, la velocidad media de los cuantos blandos, en unidades de a₁:

    v̄ = Σ_c Σ_{k < 10 m_c a₁} (k / (m_c a₁)) n_c(k)  /  Σ_c Σ_{k < 10 m_c a₁} n_c(k).

- Para cada semilla, D_s = v̄_GEN,s / v̄_O3,s − 1.
- D es la media de las tres semillas, y σ = desviación típica/√3 (emparejada).

**La lógica.**
- **IGUAL.** Pasan G0–G3 y |D| ≤ 0.15.
- **DISTINTO.** Pasan G0–G3, |D| > 0.30 y |D| > 3σ.
- **INTERMEDIO.** Pasan G0–G3 y no se cumple ninguno de los dos.
- **INCONCLUSO.** Falla alguna puerta.

**Lo esperado:** IGUAL, si vale (IND). Con (CANAL), espero más energía en el modo ligero y un D de signo incierto.

## 5. Lecturas informativas

- **f_L,** la fracción de cuantos blandos en el modo ligero (c = 2). (IND) predice, si el número de cada especie escala como m^{−1/2} (desalineamiento y colapso en radiación):
  - f_L = 1/(1 + 2^{−1/4} + 4^{−1/4}) = 0.392 en GEN;
  - f_L = 1/(2 + 4^{−1/4}) en O2 (contando los dos modos ligeros).

  (CANAL) predice más.
- **D_ρ,** la densidad de masa Σ_c m_c N_c frente a la O3 de la misma semilla. (IND) predice 3^{−1}·Σ (m_c/m_G)^{1/2} = 1.20 en GEN.

  Si la abundancia cambia, a F fijo la m_G de la banda se mueve, porque ρ ∝ m^{1/2}F². Se traduce como orden de magnitud, sin CLASS.
- **O2:** los mismos números.
- **Por modo:** c_p, r_W y el espectro final.
- **La orientación de las cuerdas:** en las cuerdas que quedan a ℓ ≈ ℓ₁, la fracción del giro en cada plano (1, c), a partir del n local en un anillo alrededor de las plaquetas atravesadas. Si es fácil de medir.

## 6. Lo que no prueba

- La red es de 192³, con ℓ₁ ≈ 3.5. Las cifras absolutas (X, Y) no se extrapolan aquí. Lo que se mide es la comparación emparejada, que cancela gran parte de los sistemáticos del retículo.
- El cociente de masas está limitado a 2 (m₄/m₂) por la ventana del retículo.
- No pasa nada por CLASS. Si sale DISTINTO, el siguiente paso sería U17e con el espectro medido.

## 7. Salidas

- `network/u68_produce.py`: comprobaciones (`checks`) y producción (`produce`, que salta las redes ya hechas). Escribe en `network/results/u68/`, con `produce.out`.
- `network/u68_analyze.py` escribe `network/results/u68/U68_RESULT.json`.
- El informe va en `docs/informes/INFORME_U68.md`.
- **Coste:** ≈ 14 min por red a 4 hilos, ≈ 2–2.5 h en total, en la nube.
