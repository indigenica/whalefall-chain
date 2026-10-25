# U5 · ¿Sobrevive el sector oscuro si las redes siguen creciendo?

*9 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U5.md`, con la enmienda 1 (R_Γ del atractor y pendiente espectral U5f).
- *Resultados:* `lya/results/u5/U5_RESULT.json` (el log, en `u5_class.log`) y las redes U5d en `network/results/u4/runs/U1_N768_s6[3-4]_xi0.5.json`.
- *Configuración:* CLASS v3.4.1 con el parche de cuadratura. Mismas puertas y umbrales que U4a: estricto δA = 0.234 y forma estricta como puerta primaria. R_Γ = 1.40 ± 0.13 (U4b fase 2).

## 0. En una frase

**El sector oscuro de la ballena sobrevive solo si, en el ℓ físico (≈ 50–85), la densidad de cuerdas deja de crecer y el espectro de emisión no es más empinado que el log-plano (q ≤ 1.02).** Los dos veredictos preregistrados son:
- **U5a: DEPENDE_DEL_ESCALADO;**
- **U5f: DEPENDE_DE_q.**

Nuestras propias redes, extrapoladas tal cual, caen del lado excluido en los dos casos: ξ crece (4.4σ) y la pendiente espectral de RP³ es q ≈ 1.1–1.2 en ℓ ≈ 5.4.

## 1. Referencia: B2 log-plano con el R_Γ del atractor

12/12 filas permitidas:
- δA = 0.104–0.167, por debajo de los 0.133–0.200 de U4a con el R_Γ de U2;
- f_em = 0.13–0.51, m_G = 1.0–5.8×10⁻¹⁵ eV, ΔN_eff(CMB) ≤ 0.005.

Con la malla ampliada (U5b), el mejor ajuste de forma es **≈ 2 % de una componente térmica de 0.2–0.33 keV**. Es decir: un pequeño componente muy caliente sobre un grueso frío. Ese límite se aplica fuera del rango de masas en que se ajustó (véase §4).

## 2. U5a · ξ creciendo con la pendiente medida (rama B4): 0/72 permitidas

| Rama | dξ/dℓ | f_em | δA |
|---|---|---|---|
| B4a (nuestra normalización, ξ(5.65) = 0.616) | 0.062 / 0.105 / 0.181 | 0.40–0.93 | 0.33–0.64 |
| B4b (normalizada a ξ(5.65) = 1.2) | 0.062 / 0.105 / 0.181 | 0.57–0.96 | 0.48–0.78 |

- Ninguna fila pasa, ni siquiera con la pendiente mínima y R a −1σ.
- **La pendiente crítica está por debajo de 0.062.** Con s = 0.062, la fila central (x₀ = 10, K_sw = 0) tiene δA = 0.40. Con ξ constante en nuestra normalización (s = 0), la emisión es la mitad que en B2, que pasa con 0.157. La extrapolación lineal entre 0.062 y 0.105 no es fiable porque δA(s) no es lineal.
- **U5d:** cuatro semillas U(1) a 768³ con ξ_i = 0.5 dan pendientes de 0.062, 0.181, 0.104 y 0.184. La media es **0.133 ± 0.030**, 4.4σ por encima de cero.
- En valor relativo, 0.133/0.63 ≈ 0.21 por unidad de ℓ, compatible con ξ ∝ ℓ (1/ℓ = 0.18). Es crecimiento logarítmico, con unas 0.11 ℓ de normalización.
- **Veredicto: DEPENDE_DEL_ESCALADO.** El sector oscuro vive solo si ξ deja de crecer entre ℓ ≈ 6 y 70, como sostienen Hindmarsh et al. (aproximación transitoria a un ξ* constante). Con la lectura de GHV y Buschmann et al. (ξ ∝ ℓ), queda excluido aunque nuestra normalización sea la mitad de la suya.

## 3. U5f · pendiente espectral q: solo pasa q ≤ 1

| B2 | q = 0.9 | q = 1 (ref.) | q = 1.1 | q = 1.2 |
|---|---|---|---|---|
| f_em | ≤ 0.003 | 0.13–0.51 | 0.52–0.89 | 0.67–0.94 |
| δA | ≤ 0.002 | 0.10–0.17 | 0.41–0.49 | 0.53–0.60 |
| permitidas | 12/12 | 12/12 | 0/12 | 0/12 |

- **q_crit(B2) = 1.025** (fila central, x₀ = 10, K_sw = 0). Con B4a pasan solo las filas con q = 0.9.
- **Por qué es tan sensible:** el peso por intervalo logarítmico es ∝ x^(1−q), y en ℓ ≈ 70 la emisión cubre unas 60 e-foldings de x. Con q = 1.1 la energía se va al infrarrojo: partículas lentas y numerosas, que hacen la componente tibia dominante. Con q = 0.9 se va al ultravioleta: la emisión apenas deja materia oscura, que pasa a ser casi toda de desalineamiento, fría.
- **Lo medido para RP³** (U4b fase 2): q = 1.10 y 1.22 en ℓ ≈ 5.4, y 0.75–0.97 en ℓ ≈ 4.5–5.0. La prueba de `u6_analyze.py` con esas dos semillas da dq/dℓ ≈ 0.7 ± 0.2. Con espectros cada tres instantáneas es una estimación pobre; U6 la mejora.
- **Veredicto: DEPENDE_DE_q.**
- La salida con q < 1 existe, pero no es gratis: con ella la emisión apenas cuenta y la predicción distintiva de la ballena (una fracción tibia medible en Lyα) desaparece. La materia oscura sería fría, de desalineamiento.

## 4. U5b · malla de forma ampliada (informativa)

- Ninguna de las 62 filas de U4a cambia de veredicto con la malla ampliada (m ≥ 0.2 keV, f ≥ 0.01).
- Las filas permitidas se ajustan ahora a 1–2.5 % de 0.2–0.33 keV.
- Siguen en el borde 4 filas: dos de B1-IR, una de B2 y una de B3.
- La cota f ≤ 0.14 m^1.1 se extrapola muy por debajo de 1 keV. Para una componente así de caliente, las cotas pertinentes son las de reliquias térmicas ligeras (CMB + estructura a gran escala) y habría que comprobarlas por separado.

## 5. U5c · isocurvatura de ruido blanco: pasa

- k_* = a H al final de la emisión es ≥ 1.0×10⁴ Mpc⁻¹.
- Con Δ²_iso = (k/k_*)³ (cota de orden uno, f_DM = 1) y el crecimiento desde la igualdad, el cociente con la potencia adiabática en z = 3 y k = 20 h/Mpc es como mucho **2.5×10⁻³**, por debajo del criterio del 1 %.

## 6. Qué queda

1. **No hay red que llegue a ℓ ≈ 70.** El destino del sector oscuro depende de la disputa de la literatura (GHV y Buschmann et al. frente a Hindmarsh et al.) sobre ξ y sobre q. Es una dependencia explícita y falsable: si se establece que ξ ∝ ℓ, o que q > 1, en el ℓ físico, la ballena 5 queda falsada en su forma actual.
2. **U6 (preregistrado):** más semillas RP³ con espectros en cada instantánea, para medir R_Γ, R_ξ y q(ℓ) con precisión. Es el resultado propio y publicable por sí solo.
3. **Cotas de reliquias calientes** para el ≈ 2 % a ~0.25 keV de §1.
4. **Auditoría independiente de U4 y U5**, como en U1–U3.
