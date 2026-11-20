# U58b · El cuscutón en la inflación, solo con el juego adoptado · Preregistro

*2026-11-15 22:21 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`. El autor adopta el juego pequeño (μ = 0.16, m = 0.02, φ∞ = 0.25) para la v4 y pide cerrar el veredicto. Escrito después de U58 y antes de U58b. Interno: no se exporta.*

## 1. Qué cambia respecto a U58, declarado

- **Solo dos juegos:** GR (el control) y el adoptado. El de C1 sale porque los datos lo excluyen (U58).
- **G2 sigue igual** (c_s² de BKG frente a c_s² = 1, ≤ 10⁻⁴), ahora solo con el juego adoptado.
- **No es una prueba ciega.** U58 ya dio, con este juego, 2.2×10⁻⁵ en G2, Δn_s = 1.7×10⁻⁶ y Δr/r = +0.26 %, y B₁ y B₂ se comprobaron después en la fuente LaTeX. U58b solo formaliza el veredicto para el juego que entra en la v4.
- **Lo que añade:** el control de la amplitud. Con el juego adoptado, la normalización de A_s pide reescalar V₀ por (1 + ε_G)³ (P_ζ ∝ V/(1 + ε_G)³), y se informa cuánto cambian n_s y r, que no deberían cambiar.

## 2. El resto, como en U58

El montaje (§3 de ESPEC_U58), las puertas G0–G3 y la cota −0.025 ≤ ε_G ≤ 0.017. Los veredictos son CUSCUTON_COMPATIBLE / CUSCUTON_EXCLUIDO / CUSCUTON_TENSO / INCONCLUSO, con las mismas reglas.

## 3. Lo que se espera

Que las puertas pasen y salga **CUSCUTON_COMPATIBLE**.

## 4. Salidas

- **Código:** `chain/u58b_cuscuton_inflacion.py`, que importa `u58_cuscuton_inflacion` sin cambiarlo.
- **Resultados:** `results/u58b_cuscuton_inflacion.json`.
- **Informe:** una sección de `INFORME_U58.md`.
