# U14 · ¿Pueden los minicúmulos de la ballena producir microlente hacia las Nubes de Magallanes? · Preregistro

*11 de octubre de 2026, 22:10 (Madrid). Tarea 4 de `docs/PLAN_NUBE_2026-10-11.md`. U13 queda reservado para m_G desde la acción común (tarea 3).*

**Declaración previa.** Antes de escribir este preregistro, al leer la literatura, el autor hizo una estimación mental de orden de magnitud para M = 1 M_⊙: δ_lens ~ 90 y F(δ > δ_lens) ~ 0.13. No se ha ejecutado ningún código de U14. Las fórmulas, las filas y los umbrales de abajo se fijan sin ajustarlos a esa estimación.

## Motivo

`docs/PREDICCIONES.md` lista los minicúmulos como predicción con la microlente «sin evaluar». La masa característica de la ballena (C35b, Fairbairn et al. 2017, arXiv:1707.03310, ec. 2 con n = 0, masa independiente de la temperatura) es

M₀ = 2.3×10⁻⁷ M_⊙ (m_G / 10⁻¹⁰ eV)^(−3/2).

Para la ventana de m_G que sobrevive (U7c: 0.8–4.9×10⁻¹⁵ eV; rama fría de U11: 2–10×10⁻¹⁵ eV), M₀ ≈ 0.23–10 M_⊙, de lleno en el rango de EROS y MACHO.

Los minicúmulos no son lentes puntuales. Solo producen microlente estándar si su radio es menor que el radio de Einstein.

## Fórmulas (fijadas)

- **Radio de esfera dura** (Fairbairn et al. 2017, ec. 43): r_mc(δ, M) = 4×10¹⁶ cm · [δ ((δ + 1) Ω_m)^(1/3)]⁻¹ · (M/M_⊙)^(1/3), con Ω_m = 0.32.
- **Radio de Einstein hacia la LMC** (su ec. 33), en el punto medio x = 1/2 y con D_S = 50 kpc: R_E(M) = 2 [G M x(1 − x) D_S / c²]^(1/2).
- **δ_lens(M):** la solución de r_mc(δ, M) = R_E(M). Por encima de δ_lens, el minicúmulo se trata como lente puntual; por debajo, no lente.
- **Distribución de δ** (ajuste de Fairbairn et al. 2017 a Kolb y Tkachev): F(δ > δ₀) = [1 + δ₀/1.023]^(−0.462). Se interpreta como fracción en masa.
- **Fracción de la materia oscura que lente como objeto compacto:** f_lens(M) = f_MC · F(δ > δ_lens(M)), con f_MC la fracción de la materia oscura en minicúmulos.
- **Cotas de lente puntual:**
  - **estricta:** f_max = 0.08 (EROS-2, Tisserand et al. 2007, astro-ph/0607207, para M ≈ 0.4 M_⊙), aplicada en todo el rango;
  - **amplia:** f_max = 0.20 (EROS-2 + MACHO, Blaineau et al. 2022, arXiv:2205.09450, de 10⁻⁷ a 200 M_⊙).
- **Cota sobre f_MC:** f_MC,max(M) = f_max / F(δ > δ_lens(M)).

## Filas

- M ∈ {0.23, 0.5, 1, 2, 5, 10} M_⊙;
- las dos cotas;
- **variantes de radio** (informativas): r_mc × 3 y r_mc / 3, por la diferencia de un factor ~3 entre la esfera dura y el radio de escala NFW (nota 8 de Fairbairn et al.).

`python3 lya/u14_microlens.py` → `lya/results/u14/U14_RESULT.json`.

## Veredicto (con el radio de esfera dura)

- **MICROLENTE_SIN_EFECTO:** f_MC,max ≥ 1 en todas las masas con la cota estricta. Ni con toda la materia oscura en minicúmulos se viola EROS.
- **MICROLENTE_ACOTA_f_MC:** con la cota estricta f_MC,max < 1 en alguna masa, pero con la amplia f_MC,max ≥ 1 en todas. Se informa el mínimo de f_MC,max.
- **MICROLENTE_EXIGE_f_MC_PEQUEÑA:** con la cota amplia f_MC,max < 1 en alguna masa. Se informa el mínimo.

Ningún veredicto falsa la ballena: f_MC no está calculada en el modelo. Si sale el tercero, la ballena necesita que la mayor parte de su materia oscura no quede en minicúmulos compactos, y eso pasa a ser una condición explícita.

## Lo que U14 no resuelve

- **La distribución de δ de la ballena.** El ajuste de Kolb y Tkachev es de axiones QCD con masa dependiente de la temperatura. Para n = 0 (masa constante) puede ser distinta, y la cola δ ≳ 50 está extrapolada.
- **La función de masas.** Se usa M₀ y no la distribución completa, ni las fusiones posteriores (que dan masas mayores y menos densas), ni la destrucción por encuentros estelares.
- **El perfil.** Esfera dura frente a NFW o perfil autosimilar; las variantes ×3 y /3 dan una idea.
- **Otras sondas** de minihalos de masa solar: cruces de cáusticas de estrellas muy magnificadas y los arreglos de púlsares. Se mencionan, pero no se evalúan.
