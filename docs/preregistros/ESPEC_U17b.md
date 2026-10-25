# U17b · La población del colapso cuerdas–paredes, por CLASS · Preregistro

*14 de octubre de 2026, 00:35 (Madrid). Escrito antes de calcular ninguna fila de U17b. Se ejecuta en la nube.*

## Motivo

U17 (`docs/informes/INFORME_U17.md`) mostró que, si el colapso de la red en H ≈ m_G/3 se cuenta con la energía que la red tiene en ese momento, K_eff = (8π/27)·ξℓ·η/c vale 3–770 en lugar de 0–1. En las ramas con ξ creciente, la materia oscura pasa a ser casi entera la población del colapso, con m_G ≈ 10⁻²⁰–5×10⁻¹⁸ eV. La comparación con Lyα se hizo con un indicador, la longitud de libre recorrido de un momento medio frente a la WDM de 5.3 keV. U17b la sustituye por las puertas de U12, con CLASS.

## Filas

- Ramas B2 (ξ constante), B4a (crecimiento medido) y B1 (ξ ∝ ℓ, GHV), con q = 0.87, x₀ = 10 y R_Γ = 1.56: la rama fría de U12.
- Tres variantes del colapso (c = momento medio/m_G en H = m_G/3; η = fracción de la energía de la red que acaba en esa población):
  - **(c, η) = (3, 1):** c/η = 3, el valor de arXiv:1202.5851 (ω̄/m ≈ 3.3 con η = 1);
  - **(30, 1):** c/η = 30, con los cuantos más calientes;
  - **(3, 0.1):** c/η = 30, con los cuantos fríos. Separa el efecto de K del de la temperatura.
- 3 × 3 = **9 filas**.

## Cálculo (`lya/u17b_class.py`)

1. **m_G autoconsistente**, como en U17: Ω_mis(m_G)·(1 + K_eff) + Ω_em(m_G) = 0.12, con K_eff = (8π/27)·ξ(ℓ)·ℓ·η/c en ℓ = ln(3m_r/m_G). Se usa `chain/u17_ksw.py` (`k_eff`, `eps_fn`) y `chain/s6_measured.source_eps`.
2. **Especie no fría:** los paquetes de la emisión (los de U12 con esa m_G) más los del colapso:
   - número comóvil tal que m_G·n = Ω_mis·K_eff·ρ_c (la parte K de Ω_mis(1 + K));
   - momento comóvil con distribución log-normal de anchura σ_ln = 0.5, y media c·m_G·a(T_osc). La anchura imita la forma ancha de la figura 6 de arXiv:1202.5851, y se declara como hipótesis.
   - El desalineamiento sin K queda como materia oscura fría (omega_cdm).
   - Fracción f_nc = (Ω_em + Ω_colapso)/0.12.
3. **Distribución en momento** con `lya/u4a_class.psd_from_packets` y parámetros con `packet_params`. CLASS v3.4.1 con el parche (C0b obligatoria), como en U12.
4. **Puertas:** las de U12 sin cambios (`lya/u5_class.gates`): δA estricto (el mínimo de las cinco referencias) y forma frente a la familia térmica ampliada; puerta primaria la de forma si el ajuste es representable (rms ≤ 0.1) y, si no, δA estricto.

## Veredicto

Por rama: **PASA** si las 3 variantes pasan la puerta primaria, **FALLA** si no pasa ninguna y **PARCIAL** en otro caso.

Global, en este orden:
1. **SOLO_XI_CONSTANTE_VIVE:** B2 PASA y B4a y B1 FALLAN en las dos variantes con η = 1.
2. **TODAS_VIVEN:** las tres ramas PASAN.
3. **XI_CONSTANTE_TAMBIEN_EN_PELIGRO:** B2 no PASA.
4. **DEPENDE_DE_ETA:** cualquier otro caso. Se dice qué variantes salvan cada rama.

**Se informa también**, sin veredicto: m_G, K_eff, T_osc, f_nc, δA, el ajuste de forma (m, f) y la masa de los minicúmulos de cada fila; si m_G < 2×10⁻²⁰ eV (cota de Lyα para materia oscura ultraligera, que no se aplica como puerta).

## Lo que U17b no resuelve

- **c y η** siguen sin medirse para RP³; eso es U19.
- **La forma del espectro del colapso** es una log-normal supuesta.
- **Solo x₀ = 10 y R_Γ = 1.56**, como la rama central de U12.


## Corrección de horas · 14 de octubre de 2026, 00:44 (Madrid)

La hora de la cabecera (00:35) está mal: el preregistro se fijó en el commit `c7ff41b` a las 00:32:50. La ejecución de U17b se lanzó después, hacia las 00:34. El orden no cambia.
