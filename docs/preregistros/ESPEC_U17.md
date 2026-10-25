# U17 · ¿Cubre K_sw ∈ {0, 1} el colapso de la red en H ≈ m_G? · Preregistro

*12 de octubre de 2026, 22:55 (Madrid). Escrito antes de calcular ninguna fila de U17. Cálculo de la cadena, sin CLASS (segundos).*

## Motivo

En la cadena, la contribución del colapso cuerdas–paredes entra como un factor constante sobre el desalineamiento: Ω_mis·(1 + K_sw), con K_sw ∈ {0, 1} en U4–U12 y {0, 1, 10} en S5 (`components/c35b_pseudo_goldstone_dm_rp3.py`). El valor «medido» (0.27, C40) salió de una rampa corta en la red, y los propios informes de cierre dicen que ese proxy no certifica la emisión cosmológica.

Pero la energía que queda en la red cuando aparecen las paredes (H ≈ m_G) **no es constante**: crece con ξ y con la tensión, ∝ ξ(ℓ)·ℓ, y en H ≈ m_G vale ℓ ≈ 85. Un árbitro de axiones preguntará por esto lo primero.

**Estimación de orden de magnitud**, hecha al escribir este preregistro (no es el cálculo):
- energía de la red en H = m_G/3: ρ_red = 4ξπf_a²ℓH²;
- desalineamiento en la misma época: ρ_mis = (3/8)m_G²F², con F = 2f_a;
- si una fracción η de ρ_red acaba en Goldstones de energía media c·m_G, el número relativo es

  **K_eff = n_colapso/n_mis = (8π/27)·η·ξℓ/c ≈ 0.93·ξℓ/(c/η).**

- Con ℓ ≈ 85: B2 (ξ = 1.2) da K_eff ≈ 94/(c/η); B1 (ξ = 0.24ℓ) da ≈ 1700/(c/η). Para que K_eff ≤ 1 en B2 haría falta c/η ≳ 100, es decir, que casi toda la energía de la red se fuera a cuantos duros.

Lo que no se sabe y se declara: c y η para una red ℤ₂ de RP³ con N_DW = 1. Para U(1) con N_DW = 1, la literatura (Hiramatsu, Kawasaki, Saikawa y Sekiguchi, arXiv:1202.5851) describe los axiones del colapso como «apenas relativistas», es decir, c de orden unidad; no se han podido leer sus números en esta sesión. Por eso c/η se recorre como un parámetro.

## Cálculo (`chain/u17_ksw.py`)

Para cada rama de U12 (B2, B4a, B4b, B1, con ξ(ℓ) de ESPEC_U12), q ∈ {0.87, 1.0}, x₀ ∈ {10, 30} y **c/η ∈ {3, 10, 30}**:
1. **m_G autoconsistente.** Se resuelve Ω_mis(m_G)·(1 + K_eff(m_G)) + Ω_em(m_G) = 0.12 con `chain/s4_dark.py` (la emisión, sin cambios) y `c35b.omega(m, F, K)`, donde K_eff se recalcula con el ℓ y el ξ de la rama en H = m_G/3 para cada m_G que pruebe el buscador de raíces.
2. **Lo que sale:** K_eff, m_G, T_osc, f_em (la fracción emitida, la de U12), la masa de los minicúmulos con la misma escala de C35b, y ΔN_eff (sin cambios: viene de la emisión).
3. **¿Son fríos los cuantos del colapso?** Longitud de libre recorrido de un Goldstone con momento c·m_G en H = m_G/3, con `s4_dark.lambda_fs`, frente a la WDM térmica de 5.3 keV (`s4_dark.lambda_wdm_reference`, 0.0718 Mpc). Es una cota del lado conservador: se toma toda la población con el momento medio.

Se dan también, como referencia, los valores con K_sw = 0 y 1 (los de U12).

## Veredicto

Con c/η = 10 como caso central y 3 y 30 como sistemática:
- **K_SW_DENTRO_DEL_RANGO:** K_eff ≤ 1 en todas las filas con c/η = 10. El corchete {0, 1} de la cadena era suficiente.
- **K_SW_FUERA_DEL_RANGO_FRIA:** K_eff > 1 en alguna fila con c/η = 10, pero, para todo c/η ∈ {3, 10, 30}, (i) los cuantos del colapso tienen λ_fs < λ_WDM(5.3 keV) y (ii) las filas de U12 con q = 0.87 siguen con f_em < 1 %. El veredicto de U12 se mantiene, pero cambian las predicciones de m_G, T_osc y minicúmulos: se dan los rangos nuevos.
- **K_SW_FUERA_DEL_RANGO_TIBIA:** K_eff > 1 y, para algún c/η ∈ {3, 10, 30}, los cuantos del colapso tienen λ_fs ≥ λ_WDM(5.3 keV). La puerta Lyα habría que rehacerla con CLASS incluyendo esa población.

**Se informa además, sin veredicto:**
- la fila B2 con q = 1 (la cola tibia de U15b): cómo cambian su f_em y su m_G con K_eff. Las tablas de U15b se hicieron con K_sw ∈ {0, 1};
- si T_osc cae por debajo de 1 MeV (en plena nucleosíntesis) y por debajo de 0.1 MeV;
- el rango nuevo de m_G frente a la superradiancia y el rango de minicúmulos frente a U14.

## Lo que U17 no resuelve

- **c y η no se miden.** Hace falta una red con la ruptura −εQ₁₁ que llegue a H ≈ m_G (con paredes), fuera del alcance de este trabajo. U17 solo dice si el corchete de la cadena es compatible con la energía que la red tiene en ese momento.
- **El espectro del colapso** se resume en un momento medio.
- **No se toca** ningún resultado anterior ni el contrato hasta que haya informe.

## Nota 1 · 13 de octubre de 2026 (Madrid): literatura leída después del cálculo

arXiv:1202.5851 se leyó después de ejecutar U17 (`docs/informes/INFORME_U17.md`, Nota 1). Da ω̄/m ≈ 3.3 con η = 1 por hipótesis, es decir, c/η ≈ 3, dentro del recorrido preregistrado {3, 10, 30}. No se cambia nada del preregistro: el caso central sigue siendo c/η = 10 y el veredicto se lee con los tres valores.
