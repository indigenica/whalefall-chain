# U17 · El colapso de la red en H ≈ m_G no cabe en K_sw ∈ {0, 1}

*12 de octubre de 2026, noche (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17.md` (commit anterior al cálculo).
- *Código:* `chain/u17_ksw.py` (segundos). *Resultado:* `results/u17_ksw.json`.

## 0. En una frase

**Veredicto preregistrado: K_SW_FUERA_DEL_RANGO_TIBIA.** Cuando aparecen las paredes (H = m_G/3), la red guarda entre 100 (ξ constante) y ~2000 (ξ creciente) veces la energía del desalineamiento. Si esa energía acaba en Goldstones blandos, K_sw no es 0–1 sino 3–770.
- En la rama de ξ constante (B2), la materia oscura sigue siendo fría, pero m_G baja dos o tres órdenes de magnitud y los minicúmulos suben otros tantos.
- En las ramas en que ξ crece (B4a, B4b, B1), m_G cae a 10⁻²⁰–5×10⁻¹⁸ eV, la oscilación empieza a T ≈ 4–80 keV, y los cuantos del colapso, que pasan a ser casi toda la materia oscura, recorren 5–13 veces más que la WDM de 5.3 keV. **La conclusión de U12 («la rama fría sobrevive aunque ξ crezca») deja de estar asegurada:** se apoyaba en K_sw ≤ 1.

## 1. La cuenta

En H = m_G/3, la densidad de la red es ρ_red = 4ξπf_a²ℓH², y la del desalineamiento, ρ_mis = (3/8)m_G²F² con F = 2f_a. Si una fracción η de ρ_red acaba en Goldstones de energía media c·m_G:

**K_eff = n_colapso/n_mis = (8π/27)·ξ(ℓ)·ℓ·η/c.**

Comprobación a mano en B2 con ℓ = 89.6: ρ_red/ρ_mis = 100.0, y K_eff = 10.0 con c/η = 10, igual que el código.

m_G se resuelve de forma autoconsistente con la emisión de U12 sin tocar: Ω_mis(1 + K_eff) + Ω_em = 0.12, con K_eff recalculado en cada paso del buscador de raíces.

## 2. Resultados (R_Γ = 1.56, x₀ ∈ {10, 30}, q ∈ {0.87, 1})

| c/η | K_eff | m_G (eV) | T_osc (MeV) | Minicúmulos (M_⊙) | λ_colapso/λ_WDM(5.3 keV) |
|---|---|---|---|---|---|
| K_sw ∈ {0, 1} (U12) | 0–1 | 1.6×10⁻¹⁷–9.6×10⁻¹⁵ | 0.14–2.7 | 0.2–10⁴ | — |
| 30 | 3.3–70 | 8×10⁻¹⁹–5×10⁻¹⁶ | 0.03–0.64 | 19–3×10⁵ | 0.65–12.8 |
| 10 | 10–220 | 1.2×10⁻¹⁹–7×10⁻¹⁷ | 0.013–0.26 | 390–5×10⁶ | 0.60–11 |
| 3 | 34–770 | 1.2×10⁻²⁰–5×10⁻¹⁸ | 0.004–0.08 | 2×10⁴–2×10⁸ | 0.66–10.5 |

**Por rama, con c/η = 10:**

| Rama | K_eff | m_G (eV) | T_osc (keV) | λ/λ_WDM | f_em (q = 0.87 / q = 1, x₀ = 10) |
|---|---|---|---|---|---|
| B2 | 10 | 6–7×10⁻¹⁷ | 240–260 | 0.60–0.66 | 2×10⁻⁵ / 0.086 |
| B4a | 109 | 5–6×10⁻¹⁹ | 25–27 | 5.5–6.0 | 1×10⁻⁵ / 0.089 |
| B4b | 219 | 1.2–1.5×10⁻¹⁹ | 13–14 | 10–11 | 9×10⁻⁶ / 0.089 |
| B1 | 204 | 1.4–1.7×10⁻¹⁹ | 14–15 | 9.6–10.4 | 1×10⁻⁵ / 0.089 |

- **Condición (i) del veredicto** (cuantos del colapso más fríos que la WDM de 5.3 keV): se cumple en B2 para los tres c/η y **falla en B4a, B4b y B1 para los tres**.
- **Condición (ii)** (f_em < 1 % con q = 0.87): se cumple en todas las filas. La emisión de las cuerdas sigue siendo despreciable con q < 1; lo que cambia es la otra fuente.

## 3. Lectura

- **Por qué es robusto en las ramas que crecen.** El cálculo toma η = 1, el lado tibio: con η < 1 los cuantos serían más blandos para el mismo c/η. Pero el problema no viene sobre todo de c, sino de que m_G baja tanto que la oscilación empieza a T ≈ 10 keV, cuando el horizonte comóvil es grande. λ ∝ c/T_osc: incluso con p ≈ m_G (c ≈ 1), en B4a λ/λ_WDM ≈ 2 y en B1 ≈ 3. Solo se salvarían si casi toda la energía de la red (c/η ≳ 10³) se fuera a cuantos duros también en el colapso, lo contrario de lo que se espera de paredes de grosor 1/m_G.
- **ξ constante (B2) sobrevive fría**, pero sus predicciones cambian:
  - m_G ≈ 5×10⁻¹⁸–5×10⁻¹⁶ eV en lugar de 2–10×10⁻¹⁵ eV;
  - minicúmulos de 20–2×10⁴ M_⊙ en lugar de 0.2–2 M_⊙ (U14 se calculó para 0.23–10 M_⊙);
  - con q = 1, la cola tibia baja de f_em = 0.15–0.52 a 0.01–0.20 (las tablas de U15b se hicieron con K_sw ∈ {0, 1}).
- **La literatura de U(1) apunta a c/η del orden de decenas** *(incorrecto: corregido en la Nota 1, con el artículo ya leído; la literatura da c/η ≈ 3)*. Para el axión QCD con N_DW = 1, las contribuciones de cuerdas y paredes suelen salir comparables al desalineamiento, o unas pocas veces mayores, con ξ ≈ 1 y ℓ ≈ 70. Eso corresponde a c/η ≈ 20–30 en esta fórmula. No se han podido verificar los números en esta sesión (sin acceso a arXiv).
- **Con q < 1 y ξ creciente, la pregunta vuelve a ser la estabilización de ξ.** U12 había mostrado que la emisión de las cuerdas no la exige. U17 dice que la energía que queda en la red al final sí la exige, salvo que c/η ≳ 10³.

## 4. Lo que U17 no resuelve

- **c y η no están medidos** para una red ℤ₂ de RP³ con la ruptura −εQ₁₁. Haría falta una red que llegue a H ≈ m_G con paredes.
- **El espectro del colapso se resume en un momento medio**, y la comparación con Lyα es la de longitudes de libre recorrido de la cadena (P5), no una puerta con CLASS. **Siguiente paso:** pasar por CLASS (la puerta δA y la de forma de U5) las filas B2, B4a y B1 con q = 0.87 y c/η ∈ {3, 30}, con la población del colapso como un paquete más. Necesita su propio preregistro (U17b).
- **m_G ≈ 10⁻²⁰–10⁻¹⁹ eV** entra en la zona de la materia oscura ultraligera, con sus propias cotas (Lyα, m ≳ 2×10⁻²⁰ eV si es toda la materia oscura). No se han aplicado.
- **No se cambia el contrato.** U17 no reclasifica U12. Señala que su veredicto depende de un corchete de K_sw que la energía de la red no respalda, y que conviene decirlo en el artículo.

## Nota 1 · 13 de octubre de 2026 (Madrid): lectura de arXiv:1202.5851, después del cálculo

El equipo aportó el PDF de Hiramatsu, Kawasaki, Saikawa y Sekiguchi, «Production of dark matter axions from collapse of string-wall systems» (PRD 85, 105020; arXiv:1202.5851v3), que no se había podido descargar. Se leyó después de calcular U17. **No cambia ni el criterio ni el veredicto**, pero corrige la frase de §3 sobre la literatura.

**Lo que dice, literalmente** (U(1), N_DW = 1, red 3D desde la transición de PQ hasta el colapso):
- ec. (77)–(78): con κ = 0.3, k̄/R = 5.76 ± 0.15 y ε_w = 3.12 ± 0.08, «the mean energy of axions» ω̄_a/m_a = √(1 + ε_w²) = 3.28 ± 0.08. En la tabla III, ε_w = 2.51–3.12 para κ = 0.4–0.3. La figura 6 muestra que la mayoría de los axiones tienen k/R ~ m_a: «not highly relativistic».
- ec. (90)–(91): «we simply assume that, after the time t₁, the whole energy stored in these defects is diluted as R⁻³», y el número de axiones del colapso es ρ_cuerdas-paredes/ω̄_a. Es decir, **η = 1 por hipótesis**.
- ec. (84) y (94): Ω_mis h² ≈ 0.58 y Ω_dec h² = 5.8 ± 2.8 (con ξ₁ ≈ 1.0 ± 0.5 y A₁ ≈ 0.50 ± 0.25): el colapso da **unas 10 veces** el desalineamiento (incluido el factor anarmónico de 1.85 en este último), y más que las cuerdas antes del colapso (2.0 ± 1.0).

**Lectura para U17:**
- Con su método, c/η ≈ 3.3: es el extremo c/η = 3 del recorrido preregistrado, el más tibio. El caso central de U17 (c/η = 10) es, por tanto, optimista respecto a la literatura.
- Con c/η = 3: B2 da K_eff ≈ 34, m_G ≈ 5×10⁻¹⁸ eV, T_osc ≈ 80 keV y λ/λ_WDM ≈ 0.66–0.68. Sigue siendo fría, pero con minicúmulos de ~2×10⁴ M_⊙. Las ramas que crecen dan K_eff ≈ 380–770, m_G ≈ 1–5×10⁻²⁰ eV y λ/λ_WDM ≈ 5.5–10.5.
- **Su cociente colapso/desalineamiento (~10 con ξ ≈ 1) confirma el orden de magnitud de K_eff** y descarta que K_sw ∈ {0, 1} sea un corchete conservador. El 0.27 «medido» de C40 queda un orden de magnitud por debajo.
- Las diferencias con nuestro caso, que se declaran: su vacío es U(1), no la RP³ de la ballena; ξ y A₁ son estimaciones con un 50 % de error; y η = 1 es una hipótesis suya, no una medida.

## Nota 2 · 14 de octubre de 2026 (Madrid): energía y momento (auditoría de Codex)

ESPEC_U17 define c como la **energía** media por cuanto: el número de cuantos es η·ρ_red/(c·m_G). El cálculo de la longitud de libre recorrido (§3 del preregistro y `solve`), en cambio, usaba **momento** c·m_G. Se corrige a √(c² − 1)·m_G (η = 1), que es lo coherente con la definición. Cambia poco: con c = 3 la longitud baja un 6 %, y con c ≥ 10, menos de un 1 %. **El veredicto (K_SW_FUERA_DEL_RANGO_TIBIA) no cambia:** las ramas que crecen siguen con λ/λ_WDM > 1. `results/u17_ksw.json` se regeneró con `make chain`.

Hay una segunda cuestión que este informe no trataba: el inventario usa ξ de U(1), sin el cociente de energías RP³/U(1) (R_ρ ≈ R_Γ = 1.56). Va como hipótesis declarada en U17d (`docs/preregistros/ESPEC_U17d.md`).
