# U17e · La población del colapso medida en la red (U19b), pasada por CLASS: el borde de la banda de F · Preregistro

*17 de octubre de 2026, ≈19:22 (Madrid). Escrito antes de calcular ninguna fila de U17e (la hora que cuenta es la de su commit). Se ejecuta en la nube, con la caché de CLASS de U12/U21.*

**Lo que ya se ha visto, y se declara:**
- U19b completo con ℓ_G = 3.2 (`docs/informes/INFORME_U19b.md`).
- Al preparar este preregistro se calculó la forma de la población blanda de las cuatro redes con la regla 2 de abajo: K_blando = 0.94–1.31 por cuanto de desalineamiento y x̄ = 5.5–6.2. No se ha calculado ninguna masa, ninguna longitud de libre recorrido ni ninguna fila de CLASS.
- U19b con ℓ_G = 2.4 no se ha visto. Su primera red está corriendo en WSL.

## Motivo

U17d, U21 y U20 describen el colapso cuerdas–paredes con una log-normal postulada: momento medio c·m_G·r, con c = 3 y σ_ln = 0.5. U19b la midió en RP³:
- Y = 1.70 ± 0.55;
- un espectro de momentos más ancho y más caliente que la log-normal (c_p ≈ 4.3, r_W ≈ 1.8);
- un número de cuantos blandos que no sigue a la densidad de cuerdas (U19b §4, enmienda 4).

ESPEC_U19b anunció que su espectro final sería la entrada de una U17e. Aquí se fija cómo.

## La población medida (`lya/u17e_class.py`, función `measured_shape`)

De cada una de las cuatro redes primarias (RP³, ℓ_G = 3.2):
1. **Exceso blando:** el espectro final del contador lineal, dN por bin de k, con k < 10·m_G·a₁. Se le resta el espectro cinético radiado antes de ℓ₁, interpolado en ℓ₁ y rebinado de forma acumulada a los mismos bins. Es la misma resta que hacen X y K en U19b (enmienda 3).
2. **El desalineamiento se quita de los bins de k más bajo,** hasta sumar N_mis cuantos. Es la atribución conservadora: deja la población del colapso tan caliente como permite el espectro. Los bins negativos se ponen a cero y se informa su peso.
3. **Forma:** la fracción de número por bin, w(x), con x = k/(m_G·a₁). Se promedia entre las cuatro redes con el mismo peso. Se informan x̄ = Σw·x y K_blando = (cuantos del exceso)/N_mis, media de las cuatro.

**En el universo,** cada bin da cuantos de momento comóvil q = x·m_G·a_osc, con a_osc en H = m_G/3. El espectro de la red es el final, así que la época del colapso ya está dentro y no se usa r. Los cuantos duros de la red (k ≥ 10·m_G·a₁) se tratan como en U17c: en el universo van a p ≫ m_G y su número es despreciable. Se informa una fila de sensibilidad con ellos.

## Filas

**Rejilla:** F/F₀ ∈ {1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03}, la de U21, con q = 0.87, x₀ = 10 y la emisión de U12.

**Tipos:**
- **Y (primaria):** la rama viva B2, con el escalado de U17, el que preregistró U19b. El número es K = (8π/27)·R_ρ·ξ_B2·ℓ·Y/x̄, con R_ρ = 1.56 y **Y = 1.70**, lo que equivale a `c_eta = x̄/Y` en `u17_ksw.solve`. m_G se resuelve de forma autoconsistente.
- **Y+ (conservadora):** igual, con Y + 1σ = 2.24.
- **X (lectura (b) de la enmienda 4 de U19b):** el número del colapso es el medido, sin escalar con ξℓ: K = K_blando, con la misma forma.
- **Informativas, sin veredicto:**
  - **A:** la rama cómoda B4a con el tipo Y;
  - **Yp:** el tipo Y con el desalineamiento quitado en proporción, no de los bins bajos, con F ∈ {0.5, 0.35};
  - **Xd:** el tipo X con los cuantos duros de la red tal como salen, con F ∈ {1, 0.5}.

**Regla de frío trivial (la de U21):** si λ_proxy/λ_WDM(5.3 keV) < 0.05, la fila pasa Lyα sin CLASS, con λ_proxy la longitud de libre recorrido con el momento x̄·m_G·a_osc. La primera fila Y que la cumpla se calcula también con CLASS como control, y debe dar δA < 0.05. Si el control falla, se calculan con CLASS todas las filas que la cumplían.

## Puertas por fila

1. **Lyα:** la puerta primaria de U12, sin cambios.
2. **Isocurvatura propia:** el estimador de U20 con C = 1 y el criterio Δ²_iso/Δ²_ad ≤ 0.1 en k = 10 h/Mpc, z = 3. k_* = a_osc·H_osc/r con **r = r_N = 2.40**, la media de U19b. Es la época más tardía medida y, por tanto, la conservadora. Se informa también con r = r_W = 1.78.
3. **Radiación oscura:** ΔN_eff(emisión) ≤ 0.3.
4. **Cuerdas en el CMB:** F ≤ 1.26×10¹⁵ GeV.

Se informan, sin puerta, LISA (f_a = F/2 ≥ 5×10¹⁴ GeV), m_G, T_osc y f_nc.

## Veredicto (filas Y)

F_max es el F más alto de la rejilla que pasa las cuatro puertas.
- **BANDA_U21_CONFIRMADA:** F_max ≥ 0.5·F₀.
- **BANDA_RECORTADA:** 0.18·F₀ ≤ F_max ≤ 0.35·F₀.
- **BANDA_MINIMA:** algún F pasa, pero F_max < 0.18·F₀.
- **SIN_BANDA:** ningún F pasa.

**Se informan también:**
- la banda conservadora (filas Y+);
- **la banda de las filas X.** Solo se aplica si la enmienda 4 de U19b da la lectura (b) con las redes de ℓ_G = 2.4. Si con (b) pasa F₀, el benchmark (y con él LISA) vuelve a estar permitido; se dirá como «BENCHMARK_PERMITIDO_SI_(b)»;
- las informativas A, Yp y Xd.

## Lo que se espera, escrito antes de calcular

- **Filas Y:** por Lyα, F hasta 0.5–0.7·F₀, entre las variantes P y H de U21. Por isocurvatura con r = 2.4, el cociente en 0.5·F₀ sube un factor ~14 respecto a la variante P de U20 (0.025), y fallaría. **Se espera BANDA_RECORTADA con F_max = 0.35·F₀.**
- **Filas X:** K ≈ 1.1 en vez de ~50. Para dar Ω = 0.12, m_G sube unas 600 veces (Ω_mis ∝ m^½·(1 + K)). Así, la isocurvatura y Lyα pasarían también en F₀. **Se espera que pasen en toda la rejilla, incluido F₀.**

## Lo que U17e no resuelve

- **La extrapolación de ℓ₁ ≈ 4.3 a ℓ ≈ 90.** El tipo Y supone el escalado de U17, y el tipo X, que el rendimiento blando no crece con ℓ. Decide U19b con ℓ_G = 2.4.
- **La normalización C de la isocurvatura** sigue tomándose como 1.
- **La forma de la población** sale de redes de ℓ₁ = 4.3. La parte de momentos altos, cerca del corte k = 10·m_G·a₁, depende de ese corte.
- **La superradiancia con autointeracción** (U23) sigue pendiente.
