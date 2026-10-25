# U12 · La rama fría sobrevive aunque la densidad de cuerdas crezca

*12 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U12.md` (commit `b0d2bc3`, 11 de octubre a las 19:35, antes de calcular), con las notas 1–3.
- *Resultado:* `lya/results/u12/U12_RESULT.json`, ejecutado en WSL en 1.26 h. CLASS v3.4.1 con el parche y las mismas puertas que U5.
- *Filas:* 64. Ramas B2 (ξ constante, control), B4a y B4b (crecimiento medido, s = 0.133) y B1 (ξ ∝ ℓ, GHV); q ∈ {0.6, 0.75, 0.87, 1.0}; R_Γ = 1.56 (U11); x₀ ∈ {10, 30}; K_sw ∈ {0, 1}.

## 0. En una frase

**Veredicto preregistrado: FRIA_SOBREVIVE_AL_CRECIMIENTO.** Con q ≤ 0.87, las tres ramas en las que la densidad crece (B4a, B4b y también B1, el crecimiento logarítmico de GHV) pasan Lyα en las 4 filas de cada casilla: la emisión es menor que el 1 % de la materia oscura, que queda fría. Con q = 1 las tres fallan en todas las filas. **El sector oscuro deja de depender de que ξ se estabilice; ahora depende solo de q a ℓ ≈ 70**, con un q crítico de 0.91–0.92 para las ramas que crecen.

## 1. Tabla (filas permitidas por la puerta primaria, de 4)

| Rama | q = 0.6 | q = 0.75 | q = 0.87 | q = 1.0 | q_crit (x₀ = 10, K_sw = 0) |
|---|---|---|---|---|---|
| B2 (ξ* = 1.20, control) | 4 | 4 | 4 | 4 | — (pasa en todo el rango) |
| B4a (s = 0.133, nuestra normalización) | 4 | 4 | 4 | **0** | 0.923 |
| B4b (s = 0.133, normalizada a 1.2) | 4 | 4 | 4 | **0** | 0.912 |
| B1 (ξ = 0.24 ℓ, GHV) | 4 | 4 | 4 | **0** | 0.913 |

Con solo δA el resultado es el mismo en las 16 casillas.

**El control B2 con q = 1** reproduce lo esperado: δA = 0.12–0.17 y f_em = 0.15–0.52, como en U5 con R_Γ = 1.40 (0.10–0.17). Los ajustes de forma caen en un 2.0–2.8 % de 0.24–0.33 keV, la cola tibia de siempre.

## 2. Lo que dejan las filas que sobreviven

| Rama, q ≤ 0.87 | f_em | m_G (eV) | Minicúmulos (M_⊙) | ΔN_eff (CMB) |
|---|---|---|---|---|
| B2 | ≤ 0.0004 | 2.4–9.6×10⁻¹⁵ | 0.24–2.0 | ≤ 0.005 |
| B4a | ≤ 0.0034 | 2.4–9.6×10⁻¹⁵ | 0.24–2.0 | ≤ 0.034 |
| B4b | ≤ 0.0066 | 2.4–9.6×10⁻¹⁵ | 0.25–2.0 | ≤ 0.067 |
| B1 | ≤ 0.0062 | 2.4–9.6×10⁻¹⁵ | 0.24–2.0 | ≤ 0.063 |

- **Materia oscura fría, de desalineamiento.** La emisión aporta como mucho un 0.7 %. m_G queda en la ventana fría de U11 (2–10×10⁻¹⁵ eV), por debajo de la superradiancia estelar.
- **Una firma nueva: radiación oscura.** Cuando la red crece y q < 1, la energía emitida se va a Goldstones duros, que no hacen materia oscura pero sí radiación. ΔN_eff llega a 0.03 (B4a) y a 0.06–0.07 (B4b, B1). Con σ(ΔN_eff) ≈ 0.03 de CMB-S4, la rama de crecimiento logarítmico daría una señal de ~2σ. Con ξ constante (B2), en cambio, ΔN_eff ≤ 0.005 es invisible.

## 3. Lectura

- **La revisión externa tenía razón.** La condición «la densidad debe estabilizarse» de U5 valía para q ≥ 1, no en general. Con q < 1, ni el crecimiento de GHV mata el sector oscuro: lo deja frío.
- **La pregunta se reduce a q a ℓ ≈ 70.** El umbral para las ramas que crecen es q ≈ 0.91–0.92; para ξ constante, q_crit = 1.025 (U5f).
- **Nuestra medida de q está justo en esa frontera.** U11 da q(ℓ = 5) = 0.87 ± 0.06 con el convenio original. La auditoría U4–U11 encuentra 0.87–0.94 según cómo se convierta la ventana a momentos de la red. Si la densidad crece, decide si q a ℓ ≈ 70 queda por debajo de ~0.91. Si se estabiliza, basta con q ≤ 1.02.
- **Las predicciones se separan en tres casos:**
  - q ≈ 1 y ξ constante: cola tibia marginal en Lyα (δA ≈ 0.12–0.17);
  - q < 0.91 y ξ creciente: materia oscura fría y ΔN_eff ≈ 0.03–0.07;
  - q < 1 y ξ constante: materia oscura fría sin firma en Lyα ni en ΔN_eff.

## 4. Lo que U12 no resuelve

- q y ξ a ℓ ≈ 70 siguen fuera del alcance de la red.
- Solo se usa R_Γ central (1.56). Con q < 1, la emisión pesa tan poco que ±0.09 no cambia nada. Con q = 1 y ramas que crecen, el fallo es enorme (δA ≈ 0.54–0.72 frente a 0.234), así que tampoco cambiaría.
- La cota de ΔN_eff no es un veredicto de U12. Planck (ΔN_eff < 0.3) permite todos los casos.
