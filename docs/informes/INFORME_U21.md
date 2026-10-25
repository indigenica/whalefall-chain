# U21 · Hay una banda de F en la que el sector oscuro es frío, y el benchmark está justo fuera

*15–16 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U21.md` (commit `a7dedad`, anterior a cualquier fila).
- *Código:* `lya/u21_fband.py`. *Resultado:* `lya/results/u21/U21_RESULT.json`. En la nube, 25 filas de CLASS en 3.4 h.
- CLASS v3.4.1 con el parche y las puertas de U12 sin cambios.

## 0. En una frase

**Veredicto preregistrado: BANDA_BAJO_EL_BENCHMARK.** Con las hipótesis primarias de U17d (rama viva B2, colapso tipo U(1)), pasan todas las puertas todos los F de la rejilla entre 0.03·F₀ y **0.7·F₀**. Es decir, F ≈ 3×10¹³–7.5×10¹⁴ GeV, o μ/H_* ≈ 0.01–0.21.

El benchmark (F₀ = 1.07×10¹⁵ GeV) solo falla por la cota de isocurvatura de la literatura: m_G = 2.6×10⁻¹⁸ eV < 10⁻¹⁷ eV.

**El control de frío trivial pasa:** la fila P con F = 0.18·F₀ da δA = 0.006 con CLASS.

## 1. Bandas por variante

| Variante | Qué supone | Banda (F/F₀) | Fila clave |
|---|---|---|---|
| **P** (primaria) | B2, η = 1, R_ρ = 1.56, r = 1 | **0.03–0.7** | F = 0.7: δA = 0.185, ya por debajo del umbral estricto |
| H (colapso caliente, como el piloto de U19b) | r = 3 | 0.03–0.5 | F = 0.5: δA = 0.231 |
| D (colapso duro) | η = 0.3 | 0.03–1 (incluye el benchmark) | F = 1: δA = 0.116 |
| A (rama cómoda, ξ creciente) | B4a, η = 1 | 0.03–0.18 | F = 0.18: δA = 0.116 |

En las cuatro variantes la banda es contigua.

## 2. La variante primaria, fila a fila

| F/F₀ | m_G (eV) | T_osc | δA | forma (m, f) | isocurvatura (lit.) | superradiancia, sin autointeracción |
|---|---|---|---|---|---|---|
| 1 | 2.6×10⁻¹⁸ | 57 keV | 0.323 | 4.3 keV, 0.64 | **falla** | SMBH (indirecta) |
| 0.7 | 1.1×10⁻¹⁷ | 0.12 MeV | 0.185 | 3.0 keV, 0.15 | pasa | SMBH |
| 0.5 | 5.5×10⁻¹⁷ | 0.23 MeV | 0.094 | 10 keV, 0.40 | pasa | SMBH |
| 0.35 | 2.6×10⁻¹⁶ | 0.47 MeV | 0.043 | 10 keV, 0.13 | pasa | SMBH |
| 0.25 | 1.1×10⁻¹⁵ | 0.91 MeV | 0.018 | — | pasa | SMBH |
| 0.18 | 4.2×10⁻¹⁵ | 1.8 MeV | 0.006 (control) | — | pasa | SMBH |
| 0.125–0.09 | 2–7×10⁻¹⁴ | 4–7 MeV | frío trivial | — | pasa | SMBH |
| 0.0625 | 3.4×10⁻¹³ | 16 MeV | frío trivial | — | pasa | **fuera de las bandas** |
| 0.045–0.03 | 1–8×10⁻¹² | 31–72 MeV | frío trivial | — | pasa | estelar |

En toda la banda, ΔN_eff(emisión) ≤ 3×10⁻³ y F está por debajo del límite de cuerdas del CMB.

## 3. Lo que significa

- **El sector oscuro no está en el filo por la física: lo estaba por el benchmark.** F₀ estaba pegado al techo de C49 y del CMB, y ese techo es justo donde la materia oscura sale más caliente (m_G ∝ F⁻⁴).
- **Con F ≲ 0.5·F₀, la rama viva es fría** incluso con un colapso caliente (H), y con margen respecto al umbral estricto de δA.
- **La rama cómoda (ξ creciente) también vive con F ≲ 0.18·F₀.** La supervivencia del sector oscuro depende ahora mucho más de F que de η o de la estabilización de ξ.

**Lo que cuesta:**
- **LISA.** Con F < F₀, f_a = F/2 queda por debajo del umbral de `s5_closure` (5×10¹⁴ GeV): la señal de ondas gravitacionales de cuerdas deja de ser una predicción detectable para LISA.
- **Superradiancia (abierta).** Sin autointeracción, casi toda la banda cae en las bandas de agujeros negros: supermasivos (indirectas, 10⁻¹⁹–10⁻¹³ eV) para F ≳ 0.09·F₀, y estelares para F ≲ 0.045·F₀. Solo F ≈ 0.06·F₀ (m_G ≈ 3×10⁻¹³ eV) cae en el hueco. Con F ~ 10¹⁴–10¹⁵ GeV la autointeracción debilita esas cotas (Baryakhtar et al. 2021), pero aquí no se ha cuantificado. **Es la siguiente puerta que hay que calcular, y puede recortar mucho la banda.**
- **La isocurvatura es la de la literatura.** Hace falta U20, el cálculo propio.

## 4. Lo que U21 no resuelve

- **Superradiancia con autointeracción** (U23, por preregistrar).
- **Isocurvatura propia** (U20).
- **El colapso sigue con las hipótesis de U17d:** U19b medirá Y, y con Y la variante que aplica.
- **El ordenamiento de Q** se retrasa con μ/H_* pequeño: H_orden ∝ μ en C49. Las condiciones iniciales de la red cambian, y no se ha comprobado.

**Nota · 18 de octubre de 2026.** La superradiancia con autointeracción (U23), la isocurvatura observacional (U20b) y la población del colapso medida (U17e) dejan la banda en F ≈ 0.09–0.5·F₀ (`INFORME_U17e_U20b_U23.md`). El veredicto preregistrado de U21 no cambia.
