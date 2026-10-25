# U20 · La isocurvatura propia confirma la cota de la literatura y recorta la banda a F ≤ 0.5·F₀

*16 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U20.md` (commit `9122623`, anterior a la ejecución).
- *Código:* `chain/u20_iso.py`. *Resultado:* `results/u20_iso.json`. En segundos.

## 0. En una frase

**Veredicto preregistrado: ISOCURVATURA_RECORTA_BANDA.** Con una estimación propia de la isocurvatura de ruido blanco (normalización C = 1, sin libre propagación, criterio de que sea ≤ 10 % de la potencia adiabática en k = 10 h/Mpc, z = 3), la rama viva pasa con **F ≤ 0.5·F₀** (U21 daba 0.7·F₀). La masa mínima, **m_min ≈ 2×10⁻¹⁷ eV**, coincide en orden de magnitud con la cota de la literatura (10⁻¹⁷ eV).

## 1. Filas (variante P, C = 1; Δ²_ad = 1.16 en k = 10 h/Mpc, z = 3; D_iso = 1149)

| F/F₀ | m_G (eV) | k_* (Mpc⁻¹) | Δ²_iso/Δ²_ad | ¿pasa? |
|---|---|---|---|---|
| 1 | 2.6×10⁻¹⁸ | 541 | 2.2 | no |
| 0.7 | 1.1×10⁻¹⁷ | 1140 | 0.23 | no |
| **0.5** | 5.5×10⁻¹⁷ | 2420 | **0.025** | **sí** |
| 0.35 | 2.6×10⁻¹⁶ | 5220 | 2.4×10⁻³ | sí |
| ≤ 0.25 | ≥ 10⁻¹⁵ | ≥ 10⁴ | ≤ 3×10⁻⁴ | sí |

## 2. Masa mínima según la normalización C

| C | variante P (r = 1) | variante H (r = 3, colapso caliente) |
|---|---|---|
| 0.3 | 9×10⁻¹⁸ eV | 9×10⁻¹⁷ eV |
| **1** | **2.1×10⁻¹⁷ eV** | **2.0×10⁻¹⁶ eV** |
| 3 | 4.4×10⁻¹⁷ eV | 4.2×10⁻¹⁶ eV |
| 10 | 1.0×10⁻¹⁶ eV | 9.4×10⁻¹⁶ eV |

La banda con C = 1 queda así:
- **variante P:** F ≤ 0.5·F₀;
- **variante H:** F ≤ 0.35·F₀. Un colapso tardío baja k_*, y la isocurvatura pesa más.

## 3. Qué significa

- **La puerta de isocurvatura de U21 ya no depende solo de una cita.** Nuestra estimación de orden de magnitud da la misma escala: m ≳ 10⁻¹⁷ eV.
- **La banda de la rama viva** queda en **F ≲ 0.5·F₀ (≈ 5×10¹⁴ GeV)**, y en ≲ 0.35·F₀ si el colapso de RP³ es tan caliente como sugiere el piloto de U19b. Dentro de ella sigue quedando mucho margen: 0.03–0.5·F₀.
- **La incertidumbre principal es C,** la amplitud del ruido blanco, que depende de la física del colapso. U19b puede medirla con el espectro de densidad de la red.

## 4. Lo que U20 no resuelve

- **La normalización C** no está calculada.
- **La libre propagación** de la población del colapso, que debilitaría la cota, no se incluye.
- **No hay espectro de flujo de Lyα ni verosimilitud:** el criterio del 10 % es aproximado.
- **La superradiancia con autointeracción** sigue abierta (U23).

**Nota · 18 de octubre de 2026.** El criterio de isocurvatura de U20 era unas 500 veces más estricto que las cotas observacionales. U20b lo sustituye como mejor estimación (ISOCURVATURA_NO_LIMITA), y la banda combinada es F ≈ 0.09–0.5·F₀ (`INFORME_U17e_U20b_U23.md`). El veredicto preregistrado de U20 no cambia.
