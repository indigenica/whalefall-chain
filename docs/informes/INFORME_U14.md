# U14 · Microlente de los minicúmulos de la ballena

*11 de octubre de 2026, noche.*
- *Preregistro:* `docs/preregistros/ESPEC_U14.md` (commit `8bbe43f`, 22:09, antes de ejecutar código).
- *Código y resultado:* `lya/u14_microlens.py` → `lya/results/u14/U14_RESULT.json`. Se ejecuta en un segundo.

## 0. En una frase

**Veredicto preregistrado: MICROLENTE_ACOTA_f_MC.**
- Los minicúmulos de la ballena (0.23–10 M_⊙) son casi todos demasiado difusos para hacer microlente: su radio típico es de miles de unidades astronómicas, frente a un radio de Einstein de 5–30 UA.
- Solo la cola más densa (δ ≳ 65–105, un 12–15 % de la masa según el ajuste de Kolb y Tkachev) actúa como lente puntual.
- Con la cota estricta de EROS-2 (f < 0.08), eso exige que **menos del ~55 % de la materia oscura esté en minicúmulos**. Con la cota combinada EROS-2 + MACHO (f < 0.20) no hay restricción.

## 1. Números (radio de esfera dura)

| M (M_⊙) | R_E (UA) | r_mc con δ = 1 (UA) | δ_lens | F(δ > δ_lens) | f_MC,max: estricta / amplia |
|---|---|---|---|---|---|
| 0.23 | 4.8 | 1 901 | 105 | 0.117 | 0.68 / 1.70 |
| 0.5 | 7.1 | 2 463 | 95 | 0.123 | 0.65 / 1.63 |
| 1 | 10.1 | 3 103 | 87 | 0.128 | 0.63 / 1.57 |
| 2 | 14.3 | 3 909 | 80 | 0.133 | 0.60 / 1.51 |
| 5 | 22.6 | 5 305 | 71 | 0.140 | 0.57 / 1.43 |
| 10 | 31.9 | 6 684 | 65 | 0.146 | **0.55** / 1.37 |

**Variantes de radio (informativas):**
- **con r × 3:** f_MC,max = 0.80–1.00 (estricta) y 2.0–2.5 (amplia);
- **con r / 3:** f_MC,max = 0.38–0.47 (estricta) y 0.95–1.17 (amplia). En este caso, a 5–10 M_⊙ la cota amplia ya morde (f_MC,max = 0.95–0.98).

## 2. Lectura

- **La microlente no falsa la ballena, pero sí la condiciona.** Si la mayor parte de la materia oscura acabara en minicúmulos (f_MC ≈ 0.7–0.8, lo típico en simulaciones de axiones QCD), la cola densa daría más eventos que los que EROS-2 permite a M ≈ 0.4 M_⊙. Esa conclusión depende de tres cosas que U14 no fija.
- **La cola de δ.** El ajuste F ∝ δ^(−0.46) está extrapolado a δ ~ 100, y es de axiones con masa dependiente de la temperatura. Para una masa constante (n = 0), como la de m_G, no hay simulaciones que conozcamos.
- **La cota de EROS-2 en masa.** 0.08 es el valor a 0.4 M_⊙. A 5–10 M_⊙ la cota de EROS-2 es más débil, así que aplicarla en todo el rango es conservador.
- **La destrucción por mareas y por encuentros estelares** reduce la fracción de minicúmulos intactos en el halo; no está incluida, y aflojaría la cota.
- **Para la predicción observable:** la ballena predice minicúmulos de masa solar. La parte densa es, en principio, visible en microlente (eventos de 1–3 meses hacia la LMC) y en cruces de cáusticas de estrellas muy magnificadas. Es una firma adicional, independiente de Lyα, que conviene listar.

## 3. Qué queda

- La distribución de δ y la función de masas para una masa constante (n = 0): necesitaría simulaciones propias o literatura específica.
- f_MC en la ballena: requiere seguir la formación de minicúmulos a partir del campo de desalineamiento inhomogéneo.
- Los cruces de cáusticas y los arreglos de púlsares, para minihalos de masa solar.
