# U12 · ¿Sobrevive la rama fría si la densidad de cuerdas sigue creciendo? · Preregistro

*11 de octubre de 2026, 19:40 (Madrid). Escrito antes de calcular ninguna fila de U12 (ni con CLASS ni con el mock).*

## Motivo

U5 dejó dos condiciones para el sector oscuro: que ξ deje de crecer (U5a) y que q ≤ 1.02 (U5f). Solo se probó parte del cruce entre las dos:
- en U5f, la rama de crecimiento medido B4a (s = 0.105) pasa Lyα con q = 0.9 y falla con q = 1.1 y 1.2;
- nunca se probó la rama B1 (ξ ∝ ℓ, normalización de GHV) con q < 1, ni B4b, ni la R_Γ limpia de U11.

U11 midió q_RP3 = 0.87 ± 0.06 en ℓ = 5, bajando con ℓ, y R_Γ = 1.56 ± 0.09 sin el baño de la condición inicial. La revisión externa (tarea 6 de `docs/PLAN_NUBE_2026-10-11.md`) pide la cuarta casilla: densidad creciente con q < 1. Si la rama fría sobrevive aunque ξ crezca, el sector oscuro deja de depender de la extrapolación más discutida de ξ.

## Filas

Todo igual que en U5 (`lya/u5_class.py`), salvo las ramas, q y R_Γ:
- la misma cadena S6 (`chain/s6_measured.node` y `source_eps`) y el mismo benchmark;
- CLASS v3.4.1 con el parche de cuadratura (`ncdm_quadrature_strategy = 4`), con la comprobación C0b obligatoria;
- las mismas referencias térmicas, el mismo umbral estricto δA (el mínimo de las referencias, ≈ 0.234) y la misma familia de forma ampliada de U5b (m ≥ 0.2 keV, f ≥ 0.01);
- la misma puerta primaria: forma estricta si el ajuste es representable (rms ≤ 0.1), si no δA estricto.

**Ramas** (ε(ℓ) como en `chain/s6_measured.py` y `lya/u5_class.py`):

| Rama | ξ(ℓ) | ε(ℓ) |
|---|---|---|
| B2 (control) | ξ* = 1.20 constante | R · 1.20 / (0.24 ℓ) |
| B4a | 0.616 + s (ℓ − 5.65), s = 0.133 | R · max(ξ, 0.3) / (0.24 ℓ) |
| B4b | (1.20/0.616) · [0.616 + s (ℓ − 5.65)], s = 0.133 | igual, con la normalización de B4b |
| B1 | 0.24 ℓ (GHV) | R |

s = 0.133 es la pendiente medida en U5d (cuatro semillas U(1) a 768³, 0.133 ± 0.030).

**Pendiente espectral:** q ∈ {0.6, 0.75, 0.87, 1.0}, con el espectro `powerlaw_q{q}` (dE/dln x ∝ x^(1−q) en [x₀, m_r/H]; q = 1 es el log-plano).

**Emisión:** R_Γ = 1.56 (U11, central), x₀ ∈ {10, 30}, K_sw ∈ {0, 1}.

**Total:** 4 ramas × 4 q × 2 x₀ × 2 K_sw = **64 filas**. Las filas con f_em < 10⁻⁴ no pasan por CLASS (T² = 1, como en U5).

Comando: `python -B lya/u12_class.py --workers 32` → `lya/results/u12/U12_RESULT.json`.

## Veredicto

Para cada rama y cada q, el estado es:
- **PASA:** las 4 filas (x₀, K_sw) permitidas por la puerta primaria;
- **PARCIAL:** entre 1 y 3;
- **FALLA:** ninguna.

El veredicto global se decide con las tres ramas que crecen (B4a, B4b, B1), en este orden:
1. **FRIA_SOBREVIVE_AL_CRECIMIENTO:** B4a, B4b y B1 están en PASA con q = 0.87.
2. **FRIA_SOBREVIVE_AL_CRECIMIENTO_MEDIDO:** B4a y B4b están en PASA con q = 0.87, pero B1 no.
3. **SOLO_CON_q_BAJA:** B4a y B4b no están ambas en PASA con q = 0.87, pero sí con q = 0.6 o con q = 0.75.
4. **DEPENDE_DE_LA_ESTABILIZACION:** en ningún q ≤ 0.75 están B4a y B4b ambas en PASA.

B2 es el control y no entra en el veredicto. Se espera que q = 1.0 con R_Γ = 1.56 se parezca a U5 y U7c (pasa, cerca del límite). Si no fuera así, se informa como discrepancia antes de leer el veredicto.

## Se informa también, sin veredicto

- La tabla completa rama × q con el número de filas permitidas por la puerta primaria y, aparte, solo por δA estricto. U7c mostró que en este régimen, por debajo de 1 keV, la puerta de forma no es monótona.
- q_crit por rama: interpolación lineal de δA frente a q en la fila central (x₀ = 10, K_sw = 0) donde cruza el umbral estricto.
- f_em, m_G y ΔN_eff de cada fila, y el ajuste de forma (m, f), señalando si cae en el borde de la malla o por debajo de 1 keV.
- Si la rama superviviente es fría (f_em < 1 %) o conserva una cola tibia medible.

## Lo que U12 no resuelve

- **q a ℓ ≈ 70.** U12 recorre q como parámetro. q = 0.87 es el valor drenado a ℓ ≈ 5, no la pendiente cosmológica; otros protocolos dan 0.6–1.0.
- **La cota de Lyα por debajo de 1 keV.** Las puertas son las de U4a y U5. La verosimilitud general (arXiv:2206.08188) es la tarea 2 del plan y tiene su propio preregistro.
- **R_Γ fija.** Solo se usa el central de U11. Con q < 1 la emisión pesa poco y se espera que el error de ±0.09 importe menos que q, pero U12 no lo comprueba.

## Nota 1 · 11 de octubre de 2026, 20:15 (Madrid)

Después de fijar este preregistro (commit `b0d2bc3`) se ejecutó `lya/u12_class.py --quick` con el mock de CLASS (`tests/mock_classy`), solo para comprobar la tubería, en la sesión de la nube, que no pudo instalar CLASS (sin acceso a PyPI). Los δA y ajustes de forma del mock no valen y no se usan. Lo único que el mock no toca son las cantidades de la cadena S6 (f_em, m_G), que ya estaban a la vista: con q ≤ 0.87, f_em < 1 % en las cuatro ramas en la fila central. No se ha cambiado nada del preregistro. Su resultado (`lya/results/u12/U12_QUICK.json`) no se versiona.

## Corrección de horas · 11 de octubre de 2026, 20:07 (Madrid)

Las horas escritas a mano en la cabecera (19:40) y en la nota 1 (20:15) estaban mal. Las reales son las de los commits: el preregistro se fijó en `b0d2bc3` a las 19:35:44 y la nota 1 en `4439218` a las 19:55:32. El orden no cambia: preregistro, después la comprobación con el mock. No se toca nada más.

## Nota 2 · 12 de octubre de 2026, 06:02 (Madrid): ejecución en la nube

- **Dónde:** la ejecución completa se hace en la sesión en la nube, con CLASS v3.4.1 y el parche compilados allí, con 2 núcleos (`--workers 2`), y no en WSL.
- **Comprobación del entorno:** las cinco referencias térmicas salen idénticas bit a bit a las de U5 en WSL (umbral estricto 0.234083460111913).
- **Reinicio:** la máquina se reinició a las 06:00 con el trabajo en las filas, después de la rejilla térmica y de 2 de las 31 filas con CLASS. Se relanzó con el mismo comando. La caché de CLASS conserva todo lo calculado, así que las ejecuciones ya hechas no se repiten y el resultado no cambia. El primer log es `lya/results/u12/u12_class.log` y el de la reanudación, `u12_class_run2.log`.
- No se ha leído ningún δA de fila antes de esta nota: los logs solo informan de que cada ejecución de CLASS terminó.

## Nota 3 · 12 de octubre de 2026, 09:25 (Madrid): resultado de WSL

U12 se ejecutó completo en WSL (`python -B lya/u12_class.py --workers 32`, 1.26 h, CLASS v3.4.1 con el parche y la caché térmica de U5) mientras la ejecución en la nube seguía con las filas. La de la nube se detuvo a mano a las 09:10, con 16 de 31 filas, sin haber calculado ningún δA (`lya/results/u12/nube/LEEME.md`). El resultado válido es el de WSL: `lya/results/u12/U12_RESULT.json`, con el log de la terminal en `u12_class_wsl.log`. Las cinco referencias térmicas coinciden bit a bit en las dos máquinas.
