# U13 · ¿Sale m_G de la acción común?

*11 de octubre de 2026, noche.*
- *Preregistro:* `docs/preregistros/ESPEC_U13.md` (commit `0fbb0b2`, 22:11, antes de ejecutar código).
- *Código y resultado:* `action/u13_mG.py` → `action/results_u13/U13_RESULT.json`. Se ejecuta en dos segundos.

## 0. En una frase

**Veredicto preregistrado: M_G_REQUIERE_RUPTURA_NUEVA.** Todos los términos de la acción común que contienen Q son invariantes bajo O(4), y el sector Q no tiene anomalía. Por eso la masa del pseudo-Goldstone es exactamente cero a todos los órdenes, con instantones incluidos. La materia oscura de la ballena necesita una ruptura explícita añadida, como la −εQ₁₁ actual; no la genera la acción.

## 1. Lo que se ha comprobado

| Vía | Resultado | Por qué |
|---|---|---|
| 1. Correcciones radiativas de ξ R S/2 | NO_SALE | S = Tr q² es invariante O(4) (desviación ≤ 2×10⁻¹⁵ con 2000 rotaciones aleatorias, reflexiones incluidas). Los lazos de una teoría simétrica solo generan operadores simétricos. |
| 2. Acoplamiento σ entre las hojas de χ | NO_SALE | σ acopla χ con el sector visible (σ f(χ) b², σ χ \|H\|²). No contiene Q, así que no puede romper su simetría. |
| 3. Instantones del sector Q | NO_SALE | Sin fermiones ni campos gauge acoplados a Q, O(4) no es anómala. Las configuraciones euclídeas (π₄(RP³) = ℤ₂) y las texturas (π₃ = ℤ) respetan la simetría: al sumar sobre sus orientaciones, el potencial efectivo sobre RP³ sigue constante. |
| 4. U(q) y término cinético (control) | NO_SALE | El vacío es un proyector de rango 1 (autovalores 0, 0, 0, 1). El hessiano tiene exactamente 3 modos nulos: los Goldstones, sin masa a nivel árbol. |
| 5. Gravedad cuántica (control) | NO_SALE | Con V ~ M_Pl² f² e^(−c M_Pl/f) y M_Pl/f = 4551, m_G ~ 10⁻⁴⁶⁷ eV para c = 0.5 y menos aún para c mayor. Para llegar a 10⁻¹⁵ eV haría falta c ≈ 0.04, muy por debajo del rango plausible [0.5, 3]. |

**Escala que habría que producir:** Λ = (m_G f)^(1/2) ≈ 23–73 keV, es decir V/f⁴ ≈ 10⁻⁷⁸–10⁻⁷⁶.

## 2. Una pista informativa (fuera del criterio)

Una simetría discreta ℤ_N gauge que proteja a Q hasta operadores de orden n daría m_G² ~ f² (f/M_Pl)^(n−4). Con coeficiente 1, **n ≈ 25** da justo m_G ≈ 10⁻¹⁵–10⁻¹⁴ eV. Es el mecanismo estándar de «calidad» de los axiones. No cumple el criterio de U13, porque N (un entero) y la simetría gauge son ingredientes nuevos. Pero convertiría el número arbitrario ε en un entero, y eso sería una mejora.

## 3. Consecuencia

- La afirmación «m_G es ruptura explícita añadida (no está en la acción común)», del contrato y de S3, queda **confirmada y explicada**: no es una carencia del cálculo, sino una consecuencia de la simetría.
- El artículo debe presentarla como un parámetro de la ballena, no como algo pendiente de derivar dentro de esta acción.
- Lo que queda abierto es la teoría UV de la ruptura (por ejemplo, ℤ_N con N ≈ 25), no un cálculo dentro de la acción actual.
