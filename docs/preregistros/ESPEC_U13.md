# U13 · ¿Sale m_G de la acción común sin parámetros nuevos? · Preregistro

*11 de octubre de 2026, 22:12 (Madrid). Tarea 3 de `docs/PLAN_NUBE_2026-10-11.md`.*

**Declaración previa.** Al leer `docs/modelo/MODELO_ACCION_COMUN.md` para escribir esto, el autor ya razonó de forma cualitativa sobre las cuatro vías de abajo. Su expectativa es que ninguna funcione: todos los términos de la acción declarada son invariantes O(4), y el sector Q no tiene anomalía. No se ha ejecutado ningún código de U13. El criterio y las cuentas se fijan aquí sin depender de esa expectativa: si una vía da el orden de magnitud correcto, se informa así.

## Objetivo

La materia oscura de la ballena necesita una masa del pseudo-Goldstone m_G ≈ 10⁻¹⁵–10⁻¹⁴ eV, con f = F/2 = 5.35×10¹⁴ GeV. Hoy esa masa entra como una ruptura explícita añadida, −εQ₁₁ (C35b), que no está en la acción común (S3: `goldstone_mass_in_action = 0`). La pregunta es si alguna pieza de la acción la genera.

**Escala que hay que producir:** V_ruptura ~ m_G² f², es decir Λ = (m_G f)^(1/2) ≈ 20–70 keV, o V_ruptura/f⁴ ~ 10⁻⁷⁸–10⁻⁷⁶.

## Vías que se examinan (las del plan, más dos de control)

1. **Correcciones radiativas de ξ R S/2**, con S = Tr q².
2. **El acoplamiento σ entre las hojas de χ** (σ f(χ) b²/2 y σ χ |H|²).
3. **Instantones del sector Q.** π₃(RP³) = ℤ da texturas; en 4D euclídeo, las configuraciones de acción finita del modelo sigma se clasifican por π₄(RP³) = π₄(S³) = ℤ₂.
4. **El potencial U(q) y el término cinético**, como control: lo que garantiza la masa cero a nivel árbol.
5. **Gravedad cuántica no perturbativa** (agujeros de gusano y agujeros negros), como control de lo que no es calculable sin hipótesis nuevas. Estimación estándar: V ~ M_Pl² f² e^(−S_wh), con S_wh ≈ c M_Pl/f y c ∈ [0.5, 3].

## Cálculo (`action/u13_mG.py`)

- **Comprobación numérica de simetría:** para Q aleatorias y rotaciones O ∈ O(4) aleatorias (incluidas reflexiones), |U(OQOᵀ) − U(Q)|, |S(OQOᵀ) − S(Q)| y |Tr(∂q∂q) rotado − original| tienen que ser 0 a precisión de máquina. Con la forma exacta del código: λ = 320, λ_P = 192, Z_Q = 400 y X₀ como raíz del gradiente.
- **Teorema que se aplica (sin cálculo):** si todos los términos de la acción que tocan Q son invariantes bajo un grupo G que actúa transitivamente sobre la variedad de vacío, el potencial efectivo sobre esa variedad es constante a todos los órdenes perturbativos, y también con instantones si G no es anómalo. Se comprueba que la acción no contiene fermiones ni campos gauge acoplados a Q (no hay anomalía posible).
- **Vía 5:** m_G(c) = [M_Pl² f² e^(−c M_Pl/f)]^(1/2)/f para c ∈ {0.5, 1, 2, 3}; y el c que haría falta para obtener 10⁻¹⁵ eV.
- **Control informativo** (fuera del criterio «sin parámetros nuevos»): el orden n de un operador c Qⁿ/M_Plⁿ⁻⁴ que respete un ℤ_N discreto y dé m_G ≈ 10⁻¹⁵–10⁻¹⁴ eV con c = 1.

## Veredicto

Para cada vía:
- **SALE:** produce, sin parámetros nuevos, m_G dentro de un factor 10 de [10⁻¹⁵, 10⁻¹⁴] eV;
- **SALE_CON_HIPOTESIS:** solo con una hipótesis o un número nuevo (un coeficiente de orden uno no calculable, un grupo discreto, una acción de agujero de gusano);
- **NO_SALE:** la simetría la prohíbe, o la estimación queda más de 10 órdenes de magnitud fuera.

**Global:**
- **M_G_DE_LA_ACCION:** alguna de las vías 1–3 da SALE;
- **M_G_REQUIERE_RUPTURA_NUEVA:** ninguna de 1–3 sale y la ruptura explícita sigue siendo un ingrediente añadido.

## Lo que U13 no resuelve

- No construye una teoría UV ni una ruptura explícita mejor motivada.
- No evalúa si −εQ₁₁ es la mejor elección entre las rupturas posibles, ni su estabilidad radiativa (ε es técnicamente natural, porque restaura una simetría).
