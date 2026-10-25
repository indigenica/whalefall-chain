# U17c · La rama con ξ creciente no está muerta, pero exige un colapso con dos poblaciones

*14 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17c.md` (commit `51762ff`, anterior a cualquier fila).
- *Código:* `lya/u17c_class.py`, ejecutado en WSL con 16 procesos (6.8 h). *Resultado:* `lya/results/u17c/U17C_RESULT.json`.
- *Lectura posterior, sin veredicto:* `lya/u17c_post.py`, que escribe `lya/results/u17c/U17C_POST.json`.
- CLASS v3.4.1 con el parche y las puertas de U12 sin cambios. Umbral estricto δA = 0.2341, el mismo valor que en U12 y U17b.

## 0. En una frase

**Veredicto preregistrado: RAMA_COMODA_FORZADA.** η_crit(B4a, c = 3, σ_ln = 0.5) = 0.1. Solo pasa Lyα si como mucho el 10 % de la energía que queda en la red al colapsar acaba en cuantos de momento ~m_G.

Las demás combinaciones dan:
- B4a: η_crit = 0.1, también con c = 1 y con σ_ln = 1;
- B1: η_crit ≈ 0.045–0.048 en las cuatro variantes.

La lectura posterior (§4) añade algo que no estaba preregistrado y que pesa más que el veredicto. **Lyα no depende de c: depende de una sola variable, ξℓη.** Por eso la rama cómoda no se salva con un colapso «más caliente». Necesita que el colapso reparta su energía en dos poblaciones muy distintas, y la única simulación de U(1) que conocemos no lo hace.

## 1. Control

Las filas que repiten U17b dan el mismo δA a la tercera cifra:
- B4a, c = 3, η = 1: 0.819;
- B4a, c = 3, η = 0.1: 0.249;
- B1, c = 3, η = 1: 0.907;
- B1, c = 3, η = 0.1: 0.378.

El control B2 (c = 1, η = 1) da δA = 0.246, con forma de 7.1 keV y f = 0.84: pasa la puerta primaria, como B2 con c = 3 en U17b (0.247). Las referencias reproducen las de U12.

## 2. Filas (σ_ln = 0.5; con σ_ln = 1, δA cambia como mucho un 4 %)

| Rama | c | η | K_eff | m_G (eV) | T_osc | δA | forma (m, f) | puerta primaria |
|---|---|---|---|---|---|---|---|---|
| B4a | 3 | 1 | 383 | 5.0×10⁻²⁰ | 8 keV | 0.819 | 1.3 keV, 1.00 | falla |
| B4a | 3 | 0.3 | 109 | 5.8×10⁻¹⁹ | 27 keV | 0.510 | 3.0 keV, 0.89 | falla |
| B4a | 3 | **0.1** | 35 | 5.3×10⁻¹⁸ | 82 keV | 0.249 | 7.5 keV, 0.93 | **pasa** (por forma) |
| B4a | 3 | 0.03 | 10 | 7.3×10⁻¹⁷ | 265 keV | 0.078 | 10 keV, 0.32 | pasa |
| B4a | 3 | 0.01 | 3.1 | 5.6×10⁻¹⁶ | 0.66 MeV | 0.023 | 6.1 keV, 0.01 | pasa |
| B4a | 1 | 1 | 1204 | 5.0×10⁻²¹ | 3 keV | 0.828 | 1.3 keV, 1.00 | falla |
| B4a | 1 | 0.3 | 343 | 6.2×10⁻²⁰ | 9 keV | 0.515 | 3.0 keV, 0.90 | falla |
| B4a | 1 | 0.1 | 109 | 5.8×10⁻¹⁹ | 27 keV | 0.252 | 4.5 keV, 0.45 | pasa (por forma) |
| B4a | 1 | 0.03 | 31 | 6.5×10⁻¹⁸ | 91 keV | 0.094 | 10 keV, 0.40 | pasa |
| B1 | 3 | 0.1 | 65 | 1.6×10⁻¹⁸ | 45 keV | 0.378 | 4.0 keV, 0.80 | falla |
| B1 | 3 | 0.03 | 18 | 2.0×10⁻¹⁷ | 150 keV | 0.144 | 10 keV, 0.73 | pasa |
| B1 | 1 | 0.1 | 204 | 1.7×10⁻¹⁹ | 15 keV | 0.382 | 4.0 keV, 0.81 | falla |
| B1 | 1 | 0.03 | 58 | 2.0×10⁻¹⁸ | 50 keV | 0.157 | 8.1 keV, 0.52 | pasa |
| B2 (control) | 1 | 1 | 105 | 6.2×10⁻¹⁹ | 28 keV | 0.246 | 7.1 keV, 0.84 | pasa (por forma) |

B1 falla con η ≥ 0.1 y pasa con η ≤ 0.03 en las cuatro variantes.

**El η_crit preregistrado** es el mayor η de la rejilla que pasa la puerta primaria. Si hay cambio, se interpola δA en ln η.
- **En B4a,** la fila de η = 0.1 pasa por forma, aunque su δA (0.249) queda un 6 % sobre el umbral estricto. La interpolación da t < 0, se acota y η_crit = 0.1.
- **Si se mirara solo δA** (no es el criterio), el cruce estaría en η ≈ 0.09.
- **Si se interpolara la masa de la forma** (7.5 keV en η = 0.1; 3.0 keV en 0.3) hasta 5.3 keV, saldría η ≈ 0.15.

En los tres casos el veredicto es el mismo: η_crit queda entre 0.03 y 0.3, es decir, **FORZADA**.

**Otros puntos que pide el preregistro:**
- **Régimen antiguo.** Ninguna fila llega a m_G > 10⁻¹⁵ eV con T_osc > 1 MeV: el máximo es 5.6×10⁻¹⁶ eV y 0.66 MeV (B4a, c = 3, η = 0.01). Ninguna vuelve a los minicúmulos de 0.2–10 M_⊙ de U14.
- **Minicúmulos.** En las filas que pasan pesan entre 18 y 5×10⁵ M_⊙; con B4a, c = 3, η = 0.1, ~2×10⁴ M_⊙.
- **Cota de la materia oscura ultraligera.** Solo la cruzan filas que ya fallan Lyα.

## 3. Radiación oscura (lo que el preregistro pedía decir)

Si la fracción 1 − η de la energía de la red es radiación oscura, su cociente con la población del colapso es ρ_DR/ρ_col = [(1 − η)c/η]·(a_osc/a). Evaluado en la igualdad materia–radiación, y tomando g*s(T_osc) = g*s(T_eq) (la elección que da una cota superior), resulta:
- **ΔN_eff ≤ 1.7×10⁻³ en B4a;**
- **≤ 3×10⁻³ en B1.**

Son cifras casi independientes de η y c, porque T_osc ∝ c/η compensa el factor (1 − η)c/η. **No se ve con CMB-S4** (σ ≈ 0.03). El temor del preregistro («con η pequeño puede ser apreciable») no se cumple.

## 4. Lectura posterior: la variable que manda es ξℓη, no c (no preregistrada)

Cuando K_eff ≫ 1 se cumple que Ω_DM ≈ Ω_mis·K_eff ∝ m_G^{1/2}·K_eff. Al ajustar Ω_DM, la masa sale m_G ∝ K_eff⁻². La longitud de libre recorrido de la población del colapso es entonces:

λ_fs ∝ (p/m_G)·a_osc ∝ c·m_G^{−1/2} ∝ c·K_eff = (8π/27)·ξℓ·η.

**c se cancela.** Un colapso más caliente (más c) da menos cuantos por unidad de energía. Eso obliga a subir m_G, que oscila antes y se enfría más, y las dos cosas se compensan.

**Las filas lo confirman:**
- con el mismo η, cambiar c de 1 a 3 mueve δA como mucho 0.016 (B4a, η = 0.03);
- con la misma m_G, cambiar c de 1 a 3 sí lo cambia: con K_eff = 109, B4a da 0.252 con η = 0.1 y c = 1, y 0.510 con η = 0.3 y c = 3;
- las 41 filas caen en una sola curva δA(A·η), con A = (8π/27)ξℓ: A ≈ 105 en B2, 1150–1200 en B4a y 2150–2250 en B1;
- **δA cruza el umbral estricto en A·η ≈ 95, es decir, ξℓη ≈ 100.**

**Consecuencias:**
- **Para la rama viva (B2, η = 1),** el resultado de U17b (equivalente a una WDM de ~7–8 keV) no depende de c mientras K_eff ≫ 1, es decir, mientras c/η ≲ 10. Es una predicción más robusta de lo que decía U17b, que la presentaba como «con c/η ≈ 3».
- **Para las ramas que crecen,** un espectro de colapso con una sola población (todo a p ~ c·m_G, η = 1) no pasa Lyα con ningún c: B4a da δA ≈ 0.81–0.83 con c de 1 a 30 (U17b y U17c). Solo pasan si η ≲ 0.1, y η < 1 significa que el colapso tiene dos poblaciones:
  - una blanda (p ~ m_G), con ≤ 10 % de la energía, que es la materia oscura;
  - otra dura, con ≥ 90 %, que tiene que comportarse como radiación.

**Cuánto de dura tiene que ser.** Los cuantos duros se vuelven no relativistas en a/a_osc ≈ p_h/m_G y entonces cuentan como materia, muy caliente. Su fracción de la materia oscura es f_h ≈ [(1 − η)c/η]·(m_G/p_h). Para que f_h < 1 %, con B4a y η = 0.1, hace falta **p_h ≳ 900 m_G (c = 1) o p_h ≳ 2600 m_G (c = 3)**. Con η = 0.03, p_h ≳ 3000–9000 m_G.

## 5. Qué significa

- **Veredicto: FORZADA, no muerta.** La rama cómoda (ξ creciente, q < 0.91) sigue siendo posible, pero ya no basta con decir «η ≲ 0.1». Hace falta:
  1. que ≥ 90 % de la energía de la red en H ≈ m_G/3 acabe en cuantos con p ≳ 10³ m_G;
  2. y que < 10 % acabe en cuantos de p ~ m_G.
- **Comparación con U(1).** En la única simulación de U(1) que tenemos (Hiramatsu et al. 2012, arXiv:1202.5851; nota 1 de INFORME_U17), la energía media de todos los axiones radiados es ω̄/m_a ≈ 3.3, y la mayoría tienen k ~ m_a.
  - Si el 90 % de la energía fuera a cuantos con p ≳ 10³ m_a, esa media saldría ≳ 30.
  - Dentro del rango que su red resuelve, **la η efectiva de U(1) es ≈ 1.**
  - Su red no resolvería una componente con p ~ 10³ m_a, así que no la excluye del todo. Pero tampoco hay nada que la sugiera: las paredes tienen grosor 1/m_G y radian a esa escala.
- **Lectura honesta.** La rama con ξ creciente solo vive si el colapso de RP³ es cualitativamente distinto del de U(1). No se puede descartar sin medirlo, pero no es la apuesta natural. **El sector oscuro de la ballena depende, en la práctica, de que ξ se estabilice**, como ya decía U17b.
- **Lo bueno:**
  - la rama viva sale más robusta: su Lyα no depende de c;
  - la radiación oscura del colapso es despreciable;
  - la pregunta pendiente para la red queda reducida a un número medible: **qué fracción de la energía del colapso de RP³ acaba en modos con k ≲ 10 m_G.**

## 6. Lo que U17c no resuelve

- **η y c de RP³ siguen sin medirse.** El §5 propone medir la fracción de energía del colapso con k ≲ 10 m_G, que es más robusta que c y η por separado. Haría falta U19b: m_G/H ≳ 30, ~1000³, en WSL.
- **Las dos poblaciones del §4 son un razonamiento,** no un cálculo con CLASS: no se han pasado filas con una componente dura explícita.
- **Solo q = 0.87, x₀ = 10, R_Γ = 1.56.** U17b mostró que, con el colapso dominando, q y R_Γ pesan poco (f_em ≤ 10⁻³ en todas las filas).
- **La curva δA(ξℓη) del §4 es una lectura posterior.** No cambia el veredicto preregistrado.

## Nota 1 · 14 de octubre de 2026 (Madrid): auditoría de Codex

- **η_crit mezcla puertas.** Es cierto: η_crit se elige con la puerta primaria (forma), pero se interpola el umbral de δA. En B4a, la fila con η = 0.1 pasa por forma con δA por encima del umbral estricto, así que la interpolación se acota y devuelve el punto de la rejilla. **η = 0.1 es el mayor punto permitido de la rejilla, no una frontera interpolada de la puerta primaria;** esa frontera está en [0.1, 0.3). El §2 ya daba el intervalo: 0.09 por δA y ≈ 0.15 por la masa de la forma. El veredicto, FORZADA, vale con cualquier valor del intervalo.
- **Energía y momento.** Las filas usaron la construcción antigua de U17b. Con el número normalizado por la energía media, la variable del §4 pasa a ser X = (8π/27)·R_ρ·ξℓη·(c/c_E): el factor c/c_E vale 0.94 con c = 3 y 0.69 con c = 1. **c ya no se cancela del todo:** sí con c ≳ 3, pero no con c ≈ 1, donde la población queda más fría que lo calculado. Las conclusiones del §5 no cambian de signo.
- **R_ρ.** Con la energía de RP³ escalada como la emisión (R_ρ ≈ 1.56), X crece un 56 % en todas las filas. Para B4a, η_crit bajaría a ≈ 0.06. U17d lo comprueba en las filas decisivas.
- **Las masas de minicúmulos** de este informe son escalas candidatas, no poblaciones (ver la Nota 2 de INFORME_U17b).
