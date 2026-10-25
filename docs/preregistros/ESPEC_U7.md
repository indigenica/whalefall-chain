# U7 · Calibración de la emisión con las redes más grandes publicadas · Preregistro

*9 de octubre de 2026, ≈12:40 (Madrid). Escrito antes de calcular nada de U7.*

## Motivo

La emisión de S4/S6 sigue la fórmula de GHV: ρ_s = ξ μ/t², con μ = π f_a² ℓ, multiplicada por la eficiencia ε. La rama B2 fija ξ* = 1.2.

Correia, Hindmarsh et al. (arXiv:2512.13653) dan, para la red U(1) que escala, dos números:
- la densidad de número de axiones, **ν_ax = n_ax/(f_a² H) = 1.66 ± 0.17**;
- en su trabajo anterior, la densidad de cuerdas (sistema en reposo), **ζ_* ≈ 1.49**.

U7 calibra nuestra emisión con ellos, antes de dar por buena la supervivencia de B2.

## U7a · ν de nuestra emisión

- **Cálculo.** Con la emisión de S6 para la U(1) en B2 (ε = 1.2/(0.24 ℓ), R = 1), con espectro log-plano y sin masa, se suma el número comóvil emitido hasta T_eval. Luego ν_S4 = (N/a³)/(f_a² H).
- **Rejilla:** x₀ ∈ {10, 30} y T_eval tal que H(T_eval) ∈ {10⁻³⁰, 10⁻³³, 10⁻³⁶} GeV. Así se comprueba que ν_S4 es constante, es decir, que la emisión escala.
- **Factor de calibración:** κ(x₀) = 1.66/ν_S4(x₀), con la media sobre T_eval. Se informa su incertidumbre, que viene del 10 % de ν_ax.

## U7b · filas calibradas por CLASS

**Familias de filas:**
- **B2κ:** ε = κ(x₀) · R · 1.2/(0.24 ℓ), con espectro log-plano. Es la emisión U(1) calibrada con ν_ax y escalada por R_Γ de RP³.
- **B2ζ:** ε = R · 1.49/(0.24 ℓ), con espectro log-plano. Solo cambia la normalización de ξ.

**Cada familia:**
- R ∈ {−1σ, central, +1σ} con R_Γ = 1.40 ± 0.13;
- x₀ ∈ {10, 30} y K_sw ∈ {0, 1};
- 24 filas en total.

Cada fila pasa por S6 (raíz de abundancia, P1–P4) y luego por la puerta Lyα de U4a/U5 sin cambios: CLASS con el parche y 45 nodos, umbral estricto 0.234, forma estricta como puerta primaria y malla ampliada de U5b.

**Veredicto U7** (sobre B2κ):
- **CALIBRADA_PERMITIDA:** alguna fila con R a +1σ permitida;
- **CALIBRADA_EXCLUIDA:** ninguna fila con R a −1σ permitida;
- **MARGINAL:** en otro caso.

Se informa también B2ζ, y κ en función de x₀.

**Advertencia.** ν_ax es un número para U(1) y para el espectro de número que miden. La traslación a RP³ usa R_Γ, que es un cociente de energías, no de números. U6 permitirá comprobar si los espectros de RP³ y de U(1) tienen la misma forma, que es lo que hace válida esa traslación.
