# U7c · Lyα con la R_Γ de las redes densas

*9 de octubre de 2026, noche.*
- *Preregistro:* `docs/preregistros/ESPEC_U7c.md` (≈22:00, antes de calcular).
- *Resultado:* `lya/results/u7/U7C_RESULT.json`. Ejecución en WSL: 1.1 h, CLASS v3.4.1 con el parche.

## 0. En una frase

**Veredicto preregistrado: CALIBRADA_SENSIBLE_A_R.** Con la emisión calibrada con Correia et al., el sector oscuro sigue pasando Lyα con R_Γ = 1.69, el valor de las redes densas de U9, en casi todas las filas. Con R_Γ = 2.09, el máximo de U9 al mover la ventana, falla una esquina: la de x₀ = 30.

## 1. Filas permitidas

Puerta primaria de forma; entre paréntesis, solo con δA ≤ 0.234:

| R_Γ | B2κ (calibrada con ν_ax) | B2ζ (ξ* = 1.49) |
|---|---|---|
| 1.39 (U6, ξ_i = 0.5) | 4/4 (4/4) | 4/4 (4/4) |
| 1.69 (U9, ξ_i = 2) | 3/4 (4/4) | 4/4 (4/4) |
| 2.09 (máximo del barrido de U9) | 2/4 (3/4) | 1/4 (3/4) |

En total, **18/24 con la puerta de forma y 22/24 solo con δA**.

**Dónde falla, en B2κ:**
- **x₀ = 10:** pasa con todas las R_Γ, con δA = 0.09–0.13, lejos del umbral.
- **x₀ = 30:**
  - con R_Γ = 2.09: δA = 0.253 (K_sw = 0), por encima del umbral, y 0.227 (K_sw = 1), que pasa δA pero no la forma;
  - con R_Γ = 1.69: la fila con K_sw = 1 (δA = 0.194) falla la forma.

**Un aviso sobre la puerta de forma en este régimen.** Los ajustes caen en mezclas de 0.2–0.34 keV con un 1–4 % de componente caliente. Ahí la cota f ≤ 0.14 m^1.1 está extrapolada muy por debajo de 1 keV, y en dos filas el ajuste cae en el borde inferior de la malla (0.2 keV). Por eso la puerta no es monótona: la fila con R_Γ = 1.69, x₀ = 30 y K_sw = 0 (δA = 0.219) pasa, y su vecina con K_sw = 1 (δA = 0.194) falla. El veredicto preregistrado usa la forma, pero en este régimen el criterio δA es más fiable. Hay que leer las dos columnas.

## 2. La región que sobrevive (18 filas)

- **Materia oscura pseudo-Goldstone:** m_G ≈ 0.8–4.9×10⁻¹⁵ eV.
- **Fracción emitida por la red:** f_em = 17–59 %.
- **Supresión en Lyα:** δA = 0.09–0.22.
- **Radiación oscura:** ΔN_eff ≤ 0.0075.
- **Minicúmulos:** 0.7–10 M_⊙.

Prácticamente la misma ventana que en U7, pero ahora **hasta un 6 % del umbral**.

## 3. Lectura

- **La sistemática de protocolo de U9 no mata al sector oscuro, pero lo lleva al borde.** Con la R_Γ de redes poco densas pasa con holgura. Con la de redes densas pasa salvo una fila. Con el extremo del barrido de ventana falla la mitad de las filas de x₀ = 30.
- **Lo que decide ahora es una combinación:** cuánto radia de verdad la red RP³ en escalado (R_Γ sin la relajación de la condición inicial, U10a), dónde empieza el espectro de emisión (x₀) y, como siempre, ξ y q a ℓ ≈ 70.
- **La predicción observable se mantiene:** una supresión de pequeña escala en el bosque de Lyα, δA ≈ 0.1–0.2, de un 1–4 % de componente muy caliente. Se ha vuelto más **arriesgada**: está más cerca del límite actual. Es una buena noticia para la falsabilidad. Unos datos de Lyα de alta resolución, o simulaciones hidrodinámicas para esta distribución, podrían cerrar esta rama.
- **La q baja de las redes densas (0.74) iría en sentido contrario:** emisión casi irrelevante y materia oscura fría. U7c no la explora; fija q = 1, como estaba preregistrado.
