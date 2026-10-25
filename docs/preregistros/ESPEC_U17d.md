# U17d · La rama viva con el balance del colapso corregido y el puente RP³ declarado · Preregistro

*14 de octubre de 2026, ≈21:05 (Madrid). Escrito antes de calcular ninguna fila de U17d (la hora que cuenta es la de su commit). Responde a los puntos 1, 2 y 6 de la auditoría de Codex (`audit/codex_20261014/`). Se ejecuta en la nube, que conserva la caché de CLASS de U12/U17b.*

## Motivo

La auditoría de Codex encontró dos fallos en la población del colapso de U17/U17b/U17c, y los dos se confirman (`audit/codex_20261014/RESPUESTA.md`).

**1. Energía y momento mezclados.**
- U17 cuenta los cuantos dividiendo la energía entre c·m_G, es decir, toma c como energía media por cuanto.
- U17b construye una log-normal cuyo **momento** medio es c·m_G y conserva aquel número.
- Por eso la energía de los cuantos supera la asignada en un 6.6 % (c = 3, σ_ln = 0.5) y hasta un 55 % (c = 1).
- Corrección (código, `lya/u17b_class.row`): el número es η·ρ_red/⟨E⟩, con ⟨E⟩ = c_E·m_G la energía media de la log-normal de momento medio c·m_G. Con c = 3 y σ_ln = 0.5, c_E = 3.199. Coincide con Hiramatsu et al.: momento 3.12 y energía 3.28.

**2. ¿Qué energía tiene la red de RP³?** El inventario de U17 usa ξ de U(1):
- B2: 1.20, de la literatura;
- B4a: la extrapolación de ξ_U1(5.65) = 0.616.

La emisión de U12, en cambio, multiplica la de U(1) por R_Γ = 1.56. Si la emisión de RP³ se escala con R_Γ, la energía de la red de RP³ debe escalarse igual.

En escalado y en radiación, la energía que pierde la red es Γ ≈ (ρ_red/t)·[1 − O(d ln ξμ/d ln t)]. El cociente de emisiones mide entonces el cociente de energías hasta correcciones O(1/ℓ): **R_ρ ≈ R_Γ**. Además, R_ξ = 1.61 en la misma red, lo que indica que las tensiones de las dos redes coinciden a un 3 %.

## Las hipótesis del puente (se declaran, no se miden)

- **R_ρ = ρ_red,RP³/ρ_red,U(1)** en el colapso:
  - **1.56 (primaria):** la misma extrapolación que ya usa la emisión;
  - **1 (alternativa):** las dos densidades convergen a ℓ ≈ 70. Si fuera así, la emisión también debería converger, pero con q = 0.87 la emisión pesa ≤ 10⁻³ de la materia oscura y no se toca.
- **r = a_decay/a_osc,** la época en la que los cuantos del colapso tienen momento c·m_G:
  - **1 (primaria):** el colapso emite en H ≈ m_G/3;
  - **2 y 4 (sensibilidad):** la receta de Hiramatsu et al. La energía se congela en H = m_G/3 y se diluye como materia, y el momento c·m_G se alcanza en a_decay (H ≈ m_G/12 y m_G/48).
  - *Por qué se toman como pesimistas.* Si la red sigue en escalado hasta a_decay, la energía que pierde entre a_osc y a_decay va a cuantos de k ~ H < m_G, que son materia fría. La población caliente baja como 1/r y su longitud de libre recorrido no cambia (c·K_eff·r es constante). Con masa constante, las paredes dominan sobre las cuerdas hacia t ≈ πℓ/(8m_G), es decir, H ≈ m_G/50. r = 4 cubre ese orden.

## Filas (`lya/u17d_class.py`)

Todas con c = 3, σ_ln = 0.5, q = 0.87, x₀ = 10, R_Γ = 1.56 y el balance corregido:
- **B2,** η ∈ {1, 0.3, 0.1}:
  - con R_ρ ∈ {1.56, 1} y r = 1;
  - con R_ρ = 1.56 y r ∈ {2, 4};
- **B4a:** η ∈ {0.1, 0.03} con R_ρ = 1.56, y η = 0.1 con R_ρ = 1;
- **Control:** la fila B2, c = 3, η = 1 de U17b tal como se ejecutó (`legacy_energy`). Debe dar δA = 0.247 ± 0.002.
- **Total: 16 filas.** Puertas de U12 sin cambios: umbral estricto 0.2341 y puerta primaria de forma si es representable.

## Veredicto (B2 con las hipótesis primarias: R_ρ = 1.56, r = 1)

- **RAMA_VIVA_CON_COLAPSO_TIPO_U1:** pasa η = 1. Un colapso como el de U(1) basta.
- **RAMA_VIVA_REQUIERE_COLAPSO_DURO:** η = 1 falla, pero pasa η = 0.3 o 0.1.
- **RAMA_VIVA_EN_PELIGRO:** η = 0.1 también falla.
- **INCONCLUSO:** el control falla o CLASS falla en alguna fila.

**Se informa también, sin veredicto:**
- η_crit de B2 con R_ρ = 1 y con r = 2 y 4;
- si B4a pasa con el balance corregido;
- m_G, T_osc y f_nc de cada fila;
- la variable de lectura X = K_eff·c·r. Con las correcciones es (8π/27)·R_ρ·ξℓη·(c/c_E)·r.

## Lo que se espera, escrito antes de calcular

En U17c, δA cruza el umbral estricto en X ≈ 95 (`lya/results/u17c/U17C_POST.json`).
- **Con las hipótesis primarias, B2 con η = 1 tiene X ≈ 152.** Por la curva de U17c, δA ≈ 0.31, y se espera que falle; η = 0.3 (X ≈ 44) debería pasar. Es decir, **se espera RAMA_VIVA_REQUIERE_COLAPSO_DURO.**
- **Con R_ρ = 1, η = 1 da X ≈ 96, justo en el umbral.** La fila U17b (X ≈ 103) pasó por forma.

Si el resultado contradice estas expectativas, se dice.

## Lo que U17d no resuelve

- **R_ρ y r** no se miden: son hipótesis declaradas, con el valor primario elegido por coherencia con la emisión.
- **La isocurvatura de la red con las masas nuevas** (punto 6 de Codex): la comprobación de U5c no vale aquí y queda pendiente. Hay que transferir el ruido blanco de la población del colapso y del desalineamiento.
- **Los minicúmulos** son escalas candidatas: la longitud de libre recorrido de la población del colapso supera su radio comóvil. No se trasladan U14 ni la fórmula fría.
- **U15c queda en espera de U17d** (ESPEC_U15c, enmienda 1). Sus casos se redefinirán con este resultado, antes de cualquier cadena.
