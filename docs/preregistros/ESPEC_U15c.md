# U15c · La rama viva (ξ constante con el colapso cuerdas–paredes) frente a Planck + BAO + Lyman-ABGD · Preregistro

*Escrito antes de calcular ninguna tabla ni cadena de U15c (la hora es la de su commit). Preparado en la nube; se ejecuta en WSL con la tubería de U15b.*

## Motivo

U17b (`docs/informes/INFORME_U17b.md`): con el colapso de la red contado, solo vive ξ constante (B2). Con c/η ≈ 3 (la literatura de U(1)), el 97 % de la materia oscura son cuantos del colapso. En Lyα equivalen a una WDM de ~8 keV: pasan la puerta de forma, pero quedan un 6 % por encima del umbral estricto de δA. Es el caso que hay que pasar por la verosimilitud completa. U15b prueba la cola tibia del escenario antiguo (K_sw ∈ {0, 1}), que ya no es el caso físico.

## Casos

Rama B2, q = 0.87, x₀ = 10, R_Γ = 1.56, con la población del colapso de U17b (log-normal, σ_ln = 0.5):
- **`U15c_B2_c3_eta1`:** c/η = 3 (literatura), f_nc = 0.97;
- **`U15c_B2_c30_eta1`:** c/η = 30, f_nc = 0.77.

La validación es la de U15b (mezclas térmicas de 1 y 2 keV): no se repite; se usa su resultado.

## Método (el de ESPEC_U15b, salvo lo que se dice)

- **Tablas** T²(k; s) = P/P_ΛCDM a z = 0 con CLASS y el parche:
  - ω_ncdm = s·f_nc·Ω_DM y ω_cdm = Ω_DM − ω_ncdm;
  - s ∈ {0, 0.25, 0.5, 0.75, 0.9, 1, s_max}, con s_max = min(1.25, 1/f_nc) (1.03 en el caso c/η = 3);
  - `python -B lya/u15b/prepare_tables.py --u15c --workers N`.
- **Cadenas** con `Lya_abgd_whale` sin cambios, prior plano de s en [0, s_max] y 4 cadenas por caso:
  - `CASES="U15c_B2_c3_eta1 U15c_B2_c30_eta1" bash lya/u15b/run_wsl.sh chains`.
- **Convergencia:** R − 1 < 0.03 en todos los parámetros, como en U15b.

## Veredicto (`lya/u15c/analyze_u15c.py`)

Como s_max es casi 1, el criterio de U15b (s_95 ≥ 1) no sirve: una posterior plana daría s_95 ≈ 0.98 y excluiría el modelo sin que los datos digan nada. La región HPD del 95 % tampoco: con una posterior casi plana, qué bin queda fuera lo decide el ruido (se comprobó con cadenas sintéticas antes de fijar este texto). Se usa el cociente de densidades posteriores, con el histograma de 30 bins de U15b:
- **PERMITIDA:** p(s = 1)/p_max ≥ e⁻² ≈ 0.135. Para una posterior gaussiana equivale a Δχ² ≤ 4, es decir, a 2σ;
- **EXCLUIDA:** p(s = 1)/p_max < e⁻².

**Global:**
- **RAMA_VIVA_PERMITIDA_LYA_COMPLETO:** los dos casos PERMITIDOS;
- **RAMA_VIVA_EXCLUIDA_LYA_COMPLETO:** los dos EXCLUIDOS;
- **RAMA_VIVA_DEPENDE_DE_C_ETA:** en otro caso.

**Se informa también:** la posterior de s, p(1)/p_max y Δχ²_min entre los puntos con s ≈ 1 (|s − 1| < 0.05) y el mínimo global.

Si la validación de U15b no pasa, U15c es **INCONCLUSO (validación)**.

## Lo que U15c no resuelve

- **c/η** sigue sin medirse para RP³.
- **T² se calcula a la cosmología fiducial,** como en U15b. La auditoría (`audit/u12_u15b`) midió el error de esa aproximación: |ΔlnL| ≈ 0.4–0.7 a 4–5σ de Planck.
- **Con f_nc ≈ 1 la materia oscura es casi toda no fría.** La rejilla de la verosimilitud (grids ABG/ABGD/térmica) la describe como cualquier otra T(k), pero el modelo queda cerca del borde de lo que las simulaciones de Hooper et al. cubren. Se dirá si los puntos de la cadena caen fuera de la rejilla (sus comprobaciones internas).

## Enmienda 1 · 14 de octubre de 2026, ≈21:10 (Madrid): en espera de U17d, y el analizador corregido (antes de generar ninguna tabla ni cadena)

**1. En espera.** La auditoría de Codex (`audit/codex_20261014`) encontró dos fallos en la población del colapso de U17b, sobre la que se definieron los dos casos:
- mezclaba energía y momento;
- no contaba la energía de la red de RP³, R_ρ.

U17d (`docs/preregistros/ESPEC_U17d.md`) rehace las filas decisivas. **Los casos de U15c se redefinirán en una enmienda 2,** con el resultado de U17d y antes de cualquier tabla o cadena. Mientras tanto, `prepare_tables.py --u15c` se niega a ejecutarse.

**2. El analizador** (`lya/u15c/analyze_u15c.py`). Lo señaló Codex, y los dos casos se reprodujeron con cadenas sintéticas.
- **Los bins.** El histograma de s se construye sobre el rango del prior, leído de `log.param`, y no sobre el de los datos. Si s = 1 queda fuera del prior, el caso da error. Antes, el índice de s = 1 se recortaba al borde: un histograma uniforme en [0, 0.5] daba p(1)/p_max = 1 y PERMITIDA. Ahora da EXCLUIDA, con p(1) = 0.
- **Las etiquetas.** PERMITIDA y EXCLUIDA no se imprimen ni se guardan hasta que el caso converja (regla de la Nota 6 de U15b) y la validación de U15b pase.
- **La estabilidad.** Si la etiqueta cambia con 20, 30 o 45 bins, o al quitar una de las cadenas, el caso es **INDECISO** y el global, **INCONCLUSO (p(1)/p_max en la frontera)**.
- **El criterio no cambia:** p(1)/p_max ≥ e⁻² con 30 bins.

## Enmienda 2 · 15 de octubre de 2026 (Madrid; la hora es la de su commit): los casos, redefinidos con U17d (antes de cualquier tabla o cadena)

U17d (`docs/informes/INFORME_U17d.md`) rehízo la población del colapso con el balance corregido y la energía de RP³ (R_ρ = 1.56, r = 1). Los dos casos de U15c pasan a ser:
- **`U15c_B2_eta1_R1.56`:** colapso tipo U(1), c = 3 y η = 1. f_nc = 0.981, m_G = 2.6×10⁻¹⁸ eV, s_max = 1.020. Es la fila que en U17d pasa solo por forma, con δA = 0.323.
- **`U15c_B2_eta0.3_R1.56`:** colapso «duro», c = 3 y η = 0.3. f_nc = 0.937, m_G = 3.2×10⁻¹⁷ eV, s_max = 1.068. Es la que pasa todas las variantes de U17d.

Los dos usan σ_ln = 0.5 y el número normalizado por la energía media (`lya/u17b_class.row`).

Ni el método, ni el criterio (p(1)/p_max ≥ e⁻², con la enmienda 1) ni la validación cambian. El global se renombra:
- RAMA_VIVA_PERMITIDA_LYA_COMPLETO: los dos casos PERMITIDOS;
- RAMA_VIVA_EXCLUIDA_LYA_COMPLETO: los dos EXCLUIDOS;
- **RAMA_VIVA_DEPENDE_DE_ETA:** en otro caso.

Un PERMITIDO prueba el sector oscuro con la cosmología marginalizada, n_s incluido, no la ballena completa (que predice n_s ≈ 0.973): se dirá en el informe. Se retira el freno de `prepare_tables.py --u15c`.

## Enmienda 3 · 16 de octubre de 2026 (Madrid; la hora es la de su commit): los casos, en los bordes de la banda de F (antes de cualquier tabla o cadena)

U21 y U20 (`docs/informes/INFORME_U21.md`, `INFORME_U20.md`) mostraron dos cosas:
- **F no está fijado.** El benchmark F₀ queda fuera por la isocurvatura.
- **La banda de la rama viva es F ≤ 0.5·F₀,** o ≤ 0.35·F₀ si el colapso es caliente.

Los casos de la enmienda 2 estaban en F₀, que ya está excluido. **Los dos casos pasan a ser los bordes de la banda,** donde la verosimilitud completa puede discrepar de las puertas de cribado. Ambos usan la población del colapso de U17d: c = 3, η = 1, R_ρ = 1.56, σ_ln = 0.5.

- **`U15c_B2_F0.5_r1`:** F = 0.5·F₀ y r = 1. Es el borde primario: m_G = 5.5×10⁻¹⁷ eV, f_nc = 0.980, δA = 0.094 en U21.
- **`U15c_B2_F0.35_r3`:** F = 0.35·F₀ y r = 3. Es el borde con colapso caliente: m_G = 2.6×10⁻¹⁶ eV, f_nc = 0.979, δA = 0.121 en U21.

El método, el criterio (p(1)/p_max ≥ e⁻², con la enmienda 1) y la validación no cambian. El global: los dos casos PERMITIDOS dan **RAMA_VIVA_PERMITIDA_LYA_COMPLETO**; los dos EXCLUIDOS, **RAMA_VIVA_EXCLUIDA_LYA_COMPLETO**; uno de cada, **RAMA_VIVA_DEPENDE_DE_ETA**, que aquí significa que depende de la época del colapso.

Se ejecuta en WSL cuando haya núcleos (8 cadenas):

```bash
python -B lya/u15b/prepare_tables.py --u15c --workers 8
CASES="U15c_B2_F0.5_r1 U15c_B2_F0.35_r3" bash lya/u15b/run_wsl.sh chains
```

## Enmienda 4 · 17 de octubre de 2026, ≈19:06 (Madrid): el mismo cambio de muestreador que ESPEC_U15b, Nota 8 (antes de leer ninguna posterior)

Las cadenas de U15c llevan corriendo desde el 16 de octubre. No se ha mirado nada de ellas, ni siquiera su convergencia. Se les aplica `boost`, como a la validación de U15b (ESPEC_U15b, Nota 8):
- la covarianza de la propuesta se rehace con las cadenas;
- el bloque de Lyα se sobremuestrea ×10;
- se usan 2 hilos de CLASS por cadena;
- cada cadena continúa desde su último punto.

La verosimilitud, los casos, las puertas y el veredicto no cambian.

## Nota · 20 de octubre de 2026, 09:57 (Madrid): resultado

- **El análisis:** el autor ejecutó `analyze_u15c.py` en WSL el 20 de octubre, a las 06:48, después de que U15b validara.
- **Los dos casos convergen:** R − 1 ≤ 0.017, y ≤ 0.025 por mitades.
- **Los dos casos son PERMITIDOS:**
  - `U15c_B2_F0.5_r1`: p(1)/p_max = 0.874;
  - `U15c_B2_F0.35_r3`: p(1)/p_max = 0.913.
- **La etiqueta es estable** con 20, 30 y 45 bins y al quitar cualquiera de las cadenas.
- **Veredicto: RAMA_VIVA_PERMITIDA_LYA_COMPLETO.**
- **Informe:** `docs/informes/INFORME_U15b_U15c.md`. Explica que la posterior de s es plana porque la supresión del modelo en k = 10 h/Mpc es < 0.15 %.
