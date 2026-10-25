# U17d · Con el balance corregido y la energía de RP³, la rama viva pasa la puerta preregistrada, pero solo por forma y sin margen

*14–15 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17d.md` (commit `bfcbaf8`, anterior a cualquier fila). Responde a la auditoría de Codex (`audit/codex_20261014/`).
- *Código:* `lya/u17d_class.py`, en la nube con 4 procesos (3.0 h). *Resultado:* `lya/results/u17d/U17D_RESULT.json`.
- CLASS v3.4.1 con el parche y las puertas de U12 sin cambios:
  - umbral estricto δA = 0.2341 (mezcla de 1 keV); umbral amplio 0.336 (WDM térmica de 5.3 keV);
  - puerta primaria = forma estricta: la fracción f del ajuste térmico (m, f) no supera f_max(m).

## 0. En una frase

**Veredicto preregistrado: RAMA_VIVA_CON_COLAPSO_TIPO_U1.** Con el balance corregido y la energía de RP³ escalada como la emisión (R_ρ = 1.56), B2 con un colapso como el de U(1) (c = 3, η = 1) **pasa la puerta primaria**.

Pero:
- pasa solo por forma: δA = 0.323 supera el umbral estricto en un 38 %;
- la expectativa preregistrada era que fallara, y no acertó (§3);
- si el colapso termina más tarde que H ≈ m_G/3 (r ≥ 2), η = 1 falla y hace falta un colapso «duro» (η ≲ 0.3);
- su masa, m_G ≈ 2.6×10⁻¹⁸ eV, cae donde la literatura reciente pone cotas de isocurvatura (§4).

## 1. Control

La fila de U17b tal como se ejecutó (B2, c = 3, η = 1, `legacy_energy`) da δA = 0.247 y forma (7.78 keV, 0.97). Son los valores de U17b: el control pasa.

## 2. Filas

| Rama | η | R_ρ | r | X | m_G (eV) | T_osc | δA | forma (m, f) | f_max(m) | puerta primaria |
|---|---|---|---|---|---|---|---|---|---|---|
| B2 | **1** | **1.56** | **1** | 152 | 2.6×10⁻¹⁸ | 57 keV | **0.323** | 4.29 keV, 0.64 | 0.72 | **pasa** (por forma) |
| B2 | 0.3 | 1.56 | 1 | 44 | 3.2×10⁻¹⁷ | 186 keV | 0.116 | 10 keV, 0.54 | 1 | pasa |
| B2 | 0.1 | 1.56 | 1 | 14 | 2.8×10⁻¹⁶ | 479 keV | 0.038 | 3.8 keV, 0.02 | 0.61 | pasa |
| B2 | 1 | 1 | 1 | 96 | 6.1×10⁻¹⁸ | 88 keV | 0.236 | 2.0 keV, 0.12 | 0.30 | pasa |
| B2 | 0.3 | 1 | 1 | 28 | 8.0×10⁻¹⁷ | 275 keV | 0.075 | 10 keV, 0.27 | 1 | pasa |
| B2 | 0.1 | 1 | 1 | 9 | 5.8×10⁻¹⁶ | 674 keV | 0.023 | 6.5 keV, 0.01 | 1 | pasa |
| B2 | 1 | 1.56 | 2 | 304 | 2.6×10⁻¹⁸ | 57 keV | 0.496 | 3.0 keV, 0.86 | 0.47 | **falla** |
| B2 | 0.3 | 1.56 | 2 | 89 | 3.2×10⁻¹⁷ | 186 keV | 0.204 | 9.1 keV, 0.93 | 1 | pasa |
| B2 | 0.1 | 1.56 | 2 | 29 | 2.8×10⁻¹⁶ | 479 keV | 0.075 | 10 keV, 0.31 | 1 | pasa |
| B2 | 1 | 1.56 | 4 | 607 | 2.6×10⁻¹⁸ | 57 keV | 0.687 | 1.9 keV, 0.93 | 0.28 | **falla** |
| B2 | 0.3 | 1.56 | 4 | 177 | 3.2×10⁻¹⁷ | 186 keV | 0.335 | 4.24 keV, 0.69 | 0.71 | pasa (por forma, al borde) |
| B2 | 0.1 | 1.56 | 4 | 58 | 2.8×10⁻¹⁶ | 479 keV | 0.137 | 4.75 keV, 0.19 | 0.83 | pasa |
| B4a | 0.1 | 1.56 | 1 | 155 | 2.5×10⁻¹⁸ | 56 keV | 0.326 | 4.0 keV, 0.61 | 0.64 | pasa (por forma) |
| B4a | 0.03 | 1.56 | 1 | 44 | 3.3×10⁻¹⁷ | 187 keV | 0.115 | 10 keV, 0.52 | 1 | pasa |
| B4a | 0.1 | 1 | 1 | 97 | 6.0×10⁻¹⁸ | 87 keV | 0.237 | 2.0 keV, 0.12 | 0.30 | pasa |

Notas:
- **X = K_eff·c·r** = (8π/27)·R_ρ·ξℓη·(c/c_E)·r. Las filas con el mismo X dan el mismo δA a ±0.01 (por ejemplo, B2 con η = 1 y R_ρ = 1, y B4a con η = 0.1 y R_ρ = 1). Confirma la lectura de U17c.
- **f_max(m)** es el límite de la forma estricta: 0.14·m^1.1 por debajo de 4 keV, que sube hasta 1 en 5.3 keV.

## 3. La expectativa no se cumplió, y por qué

El preregistro esperaba «RAMA_VIVA_REQUIERE_COLAPSO_DURO»: que B2 con η = 1 fallara, con δA ≈ 0.31.
- **δA acertó:** 0.323.
- **La conclusión no:** la puerta primaria no es δA, sino la forma. La forma acepta una mezcla de 4.3 keV con un 64 % de componente tibia, porque el límite de las mezclas en esa masa es 0.72.
- El umbral X ≈ 95 de U17c (§4 de su informe) era el de δA estricto, no el de la puerta primaria. **El error fue mío al escribir la expectativa,** no del cálculo. La regla preregistrada manda, y el veredicto es el que da.

**Lo que hay que decir al lado del veredicto:**
- **Por δA, la fila primaria está entre los dos umbrales:** 0.323 queda por encima del estricto (0.234) y por debajo del amplio (0.336), con un 4 % de margen.
- **Pasa por la puerta más permisiva de las dos que se informan.** Es el mismo tipo de «pasa por forma, falla por δA» que ya tenían la fila U17b (0.247) y B4a con η = 0.1 en U17c (0.249), pero ahora un 38 % por encima del umbral estricto, no un 6 %.
- **La equivalencia térmica está poco identificada,** como advirtió Codex:
  - filas con δA casi igual se ajustan a puntos muy distintos de la familia: X ≈ 96–103 da (7.8 keV, 0.97) en el control y (2.0 keV, 0.12) con el balance corregido;
  - **la frase «equivale a una WDM de ~8 keV» no es una predicción robusta.** Lo robusto es la supresión: δA ≈ 0.24–0.32 en la rama viva con η = 1.

## 4. Lo que cambia la lectura

- **La época del colapso decide η = 1.** Con la receta de Hiramatsu (energía congelada en H = m_G/3 y momento c·m_G al terminar), r = 2 ya hace fallar η = 1: δA = 0.50 y forma (3.0 keV, 0.86) > 0.47. Con r = 4 hace falta η ≲ 0.3, y esa fila pasa al borde. La derivación del preregistro dice que esa receta es el extremo pesimista, pero no se ha medido.
- **La isocurvatura empuja en la misma dirección.**
  - B2 con η = 1 y R_ρ = 1.56 tiene m_G ≈ 2.6×10⁻¹⁸ eV.
  - **Literatura (solo resúmenes del buscador, sin leer):** con la radiación de cuerdas dominando, las cotas de Lyα para ALPs postinflacionarios de masa constante quedan en m ≳ 10⁻¹⁷ eV en la mayoría de los casos (arXiv:2511.15790); otro análisis da m > 1.7×10⁻¹⁸ eV (arXiv:2603.04401).
  - Las filas con η ≤ 0.3 tienen m_G ≳ 3×10⁻¹⁷ eV.
  - Sin un cálculo propio (U20) no es una exclusión, pero **las dos presiones independientes, Lyα con r ≥ 2 e isocurvatura, apuntan a η ≲ 0.3.**
- **B4a (la rama cómoda)** sigue pasando la puerta primaria con η = 0.1, también por forma y con δA = 0.326, y con η = 0.03 holgadamente. Su η_crit en la rejilla no cambia, pero la fila η = 0.1 queda tan al borde como la de B2 con η = 1.

## 5. Qué significa

- **No hay exclusión nueva.** Con las correcciones de Codex, la regla preregistrada sigue dejando viva la rama de ξ constante con un colapso tipo U(1).
- **El margen ha desaparecido.**
  - Con la hipótesis coherente (R_ρ = R_Γ), la rama viva pasa solo por la puerta más permisiva.
  - Falla si el colapso termina tarde (r ≥ 2).
  - Probablemente choca con la isocurvatura.
  - **La ballena es más robusta con un colapso «duro» (η ≲ 0.3):** m_G ≈ 3×10⁻¹⁷–3×10⁻¹⁶ eV, δA ≈ 0.04–0.2, y pasa todas las variantes de r y R_ρ probadas.
- **La pregunta física decisiva es la misma para las dos ramas:** qué fracción de la energía del colapso de RP³ acaba en cuantos de p ~ m_G, y cuándo termina el colapso. Las dos cosas se miden en la red (U19b, en WSL).

## 6. Lo que U17d no resuelve

- **R_ρ y r** son hipótesis declaradas, no medidas.
- **Isocurvatura y minicúmulos:** U20 está por preregistrar. Las masas de minicúmulos (17–6×10⁴ M_⊙) son escalas candidatas.
- **Solo c = 3, σ_ln = 0.5, q = 0.87, x₀ = 10.**
- **U15c** se redefine con este resultado (ESPEC_U15c, enmienda 2).
