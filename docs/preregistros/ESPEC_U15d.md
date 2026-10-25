# U15d · El borde superior de la banda frente a la verosimilitud completa de Lyman-α, con la población del colapso medida en U19b · Preregistro

*20 de octubre de 2026, 10:41 (Madrid). Escrito antes de generar ninguna tabla ni cadena de U15d; la hora que cuenta es la de su commit. Las tablas se generan en la nube y las cadenas, en WSL, con la tubería de U15b/U15c sin cambios.*

**Lo que ya se ha visto, y se declara:**
- **U17f completo** (`docs/informes/INFORME_U17f.md`): el δA de cribado de cada fila del borde, con la población medida. Por ejemplo, en F₀ es 0.378 con Y = 1.20, 0.468 con Y = 1.70 y 0.510 con Y = 1.98; en 0.85·F₀ con Y = 1.98, 0.421; en 0.7·F₀ con Y = 1.98, 0.323.
- **U15b y U15c completos** (`docs/informes/INFORME_U15b_U15c.md`):
  - la validación da f₉₅ = 0.27 (1 keV) y 0.42 (2 keV);
  - los bordes de U17d, F ≤ 0.5·F₀, son PERMITIDOS con una posterior de s plana.
- **La calibración informativa de ese informe:** el f₉₅ de la validación equivale a δA ≈ 0.39–0.43. Por eso se espera que este cálculo decida algo en el borde actual.
- **No se ha generado ninguna tabla T² de estos casos ni ninguna cadena.**

## Motivo

U17f dejó la banda de F en ≥ 0.7·F₀ para todo el intervalo del colapso. F₀, y con él LISA, queda permitido solo si el X del universo es ≲ 194, y solo por forma. Esa conclusión sale de la puerta de cribado (δA y forma a z = 3) de U12.

La verosimilitud completa (Planck + BAO + Lyman-ABGD, con la cosmología y el IGM marginalizados) validó con mezclas térmicas y no probó nunca ese borde. U15d lo prueba con la población medida en U19b, la de U17e y U17f.

## Casos (`lya/u15b/prepare_tables.py --u15d`)

La población, el número y la rama son los de U17f, sin cambios:
- forma medida `measured_shape("low")` de U17e, con x̄ = 5.83;
- rama viva B2, q = 0.87, x₀ = 10, la emisión de U12, R_ρ = 1.56;
- K = (8π/27)·R_ρ·ξ_B2·ℓ·Y/x̄, con m_G autoconsistente (`u17e_class.row`);
- meseta de suma (`WHALEFALL_PLATEAU=sum`), como U17e y U17f.

| Caso | F/F₀ | Y | X del universo en F₀ | δA de cribado (U17f) | Qué prueba |
|---|---|---|---|---|---|
| `U15d_F1_Y1.20` | 1 | 1.20 | ≈ 195 (∝ ℓ más las cuerdas) | 0.378 | **D1:** si vuelve F₀ con la extrapolación que mejor apoyan los datos |
| `U15d_F1_Y1.70` | 1 | 1.70 | ≈ 275 (escalado de U17) | 0.468 | **D2:** F₀ con la extrapolación central de U17 |
| `U15d_F0.85_Y1.98` | 0.85 | 1.98 | ≈ 320 (la peor del intervalo) | 0.421 | **D3:** el borde de 0.85·F₀ en el peor caso |
| `U15d_F0.7_Y1.98` | 0.7 | 1.98 | ≈ 320 | 0.323 | **D4:** que el borde vigente de 0.7·F₀ resista en el peor caso |

**Tablas:** T²(k; s) = P/P_ΛCDM a z = 0 con CLASS y el parche.
- ω_ncdm = s·f_nc·Ω_DM y ω_cdm = Ω_DM − ω_ncdm;
- s ∈ {0, 0.25, 0.5, 0.75, 0.9, 1, s_max}, con s_max = min(1.25, 1/f_nc), como en U15c.

**Cadenas:** `Lya_abgd_whale` sin cambios, con prior plano de s en [0, s_max], 4 cadenas por caso, el muestreador de la Nota 8 de U15b (`boost`) y las reglas de lectura de las notas 9–12.
```bash
CASES="U15d_F1_Y1.20 U15d_F1_Y1.70 U15d_F0.85_Y1.98 U15d_F0.7_Y1.98" bash lya/u15b/run_wsl.sh chains
```

**Convergencia:** R − 1 < 0.03 en todos los parámetros, también por mitades, y 4 cadenas, como en U15b y U15c.

**Validación:** la de U15b, que pasó. No se repite.

## Criterio por caso (el de U15c, sin cambios; `lya/u15c/analyze_u15d.py`)

- **PERMITIDA:** p(s = 1)/p_max ≥ e⁻² con 30 bins sobre el prior.
- **EXCLUIDA:** p(s = 1)/p_max < e⁻².
- **INDECISO:** si la etiqueta cambia con 20, 30 o 45 bins, o al quitar una cadena.
- **No se escribe ninguna etiqueta** hasta que el caso converja.

## Veredicto

Tiene dos partes, que se informan juntas.

**F₀ (casos D1 y D2):**

| Veredicto | Condición |
|---|---|
| **F0_PERMITIDO_LYA_COMPLETO** | D1 y D2 PERMITIDOS: F₀ pasa Lyα completo con X ≲ 275 |
| **F0_PERMITIDO_SOLO_CON_X_195** | D1 PERMITIDO y D2 EXCLUIDO |
| **F0_EXCLUIDO_LYA_COMPLETO** | D1 EXCLUIDO |

**El borde en el peor caso, X ≈ 320 (casos D3 y D4):**

| Veredicto | Condición |
|---|---|
| **BORDE_0.85_SOSTENIDO** | D3 PERMITIDO |
| **BORDE_0.7_SOSTENIDO** | D3 EXCLUIDO y D4 PERMITIDO |
| **BORDE_BAJA_DE_0.7** | D4 EXCLUIDO |

**Global:** «<parte F₀>; <parte borde>».
- Si algún caso queda INDECISO, su parte es INCONCLUSA.
- Si D1 sale EXCLUIDO y D2 PERMITIDO, el resultado no es monótono en X: la parte F₀ es **INCONCLUSA (no monótona)**, y se informa.

**Se informan también:**
- la posterior de s y p(s = 0)/p_max;
- Δχ²_min con s ≈ 1, que es informativo;
- T² en k = 10 h/Mpc con s = 1.

**U15d solo prueba Lyα.** Las demás puertas de la banda no cambian:
- la isocurvatura de U20b, que no recorta en el borde;
- la superradiancia estelar, que no recorta en el borde;
- la cota secundaria de supermasivos, que excluye F₀ en todo el intervalo.

Si F₀ pasara aquí, la banda principal llegaría a F₀ con esa extrapolación, pero con la cota de supermasivos seguiría en 0.5–0.7·F₀.

## Lo que se espera, escrito antes de calcular

**Con la calibración de U15b** (el f₉₅ de la validación equivale a δA ≈ 0.39–0.43), y sabiendo que la comparación es solo orientativa:
- **D1** (δA 0.378): PERMITIDO;
- **D2** (0.468): EXCLUIDO;
- **D3** (0.421): incierto;
- **D4** (0.323): PERMITIDO.

**Veredicto esperado:** F0_PERMITIDO_SOLO_CON_X_195, con BORDE_0.7_SOSTENIDO o BORDE_0.85_SOSTENIDO.

## Lo que U15d no resuelve

- **Cuál es la extrapolación correcta del colapso.** Eso lo dirá U19c, la tercera masa.
- **La cota de supermasivos y la normalización C de la isocurvatura.**
- **T² en la cosmología fiducial,** como en U15b. La auditoría acotó ese error en |ΔlnL| ≈ 0.4–0.7.
- **Con f_nc ≈ 0.98, el modelo queda cerca del borde de lo que cubren las simulaciones de Hooper et al.** Se dirá si las comprobaciones internas de la verosimilitud avisan.

## Nota · 20 de octubre de 2026, 12:25 (Madrid): las tablas, hechas; ninguna cadena todavía

- **Las 28 tablas de CLASS** están en `lya/u15b/tables/U15d_*.npz` y `manifest_u15d.json`. Se generaron en la nube con `prepare_tables.py --u15d --workers 4`, en ≈ 1.6 h.
- **Coinciden con las filas de U17f:** m_G, X y f_nc reproducen esas filas exactamente.
- **T² con s = 1 en k = 10 h/Mpc,** visto al hacer las tablas; son las entradas, no un resultado:

  | Caso | T² |
  |---|---|
  | D1 | 0.972 |
  | D2 | 0.946 |
  | D3 | 0.961 |
  | D4 | 0.982 |

  Son supresiones del 2–5 %, frente a < 0.15 % en U15c y al 6–15 % que separan los límites del 95 % de la validación. Esta vez la verosimilitud debería poder discriminar.
- **En WSL:** `git pull` de esta rama, o de `main` cuando se fusione, y:
  ```bash
  CASES="U15d_F1_Y1.20 U15d_F1_Y1.70 U15d_F0.85_Y1.98 U15d_F0.7_Y1.98" bash lya/u15b/run_wsl.sh chains
  bash lya/u15b/run_wsl.sh status U15d                    # convergencia
  python -B lya/u15c/analyze_u15d.py ~/u15b/chains        # cuando converjan
  ```

## Nota · 24 de octubre de 2026, 09:24 (Madrid): cómo corrieron las cadenas y el resultado

**Lo que pasó con las cadenas en WSL,** sin efecto sobre el resultado:
- **20 de octubre:** se lanzaron las cuatro cadenas de cada caso con `run_wsl.sh chains`.
- **21 de octubre: `boost`** (el muestreador de la Nota 8 de U15b, como decía este ESPEC). La primera vez no hizo nada, porque `boost U15d` no encontraba los casos fuera de la lista por defecto. Se arregló (`CHANGELOG`, 21 oct) y se aplicó ese mismo día. Antes del `boost`, R − 1 iba de 0.15 a 10.5; el 22 de octubre ya estaba por debajo de 0.05 en los cuatro casos.
- **22 de octubre: un `chains` sin `CASES`** reanudó por error U15b y U15c, que ya estaban cerrados. Se paró enseguida. Las copias parciales que llegó a escribir MontePython se comprobaron byte a byte contra sus originales y se apartaron a `restart_accidental_20261022/`, sin borrar nada. **U15d no se tocó**, y U15b y U15c volvieron a dar exactamente los números del 21 de octubre. Desde entonces, `chains` sin `CASES` solo lanza U15d.
- **23 de octubre:**
  - Tras un reinicio de la máquina, se relanzó `chains` sobre U15d, y se paró a los pocos minutos.
  - **`status` dio un «converged» falso.** El lector contaba como cadenas las copias parciales de esa reanudación cortada: 5–7 en lugar de 4. Al repetir puntos, R − 1 bajaba.
  - **Lo que se hizo:**
    - las copias se apartaron a `restart_accidental_20261023/`, otra vez comprobadas byte a byte y sin borrar;
    - las originales volvieron a su sitio;
    - la convergencia pasó a exigir exactamente 4 cadenas (PR 28).
  - **Con las 4 cadenas reales, los casos no habían convergido:** R − 1 por mitades = 0.028–0.031. Por eso se reanudó `chains` de forma normal. La decisión se tomó solo con `status`, que únicamente informa de la convergencia; no se miró ninguna posterior de s.
- **24 de octubre:** con los cuatro casos convergidos, se ejecutó `python -B lya/u15c/analyze_u15d.py ~/u15b/chains` con el código de `main`. **Es el único análisis con veredicto que se ha hecho.**

**La lectura de las cadenas:**
- 4 cadenas por caso, sin copias de reanudación ni bifurcaciones;
- entre 178 720 y 195 821 puntos por caso tras los cortes;
- las reglas de las notas 9–12 de U15b quitan entre 0 y 8 712 puntos no markovianos por caso, los del cambio de propuesta;
- R − 1 ≤ 0.017, y ≤ 0.026 por mitades.

**El resultado (`lya/results/u15d/U15D_RESULT.json`): F0_EXCLUIDO_LYA_COMPLETO; BORDE_0.7_SOSTENIDO.**

| Caso | p(s = 1)/p_max |
|---|---|
| D1 | 0.111 |
| D2 | 0.047 |
| D3 | 0.114 |
| D4 | 0.238 |

**Ninguna etiqueta cambia con 20 o 45 bins ni al quitar una cadena.** El informe es `docs/informes/INFORME_U15d.md`.

**Las comprobaciones internas de la verosimilitud.** Este ESPEC pedía decir si avisan con f_nc ≈ 0.98. El JSON del analizador no recoge esos avisos, y los registros de MontePython siguen en WSL, sin revisar. **Queda pendiente;** se añadirá en otra nota.

## Nota · 24 de octubre de 2026, 13:11 (Madrid): las comprobaciones internas de la verosimilitud no saltaron

Esto cierra lo que la nota anterior dejó pendiente.

**Dónde avisa la verosimilitud.** Se leyó el código de `Lya_abgd` (MontePython 5b365b0). Cuando un modelo se sale de lo que cubre, no escribe ningún «warning». Hace tres cosas:
- devuelve el valor límite, de modo que el punto se rechaza;
- escribe `#Error_cosmo`, `#Error_equiv`, `#Error_abgd_fit` o `#Error_smooth` en la salida de error;
- deja una línea con los parámetros en `Lya_abgd_bin_file.txt`, en la carpeta de cada caso.

Las cuatro comprobaciones son:
- σ₈ y n_eff dentro de la rejilla de simulaciones;
- T²(k) → 1 a escalas grandes;
- que el ajuste α–β–γ–δ converja;
- que el suavizado no cambie T² en más de 0.1.

Si saltaran con preferencia cerca de s = 1, rechazarían esos puntos y empujarían artificialmente hacia la exclusión.

**Lo que se comprobó, en WSL y con las cadenas del autor:**
- se contaron las líneas `#Error_*` de `Lya_abgd_bin_file.txt` por tipo y por tramo de s (s < 0.5, 0.5–0.9, ≥ 0.9), leyendo el orden de los parámetros de cada `log.param`;
- en los registros de `~/u15b/logs/` no aparece ningún «warning» ni «boundary»; esa búsqueda no incluía `#Error`, así que la que cuenta es la del fichero de errores.

**Resultado: ningún error, en ningún caso.** Ni en los cuatro de U15d, ni en U15c, ni en la validación o las filas de U15b. Las comprobaciones internas no saltaron nunca, ni siquiera con f_nc ≈ 0.98, y no hay puntos rechazados por ellas que puedan sesgar p(s = 1). **El veredicto no cambia.**

**De paso, en los registros apareció un error de MontePython que no afecta a los resultados.** Al terminar sus 100 000 pasos, cada cadena intenta borrar el fichero desde el que se reanudó (`mcmc.py`, línea 878). `move_restarted.sh` ya lo había movido a `restarted_from/`, así que el borrado falla con `FileNotFoundError` y `mpirun` corta las demás cadenas del caso.
- Todos los puntos ya estaban escritos.
- Como mucho se pierden los pasos que les faltaban a las otras cadenas.
- Explica por qué `chains` «terminaba» y había que relanzarlo.
