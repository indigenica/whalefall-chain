# U26 · La inflación de polo como potencial por defecto de la cadena · Preregistro

*18 de octubre de 2026, 16:55 (Madrid; la hora que cuenta es la de su commit). Escrito antes de cambiar el código y de ejecutar nada con el polo por defecto. Se ejecuta en la nube (minutos).*

## Motivo

- **U22** (POLO_CAMBIA_NACIMIENTO) y **U22b** (NACIMIENTO_ROBUSTO) mostraron que la meseta puede ser la de la **inflación de polo de orden 4**, V = V₀·[1 − C/√(χ² + C²)] con C = 0.3516. Es un lagrangiano estándar, en lugar de la suma postulada de 10 niveles (el parche D2).
  - El CMB es el mismo: Δn_s = −1.2×10⁻⁴ y Δr/r = −0.2 %.
  - La estadística de nacimientos también: f_ES = 0.153 frente a 0.163.
- **Hoy la cadena calcula con la suma** (S1, S2, la acción común), mientras que los textos dicen que la meseta «puede ser» de polo. El autor decide, el 18 de octubre, que el polo pase a ser el potencial por defecto y que la suma se quede como opción para reproducir lo anterior. U26 comprueba que el cambio es coherente y fija qué pasa si no lo es.
- **Lo que el cambio mueve,** visto antes de escribir este documento con `chain/u22_pole.py` (sin veredicto):
  - casi nada en F, H_* y μ (−0.09 %);
  - **m_χ (+33 %), Γ_χ (+33 %) y T_rh (+15 %),** porque la curvatura de la meseta en el mínimo es otra (1/C² frente a la de la suma);
  - χ_end (−18 %).
  - Una prueba rápida de la raíz de materia oscura (`s4_dark.dm_root`, dos filas) da Δm_G/m_G ≈ +0.4 %. Por eso se espera que la parte oscura no cambie.

## Qué se cambia (código, sin veredicto)

1. **Un solo interruptor:** `WHALEFALL_PLATEAU = pole | sum` (variable de entorno; por defecto `pole`). Lo leen `chain/ucommon.py` (`benchmark()`, `provenance()`, que lo guarda en cada resultado) y `action/src/common.py`. Con `sum`, todo debe reproducir lo guardado hasta hoy.
2. **El potencial:** `chain/vendor_potential.py` exporta `Potential` según el interruptor. La clase de polo es la de U22 (`chain/u22_pole.py`), con `eternal()` y `pzeta()` añadidos.
3. **Entrada nueva congelada:** `data/cadena_cierre/thermal_pole.json`. Es el emparejamiento continuo del pivote de U22 con el polo, a rtol = 2×10⁻¹⁰ como `thermal.json`. Lo escribe `chain/u26_thermal_pole.py`, y se registra en `data/MANIFEST.json`. `benchmark()` lee `thermal_pole.json` o `thermal.json` según el interruptor.
4. **La acción común** (`action/`): con el polo, `common.Potential` es el polo y `model.py` lee `vendor/thermal_pole.json`, que es una copia de la entrada del punto 3. Hay dos cambios más:
   - `portal.py` hace la cota de no negatividad con el polo (punto C4);
   - los resultados del polo van a `action/results_pole/`, para no pisar `action/results/`, y se congelan en `data/cierre_accion_pole/`.
5. **Las etapas y las pruebas:** S1 y S3 se reformulan como dice el apartado siguiente. Las pruebas que congelan números de la suma (χ_et = 415.6, el pivote de `thermal.json`, las raíces de `dark.json`) fijan `sum` de forma explícita, y se añaden sus equivalentes con el polo.
6. **Lo que se queda con la suma, declarado:**
   - los resultados de `lya/`, `audit/` y `network/` ya guardados;
   - para reproducirlos se ejecuta con `WHALEFALL_PLATEAU=sum`, y así se dice en sus README.

## Las filas (con el polo por defecto)

**C1 · Pivote (S2).**
- El emparejamiento con el polo reproduce `thermal_pole.json` (máx. rel. < 10⁻⁵, la puerta actual de S2).
- Con `sum`, S2 sigue reproduciendo `thermal.json`.
- **Pasa** si además |Δn_s| ≤ 0.002 y |Δr/r| ≤ 0.10 frente a la suma (los criterios de CMB de U22).

**C2 · Nacimiento (S1), reformulado.**
- **Control:** el modo 0 (la suma) reproduce bit a bit los 20 hijos viables de C45j y los 70 controles, como hasta hoy. Es la garantía del puerto C, no depende del potencial por defecto.
- **Fila primaria:** la rejilla completa de C45h (1100 filas, frac_rad = 0) en el modo 4 (polo), con el χ_et del V₀ del polo. Se calcula f_ES(polo), con la definición de S1 y U22b, y R = f_ES(polo)/f_ES(suma), con la suma en el modo 0 y su χ_et.
- **Pasa** si el control pasa y R ≥ 0.5 (el umbral de NACIMIENTO_ROBUSTO de U22b).
- **Cambio de enunciado:** el eslabón 2 del contrato deja de decir «los 20 hijos viables son eternos» y dice f_ES. Con el polo solo 10 de esos 20 son eternos. El filtro de perturbaciones de C45j (ζ₀) está calculado con la suma y no se rehace; se declara.

**C3 · Fase eterna (S2).**
- Se calculan con el polo y su V₀: χ_et, el umbral de volumen (P_ζ = 1/6) y la cota de planitud D3, sqrt(V₀/3)/(2π).
- **Pasa** si la meseta sin deformar cumple la cota por encima de χ_et (el mismo indicador que hoy, `undeformed_plateau_meets_bound`).
- El veredicto C45i (deformación FRÁGIL) es de la suma y se arrastra declarado.

**C4 · Acción común (S3).** `make action` con el polo.
- **Pasa** si se cumplen todas estas condiciones:
  - los residuos simbólicos son nulos;
  - las cuatro integraciones de `experiments` cumplen `pass_numerics`;
  - `late_fluid` cumple `pass_numerics`;
  - todos los modos de Q tienen canal de desintegración, con min Γ_Q/Γ_χ > 1;
  - la cota de no negatividad del portal se cumple con el polo.
- **La cota con el polo:** V_polo(χ) ≥ σ²M²·tanh²(χ/M)/(4λ) para todo χ, comprobada en una rejilla logarítmica de |χ| entre 10⁻⁶ y 10⁶ M_Pl y en el límite χ → 0, con V_polo ≈ V₀χ²/(2C²).
- **S3 pasa** con su criterio de hoy, T_tardía/T_rh − 1 < 10⁻³ y la identidad de F < 10⁻⁹. Con el polo, T_tardía se compara con la acción del polo y F con `thermal_pole.json`. La comparación de F con el contrato de `cadena_cierre` (de la suma) se queda como control con `sum`.

**C5 · Parte oscura (S4–S6, U17).**
- **S4:** sus puertas de reproducción comparan el código vectorizado con el original y las raíces con `dark.json`, que es de la suma. Son puertas de código y se evalúan con la suma.
  - Con el polo se calculan las mismas raíces y se informa max |Δlog₁₀ m_G|.
- **S5, S6 y U17:** se rehacen con el polo.
  - **Pasan** si mantienen la etiqueta de veredicto (CIERRE_ESTRECHO con 309 nodos; FALSADA_EN_FORMA_ACTUAL; K_SW_FUERA_DEL_RANGO_TIBIA).
  - **Se arrastran los veredictos de Lyα con CLASS** (U3a, U4a, U5, U7, U12, U15, U17b–e, U21, U20b, U23), sin rehacerlos, si en todas las filas de S5, S6 y U17 que dan m_G se cumple |Δm_G/m_G| ≤ 2 %.

## Veredicto

- **POLO_ADOPTADO:** pasan C1–C5. El polo queda por defecto y la suma como opción.
- **POLO_ADOPTADO_CON_CAMBIOS:** pasan C1 y C2, pero falla alguna de C3–C5.
  - El polo queda por defecto, porque las condiciones de equivalencia de U22 y U22b se cumplen.
  - Cada etapa que cambia se informa con su veredicto nuevo como consecuencia del potencial, y el contrato lo refleja.
  - Si falla el arrastre de Lyα (> 2 %), se listan las pruebas que hay que rehacer, y sus veredictos se marcan «calculados con la suma».
- **POLO_NO_ADOPTADO:** falla C1 o C2. Se deja la suma por defecto, el interruptor se queda y se informa por qué.

## Lo que se espera, escrito antes de calcular

**POLO_ADOPTADO.**
- C1 y C2 son, por construcción, U22 y U22b.
- En C4, los márgenes actuales son enormes (min Γ_Q/Γ_χ ≈ 4×10⁵; cota del portal 1.5×10⁻¹⁷ frente a 8.6×10⁻¹¹). Con +33 % en σ y Γ_χ siguen siéndolo.
- En C5, Δm_G ≈ 0.4 %.
- **El riesgo está en C4:** χ_end y m_χ cambian, y la integración homogénea con Q y NPQT no se ha hecho nunca con el polo.

## Lo que U26 no resuelve

- **El filtro de perturbaciones de C45j** (la «fracción viable» de 3.6 %) sigue siendo el de la suma.
- **El veredicto de deformación C45i** (D3 FRÁGIL) tampoco se rehace con el polo.
- **La conjetura de la distancia de la Ciénaga** sigue sin evitarse: el polo protege la meseta frente a correcciones radiativas, no frente al campo transplanckiano.
- **No se rehacen las redes ni las pruebas con CLASS;** solo se arrastran si C5 lo permite.

## Enmienda 1 · 18 de octubre de 2026, 18:11 (Madrid): dos errores del verificador, corregidos después del primer cálculo

La revisión de Codex en la PR #11 encontró dos errores en `chain/u26_pole_default.py`. Los dos hacían el verificador menos estricto de lo que pide este ESPEC. **Se corrigen después de haber calculado el veredicto** (POLO_ADOPTADO, commit `1656a97`), y el veredicto se recalcula. **Los criterios no cambian.**

1. **C4 se evaluaba sobre la copia congelada** `data/cierre_accion_pole/`, no sobre la acción que acaba de ejecutar `make u26` (`action/results_pole/`). Ahora C4 se evalúa sobre la ejecución nueva, y exige además que la copia congelada, que es la que lee S3, coincida con ella: los mismos indicadores y una diferencia relativa ≤ 10⁻⁶ en la cota del portal, en Γ_Q/Γ_χ y en las temperaturas tardías.
2. **C5 comparaba solo las masas presentes en las dos cadenas.** Una fila que perdiera su m_G con el polo se habría saltado sin avisar. Ahora C5 exige que los conjuntos de m_G de la suma y del polo sean idénticos.

**Resultado recalculado:** **POLO_ADOPTADO**, igual.
- La copia congelada coincide con la ejecución nueva.
- Los conjuntos de m_G son idénticos: 373, 72 y 98 valores en S5, S6 y U17.

## Enmienda 2 · 18 de octubre de 2026, 20:57 (Madrid): una puerta más estricta, una comprobación de F con sentido y lo que reproduce la suma, tras una auditoría independiente

Una auditoría independiente reimplementó la física del polo y comprobó la reproducción de la suma, sin encontrar nada que cambie el veredicto. Se corrige:

1. **C2:** además del control, la rejilla en modo 0 tiene que reproducir lo guardado (diferencia relativa ≤ 10⁻⁸, el criterio INVALIDO de U22b). Es más estricta que el texto original. Su valor ya visto es 0.
2. **C4 / S3:** con el polo, la identidad de F comparaba `thermal_pole.json` consigo mismo. Ahora se compara con la copia que lee la acción común (`action/vendor/thermal_pole.json`).
3. **La cota del portal** se refina con un minimizador alrededor del mínimo de la rejilla, porque un mínimo de rejilla no es una cota inferior. El margen es de 6×10⁶.
4. **El contrato con la suma no lleva la fila U26.** El veredicto de U26 es del polo, y con eso el contrato de la suma vuelve a ser el de antes de U26.
5. **C3 no puede fallar por construcción.** P_ζ > 1 por encima de χ_et equivale a la cota D3, y viene así de antes de U26. Se declara.
6. **El verificador** (`chain/u26_pole_default.py`) se escribió después de la primera ejecución de la cadena con el polo. Su lógica era la de este ESPEC, fijada antes.
7. **Se rehace `make u26`** con el código corregido. Lo que vuelva a salir se informa en `INFORME_U26.md`.
