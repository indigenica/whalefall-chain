# U26 · La inflación de polo pasa a ser el potencial de la cadena: nada de lo medido cambia

*18 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U26.md` (commit `d70a495`, anterior a cualquier cambio de código).
- *Código:*
  - el interruptor `WHALEFALL_PLATEAU` (`chain/ucommon.py`, `action/src/common.py`);
  - `chain/pole_potential.py`, `chain/u26_thermal_pole.py` y `chain/u26_pole_default.py`;
  - `make u26`.
- *Resultados:*
  - `results/u26_pole_default.json`;
  - la cadena con el polo en `results/`;
  - la misma cadena con la suma en `results/plateau_sum/`;
  - la acción con el polo en `data/cierre_accion_pole/`.
- *Ejecución:* en la nube, unos 12 min.

## 0. En una frase

**Veredicto preregistrado: POLO_ADOPTADO,** como se esperaba. La meseta de la cadena es desde hoy la **inflación de polo de orden 4**, V = V₀·[1 − C/√(χ² + C²)] con C = 0.3516, en lugar de la suma postulada de 10 niveles (el parche D2).
- Pasan las cinco filas.
- Ningún veredicto de la cadena cambia.
- Las masas de la materia oscura se mueven un 0.4 % como mucho.
- La suma queda como opción y **reproduce exactamente todo lo guardado antes de U26** (`WHALEFALL_PLATEAU=sum`).

## 1. Las cinco filas

| Fila | Qué se mide | Polo | Suma | Puerta | |
|---|---|---|---|---|---|
| **C1 · Pivote** (S2) | n_s, r | 0.97277, 0.00502 | 0.97289, 0.00502 | \|Δn_s\| ≤ 0.002, \|Δr/r\| ≤ 10 % | pasa (Δn_s = −1.2×10⁻⁴, Δr/r = −0.2 %) |
| **C2 · Nacimiento** (S1) | f_ES en las 1100 filas de C45h | 0.153 | 0.163 | R ≥ 0.5 | pasa (R = 0.94) |
| **C3 · Fase eterna** (S2) | χ_et; planitud por encima de χ_et | 541 M_Pl; cumple | 416 M_Pl; cumple | cumple la cota D3 | pasa |
| **C4 · Acción común** (S3) | las integraciones con Q y NPQT, los anchos de Q, el portal | todo pasa | todo pasa | las de la acción, más min Γ_Q/Γ_χ > 1 | pasa (Γ_Q/Γ_χ ≥ 3×10⁵; cota del portal 2.7×10⁻¹⁷ frente a 1.7×10⁻¹⁰) |
| **C5 · Parte oscura** (S4–S6, U17) | etiquetas; Δm_G | las mismas; ≤ 0.42 % | — | mismas etiquetas; \|Δm_G\| ≤ 2 % | pasa: se arrastran los veredictos de CLASS |

**Lo que se mueve al cambiar de meseta.** Casi nada de lo observable; sí lo que depende de la curvatura en el mínimo:

| | Polo / suma − 1 |
|---|---|
| F, H_*, μ | −0.09 % |
| V₀ | −0.18 % |
| T_rh | **+15 %** |
| m_χ, Γ_χ | **+33 %** (m_χ² = V₀/C²) |
| χ al final de la inflación | −18 % |
| m_G (S5, S6, U17; 543 valores) | +0.40 a +0.42 % |

T_rh entra en la materia oscura solo como la temperatura en que empieza la emisión de la red, y por eso m_G casi no cambia.

## 2. Qué cambia en cada etapa

- **S1 se enuncia de otra forma.**
  - Antes decía que los 20 hijos viables de C45j son eternos. Con el polo, solo 10 de esos 20 lo son: χ_et sube de 416 a 541 y los aterrizajes concretos cambian.
  - La fracción de nacimientos eternos es la misma (U22b). Por eso **S1 se mide ahora con f_ES**, y el contrato lo dice así.
  - La reproducción bit a bit de los 20 hijos con la suma se queda como control del puerto C.
- **S2:** el pivote con el polo reproduce la entrada congelada nueva, `data/cadena_cierre/thermal_pole.json` (máx. rel. 1.1×10⁻⁸), y la suma sigue reproduciendo `thermal.json`.
- **S3:** la acción común se rehace con el polo (`make u26`).
  - La temperatura tardía coincide con la T_rh del polo a 4×10⁻⁷.
  - ΔN_eff sigue siendo despreciable (3.4×10⁻⁹).
  - Las cuatro integraciones con Q y NPQT pasan. La de estrés con η_Q = 1000 lleva la fracción de energía de Q al 0.53 % (0.53 % con la suma).
- **El portal visible:** con la suma, la prueba de no negatividad usaba el primer nivel de la suma (V₀/Z). Con el polo, se usa el ínfimo de V_polo(χ)/tanh²(χ/M), que es 1.7×10⁻¹⁰ en χ ≈ 41 M_Pl. La corrección negativa es 2.7×10⁻¹⁷, siete órdenes por debajo.
- **S4:** sus puertas comparan el código con archivos de la suma (`dark.json`), así que se evalúan con la suma. Con el polo, las mismas raíces de m_G se desplazan +0.0015 a +0.0017 dex.
- **S5, S6 y U17:**
  - las etiquetas son las mismas (CIERRE_ESTRECHO con 309 nodos, 178 conservadores y 112 estrictos; FALSADA_EN_FORMA_ACTUAL; K_SW_FUERA_DEL_RANGO_TIBIA);
  - m_G sube un 0.4 % en todas las filas.

## 3. Una diferencia física que conviene saber

**Lejos de la ventana del CMB, el polo es mucho menos plano que la suma.**
- V′/V cae como C/χ² en el polo y de forma exponencial en la suma. Con χ = 10³ M_Pl, el cociente frente a la cota D3 es 0.29 con el polo y 0.010 con la suma. Con χ = 10⁴ M_Pl, es 3×10⁻³ frente a 3×10⁻³³.
- **Las dos cumplen la cota D3** por encima de χ_et, y no cambia ni el nacimiento ni el CMB.
- Pero la «meseta jerárquica» de la suma era una propiedad postulada a medida; el polo la sustituye por la cola estándar 1 − C/χ.
- **El veredicto de deformación C45i (D3 FRÁGIL) es de la suma** y no se ha rehecho con el polo.

## 4. Reproducibilidad

- `WHALEFALL_PLATEAU=sum make chain` reproduce todos los resultados de S1–S6, U17 y el contrato de antes de U26: diferencia relativa máxima 0, salvo el tiempo de ejecución. Están en `results/plateau_sum/`.
- `WHALEFALL_PLATEAU=sum make action` reproduce `action/results_reference/` con diferencias solo en los diagnósticos de error de redondeo, como ya decía su README.
- `make test` pasa con los dos valores del interruptor (14 pruebas). Las pruebas que congelan números de la suma (χ_et = 415.6, `thermal.json`, `dark.json`, U17d) la fijan explícitamente, y hay equivalentes con el polo (χ_et = 541.0, `thermal_pole.json`).
- **Los resultados ya guardados de `lya/`, `network/` y `audit/` se calcularon con la suma.** Para reproducirlos exactamente, ejecútalos con `WHALEFALL_PLATEAU=sum`. C5 muestra que con el polo cambiarían como mucho un 0.4 % en m_G, así que sus veredictos se arrastran.

## 5. Lo que U26 no resuelve

- **El filtro de perturbaciones de C45j** (la fracción viable de 3.6 %) sigue siendo el de la suma.
- **El veredicto de deformación C45i** (D3 FRÁGIL) sigue siendo el de la suma.
- **La conjetura de la distancia de la Ciénaga** sigue sin evitarse. El polo protege la meseta frente a correcciones radiativas, no frente a un recorrido transplanckiano del campo (χ_et ≈ 541 M_Pl).
- **No se rehacen las redes ni las pruebas con CLASS.** Se arrastran porque m_G cambia un 0.4 %.

## Nota · corrección del verificador (enmienda 1 del ESPEC)

Codex encontró dos errores en el verificador:
- C4 leía la copia congelada de la acción en vez de la ejecución nueva;
- C5 no exigía que las masas de la suma y del polo fueran las mismas filas.

Corregidos después del primer cálculo, **el veredicto no cambia:** la copia coincide con la ejecución nueva y los conjuntos de m_G son idénticos.

## Corrección · 18 de octubre de 2026, 20:57 (Madrid), tras una auditoría independiente

**El veredicto no cambia.** La auditoría reimplementó la física del polo: n_s, r, V₀, χ_end, m_χ, χ_et = 540.95 y la cota del portal coinciden. Tampoco encontró diferencias numéricas entre la suma y los resultados de antes de U26. Lo que se corrige:

1. **§4 se pasaba al decir «con `sum` se reproduce todo, incluidos `lya/`, `network/` y `audit/`».**
   - `network/` no depende de la meseta.
   - **`lya/`:** varios scripts leían `results/` con una ruta fija, y desde U26 esa carpeta es la del polo. Ahora leen `results/plateau_sum/` con la suma, y `make lya` y `make u4a` fijan la suma.
   - **`audit/`:** los scripts de las auditorías no se tocan, porque son código de los auditores. Leen `results/` con la ruta fija, así que para reproducir exactamente sus salidas hay que usar la versión anterior a U26 (`git checkout 3d38b3b`).
2. **Errata de §1:** m_G sube entre +0.35 % y +0.42 %, no entre +0.40 % y +0.42 %.
3. **La reproducción la comprueba ahora un script, también en el CI:** `chain/check_reproduction.py`, en `results/` y en `results/plateau_sum/`. Antes no la comprobaba nada.
4. **Enmienda 2 del ESPEC:**
   - la rejilla en modo 0 es ahora una puerta;
   - la identidad de F se comprueba contra la copia de la acción;
   - la cota del portal está refinada;
   - el contrato de la suma no lleva U26;
   - se declara que C3 no puede fallar por construcción y que el verificador se escribió después de la primera ejecución.
