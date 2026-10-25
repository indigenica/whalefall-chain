# U15d · El borde superior frente a la verosimilitud completa de Lyman-α: F0_EXCLUIDO_LYA_COMPLETO; BORDE_0.7_SOSTENIDO

*24 de octubre de 2026, 09:24 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U15d.md`, con la nota del 24 de octubre sobre cómo corrieron las cadenas.
- *Código:*
  - `lya/u15b/prepare_tables.py --u15d` (las tablas);
  - `lya/u15b/run_wsl.sh` (las cadenas);
  - `lya/u15c/analyze_u15d.py`, con el lector y la convergencia de `lya/u15b/analyze.py` (exactamente 4 cadenas desde la PR 28).
- *Resultado:* `lya/results/u15d/U15D_RESULT.json` y las tablas T²(k; s) de `lya/u15b/tables/U15d_*.npz`.
- *Ejecución:*
  - **las tablas:** 28 filas de CLASS en la nube, el 20 de octubre;
  - **las cadenas:** MCMC con MontePython (Planck lite + BAO + Lyman-ABGD) en WSL, del 20 al 24 de octubre, con 4 cadenas por caso y el `boost` desde el 21 de octubre. Las cadenas no están en el repositorio: el autor las conserva en WSL (`~/u15b/chains`).
  - **el análisis:** lo ejecutó el autor el 24 de octubre con el código de `main`, y aquí se sube su JSON.

## 0. En una frase

**Veredicto preregistrado: F0_EXCLUIDO_LYA_COMPLETO; BORDE_0.7_SOSTENIDO.**

Con la población del colapso medida en U19b, la verosimilitud completa de Lyα:
- **excluye F₀ incluso con X ≈ 193,** el caso más favorable de los probados;
- **excluye 0.85·F₀ con X ≈ 320;**
- **permite 0.7·F₀ con X ≈ 320.**

**F₀ queda excluido por dos caminos independientes:** por la extrapolación del colapso (U19d: X ≥ 244) y ahora por Lyα completo aunque X fuera 193. Con él quedan fuera el benchmark y LISA.

**El borde de 0.7·F₀ con X ≈ 320 resiste.** La puerta de cribado de U17f acertó en ese borde, pero en F₀ con X ≈ 193 era más permisiva que la verosimilitud completa. Ahí dejaba pasar F₀ solo por forma, y Lyα completo lo excluye.

## 1. Los casos

La población, el número y la rama son los de U17f:
- la forma medida, con x̄ = 5.83;
- B2, q = 0.87, x₀ = 10, R_ρ = 1.56;
- meseta de suma.

Las tablas reproducen las filas de U17f exactamente.

| Caso | F/F₀ | Y | X del universo | m_G (eV) | f_nc | δA de cribado (U17f) | T² con s = 1 en k = 10 h/Mpc |
|---|---|---|---|---|---|---|---|
| D1 `U15d_F1_Y1.20` | 1 | 1.20 | 193 | 5.8×10⁻¹⁸ | 0.971 | 0.378 | 0.972 |
| D2 `U15d_F1_Y1.70` | 1 | 1.70 | 275 | 2.9×10⁻¹⁸ | 0.979 | 0.468 | 0.946 |
| D3 `U15d_F0.85_Y1.98` | 0.85 | 1.98 | 318 | 4.2×10⁻¹⁸ | 0.982 | 0.421 | 0.961 |
| D4 `U15d_F0.7_Y1.98` | 0.7 | 1.98 | 315 | 9.4×10⁻¹⁸ | 0.982 | 0.323 | 0.982 |

El prior de s es plano en [0, s_max], con s_max = 1/f_nc ≈ 1.02–1.03. s = 1 es el modelo, y s = 0 es CDM pura.

## 2. Convergencia

| Caso | Cadenas | Puntos (tras los cortes) | R − 1 máx. | Por mitades | ESS de s |
|---|---|---|---|---|---|
| D1 | 4 | 188 880 | 0.017 | 0.026 | 3 060 |
| D2 | 4 | 178 720 | 0.014 | 0.018 | 1 380 |
| D3 | 4 | 186 309 | 0.016 | 0.018 | 2 645 |
| D4 | 4 | 195 821 | 0.010 | 0.014 | 3 613 |

- **En los cuatro casos:** ninguna copia de reanudación ni bifurcación, y exactamente 4 cadenas, como exige la PR 28.
- **Las reglas de lectura de U15b** (notas 9–12) quitan 7 606, 0, 8 712 y 5 397 puntos no markovianos, los del cambio de propuesta del `boost`.
- **Los incidentes de las cadenas** (dos `chains` accidentales y un «converged» falso el 23 de octubre) están en la nota del 24 de octubre del ESPEC. No afectan a estas cadenas: las copias parciales se apartaron tras comprobarlas byte a byte, y la única lectura con veredicto es la del 24 de octubre.

## 3. El resultado por caso

| Caso | p(s = 1)/p_max (30 bins) | 20 / 45 bins | Quitando una cadena | Estado |
|---|---|---|---|---|
| D1 (F₀, X ≈ 193) | **0.111** | 0.116 / 0.126 | 0.106–0.117 | **EXCLUIDA** |
| D2 (F₀, X ≈ 275) | **0.047** | 0.045 / 0.044 | 0.038–0.057 | **EXCLUIDA** |
| D3 (0.85·F₀, X ≈ 320) | **0.114** | 0.112 / 0.113 | 0.110–0.120 | **EXCLUIDA** |
| D4 (0.7·F₀, X ≈ 320) | **0.238** | 0.239 / 0.242 | 0.225–0.246 | **PERMITIDA** |

El umbral preregistrado es e⁻² = 0.135. **Ninguna etiqueta cambia con los bins ni con el jackknife,** así que no hay INDECISO.

**D1 y D3 están excluidos por poco.** −2 ln(p(1)/p_max) vale 4.39 y 4.34, frente a 4 en el umbral. El valor más alto de las variantes es 0.126 (D1 con 45 bins), todavía por debajo de 0.135. El criterio preregistrado los excluye sin ambigüedad, pero el margen es pequeño.

**Se informan también, sin veredicto:**

| Caso | s₉₅ (cota superior al 95 %) | ⟨s⟩ | p(s = 0)/p_max | Δχ²_min con s ≈ 1 | Puntos con s ≈ 1 |
|---|---|---|---|---|---|
| D1 | 0.85 | 0.37 | 0.94 | 4.7 | 2 388 |
| D2 | 0.63 | 0.27 | 0.95 | 2.7 | 782 |
| D3 | 0.83 | 0.34 | 0.94 | 4.9 | 2 058 |
| D4 | 0.89 | 0.41 | 0.92 | 2.7 | 3 093 |

- **La posterior de s ya no es plana,** como lo era en U15c. Baja de forma monótona de s ≈ 0 a s_max, y al final del prior la densidad es entre 0.12 y 0.34 veces la media. **Los datos prefieren CDM a la rama viva**, con más fuerza cuanto mayor es la supresión.
- **s₉₅ < 1 en los cuatro casos, también en D4.** No es el criterio preregistrado, y no se puede usar como tal: s = 1 está a un 2–3 % del borde del prior. Por eso, una cota superior unilateral queda por debajo de 1 incluso con una posterior plana (con una plana, s₉₅ ≈ 0.97). El cociente de densidades, que es el criterio de U15c, no tiene ese sesgo.
- **Δχ²_min sigue siendo ruidoso,** como en U15c. D2, el caso más excluido, da uno de los Δχ² menores (2.7). Es el mínimo de 800–3 000 muestras frente al de ≈ 190 000, en 25 dimensiones.

## 4. Lo que se esperaba

| Caso | Esperado (calibración de U15b, δA ≈ 0.39–0.43) | Obtenido |
|---|---|---|
| D1 | PERMITIDO | **EXCLUIDO** |
| D2 | EXCLUIDO | EXCLUIDO |
| D3 | incierto | EXCLUIDO |
| D4 | PERMITIDO | PERMITIDO |

**Veredicto esperado:** F0_PERMITIDO_SOLO_CON_X_195, con BORDE_0.7_SOSTENIDO o BORDE_0.85_SOSTENIDO. **Obtenido:** F0_EXCLUIDO_LYA_COMPLETO; BORDE_0.7_SOSTENIDO.

**El fallo de la expectativa está en D1.** La calibración de U15b, que traducía el f₉₅ de las mezclas térmicas a δA ≈ 0.39–0.43, era demasiado permisiva para esta población.

## 5. Lectura informativa: dónde corta Lyα completo en δA

*Esta lectura se hace después de ver los datos y no tiene estatus de veredicto.*

**El umbral, en δA de cribado:**
- D4 (δA = 0.323) pasa, y D1 (0.378) no. **Si p(1)/p_max baja con δA, el umbral está entre 0.32 y 0.38;** interpolando linealmente entre D4 y D1, ≈ 0.37.
- **δA no ordena los casos del todo:** D1 (0.378) y D3 (0.421) dan casi el mismo p(1)/p_max (0.111 y 0.114). La verosimilitud ve toda la forma de T²(k) y su dependencia de z, no solo el δA a z = 3. Interpolando entre D4 y D3 el umbral saldría ≈ 0.40.
- **La lectura prudente es la de los casos:** δA ≤ 0.32 pasa y δA ≥ 0.38 no. Es más estricto que la calibración de U15b (0.39–0.43), y la referencia ancha de la puerta de cribado (0.336) cae dentro de ese intervalo.

**Qué implica para las filas de U17f y U17g que pasaban «solo por forma»** (δA entre 0.336 y la puerta de forma):

| Fila | X | δA | Con el umbral de U15d |
|---|---|---|---|
| U17f F₀, Y = 1.20 | 193 | 0.378 | excluida (es D1) |
| U17f 0.85·F₀, Y = 1.70 | 275 | 0.383 | probablemente excluida (δA ≥ D1) |
| U17f F₀, Y = 1.05 | 170 | 0.347 | en la zona 0.32–0.38, sin decidir |
| U17g 0.6·F₀, Y = 3.0 | 486 | 0.350 | en la zona 0.32–0.38, sin decidir |
| U17g 0.5·F₀, Y = 4.25 | 689 | 0.346 | en la zona 0.32–0.38, sin decidir |
| U17g 0.35·F₀, Y = 8.44 | 1 367 | 0.336 | en la zona 0.32–0.38, sin decidir |

**Las filas del borde que pasan con δA ≤ 0.32 no cambian:** 0.85·F₀ con X ≈ 220 (0.328, en el límite), 0.6·F₀ con X ≈ 405 (0.301), **0.5·F₀ con el X central de U19d (556; δA = 0.288)** y 0.35·F₀ con X ≈ 890 (0.231).

**El borde superior de la banda con esta lectura,** según el X del universo:

| X | Puerta de cribado (U17f + U17g) | Con el umbral de U15d (informativo) |
|---|---|---|
| ≈ 193 | F₀ (solo por forma) | **< F₀** (D1, preregistrado) |
| ≈ 220 | 0.85·F₀ | 0.85·F₀ (δA 0.328, en el límite) |
| ≈ 275 | 0.85·F₀ | probablemente 0.7·F₀ |
| ≈ 320 | 0.7·F₀ | **0.7·F₀** (D3 y D4, preregistrado) |
| ≈ 400–490 | 0.6·F₀ | 0.6·F₀ (X ≈ 405); 0.5–0.6·F₀ (X ≈ 486) |
| **≈ 556** | **0.5·F₀** | **0.5·F₀** |
| ≈ 690 | 0.5·F₀ | 0.35–0.5·F₀ |
| ≈ 890–1 370 | 0.35·F₀ | 0.35·F₀ (X ≈ 890); 0.25–0.35·F₀ (X ≈ 1 370) |

**En el valor central de U19d, el borde sigue en 0.5·F₀,** lejos del umbral. Con las extrapolaciones ∝ ℓ₁ (X ≈ 218) y lineal (X ≈ 273), el borde sería 0.85 y probablemente 0.7·F₀, en lugar de 0.85·F₀ en las dos.

## 6. Qué significa

1. **F₀ está excluido por Lyα completo** con la población medida, incluso en el caso más favorable de los probados (X ≈ 193). U19d ya lo dejaba fuera por la extrapolación (X ≥ 244). **Ahora hay dos argumentos independientes.** Para que F₀ volviera haría falta X < 193, y eso solo lo da el modelo de X₀ constante, que U19d descarta. Con F₀ quedan fuera el benchmark y la señal de cuerdas en LISA.
2. **El borde de 0.7·F₀ con X ≈ 320 resiste la verosimilitud completa,** y 0.85·F₀ no. En ese punto, la puerta de cribado y la verosimilitud completa coinciden.
3. **La puerta de cribado no es conservadora en las filas que pasan solo por forma.** Es lo contrario de lo que sugería U15c, que solo probó filas con δA ≈ 0.1, muy lejos del borde. Por eso las filas marcadas «por forma» de U17f y U17g pasan a ser dudosas (§5). Las que pasan con δA ≤ 0.32 no cambian.
4. **La mejor estimación de la banda no cambia en el X central:** F ≈ 0.06–0.5·F₀, o ≈ 0.06–0.35·F₀ con la cota de supermasivos.
   - El extremo de 0.85·F₀ que daban las extrapolaciones ∝ ℓ₁ y lineal se sostiene solo con X ≈ 220, y en el límite.
   - **La estimación más prudente sigue siendo F ≈ 0.03–0.85·F₀** en todo el intervalo de U19d, ahora con F₀ fuera por dos caminos.
5. **Una U15e,** con Lyα completo en las filas «por forma» de U17g (sobre todo 0.6·F₀ con X ≈ 486 y 0.5·F₀ con X ≈ 689), decidiría la zona de δA ≈ 0.32–0.38. **No toca el valor central,** que ya pasa con margen. Tendría que preregistrarse, y lanzarla en WSL lo decide el autor.

## 7. Lo que U15d no resuelve

- **Cuál es el X del universo.** U19d lo deja entre 244 y 1 367; la forma de la extrapolación sigue abierta.
- **Lyα completo con X > 320.** La zona dudosa de §5 necesitaría una U15e.
- **T² en la cosmología fiducial,** como en U15b. La auditoría acotó ese error en |ΔlnL| ≈ 0.4–0.7. Es del orden del margen de D1 y D3 (−2 ln p = 4.3–4.4 frente a 4), así que su exclusión es la parte más sensible a esa aproximación.
- **Las comprobaciones internas de la verosimilitud con f_nc ≈ 0.98,** cerca del borde de las simulaciones de Hooper et al.: los registros de MontePython siguen en WSL, sin revisar (nota del ESPEC).
- **Lo de siempre:** la normalización C de la isocurvatura, la forma de la población (ℓ₁ = 4.3) y la cota de supermasivos.

## Nota · 24 de octubre de 2026, 13:11 (Madrid): las comprobaciones internas de la verosimilitud

El punto de §7 sobre las comprobaciones internas con f_nc ≈ 0.98 ya tiene respuesta: **no saltaron nunca.**
- No hay ninguna línea `#Error_*` en `Lya_abgd_bin_file.txt` en ninguno de los casos de U15b, U15c y U15d (nota del 24 de octubre, 13:11, en el ESPEC).
- Ningún punto se rechazó por salirse de la rejilla de simulaciones, así que eso no puede haber empujado hacia la exclusión.
- **Sigue abierta** la aproximación de T² en la cosmología fiducial, que es del orden del margen de D1 y D3.

## Corrección · 25 de octubre de 2026, 12:22 (Madrid): los dos argumentos contra F₀ no son independientes

§0 y §6 dicen que F₀ queda excluido «por dos caminos independientes». No lo son: se complementan.
- **Lyα completo** excluye F₀ en X ≈ 193 y por encima. Por debajo no se probó.
- **La extrapolación del colapso** da X ≥ 218 en todos los modelos que ajustan.

Los dos argumentos suponen que el rendimiento blando crece con ℓ. Lo correcto es decir que **F₀ queda excluido salvo que el rendimiento crezca más despacio que cualquiera de los modelos que ajustan**: un X₀ constante, que está desfavorecido por Δχ² ≈ 19, lo devolvería por debajo de X ≈ 190.

El veredicto preregistrado no cambia. La v2 del preprint ya lo dice así.
