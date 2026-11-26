# INFORME U70 · Los agujeros de gusano de Q en la realización NPQT: PROTEGE

*Escrito el 25 de noviembre de 2026, 16:07 (Madrid).*

**Preregistro y código:**
- el preregistro es `docs/preregistros/ESPEC_U70.md` (14:23);
- el solver, `chain/u70_agujeros.py`, se escribió después del preregistro y sin reutilizar el código de la servilleta;
- el análisis, `chain/u70_analyze.py`, se subió antes de que terminara el barrido.

## En una frase

En la realización NPQT, con ℓ_QT = 300 ℓ_P, el agujero de gusano con carga de Q de menor acción tiene S_min = 9.3 × 10⁶, frente a los 221 que exige la banda. La saturación de la curvatura obliga a que toda garganta mida al menos ℓ_QT, y con el modo radial dinámico no existe ninguna solución. **Veredicto: PROTEGE**, con los límites del §5.

## Qué se calculó

- **Las soluciones.** Agujeros de gusano euclídeos O(4)-simétricos con la carga de un U(1) ⊂ SO(4) y el valle radial de 10 componentes del potencial de la acción común. Hay dos versiones: el modo radial congelado (C) y el dinámico (D), que se resuelve por tiro.
- **La gravedad.** La relatividad general y la torre NPQT de la cadena, 3h(ψ) = ρ, con ℓ_QT = 300, 1000 y 2000 ℓ_P.
- **El barrido.** Tres valores de F de la banda, 40 cargas de n = 1 a 10⁶ y la carga mínima con garganta: 954 problemas, unos 50 min con 4 procesos en la nube.
- **La acción.** Sale del lagrangiano reducido de primer orden, que es la única función analítica cuya ligadura da exactamente la ecuación de Friedmann de la cadena, menos el espacio plano con la misma frontera.

## Las puertas

| Puerta | Resultado | |
|---|---|---|
| G1: Giddings–Strominger con α = 0 | S = 3π³a₀² a 7 × 10⁻¹² | pasa |
| G2: Alvey–Escudero con α = 0, (D), n = 1 | S = 10.46, 11.38, 12.54 (≤ 6 × 10⁻⁴ de (C)); pendiente dS/d ln(1/f) = 1.000 | pasa |
| G3: la ligadura de NPQT desde F | 2.5 × 10⁻¹² | pasa |
| G4: residuo de la ligadura; tolerancias × 1/10 | 7.6 × 10⁻⁷; cambios ≤ 4 × 10⁻⁶ | pasa |
| G5: a₀ ≥ ℓ_QT en toda garganta de NPQT | ninguna violación | pasa |

**Dos arreglos de código, antes del barrido y documentados en el propio código:**
1. La primera búsqueda de la garganta tomaba la raíz más externa. Con potencial, esa raíz es la región de tipo de Sitter. Ahora toma el primer cruce hacia fuera.
2. G3 usaba diferencias finitas, con un error de 2 × 10⁻⁵. Ahora usa la derivada analítica, ∂F/∂Y = −3∫h′ dt.

Ninguno de los dos cambió una puerta después de ver resultados del barrido.

## Los resultados

| | S_min | dónde | soluciones (C) / (D) |
|---|---|---|---|
| Relatividad general | 10.5 | (D), n = 1, F = 0.5 F₀ | 117 / 109 |
| NPQT, ℓ_QT = 300 ℓ_P | **9.3 × 10⁶** | (C), n = 308, F = 0.2 F₀; a₀ = 333 ℓ_P | 74 / **0** |
| NPQT, ℓ_QT = 1000 ℓ_P | 1.0 × 10⁸ | (C) | 54 / 0 |
| NPQT, ℓ_QT = 2000 ℓ_P | 4.1 × 10⁸ | (C) | 42 / 0 |

S_req va de 208 a 221 según el borde de la banda: S_min(300) = 9.3 × 10⁶ ≥ 2 × 221, así que el veredicto es **PROTEGE**.

## Lecturas informativas

**En relatividad general, la ruptura no está suprimida.** Con el modo radial dinámico, S ≈ 10.5 n para cargas pequeñas, la escala de Alvey y Escudero con su pendiente uno en ln(1/f). La garganta mide 0.10 ℓ_P y el campo radial sube hasta r₀ ≈ M_Pl/f_a. Esta es la rama del cuscutón, que no tiene NPQT, y confirma la corrección hecha hoy en el texto de la v4.1: allí la supresión exponencial es un supuesto, no un resultado.

**En NPQT no existe solución con el modo radial dinámico.**
- Para salir de la garganta y volver al vacío, el campo radial tendría que subir hasta que el gradiente del potencial venciera al flujo, con r₀ ∼ 10⁴.
- Con r₀ tan grande, el flujo baja como 1/r₀², y la garganta tendría que encogerse por debajo de ℓ_QT, cosa que la saturación no permite.
- En todo el barrido, cada tiro sobrepasa el vacío o se queda sin garganta.

**La solución congelada tampoco es una silla consistente.** El modo radial es mucho más ligero que la curvatura de la garganta (m_r a₀ ≈ 10⁻⁴–10⁻³), así que congelarlo no está justificado.

**La lectura honesta:** en este truncamiento, con NPQT no hay ninguna silla O(4)-simétrica cargada con acción pequeña. O no hay silla, o la congelada tiene S ∼ 10⁷.

**Las cargas mínimas con garganta** en (C), con ℓ_QT = 300 ℓ_P, son n = 770, 308 y 96 para F = 0.5, 0.2 y 0.0625 F₀. Por debajo, no hay garganta.

## Lo que esta prueba no decide (§6 del ESPEC)

1. **Los operadores locales.** Si NPQT señala física gravitatoria nueva a curvaturas ∼ 1/ℓ_QT², los operadores locales que rompen O(4) podrían estar suprimidos por la escala ℓ_QT⁻¹ ≈ 10¹⁶ GeV en lugar de M_Pl, y el problema de calidad sería peor. U70 no dice nada de ellos. Es el arma de doble filo de la idea.
2. **Otras contribuciones de la gravedad cuántica,** fuera de las sillas semiclásicas O(4)-simétricas.
3. **El truncamiento.**
   - Minisuperespacio O(4), un U(1) de la SO(4), las otras direcciones de Q en su valle y acoplo mínimo.
   - La continuación euclídea de las ecuaciones reducidas de NPQT, que la cadena solo ha usado en el rebote lorentziano.
4. **La rama del cuscutón,** que no tiene NPQT. Allí vale la fila de relatividad general, sin protección.
5. **El valor de m_G.** Aunque proteja, ε sigue siendo un parámetro desnudo. U70 explica por qué los agujeros de gusano no estropearían la masa en la realización NPQT, pero no de dónde sale.

## Para el preprint (lo decide el autor)

Si se adopta, bastaría una frase en el párrafo de m_G, tras la cita de Alvey–Escudero:

> In the NPQT realisation the curvature bound forces every Euclidean throat to be wider than $\ell_{\rm QT}$; a preregistered calculation with the radial mode of $Q$ finds no charged wormhole with $S<9\times10^{6}$ for $\ell_{\rm QT}\ge300\,\ell_{\rm P}$, against $S\approx210$--$220$ needed for the band, while local operators suppressed by $\ell_{\rm QT}$, and the cuscuton realisation, are not covered.

## Ficheros

- Las filas: `results/u70/rows.json`.
- El resultado: `results/u70/U70_RESULT.json`.
- Los registros: `results/u70/run.log` y `produce.out`.
- El código: `chain/u70_agujeros.py` y `chain/u70_analyze.py`.

---

## Nota del 25 de noviembre de 2026, 20:27 (Madrid): un error del cálculo de la acción en filas de carga grande

**El error, que encontró la revisión de Codex de la PR #106.** En las filas de relatividad general con el modo radial dinámico (D) y carga grande, el tiro terminaba antes de llegar al vacío asintótico. Aun así se le sumaba una cola de la acción que supone r = 1 y un exterior casi plano, y en algunas filas esa cola era casi toda la acción.

**El arreglo.** El código marca ahora como *no resueltas* las filas en las que el término de flujo pesa ≥ 10⁻³ en a′² al final, o la cola es ≥ 1 % de S, y las deja fuera de S_min. El barrido se repitió entero. El primero está en `results/u70/_superseded/barrido1_cola_sin_validar/`, con su nota.

**Lo que cambia:**
- 41 filas quedan sin resolver, todas de relatividad general en (D) con n ≥ 1701;
- la fila de relatividad general de la tabla de resultados pasa de «117 / 109» a **«117 / 68»**, con las 41 no resueltas aparte.

**Lo que no cambia:** las puertas, los S_min (10.5; 9.3 × 10⁶; 1.0 × 10⁸; 4.1 × 10⁸), el veredicto **PROTEGE** y las cifras que cita la v4.1 (`paper_numbers_v4.py` 198/0). La escala S ≈ 10.5 n se lee en cargas pequeñas, que estaban resueltas.
