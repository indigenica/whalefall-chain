# ESPEC_U70 · ¿Suprime la saturación de NPQT los agujeros de gusano que romperían la O(4) de Q?

*Escrito el 2026-11-25 a las 14:23 (Madrid), antes de escribir el código de esta prueba.*

**Lo que ya se ha visto, y se declara:** la «servilleta 2» (24 de noviembre, código externo en `audit/servilletas_20261124/`, revisada en `docs/informes/NOTA_SERVILLETAS_20261125.md` §3).
- En relatividad general, con el modo radial de Q dinámico, da una acción S ≈ 10–13, con la escala S ∼ ln(M_Pl/f) de Alvey y Escudero ([arXiv:2009.03917](https://arxiv.org/abs/2009.03917)).
- Con NPQT y ℓ_QT = 300–2000 ℓ_P, dice que no existe garganta para la carga mínima.
- No minimizó sobre la carga, y su acción para el caso con corte UV era una estimación.

Esta prueba no reutiliza su código.

## 1. La pregunta

**Por qué importa.** El preprint supone que la gravedad cuántica rompe las simetrías globales solo con efectos exponencialmente suprimidos (§F de `main-v4.tex`). En relatividad general, los agujeros de gusano con el modo radial dinámico contradicen ese supuesto.

**La observación.** En la realización NPQT, la ecuación de Friedmann de la cadena, 3h(ψ) = ρ con h = ψ/√(1 − α²ψ²) y α = ℓ_QT², acota ψ ≤ 1/α. La garganta euclídea de un agujero de gusano O(4)-simétrico tiene a′ = 0, es decir, ψ_E = (1 − a′²)/a² = 1/a₀². Así que exige a₀ ≥ ℓ_QT, sea cual sea la materia.

**La pregunta con veredicto:** con NPQT, ¿cuál es la acción mínima S_min de un agujero de gusano con carga de Q, y basta para que la masa inducida quede por debajo de la banda?

## 2. El modelo

**La gravedad.** El lagrangiano reducido de minisuperespacio es L = N a³ [F(Y, κ) − ρ], con Y = ȧ²/(N²a²) y κ = 1/a² (k = +1).
- F es la única función analítica en Y cuya ligadura, ∂L/∂N = 0, da exactamente 3h(Y + κ) = ρ:
  F(Y, κ) = 3h(κ) − 3∫₀¹ dt t⁻² [h(κ + Y t²) − h(κ)].
- Con α = 0 se reduce a Einstein–Hilbert. El paso a la signatura euclídea es Y → −a′²/a².
- Que la continuación euclídea de las ecuaciones reducidas valga es un supuesto. Ya lo es en el nacimiento, para la versión lorentziana.

**La materia.** Q en el valle radial de 10 componentes del potencial de la acción común (λ = 320, λ_P = 192, Z_Q = 400), con la carga de un U(1) ⊂ SO(4) que gira la dirección de n en un plano.
- La carga entra como flujo conservado: es la dualización estándar del axión.
- El valle se obtiene minimizando el potencial en las otras direcciones, para cada amplitud.
- Hay dos versiones:
  - **(C)** el modo radial congelado en el vacío;
  - **(D)** el modo radial dinámico.
- El acoplo es mínimo (ξ_Q = 0). El acoplo ξ_Q = 1/4 de la acción común entra solo como lectura informativa, en relatividad general.

**La acción.** S es la acción euclídea de la solución menos la de la geometría sin agujero de gusano, con los términos de frontera de Gibbons–Hawking correspondientes, integrada hasta que el integrando cae por debajo de 10⁻¹⁰ de su máximo.

## 3. Las filas

- **ℓ_QT:** 300, 1000 y 2000 ℓ_P, más la relatividad general (α = 0).
- **F:** 0.5, 0.2 y 0.0625 F₀, es decir, el borde superior, el centro y el borde inferior de la banda.
- **La carga:** n entre 1 y 10⁶, en 40 valores espaciados logarítmicamente, más la carga mínima con garganta, si existe, refinada por bisección.
- Las versiones (C) y (D).

**S_min(ℓ_QT)** es el mínimo de S sobre n, F y las dos versiones, entre las soluciones con garganta.

**La acción necesaria.** El operador inducido es ∼ M_Pl⁴ e^{−S}, el prefactor más grande, que es la elección conservadora. Una masa no mayor que la de la banda exige S ≥ S_req = ln(M_Pl⁴/ε_req), con ε_req = ½ m_G² F². Para los bordes de la banda, S_req ≈ 208–221. El guion lo recalcula, y S_req,max es el mayor de los tres valores.

## 4. Las puertas

Son comprobaciones del código. Si una falla por un error de programación, se corrige y se documenta.

- **G1 (Giddings–Strominger):** con α = 0 y (C), S reproduce la fórmula analítica 3π³a₀², con a₀⁴ el de Giddings–Strominger, a menos de un 1 % en los tres F.
- **G2 (Alvey–Escudero):** con α = 0, (D) y n = 1, S es menor que un 10 % de la de (C), y crece con ln(1/f). La pendiente dS/d ln(1/f) entre los tres F debe caer en [0.3, 3]. Es una ventana ancha a propósito, porque el prefactor del logaritmo no lo he verificado en el artículo.
- **G3 (la ligadura de NPQT):** en 20 puntos al azar con αψ ∈ [0, 0.999], la ligadura obtenida de F coincide con 3h(ψ) a 10⁻⁸ relativo.
- **G4 (numérica):**
  - el residuo de la ligadura de Friedmann a lo largo de cada solución es < 10⁻⁶ relativo;
  - S cambia menos de un 1 % al apretar 10 veces las tolerancias, en tres filas elegidas antes de correr: la de S_min con ℓ_QT = 300, una con α = 0 en (C) y otra con α = 0 en (D).
- **G5 (la cota):** toda garganta encontrada con NPQT tiene a₀ ≥ ℓ_QT. Una violación sería un error.

Si alguna puerta falla sin un error de código que la explique, el veredicto es **INCONCLUSO**.

**Lección de U68:** las puertas se comprueban primero sobre casos analíticos (G1, G3) y con tolerancias realistas, para no acabar en INCONCLUSO por una puerta demasiado estricta.

## 5. El veredicto

Lo decide ℓ_QT = 300 ℓ_P, el valor más bajo que usa la cadena y el menos protector.
- **PROTEGE** si S_min(300) ≥ 2 S_req,max. El factor 2 cubre los prefactores O(1) de la acción.
- **NO PROTEGE** si S_min(300) < S_req,max.
- **INTERMEDIA** en los demás casos.

Si no existe garganta para ninguna carga hasta n = 10⁶, se informa la cota inferior que da la carga mayor del barrido, con su a₀ mínimo posible, y PROTEGE se decide sobre esa cota.

**Lo esperado:** PROTEGE con mucho margen. En (C), la carga mínima con garganta tiene a₀^GS ≈ 1.27 ℓ_QT, y su acción es ∼ 3π³(1.27 ℓ_QT)² ∼ 10⁷ para ℓ_QT = 300. Lo que no sé es cuánto la rebaja el modo radial dinámico.

## 6. Lo que esta prueba no decide (se declara en el informe, salga lo que salga)

1. **Los operadores locales.** Si NPQT señala física gravitatoria nueva a curvaturas ∼ 1/ℓ_QT², los operadores locales que rompen O(4) podrían estar suprimidos por esa escala, no por M_Pl, y el problema de calidad sería peor. U70 no dice nada de ellos.
2. **Otras contribuciones de la gravedad cuántica** fuera de las sillas semiclásicas O(4)-simétricas.
3. **El truncamiento.** Minisuperespacio O(4), un U(1) de la SO(4) y acoplo mínimo.
4. **La rama del cuscutón,** que no tiene NPQT. Allí vale la fila α = 0.
5. **El valor de m_G.** Aunque proteja, ε sigue siendo un parámetro desnudo.

## 7. Salidas

- El cálculo es `chain/u70_agujeros.py`, que escribe `results/u70/U70_RESULT.json`.
- El informe va en `docs/informes/INFORME_U70.md`.
- Son minutos u horas en la nube, sin redes.

---

## Nota del 2026-11-25 a las 17:55 (Madrid): el material interno

**Por qué esta nota.** Este preregistro se exporta a `whalefall-chain` porque la v4.1 cita su resultado. Las «servilletas» del 24 de noviembre y su revisión son material interno del proyecto, por decisión del autor. Por eso no se exportan:
- `audit/servilletas_20261124/`;
- `docs/informes/NOTA_SERVILLETAS_20261125.md`.

Las referencias a ellas de arriba no tienen destino en la exportación.

**Lo que se usó de ese material:** solo lo que aquí se declara como ya visto. No se reutilizó su código. Esta prueba, su código y sus resultados están completos en la exportación y se reproducen sin él.

**Lo que no cambia:** las filas, los umbrales y el veredicto.
