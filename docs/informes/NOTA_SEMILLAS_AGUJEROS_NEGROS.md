# NOTA · Los agujeros negros como fábricas de semillas: órdenes de magnitud

*27 de octubre de 2026, 15:01 (Madrid). Claude, en la nube, a petición del autor.*

- **Qué es:** una exploración sin veredicto ni preregistro. Recoge condiciones necesarias y estimaciones dimensionales para el supuesto abierto «Seed» del preprint (§8 de `main-v3.tex`), y no cambia la cadena, el contrato ni el preprint.
- **Código:** `chain/explora_semillas_bh.py`. **Resultado:** `results/explora_semillas_bh.json`.

## 1. Qué es una semilla en whalefall

Así son los padres de los 30 hijos de U28 en su punto de máxima expansión (el estado de C45 con la meseta de polo):

| Magnitud | Rango |
|---|---|
| χ₀ | **0 en los 30**: el inflatón está en el mínimo de su potencial |
| energía cinética de χ / polvo | 5.0×10⁻¹¹ – 4.7×10⁻² |
| A·H (radio de curvatura en radios de Hubble de la meseta) | 0.42 – 21 |
| densidad del polvo / V₀ | 0.002 – 5.7 |

**Una semilla no es un parche de inflatón en la meseta.** Es una bola de polvo cerrada:
- con un tamaño de unos pocos radios de Hubble de la inflación;
- con una densidad del orden de la de la meseta;
- con χ en reposo en el mínimo, salvo una pequeña velocidad.

La antifricción de la contracción sube χ a la meseta, y el latigazo del rebote, Δχ = (π/2)·p_χ/α, lo lleva hasta χ_land. Sin p_χ no hay latigazo: el hijo es estéril.

## 2. Tres maneras de producir semillas en un universo como el nuestro

**(a) Regiones que se recolapsan.**
- Una región cerrada de polvo con Λ se recolapsa si y solo si Λ ≤ 4/(9A²) (cálculo exacto con sympy, 8πG = 1, A = ρ_d a³/3). Es decir, si en su máximo su densidad de polvo es al menos 2ρ_Λ. Nuestra Λ permite que se recolapsen los cúmulos.
- Pero al darse la vuelta su densidad es ~10¹¹⁰ veces menor que la de una semilla.
- Y χ está muerto: el inflatón se desintegró en el recalentamiento, con Γ_χ ~ 6.5 GeV y Γ_χ t₀ ~ 4×10⁴². Su energía residual es ~e^(−10⁴²).
- Además choca con el no-go del colapso local dentro de un padre que se expande. **Descartada.**

**(b) Recolapso global del universo** (energía oscura que acaba negativa).
- Para llegar del radio observable actual a ℓ_QT ~ 10³ ℓ_P, el universo tendría que contraerse un factor ~3×10⁵⁸.
- La cizalla crece frente al polvo como a⁻³, es decir ~10¹⁷⁵ veces. Ningún rebote coherente lo resiste, salvo con una fase ecpirótica (w ≫ 1), como en el modelo cíclico, y χ seguiría muerto. **Descartada** dentro de este modelo.

**(c) Interiores de agujeros negros.**
- Es la única de las tres que no queda descartada por órdenes de magnitud.
- La materia que colapsa pasa por densidades del orden de la meseta, que es el régimen de las semillas.
- La única fuente de p_χ sería cuántica. Se compara una estimación dimensional, f_q ~ (H⁴/8π²)/ρ_d, con la fracción cinética de cada padre. Se evalúa en unidades físicas: los estados de C45 están reescalados por c₂ = V₀/V₀_bench, como en el criterio sub-Planck de U28, y f_q no es invariante bajo ese reescalado, a diferencia de la fracción cinética.
  - **En los 20 padres viables, f_q se queda corta** entre 10²·⁴ y 10¹¹·⁰ veces (mediana, 10⁷·⁹).
  - **En ninguno de los 30 padres alcanza la fracción necesaria:** al más cercano le falta un factor 10⁰·⁷.

  *Corrección:* en la conversación del 26 de octubre se dijo que esta estimación caía «en la parte baja del rango de los padres fértiles». Comparada padre a padre, no es así.

  *Nota de la misma tarde, antes de fusionar:* la primera versión de este cálculo usaba las densidades reescaladas del solver. Daba una carencia mediana de 10⁵·⁵ y dos padres no viables que sí llegaban. Codex lo señaló en la PR 47; corregido, los números son los de arriba.

## 3. Lo que tendría que pasar para que (c) funcionara

1. **Un p_χ suficiente.** Por (c), las fluctuaciones cuánticas a la escala de la meseta no bastan en ninguno de los 30 padres. Haría falta otra fuente: producción de partículas durante el colapso, un acoplamiento de χ con la materia que colapsa, o un campo distinto que haga de χ en el hijo.
2. **Un núcleo homogéneo de unos pocos radios de Hubble de la meseta.** El interior de un agujero negro estelar es muchísimo mayor cuando llega a esa densidad.
3. **Aislarse del exterior.** El horizonte separa causalmente el interior del universo en expansión, pero no homogeneiza nada: dentro de una superficie atrapada las capas siguen siendo distintas si lo eran al formarse.

## 4. Sobre la propuesta de una prueba U29 (revisión externa del 27 de octubre)

Una revisión externa propone evolucionar un perfil de sombrero de copa suavizado, ρ(r) = ρ₀[1 − tanh((r − r₀)/Δr)]/2, con el solver de U28, y comprobar que el núcleo central rebota en fase dentro del horizonte aparente. Antes de preregistrarla hay que ajustar tres cosas:
- **El núcleo homogéneo rebota en fase por ser homogéneo, no por el horizonte.** Las capas de un núcleo LTB con K y M/R³ iguales llegan juntas al rebote con o sin superficie atrapada. Lo que habría que medir es hasta dónde llega la perturbación del borde, por las ondas de rarefacción y los gradientes de χ, antes del rebote del centro.
- **El solver de U28 no sirve tal cual.** Evoluciona un universo cerrado entero (χ ∈ [0, π/2], con espejo en el ecuador), no un parche empalmado con un exterior asintótico. Una U29 necesitaría un dominio nuevo, con la frontera exterior y su condición de contorno.
- **No resolvería el punto 1.** Aunque el núcleo rebotara en fase, sin p_χ el hijo sería estéril. Una U29 útil tendría que incluir también el origen de p_χ.

## 5. Lectura

- La vía de los agujeros negros cerraría el ciclo de whalefall de forma concreta: universo, estrellas, agujeros negros, rebote NPQT en el interior y un hijo.
- Con lo que hay hoy, le faltan dos piezas: un mecanismo para p_χ y un núcleo de tamaño de semilla.
- La Λ actual no interviene en la producción de semillas. Como mucho, la vía de los agujeros negros daría una cota de tipo Weinberg: si Λ es demasiado grande no hay estrellas, sin estrellas no hay agujeros negros y sin agujeros negros no hay hijos.
- **Recomendación:** mantener «Seed» como supuesto abierto. En una próxima versión, añadir una línea que diga qué tendría que cumplir un interior de agujero negro (p_χ, tamaño y homogeneidad del núcleo).
