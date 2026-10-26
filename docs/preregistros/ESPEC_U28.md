# ESPEC_U28 · C45k: la coherencia de los rebotes del padre con χ inhomogéneo, para los hijos eternos y los directos

*25 de octubre de 2026, 21:38 (Madrid). Claude, en la nube, sobre `main` con U27 fusionado. Preregistro escrito antes de cualquier código de U28 y antes de ver el solver que se va a usar.*

*Sobre el nombre: en `audit/revisor_20261024/` aparece una «U28» que proponía un revisor (medir C de la isocurvatura). No se adoptó, y no tiene nada que ver con esta prueba.*

## 1. Pregunta

U27 dio RAMA_DIRECTA_VIABLE: el 14.9 % de los fértiles son hijos directos que aguantan los rebotes del padre, frente al 3.6 % de los eternos. Pero lo dio en el modelo de C45j, que trata el padre como polvo. En los hijos viables, el rebote está dominado por la energía cinética de χ, entre 46 y 4 200 veces el polvo. Un campo no se desincroniza como el polvo: sus gradientes tienen presión.

`docs/evidencia/nacimiento_c45/ESPEC_C45K.md` (4 de octubre) planteó esta prueba, pero nunca se construyó. U28 la retoma con la meseta de polo y con las dos ramas.

**U28 pregunta:** ¿aguanta un padre dominado por χ, con perturbaciones del campo, mejor o peor que el polvo? ¿Qué queda entonces de la fracción viable de los hijos directos y de los eternos?

## 2. Herramienta

- **El solver** es `c58c/regular_center`, del proyecto local del autor: un escalar canónico con geometría QT esférica, en coordenada areal y con centro regular. Ver `C58C_CENTRO_STATUS_20261003.md`.
  - Se copia sin cambios a `components/c58c_regular_center/`, con los ficheros de `c58c/regular_center/`: `model.py`, `evolution.py`, `diagnostics.py`, `run_bank.py`, `__init__.py`, `__main__.py`, `README.md` y `tests/test_center.py`.
  - Se copia también lo que importen de `c58c/`.
  - Los hashes se comprueban contra `data/cadena_cierre/originals_manifest.json`. Si alguno no coincide, se para y se pregunta al autor.
- **Lo nuevo,** en un módulo aparte (`chain/u28_field.py` o `components/c58c_regular_center/` sin tocar sus ficheros):
  - un fluido de polvo sin presión, con masa de Misner–Sharp y la fuente ∂_a M de C45b;
  - el potencial de polo (C = 0.3516) para χ;
  - los datos iniciales de §4.
- **Antes de escribir el código** se lee el solver. Lo que haya que concretar, como las unidades, la malla, la condición de parada o cómo se detecta el rebote, se fija en una enmienda fechada, antes de ejecutar nada.

## 3. Validaciones (puertas; si alguna falla, U28 es INVÁLIDO)

- **V0:** las pruebas propias del solver (`tests/test_center.py`) pasan sin cambios.
- **V1, el límite homogéneo:** sin perturbación, el solver con polvo y χ reproduce `components/c45b_hom/frw_core.c`.
  - Se comprueba en 5 hijos: los cuatro de los extremos de A·H en E y en D, más uno central.
  - Exige χ_land con una diferencia relativa ≤ 10⁻³ y el mismo número de rebotes.
- **V2, el límite sin χ:** con solo polvo y la perturbación de curvatura de C45d, el solver reproduce el tiempo de rebote t_B de LTB-QT de `components/c45j_ltb`, con una diferencia relativa ≤ 10⁻³. Se comprueba en los mismos 5 padres.
- **V3, la convergencia:** en 4 corridas (2 de E y 2 de D, con ε = 10⁻³ y L·H = 1), doblar la resolución cambia Δt_B en menos de un 10 % y no cambia ninguna etiqueta COHERENTE.
- **Corridas inválidas:** si se forma un horizonte aparente antes del rebote (el límite de *polar slicing* del solver), la corrida se marca INVÁLIDA y no cuenta. Si más del 20 % de las corridas de una clase son inválidas, esa clase es INCONCLUSA.

## 4. Filas

**Hijos,** con la meseta de polo, de U27 (`results/u27_direct_children.json` y su rejilla). La selección es determinista: se ordenan por ζ0_max del modelo de polvo (L·H = 1) y se toman a intervalos iguales.

| Grupo | Hijos | Papel |
|---|---|---|
| E-v | 10 de los 20 eternos viables | medida |
| E-n | 5 eternos no viables | control |
| D-v | 10 de los 83 directos viables, de ellos al menos 4 de la banda 10–50 M_Pl | medida |
| D-n | 5 directos no viables | control |

**Datos iniciales,** como en ESPEC_C45K:
- un parche cerrado en el turnaround, con los (λ, sep, σ) de cada hijo;
- dos familias de perturbación con el pico monocromático g(r) de C45d:
  - **(a)** δK/K = ε·g, solo curvatura, la misma que en el modelo de polvo;
  - **(b)** δχ̇/χ̇ = ε_χ·g, solo el campo;
- **amplitudes:** ε, ε_χ ∈ {10⁻⁵, 3×10⁻⁵, 10⁻⁴, 3×10⁻⁴, 10⁻³, 3×10⁻³, 10⁻²};
- **escala:** L·H = 1 para el veredicto, y L·H = 0.3 para los grupos E-v y D-v (informativa).

**Coste:** se fijará con el tiempo medido en V3. Si una corrida tarda de media más de 5 minutos, una enmienda puede reducir el número de hijos o de amplitudes antes de la producción, nunca después. Si las corridas son largas, la producción se lanza con `nohup` y un productor que salte las ya hechas.

## 5. Medidas y criterio por corrida

**Medidas,** en la región r ≤ 0.8·R_parche:
- **Δt_B:** la dispersión de los tiempos de rebote, con el rebote como máximo local de ψ;
- **Δχ_land/χ_land:** la dispersión del aterrizaje;
- **f_loc:** la fracción de volumen que aterriza en la misma clase (E o D) que el hijo homogéneo y por debajo de la densidad de Planck.

**Una corrida es COHERENTE** si Δt_B < 0.3·√α y f_loc ≥ 0.8. Es el criterio de ESPEC_C45K. Solo cambia que «eterno y sub-Planck» pasa a «la misma clase y sub-Planck», para que valga también para los directos.

**ε\*_χ** es la mayor amplitud de la rejilla que sigue COHERENTE, en las dos familias, y se toma la menor de las dos. Con esa misma rejilla de amplitudes, ε\*_polvo es la mayor que no supera ε_est/n_rebotes del código de C45j en ese hijo.

## 6. Veredicto

**Por clase (E y D, con los grupos -v y L·H = 1),** a partir de la mediana de ρ = ε\*_χ/ε\*_polvo:
- **MEJOR_QUE_POLVO** si la mediana es ≥ 3;
- **IGUAL** si está entre 1/3 y 3;
- **PEOR** si es < 1/3.

**La fracción viable corregida** se calcula con ζ0_max,χ = ζ0_max,polvo × ρ:
- en los hijos simulados, con su propio ρ;
- en los demás de su clase, con la mediana de ρ de la clase.

Con ella, f_v,χ(X) = viables / 556 fértiles. **Veredicto D11 para cada rama,** con los umbrales de U27:
- **VIABLE** si f_v,χ ≥ 5 %;
- **MARGINAL** entre el 1 % y el 5 %;
- **NO_COHERENTE** por debajo del 1 %.

**Controles:** si algún hijo de E-n o D-n sale COHERENTE a amplitud inflacionaria y su ρ es ≥ 3, se informa como lectura (el campo rescata padres que el polvo daba por perdidos). No cambia el veredicto.

**Qué cambiaría:**
- **Si la rama directa sale VIABLE,** D11 pasa de «estimación de polvo» a «calculada con χ». Junto con U27, eso sostiene la rama directa en la v3.
- **Si sale MARGINAL o NO_COHERENTE,** el texto de la v3 lo dice, y la rama eterna tampoco quedaría mejor si su clase sale igual.

## 7. Expectativa escrita antes de calcular

La presión de los gradientes de χ debería frenar la desincronización que en el polvo crece sin freno, como cuando la presión estabiliza un fluido frente al colapso. Por eso espero **IGUAL o MEJOR_QUE_POLVO** en las dos clases, y que la rama directa siga **VIABLE**.

Con incertidumbre real: en el rebote QT, el campo puede transferir desfases entre capas a través de los gradientes y empeorar la coherencia. No he visto el solver.

## 8. Salidas

- **Código:** `chain/u28_c45k.py` y el módulo de §2.
- **Resultados:** `results/u28_c45k.json`, con una corrida por fichero en `results/u28_runs/`.
- **Informe:** `docs/informes/INFORME_U28.md`.
- **Contrato:** U28 no cambia la cadena. La entrada D11 del contrato se actualiza con el veredicto si el autor lo lleva a la v3.

## Enmienda 1 (25 de octubre de 2026, 21:42, Madrid): umbrales como intervalos, y umbrales censurados

*Antes de cualquier código o cálculo de U28. Las señaló Codex en la PR 43. Sustituye la definición de ρ y la corrección de §5–§6. Las puertas, las filas, los criterios por corrida y los umbrales del veredicto no cambian.*

**1. Problema: la rejilla redondea hacia abajo la tolerancia del polvo.** §5 llevaba ε\*_polvo a la rejilla, redondeado hacia abajo, y §6 multiplicaba la tolerancia continua de C45j por ε\*_χ/ε\*_polvo. Las dos operaciones no se cancelan. Por ejemplo, ε_est = 0.0062 pasa a 0.003, y un campo coherente solo hasta 0.003 daría ρ = 1 aunque su umbral real fuera 0.004. Con eso, contaría como viables hijos que el campo no tolera.

**Nueva definición:**
- **El umbral del campo es un intervalo.** En cada hijo y familia, ε_lo es la mayor amplitud de la rejilla que sigue COHERENTE y ε_hi es la siguiente, la primera incoherente. Se toma la familia más desfavorable, es decir, el menor ε_lo, con su ε_hi.
- **El campo, como el polvo, se compara a un solo rebote.** El solver simula un rebote del padre, como C45d. El de varios rebotes se obtiene como en C45j: ε\*_n = ε/n_rebotes.
- **ρ se calcula con la tolerancia continua del polvo** de ese hijo, ε_est(A·H, αH², L·H), del código de C45j, sin pasar por la rejilla:
  - ρ_lo = ε_lo/ε_est;
  - ρ_hi = ε_hi/ε_est.
- **Veredicto por clase:**
  - **MEJOR_QUE_POLVO** si la mediana de ρ_lo es ≥ 3;
  - **PEOR** si la mediana de ρ_hi es ≤ 1/3;
  - **IGUAL** en cualquier otro caso, incluidos los intervalos que caen a los dos lados.
- **Fracción viable corregida:** se calcula dos veces.
  - **f_lo:** en los hijos simulados, con su propio ε_lo y ζ0_max,χ = (ε_lo/n_rebotes)/[(2/3)π²(A·H/L·H)²]; en los demás de su clase, con ε_est × la mediana de ρ_lo de la clase.
  - **f_hi:** igual, con ε_hi y la mediana de ρ_hi.
- **El veredicto D11 de cada rama usa la cota inferior f_lo,** con los umbrales de §6: VIABLE ≥ 5 %, MARGINAL entre el 1 % y el 5 %, NO_COHERENTE por debajo del 1 %. También se informa el veredicto con f_hi. Si los dos no coinciden, el informe lo dice.

**2. Problema: los umbrales censurados no tenían regla.** Puede que un hijo no sea coherente a ninguna amplitud, o que lo sea a todas.
- **Incoherente a todas las amplitudes:** se amplía la rejilla hacia abajo, en medias décadas, hasta 10⁻⁷ (3×10⁻⁶, 10⁻⁶, 3×10⁻⁷ y 10⁻⁷).
  - Si sigue incoherente, ε_lo = 0 y ε_hi = 10⁻⁷. El hijo no es viable en ninguna cota, porque ε_phys a amplitud inflacionaria es ≥ 3.6×10⁻⁵ en toda la rejilla de U27 (con el padre más pequeño, A·H = 0.34, y L·H = 1).
- **Coherente a todas las amplitudes:** se amplía la rejilla hacia arriba, hasta 3×10⁻² y 10⁻¹. Es el techo de C45d, a partir del cual la perturbación deja de ser pequeña.
  - Si sigue coherente, ε_lo = 10⁻¹ y ε_hi = ∞, con ρ_hi = ∞.
- **Las medianas** se calculan por orden con esos valores. Un valor censurado está en el lado correcto de la cota que representa, así que la mediana de ρ_lo y la de ρ_hi siguen siendo cotas válidas de la mediana verdadera.
- **Si más de la mitad de los hijos de una clase** están censurados en el mismo extremo, el informe da para esa clase la cota (por ejemplo, «ρ ≥ 3, MEJOR al menos») y el veredicto que se deduce de ella. Si no se deduce ninguno, la clase es INCONCLUSA.
- **Las ampliaciones son parte del protocolo,** no una decisión tomada al ver los datos: se aplican a todo hijo que caiga en esos casos. Las corridas extra cuentan en el presupuesto de §4.

## Enmienda 2 (25 de octubre de 2026, 21:50, Madrid): el solver preregistrado no puede cubrir el parche; producción suspendida

*Antes de cualquier código o cálculo de U28.*

**El código.** El autor envió `c58c_regular_center.zip`: 32 ficheros, cuyos hashes coinciden todos con `data/cadena_cierre/originals_manifest.json`. Contiene `regular_center`, sus dependencias reales (`c58c/matter_qt`, `c58b/src` y dos paquetes de `sidequests`) y `C58C_CENTRO_STATUS_20261003.md`. Se leyó antes de construir nada.

**Lo que es el solver.** Es una carta polar-areal, ds² = −N²f dt² + dr²/f + r²dΩ², con f = 1 − r²ψ. Evoluciona un escalar sin masa ni potencial, con centro regular y un muro reflectante exterior, en un espacio asintóticamente plano.
- **La carta exige f > 0 en todo punto,** y se detiene si no se cumple: «La carta requiere f>0; no cruza horizontes».
- **No tiene** polvo, potencial ni geometría cerrada. Se diseñó para pulsos que enfocan y salen, no para un universo cerrado en contracción.

**Comprobación de viabilidad.** No es un cálculo de U28: es el fondo homogéneo de C45, sin perturbar. En un FRW cerrado con r = a·sen χ, la carta polar-areal da f = cos²χ − sen²χ·ȧ². Las capas con |ȧ| > cot χ quedan atrapadas (f ≤ 0). Para medir |ȧ| se integró `components/c45b_hom/frw.py` en 6 de los 20 padres viables de C45j, los más favorables (sep = 10⁻⁴ y A·H = 0.51–0.62).
- El máximo de |ȧ| antes del primer rebote es 3.8–6.7.
- Por tanto, quedan atrapadas todas las capas con sen χ ≳ 0.15–0.25.
- Con separaciones más pequeñas el rebote es más violento, y la región cubierta, aún menor.

La región que mide U28 es r ≤ 0.8·R_parche, y la perturbación de C45d tiene su pico en x = sen χ ~ L·H/A·H, de orden 1. Las dos quedan fuera de la carta. Con §3, todas las corridas serían inválidas y U28 saldría INCONCLUSA por construcción.

**Lo que se decide:**
- **La producción de U28 con `c58c/regular_center` queda suspendida.** No se ha escrito código de U28 ni se ha ejecutado ninguna corrida.
- **Esto se informa como lo que es:** el solver que preregistraron ESPEC_C45K y este ESPEC no sirve para esta pregunta. No es un resultado de C45k.
- **Las vías posibles,** para que el autor decida antes de seguir:
  - **(A)** un solver nuevo en coordenadas comóviles con el polvo (lagrangianas, como LTB-QT), que no necesita f > 0. Sería la relación de masa QT, h(ψ) = 2Gm/r³, con un escalar con potencial y polvo. Las puertas V1 (el límite homogéneo frente a `frw_core.c`) y V2 (el límite de polvo frente a LTB-QT de `components/c45j_ltb`) siguen igual, y son justo las que comprobarían la derivación.
  - **(B)** la teoría lineal de perturbaciones esféricas alrededor del rebote homogéneo. El criterio de coherencia es lineal en ε (C45d), pero necesita las mismas ecuaciones QT, linealizadas.
  - **(C)** dejar D11 con χ abierto, y U27 como la estimación de polvo.
- Cualquiera de las vías exige una enmienda 3, con las ecuaciones y la malla, antes de calcular.

## Enmienda 3 (25 de octubre de 2026, 22:09, Madrid): vía A, un solver comóvil con el polvo; ecuaciones, datos, puertas y calibración

*Antes de cualquier código de U28 y antes de cualquier corrida. El autor eligió la vía A (enmienda 2). Esta enmienda sustituye §2, §3 y la parte de §5 que la contradiga. Las filas de §4, los veredictos de §6 y la enmienda 1 siguen vigentes, salvo lo que se corrige en el punto 6.*

**1. De dónde salen las ecuaciones.** El autor envió también `c45b_dn.zip`: 38 ficheros, cuyos hashes coinciden con `originals_manifest.json`. Contiene el solver en doble nulo de C45b (3 de octubre), con NPQT esférica, un escalar con potencial y masa generalizada, ya auditado (`INFORME_C45B.md`, `AUDIT_C45B.md`).
- **El solver no sirve como base.** Su gauge exige r_v > 0 en el rayo inicial, y en un universo cerrado que se contrae los rayos salientes acaban atrapados. Además, no tiene polvo, y en estos padres el polvo domina en el turnaround: la energía cinética inicial de χ es una fracción ~σ² ≲ 10⁻⁸.
- **Sus ecuaciones** (`dn/equations.py`, ec. 33–34 de Bueno, Cano, Hennigar y Murcia) se reescriben en forma covariante en 2D:
  - (E1) ∇_a∇_b r = −(κ/2)·r·m^{3/2}·(T_ab − g_ab T) + g_ab·r·ψ(3m − 2)/2;
  - (E2) R⁽²⁾/2 = □r/r + 3α²ψ³ − 6α²ψ m J − κ m^{3/2} p_⊥.
  
  Aquí m = 1 − α²ψ², T es la traza 2D, J = b² − κ²m[T_abT^{ab}/2 − T²/4]/4 y b = −κT√m/4 − 3ψ/2.
- **Comprobación antes de escribir el solver:** en el límite homogéneo, (E1) reproduce la ecuación de `frw_core.c`, y (E2) es idéntica a ella con presión arbitraria. El residuo es ≤ 4×10⁻¹⁶ en 5 estados al azar.

**2. Gauge y variables.** El polvo sigue geodésicas irrotacionales, así que se usa la gauge síncrona comóvil, ds² = −dt² + e^{2Λ}dχ² + r²dΩ², en todo el universo cerrado, χ ∈ [0, π].
- **Malla:** centrada en celdas, sin nodos en los polos, con paridad en los dos polos (r y U impares; Λ, L, φ y Π pares).
- **Variables:** r, U = ṙ, Λ, L = Λ̇, φ, Π = φ̇ e Y = r³h(ψ), con κY = 2·(8πG)·M_gen/2.
- **Evolución:**
  - ṙ = U, con U̇ = −(κ/2) r m^{3/2} p_r − rψ(3m − 2)/2 (la componente tt de E1);
  - L̇ = −L² + □r/r + 3α²ψ³ − 6α²ψ m J − κ m^{3/2} p_⊥ (E2), con □r = −U̇ − LU + e^{−2Λ}(r'' − Λ'r');
  - Π̇ = −(L + 2U/r)Π + e^{−Λ}r^{−2}(e^{−Λ}r²φ')' − V'(φ);
  - Ẏ = −κr²(p_r U − e^{−2Λ}Πφ'r'), la ley de flujo de la masa generalizada, como en `dn`;
  - ψ = ψ(Y/r³), que cumple |αψ| < 1 sin recortes.
- **El polvo** es exacto: ρ_d = D(χ)/(e^Λ r²), con D fijo.
- **Ligaduras vigiladas,** no impuestas:
  - la hamiltoniana, r²ψ − (1 + U² − e^{−2Λ}r'²);
  - la de momento (la componente tR de E1), U' − Lr' + (κ/2) r m^{3/2} Πφ'.
- **Integración:** RK4 con paso CFL en tiempo propio, dt = c·Δχ·min e^Λ, con c = 0.25.

**3. Datos iniciales:** en la rodaja del rebote anterior del padre, U = 0 y L = 0.
- **El fondo** es el padre homogéneo de C45 de cada hijo (`c45_state`, con el polo, C = 0.3516), integrado hacia atrás desde el turnaround hasta su rebote anterior.
- **Es la sincronización de C45d y C45j,** «tiempo de rebote previo uniforme»: cada capa parte de su rebote a la vez.
- **La ligadura de momento** se cumple exactamente con φ uniforme en la rodaja.
- **Familia (a):** la curvatura de C45j, K(x) = 1 + ε·sinc(x/(L·H/A·H)), con x = sen χ y una masa total Y ∝ x³, como M = μx³ en C45j. Se atenúa suavemente a cero entre x = 0.95 y el ecuador.
  - El atenuador es un C∞ tipo «bump», y se fija en el código antes de ejecutar.
  - Como es función de x, es simétrica respecto al ecuador y la geometría cierra sin singularidad.
  - r sale de la hamiltoniana con U = 0, r²ψ(Y/r³) = Kx². El polvo es el resto, ρ_d = ρ_total − Π²/2 − V.
  - Una corrida con ρ_d < 0 en algún punto es INVÁLIDA.
- **Familia (b):** δΠ/Π = ε_χ·g(x) con el mismo g, la misma atenuación y K = 1, y el polvo sin perturbar. Y se integra con la ley de masa en χ.

**4. Puertas** (sustituyen a las de §3):
- **V0, autoconsistencia:** las ligaduras convergen a orden ≈ 2 al doblar N.
- **V1, límite homogéneo:** sin perturbación, el solver reproduce `frw_core.c` en 5 padres, con t_B y a_min relativos ≤ 10⁻³ y la misma rama.
- **V2, límite de polvo:** sin χ y con los mismos (μ, K, α) que `ltb/growing_qt.build_qt` de C45j, el t_B de cada capa con x ≤ x_hi coincide con 2τ_B de `ltb/shells.Shell`, con diferencia relativa ≤ 10⁻³. Se comprueba en 3 padres y ε = 10⁻⁴.
- **V3, convergencia:** como en §3. Doblar N cambia Δt_B menos de un 10 % y ninguna etiqueta.
- **V4, independiente:** en el límite GR (α → 0) con solo polvo, el colapso homogéneo reproduce Oppenheimer–Snyder, τ = (π/2)√(R0³/2GM), a 10⁻³.
- **Corrida inválida:** si una ligadura relativa supera 10⁻² antes del rebote, o si aparece cruce de capas (r' cambia de signo fuera del ecuador). Mantiene la regla del 20 % de §3.

**5. La medida:**
- **Δt_B:** la dispersión de los tiempos de rebote (U pasa de − a +) en las capas con x ≤ x_hi, la misma región que `eps_est` de C45j.
- **f_loc:** la fracción de volumen de esa región que aterriza en la clase del hijo homogéneo y por debajo de la densidad de Planck.
  - Se mide con universos separados: tras el rebote de la última capa, cada capa se continúa como un FRW homogéneo con su estado local (H = U/r, ψ, χ, χ̇ y polvo), con la misma integración que `frw_core.c` pero con H ≠ 0 inicial.
  - Se continúa hasta 60 e-folds o un nuevo colapso.

**6. Calibración frente al polvo.** Corrige la enmienda 1 por un factor de criterio que no tenía en cuenta.
- **El problema:** ε_est de C45j es la amplitud a la que Δt_B = √α. U28 llama COHERENTE a Δt_B < 0.3·√α. En régimen lineal, el umbral del polvo con el criterio de U28 es 0.3·ε_est.
- **Por tanto:**
  - ρ_lo = ε_lo/(0.3·ε_est);
  - ρ_hi = ε_hi/(0.3·ε_est).
- **La tolerancia corregida** que entra en f_lo y f_hi es ε_est·ρ (la de C45j multiplicada por el efecto del campo), dividida por n_rebotes y comparada con ε_phys(ζ0), como en C45j.
- **El resto de la enmienda 1 no cambia:** los intervalos, la censura, las medianas y el veredicto con f_lo.

**7. Expectativa (§7):** no cambia.

**8. Coste:** pruebas y desarrollo, del orden de horas. Cada corrida es 1D con N = 256–512, y su duración se medirá en V3. Si hace falta recortar, se hará con una enmienda antes de la producción, como en §4.

## Enmienda 4 (25 de octubre de 2026, 22:11, Madrid): reparto adiabático en la familia (a), y dos precisiones

*Antes de cualquier código o corrida de U28.*

1. **La familia (a) es adiabática.** La enmienda 3 hacía que el polvo absorbiera toda la perturbación de densidad (ρ_d = ρ_total − Π²/2 − V). En la rodaja del rebote, la energía cinética de χ es entre 46 y 4 200 veces el polvo (C45j). Con ese reparto, ρ_d sería negativo a partir de ε ~ 10⁻⁴–10⁻², y casi todas las corridas útiles serían inválidas.
   - **El reparto nuevo:** el campo y el polvo se reparten la densidad perturbada en la proporción del fondo. Con f_k = (Π_b²/2)/(ρ̄ − V(φ_b)), queda Π²/2 = f_k·(ρ − V(φ_b)) y ρ_d = (1 − f_k)·(ρ − V(φ_b)).
   - **φ sigue uniforme** en la rodaja, así que la ligadura de momento se cumple exactamente.
   - **El límite de solo polvo (V2) no cambia,** porque ahí f_k = 0.
   - La regla de invalidez por ρ_d < 0 se mantiene, pero ya no debería activarse.
2. **Precisión numérica.** m se calcula siempre como 1/(1 + α²y²), con y = Y/r³, y nunca como 1 − α²ψ². En los padres viables, 1 − αψ llega a 10⁻¹⁰ en el rebote (`z_max` de C45h), y la resta perdería toda la precisión.
3. **Errata de la enmienda 3, punto 2.** Donde dice «κY = 2·(8πG)·M_gen/2», debe decir: Y = r³h(ψ) = 2G·M_gen = κ·M_gen/(4π), con M_gen la masa generalizada de `dn`.

## Enmienda 5 (25 de octubre de 2026, 22:15, Madrid): formulación de Misner–Sharp en lugar de evolucionar Λ con E2

*Durante el desarrollo de la puerta V4, antes de cualquier corrida de U28 con perturbaciones o con χ.*

- **Lo que pasó:**
  - Con el esquema de la enmienda 3, el colapso homogéneo en relatividad general (V4) reproduce la cicloide de Oppenheimer–Snyder a 3×10⁻¹¹ con N = 64.
  - Pero las ligaduras crecen hasta O(1), y con N = 128 la corrida revienta a mitad del colapso.
- **La causa:** la ecuación E2 tomada como evolución de Λ contiene un término Λ_tt ∝ −e^{−2Λ}(r'/r)·Λ'. Es un problema mal planteado en el sentido de Hadamard: los modos de número de onda k crecen a un ritmo ∝ √k, así que cuanto más fina es la malla, antes revienta.
- **El cambio:**
  - Se evoluciona Γ = e^{−Λ}r' con la componente tR de E1, Γ̇ = −(κ/2)·r·m^{3/2}·Π·q, con q = e^{−Λ}φ'.
  - e^Λ se obtiene como r'/Γ.
  - Klein–Gordon se escribe en forma conservativa, con P = e^Λ r²Π: Ṗ = (e^{−Λ}r²φ')' − e^Λ r² V'(φ).
  - Con esto, la ligadura de momento se cumple de forma idéntica, y la vigilada pasa a ser la hamiltoniana.
  - E2 queda como diagnóstico, no como ecuación de evolución.
- **El dominio pasa a ser χ ∈ [0, π/2],** con simetría especular en el ecuador: r, U, φ, P e Y son pares, y Γ es impar.
  - Los datos de las dos familias son funciones de x = sen χ, así que son especulares.
  - En el ecuador, r' y Γ se anulan a la vez; las celdas están centradas y ninguna cae justo en el ecuador, de modo que el cociente r'/Γ queda bien definido.
- **Los datos iniciales** usan las mismas derivadas discretas que la evolución (e^Λ = r'/Γ y ρ = Y'/(κr²r')). Así, la ligadura hamiltoniana se cumple exactamente en t = 0.
- **Lo demás no cambia:** las puertas V0–V4, las filas, la medida y el veredicto.

## Enmienda 6 (25 de octubre de 2026, 23:16, Madrid): cierre sin atenuador, reglas de umbral y de clasificación, y lo ya ejecutado en desarrollo

*Antes de las puertas registradas y de la producción.*

**Lo ejecutado en desarrollo, que no entra en el veredicto.** Se declara aquí, y la producción lo repite todo.
- **V4:** pasa a 10⁻¹² con N = 64, 128 y 256.
- **V2:** pasa en 3 padres de C45j, con ≤ 2.3×10⁻⁷ capa a capa.
- **V1:** pasa en 5 padres, con t_B ≤ 4×10⁻⁵ y a_rebote ≤ 1.1×10⁻⁵.
- **Dos corridas perturbadas,** familia (a), L·H = 1 y ε = 10⁻⁴:
  - el primer hijo E-v da Δt_B/√α = 0.0576, frente a 0.0575 del modelo de polvo, y f_loc = 0.88;
  - el primer hijo D-v da una caústica (punto 1);
  - además, el clasificador con ε = 0 da f_loc = 1 en 3 hijos y reproduce el χ_land de la rejilla de polo (por ejemplo, 2627.81 frente a 2627.83).
- **La fracción cinética del campo en la rodaja del rebote** es f_k = 2×10⁻⁴ en el primer E-v y 0.049 en el primer D-v. En estos padres el polvo domina en el rebote; el informe dará la distribución de f_k.

**1. Sin atenuador; cierre en la «panza».**
- **El problema:** el atenuador de la enmienda 3 llevaba la perturbación de curvatura a cero en x ∈ [0.95, 1]. Eso exige un dK/dx unas 30 veces mayor que el del sinc, y justo ahí r' ∝ cos χ ya es pequeño. En el primer hijo D-v con ε = 10⁻⁴, eso produce un cruce de capas en x ≈ 0.977 antes del rebote, en el mismo sitio y con la misma profundidad con N = 128, 256 y 512.
- **Diagnóstico:** es un artefacto de la construcción, no del perfil de C45j.
- **El cambio:** la etiqueta de capa pasa a ser x = x_e·sen χ, con K(x_e)·x_e² = 1.
  - Es el máximo del radio areal en la rodaja, la «panza» del universo cerrado.
  - En ella Γ = 0, y la geometría cierra con simetría especular sin modificar K.
  - El perfil de C45j se aplica entero hasta x_e ≈ 1 − εg/2, y las capas con x ≤ 0.95 son las mismas de C45j.
  - La familia (b) tampoco lleva atenuador: con K = 1, ya cierra.
- **V2 se repite en las puertas registradas** con esta construcción.

**2. Clasificación por universos separados (f_loc).**
- Se toman 24 celdas de la región x ≤ 0.95, equiespaciadas en índice. Cada una pesa el volumen e^Λ r² dχ de las celdas hasta el punto medio con sus vecinas.
- **Cada celda se clasifica** al final de la corrida, cuando todas han rebotado, con la integración de `frw_core.c` en H ≠ 0:
  - **fértil** si llega a 60 e-folds desde a₀ y max|χ| > 5;
  - **sub-Planck** si max(ρ en la corrida, E_max de la continuación)/c₂ < 1;
  - su clase es E si |χ_fin| > χ_et = 540.95, y D si no.
- Una celda con ψ_loc ≤ H² (localmente abierta) cuenta como de otra clase.

**3. Umbrales con corridas inválidas.**
- ε_lo es la mayor amplitud válida y COHERENTE.
- ε_hi es la menor amplitud válida e incoherente por encima de ε_lo. Si no hay ninguna, es la siguiente amplitud de la rejilla, que es la cota más estrecha.
- Si el conjunto coherente no es monótono, se informa, y ε_lo sigue siendo la mayor coherente.
- La ampliación de la rejilla de la enmienda 1 se aplica si todas las amplitudes válidas son coherentes, o si ninguna lo es.

**4. Resolución.** N = 256 en producción si V3 pasa con N = 256 frente a 512. Si no pasa, se hará otra enmienda antes de producir.

## Enmienda 7 (25 de octubre de 2026, 23:50, Madrid): error de código en Klein–Gordon cerca del polo, corregido antes de producir

*Antes de la producción. Las puertas registradas se repiten enteras con el código corregido.*

**El síntoma.** En los hijos donde el campo domina en el rebote (f_k = 0.97–0.996), las corridas perturbadas (ε = 10⁻³) se atascaban en t/T ≈ 0.98. e^Λ caía a cero en la cuarta celda desde el polo, con N = 256 y con N = 512: siempre la misma celda, no el mismo x, así que el fallo era numérico. Ya en t/T = 0.9, las celdas junto al polo tenían Y/r³ un 30–70 % por debajo de κρ/3 y una U/r del centro de 0.3 veces el fondo.

**La causa.** La divergencia de Klein–Gordon se dividía por e^Λ·r²·h, no por el volumen exacto de la celda esférica.
- En la celda del polo, eso es un factor 4/3; en las siguientes, un error relativo de orden h²/r².
- Cuando el campo domina, esa fuente espuria se concentra en el centro.

**El arreglo.** Klein–Gordon pasa a volúmenes finitos con el volumen exacto de cada celda, Vol_i = (r³_{i+½} − r³_{i−½})/(3Γ_i), donde d(r³)/(3Γ) = e^Λ r² dχ. Se usa P = Π·Vol/h como variable conservada. `regular_center` hace lo mismo.

**Comprobación.** En el mismo hijo y con ε = 10⁻³, las desviaciones junto al polo pasan a ~10⁻⁴, de acuerdo con una perturbación suave de 10⁻³.

**Lo demás no cambia:** ni las ecuaciones, ni los datos, ni las filas, ni los criterios. Las corridas de V1 a V3 hechas antes de esta enmienda quedan superadas, y las puertas se vuelven a ejecutar.
