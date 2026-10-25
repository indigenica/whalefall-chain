# U17f · El borde superior de la banda de F con el colapso extrapolado de U19b (X del universo ≈ 170–320) · Preregistro

*19 de octubre de 2026, 16:45 (Madrid). Escrito antes de calcular ninguna fila de U17f; la hora que cuenta es la de su commit. Se ejecuta en la nube, con la caché de CLASS de U12, U21 y U17e.*

**Lo que ya se ha visto, y se declara:**
- **Todas las filas de U17e** (`lya/results/u17e/U17E_RESULT.json`):
  - filas Y (Y = 1.70, X ≈ 275 en F₀): Lyα falla en F₀ (δA = 0.47) y pasa en 0.7·F₀ (δA = 0.28);
  - filas Y+ (Y = 2.24): pasan Lyα en 0.7·F₀ y por debajo;
  - filas X (K = 1.11): pasan en toda la rejilla.
- **U20b, U23 y la banda combinada** (`INFORME_U17e_U20b_U23.md`): F ≈ 0.09–0.7·F₀, o 0.09–0.5·F₀ con la cota secundaria de supermasivos.
- **U19b completo** (12 redes):
  - el X del universo (B2, ℓ ≈ 90) cae entre ≈ 170 y ≈ 320 con las extrapolaciones que crecen con ℓ (§7, corregido);
  - la lectura (ii) del control U(1), calibrada con un U(1) de Y ≈ 0.94, da Y_RP³ ≈ 1.0 ± 0.5 (§8).
- **No se ha calculado** ninguna masa, ninguna longitud de libre recorrido ni ninguna fila de CLASS con los Y de abajo, ni con F = 0.85·F₀.

## Motivo

- La corrección de §7 de U19b deja el X del universo entre ≈ 170 y ≈ 320.
- U17e pasó por CLASS un solo valor, X ≈ 275 (Y = 1.70), y una rejilla de F con un salto entre 0.7·F₀ y F₀.
- La enmienda 5 de U19b y su informe (§7 y §8) piden una U17f preregistrada antes de usar esas extrapolaciones en filas de CLASS.

**U17f pregunta dos cosas:**
1. ¿Se mueve el borde superior de la banda (0.7·F₀) en todo el intervalo de X?
2. ¿Vuelve F₀, y con él LISA, en el extremo favorable del intervalo?

Para no elegir una extrapolación después de ver los datos, **el veredicto se da sobre todo el intervalo**, no sobre una de sus extrapolaciones.

## Filas (`lya/u17f_class.py`)

**La población** es la de U17e, sin cambios: la forma medida `measured_shape("low")` de las cuatro redes primarias de U19b (x̄ = 5.83), con la rama viva B2, q = 0.87, x₀ = 10 y la emisión de U12.

**El número** es el de las filas Y de U17e, K = (8π/27)·R_ρ·ξ_B2·ℓ·Y/x̄, con R_ρ = 1.56 y m_G autoconsistente (`u17_ksw.solve` con `c_eta = x̄/Y`). Es decir, X ∝ ξℓ sobre la rejilla de F, y en F₀ vale X ≈ 162·Y.

**Familias primarias, el intervalo de U19b §7:**

| Y | 1.05 | 1.20 | 1.36 | 1.70 | 1.98 |
|---|---|---|---|---|---|
| X en F₀ | ≈ 170 | ≈ 195 (∝ ℓ + cuerdas) | ≈ 221 (lineal en ℓ₁ + cuerdas) | ≈ 275 (U17; = U17e Y) | ≈ 320 (potencia ajustada + cuerdas) |

**Familias informativas, sin veredicto:** la lectura (ii) de U19b §8, con Y = 1.0 (X ≈ 162) e Y = 0.5 (X ≈ 81, su −1σ).

**Lo que no se modela:** la dependencia con ℓ de cada extrapolación sobre la rejilla de F. Todas las familias escalan como ξℓ, igual que las filas Y de U17e. Entre F₀ y 0.5·F₀, K cambia un 4 %, así que no importa para el borde superior; para el inferior, se dice en el informe.

**Rejilla de F:** F/F₀ ∈ {1, **0.85**, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03}. Es la de U21 y U17e, más 0.85 para resolver el borde.

## Qué filas pasan por CLASS

1. **F ∈ {1, 0.85, 0.7, 0.5}: todas, en las siete familias.** Si alguna cumple la regla de frío trivial, se calcula igual.
2. **F ≤ 0.35, por dominancia:** una fila con Y ≤ 2.24 pasa Lyα si la fila Y+ de U17e con el mismo F la pasó, y las once lo hicieron. El argumento: con menos cuantos del colapso, m_G sube, la oscilación empieza antes y la fracción no fría baja.
   - **Control:** la fila Y = 1.98, F = 0.35 se calcula también con CLASS y debe pasar la puerta primaria.
   - **Si el control falla,** se calculan con CLASS todas las filas F ≤ 0.35 que no cumplan la regla de frío trivial de U21.
3. **Regla de frío trivial de U21:** λ_proxy/λ_WDM(5.3 keV) < 0.05. Se informa en todas las filas; las de F ≤ 0.5 que la cumplen pasan Lyα sin CLASS.

Son 28 filas de CLASS, más el control.

## Puertas por fila (las de la banda combinada de U17e, U20b y U23)

1. **Lyα:** la puerta primaria de U12, sin cambios (`U5.gates`).
2. **Isocurvatura, con la cota observacional de U20b:** f_iso = C·(k_CMB/k_*)³/A_s < 0.0084, con C = 1 y r = r_N = 2.40.
   - Se informan también r = r_W = 1.78 y el criterio de U20.
3. **Radiación oscura:** ΔN_eff(emisión) ≤ 0.3.
4. **Cuerdas en el CMB:** F ≤ 1.26×10¹⁵ GeV.
5. **Superradiancia estelar (U23, puerta primaria):** la fila no puede estar «excluida» por Witte y Mummery 2025 con f_eff ∈ [F/4, F].
   - Se informan sin puerta: «frontera», «probable», la cota secundaria de supermasivos (m_G < 1.7×10⁻¹⁷ eV) y LISA (f_a = F/2 ≥ 5×10¹⁴ GeV).

## Veredicto (familias primarias)

Para cada familia, F_max es el F más alto de la rejilla que pasa las cinco puertas. Sean:
- **F_max,bajo:** el menor F_max de las cinco familias primarias (en principio, Y = 1.98);
- **F_max,alto:** el mayor (en principio, Y = 1.05).

| Veredicto | Condición |
|---|---|
| **BENCHMARK_PERMITIDO** | F_max,bajo = 1 |
| **BORDE_SUBE** | F_max,bajo = 0.85 |
| **BORDE_0.7_CONFIRMADO** | F_max,bajo = 0.7 |
| **BORDE_BAJA** | F_max,bajo ≤ 0.5 |
| **SIN_BANDA** | alguna familia primaria no pasa en ningún F |

- **Calificador:** si F_max,alto = 1 > F_max,bajo, se añade «; F₀ PERMITIDO CON X ≤ X*», con X* el de la familia más alta que deja F₀. Es decir, el benchmark depende de la extrapolación.
- **Se informan también:**
  - la banda completa de cada familia (F_min, F_max, contigüidad) y cómo la parte la superradiancia estelar;
  - la banda con la cota secundaria de supermasivos;
  - el δA de cada fila;
  - las familias informativas (ii).

## Lo que se espera, escrito antes de calcular

- **El borde superior:** F₀ falla en Y = 1.70 (ya visto) y, por monotonía, en Y = 1.98. En Y = 1.05, K baja de ≈ 47 a ≈ 29 en F₀.
  - **Visto al escribir este ESPEC:** U21 P, con Y ≈ 1 y la log-normal, pasaba Lyα en F₀ por forma (δA = 0.32), aunque no por δA estricto. La cortaba su propia isocurvatura, que U20b retiró después.
  - Con la forma medida, más caliente, en Y = 1.70 falla también la forma (δA = 0.47). En Y = 1.05 es incierto.
  - **Expectativa: BORDE_0.7_CONFIRMADO o BORDE_SUBE.** Si en el extremo favorable F₀ pasa por forma, con el calificador.
- **El borde inferior:** con menos K, m_G sube. Ω_mis ∝ m^½·(1 + K), así que en Y = 1.05 la masa es ≈ 2.6 veces mayor. La ventana estelar se desplaza hacia F mayores, un factor ≈ 1.3 en F. **Se espera que el borde inferior suba de 0.09 a ≈ 0.125·F₀ en el extremo favorable.**
- **Supermasivos:** con m_G mayor, la cota secundaria (m_G < 1.7×10⁻¹⁷ eV) recorta menos. Su borde de 0.5·F₀ podría subir a 0.7·F₀ en Y ≈ 1.

## Lo que U17f no resuelve

- **Cuál de las extrapolaciones es la buena:** U17f solo acota el intervalo. Lo decidiría una tercera masa de U19b (ℓ_G ≈ 4.0), con su propio preregistro.
- **La normalización C de la isocurvatura,** que sigue siendo 1.
- **La forma de la población,** que sigue siendo la de ℓ₁ = 4.3.
- **La superradiancia con autointeracción para supermasivos,** que sigue siendo secundaria.

## Enmienda 1 · 19 de octubre de 2026, 16:47 (Madrid): la meseta y un control de la tubería, antes de calcular ninguna fila

- **La meseta:** U17f se ejecuta con `WHALEFALL_PLATEAU=sum`, igual que U17e, U21 y los objetivos de `lya/` del `Makefile`.
  - U26 mostró que el polo mueve m_G ≤ 0.4 % y arrastra los veredictos de CLASS.
  - Con la suma, la dominancia sobre las filas Y+ de U17e compara filas calculadas igual.
  - Visto al escribir esta enmienda: F₀ vale 1.0701×10¹⁵ GeV con la suma y 1.0691×10¹⁵ GeV con el polo. No se ha calculado ninguna fila.
- **Control de la tubería, sin veredicto:** la familia Y = 1.70 es la de las filas Y de U17e. En F ∈ {1, 0.7, 0.5}, su δA tiene que coincidir con el de U17e a 10⁻³. Si no coincide, se busca el error antes de dar el veredicto, y se informa.

## Nota · 19 de octubre de 2026, 22:44 (Madrid): un error al escribir el resultado, corregido sin tocar ningún criterio

- **El fallo:** la primera pasada completa (29 filas de CLASS) falló al escribir `U17F_RESULT.json`. `u23_superradiance.f_cd` devuelve `inf` cuando el modo 211 del supermasivo no es superradiante, y el JSON se escribía en modo estricto.
- **La corrección:** el código de salida escribe `inf` como `null` (`_finite`). Se volvió a ejecutar desde la caché de CLASS.
- **Lo que no cambia:** ninguna fila, ninguna puerta ni la lógica del veredicto.
- **Los dos controles pasan:** la tubería (Y = 1.70 reproduce a U17e exactamente) y la dominancia (δA = 0.090).
- **El resultado** está en `docs/informes/INFORME_U17f.md`.
