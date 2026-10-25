# U15b · La cola tibia frente a Planck + BAO + Lyman-ABGD, con cadenas MCMC · Preregistro

*12 de octubre de 2026, 10:45 (Madrid). Sustituye al método de U15, cuyo arnés a cosmología fija no pasó su validación (`docs/informes/INFORME_U15.md`). Se ejecuta en WSL.*

**Lo que se sabe antes.** Lo de U15: con la cosmología fija y perfilando los parámetros astrofísicos, MIKE/HIRES no distingue la WDM de 2 keV de ΛCDM. No se ha evaluado ninguna fila del modelo con ninguna versión de la verosimilitud.

## Idea

Se reproduce el procedimiento de Hooper et al. (2022): una cadena de MontePython con Planck + BAO + `Lya_abgd`, con la cosmología y los parámetros astrofísicos libres. Lo que cambia es cómo se hace la pregunta. En lugar de comparar χ² entre modelos fijos, cada distribución de la ballena entra con una amplitud libre s, la fracción de la componente emitida respecto a la que predice el modelo:
- s = 0 es ΛCDM;
- s = 1 es la fila tal cual.

La cadena da la posterior de s, y la pregunta es si s = 1 está permitido.

**Aproximación (declarada).** La supresión del modelo se toma como una función de transferencia precalculada, T²(k; s) = P_modelo/P_ΛCDM, a la cosmología fiducial de U4. Se aplica a la P(k) ΛCDM de cada punto de la cadena. Es decir, se supone que T² no depende de los parámetros cosmológicos dentro de la región de Planck. Así cada paso solo ejecuta CLASS en ΛCDM; con la especie ncdm completa, cada paso tardaría decenas de minutos.

## Componentes (`lya/u15b/`)

- **`prepare_tables.py`:** calcula T²(k; s) con CLASS y el parche, a z = 0, con k ≤ 300 h/Mpc:
  - **filas:** `packet_params` de cada fila, con ω_ncdm = s·f_em·Ω_DM y ω_cdm = Ω_DM − ω_ncdm, para s ∈ {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3}, recortado a s·f_em ≤ 1;
  - **validación:** mezclas térmicas (`thermal_params`) de m = 1 y 2 keV con f ∈ {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6}.
- **`Lya_abgd_whale/`:** verosimilitud de MontePython que hereda de `Lya_abgd` sin tocar su código. Multiplica la P(k) lineal ΛCDM de cada punto por T²(k; s), interpolada linealmente en s. La pendiente n_eff se recalcula sobre esa P(k), y s es un parámetro nuisance de la cadena.
- **`params/*.param`:**
  - Planck 2018 high-ℓ TTTEEE *lite* + low-ℓ TT + low-ℓ EE;
  - BAO: `bao_boss_dr12` + `bao_smallz_2014`;
  - `Lya_abgd_whale`;
  - los 6 parámetros ΛCDM libres con los priors del `.param` de ejemplo de MontePython; los 9 astrofísicos con sus rangos y las comprobaciones de rejilla de la propia verosimilitud;
  - s con prior plano en [0, s_max].
- **`run_wsl.sh`:** instalación (MontePython en el commit `5b365b0`, el código y los datos de la verosimilitud de Planck 2018, la nueva verosimilitud) y una cadena por caso, con `mpirun`.
- **Convergencia:** R − 1 < 0.03 en todos los parámetros, con 4 cadenas por caso.
- **`analyze.py`:** posterior de s (o f) y su límite superior al 95 %, s_95.

## V · Validación (antes de leer ninguna fila)

Con las mezclas térmicas, el límite f_95 debe caer cerca del ajuste de Hooper et al. (ec. 4.1), f ≈ 0.1 + m/7.2 keV:
- m = 1 keV: f_95 ∈ [0.14, 0.38] (es decir, 0.24 multiplicado por [0.6, 1.6]);
- m = 2 keV: f_95 ∈ [0.23, 0.61] (0.38 multiplicado por [0.6, 1.6]);
- además, f_95(1 keV) < f_95(2 keV).

Si falla, U15b se informa como **INCONCLUSO (validación)** y no se leen las filas.

## Filas

Las 4 filas B2 de U12 con q = 1 y R_Γ = 1.56: x₀ ∈ {10, 30} y K_sw ∈ {0, 1}. Son las que llevan cola tibia; las de q < 1 tienen f_em < 1 % y no necesitan esta prueba.

## Veredicto

Por fila:
- **PERMITIDA:** s_95 ≥ 1;
- **EXCLUIDA:** s_95 < 1.

Global:
- **COLA_TIBIA_PERMITIDA_LYA_COMPLETO:** las 4 filas, PERMITIDAS;
- **COLA_TIBIA_EXCLUIDA_LYA_COMPLETO:** las 4, EXCLUIDAS;
- **COLA_TIBIA_PARCIAL_LYA_COMPLETO:** cualquier otra combinación.

Se informa también la posterior de s completa y si s = 0 queda fuera del 95 %, es decir, si los datos prefieren algo de cola.

## Lo que U15b no resuelve

- **T² fija:** el error de la aproximación no se mide. Se puede comprobar después repitiendo un punto de la cadena con la cosmología del mejor ajuste.
- **Planck lite** en lugar de la versión completa, para acelerar; es la compresión validada por la colaboración.
- **Solo R_Γ = 1.56 y q = 1.**

## Nota 1 · 12 de octubre de 2026, 11:20 (Madrid): preparación

- **Trabajo preparado en `lya/u15b/`** (README, `run_wsl.sh`), para ejecutarlo en WSL.
- **Prueba en la nube:** solo de la tubería, sin Planck (no hay clik en la nube): BAO + `Lya_abgd_whale` durante 6 pasos de MontePython 3.6.1 con CLASS v3.4.1, con una tabla de prueba de la mezcla de 1 keV (3 valores de f, k ≤ 30 h/Mpc). Funciona y da −lnL ≈ 30–32. No se ha leído ninguna posterior ni se ha evaluado ninguna fila.
- **Dos ajustes técnicos que no cambian el método:**
  - un parche de dos líneas a MontePython para Python 3.13 (`patch_montepython_py313.py`);
  - s se declara nuisance de `Lya_abgd_whale` en su fichero `.data`.
- **Las cadenas usan neutrinos sin masa (N_ur = 3.044)**, igual que las tablas de T². Se dice porque el análisis de referencia de Planck incluye un neutrino de 0.06 eV.

## Nota 2 · 12 de octubre de 2026, 18:25 (Madrid): pausas, reanudación y parte no markoviana

Las cadenas llevan corriendo unas horas en WSL. Hasta esta nota solo se ha mirado lo siguiente:
- el R − 1 máximo de cada caso, con `MontePython.py info`, a través de `status`;
- en los registros, el rango de s visitado, para comprobar que la cadena se mueve.

No se ha calculado ninguna posterior ni ningún s_95, ni se ha ejecutado `analyze.py` sobre cadenas reales.

**1. `analyze.py` quita la parte no markoviana.** Con `--update`, MontePython reajusta la propuesta mientras corre y marca cada reajuste en el fichero de la cadena. Los puntos anteriores al último reajuste no forman una cadena de Markov. `MontePython.py info` los quita por defecto, y por eso el R − 1 que muestra `status` ya está calculado sin ellos. `analyze.py` no lo hacía: solo descartaba el 30 % inicial. Es una corrección de código. Ahora, en cada fichero:
- se quitan los puntos anteriores a la última línea `update proposal`;
- de lo que queda, se descarta el 30 % inicial, como antes.

La regla de convergencia (R − 1 < 0.03 en todos los parámetros, calculado por `analyze.py`) y la del veredicto no cambian.

**2. Pausas y reanudación.** El equipo no puede tener la máquina ocupada varios días seguidos. Hay dos formas de liberarla:
- **`pause` y `resume`:** congelan los procesos y los reanudan tal cual. La cadena no cambia.
- **`stop` y `chains`:** paran las cadenas, y `chains` las continúa más tarde, cada una desde su último punto, con `-r` de MontePython. MontePython copia la cadena vieja al principio de la nueva. Una vez copiada, el fichero viejo se mueve a `restarted_from/` sin borrarlo, para que ningún punto cuente dos veces.

Ninguna de las dos formas cambia lo que se calcula.

**3. Parar un caso cuando converge.** Un caso se puede parar cuando `status` dé R − 1 < 0.03 en dos comprobaciones separadas por unas horas. Si después `analyze.py` da R − 1 ≥ 0.03 para ese caso, el resultado queda PENDIENTE y la cadena se continúa. Parar antes no cambia la regla del veredicto.

## Nota 3 · 13 de octubre de 2026, 09:50 (Madrid): fallo en el traslado de las cadenas reanudadas

La instancia de WSL se quedó colgada durante la noche y las cadenas se relanzaron con `chains`. No apareció la carpeta `restarted_from/`.

**El fallo.** `move_restarted.sh` (Nota 2) solo movía los ficheros viejos cuando cada cadena nueva superaba en puntos a la cadena vieja *más larga*, no a la suya. Con cadenas de longitudes distintas eso podía tardar más que su límite de 2 h, y entonces no movía nada. Mientras tanto, cada punto viejo estaba dos veces en la carpeta: en el fichero viejo y en su copia.

**La corrección.** `restarted.py` decide por contenido: un fichero sobra cuando todos sus puntos son exactamente los primeros de otro fichero más largo de la misma carpeta. Tolera una última línea cortada si la cadena murió a mitad de escritura.
- `analyze.py` usa esa misma comprobación para no contar nunca esos ficheros, aunque sigan en la carpeta.
- `status` avisa si hay copias sin mover, y `fixrestart` las mueve a `restarted_from/` sin borrarlas.

Esto no cambia el método ni el veredicto. Lo que pasó en concreto (si `chains` reanudó desde el último punto o empezó cadenas nuevas, según la versión del script que se usó) se anota cuando se revise la carpeta.

## Nota 4 · 12 de octubre de 2026, 23:20 (Madrid): corrección de `tail_preferred` (auditoría, antes de leer ninguna posterior)

*Escrita en la rama de auditoría como «Nota 3»; se numera 4 al fusionarla con la Nota 3 de `main` (13 de octubre, 09:50), que se escribió en paralelo. El contenido no cambia.*

La auditoría `audit/u12_u15b/` probó `analyze.py` con cadenas sintéticas de posterior conocida (`a1_analyze_synthetic.py`). s_95, el R − 1 y la eliminación de la parte no markoviana son correctos. Pero la etiqueta informativa `tail_preferred` («si s = 0 queda fuera del 95 %») estaba programada como s_05 > 0, que se cumple siempre que s ≥ 0: salía verdadera con una posterior plana y con una posterior con el máximo en s = 0.

Es un error de código y se corrige: ahora `tail_preferred` es verdadero si s = 0 queda fuera de la región de máxima densidad posterior que contiene el 95 % (histograma de 30 bins). También se guarda p(0)/p_max. Con una posterior casi plana la etiqueta es frágil y hay que leerla junto a p(0)/p_max. **No cambia s_95 ni el veredicto.** En el momento de esta nota no hay en el repositorio ninguna cadena ni ninguna posterior, y el auditor no ha visto ninguna.

## Nota 5 · 14 de octubre de 2026, 00:20 (Madrid): `status` solo leía omega_b; valores provisionales de s vistos antes de converger

**El fallo de `status`.** `MontePython.py info` escribe la etiqueta «-> R-1 is» solo en la línea del primer parámetro, omega_b. Los demás R − 1 salen en líneas sin etiqueta. `status` buscaba esa etiqueta, así que el «max R-1» que mostraba era en realidad el R − 1 de omega_b, nunca el de los otros 16 parámetros ni el de los derivados.

El 14 de octubre, a las 00:14, `status` daba entre 0.008 y 0.032 en los seis casos. En cambio, `analyze.py`, que calcula R − 1 en todos los parámetros con los cortes de la Nota 2, daba entre 0.19 y 5.1. Las cadenas **no han convergido**.

**Lo que se vio.** `analyze.py` imprimía s_95 y s_05 aunque las cadenas no hubieran convergido. Así, el equipo vio en esa ejecución los valores provisionales de s_95 de los seis casos, que se dejan anotados aquí:

| Caso | s_95 provisional |
|---|---|
| validación 1 keV | 0.26 |
| validación 2 keV | 0.34 |
| x10, K_sw = 0 | 1.41 |
| x10, K_sw = 1 | 1.92 |
| x30, K_sw = 0 | 2.69 |
| x30, K_sw = 1 | 2.80 |

Con R − 1 entre 0.19 y 5.1, estos valores no son un resultado.

**Qué cambia.**
1. **`status` ya no usa `info`.** Llama a `analyze.py --check`: el mismo R − 1, en todos los parámetros y con los mismos cortes, que decide el veredicto. Muestra los tres parámetros con el R − 1 más alto y no muestra nada sobre s.
2. **Nada sobre s antes de converger.** Mientras algún caso tenga R − 1 ≥ 0.03, `analyze.py` no imprime ni guarda ningún dato sobre s, y el resultado queda PENDIENTE.

La regla de convergencia (R − 1 < 0.03 en todos los parámetros) y la del veredicto no cambian.

**La regla de parada de la Nota 2 (punto 3) se precisa:** un caso se puede parar cuando `status` (es decir, `analyze.py --check`) dé R − 1 < 0.03 en dos comprobaciones separadas.

## Nota 6 · 14 de octubre de 2026, ≈21:10 (Madrid): la puerta de convergencia no puede aceptar NaN (auditoría de Codex, con las cadenas sin converger)

**El fallo,** señalado por Codex (`audit/codex_20261014`) y reproducido con cadenas sintéticas:
- `analyze.py` daba R − 1 = NaN a un parámetro sin varianza dentro de las cadenas, aunque difiriera entre ellas;
- `np.nanmax` lo ignoraba, así que dos cadenas congeladas en s = 0.2 y s = 0.8 salían con R − 1 máximo = 0 y convergidas.

El resultado local de U15b seguía PENDIENTE, con un máximo de 0.13–0.28: el fallo no aprobó nada. **Se corrige antes de que ninguna cadena converja.**

**La regla de convergencia pasa a exigir tres cosas** (`converged()`):
1. R − 1 < 0.03 y finito en todos los parámetros que varían. Si no hay varianza dentro de las cadenas pero sí entre ellas, R − 1 = ∞;
2. R − 1 < 0.03 también **por mitades**: cada cadena partida en dos (Gelman et al., BDA3). Así se detecta una deriva común a todas las cadenas;
3. al menos las **4 cadenas** por caso del preregistro.

**Se informa** el tamaño efectivo de muestra de s (`ESS_s`). El JSON se escribe en JSON estricto, con null en lugar de NaN o ∞.

El criterio del veredicto (s_95 ≥ 1) y el umbral de 0.03 no cambian. **La regla es más estricta,** no más laxa.

## Nota 7 · 16 de octubre de 2026 (Madrid): se detienen las cuatro filas físicas; sigue la validación

**Qué se detiene.** Las cuatro filas `U15b_B2_R1.56_q1_x{10,30}_K{0,1}` prueban la cola tibia del escenario anterior a U17: K_sw ∈ {0, 1}, q = 1, F = F₀. Ese ya no es el caso físico:
- **U17–U17d:** con el colapso cuerdas–paredes contado, K_eff ≫ 1;
- **U21–U20:** F₀ queda fuera de la banda de F.

**Por qué ahora.** Tras días de cadenas, sus R − 1 siguen entre 0.04 y 0.15, y por mitades entre 0.12 y 0.28. Se estima una semana más. **Se detienen sin haber leído nada sobre s** (la Nota 5 lo impide hasta la convergencia): la decisión no depende de los datos.

**Qué sigue.** Siguen las dos cadenas de **validación** (mezclas térmicas de 1 y 2 keV), que U15c necesita (enmienda 1 de ESPEC_U15c). `analyze.py` (`ROWS_STOPPED = True`) decide la convergencia solo con ellas y da:
- **DETENIDO (filas físicas superadas por U17–U21; validación PASA/FALLA),**
- con `validation_passed`, y sin ningún dato de s de las filas detenidas.

**El criterio de validación no cambia.** Los núcleos liberados se usan para U15c.

**Órdenes:**

```bash
bash lya/u15b/run_wsl.sh stop B2_R1.56           # para las cuatro filas físicas (las de validación siguen)
CASES="U15c_B2_F0.5_r1 U15c_B2_F0.35_r3" bash lya/u15b/run_wsl.sh chains   # tras generar sus tablas
```

## Nota 8 · 17 de octubre de 2026, ≈19:06 (Madrid): el muestreador, no la verosimilitud: covarianza desde las cadenas y sobremuestreo de Lyα

**Qué se ha visto.** Solo la convergencia: R − 1 por parámetro, como en la salida de `analyze` del 17 de octubre. La validación sigue lejos de la puerta (1 keV: 0.160, por mitades 0.254; 2 keV: 0.135 y 0.171), tras varios días. **No se ha leído ninguna posterior de s.**

**Causa.** Los parámetros que no convergen son los de molestia de Lyα (F_UV, gamma_a, T0a, Fz3).
- La plantilla dejó `data.over_sampling` en el valor por defecto de MontePython: un paso por bloque. Así, los parámetros de Lyα solo se mueven a la velocidad de CLASS, aunque moverlos no exige recalcular CLASS.
- La entrada de referencia de MontePython, cuyos priors y pasos dice seguir la plantilla (`input/base2018TTTEEE_lya_abgd.param`), usa `data.over_sampling=[1, 5, 4]`: sobremuestrea el bloque de Lyα ×4.
- Además, la máquina estaba sobresuscrita: 8 cadenas de U15b con 1 hilo, 8 de U15c con 5 y U19b con 34, unos 82 hilos en 44 núcleos.

**Cambio** (`bash lya/u15b/run_wsl.sh boost`, por caso):
1. Se detienen las cadenas.
2. Se comprueba la disposición de bloques con una ejecución de un paso en una carpeta aparte. Planck lite, low-ℓ TT y low-ℓ EE comparten A_planck, y MontePython forma 5 bloques: [cosmología, A_planck, —, —, Lyα].
3. Se rehace la matriz de covarianza de la propuesta con las propias cadenas (`info --want-covmat --keep-fraction 0.5`). La anterior se guarda al lado.
4. Se escribe `data.over_sampling=[1, 1, 1, 1, 10]` en `log.param`. La copia anterior se guarda al lado.
5. Cada cadena continúa desde su último punto (`-r`), con 2 hilos de CLASS por cadena.

**Lo que no cambia:**
- la verosimilitud, los priors, las tablas, los cortes del analizador, la puerta de convergencia y el veredicto;
- **la posterior objetivo es la misma.** Los puntos ya generados siguen contando: son muestras de la misma distribución, con una propuesta fija en cada tramo;
- la nueva covarianza la adopta cada cadena en su primer paso sin escribir «update proposal», así que la regla de corte del analizador no se dispara por este cambio.

Se aplica igual a U15c (enmienda 4 de su ESPEC).

## Nota 9 · 17 de octubre de 2026, ≈20:33 (Madrid): lo que pasó con las cadenas tras `boost`, y cómo se leen ahora las cadenas reanudadas

**Qué se ha visto.** Solo la salida de `analyze` y `status` que envió el equipo; nada de s.
- Tras `boost`, los dos casos de validación aparecen con 8 cadenas de 228 a 6165 puntos tras los cortes.
- El R − 1 subió (1 keV: 0.758, por mitades 0.898).
- `status` no avisaba de copias y no mostraba los casos de U15c.

**La causa más probable.** `-r` copia la cadena vieja, con sus comentarios, al principio de la nueva. Al seguir con `--update`, MontePython vuelve a calcular su R − 1 sobre toda la carpeta, y con cadenas recién empezadas puede reajustar la propuesta y escribir otra línea «update proposal».
- La regla de la Nota 2 quitaba **todo** lo anterior a la última de esas líneas, también los ~6000 puntos copiados, que la propia cadena vieja había generado con una propuesta fija.
- Las cadenas nuevas quedaban con ~250 puntos, y las viejas seguían contando aparte. Por qué no se reconocieron como copiadas lo dirá `run_wsl.sh diagnose` sobre las carpetas reales.

**Cambios** (código y lectura; ni la puerta de convergencia ni el veredicto cambian):
1. **El corte no markoviano se aplica por ejecución** (`lya/u15b/chainread.py`, que usan `analyze.py` y `analyze_u15c.py`).
   - Un fichero reanudado se parte en tramos: cada copia de una ejecución anterior y la ejecución nueva. Los límites salen de los ficheros viejos que contiene, también los de `restarted_from/`.
   - En cada tramo se quitan los puntos anteriores a su propia última línea «update proposal».
   - Motivo: un reajuste de la propuesta en un proceso posterior no hace no markovianos los puntos de un proceso anterior. Es lo que hace MontePython dentro de una sola ejecución, aplicado a cada una.
2. **Bifurcaciones.** Puede pasar que un fichero viejo siguiera creciendo después de que otro más nuevo copiara su principio, porque el proceso viejo no se paró.
   - Se avisa de ello.
   - El tramo común cuenta una sola vez, en el viejo; el nuevo cuenta desde el final del tramo común.
3. **`status` y `analyze --check`** muestran todas las carpetas de casos que encuentren, también U15c, con las copias saltadas y las bifurcaciones.
4. **`pause`, `resume` y `stop`** actúan sobre todas las cadenas que estén corriendo de verdad, no solo sobre la lista por defecto.
5. **`chains`** sin `CASES` lanza solo los casos activos: la validación de U15b y los dos casos de U15c. Las cuatro filas físicas (Nota 7) ya no se relanzan por accidente.
6. **Nueva orden `diagnose`:** por fichero, puntos, líneas «update proposal», tramos, copias y bifurcaciones.

Probado con carpetas sintéticas (tres generaciones de reanudación, una bifurcación) y con procesos simulados para `pause` y `stop`.

## Nota 10 · 17 de octubre de 2026, ≈20:39 (Madrid): `diagnose` encuentra la causa, una carrera de MontePython al reanudar

**Qué dice `diagnose` sobre las carpetas reales** (no se ha leído nada de s):
- En los dos casos de validación, las cuatro cadenas reanudadas por `boost` (`2026-10-17_300000__1–4`, unos 500 puntos cada una) **no llevan la copia** de las viejas (`2026-10-13_200000__1–4`, más de 10 000 puntos), que siguen arriba.
- En U15c con F = 0.35, dos de las cuatro reanudadas sí la llevan y dos no.

**La causa** está en `montepython/io_mp.py`. Al arrancar con MPI, cada proceso abre su fichero de cadena con el modo `'w'`, que lo vacía, y solo después intenta bloquearlo. Si dos procesos eligen el mismo nombre, el que pierde vacía el fichero del que ganó, que ya podía tener copiada la cadena vieja (`-r`). El ganador sigue escribiendo sus puntos nuevos a partir de un principio vacío. Las cadenas nuevas son, por tanto, continuaciones de las viejas: cada una arranca en el último punto de la vieja con el mismo número.

**Cambios** (lectura y lanzamiento; la puerta y el veredicto no cambian):
1. **`chainread.py` reconoce las continuaciones** y las lee junto con su cadena vieja, como una sola cadena.
   - Criterio principal: el primer punto del fichero nuevo es el último del viejo, con los mismos −lnL y parámetros; solo cambia la multiplicidad.
   - Si eso falla, por nombre: el mismo número de cadena y la base anterior más reciente. Se marca «(por nombre)».
   - Cada tramo conserva su propio corte no markoviano (Nota 9).
2. **`patch_montepython_restart.py`** corrige la carrera: abre con `'a'`, bloquea y solo el dueño vacía el fichero. `run_wsl.sh` lo aplica, de forma idempotente, antes de cada lanzamiento o reanudación, y también en `install`.

Las cadenas ya lanzadas no hay que tocarlas: la lectura nueva las une.

**Comprobación de la Nota 10 · 17 de octubre, ≈20:46 (Madrid).** Con las carpetas reales, `diagnose` une las 4 continuaciones de cada caso de validación, y 2 de U15c con F = 0.35, **por nombre**: el primer punto nuevo no coincide exactamente con el último viejo. El emparejamiento por nombre sigue la regla de MontePython (`-r` reanuda la cadena __i desde la vieja __i). `status` da ya 4 cadenas por caso. R − 1 (por mitades): 1 keV, 0.107 (0.224); 2 keV, 0.091 (0.151). Nada de s.

## Nota 11 · 18 de octubre de 2026, 21:55 (Madrid): la lectura de cadenas reanudadas, corregida tras una auditoría independiente; antes de cualquier veredicto con s

Una auditoría independiente construyó seis carpetas sintéticas adversas para `lya/u15b/chainread.py` y encontró tres errores de lectura. Ninguno afecta a lo ya informado: la validación (R − 1) no da veredicto y no se ha mirado nada de s. Pero podían contar puntos dos veces o unir cadenas ajenas.

1. **Una continuación que copia otra reanudación se contaba dos veces.** Es justo lo que pasa con la siguiente `chains` o `boost` tras la Nota 10. Ahora el linaje sigue a la copia: la vieja, su continuación y la reanudación que la copia son una sola cadena.
2. **La carrera de truncado puede dejar un hueco de bytes NUL seguido de la cola copiada.** Esa cola se contaba otra vez. Ahora se quitan los NUL iniciales y la línea parcial, y la continuación cuenta desde el final del tramo que repite el final de la vieja.
3. **El emparejamiento por nombre**, que usó la comprobación de la Nota 10, no es fiable. En la carrera, el sufijo __i de cada proceso depende de quién gana el cerrojo, no de qué cadena reanuda, así que puede unir dos cadenas independientes. Eso baja R − 1 y hace parecer convergido lo que no lo está.
   - Ahora una continuación se reconoce **solo por contenido**: sus primeros puntos, tras como mucho tres líneas dañadas, repiten el final de la vieja, sin mirar la multiplicidad.
   - Un fichero que solo comparte el sufijo cuenta como **cadena aparte**, con un aviso. Partir una cadena en dos es conservador para R − 1; unir dos cadenas no lo es.
4. **Además:**
   - las bifurcaciones se ordenan por (fecha, N, fecha de modificación), no solo por N;
   - `status` de `run_wsl.sh` ya no se aborta cuando pgrep no encuentra nada, y detecta «en pausa» con todos los procesos del caso;
   - `analyze_u15c.py` ya no guarda `n_points_near_1` (un recuento sobre s) antes de converger, y distingue la validación no disponible (PENDIENTE) de la fallida (INCONCLUSO).

**Qué hacer en WSL antes de leer nada:**
- Ejecutar `bash lya/u15b/run_wsl.sh diagnose` con este código.
- Si las continuaciones de la Nota 10 no repiten el final de la vieja, cada caso puede pasar de 4 a 8 cadenas. Eso es más conservador; las R − 1 de la comprobación de la Nota 10 (0.107 y 0.091) pueden cambiar.
- Mirar el principio de esos ficheros (`head -c 64 <fichero> | od -c`) para saber si llevan un hueco de NUL.

La prueba `test_chainread_restart_segments` cubre ahora el hueco de NUL con la cola copiada y el fichero que solo comparte el sufijo.


## Nota 12 · 19 de octubre de 2026, ≈11:40 (Madrid): una salida partida en dos por un movimiento prematuro, en U15c_B2_F0.5_r1; antes de leer nada de s

**Qué se vio.**
- Al continuar las cadenas el 19 de octubre (`run_wsl.sh chains`), tres y luego cuatro casos quedaron bien. En `U15c_B2_F0.5_r1`, en cambio, las cadenas `__1`, `__2` y `__4` nuevas (`2026-10-19_300000__k`) arrancaron sin la copia de sus cadenas antiguas: unos 1000 puntos, sin hueco de NUL.
- `chainread` las contó, con razón, como cadenas aparte (Nota 11): 7 cadenas y R − 1 = 0.27–0.42.
- En `restarted_from/` aparecían ficheros con **el mismo nombre** que los nuevos.

**Causa.**
- `restarted.py` decidía que un fichero estaba «contenido» en otro sin mirar cuál era el más reciente. Justo después de reanudar, un fichero nuevo que todavía está copiando su cadena antigua es un **prefijo** del fichero viejo, así que se tomaba por contenido en él.
- `move_restarted.sh` movió así los tres ficheros nuevos a medio copiar, junto con la `__3` vieja, que sí estaba copiada. MontePython siguió escribiendo en los ficheros movidos hasta que reabrió su salida **por nombre** (`io_mp.refresh_file`). Desde ese momento escribió en un fichero nuevo con el mismo nombre arriba.
- Resultado: la salida de un solo proceso quedó en dos ficheros. `restarted_from/<nombre>` tiene la copia de la cadena vieja y los primeros puntos nuevos; `<nombre>`, arriba, tiene los puntos siguientes. No se perdió ni se duplicó ningún punto.
- No se ha visto nada de s.

**Corrección** (código, sin tocar la regla de convergencia ni la del veredicto):
1. `restarted.py`: un fichero solo puede contener a otro si su N en el nombre es mayor. Una reanudación siempre escribe un N mayor que el del fichero que copia.
2. `rejoin_moved.py` (nuevo), solo con las cadenas del caso paradas:
   - escribe `<nombre> = restarted_from/<nombre> + <nombre>`;
   - guarda los dos originales y un registro con sus sha256 en `restarted_from/_rejoined_<fecha>/`, sin borrar nada;
   - después mueve los ficheros viejos que ahora quedan contenidos.
3. Prueba nueva: `test_restart_move_and_rejoin` reproduce la secuencia y comprueba que la cadena reunida es una sola, con todos sus puntos.

**Qué se hace en WSL:** parar `U15c_B2_F0.5_r1`, reunir sus ficheros, comprobar con `diagnose` que vuelven a ser 4 cadenas, y continuar. Cada cadena reunida es la salida de un único proceso de MontePython, que continúa su cadena antigua con el mismo sufijo (`io_mp.create_output_files` toma como fichero de reanudación el del sufijo que obtiene). Por eso reunirlas no es elegir cómo emparejar cadenas.


## Nota 13 · 20 de octubre de 2026, 09:57 (Madrid): resultado

- **El análisis:** el autor ejecutó `analyze.py` en WSL el 20 de octubre, a las 06:37, sobre las cadenas finales. Los JSON y las tablas se suben al repositorio; las cadenas quedan en WSL.
- **La validación PASA:**
  - 1 keV: f₉₅ = 0.271, en [0.14, 0.38];
  - 2 keV: f₉₅ = 0.422, en [0.23, 0.61];
  - y f₉₅(1 keV) < f₉₅(2 keV);
  - R − 1 ≤ 0.021.
- **Veredicto: DETENIDO** (filas físicas superadas por U17–U21; validación PASA).
- **De las filas físicas detenidas (Nota 7)** no se calcula nada de s.
- **Informe:** `docs/informes/INFORME_U15b_U15c.md`.
