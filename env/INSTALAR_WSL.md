# U4 en `lighthouse-sage` (WSL Ubuntu-24.04 de indi-xeon)

Esta guía usa **vuestro entorno de ciencia `lighthouse-sage`**: Python 3.13, numpy 2.5, scipy 1.16 y JAX con GPU. Es el mismo con el que la cola (`WHALEFALL/claude/queue`) ejecuta los jobs. Solo hay que **añadirle CLASS** y, si faltan, sympy y cython; no se cambia nada de lo que ya tiene.

**Máquina** (según `claude/logs/J000_envcheck.log`): 44 núcleos, 31 GB de RAM visibles en WSL y una RTX 3060 Ti de 8 GB.
- Para U4, **la CPU es lo que manda**: 44 núcleos para CLASS en paralelo y para el kernel OpenMP de las redes.
- Los 8 GB de la GPU no alcanzan para las redes que hacen falta (una RP³ de 384³ ya necesita ~7 GB solo para los campos), así que U4 no usa JAX.

## 0. Dar más memoria a WSL (muy recomendable para U4b)

WSL ve 31 GB, probablemente el 50 % de los 64 GB de Windows. El tamaño de red que cabe depende de eso:

| RAM de WSL | U(1), fase 1 | RP³, fase 2 |
|---|---|---|
| 31 GB (ahora) | N = 576, ℓ_fin ≈ 5.66 | N = 448, ℓ_fin ≈ 5.41 |
| ~52 GB | N = 768, ℓ_fin ≈ 5.95 | N = 512, ℓ_fin ≈ 5.55 |

En Windows, cread o editad `C:\Users\<usuario>\.wslconfig`:

```
[wsl2]
memory=52GB
processors=44
```

Después, en PowerShell: `wsl --shutdown`. Al volver a WSL, `free -g` debe mostrar ~52 GB. Ajustad el valor si Windows tiene otra RAM; dejad unos 10 GB para Windows.

## 1. Comprobar qué le falta a `lighthouse-sage`

```bash
conda activate lighthouse-sage
REPO="/mnt/f/Proyectos/Cencia/Akashic alpha engine/whalefall/WHALEFALL/claude_dual/unification"
cd "$REPO"
python env/check_env.py
```

Lo normal es que salgan en FALLA solo `classy`, y quizá `sympy`. El resto (numpy, scipy, matplotlib, gcc + OpenMP) debería estar en OK.

## 2. Instalar lo que falte, sin tocar lo demás

```bash
conda activate lighthouse-sage
conda install -c conda-forge cython sympy        # solo si check_env los marca; conda no cambia numpy/scipy si ya cumplen
which gcc || sudo apt install -y build-essential  # gcc del sistema para CLASS y para el kernel de red
```

## 3. Compilar e instalar CLASS (classy), con el parche de cuadratura, dentro de `lighthouse-sage`

U4a necesita un parche de 18 líneas en CLASS (`env/class_qm_trapz_log.patch`). Añade una regla de cuadratura en momentos uniforme en ln q (`ncdm_quadrature_strategy = 4`) y no toca nada más; el motivo está en la enmienda 1 de `docs/preregistros/ESPEC_U4.md`. Sin él, las filas de U4a tardan días.

```bash
conda activate lighthouse-sage
mkdir -p ~/src && cd ~/src
git clone https://github.com/lesgourg/class_public.git       # (si ya lo teníais clonado, saltad esta línea)
cd ~/src/class_public
git log -1 --format="%H %cd" > "$REPO/env/CLASS_VERSION.txt"   # versión exacta usada (queda en el repositorio)
git apply "$REPO/env/class_qm_trapz_log.patch"              # el parche (se comprueba con: git diff --stat → 2 ficheros)
make clean
make -j 40                                        # class + libclass.a, con el gcc del sistema
CC=gcc pip install --no-build-isolation --no-deps .   # instala classy usando el numpy/cython del entorno, sin tocar dependencias
python -c "import classy; print('classy OK', classy.__file__)"
```

**Si ya instalasteis CLASS sin el parche** (6–7 de octubre): basta con las líneas desde `cd ~/src/class_public` hasta `pip install`, sin el `git clone`. Si `git apply` dice que ya está aplicado, no hace falta repetirlo.

- `--no-deps` y `--no-build-isolation` garantizan que pip **no actualiza ni cambia** numpy, scipy ni nada más del entorno.
- Si algún día queréis quitarlo: `pip uninstall classy`.

## 4. Comprobar todo

```bash
cd "$REPO"
python env/check_env.py     # tiene que terminar en "ENTORNO LISTO"
make test                   # 9 pruebas, ~2 min
```

`check_env.py` comprueba en concreto dos cosas:
- que **CLASS lee bien una distribución de momentos desde fichero** (WDM de 3.5 keV leído de fichero frente a la versión analítica, con < 0.5 % de diferencia);
- que **el parche de cuadratura está instalado** (línea «CLASS con el parche de cuadratura logarítmica»). Un CLASS sin parche no da error por sí solo, pero deja los pesos sin calcular; por eso se comprueba aquí y otra vez al principio de U4a (C0b).

Si alguna falla, no lancéis U4a y pasadme la salida.

## 5. Lanzar U4

**A mano (recomendado ahora).** Mirad el bloque «A mano» más abajo.

**Con la cola.** El runner (`WHALEFALL\claude\runner\START_RUNNER.bat`, doble clic en Windows) abre WSL, activa `lighthouse-sage`, se sitúa en `WHALEFALL` (por eso los jobs hacen `cd claude_dual/unification`) y ejecuta por orden de nombre los `.sh` de `WHALEFALL\claude\queue`. Cada log queda en `claude\logs\<job>.log` y el job pasa a `claude\done\`. Para pararlo se crea un fichero vacío `claude\runner\STOP`.
- **Ojo:** ahora mismo la cola tiene también 13 jobs antiguos (J078–J090) que se ejecutaron en otro sitio. Si arrancáis el runner tal cual, los ejecutaría primero. Antes de usarlo, movedlos a otra carpeta (por ejemplo `claude\queue_aparcados\`).

Los jobs de U4 son estos:

| Job | Qué hace | Tiempo aprox. (44 núcleos) |
|---|---|---|
| `J094_u4_envcheck.sh` | `check_env` + `make test` + plan de redes (lanzadlo primero) | 5 min |
| `J091_u4a_class.sh` | U4a: Lyα con CLASS (comprobación rápida + prueba completa; se puede reanudar) | 3 h + C1 (6–12 h) |
| `J092_u4b_phase1.sh` | U4b fase 1: 4 redes U(1) al N que quepa + puertas | 3–8 h (N = 576–768; las medidas son en parte en serie) |
| `J093_u4b_phase2.sh` | U4b fase 2: RP³ + U(1), solo si la fase 1 pasa; luego `make chain` | 2–4 h |

**A mano**, en `tmux` (es lo mismo que hacen los jobs):

```bash
conda activate lighthouse-sage && cd "$REPO"
make u4a-quick      # ~1 h: C0, C0b, referencias, familia pequeña y 10 filas (sin C1); va imprimiendo [n/N] por ejecución
make u4a            # ~2 h de filas (62, de 20–40 min cada una, una por núcleo) + C1 (dos ejecuciones de 6–12 h en paralelo)
make u4b-plan       # revisad el N y las horas
make u4b-phase1 && make u4b-analyze
make u4b-phase2 && make u4b-analyze    # solo si u4b-analyze dice GATES_PASSED
make chain
```

**Hilos.**
- U4a reparte las ejecuciones de CLASS en procesos de un hilo cada uno. Por defecto usa min(núcleos, RAM libre/1.2 GB, 32); con `WORKERS=` se cambia.
- U4b usa todos los núcleos en el kernel OpenMP; con `U4_THREADS=` se cambia.
- No lancéis U4a y U4b a la vez si la RAM es de 31 GB.

**Cómo saber que va bien.** Cada ejecución de CLASS que termina escribe una línea `[n/N] … ok (x min)`. Entre línea y línea puede pasar media hora: cada fila es una ejecución de CLASS de un solo núcleo. En otra terminal:
- `top`: tiene que haber tantos procesos `python` al ~100 % como ejecuciones en marcha;
- `ls lya/results/u4a/cache | wc -l`: crece con cada ejecución terminada;
- `cat lya/results/u4a/class_errors.log`: si existe, CLASS ha fallado en alguna ejecución.

Se puede parar con Ctrl+C y relanzar: lo ya calculado está en la caché y no se repite.

## 6. Resultados

Se escriben directamente en el repositorio:
- `lya/results/u4a/U4A_RESULT.json`
- `network/results/u4/U4B_RESULT.json` y las redes en `network/results/u4/runs/`
- los logs en `logs/`

La caché de CLASS (`lya/results/u4a/cache`) se puede borrar al final. Cuando terminen, decídmelo y hago el análisis, la auditoría y la actualización de los documentos.

## Problemas frecuentes

- **`make` de CLASS falla con un error de OpenMP:** comprobad que `which gcc` da `/usr/bin/gcc` y que está `build-essential`.
- **`pip install .` falla con `-lmvec` o con el compilador de conda:** `CC=gcc pip install --no-build-isolation --no-deps .`, tal como está en el paso 3.
- **Python 3.13 y Cython:** CLASS necesita Cython ≥ 3.0 con Python 3.13. Si `conda list cython` da una versión 0.29, actualizadlo: `conda install -c conda-forge "cython>=3"`.
- **Una red se queda sin memoria:** relanzad con un N menor, por ejemplo `make u4b-phase1 N=512`. Las redes ya terminadas se saltan.

---

### Alternativa: un entorno aparte (`whalefall`)

Si preferís no tocar `lighthouse-sage`:

```bash
conda env create -f env/environment.yml && conda activate whalefall
```

Después, el paso 3 igual. En los jobs, exportad `WF_PYTHON=$HOME/miniforge3/envs/whalefall/bin/python`.
