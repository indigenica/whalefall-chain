# U15b · Planck + BAO + Lyman-ABGD para la cola tibia (WSL)

Preregistro: `docs/preregistros/ESPEC_U15b.md`. Todo se lanza con `lya/u15b/run_wsl.sh`, desde la carpeta del repositorio y en el entorno `lighthouse-sage`, que ya tiene el CLASS con el parche.

```bash
conda activate lighthouse-sage
cd "/mnt/f/Proyectos/Cencia/Akashic alpha engine/whalefall/WHALEFALL/claude_dual/unification"
git pull
bash lya/u15b/run_wsl.sh install    # MontePython, clik 16 (código de Planck, de GitHub), los datos de Planck 2018 (PLA) y Lya_abgd_whale
bash lya/u15b/run_wsl.sh tables     # 51 ejecuciones de CLASS → lya/u15b/tables/*.npz   (~1 h con 40 procesos)
bash lya/u15b/run_wsl.sh test       # 30 pasos del caso de validación de 1 keV
bash lya/u15b/run_wsl.sh chains     # 6 casos × 4 cadenas MPI, en segundo plano (1–3 días)
bash lya/u15b/run_wsl.sh status     # R − 1 de cada caso; se para cuando todos estén por debajo de 0.03
bash lya/u15b/run_wsl.sh pause [x]  # congela las cadenas (no usan CPU; la memoria sigue ocupada); `resume [x]` las sigue tal cual
bash lya/u15b/run_wsl.sh stop [x]   # las termina; `chains [x]` las continúa después desde su último punto (sirve tras reiniciar)
                                    # [x], opcional: parte del nombre del caso (x30, x10_K0, val_thermal…); sin él, todos
bash lya/u15b/run_wsl.sh analyze    # → lya/results/u15b/U15B_RESULT.json
```

## Variables de entorno

| Variable | Por defecto | Qué controla |
|---|---|---|
| `U15B_WORK` | `~/u15b` | MontePython, Planck y las cadenas (fuera del repositorio) |
| `CLASS_DIR` | `~/src/class_public` | El CLASS con el parche (`env/INSTALAR_WSL.md`) |
| `NPROC` | 40 | Procesos para las tablas |
| `CHAINS_PER_CASE` | 4 | Cadenas MPI por caso; 6 casos × 4 = 24 núcleos |
| `STEPS` | 100000 | Pasos máximos por cadena; se puede parar antes, cuando R − 1 < 0.03 |

## Qué hay aquí

- **`prepare_tables.py`:** T²(k; s) = P_modelo/P_ΛCDM a z = 0, para las 4 filas B2 con q = 1 y R_Γ = 1.56 (las de U12) y para las dos mezclas térmicas de validación.
- **`Lya_abgd_whale/`:** verosimilitud de MontePython. Hereda de `Lya_abgd` (Hooper et al. 2022) sin modificarla y aplica T²(k; s) a la P(k) ΛCDM de cada paso. La amplitud s es un parámetro nuisance.
- **`params/template.param`:**
  - Planck 2018 high-ℓ TTTEEE lite + low-ℓ TT + low-ℓ EE;
  - BAO: BOSS DR12 + 6dF/MGS;
  - `Lya_abgd_whale`.
- **`patch_montepython_py313.py`:** MontePython importa las verosimilitudes con dos `exec` encadenados, y eso falla en Python 3.13 (PEP 667). El parche cambia solo esas dos líneas.
- **`analyze.py`:** posterior de s, límite al 95 %, R − 1, validación y veredicto.

## Si algo falla

- **La verosimilitud de Planck no compila** (`./waf configure --install_all_deps`): instaladla con Cobaya, que la parchea y compila:
  ```bash
  pip install cobaya
  cobaya-install planck_2018_highl_plik.TTTEEE_lite planck_2018_lowl.TT planck_2018_lowl.EE --packages-path ~/u15b/cobaya
  ```
  Después, en `~/u15b/montepython_public/default.conf`, apuntad `path['clik']` a la carpeta `plc-3.1` que crea Cobaya, y exportad `CLIK_DIR=…/plc-3.1` y `PLANCK=…/baseline/plc_3.0` antes de `test` y `chains`.
- **La URL de datos del Planck Legacy Archive no responde:** se puede descargar a mano desde pla.esac.esa.int → Cosmology → Likelihood (`COM_Likelihood_Data-baseline_R3.00.tar.gz`) y dejarla en `~/u15b/planck/data.tar.gz`.
- **Cadenas en marcha:**
  - `chains` fija un hilo de CLASS por cadena (`OMP_NUM_THREADS`) y comprueba antes que `mpirun` + `mpi4py` funcionan;
  - los registros quedan en `~/u15b/logs/<caso>.log` y las cadenas en `~/u15b/chains/<caso>/`;
  - hay que dejar abierta una ventana de WSL, porque al cerrarlas todas Windows puede apagar la distribución y, con ella, las cadenas. Si se paran, `chains` las continúa desde su último punto (con `-r` de MontePython). Los ficheros viejos se mueven a `restarted_from/` una vez copiados;
  - **para liberar la máquina:**
    - `pause` es lo más limpio si WSL no se reinicia;
    - si se va a reiniciar o apagar, `stop` y, más tarde, `chains`;
    - los casos que ya han convergido (R − 1 < 0.03 en dos `status` seguidos) se pueden parar con `stop <caso>`;
    - ESPEC_U15b, Nota 2;
  - `status` muestra el R − 1 de todos los parámetros (los tres peores), calculado igual que en el veredicto (`analyze.py --check`). No muestra nada sobre s: hasta que todo converja, `analyze` tampoco lo hace (ESPEC_U15b, Nota 5).
- **`mpirun` falla con Intel MPI** (el de conda): ya se exporta `I_MPI_FABRICS=shm`. Si aun así falla, `conda install -c conda-forge "mpi4py=*=*openmpi*"` cambia a Open MPI.
- **Requisitos del sistema para compilar clik:** `sudo apt install -y build-essential gfortran liblapack-dev libblas-dev libcfitsio-dev`. Se compila con el gcc y el gfortran del sistema, no con los compiladores de conda.
- **Código de Planck:** se usa clik 16 (github.com/benabed/clik, commit `9c6a194`), sucesor mantenido de PLC 3.x que lee los mismos ficheros `.clik`. El código PLC 3.01 del archivo no construía su envoltorio de Python con Python 3.13. Si el envoltorio no se construye con `waf`, `build_clik_wrapper.py` lo compila con Cython.
- **Comprobado en la nube:** el arnés funciona sin Planck (BAO + `Lya_abgd_whale`, 6 pasos con MontePython 3.6.1 y CLASS v3.4.1). La parte de Planck no se ha podido probar allí.

## U15c (ESPEC_U15c): la rama viva con la misma tubería

```bash
python -B lya/u15b/prepare_tables.py --u15c --workers 8                 # 14 CLASS runs -> tables/U15c_*.npz, manifest_u15c.json
CASES="U15c_B2_F0.5_r1 U15c_B2_F0.35_r3" bash lya/u15b/run_wsl.sh chains   # 2 cases x 4 chains = 8 cores
CASES="U15c_B2_F0.5_r1 U15c_B2_F0.35_r3" bash lya/u15b/run_wsl.sh status
python -B lya/u15c/analyze_u15c.py ~/u15b/chains                           # -> lya/results/u15c/U15C_RESULT.json
```
Sin `CASES`, `run_wsl.sh` sigue con los seis casos de U15b, como siempre.

## Acelerar la convergencia sin cambiar nada de la física (ESPEC_U15b, Nota 8)

```bash
export OMP_NUM_THREADS=2                                                     # hilos de CLASS por cadena
bash lya/u15b/run_wsl.sh boost val                                           # validación de U15b
CASES="U15c_B2_F0.5_r1 U15c_B2_F0.35_r3" bash lya/u15b/run_wsl.sh boost      # U15c
```

`boost` hace cinco cosas por caso:
1. para las cadenas;
2. prueba en una carpeta aparte (`~/u15b/tests/`) que MontePython acepta `data.over_sampling=[1, 1, 1, 1, 10]`;
3. rehace la covarianza de la propuesta con las propias cadenas (la anterior se guarda al lado);
4. escribe el sobremuestreo en `log.param` (la copia anterior se guarda al lado);
5. continúa cada cadena desde su último punto.

Los puntos ya generados siguen contando. Si la prueba del paso 2 falla, las cadenas quedan paradas sin tocar y `chains` las reanuda como estaban.
