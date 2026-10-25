# chain/ · la cadena unificada

| Etapa | Qué hace | Salida |
|---|---|---|
| `s1_birth.py` | Nacimiento: rehace los hijos viables de C45j con `components/c45b_hom/frw_core.c` (necesita gcc) | `results/s1_birth.json` |
| `s2_pocket.py` | Bolsillo: ajuste del pivote, umbral eterno y planitud D3 | `results/s2_pocket.json` |
| `s3_sectors.py` | Especies, desintegraciones, fluido tardío y estado de la red | `results/s3_sectors.json` |
| `s4_dark.py` | Libro oscuro vectorizado, flujo libre y referencia WDM | `results/s4_dark_checks.json` |
| `s5_closure.py` | Cierre U1: barrido preregistrado de ε | `results/s5_closure.json` |
| `s6_measured.py` | Cierre con la ε medida en `network/` (ramas B1/B2/B3) | `results/s6_measured_closure.json` |
| `u_contract.py` | Contrato de cierre | `results/CONTRATO_CIERRE.json` |

`run_all.py` lo ejecuta todo (`make chain`). `ucommon.py` fija las rutas y el único conjunto de parámetros (`benchmark()`).
