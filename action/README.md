# action/ · la acción común (Q + χ + radiación + NPQT)

Copia, sin cambios, de la campaña `cierre_accion` (4 de octubre); el modelo está en `docs/modelo/MODELO_ACCION_COMUN.md`.
- `symbolic.py`: residuos simbólicos exactos de la conservación y la restricción (necesita sympy).
- `portal.py`: compleción acotada del portal visible.
- `experiments.py`: integración homogénea con NPQT.
- `decays.py`: desintegraciones a nivel árbol (sympy).
- `late_fluid.py`: fluido tardío y ΔN_eff homogéneo.

**Desde U26 (18 oct)** la meseta es la inflación de polo: `make action` usa `vendor/thermal_pole.json` y escribe en `results_pole/` (congelado en `data/cierre_accion_pole/`), y `portal.py` hace la cota de no negatividad con el polo. Con `WHALEFALL_PLATEAU=sum` la carpeta es la copia sin cambios de `cierre_accion` y escribe en `results/`.

Con `WHALEFALL_PLATEAU=sum`, `make action` lo ejecuta todo y escribe en `results/`; `results_reference/` guarda los resultados originales. Sin sympy, `portal`, `experiments` y `late_fluid` (este último con `results_reference/decays.json` copiado a `results/`) reproducen los originales con diferencias ≤ 10⁻⁵, salvo diagnósticos de error de redondeo.
