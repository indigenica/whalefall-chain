# c45b — colapso con horizonte, rebote QT y lanzamiento (ballena 4)

- `INFORME_C45B.md`: informe (español). Auditoría independiente: `audit/AUDIT_C45B.md`.
- `dn/`: código de doble nulo NPQT (ecuaciones, solver con centro regular y enmascarado causal, análisis, convergencia, perfiles).
- `dn_bank/run_c45b_dn.py`: banco preregistrado (PROTOCOL en el propio fichero y en results/*/PROTOCOL.json).
- `hom/`: parche homogéneo cerrado solo de χ (port de C45), banco C45b-0.
- `dn_runs/`: sondas exploratorias declaradas.
- `results/`, `figures/`, `_superseded/`.
- `INFORME_C45B_LTB.md`: prueba LTB-QT de coherencia del rebote (FAIL para la estrella de χ) y proyección a la cuenca de C45. Código en `ltb/`, `hom/projection*.py`; auditoría en `audit_ltb/`.
- `INFORME_C45B_PBH.md`: test decisivo con perfiles tipo PBH (modo creciente en padre de polvo). Rama local cerrada. Código `ltb/growing.py`, `ltb/bank_pbh.py`, `ltb/coherence.py`; auditoría `audit_pbh/`.
- `INFORME_C45B_RESCATE.md`: búsqueda exploratoria de rescates (QT blanda, borde nítido). Código `ltb/growing_qt.py`, `ltb/soft_qt.py`, `ltb/joint_soft.py`, `ltb/edge_scan.py`, `hom/c45_soft.py`, `hom/qt_inflation.py`.
