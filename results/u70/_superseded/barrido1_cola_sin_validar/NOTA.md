# Primer barrido de U70 (25 de noviembre de 2026), superado

La revisión de Codex de la PR #106 encontró un error. En las filas de relatividad general con el modo radial dinámico y carga grande, el tiro terminaba antes de llegar al vacío asintótico, y la cola de la acción, que supone r = 1 y un exterior casi plano, se sumaba igualmente.

Afectaba a 41 filas, todas de relatividad general en (D) con n ≥ 1701. En algunas, la cola era casi toda la acción.

No afectaba ni al veredicto (PROTEGE, decidido por las filas de NPQT en (C), que llegan a 10⁴ a₀) ni a las puertas (G2 usa n = 1).

El código corregido marca esas filas como no resueltas y las deja fuera de S_min. El barrido se repitió con él.
