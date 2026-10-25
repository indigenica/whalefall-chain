# Ejecución parcial de U12 en la nube

Logs de la ejecución en la sesión de la nube (2 núcleos):
- **run1:** rejilla térmica completa y 2 filas; se cortó por un reinicio de la máquina a las 06:00.
- **run2:** reanudación desde la caché de CLASS; se detuvo a mano el 12 de octubre a las 09:10, con 16 de 31 filas hechas, porque U12 ya había terminado en WSL.

No se calculó ningún δA ni veredicto en la nube: los logs solo registran que cada ejecución de CLASS terminó. El resultado válido de U12 es el de WSL (`lya/results/u12/U12_RESULT.json`). Las cinco referencias térmicas de la nube coinciden bit a bit con las de WSL.
