# lya/ · Lyman-α con una transferencia de Boltzmann lineal (U3a)

> **Meseta (U26, 18 de octubre de 2026).** Los resultados guardados en esta carpeta se calcularon con la suma de 10 niveles, antes de que la meseta por defecto pasara a ser la inflación de polo.
> - **Para reproducirlos**, ejecuta con `WHALEFALL_PLATEAU=sum`. Así leen `benchmark()` y `results/plateau_sum/`. Los objetivos `make lya`, `make u4a` y `make u4a-quick` ya lo fijan.
> - **Con el polo** (el valor por defecto), m_G cambia un +0.4 % y T_Γ un +15 %, lo que no cambia ningún veredicto (`docs/informes/INFORME_U26.md`).


- `boltz.py`: resolvedor lineal de Einstein–Boltzmann en gauge síncrono (Ma & Bertschinger 1995).
  - CDM; bariones que no se agrupan antes de a = 1/1090; radiación sin masa con RSA.
  - Especie no fría con distribución de momentos arbitraria (jerarquía Ψ_l(q)).
  - δA de Murgia et al. 2017.
- `run_u3a.py`: validaciones W1/W2/W3, referencias (WDM de 5.3 keV y mezclas en el borde de Garcia-Gallego et al. 2025) y puerta sobre las filas de S6.
- `w3_decisive.py`, `figure_u3a.py`.
- `results/U3A_RESULT.json` (k_max = 1000 h/Mpc, final); `U3A_RESULT_kmax120.json` (primera pasada, no convergida en k_max).

**Estado:**
- No validado según el preregistro: W1 frente a Viel et al. falla en un 11–18 % en k_{1/2}.
- Auditado y reproducido de forma independiente (`audit/u3/`).
- Pendiente: confirmación con `classy`.
