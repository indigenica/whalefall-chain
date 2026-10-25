# U22 · La meseta como inflación de polo (orden 4): ¿se puede quitar el parche D2 sin cambiar las predicciones? · Preregistro

*15 de octubre de 2026, ≈20:00 (Madrid). Escrito antes de ejecutar U22 (la hora que cuenta es la de su commit). Lo que ya se calculó antes, y se declara: el diagnóstico en slow-roll de `audit/referee_20261015/RESPUESTA.md` §1, con N fijo = 51.55: n_s 0.97357 → 0.97343 y r 0.004858 → 0.004856. Se ejecuta en la nube.*

## Motivo

El potencial del inflatón es una suma postulada de 10 niveles, Σ 2^−j tanh²(χ/(0.5·2^j)) (hipótesis D2, «q·s = 1 postulado»). Es una versión discreta y truncada de la cola V₀(1 − C/χ). Esa cola sale exactamente de un lagrangiano estándar, la **inflación de polo de orden 4**:
- −C²/(2φ⁴)(∂φ)² − V₀(1 − φ), con χ = C/φ;
- Galante, Kallosh, Linde y Roest 2014; Broy et al. 2015;
- protege la meseta de las correcciones radiativas.

Si las predicciones no cambian, **D2 deja de ser un parche.**

## Modelo

**Potencial:** V = V₀·[1 − C/√(χ² + C²)] (polo regularizado), con V(0) = 0, m_χ² = V₀/C² y cola exacta V₀(1 − C/χ).
- **C = 0.3516** es el ajuste geométrico de la cola de la suma en χ ∈ [5, 200]. Es un ajuste de forma, no a datos.
- **Sensibilidad,** sin veredicto: C × {0.8, 1.25}.

## Qué se ejecuta (`chain/u22_pole.py`)

1. **S2 (CMB):** el mismo flujo que `cadena_cierre/thermal.normalize_continuous`, sin tocar su código:
   - final exacto en ε_H = 1;
   - recalentamiento continuo con σ/m = 1e-5;
   - empalme del pivote.

   Se sustituyen el potencial, la masa m_χ y `predict` por los del polo. Salen N_*, n_s, r, H_*, T_rh y F (con μ/H_* = 0.30, como el benchmark).
2. **Umbral eterno:** χ_et con P_ζ = 1, con el V₀ común de S2.
3. **S1 (nacimiento):** las 90 ejecuciones de S1 (los 20 hijos viables de C45j y 70 de control) con el puerto C de C45 en el modo nuevo 4, el polo de `frw_core.c` (los modos 0–3 no cambian; S1 sigue reproduciendo bit a bit, comprobado). Se cuentan:
   - los hijos viables que siguen siendo eternos (χ_land > χ_et y E_phys,max < 1);
   - f_ES en la ventana fértil.

## Veredicto

- **POLO_EQUIVALENTE:** se cumplen las tres condiciones:
  - |Δn_s| ≤ 0.002 (la mitad de σ de Planck);
  - |Δr|/r ≤ 10 %;
  - ≥ 15 de los 20 hijos viables siguen siendo eternos.

  D2 se puede sustituir por un lagrangiano estándar.
- **POLO_CAMBIA_CMB:** falla la condición de n_s o la de r.
- **POLO_CAMBIA_NACIMIENTO:** el CMB es equivalente, pero quedan < 15 hijos eternos.

## Lo que se espera, escrito antes de calcular

- Δn_s ~ 10⁻⁴ y r igual.
- χ_et sube de ~416 a ~550 M_Pl (HIPOTESIS_BALLENA_5: 551 sin truncar), así que ~18 de los 20 hijos siguen siendo eternos (χ_land de 418 y 463 caen).
- **Se espera POLO_EQUIVALENTE.**

## Lo que U22 no resuelve

- **La conjetura de la distancia de la Ciénaga** no se evita: Δχ ~ 10³ M_Pl en el espacio de campos. Se gana protección radiativa, no el criterio de la distancia.
- **La regularización del polo cerca de χ = 0** cambia m_χ, y con ella el recalentamiento. Se informa su efecto.
- **ζ₀,max de C45j** se toma tal como se calculó con la suma: el potencial apenas interviene en el rebote, pero no se rehace.
