# U1 · Cierre cosmológico unificado de WHALEFALL (ballena 5 + campañas de cierre)

*6 de octubre de 2026 · este repositorio*
- *Preregistro:* `ESPEC_CIERRE.md` (06:55 UTC, antes del código; erratas de la auditoría al final).
- *Resultados:* `results/`. *Figura:* `figures/u1_cierre.png`. *Contrato:* `results/CONTRATO_CIERRE.json`.
- *Auditoría independiente:* `audit/u1/AUDIT_U1.md`.

## 0. En una frase

**Con un único conjunto de parámetros, la cadena entera, desde la semilla hasta la materia oscura, se cierra de forma condicional (veredicto preregistrado CIERRE_ESTRECHO).** Hay una región no vacía en la que se cumplen a la vez:
- Ω_DM h² = 0.12;
- ΔN_eff (CMB y BBN);
- Lyman-α (materia oscura ondulatoria y flujo libre);
- la superradiancia.

La condición es que **la red RP³ emita en Goldstones libres como mucho entre el 0.1 % y el 6 % de lo que emitiría una red de cuerdas U(1) con la normalización asintótica de GHN.** Ese número, ε, es ahora la cantidad que tiene que medir una simulación representativa de la red. Por encima del 1–6 %, la componente emitida es demasiado tibia y Lyman-α excluye el modelo.

**Esto no es una cosmología cerrada:** la fuente sigue NOT_CALIBRATED. Es la primera vez que la cadena completa se evalúa contra los datos con una sola incógnita explícita.

## 1. Qué se ha unido

| Etapa | Recalculado aquí | Resultado | Estado |
|---|---|---|---|
| S1 Nacimiento (C45, ballena 5) | 90 hijos con la copia de `frw_core.c`: los 20 viables de C45j y 70 de control | Reproducción exacta (diferencia 0, mismos rebotes). Con el V0 común, χ_et pasa de 417.28 a **415.58**, y los 20 viables siguen siendo eternos y sub-Planck | PASS |
| S2 Bolsillo y fase eterna | Ajuste del pivote continuo (`thermal.py`, copia) | n_s = 0.97289, r = 0.00502, N\* = 51.55, T_rh = 1.98×10⁹ GeV, F = 1.07×10¹⁵ GeV; diferencia ≤ 3×10⁻¹⁰ con el archivo | PASS |
| S2b D3 | Cota de planitud √(V0/3)/2π = 1.20×10⁻⁶ /M_Pl | La meseta sin deformar la cumple: en χ = 10³ está 100 veces por debajo del límite. Cualquier corrección mayor la rompe (C45i FRÁGIL). UV: lectura (b) | requisito |
| S3 Sectores | Lectura de especies, desintegraciones y fluido tardío | T_tardía/T_rh − 1 = 2.5×10⁻⁷; F idéntica; ΔN_eff homogéneo 3×10⁻⁹ | PASS |
| S4 Libro de materia oscura | Versión vectorizada de `dark.py`, más flujo libre | Igual a la copia del original a 1×10⁻¹⁵ por paquete; raíces de `dark.json` a 1.2×10⁻⁷ dex; horizonte 14 145 Mpc | PASS |
| S5 Cierre | 309 nodos (espectro × x0 × K_sw × ε) | 178 conservadores, 112 estrictos | **CIERRE_ESTRECHO** |

**Unidad de parámetros.** S1 a S5 usan el mismo potencial (la clase `Potential` de cadena_cierre, de 10 niveles), el mismo V0, la misma F y la misma Γ. Solo hay una excepción declarada: los barridos de nacimiento de C45 se hicieron con V0_BENCH = 1.67×10⁻¹⁰, un 2.8 % menos. S1 los reclasifica con el V0 común y nada cambia.

## 2. El cierre (S5)

**Montaje.** Ω_mis(m_G; K_sw) + Ω_em(m_G; ε, espectro, x0) = 0.12. Se busca m_G para:
- 17 eficiencias ε entre 10⁻⁴ y 1, más ε = 0;
- dos espectros de emisión (delta infrarroja en x0·H, y energía plana en log p);
- x0 ∈ {10, 30, 100};
- K_sw ∈ {0, 1, 10} (horquilla de muros de C35b).

**Lo que decide es la tibieza de lo emitido (P5), no ΔN_eff ni Lyα.**
- Los Goldstones emitidos cerca del colapso de la red (T_stop ≈ 0.4–2.7 MeV) nacen con p ≈ 10–100 H ≈ 3–30 m_G. Recorren hasta ~0.4 Mpc comóviles antes de enfriarse.
- La referencia es un WDM térmico de 5.3 keV, con λ_fs = 0.072 Mpc calculada con el mismo integrador.
- En las 17 configuraciones con borde, el primer nodo que falla lo hace solo por P5 (auditoría, 2b).

**Eficiencia máxima permitida** (contigua desde ε = 0, en la malla de 0.25 dex, K_sw = 0):

| Espectro | x0 = 10 | x0 = 30 | x0 = 100 |
|---|---|---|---|
| ir_delta | 5.6×10⁻³ | 5.6×10⁻³ | 1.0×10⁻³ |
| energy_logflat | 1.8×10⁻² | 1.8×10⁻² | 5.6×10⁻² (\*) |

(\*) Este borde lo fija la franquicia f_em ≤ 0.1, no el flujo libre. Sin ella, el máximo con K_sw = 0 es 1.8×10⁻².
- Interpolando entre nodos, el rango es ~1.3×10⁻³ a 9.5×10⁻².
- Con WDM de 3.5 keV (historia térmica no suave), ε_max sube hasta 10⁻²–10⁻¹. Con 6 keV baja hasta 5.6×10⁻³ en el caso logflat x0 = 10 (`results/s5_sensitivity_wdm.json`).

**Lectura física.** Para ε ≪ 1, la materia oscura es casi toda desalineamiento (m_G ≈ 10⁻¹⁴ eV, como en C35b). Al subir ε:
1. la emisión toma el relevo;
2. la masa necesaria baja;
3. la red colapsa más tarde;
4. lo emitido sale más caliente.

La cadena solo se cierra si la red RP³ es **mucho menos eficiente que una red U(1)** convirtiendo energía en Goldstones libres. Las campañas de defectos dan una pista compatible: el flujo neto medido incluye mucha energía que se queda ligada o vuelve al núcleo (`cierre_red`, bloqueadores). Pero **no es una medida.**

## 3. Predicciones dentro de la región permitida (K_sw ∈ {0, 1})

| Observable | Predicción | Comparación |
|---|---|---|
| n_s, r | 0.9729, 0.0050 | Planck 2018 a +1.9σ; ACT DR6 (P-ACT) a −0.4σ; r < 0.036 (BK18). Un CMB-S4 con σ(r) ~ 10⁻³ lo vería a ~5σ |
| m_G | **1.6×10⁻¹⁶ – 9.6×10⁻¹⁵ eV** | En el hueco entre la ventana de superradiancia de agujeros negros supermasivos (≤ 1.3×10⁻¹⁶, inactiva para f_a < 2×10¹⁶ GeV) y la estelar (≥ 3.8×10⁻¹⁴). El extremo inferior está a solo un factor 1.2 del borde SMBH |
| ΔN_eff | **≤ 7×10⁻³** (≤ 2×10⁻³ con K_sw = 0) | **Indetectable para CMB-S4.** Corrige la ficha anterior («0.002–0.15»): con masa, emisión finita y ε ≲ 10⁻², el momento de radiación es mucho menor |
| Fracción emitida f_em | 0 – 0.85 | Parte de la materia oscura puede ser «tibia ligera», con λ_fs hasta ~0.4 Mpc |
| Minicúmulos | 0.24 – 115 M_⊙ | Escala de Fairbairn (C35b). Microlente: abierto |
| Ondas gravitacionales de la red | f_a = 5.35×10¹⁴ GeV = 1.07 × umbral LISA de GHN | En el umbral, sin SNR calculada |
| T_osc | 0.4 – 2.7 MeV | Coincide con la BBN, pero la densidad de energía del campo es entonces despreciable (ΔN_eff de BBN ≤ 7×10⁻³). Es también cuando colapsa la red |

## 4. Qué sigue abierto (estado heredado, sin maquillar)

1. **ε (eslabón 8): la incógnita decisiva.** El paso siguiente es medir la eficiencia de emisión libre de una red RP³ representativa, con núcleos resueltos y AMR bidireccional (`docs/evidencia/cierre_cosmologia_20261005/SIGUIENTE_PASO.md`), en la normalización de `s4_dark.source`.
   - Si sale ε ≲ 10⁻², la ballena se cierra.
   - Si sale ε ~ 1, **queda falsada** en su versión actual: materia oscura demasiado tibia.
2. **m_G** es un término de ruptura explícita añadido. En la acción común (`cierre_accion`), los Goldstone no tienen masa.
3. **D3:** la meseta tiene que ser plana a 10⁻⁶ hasta χ ~ 10³ M_Pl, fuera de la ciénaga (`docs/literatura/NOTA_D3_UV.md`).
4. **D11:** coherencia de varios rebotes, MARGINAL (3.6 %). La prueba con campo (C45k) está montada.
5. **Semilla (D1), medida (D9), islas (D12):** como en la ballena 5.
6. **Método:** el momento 3P no es una verosimilitud de N_eff, y P5 es un indicador de orden de magnitud. La auditoría lo acepta, con un factor de pocos en el borde. Faltan Boltzmann completo, isocurvatura y la SNR de LISA.

## 5. Reproducir

```
# desde la raíz del repositorio
make chain        # ~2 min; S1 necesita gcc (WSL). Sin gcc, S1 = SKIPPED
```

Cola WSL: `make chain`.

---

## 6. Actualización (6 de octubre, tarde): la ε medida en la red (U2, `docs/informes/INFORME_U2.md`)

**La incógnita de U1 ya está medida, en cociente.** Simulación 3D del campo Q de la acción común (vacío RP³, cuerdas ℤ₂ con núcleo físico, era de radiación), frente a una U(1) de control con la misma f y la misma masa de núcleo, en el mismo código:
- **R_Γ = Γ_RP³/Γ_U1 = 1.69 ± 0.29**, para ℓ ≈ 3–5. La red RP³ emite más que la U(1), no menos.
- 20 redes de hasta 256³; validaciones V2–V4 aprobadas; auditoría independiente sin errores de normalización.

**Veredicto preregistrado: INCONCLUSO.** El control U(1) da ξ = 0.45 en ℓ = 4.5, justo por debajo del rango de la literatura (≥ 0.5): las redes aún no han llegado al atractor de escalado.

**Lo que indica.** Al ℓ físico (≈ 50–85), según la extrapolación:

| Rama | ε |
|---|---|
| GHV (ξ ≈ 0.24ℓ) | ≈ 1.7 |
| Hindmarsh et al. (ξ\* = 1.2) | ≈ 0.1–0.3 |
| Tendencia medida | ≈ 0.05–0.3 |

- **Con las puertas preregistradas de §2, las tres quedan excluidas incluso a −1σ.** La materia oscura emitida sale demasiado tibia.
- **Con una cota de Lyα para materia oscura mixta fría + tibia** (Garcia-Gallego et al. 2025; exploratorio) sobreviven la rama de Hindmarsh y la de la tendencia, siempre con el espectro log-plano. La predicción pasa a ser:
  - m_G ≈ 2–7×10⁻¹⁵ eV;
  - un 15–54 % de la materia oscura tibia, con un λ_fs como el de un WDM de 2–4 keV;
  - ΔN_eff ≈ 0.002–0.006.
- **La región «ε ≪ 1» de §2–§3 queda descartada por la medida** si se confirma en el atractor.

**Estado:** el sector oscuro de la ballena 5 queda en **tensión fuerte**. Depende de tres cosas:
1. que la literatura resuelva a favor del escalado con ξ constante;
2. que el espectro sea log-plano;
3. que Lyα admita una fracción tibia.

**Siguiente prueba:**
- redes que lleguen al atractor (preparación con cuerdas gruesas y N ≥ 384);
- la puerta Lyα mixta preregistrada, con Boltzmann.

---

## 7. Actualización (6 de octubre, noche): U3 (`docs/informes/INFORME_U3.md`)

- **La red no llega al atractor con N = 256** (U3b, ni con preparación de núcleo). La medida R_Γ = 1.69 ± 0.29 sigue en INCONCLUSO.
- **Lyα con transferencia de Boltzmann** (U3a, que sustituye al indicador λ_fs de §2). Con un espectro log-plano, las ramas de escalado constante (B2) y de la tendencia medida (B3) pasan; la de crecimiento logarítmico (B1) no.
  - Predicción en ese caso: m_G ≈ 1.7–6.9×10⁻¹⁵ eV, un 15–58 % de materia oscura tibia y δA ≈ 0.09–0.22.
  - Es una lectura condicional a la red, y el resolvedor no está validado contra CLASS (W1 falla).
