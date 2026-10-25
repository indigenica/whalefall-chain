# U22 · La inflación de polo reproduce el CMB, pero no el nacimiento: el parche D2 no se puede quitar sin pagar en S1

*15 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U22.md` (commit `71bcaa6`, anterior a la ejecución).
- *Código:* `chain/u22_pole.py`. *Resultado:* `results/u22_pole.json`. En la nube, en ~2 min.

## 0. En una frase

**Veredicto preregistrado: POLO_CAMBIA_NACIMIENTO.**
- En el CMB, la inflación de polo de orden 4 es equivalente a la suma de 10 niveles: Δn_s = −1.2×10⁻⁴ y Δr = −0.2 %.
- En el nacimiento, solo 10 de los 20 hijos viables siguen siendo eternos, cuando el umbral preregistrado era 15.

## 1. CMB (S2): equivalente

Mismo flujo que S2: final exacto en ε_H = 1, recalentamiento continuo con σ/m = 1e-5 y empalme del pivote.

| | N_* | n_s | r | T_rh (GeV) | F (GeV) | χ_et (M_Pl) |
|---|---|---|---|---|---|---|
| suma de 10 niveles (S2) | 51.55 | 0.97289 | 0.005025 | 1.98×10⁹ | 1.070×10¹⁵ | 416 |
| polo, C = 0.3516 | 51.59 | 0.97277 | 0.005016 | 2.29×10⁹ | 1.069×10¹⁵ | 541 |
| polo, C × 0.8 | 51.56 | 0.97288 | 0.004313 | 2.45×10⁹ | 0.99×10¹⁵ | 504 |
| polo, C × 1.25 | 51.63 | 0.97265 | 0.005839 | 2.13×10⁹ | 1.15×10¹⁵ | 580 |

- **n_s no depende de la forma concreta de la cola.**
- **r va con C,** porque C fija la amplitud de la cola: con C de ±25 %, r cambia un ±15 %. El C primario sale de la forma de la suma, no de los datos.

## 2. Nacimiento (S1): cambia

Son las 90 ejecuciones de S1 con el puerto C de C45 en el modo 4 (polo). Los modos 0–3 no cambian, y S1 sigue reproduciendo bit a bit.

| | eternos con la suma | eternos con el polo | fértiles con el polo |
|---|---|---|---|
| 20 hijos viables de C45j | 20 | **10** | 20 |
| 70 controles | ~35 (χ_land > 416) | **7** | 64 |

- **χ_land cambia mucho y sin un patrón sencillo:** la mediana de |ln(χ_land,polo/χ_land,suma)| es 0.68, y la media de ln es −1.16, un factor ~3 a la baja. Los números de rebotes casi no cambian (8 de 90).
- **El umbral eterno sube de 416 a 541 M_Pl.**
- **Por qué.** El aterrizaje depende de cómo χ atraviesa el potencial a χ ≲ 10–10³ tras el salto, y ahí la suma y el polo difieren: m_χ² = 4.6·V₀ con la suma y 8.1·V₀ con el polo, y la suma tiene ondulaciones entre niveles. **El CMB solo ve χ ≈ 4; el nacimiento ve todo el potencial.**

## 3. Qué significa

- **Lo bueno.** La meseta del CMB sí sale de un lagrangiano estándar (polo de orden 4) con las mismas n_s y r. Para S2, D2 deja de ser un parche.
- **Lo que no esperábamos.**
  - La eficiencia del nacimiento depende de la forma del potencial lejos de la ventana del CMB, que nada mide.
  - Con el polo, los nacimientos eternos son unas 5 veces más raros entre los fértiles de control.
  - **El eslabón S1 es más frágil de lo que parecía: su resultado depende de un detalle del potencial que el modelo no fija.**
- **La expectativa preregistrada** (~18 de 20, por el cambio de umbral) **no se cumplió.** Suponía que χ_land no cambiaría, y sí cambia.
- **Lectura honesta.** Hay que elegir:
  - (a) la suma de 10 niveles, que es un parche pero da el nacimiento calculado;
  - (b) el polo, que es un lagrangiano estándar pero rehace las estadísticas del nacimiento (C45h/C45j) con nacimientos eternos más raros.

  Con (b), el nacimiento sigue existiendo (10 de los 20 hijos y 7 de los 70 controles siguen siendo eternos), pero la fracción viable baja de ~3.6 % a algo del orden de ~1 %. Hay que recalcularla con la estadística de C45j: el criterio ζ₀ se calculó con la suma.

## 4. Lo que U22 no resuelve

- **La fracción viable con el polo** no se ha recalculado: haría falta rehacer C45h/C45j con el potencial nuevo, que cabe en la nube.
- **La conjetura de la distancia de la Ciénaga** no se evita: Δχ ~ 10³ M_Pl.
- **La regularización del polo cerca de χ = 0** es una elección. Afecta a m_χ, al recalentamiento y, como se ve aquí, al nacimiento.

## Nota 1 · 16 de octubre de 2026 (Madrid): la lectura «5 veces más raros» era errónea (U22b)

En el §2 conté los controles eternos de la suma sin la condición E_phys,max < 1 ni la de fertilidad, que sí apliqué al polo. Con la definición completa de S1 son **9 de 70 con la suma, frente a 7 de 70 con el polo**, no ~35 frente a 7. Fue un error de cuenta en la lectura posterior, no del cálculo.

**U22b** (`docs/informes/INFORME_U22b.md`) lo confirma con toda la rejilla de C45h: la fracción de fértiles que acaban eternos es 0.163 con la suma y 0.153 con el polo (R = 0.94, NACIMIENTO_ROBUSTO).

**Lo que sigue en pie:**
- el veredicto preregistrado de U22 (POLO_CAMBIA_NACIMIENTO), que cuenta hijos concretos;
- y el hecho de que **qué hijo concreto nace eterno sí depende del potencial**: el aterrizaje cambia mucho de un caso a otro.

**Lo que se retira:** la conclusión de que «los nacimientos eternos son ~5 veces más raros» y de que «S1 es más frágil». Estadísticamente, el nacimiento no cambia.
