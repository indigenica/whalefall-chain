# U23 · Superradiancia de agujeros negros sobre la banda de F, con autointeracción · Preregistro

*17 de octubre de 2026, ≈20:02 (Madrid). Escrito antes de ejecutar U23 (la hora que cuenta es la de su commit). Se ejecuta en la nube, en segundos, cuando termine U17e.*

**Lo que ya se ha visto, y se declara:**
- Los artículos que aportó el equipo el 17 de octubre:
  - Baryakhtar, Galanis, Lasenby y Simon (BGLS), arXiv:2011.11646;
  - Witte y Mummery (WM), arXiv:2412.03655;
  - Mehta et al., arXiv:2011.08693;
  - Stott, arXiv:2009.07206;
  - Cardoso et al., arXiv:1801.01420;
  - Bai et al. (Fermilab), arXiv:2507.07482.
- Al leer las figuras se vio que la ventana estelar cae sobre el extremo de F bajo de las filas de U17e (en su ejecución en seco) y de U21. No se ha calculado ninguna puerta fila a fila.

## Corrección previa: la nota U23 del 16 de octubre estaba mal

La nota (`docs/informes/NOTA_U23_SUPERRADIANCIA.md`) decía que la autointeracción apaga la superradiancia por debajo de f ~ 10¹⁵–10¹⁷ GeV para agujeros negros estelares. Venía de fragmentos de un buscador, sin leer los artículos, y **es incorrecto**.
- **BGLS** (Fig. 10, ec. 62) y **WM** (Fig. 1) muestran que, con agujeros negros estelares (Cyg X-1, GRS 1915+105, M33 X-7 y otros), **la exclusión por frenado del giro sigue en pie hasta f ~ 10¹²–10¹³ GeV** con m ≈ 2×10⁻¹³–6×10⁻¹² eV. Cuanto mayor es f, más fuerte es la exclusión, no más débil.
- **La autointeracción sí apaga la superradiancia de los agujeros negros supermasivos** para f ≲ 10¹⁵ GeV (BGLS, Fig. 11; WM §V; Mehta et al.), y esas cotas no se consideran robustas.
- La nota recibe una corrección fechada; no se reescribe.

## Qué se calcula (`chain/u23_superradiance.py`)

Para cada fila de U17e (Y, Y+ y X) y de U21 (P y H), el par (m_G, f_eff).

**La constante de la autointeracción.** El potencial del Goldstone de RP³ es V = ε(Λ − μ)·sin²(χ/F_c), es decir, un coseno en 2χ/F_c. Eso da f_eff = F_c/2, la f_a = F/2 que ya usa la cadena.
- Los tres Goldstones de RP³ tienen además interacciones derivativas del modelo sigma. Se anulan para una nube con polarización fija, pero no en general.
- Por eso f_eff se toma con una incertidumbre de O(1): **f_eff ∈ [F/4, F]**.

**Puerta primaria: agujeros negros estelares.** Es la región excluida a 2σ de WM (Fig. 1, n ≤ 5, m ≤ 2; Cyg X-1 conservador y GRS 1915+105), digitalizada del PDF a 300 ppp en `data/literatura/witte_mummery_2412.03655_fig1.json`.
- La frontera superior se suaviza con una mediana móvil de 9 columnas y se recorta en 1/f = 10^−11.4 GeV⁻¹, para quitar los trazos de texto.
- Una fila **queda excluida** si m_G cae en el intervalo de masas de la región oscura y 1/f_eff queda por debajo de la frontera **en todo** el intervalo f_eff ∈ [F/4, F].
- **Frontera:** solo en una parte del intervalo.
- **Probable:** fuera de la región oscura, pero dentro de la región clara (m > 2, menos robusta).

**Informativo: agujeros negros supermasivos.** Sin veredicto, porque BGLS y WM no las consideran robustas.
- La ventana sin autointeracción de Stott: 5.5×10⁻²⁰ ≤ m ≤ 1.3×10⁻¹⁶ eV.
- El umbral de apagado de BGLS (ec. 62), f_CD, para un agujero negro de 10⁷ M_⊙ con a* = 0.9 y T = t_Edd = 4×10⁸ yr, con α = G·M·m.
- La fila se marca «SMBH posible» si m está en la ventana y f_eff(F) > f_CD.

**No se aplica:** Bai et al. (estrellas S) necesita un acoplamiento cuadrático del axión a fotones, y los Goldstones de la ballena no lo tienen.

## Veredicto (sobre la banda de las filas Y de U17e)

- **SUPERRADIANCIA_NO_RECORTA:** ninguna fila de la banda queda excluida.
- **SUPERRADIANCIA_RECORTA_EL_BORDE_INFERIOR:** solo quedan excluidas filas del extremo de F bajo, y la banda sigue contigua con el mismo F_max.
- **SUPERRADIANCIA_PARTE_LA_BANDA:** alguna fila excluida queda dentro de la banda, que deja de ser contigua.
- **SUPERRADIANCIA_CIERRA:** todas las filas de la banda quedan excluidas.

Se informa lo mismo para Y+, X (que solo se aplica si U19b da la lectura (b)), P y H.

## Lo que se espera, escrito antes de calcular

- **Filas Y:** quedan excluidas F = 0.0625 y 0.045·F₀ (m_G ≈ 4×10⁻¹³ y 1.5×10⁻¹² eV). F = 0.03·F₀ (9×10⁻¹² eV) sale probable. La banda quedaría en 0.09–0.35·F₀. **Se espera SUPERRADIANCIA_RECORTA_EL_BORDE_INFERIOR.**
- **Filas X:** F = 0.25 y 0.18·F₀ quedan excluidas, y 0.125 sale probable. La banda X se partiría: F ≥ 0.35 por un lado y F ≤ 0.09 por otro.

## Lo que U23 no resuelve

- La región se lee de una figura, con una precisión de unas ±0.1 dex.
- No se rehace el cálculo de BGLS ni el de WM, ni se usan agujeros negros concretos.
- Las ondas gravitacionales y las ondas de axión de las nubes saturadas no se tratan aquí; podrían ser una predicción, no una exclusión.

## Enmienda 1 · 17 de octubre de 2026, ≈20:03 (Madrid): la estimación informativa de los supermasivos, antes de ver las filas de U17e

Al probar el código con las filas de U21 (las de U17e todavía no existen) se vio que el agujero negro fijo de 10⁷ M_⊙ da α = G·M·m ≈ 4 para m = 5×10⁻¹⁷ eV. Ahí la ec. 62 de BGLS no vale, y el nivel 211 ni siquiera es superradiante.
- **Cambio:** se toma la masa del agujero negro que da α = 0.1, el crecimiento más rápido de 211, limitada a [10⁶, 10¹⁰] M_⊙.
- Si con esa masa α > 0.45, la fila no tiene superradiancia de supermasivos.
- Es solo la parte informativa; la puerta estelar y el veredicto no cambian.
- **Visto antes de la enmienda:** la clasificación estelar de las filas P y H de U21. Coincide con lo esperado: F = 0.0625 y 0.045 quedan excluidas y 0.03 sale probable. Con la regla escrita, eso parte la banda de P y H, porque la fila de 0.03 no cuenta como excluida.

## Enmienda 2 · 17 de octubre de 2026, ≈20:05 (Madrid): una cota de supermasivos más, informativa, antes de ver las filas de U17e

Chathirathas y Schwetz (arXiv:2511.15790, §2 y §7) usan las regiones de Ünal, Pacucci y Loeb (arXiv:2012.12790), que incluyen la autointeracción. Para ALPs postinflacionarios con masa constante, como los Goldstones de la ballena, concluyen **m ≳ 1.7×10⁻¹⁷ eV**.
- Eso contradice en parte la lectura de BGLS de que la autointeracción apaga la superradiancia de los supermasivos para f ≲ 10¹⁵ GeV.
- No se ha podido leer Ünal et al. Se pide al equipo.
- **Se añade como columna informativa** (`SMBH_CS_excluded`: m_G < 1.7×10⁻¹⁷ eV), sin veredicto. Si Ünal et al. confirma que la región cubre nuestros f, se convertirá en puerta en una enmienda posterior, antes de usarla.
- En las filas de U21 afecta solo a F ≥ 0.7·F₀, que la isocurvatura de U20 ya había excluido.

## Enmienda 3 · 17 de octubre de 2026, ≈20:12 (Madrid): Ünal et al. confirma la región de los supermasivos; pasa a puerta secundaria, antes de ver las filas de U17e

El equipo aportó Ünal, Pacucci y Loeb (arXiv:2012.12790). Su figura 6 excluye, con autointeracción, todo m entre ~3×10⁻²¹ y 1.7×10⁻¹⁷ eV para f por encima de una frontera:
- entre 10⁻¹⁹ y 1.7×10⁻¹⁷ eV, la frontera está entre 2×10¹⁴ y 3×10¹³ GeV;
- por debajo de 10⁻¹⁹ eV, entre 10¹⁵ y 5×10¹⁵ GeV.

Las filas con m < 1.7×10⁻¹⁷ eV tienen f_a ≥ 3.7×10¹⁴ GeV, y todas caen dentro. Como dejó escrito la enmienda 2:
- **se añade una puerta secundaria:** una fila queda excluida por supermasivos si m_G < 1.7×10⁻¹⁷ eV;
- **dos veredictos:** el preregistrado (solo estelar) y uno «con supermasivos», con la misma lógica sobre la unión de las dos exclusiones;
- **matiz:** BGLS y WM consideran estas cotas poco robustas, por el entorno de los núcleos galácticos y los espines inferidos;
- en U20b esta exclusión sigue siendo una marca, como dice su preregistro.
