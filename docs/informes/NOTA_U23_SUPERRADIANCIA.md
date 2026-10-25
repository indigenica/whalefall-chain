# Nota U23 · Superradiancia en la banda de F: lectura de la literatura (sin veredicto preregistrado)

*16 de octubre de 2026 (Madrid). Es una lectura de la literatura, no un cálculo. Desde la nube no se pudieron leer los artículos completos (arXiv y sus espejos están bloqueados); solo los resúmenes y fragmentos del buscador. Cuando se puedan leer, U23 debe hacerse como cálculo preregistrado.*

## La pregunta

U21 marcó que, **sin autointeracción**, las bandas excluidas por superradiancia de agujeros negros cubrirían casi toda la banda de F:
- supermasivos, cota indirecta, ~10⁻¹⁹–10⁻¹³ eV;
- estelares, 6×10⁻¹³–10⁻¹¹ eV (Cardoso et al. 2018).

Los pseudo-Goldstones de la ballena tienen autointeracción cuártica λ ~ m_G²/f², con f = F/2 ≈ 2.7×10¹⁴ GeV en el borde de la banda (F = 0.5·F₀) y menos por debajo. ¿Siguen valiendo esas cotas?

## Lo que dice la literatura (fragmentos verificados por el buscador)

- **Baryakhtar, Galanis, Lasenby y Simon (PRD 103, 095019; arXiv:2011.11646):** con autointeracción grande, el crecimiento superradiante se satura en un estado de equilibrio con ocupación reducida. La extracción de giro queda suprimida paramétricamente, y no hay «bosenova» en la mayor parte del espacio astrofísico. Los umbrales en f por debajo de los cuales el giro ya no se frena son **~3×10¹⁶ GeV y ~2×10¹⁷ GeV** para m = 8×10⁻¹³ y 2.5×10⁻¹² eV, agujeros estelares.
- **Fermilab, S-stars (2025):** para f por debajo de ~**10¹⁶ GeV**, la autointeracción detiene el crecimiento antes de que el agujero negro supermasivo pierda giro.
- **«Stepping up superradiance constraints on axions» (PRD 111, 083044; arXiv:2412.03655):** un reanálisis bayesiano dice que la sensibilidad en f estaba subestimada en un orden de magnitud. Aun así, el umbral quedaría en ~10¹⁵ GeV.

## Lectura para la ballena

- **En toda la banda de U20/U21,** f ≤ 2.7×10¹⁴ GeV: entre uno y tres órdenes de magnitud por debajo de esos umbrales. **Por la literatura, la superradiancia no excluiría la banda:** la autointeracción saturaría las nubes antes de frenar el giro de los agujeros negros.
- **El benchmark F₀** (f = 5.35×10¹⁴ GeV) también quedaría por debajo.

**Matices:**
1. Los umbrales citados son para otras masas. Su dependencia en α y en la masa del agujero negro no se ha comprobado aquí.
2. Los Goldstones de RP³ son tres campos con la autointeracción del sigma no lineal, no un ALP con potencial coseno. La estructura de λ es parecida, pero no idéntica.
3. **Las señales de ondas continuas** («transitions») pueden sobrevivir con interacción moderada (Baryakhtar et al.). No excluyen, pero podrían ser una predicción.

**Conclusión provisional:** la superradiancia probablemente no recorta la banda. **Pendiente:** U23 como cálculo preregistrado, con la tasa de saturación de Baryakhtar et al. aplicada a nuestros (m_G, f) y a agujeros negros concretos (M33 X-7, IRAS 09149-6206, M87*), cuando se puedan leer los artículos.

## Fuentes

- https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.095019 (arXiv:2011.11646)
- https://inspirehep.net/files/f997527fa5ae22cbec923fb868eb8f4e
- https://lss.fnal.gov/archive/2025/pub/fermilab-pub-25-0150-t.pdf
- https://arxiv.org/html/2412.03655 (PRD 111, 083044)
- https://arxiv.org/pdf/1801.01420 (Cardoso et al. 2018)

## Corrección · 17 de octubre de 2026 (Madrid): la conclusión de esta nota era incorrecta

Con los artículos completos, que aportó el equipo, se vio que la lectura de arriba estaba mal. Esta nota no se reescribe; queda superada por `docs/preregistros/ESPEC_U23.md` y su informe.
- **Los umbrales «~3×10¹⁶ y ~2×10¹⁷ GeV»** no son umbrales de apagado. En BGLS (arXiv:2011.11646, Fig. 9 y 10, ec. 62), con agujeros negros estelares el frenado del giro sigue ocurriendo hasta f ~ 10¹²–10¹³ GeV. La región excluida por Cyg X-1, GRS 1915+105 y otros cubre m ≈ 2×10⁻¹³–6×10⁻¹² eV para **todo f ≳ 10¹³ GeV**. WM (arXiv:2412.03655, Fig. 1) la refuerza.
- **Para la ballena:** las filas de la banda con m_G en esa ventana (F ≲ 0.06·F₀ en la rama viva) **sí quedan excluidas**. Las de F mayor tienen m_G por debajo de la ventana estelar.
- **Lo de los agujeros negros supermasivos** sí era correcto en lo esencial: la autointeracción los apaga para f ≲ 10¹⁵ GeV (BGLS, Fig. 11), y sus cotas no son robustas.
