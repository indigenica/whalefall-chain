# U17e, U23 y U20b · La población del colapso medida, la superradiancia y la isocurvatura observacional: la banda de F queda en 0.09–0.5·F₀

*17 de octubre de 2026 (Madrid).*
- *Preregistros:* `docs/preregistros/ESPEC_U17e.md`, `ESPEC_U23.md` (enmiendas 1–3) y `ESPEC_U20b.md`, todos anteriores a su ejecución.
- *Código:* `lya/u17e_class.py`, `chain/u23_superradiance.py`, `chain/u20b_iso.py`.
- *Resultados:* `lya/results/u17e/U17E_RESULT.json`, `results/u23_superradiance.json`, `results/u20b_iso.json`.
- *Ejecución:* en la nube, 27 filas de CLASS en ~3 h; U23 y U20b en segundos.

## 0. En una frase

| Prueba | Veredicto preregistrado | Qué dice |
|---|---|---|
| **U17e** (población medida en U19b, por CLASS) | **BANDA_RECORTADA** (F_max = 0.35·F₀), como se esperaba | Lyα solo deja pasar hasta 0.7·F₀; el corte en 0.35 lo pone el criterio de isocurvatura de U20. |
| **U23** (superradiancia, con autointeracción) | **SUPERRADIANCIA_PARTE_LA_BANDA** (se esperaba «recorta el borde inferior») | Los agujeros negros estelares excluyen F = 0.0625 y 0.045·F₀. La fila de 0.03·F₀ cae solo en la región menos robusta, y por la regla escrita la banda queda partida. |
| **U20b** (isocurvatura con las cotas observacionales) | **ISOCURVATURA_NO_LIMITA**, como se esperaba | Con f_iso < 0.0084 (Lyα, GG 2026), toda la rejilla pasa, incluido F₀. El criterio de U20 era unas 500 veces más estricto que los datos. |

**La mejor estimación combinada** (filas Y: Lyα de U17e, isocurvatura de U20b, superradiancia estelar de U23, ΔN_eff y cuerdas del CMB) es **F ∈ [0.09, 0.5]·F₀**:
- f_a ≈ 0.5–2.7×10¹⁴ GeV, m_G ≈ 6×10⁻¹⁷–8×10⁻¹⁴ eV;
- el colapso es el 98 % de la materia oscura, que es fría (δA ≤ 0.15);
- más una fila suelta en 0.03·F₀, en la zona probable de superradiancia;
- **F = 0.7·F₀** pasa Lyα (por forma) y la isocurvatura observacional, pero cae bajo la cota de supermasivos (m_G < 1.7×10⁻¹⁷ eV), que es poco robusta.

## 1. U17e: la población medida, pasada por CLASS

**La forma.** El exceso blando de las cuatro redes de U19b, con el desalineamiento quitado de los bins de k bajo (la atribución conservadora), tiene un momento medio x̄ = 5.83·m_G·a_osc. La log-normal de U17d tenía 3: la población medida es más caliente.

**El número.** En las filas Y, K = (8π/27)·R_ρ·ξ_B2·ℓ·Y/x̄ ≈ 38–47 por cuanto de desalineamiento; en las filas X, el número medido sin escalar, K = 1.12.

| F/F₀ | Y: m_G (eV) | Y: δA | Y: Lyα | Y: isocurvatura de U20 (r = 2.4) | X: m_G (eV) | X: δA |
|---|---|---|---|---|---|---|
| 1 | 2.9×10⁻¹⁸ | 0.468 | falla | 25 | 2.1×10⁻¹⁵ | 0.017 |
| 0.7 | 1.35×10⁻¹⁷ | 0.283 | pasa por forma | 2.6 | 9×10⁻¹⁵ | frío trivial |
| 0.5 | 6.4×10⁻¹⁷ | 0.154 | pasa | 0.27 | 3.5×10⁻¹⁴ | frío trivial |
| **0.35** | 3.1×10⁻¹⁶ | 0.077 | pasa | **0.027** | 1.4×10⁻¹³ | frío trivial |
| 0.25 | 1.2×10⁻¹⁵ | 0.039 | pasa | 0.003 | 5.7×10⁻¹³ | frío trivial |
| 0.18 | 4.8×10⁻¹⁵ | 0.019 | pasa | ~0 | 2.2×10⁻¹² | frío trivial |
| ≤ 0.125 | ≥ 2×10⁻¹⁴ | ≤ 0.005 | pasa | ~0 | ≥ 10⁻¹¹ | frío trivial |

- **Control de frío trivial:** la fila Y con F = 0.125 da δA = 0.005 < 0.05 con CLASS; pasa.
- **Banda de U17e** (Lyα y la isocurvatura de U20): Y, 0.03–0.35; Y+, 0.03–0.35; X, toda la rejilla; la rama cómoda (A), ≤ 0.125.
- **Sensibilidades:** con el desalineamiento quitado en proporción (Yp), δA baja un poco (0.168 y 0.083 en 0.5 y 0.35), pero la isocurvatura de U20, con su k_* algo más bajo, las tumba. Con los cuantos duros de la red tal como salen (Xd), F₀ sigue pasando (δA = 0.057).
- **La expectativa preregistrada** (BANDA_RECORTADA con F_max = 0.35 en Y; X pasa en toda la rejilla) **se cumplió**.

## 2. U23: superradiancia de agujeros negros

Puerta primaria: la región excluida a 2σ de Witte y Mummery (2025), digitalizada, con f_eff entre F/4 y F.

| Filas | Excluidas por agujeros estelares | Probables (región menos robusta) | Excluidas por supermasivos (secundaria, poco robusta) |
|---|---|---|---|
| Y | 0.0625, 0.045 | 0.03 | 1, 0.7 (ya fuera por Lyα o isocurvatura) |
| Y+ | 0.045, 0.03 | — | 1, 0.7 |
| X | 0.25, 0.18 | 0.125 | ninguna |
| P y H de U21 | 0.0625, 0.045 | 0.03 | 1, 0.7 |

- **Veredicto: SUPERRADIANCIA_PARTE_LA_BANDA.** La fila de 0.03·F₀ (m_G = 9×10⁻¹² eV) no cuenta como excluida porque cae en la región menos robusta (m > 2), y por la regla escrita la banda queda partida. Se esperaba «recorta el borde inferior»; **la expectativa no se cumplió por esa fila.** Que esa fila sobreviva no tiene peso físico: es el borde de la región de exclusión.
- **Supermasivos** (Ünal, Pacucci y Loeb): excluyen m_G < 1.7×10⁻¹⁷ eV. En las filas Y eso es F ≥ 0.7·F₀.
- **Filas X:** la ventana estelar cae en medio (0.18–0.25·F₀). Si se confirma la lectura (b) de U19b, la banda serían dos tramos: F ≥ 0.35·F₀, F₀ incluido, y F ≤ 0.09·F₀.

## 3. U20b: la isocurvatura con las cotas de los datos

f_iso = C·(k_CMB/k_*)³/A_s, con C = 1, r = 2.4 y la cota de Lyα de Garcia-Gallego et al. (2026), f_iso < 0.0084:

| F/F₀ (filas Y) | 1 | 0.7 | 0.5 | 0.35 | ≤ 0.25 |
|---|---|---|---|---|---|
| f_iso | 4.2×10⁻³ | 4.4×10⁻⁴ | 4.7×10⁻⁵ | 4.7×10⁻⁶ | < 10⁻⁶ |
| Lyα, GG 2026 (0.0084) | pasa | pasa | pasa | pasa | pasa |
| Lyα anterior (4×10⁻³) | falla por poco | pasa | pasa | pasa | pasa |
| Enanas ultradébiles (6×10⁻⁴) | falla | pasa | pasa | pasa | pasa |

- **Veredicto: ISOCURVATURA_NO_LIMITA.** La isocurvatura deja de recortar la banda.
- **El corte en 0.35·F₀ de U17e era un artefacto del criterio de U20,** no de los datos. U20 y U17e siguen con sus veredictos preregistrados; U20b los sustituye como mejor estimación.
- **La detección tentativa de GG** (f_iso ≈ 0.005–0.0076) no encaja con ninguna fila con C = 1. F₀ da 0.0042, justo debajo, y además falla Lyα. Encajaría con C ≈ 1.5 en F₀, pero F₀ falla Lyα en las filas Y.

## 4. Qué cambia en el modelo

1. **La banda de la rama viva** pasa de «F ≲ 0.35·F₀» a **F ≈ 0.09–0.5·F₀** (f_a ≈ 0.5–2.7×10¹⁴ GeV).
   - El techo lo ponen Lyα (δA = 0.154 en 0.5·F₀; 0.7·F₀ pasa solo por forma) y la cota de supermasivos.
   - El suelo lo pone la superradiancia estelar.
2. **La materia oscura de la ballena es fría,** con m_G ≈ 6×10⁻¹⁷–8×10⁻¹⁴ eV. Casi toda (98 %) viene del colapso cuerdas–paredes, más caliente que la log-normal postulada, pero fría en toda la banda.
3. **La lectura (b) de U19b,** que se decide con las redes de ℓ_G = 2.4, cambiaría mucho el cuadro: F₀ pasaría todas las puertas (m_G = 2×10⁻¹⁵ eV, δA = 0.017) y **la señal de cuerdas en LISA volvería a ser una predicción.**
4. **La rama cómoda (ξ creciente, B4a)** solo vive con F ≲ 0.125·F₀ con el criterio de U20, y hasta 0.25 con el de U20b. Sigue forzada.

## 5. Lo que no resuelven

- La extrapolación del colapso de ℓ₁ ≈ 4.3 a ℓ ≈ 90 (filas Y frente a X); la decide U19b con ℓ_G = 2.4.
- La normalización C de la isocurvatura sigue tomándose como 1. Medir el espectro de densidad en una red lo cerraría.
- La región de superradiancia es la de una figura digitalizada, con ±0.1 dex. La cota de supermasivos es la de un solo trabajo, considerado poco robusto por otros dos.
- La verosimilitud completa (U15c, en WSL) sigue en curso, con los casos 0.5·F₀ y 0.35·F₀: ahora caen justo dentro de la banda combinada.

## Nota · 18 de octubre de 2026 (Madrid): la lectura (b) de U19b no se confirmó

U19b con ℓ_G = 2.4 da la lectura (c), mixta (`INFORME_U19b.md` §7 y su corrección).
- **Las filas X y la vuelta de F₀ y de LISA** (§2, punto 3 de §4) quedan como variante no apoyada por los datos.
- **El X del universo sería ≈ 170–320,** el entorno de las filas Y (275).
- **La banda combinada de §0, F ≈ 0.09–0.5·F₀,** sigue siendo la mejor estimación.
- **Errata:** en §1, «K = 1.12» en las filas X es 1.11 (1.115 en `U17E_RESULT.json`).

## Corrección · 18 de octubre de 2026, 21:46 (Madrid), tras una segunda auditoría independiente: hay dos bandas, y la de 0.5·F₀ lleva la cota de supermasivos

La auditoría reprodujo las 27 etiquetas de superradiancia, los f_iso de U20b y las filas de U17e. Los veredictos no cambian. Pero §0 y §4 presentaban como «mejor estimación combinada» una banda que no es la combinada preregistrada:

- **La banda combinada preregistrada (ESPEC_U20b)** es **F ∈ [0.09, 0.7]·F₀**, más la fila de 0.03·F₀ (`results/u20b_iso.json`, `bands.U17e_Y.combined`). Combina:
  - la puerta primaria de Lyα;
  - la isocurvatura observacional (U20b);
  - la superradiancia estelar no excluida (U23);
  - ΔN_eff;
  - las cuerdas del CMB.

  La cota de supermasivos, en U20b, solo **marca** la fila de 0.7.
- **Con la cota secundaria de supermasivos** (U23, enmienda 3), la banda es **[0.09, 0.5]·F₀**.
  - La fila de 0.7 tiene m_G = 1.35×10⁻¹⁷ eV, frente al umbral de 1.7×10⁻¹⁷ eV. El margen es de 0.1 dex, del orden de la incertidumbre de la digitalización.
  - La cota viene de un solo trabajo, que otros dos consideran poco robusto.
- **§4, punto 1, atribuía el techo a Lyα. Es incorrecto:** Lyα no corta entre 0.5 y 0.7·F₀. La fila de 0.7 pasa la puerta primaria por forma y solo falla el umbral estricto de δA (0.283 > 0.234). Entre 0.5 y 0.7 corta únicamente la cota de supermasivos.
- **Cómo hay que decirlo:** banda combinada **0.09–0.7·F₀**, o **0.09–0.5·F₀** si se acepta la cota de supermasivos. Con 0.7, m_G baja hasta 1.35×10⁻¹⁷ eV y δA llega a 0.28, que pasa solo por forma.

**Detalles menores:**
- «fría (δA ≤ 0.15)» es δA ≤ 0.154 en 0.5·F₀.
- La rejilla termina en 0.03·F₀. Por debajo no se ha evaluado, así que puede haber un segundo tramo abierto.
- El suelo de 0.09 es conservador: interpolando, el borde estelar está en ≈ 0.075·F₀.
- En la tabla de §3, las marcas de UFD con F ≥ 0.7 caen fuera de su rango de validez (k_* < 1000 Mpc⁻¹).
- En la digitalización de Witte y Mummery hay una franja espuria (log m ≈ −12.74 a −12.67, trazos de texto). Ninguna fila cae en ella, pero conviene enmascararla.
