# U20b · La isocurvatura con las cotas observacionales de f_iso, en vez del criterio aproximado de U20 · Preregistro

*17 de octubre de 2026, ≈20:11 (Madrid). Escrito antes de ejecutar U20b (la hora que cuenta es la de su commit), con U17e todavía corriendo. Se ejecuta en la nube, en segundos.*

**Lo que ya se ha visto, y se declara:**
- La lectura de Chathirathas y Schwetz (C&S, arXiv:2511.15790, §4 y §7) y de Garcia-Gallego et al. (GG, arXiv:2603.04401), aportados por el equipo el 17 de octubre.
- La cuenta de una línea de abajo (f_iso en F₀ para la variante P de U20). No se ha calculado ninguna fila de U20b.

## Motivo: dos cosas que la lectura cambia

**1. La «cota estricta de isocurvatura» de U21 no era de isocurvatura.** La cifra m ≳ 10⁻¹⁷ eV del resumen de C&S es una combinación. Para ALPs postinflacionarios con masa constante, la isocurvatura sola da una cota débil, y el 1.7×10⁻¹⁷ eV lo pone la **superradiancia de agujeros negros supermasivos** (Ünal, Pacucci y Loeb, arXiv:2012.12790; C&S §7). Esa cota se trata en U23, como información.

**2. El criterio de U20 era mucho más estricto que los datos.** U20 pedía Δ²_iso/Δ²_ad ≤ 0.1 en k = 10 h/Mpc y z = 3. Las cotas observacionales están en términos de f_iso, el cociente entre la potencia primordial de isocurvatura y la adiabática en k_CMB = 0.05 Mpc⁻¹, con Δ²_S(k) = f_iso·A_s·(k/k_CMB)³ (GG, ec. 1; C&S, ec. 4.11):
- **Lyα, GG 2026:** f_iso < 0.0084 a 95 %, con el modelo de ruido conservador. Con el modelo de ruido publicado hay una detección tentativa: f_iso = 0.0064 (+0.0012/−0.0014).
- **Lyα, la cota anterior** usada por C&S: f_iso ≤ 4×10⁻³.
- **Galaxias enanas ultradébiles (UFD):** f_iso ≲ 6×10⁻⁴ en k = 10–1000 Mpc⁻¹. Está discutida por los efectos de marea.

En F₀ (variante P de U20, C = 1, r = 1), f_iso = (k_CMB/k_*)³/A_s ≈ 4×10⁻⁴. U20 la daba como fallo por un factor 22, y con GG pasa por un factor 20: los dos criterios difieren en unas 500 veces.

## Qué se calcula (`chain/u20b_iso.py`)

Para cada fila de U17e (Y, Y+ y X) y de U21 (P y H):
- f_iso = C·(k_CMB/k_*)³/A_s, con k_* = a_osc·H_osc/r como en U20;
- A_s = 2.10×10⁻⁹;
- C ∈ {0.3, 1, 3, 10}, con **C = 1 como primaria**;
- **r = r_N = 2.40** como primaria en las filas de U17e (la época de U19b); r = 1 y 3 en P y H, como en U20; r_W = 1.78, informativa.

**Validez:** las cotas se aplican si su intervalo de k queda por debajo de k_* (ruido blanco). Para Lyα, k ≲ 10 Mpc⁻¹; para UFD, hasta 1000 Mpc⁻¹. Si no, la fila se marca.

**Puertas:**
- **Primaria:** GG, f_iso < 0.0084.
- **Informativas:** Lyα de C&S (4×10⁻³), UFD (6×10⁻⁴) y la compatibilidad con la detección tentativa de GG (f_iso entre 0.005 y 0.0076).

## Veredicto (filas Y de U17e, C = 1, r = r_N)

F_max,iso es el F más alto de la rejilla que pasa la puerta primaria.
- **ISOCURVATURA_NO_LIMITA:** pasa toda la rejilla, F₀ incluido.
- **ISOCURVATURA_LIMITA:** F_max,iso < F₀. Se da su valor.

**Banda combinada** (se informa; es la mejor estimación actual):
- puerta primaria de Lyα de U17e;
- isocurvatura de U20b;
- superradiancia estelar de U23 (fila no «excluida»);
- ΔN_eff y cuerdas del CMB.

Se dice aparte qué filas de la banda combinada caen bajo la cota informativa de supermasivos (m < 1.7×10⁻¹⁷ eV) y bajo la de UFD.

## Lo que se espera, escrito antes de calcular

- **En F₀,** con r = 2.4 y C = 1, f_iso ≈ 5×10⁻³ < 0.0084. **Se espera ISOCURVATURA_NO_LIMITA.**
- **El techo de la banda combinada** pasaría a ser Lyα (U17e). Por encima de 0.5·F₀ quedaría además la cota informativa de supermasivos.
- **La detección tentativa de GG** solo encaja cerca de F₀ con r ≈ 2.4.

## Lo que U20b no resuelve

- **La normalización C** sigue sin medirse. La medida directa sería el espectro de densidad al final de una red de U19b; queda propuesta.
- **Lo no lineal en escalas pequeñas** (minicúmulos) no se trata.
- **Los veredictos de U20 y U21 no se cambian:** se calcularon con sus criterios preregistrados, y U20b los sustituye como mejor estimación.
