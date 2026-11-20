# U58 · Lo que el cuscutón deja en la inflación y en el universo tardío · Preregistro

*15 de noviembre de 2026, 21:18 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U58»). Es la prueba B de `paper/zenodo_whale/PLAN_V4_CUSCUTON.md`. Escrito antes del código y de cualquier espectro de U58. La hora que cuenta es la de su commit. Interno: no se exporta.*

## 1. Las preguntas

1. **¿Qué cota observacional hay sobre ε_G?** El cuscutón con el potencial de Boruah cambia la ecuación de Friedmann a tiempos tardíos, 3H²(1 + ε_G) = ρ, con ε_G = 3μ⁴/(8φ∞²). La pregunta es qué juegos de parámetros la cumplen.
2. **¿Cuánto mueve el cuscutón n_s y r** durante la inflación del bolsillo, calculado con su acción cuadrática y no con una estimación?

## 2. Lo que ya se sabe, y se declara

- **La cota, leída hoy en el cuerpo de los artículos** (nota final de `INFORME_U57.md`):
  - Afshordi (2009, arXiv:0907.5201), ec. 13: ε_G = (3/2)(1 − λ);
  - ec. 14: λ − 1 = 0.003 ± 0.014 al 95 %, tomada de Robbers, Afshordi y Doran (2008, PRL 100, 111101), que es una medida de dos lados del desajuste entre la masa de Planck dentro y fuera del horizonte.
  - En ε_G: **−0.025 ≤ ε_G ≤ 0.017**. Esta parte de U58 no es ciega: la cuenta ya está hecha.
- **La acción cuadrática,** leída hoy en Boruah, Kim y Geshnizjani (2017, arXiv:1704.01131), ecs. 3.16–3.24:
  - S = ∫ d⁴x z²[ζ′² − c_s²k²ζ²], en tiempo conforme, en la gauge δχ = 0 y con ζ la curvatura comóvil con χ;
  - z² = a²α(k² + 3α𝓗²)/(k² + α𝓗²(3 − σ));
  - c_s² = (k⁴ + k²𝓗²B₁ + 𝓗⁴B₂)/(k⁴ + k²𝓗²A₁ + 𝓗⁴A₂);
  - con α = χ̇²/(2H²), ε = −Ḣ/H² y σ = ε − α (negativo con nuestro signo, −μ² en la acción).

  El texto extraído del PDF tiene mal B₁ y B₂ (se pierden símbolos). Uso B₁ = A₁ + σ(6 + η + β − 2ε) + α(η − β) y B₂ = A₂ + σα(12 − 4σ + 3η) + 3α²(η − β), con η = ε′/(𝓗ε) y β = α′/(𝓗α). **Por eso G2 comprueba que c_s² no importa en el cruce del horizonte.** z² no tiene esa duda: se obtiene de la factorización de su ec. 4.4.
- **Una estimación de servilleta, de la nota:**
  - en slow-roll con 3H²(1 + ε_G) = V: α = ε_V(1 + ε_G)², P_ζ ∝ V/[ε_V(1 + ε_G)³] y P_t ∝ V/(1 + ε_G);
  - a N_* fijo, r ∝ (1 + ε_G)^(2/3), es decir **Δr/r ≈ +⅔ε_G**, y **n_s no cambia a primer orden.**
- Ningún espectro con cuscutón calculado.

## 3. El montaje

**El fondo.** FLRW plano en la inflación (la curvatura y la radiación del hijo ya son despreciables), con la meseta de polo de la cadena (V₀ y C de `benchmark()`) y el cuscutón de Boruah en su atractor tardío (φ → φ∞):
- **H y φ** salen de resolver a la vez U′(φ) = 3μ²H y 3H² = ½χ̇² + V + U(φ);
- **φ̇** = 3μ²·½χ̇²/D, con D = (3/2)μ⁴ − U″;
- **Ḣ** = −½χ̇² + ½μ²φ̇;
- **χ** sigue χ̈ + 3Hχ̇ + V′ = 0;
- se integra en e-folds desde χ = 6 hasta el final (ε = 1).

**Los juegos:**
- **GR**, sin cuscutón (ε_G = 0, el control);
- **«pequeño»**: μ = 0.16, m = 0.02, φ∞ = 0.25 (ε_G = 0.0039);
- **C1**: μ = 0.3, m = 0.05, φ∞ = 0.25 (ε_G = 0.049).

**El pivote:** k_* = aH cuando faltan N_* = 51.59 e-folds (`benchmark()`) para el final.

**Los escalares:**
- se integra la acción de BKG (z² y c_s²) para k_*·e^(−0.5), k_* y k_*·e^(+0.5);
- vacío de Bunch–Davies con v = √2·zζ, puesto cuando k = 100·aH;
- hasta k = aH/1 000;
- P_ζ = k³|ζ|²/(2π²) y n_s − 1 por diferencia centrada en ln k.

**Los tensores:** como en relatividad general, con la masa de Planck desnuda (el cuscutón no tiene tensores). Con v = ah/2, P_t = 4k³|v|²/(π²a²) y r = P_t/P_ζ en k_*.

**Informativo:** el ΔN_eff equivalente en BBN, −(43/7)·ε_G, por la cuenta de g_* = 10.75.

## 4. Puertas (si falla alguna, U58 es INCONCLUSO)

- **G0, el control:** con GR, n_s a ≤ 0.002 del de la cadena (0.97277) y r a ≤ 10 % (0.005016). Son las tolerancias de la puerta de equivalencia CMB de S2. Además, P_ζ(k_*) a ≤ 5 % de H²/(8π²α) en el cruce.
- **G1, la convergencia:** con rtol de 10⁻⁹ y de 10⁻¹⁰, n_s y r cambian a ≤ 10⁻⁴ relativo, en los tres juegos.
- **G2, c_s²:** con c_s² de BKG y con c_s² = 1, P_ζ(k_*) cambia a ≤ 10⁻⁴ en los dos juegos con cuscutón.
- **G3, el fondo:**
  - residuo de Friedmann ≤ 10⁻¹⁰ relativo;
  - D > 0 y φ̇ > 0 en toda la integración;
  - z² > 0 en todos los modos.

## 5. Veredicto

- Un juego está **permitido por los datos** si −0.025 ≤ ε_G ≤ 0.017.
- **INCONCLUSO:** falla una puerta.
- **CUSCUTON_COMPATIBLE:** algún juego con cuscutón está permitido por los datos, tiene hijos en U57b y cumple |Δn_s| ≤ 0.001 frente a GR.
- **CUSCUTON_EXCLUIDO:** ningún juego con hijos está permitido por los datos.
- **CUSCUTON_TENSO:** el juego permitido mueve n_s más de 0.001.

**Se informan siempre:**
- Δn_s y Δr/r de cada juego, frente a la servilleta (⅔ε_G);
- el ΔN_eff equivalente;
- la posición de cada juego frente a la cota.

## 6. Lo que no cubre

- El rebote y la perturbación fuera de la inflación (U59).
- La curvatura del hijo durante la inflación.
- **El contraste de un Ω_Q < 0** («early dark energy» negativa) con datos posteriores a 2008: la cota usada es la de Robbers et al. (2008); una actualización con Planck sería un trabajo aparte.
- La inestabilidad ultravioleta con λ < 1 que menciona Afshordi (2009).

## 7. Lo que se espera

- **G0–G3 pasan.**
- **Δr/r** ≈ +0.26 % con el juego pequeño y ≈ +3.2 % con el de C1; |Δn_s| < 10⁻⁴.
- **CUSCUTON_COMPATIBLE** con el juego pequeño. El de C1 queda excluido por la cota.

## 8. Coste y salidas

- **Coste:** minutos.
- **Código:** `chain/u58_cuscuton_inflacion.py`.
- **Resultados:** `results/u58_cuscuton_inflacion.json`.
- **Informe:** `docs/informes/INFORME_U58.md`.
- Interno.
