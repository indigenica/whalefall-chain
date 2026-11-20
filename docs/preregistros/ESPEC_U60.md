# U60 · La anisotropía con cuscutón: un Bianchi IX con Einstein + cuscutón · Preregistro

*2026-11-17 14:24 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U60»). Es la prueba D de `paper/zenodo_whale/PLAN_V4_CUSCUTON.md`, la laguna más importante de la completitud 4D del nacimiento con cuscutón. Escrito antes del código. Interno: no se exporta.*

## 1. La pregunta

¿Sobreviven el rebote con cuscutón y la inflación del hijo a una anisotropía en la madre?
- En una contracción, la cizalla crece frente a la radiación, el problema clásico de los rebotes.
- Con NPQT, U37 encontró que los hijos aguantan anisotropías de forma de 3×10⁻³ a 3×10⁻² (ANISOTROPIA_INOCUA).
- Aquí no hay prescripciones: el cuscutón es una acción, y la ligadura fija la expansión media también con cizalla.

## 2. Las ecuaciones

**El modelo:** un Bianchi IX diagonal, con las formas invariantes de SU(2) normalizadas como en REVENGE-R2 del explorador, de modo que la esfera redonda de radio a tiene ³R = 6/a².
- q_i = ln a_i;
- c_i = 2 exp(q_i − q_j − q_k);
- R_i = ½[c_i² − (c_j − c_k)²], y R3 = ΣR_i.

**La materia:** radiación isótropa (ρ_r = C_rad·e^{−4α}, con α = Σq_i/3), χ en la meseta de polo y el cuscutón (juego adoptado μ = 0.16, m = 0.02, φ∞ = 0.25).

**La cinemática:** H_i = q̇_i, H = ΣH_i/3, s_i = H_i − H (Σs_i = 0) y σ² = ½Σs_i².

**El cuscutón,** con φ(t) homogéneo:
- ρ = U y p = −U − μ²|φ̇|;
- su ecuación da U′(φ) = 3μ²H (la rama φ̇ > 0), independiente de la cizalla;
- derivándola, φ̇ = 3μ²W/D, con D = 1.5μ⁴ − U″ y W = ⅔ρ_r + ½χ̇² + σ² − R3/6. En el límite isótropo es la W de `explora_cuscuton`;
- **la rama se rompe si W ≤ 0.**

**La evolución:**
- q̇_i = H + s_i;
- ṡ_i = −3Hs_i − (R_i − R3/3);
- χ̈ = −3Hχ̇ − V′(χ);
- φ con los modos −1, 0 y +1 de `explora_cuscuton` (Cusc.all).

**La ligadura hamiltoniana,** que se vigila y no se impone en la evolución: 3H² = ρ_r + ½χ̇² + V + U + σ² − R3/2.
- Se proyecta (se reajusta φ) al empezar los tramos de inicio de la aceleración, de fin y de cola, como en U61b.

**Los tramos y el resultado:** como `u61b_cuscuton_huella.host`:
- **los tramos:** half−, half+, deceleración, inicio de la aceleración, fin de la inflación, y la cola hasta χ_end de la cadena;
- **la aceleración** se mide con el factor de escala medio;
- **los eventos terminales añadidos:** W = 0 (RAMA_ROTA), máx q_i − mín q_i > 10 (ANISOTROPO) y una a_i < 10⁻³·a_b,iso (COLAPSO);
- **N_acc** se mide entre el inicio de la aceleración y su fin, como en U61b. **Ω_σ = σ²/(3H²)** se toma en χ_end.

## 3. Los datos

- **Las madres:** los 52 anfitriones de U61b (`results/u61b_cuscuton_huella.json`), con sus (A·H, f).
- **La rodaja del giro de la madre:**
  - en H = 0 y φ = −φ∞, con χ = 0 y s_i = 0;
  - q_i = ln ā + ε·d_i;
  - se conserva la materia, como en la enmienda 1 de U37: C_rad = 3a_t²/(1 + f) (el valor isótropo, con a_t = A·H/H_meseta), ½χ̇² = f·ρ_r(ā), y ā se despeja de la ligadura: ā² = 2(1 + f)C_rad/r3(ε), con R3 = r3(ε)/ā²;
  - si r3 ≤ 0 no hay rodaja, y la amplitud se excluye.
- **Las direcciones** (de Misner):
  - +β₊: d = (1, 1, −2);
  - −β₊: d = (−1, −1, 2);
  - +β₋: d = (√3, −√3, 0).
- **Las amplitudes:** 12 valores de ε, espaciados logarítmicamente de 10⁻⁶ a 10^−0.5, como en U37, más ε = 0.

## 4. Las corridas y su clasificación

**Una corrida es ROBUSTA** si se cumplen las cuatro condiciones:
1. el hijo nace (llega a χ_end tras aterrizar, sin RAMA_ROTA, ANISOTROPO, COLAPSO ni recolapso);
2. |N_acc − N_acc(ε = 0)| ≤ 1;
3. Ω_σ(χ_end) < 10⁻²;
4. el residuo de la ligadura es ≤ 10⁻⁷ (B2).

**ε\* de un hijo** es la mayor amplitud tal que todas las corridas con ε ≤ ε\*, en las tres direcciones, son robustas. Una corrida inválida por B2 corta el intervalo, que es la lectura conservadora.

**El hijo es robusto** si ε\* ≥ ζ₀ = 4.6×10⁻⁵, como en U37.

## 5. Puertas (si falla alguna, U60 es INCONCLUSO)

- **B0, el límite isótropo** (ε = 0 en los 52): χ_land y N_acc frente al `host` de U61b, a ≤ 10⁻⁵ relativo y ≤ 10⁻⁴ absoluto.
- **B1, el modo de forma lineal:**
  - con ε = 10⁻⁷ en +β₋, en los tres primeros hijos, la deformación δ = (q₁ − q₂)/(2√3) se compara con la ecuación lineal δ̈ + 3Hδ̇ + 8δ/a² = 0 (la de R2), integrada sobre el fondo isótropo, desde el giro hasta φ = 0;
  - máx|δ − δ_lin|/máx|δ_lin| ≤ 10⁻⁴.
- **B2, la ligadura:** se comprueba en cada corrida. Si más del 20 % de las corridas con ε ≤ 10⁻² son inválidas, U60 es INCONCLUSO.
- **B3, la tolerancia:** rtol 10⁻¹¹ frente a 10⁻¹² en 8 corridas (dos hijos, dos amplitudes y dos direcciones):
  - la misma clasificación;
  - |ΔN_acc| ≤ 10⁻³;
  - |ΔΩ_σ| ≤ 10⁻⁴.

## 6. Veredicto

Con f_rob, la fracción de los hijos robustos:
- **INCONCLUSO:** falla una puerta.
- **ANISOTROPIA_INOCUA:** f_rob ≥ 0.9.
- **ANISOTROPIA_FATAL:** f_rob ≤ 0.1.
- **ANISOTROPIA_ACOTADA:** cualquier otro caso.

**Se informan siempre:**
- la distribución de ε\* frente a A·H, y su comparación con U37;
- las causas de fallo por encima de ε\*;
- la cizalla en el rebote;
- Ω_σ al final.

## 7. Lo que no cubre

- **Las anisotropías no homogéneas,** que serían los modos de tensor y de vector de U59 y C4.
- **Los Bianchi con inclinación** (radiación que se mueve respecto a la foliación).
- **El estado de la madre antes del giro.**
- **La cizalla inicial con velocidad** (s_i ≠ 0 en la rodaja). Se informa como fila adicional, sin veredicto, en un hijo.

## 8. Lo que se espera

- **B0–B3 pasan.**
- **La servilleta** (NOTA_CUSCUTON_V4 §3.5, con la cizalla como energía rígida) daba una tolerancia de energía ≈ 3×10⁻⁵, es decir, un ritmo de ≈ 5×10⁻³, que empeora con el tamaño de la madre.
- **Con la curvatura de Bianchi IX,** la forma oscila y su energía crece como la radiación, no como a⁻⁶, mientras la curvatura domine sobre la fricción de Hubble. Por eso espero ε\* ~ 10⁻³–10⁻², como en U37.
- **Veredicto: ANISOTROPIA_INOCUA,** con confianza media.

## 9. Salidas

- **Código:** `chain/u60_bianchi_cuscuton.py`.
- **Resultados:** `results/u60_bianchi_cuscuton.json` y `results/u60_runs/`.
- **Informe:** `docs/informes/INFORME_U60.md`.
- Interno.
