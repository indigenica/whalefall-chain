# INFORME U58 · Lo que el cuscutón deja en la inflación y en el universo tardío

**Preregistro:** `docs/preregistros/ESPEC_U58.md` (15 nov 2026, 21:18, Madrid).
**Código:** `chain/u58_cuscuton_inflacion.py`. **Resultados:** `results/u58_cuscuton_inflacion.json`. Tardó 36 s en la nube.

Interno: no se exporta.

---

## 2026-11-15 21:21 (Madrid)

### 1. Veredicto preregistrado: INCONCLUSO (falla G2, solo con el juego de C1)

| puerta | resultado |
|---|---|
| G0, el control con GR | **pasa:** n_s = 0.97343 (la cadena da 0.97277, Δ = 6.6×10⁻⁴); r = 0.004877 (−2.8 %); P_ζ/P_slow-roll = 1.019 |
| G1, la convergencia (rtol 10⁻⁹ frente a 10⁻¹⁰) | **pasa:** n_s y r cambian ≤ 10⁻⁹ relativo |
| G2, c_s² de BKG frente a c_s² = 1 (≤ 10⁻⁴ en P_ζ) | **falla:** 2.2×10⁻⁵ con el juego pequeño (pasa) y **2.7×10⁻⁴ con el de C1** (falla) |
| G3, el fondo | **pasa:** Friedmann ≤ 1.6×10⁻¹⁴; D, φ̇ y z² > 0 |

G2 existía porque B₁ y B₂ de BKG salieron mal en la extracción del PDF. Con el juego de C1, la parte de c_s² que depende de esos coeficientes mueve P_ζ un 0.03 %. Es pequeño para cualquier observable, pero supera el umbral preregistrado. El juego de C1 ya está excluido por la cota, y el pequeño pasa G2, pero la regla era sobre los dos.

### 2. Lo que sale (informativo con el veredicto INCONCLUSO)

| juego | ε_G | λ − 1 | ¿permitido? (−0.025 ≤ ε_G ≤ 0.017) | Δn_s | Δr/r | servilleta ⅔ε_G | ΔN_eff equivalente en BBN |
|---|---|---|---|---|---|---|---|
| pequeño | 0.0039 | −0.0026 | **sí**, a 0.8σ del valor central | +1.7×10⁻⁶ | **+0.26 %** | +0.26 % | −0.024 |
| C1 | 0.049 | −0.032 | **no**, a unas 5σ | +2.0×10⁻⁵ | +3.2 % | +3.2 % | −0.30 |

- **La servilleta de la nota se confirma con la acción cuadrática:** el cuscutón no mueve n_s (Δn_s ~ 10⁻⁶–10⁻⁵) y sube r en ≈ ⅔ε_G. El cálculo y la servilleta coinciden al 2 %.
- **En la inflación, σ/ε = −ε_G** exactamente (−0.0039 y −0.0487 en el pivote), como se deduce del atractor φ∞ − φ ∝ H.
- **Con las puertas en verde** habría sido **CUSCUTON_COMPATIBLE**: el juego pequeño tiene hijos (U57b), está dentro de la cota y deja n_s intacto.

**La predicción que deja el cuscutón con este potencial:** G_cos/G_N = 1/(1 + ε_G) < 1, es decir λ < 1. **Una medida de λ > 1, o de G_cos > G_N, lo falsaría.** El valor de ε_G es libre entre m²/(2φ∞²), que exige el rebote (D(0) > 0), y 0.017, que exigen los datos.

### 3. Frente a lo que se esperaba (§7 del ESPEC)

- **«G0–G3 pasan»: fallado en G2,** por el juego de C1.
- **«Δr/r ≈ +0.26 % y +3.2 %; |Δn_s| < 10⁻⁴»: acertado.**
- **«CUSCUTON_COMPATIBLE con el juego pequeño»:** sería el resultado con las puertas en verde.

### 4. Siguiente

- **Si se quiere un veredicto,** una U58b solo con los juegos que puedan entrar en la v4 (el pequeño u otro que elija el autor), con la lectura de B₁ y B₂ comprobada en la fuente LaTeX de arXiv.
- **U59** (los escalares a través del rebote) se haría con el juego que se adopte. Las matrices de C2–C3e del explorador son del juego de C1, que los datos excluyen.

---

## Nota · 2026-11-15 21:22 (Madrid): B₁ y B₂ comprobados en la fuente LaTeX

He descargado la fuente LaTeX del artículo de arXiv (`cuscutonPerturbation.tex`, líneas 270–280). Las definiciones de z², c_s², η, β, A₁, A₂, B₁ y B₂ son **exactamente** las que usa el código:
- η = ε′/(𝓗ε) y β = α′/(𝓗α);
- B₁ = A₁ + σ(6 + η + β − 2ε) + α(η − β);
- B₂ = A₂ + σα(12 − 4σ + 3η) + 3α²(η − β).

**Consecuencias:**
- **La duda que justificaba G2 desaparece.** Los 2.7×10⁻⁴ que separan c_s² de BKG de c_s² = 1 con el juego de C1 no son incertidumbre de transcripción, sino la contribución real de c_s² en el cruce del horizonte. El resultado correcto es el que la incluye (las filas «main»).
- **El veredicto preregistrado no cambia:** sigue siendo INCONCLUSO, porque la regla de G2 estaba escrita así.
- **Lo informativo de §2 queda más firme:** el juego pequeño da Δn_s = 1.7×10⁻⁶ y Δr/r = +0.26 %, con c_s² ya verificado.

---

## U58b · 2026-11-15 22:21 (Madrid): el veredicto con el juego adoptado

**Preregistro:** `docs/preregistros/ESPEC_U58b.md`. No es ciega: lo declara. **Código:** `chain/u58b_cuscuton_inflacion.py`. **Resultados:** `results/u58b_cuscuton_inflacion.json`.

**Veredicto preregistrado: CUSCUTON_COMPATIBLE.** Pasan las cuatro puertas: G2 da 2.2×10⁻⁵ y G1 ≤ 10⁻⁹.

| | juego adoptado (μ = 0.16, m = 0.02, φ∞ = 0.25) |
|---|---|
| ε_G | 0.0039 (λ − 1 = −0.0026, dentro de la cota a 0.8σ) |
| Δn_s frente a GR | +1.7×10⁻⁶ |
| Δr/r | +0.257 % (servilleta: +0.262 %) |
| ΔN_eff equivalente en BBN | −0.024 |
| control de amplitud, V₀ → V₀(1 + ε_G)³ | P_ζ/P_ζ(GR) = 1.005; n_s y r no cambian (Δn_s = 1.7×10⁻⁶, Δr/r = +0.257 %) |

**Para la v4:** con el juego adoptado, el bolsillo predice el mismo n_s que con NPQT y un r un 0.26 % mayor. La única huella tardía es G_cos/G_N = 1/(1 + ε_G) = 0.996.
