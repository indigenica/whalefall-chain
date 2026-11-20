# U62d · ¿Excluye el CMB a los hijos visibles del cuscutón? El espectro desde el giro de la madre, proyectado a C_ℓ · Preregistro

*2026-11-16 18:29 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U62d»). Escrito antes del código de U62d. Interno: no se exporta.*

## 1. La pregunta

¿Son compatibles con la potencia del CMB en ℓ = 2–30 los hijos de cuscutón que dejan su rebote a la vista (U61b)?
- U62c encontró excesos de potencia de ×1.7–18 en ℓ = 2–5, con el estado A2, pero quedó INCONCLUSO.
- El motivo: el fondo integrado hacia atrás hasta el giro de la madre perdía la ligadura, así que el estado A1 no era fiable.

## 2. Comprobación previa (declarada)

`results/u62d_runs/previo/viabilidad.{py,out}`, solo fondo, sin perturbaciones:
- **el fondo integrado hacia delante desde el giro de la madre** (el estado inicial de U61b: χ = 0, χ̇ = √(2fρ_t), φ = −φ∞) conserva la ligadura:
  - residuo ≤ 7×10⁻¹¹ en la contracción;
  - ≤ 4×10⁻¹³ en la ventana;
  - ≤ 9×10⁻⁸ en la expansión, hasta ρ_r/ρ = 10⁻⁵.
- **El quinto hijo** (A·H = 1 000, f = 0.898) acaba en la hoja χ < 0.

## 3. El montaje

**Los hijos:** los cuatro visibles de U61b en la hoja χ > 0. El de la hoja χ < 0 se excluye: según el supuesto de medida del preprint, allí el Higgs es taquiónico y no hay observadores. Se informa solo con lo de U62c.

**El fondo:** hacia delante desde el giro, sin la trayectoria de U59, en tres tramos:
- el modo −1 hasta φ = −0.2;
- el modo 0 hasta φ = +0.2;
- el modo +1 hasta ρ_r/ρ = 0.9×10⁻⁵.

Las matrices son las de U59, con el fondo de cada tramo (u62_diagnostico.bg_mode).

**Los estados:** el vacío instantáneo del marco sin B_s (U62b), como en U62c.
- **A1:** a Δln a = 0.1 del giro (allí H ≈ 0). Si E′ no es positiva, se acerca al rebote en pasos de 0.1 e-folds.
- **A2:** a Δln a = 2 antes de la ventana.
- **B:** a Δln a = 2 después de la ventana. Es informativo.

**La propagación y la referencia:** como en U62c, con el código de U62c importado:
- el cambio de etapa en ρ_r/ρ = 10⁻⁴;
- la acción de Boruah, Kim y Geshnizjani hasta k/aH = 10⁻³;
- el atractor de U58 en |χ_×|.

**La rejilla:** por hijo, 100 n espaciados logarítmicamente entre n_lo = 0.25·n(ℓ = 2) y n_hi = mín(máx(2.5·n(ℓ = 30), 3·n_mem), 3×10⁴), con n(ℓ) = ℓ/sen(D√|Ω_K,0|). Más un control n_c = 10·n_mem con A2.

**La proyección a C_ℓ.**
- **La proyección** es de Sachs–Wolfe en el límite de ángulo pequeño, con el ángulo del último scattering sobre S³, χ_LS = D√|Ω_K,0|, entre 2×10⁻³ y 2×10⁻²:
  - C_ℓ[P] ∝ Σ_{n=2}^{n_max} (1/n)·P(n)·j_ℓ²(n·χ_LS);
  - se suma sobre los n enteros, con n_max = 50·n(ℓ = 30).
- **Los ingredientes:** P(n) = F(n)·P_att(n).
  - F se interpola linealmente en (ln n, ln F) sobre la rejilla, y vale 1 fuera de [n_lo, n_hi];
  - P_att se interpola en (ln n, ln P), y fuera de la rejilla se extiende como ley de potencias.
- **El cociente:** R_ℓ = C_ℓ[F·P_att]/C_ℓ[P_att], para ℓ = 2–30, con A1 y con A2.

**La comparación con el CMB:**
- R̃_ℓ = 1 + w·(R_ℓ − 1), con w = 0.6 para el veredicto. Es una dilución conservadora por los términos que no son de Sachs–Wolfe (ISW, Doppler). Se informa también w = 1.
- Δχ² = f_sky·Σ_{ℓ=2}^{30} (2ℓ+1)·[1/R̃_ℓ + ln R̃_ℓ − 1], con f_sky = 0.8.
- Es la verosimilitud de Wishart con el cielo ΛCDM como dato. No usa el déficit observado en los ℓ bajos, que haría la exclusión más fuerte para los excesos.

## 4. Puertas (si falla alguna, U62d es INCONCLUSO)

- **G0, la normalización:** ≤ 10⁻¹⁰.
- **G1, el wronskiano al final:** ≤ 10⁻³.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en tres n del primer hijo. P_A1 cambia ≤ 10⁻³.
- **G3, el control:** |F_tot(n_c) − 1| ≤ 0.03 con A2.
- **G4, el cambio de etapa:** con 10⁻⁵ en lugar de 10⁻⁴, P_A1 cambia ≤ 1 % en los tres n.
- **G5, el fondo:** el residuo de Friedmann es ≤ 3×10⁻⁵ en todos los tramos.
- **G6, el muestreo:** con la rejilla a la mitad (los puntos pares), Δχ² cambia ≤ 20 % relativo, o ≤ 2 absoluto, en cada hijo y estado.

## 5. Veredicto

Un hijo queda **excluido** si Δχ² ≥ 25 con A1 y con A2, y **compatible** si Δχ² < 9 con los dos.
- **INCONCLUSO:** falla una puerta.
- **VISIBLES_EXCLUIDOS:** los cuatro hijos de la hoja χ > 0, excluidos.
- **VISIBLES_COMPATIBLES:** los cuatro, compatibles.
- **MIXTO:** cualquier otro caso.

**Se informan siempre:**
- F_tot(n) para A1, A2 y B;
- R_ℓ;
- Δχ² con w = 0.6 y con w = 1;
- la dependencia del estado;
- la fracción del núcleo de Sachs–Wolfe de C_2 que viene de n < n_lo.

## 6. Lo que no cubre

- **El cálculo completo de C_ℓ:** CLASS con la geometría cerrada y el espectro externo, con el ISW tardío.
- **Los datos reales de Planck:** su espectro de los ℓ bajos y su verosimilitud.
- **La polarización.**
- **Otros estados de la madre,** antes del giro.
- **El hijo de la hoja χ < 0.**

## 7. Lo que se espera

- **G0–G6 pasan.**
- **Con A·H = 1 000, excluidos.**
- **Con A·H = 30–100, no está claro.**
- **Veredicto: MIXTO o VISIBLES_EXCLUIDOS,** con confianza baja en cuál de los dos.

## 8. Salidas

- **Código:** `chain/u62d_cl_rebote.py`.
- **Resultados:** `results/u62d_cl_rebote.json` y `results/u62d_runs/`.
- **Informe:** sección U62d de `docs/informes/INFORME_U62.md`.
- Interno.
