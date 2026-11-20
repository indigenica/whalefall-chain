# U63 · El orden cúbico y el acoplamiento fuerte del cuscutón cerca de φ̇ → 0 · Preregistro

*2026-11-17 16:57 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con la acción cúbica cerca de φ̇ → 0»). Es la prueba F del plan de la v4, el punto P02 de la revisión de Codex del explorador, y la prioridad que señala E4c del explorador (`whalefall-birth-explorer` PR #13). Escrito antes del código. Interno: no se exporta.*

## 1. La pregunta

¿Sigue siendo perturbativa la teoría (Einstein + cuscutón + χ + radiación, en S³, con el juego adoptado) a lo largo del nacimiento?
- **Dónde mirar:** en la contracción, en el rebote, y al empezar la inflación, donde φ̇₀ → 0 (φ̇₀/H cae varios órdenes).
- **Qué preguntar:** ¿hasta qué escala son perturbativos los modos que lleva el rebote, los de n ≲ n_mem, en los que U62b–U62e encontraron memoria?
- **La referencia** es Dehghani, Geshnizjani y Quintin (2025, arXiv:2503.01992), en un rebote con cuscutón plano, ecpirótico y sin radiación:
  - débilmente acoplado en el infrarrojo;
  - en el ultravioleta, una escala fuerte Λ ~ (H_b²·M_Pl)^{1/3} a través del rebote, por encima de la escala del fondo.

  Su modelo no es el nuestro, y su éxito no certifica el de la ballena (E4c del explorador lo advierte).

## 2. Comprobación previa (declarada)

`results/u63_runs/previo/viabilidad_L3.py`: el coeficiente ε³ del lagrangiano ADM de U59 (el mismo Ltot, con un armónico zonal de S³) se deriva en 3 s y se convierte en función en otros 0.1 s, con 2 982 operaciones. Queda en caché en `results/u59_runs/L3_symbolic.pkl`. No se evaluó ningún valor físico.

## 3. El montaje

- **Los hijos (8):**
  - los cinco visibles de U61b;
  - los de f más bajo con A·H = 48 y A·H = 300;
  - el de f más alto con A·H = 100.
- **El fondo:** el de U62d, integrado hacia delante desde el giro, desde A1 (0.1 e-fold antes del giro) hasta ρ_r/ρ = 10⁻⁵, ya en inflación. Tiene la ligadura a ≤ 10⁻⁷.
- **Los modos:**
  - n ∈ {2, 10, 30, 100, 300, 10³, 3×10³, 10⁴, 3×10⁴};
  - los dos modos de entrada en el vacío del marco sin B_s, en A1 (como en U62d), propagados con las matrices de U59, guardando la solución densa.
- **Los campos sobre la solución lineal:**
  - q = (ζ, x) y q̇ = K⁻¹(P − Bq);
  - las auxiliares u = (α, β, f) = −S⁻¹(Tᵀq̇ + Rq), la solución lineal de las ligaduras;
  - ḟ por diferencias centradas a lo largo de la solución.
- **La amplitud:** la de una banda logarítmica de vacío. Cada modo se multiplica por n^{3/2}, de modo que su rms espacial sea la del campo en una banda Δln n = 1. Cada cuadratura (parte real e imaginaria del modo de vacío) se evalúa por separado, y se toma la mayor.
- **Tres medidas, en cada instante:**
  1. **ρ₃** = ⟨|L₃|⟩/⟨|L₂|⟩: la densidad cúbica frente a la cuadrática del lagrangiano ADM completo, en valor absoluto y promediadas sobre S³ (∫ dr con el peso de U59), evaluadas sobre la solución lineal;
  2. **Q_u** = |2ḟ/φ̇₀ − 2α − 2fφ̈₀/φ̇₀²|: la perturbación de 2X sobre las hojas de φ constante (el gauge unitario del cuscutón). Es la variable de expansión de la raíz √(2X), invariante de gauge a primer orden;
  3. **Q_t** = k_fís·|f|/φ̇₀: la inclinación de las hojas de φ constante. Q_t ≥ 1 significa que ∂φ deja de ser temporal y la acción del cuscutón deja de estar definida: el fallo no perturbativo.
- **La escala fuerte de un hijo:** n_sc, el menor n de la rejilla en el que máx_t máx(ρ₃, Q_u, Q_t) ≥ 1, o «> 3×10⁴» si no se alcanza. En momento físico en el rebote, k_sc = n_sc/a_b.
- **Las referencias:**
  - k_mem = n_mem/a_b (U61b);
  - H_b, la |H| máxima de la ventana;
  - la de DGQ, (H_b²)^{1/3}.

## 4. Puertas (si falla alguna, U63 es INCONCLUSO)

- **G0, el L₂ de este código frente al L₂ en caché de U59:** en 20 puntos aleatorios de (fondo, perturbaciones, r, Y), relativo ≤ 10⁻¹².
- **G1, la integral de L₂ frente a la forma de U59:** ∫L₂ dr (con el 4π) frente a Σ items·integrales de U59, en 3 hijos × 3 n × 5 instantes, relativo ≤ 10⁻⁸.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en el primer hijo y tres n. El máximo de máx(ρ₃, Q_u, Q_t) cambia ≤ 10⁻³ relativo.
- **G3, la cuadratura en r:** doblar los puntos cambia el máximo de ρ₃ ≤ 10⁻³ relativo.

## 5. Veredicto

- **INCONCLUSO:** falla una puerta.
- **ACOPLAMIENTO_FUERTE:** en algún hijo, n_sc ≤ n_mem, o máx(ρ₃, Q_u, Q_t) ≥ 1 en algún n ≤ n(ℓ = 30).
- **PERTURBATIVO:** en todos los hijos, n_sc ≥ 10·n_mem, y en todos los n ≤ n(ℓ = 30), máx(ρ₃, Q_u, Q_t) ≤ 0.1.
- **LIMITE:** cualquier otro caso.

**Se informan siempre:**
- los máximos de cada medida frente a n y frente al tiempo: contracción, ventana e inflación;
- φ̇₀/H a lo largo del camino;
- k_sc frente a H_b, k_mem y la escala de DGQ;
- en qué fase se alcanza la escala fuerte.

## 6. Lo que no cubre

- **La inflación profunda,** más allá de ρ_r/ρ = 10⁻⁵. Se informa la tendencia de las medidas con φ̇₀/H.
- **El acoplamiento de modos distintos:** las medidas usan un solo armónico, la autointeracción. La integral completa del cúbico con tres armónicos de S³ y el biespectro quedan fuera.
- **Los vectores y los tensores.**
- **Los bucles,** y una escala de corte en sentido estricto (unitariedad de amplitudes).

## 7. Lo que se espera

- **G0–G3 pasan.**
- **En el infrarrojo, perturbativo:** los modos de ℓ ≤ 30 tienen amplitudes ~10⁻⁵.
- **n_sc** entre ~10·n_mem y la escala de DGQ (≈ 0.02 M_Pl en el rebote, frente a k_mem ≈ H_b ≈ 3×10⁻³).
- **Q_t puede crecer al empezar la inflación,** porque φ̇₀ → 0. Si lo hace, también los modos que salen del horizonte serían sensibles: es el riesgo de esta prueba.
- **Veredicto: PERTURBATIVO o LIMITE.** Confianza baja.

## 8. Salidas

- **Código:** `chain/u63_cubico_cuscuton.py`.
- **Resultados:** `results/u63_cubico_cuscuton.json` y `results/u63_runs/`.
- **Informe:** `docs/informes/INFORME_U63.md`.
- Interno.

---

## Enmienda 1 · 2026-11-17 16:58 (Madrid), antes de escribir el código

Cambio de método por coste. No se ha calculado nada de U63.

1. **ρ₃** ya no se promedia sobre todo S³.
   - Para n hasta 3×10⁴, la cuadratura de U59 necesita ~8n puntos por instante, del orden de 10⁹ evaluaciones por fila.
   - ⟨|L|⟩ se toma sobre una longitud de onda local alrededor de r = π/2, donde está la mayor parte del volumen de S³: 64 puntos uniformes en r ∈ [π/2 − π/(n+1), π/2 + π/(n+1)].
   - **G3 pasa a ser:** con 128 puntos en lugar de 64, el máximo de ρ₃ cambia ≤ 10⁻³ relativo.
2. **G1** (la integral de L₂ sobre S³ frente a la forma de U59) se hace con la cuadratura completa de U59, en n = 2, 10 y 100.
3. **Las amplitudes de banda** de Q_u y Q_t usan las rms del armónico normalizado (∫Y² = 1 en el S³ de volumen 2π²):
   - rms(Y) = 1/(π√2);
   - rms(|∇Y|) = √(n(n+2))/(π√2);
   - Q_u = n^{3/2}·|2ḟ/φ̇₀ − 2α − 2fφ̈₀/φ̇₀²|/(π√2);
   - Q_t = n^{3/2}·√(n(n+2))·|f|/(a·φ̇₀·π√2).
   - ρ₃ se multiplica por n^{3/2}; el perfil de Y ya lleva su rms.
4. **Los instantes de evaluación:** 120 en la contracción (uniformes en ln a), 120 en la ventana (uniformes en t) y 160 en la expansión (uniformes en ln a).
5. **El paso de las derivadas** (ḟ y φ̈₀) es h = 10⁻⁴·mín(a/k, 1/|H|, 46) en unidades de Planck (46 es la duración de la ventana). No es relativo a |t|, por la lección de U62d.
