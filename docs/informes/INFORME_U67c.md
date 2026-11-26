# INFORME U67c · Una red que nace muy diluida olvida su densidad: OLVIDA

*Escrito el 25 de noviembre de 2026, 15:47 (Madrid). Preregistro: `docs/preregistros/ESPEC_U67c.md` (13:27), con el análisis `network/u67c_analyze.py` subido antes de que terminara ninguna red.*

## En una frase

Una red RP³ que empieza diez veces más diluida que la de producción cierra su déficit de densidad con una escala τ_ℓ = 2.46 ± 0.34 unidades de ℓ. Con las dos extrapolaciones declaradas, una red que nazca en ℓ ≈ 3.4 con ξ = 0.05, como da la estimación lineal, llega al colapso (ℓ ≈ 80) con un déficit de como mucho un 2 %. **Veredicto: OLVIDA.**

## Qué se hizo

- **Las redes:** ocho redes RP³ de 256³ con el protocolo RD_SPARSE de U67 sin cambios (arranque en ℓ_i = 1, drenaje hasta ℓ = 3.8, final en ℓ_end = ln 128 = 4.85), semillas 601–604:
  - fila A: ξ_i = 0.05;
  - fila B: ξ_i = 0.5, la de producción.
- **El coste:** cada red tardó 16–18 min con 4 hilos en la nube; en total, 2 h 20 min.
- **La medida:** el cociente R(ℓ) = ⟨ξ_A⟩/⟨ξ_B⟩, el déficit D = 1 − R, y dos ajustes sobre 2.5 ≤ ℓ ≤ 4.85 con errores de jackknife sobre los pares de semillas:
  - (exp) ln D lineal en ℓ;
  - (pot) ln D lineal en ln ℓ.

## Las puertas

| Puerta | Resultado | |
|---|---|---|
| G1: campos finitos, ξ_B(ℓ_end) ∈ [0.5, 1.2] | ξ_B(4.85) = 0.92 | pasa |
| G2: ξ_B(4.6) a menos del 15 % de U67 a 192³ (0.833) | 0.835 (+0.2 %) | pasa |

G2 dice que el protocolo de referencia no cambia de 192³ a 256³.

## Los resultados

| ℓ | ξ_A (0.05) | ξ_B (0.5) | R |
|---|---|---|---|
| 2.56 | 0.155 | 0.859 | 0.18 |
| 3.07 | 0.219 | 0.533 | 0.41 |
| 3.45 | 0.279 | 0.598 | 0.47 |
| 3.96 | 0.386 | 0.697 | 0.56 |
| 4.47 | 0.519 | 0.798 | 0.65 |
| 4.85 | 0.638 | 0.920 | 0.69 |

**El ajuste.** La pendiente de (exp) es −0.406 ± 0.058, negativa a 7σ.
- τ_ℓ = 2.46 ± 0.34;
- k = 1.47 ± 0.21.

**La extrapolación declarada.** La red física nace en ℓ_b = 3.4, donde ξ_B = 0.588, con un déficit D_b = 0.915. En el lado pesimista (τ_ℓ + σ, k − σ):
- D₈₀^exp = 1.2 × 10⁻¹²;
- D₈₀^pot = 0.017.

Los dos quedan por debajo de 0.10, así que el veredicto es OLVIDA. Con los valores centrales, D₈₀^pot = 0.009.

**Lo esperado,** τ_ℓ ≈ 4–5 y k ≈ 0.9, era una estimación a mano sobre la exploración de 192³. Con el mismo análisis, esa exploración da τ_ℓ = 3.3 y k = 1.07: 256³ cierra algo más deprisa.

## Lecturas informativas

- **Sin la última instantánea:** τ_ℓ = 2.45, k = 1.45. El final de la caja no decide.
- **Empezando en ℓ ≥ 3.0**, para dejar fuera el transitorio del drenaje de la fila B (que cae de ξ ≈ 1.05 a 0.53 entre ℓ = 2.4 y 3.1): τ_ℓ = 2.56, k = 1.51, y D₈₀^pot = 0.008. El veredicto no depende de ese transitorio.
- **La fila A crece de forma monótona** desde ξ ≈ 0.14 en ℓ ≈ 2. El arranque con ξ_i = 0.05 no queda en 0.05: el arranque disperso produce cuerdas en el transitorio inicial. Lo que se mide es el acercamiento desde un ξ unas cuatro veces menor que el de producción a la misma ℓ.

## Qué significa y qué no

**Significa** que el supuesto «la red llega al colapso con la densidad del atractor» se apoya ahora en un ritmo de cierre medido. Con él, nacer diluida, como dice la estimación lineal (`NOTA_U67C_LINEAL.md`), cambia la densidad en el colapso como mucho un 2 % en la forma más lenta. Es mucho menos que la incertidumbre de la banda, que llega a un factor de varias veces en X_univ.

**No significa:**
- **No es una medida hasta ℓ ≈ 80.** El paso de ℓ ≤ 4.85 a ℓ ≈ 80 es una extrapolación con dos formas declaradas, como la del colapso.
- **No mide el nacimiento físico.** La red no arranca en ℓ ≈ 3.4 con el espectro lineal, porque eso no cabe en una caja de 256³ (lo avisó la servilleta 3). Se supone que el ritmo de cierre no depende de cómo se llegó al déficit.
- **No cubre la era de materia.** El ordenamiento en la era de materia sigue siendo el de U67 y U67b.

## Para el preprint (lo decide el autor)

En «Network before radiation domination» bastaría añadir una frase como:

> A network ten times sparser than the production one closes its deficit with a scale of $2.5\pm0.3$ units of $\ell$ at $256^3$, which extrapolated leaves at most $2\%$ at the collapse (preregistered).

Si se adopta, U67c pasa a exportarse a `whalefall-chain`, y `paper_numbers_v4.py` comprobaría los números.

## Ficheros

- Las redes: `network/results/u67c/runs/RP3_N256_{A,B}_s60{1..4}.json`.
- El registro: `produce.out`.
- El resultado: `network/results/u67c/U67C_RESULT.json`.
- El código: `network/u67c_produce.py` y `network/u67c_analyze.py`.
