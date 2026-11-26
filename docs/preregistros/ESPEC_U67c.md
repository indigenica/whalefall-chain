# ESPEC_U67c · ¿Olvida su densidad una red que nace muy diluida? El acercamiento al atractor desde abajo, a 256³

*Escrito el 2026-11-25 a las 13:27 (Madrid), antes de lanzar ninguna de las redes de esta prueba.*

**Lo que ya se ha visto, y se declara:**
- `NOTA_U67C_LINEAL.md`: la red física nacería diluida, con ξ ≈ 0.05 a ℓ ≈ 3.4.
- `network/explora_disperso.py` (192³, semillas 101–102, exploratorio): con ξ_i = 0.05 la red va a ξ ≈ 0.51 en ℓ = 4.56, frente a 0.83 de RD_SPARSE de U67. El déficit 1 − R baja de ≈ 0.49 en ℓ = 3.4 a ≈ 0.39 en ℓ = 4.56.
- La «servilleta 3» que el autor recibió de otra sesión (24 de noviembre, 256³, semillas 201–203, revisada en `docs/informes/NOTA_SERVILLETAS_20261125.md`). Sus redes de radiación empiezan en ℓ ≈ −2.8, no en ℓ = 1, porque toma ℓ = ln τ en lugar de ℓ = 2 ln τ − ln τ_end. No sirven como arranque disperso, pero sus números se han visto.

Esta prueba usa el protocolo RD_SPARSE de U67 sin cambios y semillas nuevas.

## 1. La pregunta

El preprint supone que la red llega al colapso, en ℓ ≈ 80–90, con la densidad del atractor («Network before radiation domination» y la memoria de U67b). La estimación lineal dice que nace unas diez veces más diluida. ¿Se cierra ese déficit, y a qué ritmo?

## 2. Las redes

Son ocho redes RP³ de 256³, con `lattice.run(1, 256, 128.0, semilla, sfac=1.0, dtau=0.25, nsnap=30, damp_ell=3.8)`. Es el protocolo RD_SPARSE de U67, con el arranque en ℓ_i = 1 y el drenaje hasta ℓ = 3.8, y llega a ℓ_end = ln 128 = 4.85.

| Fila | ξ_i | Semillas |
|---|---|---|
| A | 0.05 | 601, 602, 603, 604 |
| B | 0.5 | 601, 602, 603, 604 |

- El productor es `network/u67c_produce.py`, que salta las redes ya terminadas. El análisis es `network/u67c_analyze.py`.
- Cada red tarda ≈ 10 min con 4 hilos (las de 192³ tardaron ≈ 3 min).

## 3. La medida

- **El cociente.** En cada instantánea, R(ℓ) = ⟨ξ_A⟩/⟨ξ_B⟩, con la media sobre las cuatro semillas. El déficit es D(ℓ) = 1 − R(ℓ). Las dos filas comparten la malla de instantáneas.
- **Dos ajustes** sobre las instantáneas con 2.5 ≤ ℓ ≤ ℓ_end y D > 0:
  - **(exp)** ln D = a − ℓ/τ_ℓ;
  - **(pot)** ln D = b − k ln ℓ.
  - Los errores salen de un jackknife sobre los cuatro pares de semillas.
- **La extrapolación (un modelo declarado, no una medida).**
  - La red física nace en ℓ_b = 3.4 con ξ = 0.05. Su déficit inicial es D_b = 1 − 0.05/ξ_B(3.4), con ξ_B(3.4) interpolado en la fila B.
  - Se supone que el déficit evoluciona con el ritmo medido, sea cual sea el arranque.
  - En el lado pesimista (τ_ℓ + σ, k − σ):
    - D₈₀^exp = D_b·exp(−(80 − 3.4)/τ_ℓ);
    - D₈₀^pot = D_b·(3.4/80)^k.

## 4. Las puertas

- **G1:** las ocho redes terminan con campos finitos, y ξ_B(ℓ_end) cae en [0.5, 1.2].
- **G2:** ξ_B a ℓ ≈ 4.56 (la instantánea más próxima) queda a menos de un 15 % de la media de RD_SPARSE de U67 a 192³ en ℓ = 4.56 (0.833). Comprueba que la caja mayor no cambia el protocolo de referencia.
- Si falla una puerta, el veredicto es **INCONCLUSO**.

## 5. El veredicto

- **OLVIDA** si el déficit disminuye de forma significativa (la pendiente de (exp) es negativa a 2σ) y los dos D₈₀ pesimistas son menores que 0.10.
- **NO OLVIDA** si la pendiente de (exp) no es negativa a 1σ, o si alguno de los dos D₈₀ pesimistas es al menos 0.30.
- **INTERMEDIA** en los demás casos.

**Lo esperado,** por la exploración de 192³: τ_ℓ ≈ 4–5 y k ≈ 0.9. Eso da D₈₀^exp ≪ 10⁻³ y D₈₀^pot ≈ 0.05, así que OLVIDA, aunque la forma potencial queda cerca del umbral si k − σ baja de ≈ 0.6.

**Lecturas informativas:**
- R en cada instantánea;
- los dos ajustes sin el lado pesimista;
- las dos formas sin la instantánea final;
- la comparación con las redes de 192³.

## 6. Salidas

- `network/results/u67c/`: las redes en `runs/`, `produce.out` y `U67C_RESULT.json`.
- El informe va en `docs/informes/INFORME_U67c.md`.

---

## Nota del 2026-11-25 a las 17:55 (Madrid): el material interno

**Por qué esta nota.** Este preregistro se exporta a `whalefall-chain` porque la v4.1 cita su resultado. Las «servilletas» del 24 de noviembre y su revisión son material interno del proyecto, por decisión del autor. Por eso no se exportan:
- `audit/servilletas_20261124/`;
- `docs/informes/NOTA_SERVILLETAS_20261125.md`.

Las referencias a ellas de arriba no tienen destino en la exportación.

**Lo que se usó de ese material:** solo lo que aquí se declara como ya visto. No se reutilizó su código. Esta prueba, su código y sus resultados están completos en la exportación y se reproducen sin él.

**Lo que no cambia:** las filas, los umbrales y el veredicto.
