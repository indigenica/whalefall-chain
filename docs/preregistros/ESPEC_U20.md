# U20 · Isocurvatura de ruido blanco: estimación propia sobre la banda de U21 · Preregistro

*16 de octubre de 2026 (Madrid; la hora que cuenta es la de su commit). Escrito antes de ejecutar U20. Se ejecuta en la nube, en segundos.*

## Motivo

U21 usó como puerta la cota de isocurvatura de la literatura, m_G ≥ 10⁻¹⁷ eV (arXiv:2511.15790, del que solo se leyó el resumen). Antes de apoyarse en esa cota hace falta una estimación propia, aunque sea de orden de magnitud, con nuestro T_osc, nuestra época del colapso (r) y nuestro espectro adiabático.

## Modelo (declarado)

**El campo.** La materia oscura es postinflacionaria: tanto el desalineamiento como la población del colapso están sin correlacionar más allá del horizonte en el colapso. Su contraste de densidad es ruido blanco por debajo de k_* = a_d·H_d = a_osc·H_osc/r (comóvil):
- Δ²_iso(k) = C·(k/k_*)³ para k ≪ k_*;
- **C** es la normalización del ruido blanco, O(1) e incierta: **primaria C = 1**, y se informan C ∈ {0.3, 3, 10}.

**El crecimiento.** El modo de isocurvatura de CDM, dentro del horizonte, crece como en Mészáros. En z = 3, D_iso = 1 + 1.5·g_b·(1 + z_eq)/(1 + z), con z_eq = 3400 y g_b = 0.9 por el retraso de los bariones. Es una aproximación: su incertidumbre, un factor ~1.2 en Δ², queda dentro del rango de C.

**La libre propagación** de la población del colapso, que suprime la isocurvatura a k > k_fs, **no se incluye.** Es el lado conservador.

**La comparación.** Δ²_ad(k, z = 3) sale del espectro lineal de CLASS (caché de U4a, ΛCDM).

## Criterio

En k = 10 h/Mpc, z = 3, una escala del bosque Lyα: **Δ²_iso/Δ²_ad ≤ 0.10.** Cada fila de U21 de las variantes P (r = 1) y H (r = 3) pasa o falla con C = 1.

## Veredicto (variante P, C = 1)

- **ISOCURVATURA_CONFIRMA_BANDA:** pasan todos los F de la banda de U21 (≤ 0.7·F₀).
- **ISOCURVATURA_RECORTA_BANDA:** el mayor F que pasa es menor que 0.7·F₀, pero pasa alguno. Se da el nuevo F_max.
- **ISOCURVATURA_CIERRA_BANDA:** no pasa ningún F de la banda.

**Se informa también:** la masa mínima m_min(C) que satisface el criterio, comparada con 10⁻¹⁷ eV (la literatura), y la variante H.

## Lo que se espera, escrito antes de calcular

Con C = 1, la cuenta a mano da Δ²_iso/Δ²_ad ≈ 1.8 en F₀ (falla), ≈ 0.2 en 0.7·F₀ (falla) y ≈ 0.03 en 0.5·F₀ (pasa). **Se espera ISOCURVATURA_RECORTA_BANDA, con F_max ≈ 0.5·F₀,** y m_min ≈ 2–4×10⁻¹⁷ eV, cerca de la cota de la literatura.

## Lo que U20 no resuelve

- **La normalización C del ruido blanco** viene de la física del colapso: U19b podría medirla con el espectro de densidad de la red.
- **No se calcula el espectro de flujo de Lyα** ni se hace una verosimilitud: el criterio del 10 % en la potencia lineal es una aproximación.
- **No hay libre propagación,** lo que favorece al lado conservador.
