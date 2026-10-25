# U19 · c/η en la red: la caja de la nube no llega a ver terminar el colapso

*14 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U19.md`, con las enmiendas 1–3 (todas fechadas).
- *Código:* `network/kern.c` (`set_tilt`, desactivado por defecto), `network/u19_produce.py` y `network/u19_analyze.py`.
- *Resultados:* `network/results/u19/` (validación, 4 redes U(1), 1 RP³ y `U19_RESULT.json`).

## 0. En una frase

**Veredicto preregistrado: INDETERMINADO (contador no válido).** Las herramientas funcionan:
- la ruptura explícita da la masa correcta a 0.5 %;
- el contador de partículas reproduce la ecuación homogénea a ≤ 1.6 %.

Pero en una caja de 320³ el final llega cuando m_G/H ≈ 10. Para entonces el desalineamiento sigue oscilando con amplitud de ~1 rad (12–31 % de los sitios lejos del vacío) y las cuerdas no han desaparecido (ξ = 0.2–0.7). **La prueba no puede medir c/η con esta dinámica.** Para RP³ se produjo una sola red de confirmación (enmienda 3).

## 1. Validación (PASA)

| | RP³ | U(1) |
|---|---|---|
| V1: ω²/m_G² − 1 (Minkowski) | −0.49 % | −0.49 % |
| V2: R/R_EDO − 1, arranque pronto | +0.13 % | −0.0006 % |
| V2: R/R_EDO − 1, arranque denso (ℓ_i = 3.5) | +1.6 % | +0.6 % |
| Anarmonicidad R(2.5 rad)/R(0.3 rad), pronto / denso | 3.84 / 2.62 | 3.87 / 2.72 |

Con ε = 0 la evolución es idéntica bit a bit a la anterior (`make test`).

**Lo que dice la anarmonicidad.** Un ángulo inicial grande multiplica por 2.6–3.9 el número de partículas respecto al caso armónico, según cuándo se suelte el campo. En una red los ángulos iniciales son aleatorios, así que el factor anarmónico pesa tanto como el colapso que se quería medir.

## 2. Las redes

| Red | ξ₁ℓ₁ | K_lat | ξ al final / máximo | fracción lejos del vacío al final |
|---|---|---|---|---|
| U(1) poco densa, s81 | 1.57 | 2.9 | 0.21 / 0.40 | 12 % |
| U(1) poco densa, s82 | 1.92 | 8.0 | 0.48 / 0.54 | 19 % |
| U(1) densa, s81 | 10.7 | 3.9 | 0.58 / 6.3 | 29 % |
| U(1) densa, s82 | 10.9 | 4.1 | 0.74 / 6.3 | 26 % |
| RP³ poco densa, s81 | 2.5 | 3.1 | 0.52 / 0.74 | 31 % |

- **El contador lineal no es válido al final** (la enmienda 1 exige < 1 % de sitios lejos del vacío). En RP³, además, el número «conservado» baja un 27 % en las tres últimas instantáneas, cuando debería ser constante.
- **El colapso no termina.** Las paredes bajan a ≤ 4 % de su volumen máximo, pero eso pasa porque se adelgazan frente al horizonte, no porque desaparezcan: las cuerdas siguen ahí. El criterio de ESPEC_U19 (volumen de pared ≤ 5 % del máximo) resultó ser mal indicador del final del colapso, y se dice.
- **Ajuste U(1):** B = −0.16 ± 0.30, sin señal. K_lat no depende de ξℓ, y su dispersión entre semillas (2.9 frente a 8.0) es del orden de su valor.

## 3. Qué hace falta para medir c/η

- **Más rango dinámico:** m_G/H ≳ 30–50 al final, para que el colapso termine y el campo sea casi lineal. Con m_r·a·dx = 1, eso pide N ≳ 1000 (WSL, solo con U(1)) o técnicas de rango dinámico. Hiramatsu et al. usan una masa que crece con T y un parámetro de grosor de cuerda; con el truco de Press–Ryden–Spergel, el núcleo de anchura comóvil fija daría ~e² más de rango.
- **Un contador para el régimen no lineal:** la acción adiabática de la oscilación (energía/frecuencia, promediada sobre varios periodos), o medir cuando la amplitud ya sea pequeña.
- **Separar la anarmonicidad:** comparar con la misma red sin cuerdas (la misma distribución de ángulos, imposible en una red real) o modelarla.

**Para la ballena:** U19 no cambia nada. c/η de RP³ sigue sin medir, y U17 se lee con el recorrido c/η ∈ [3, 30], con c/η ≈ 3 como valor de la literatura de U(1).
