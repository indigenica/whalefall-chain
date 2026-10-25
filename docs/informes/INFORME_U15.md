# U15 · La cola tibia frente a Lyman-ABGD: arnés no validado

*11 de octubre de 2026, noche.*
- *Preregistro:* `docs/preregistros/ESPEC_U15.md`, con la enmienda 1.
- *Código:* `lya/u15_lya_abgd.py`. Usa la verosimilitud pública `Lya_abgd` de MontePython sin modificarla (commit `5b365b0`).
- *Resultado:* `lya/results/u15/U15_RESULT.json`.

## 0. En una frase

**Veredicto preregistrado: INCONCLUSO (arnés no validado).** Con la cosmología fija y perfilando los nueve parámetros astrofísicos, la verosimilitud de MIKE/HIRES no distingue la WDM de 2 keV de ΛCDM. Ni siquiera la penaliza: Δχ² ≈ −2. Por tanto, el arnés no reproduce las cotas publicadas, y no se ha evaluado ninguna fila del modelo.

## 1. Validación

Δχ² = χ²_min(modelo) − χ²_min(ΛCDM). Entre paréntesis, el criterio de área de la propia verosimilitud.

| Modelo | Sin enmienda | Con la enmienda 1 (rejilla térmica) | Criterio de V |
|---|---|---|---|
| WDM 2 keV | −0.57 | −2.36 (0.66) | > 4 |
| WDM 3.5 keV | −2.09 | −3.88 (0.44) | decreciente con m |
| WDM 5.3 keV | +0.37 | −1.38 (0.29) | decreciente con m |
| WDM 10 keV | −0.10 | −0.11 (0.12) | < 4 |
| mezcla 1 keV, f = 0.5 | −2.68 | −4.46 (0.64) | > 4 |
| mezcla 2 keV, f = 0.2 | — | −2.09 (0.16) | < 4 |

El χ² mínimo de ΛCDM es 35.4, con 70 puntos de MIKE/HIRES.

## 2. Lectura

- **No es un fallo de la ballena, es del arnés.** El perfil a cosmología fija deja libres F(z) en cada redshift, T0(z) y γ(z). Con esa libertad, los datos de alto z de MIKE/HIRES prefieren algo de supresión. Es una preferencia conocida en esos datos cuando la historia térmica es libre.
- **Las cotas publicadas no salen de este procedimiento.** Hooper et al. combinan Lyα con Planck + BAO en una cadena MCMC y dan cotas marginales (m > 5.9 keV para WDM pura, con un prior plano en m). Una cadena así necesita las verosimilitudes de Planck (clik), que no hay aquí.
- **La enmienda 1** (rechazar historias térmicas fuera de la rejilla, como dice el artículo) era necesaria, pero no basta.

## 3. Qué haría falta

- **La cadena completa en WSL:** MontePython con Planck 2018 + BAO + `Lya_abgd`, una cadena por fila representativa (por ejemplo, B2 con q = 1, R_Γ = 1.56, x₀ = 10, K_sw = 0), comparando la posterior con la de ΛCDM. Hay que preregistrarla aparte.
- Mientras tanto, la cola tibia sigue siendo **marginal** según la extrapolación de Garcia-Gallego et al. La fórmula de Hooper et al. para mezclas (m ≳ 7.2 keV (f − 0.1)) dejaría f ≲ 0.1 sin restringir, pero no se puede usar como veredicto: es un ajuste a sus posteriores y está fuera de su rango de masas.
