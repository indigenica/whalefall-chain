# U62 · ¿Excita el rebote con cuscutón los modos observables? Bogoliubov a través del rebote · Preregistro

*2026-11-16 08:57 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («luego U62»). Es la versión acotada de la prueba que proponía `INFORME_U61.md` §4. Escrito antes del código de U62. Interno: no se exporta.*

## 1. La pregunta

Con el juego adoptado, ¿en qué estado salen del rebote los modos escalares que serían observables en los hijos visibles de U61b?
- Si salen en el vacío, el espectro primordial es el de siempre.
- Si salen excitados (número de ocupación N_k de orden 1), su potencia se multiplica por un factor entre (|α| − |β|)² y (|α| + |β|)², que oscila con n. Esos serían rasgos de orden uno en P_ζ, y el espectro observado no los tiene.

## 2. Lo que ya se sabe, y se declara

- **U59 y U59b:** la ganancia de energía en φ ∈ [−0.2, 0.2] llega a ≈ 14 y no decae con n hasta 3×10⁵. Eso sugiere N_k ~ (E − 1)/2 de orden 1 en todos los n.
- **El explorador (C3c):** justo en H = 0, las frecuencias congeladas de los dos modos son ≈ k, y no k y k/√3. Hay un cambio brusco, independiente de k, que encaja con una mezcla que no decae con n.
- **U61b:** los hijos visibles tienen n_mem ≈ 400–2 200 y Ω_K,0 entre −5×10⁻⁷ y −5×10⁻⁵. Los modos de ℓ = 2–30 tienen n_obs = ℓ/sen(D√|Ω_K,0|), unos 10²–10⁴.

## 3. El montaje

- **Los fondos:** los cinco hijos visibles de U61b, la ventana φ ∈ [−0.2, 0.2] de cada uno (U59.Trajectory, con el estado de su rebote de `explora_cuscuton.run`) y las matrices de U59, en el gauge de radiación uniforme, q = (ζ, x).
- **Los modos:** n ∈ {100, 200, 500, 10³, 2×10³, 5×10³, 10⁴, 3×10⁴, 10⁵} por hijo, más los n_obs de ℓ = 2, 10 y 30 de cada hijo.
- **Los estados:**
  - en φ = −0.2, el vacío instantáneo de la energía E (definida positiva en los extremos, G4 de U59). Los modos normales se obtienen por la diagonalización simpléctica de J·E y se normalizan con Ω(f, f*) = i;
  - se evolucionan con la matriz de transferencia;
  - en φ = +0.2 se descomponen en los modos normales instantáneos: f_j = Σ_k (α_jk g_k + β_jk g_k*).
- **Las cantidades:**
  - N_k = Σ_j |β_jk|²;
  - y para el modo de χ, el que ocupa el slot del inflatón (el de mayor frecuencia, de velocidad ≈ 1), los factores de potencia R_min = Σ_j(|α_jk| − |β_jk|)² y R_max = Σ_j(|α_jk| + |β_jk|)².
- **El supuesto** (declarado, no comprobado aquí): que después de φ = +0.2 la evolución es adiabática hasta que el modo sale del horizonte en la inflación, de modo que N_k se conserva.

## 4. Puertas (si falla alguna, U62 es INCONCLUSO)

- **G0, la normalización:** en φ = −0.2, los modos construidos cumplen Ω(f_j, f_k*) = iδ_jk y Ω(f_j, f_k) = 0 a ≤ 10⁻¹⁰.
- **G1, Bogoliubov:** Σ_k(|α_jk|² − |β_jk|²) = 1 para cada j, a ≤ 10⁻⁶ (la simplecticidad de la transferencia).
- **G2, la convergencia:** con rtol de 10⁻⁹ y de 10⁻¹⁰, N_k cambia ≤ 10⁻⁴ absoluto + 10⁻⁴ relativo, en tres n.

## 5. Veredicto

Sobre el modo de χ, con los n_obs de ℓ = 2–30 de los cinco hijos visibles:
- **INCONCLUSO:** falla una puerta.
- **REBOTE_ADIABATICO:** N_χ ≤ 0.01 en todos (rasgos de potencia ≤ 20 %).
- **REBOTE_EXCITA:** N_χ ≥ 0.1 en al menos 3 de los 5 hijos (rasgos de potencia de al menos ±60 %).
- **INTERMEDIO:** cualquier otro caso.

**Se informan siempre:**
- N_k frente a n (los dos modos);
- R_min y R_max;
- la relación con la ganancia de energía de U59b.

## 6. Lo que no cubre

- La evolución tras φ = +0.2 (el supuesto de adiabaticidad).
- La proyección a C_ℓ.
- Los modos con k_fís ≫ M_Pl en el rebote, los de los hijos con muchos e-folds, fuera de la EFT.
- Otros estados iniciales: vacíos de orden superior, o el estado puesto en el giro de la madre.

## 7. Lo que se espera

- **G0–G2 pasan.**
- **REBOTE_EXCITA,** con confianza media: N ~ 1–6, por la ganancia de energía de U59b.
- **Si sale así, la consecuencia** es que con este potencial los hijos visibles del cuscutón tendrían rasgos de orden uno en P_ζ a ℓ ≲ 30, que los datos excluyen. Faltaría ver si el supuesto de §3 se sostiene.

## 8. Salidas

- **Código:** `chain/u62_bogoliubov_rebote.py`.
- **Resultados:** `results/u62_bogoliubov_rebote.json`.
- **Informe:** `docs/informes/INFORME_U62.md`.
- Interno.

---

## Enmienda 1 · 2026-11-16 09:00 (Madrid), antes de escribir el código

Aclara tres puntos que el texto deja abiertos. No se ha calculado nada de U62.

1. **Qué hijo cuenta en el veredicto:**
   - un hijo cuenta como **excitado** si N_χ ≥ 0.1 en **alguno** de sus tres n_obs (ℓ = 2, 10 y 30);
   - para **REBOTE_ADIABATICO**, N_χ ≤ 0.01 en los quince valores (5 hijos × 3 ℓ).
   - Se informa también la lectura estricta (N_χ ≥ 0.1 en los tres ℓ).
2. **N_χ** es el número de cuantos que acaban en el modo de salida de mayor frecuencia: N_χ = Σ_j |β_jχ|², sumando sobre los dos modos de entrada.
3. **Las coordenadas.** Para el vacío instantáneo se usan las variables canónicas reescaladas (a^{3/2}q, a^{−3/2}p), una transformación simpléctica en cada instante en la que la matriz del hamiltoniano congelado es exactamente la E de U59. Los n_obs se redondean al entero más próximo. Los hijos se toman de `results/u61b_cuscuton_huella.json` con la misma condición de visibilidad de U61b.
