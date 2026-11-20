# U62c · El espectro primordial de los hijos visibles del cuscutón, a través del rebote · Preregistro

*2026-11-16 17:35 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con U62c»). Escrito antes del código de U62c. Interno: no se exporta.*

## 1. La pregunta

¿Tiene P_ζ rasgos en ℓ = 2–40 en los cinco hijos visibles de U61b?
- U62b mostró que la memoria del rebote se borra por encima de unas pocas n_mem.
- Pero los modos con n ≲ n_mem nunca son adiabáticos cerca del rebote, así que el número de partículas no responde por ellos.
- Aquí se calcula P_ζ directamente: los modos se propagan desde la madre que se contrae hasta que se congelan fuera del horizonte, en la inflación.

## 2. Comprobación previa (declarada)

Hecha antes de este ESPEC: `results/u62c_runs/previo/viabilidad.{py,out}`. Solo fondo y matrices, sin propagar ningún modo:
- en el tramo de expansión, la radiación cae a 10⁻⁴ de la densidad total a Δln a ≈ 5.1 de la ventana, y a 10⁻⁵ a ≈ 5.7, ya en inflación;
- el número de condición de K del gauge de radiación uniforme crece como 1/(ρ_r/ρ): ≈ 10³ y ≈ 10⁴;
- el residuo de Friedmann se queda en ≤ 1.1×10⁻⁵ sin proyección;
- la contracción, desde φ = −0.2 hasta el giro de la madre, dura 4.5–6.4 e-folds.

## 3. El montaje

**Etapa 1 (dos campos, como U62b).**
- **Las matrices:** las de U59 en el gauge de radiación uniforme, con el fondo de `explora_cuscuton`:
  - un tramo de contracción en el modo −1, integrado hacia atrás desde φ = −0.2;
  - la ventana de U59.Trajectory;
  - un tramo de expansión en el modo +1, desde φ = +0.2 hasta que ρ_r/ρ = 10⁻⁵.
- **Tres estados iniciales.** Cada uno es el vacío instantáneo del marco sin B_s (el de U62b), con los dos modos de entrada.
  - **A1:** en la contracción, a 0.5 e-folds del giro de la madre. Allí H ≈ 0 y los modos son lo más adiabáticos posible. Si E′ no es positiva allí, se acerca al rebote en pasos de 0.1 e-folds hasta que lo sea, y se anota.
  - **A2:** a Δln a = 2 antes de la ventana, como en U62b.
  - **B:** a Δln a = 2 después de la ventana. Es un hijo que «empieza de cero» tras el rebote, sin memoria de la contracción.
- **El cambio de etapa** se hace cuando ρ_r/ρ = 10⁻⁴. Se pasa a la curvatura comóvil con χ, la del gauge C de C3d del explorador: ζ_C = ζ_R − (H/χ̇)x_R, con su derivada temporal completa. Se descarta la radiación.

**Etapa 2 (un campo).**
- **La acción** de Boruah, Kim y Geshnizjani, con z² y c_s² de `u58_cuscuton_inflacion.z2_cs2`.
- **El fondo:** U58.Model, plano, con el juego adoptado, desde (χ, dχ/dN) del hijo en el cambio de etapa. Se descartan la curvatura (|Ω_K| ≲ 4×10⁻⁴ allí) y la radiación.
- **El modo** se sigue hasta k/aH = 10⁻³, con k = √(n(n+2)) y a en las unidades del hijo.
- **El espectro:** P_ζ(n) = n³ Σ_j |ζ_j|²/(2π²), sumando los dos modos de entrada. Con el armónico normalizado en S³ (∫Y² = 1), este P coincide con la convención plana para n ≫ 1.

**La referencia (el atractor).**
- **Dónde se compara:** para cada n, el valor χ_× del inflatón cuando el modo sale del horizonte por última vez (k = aH) en el hijo.
- **Qué se compara:** el P_ζ de U58.mode, con el vacío de Bunch–Davies, en el fondo plano de U58 que arranca en χ = 6, para el k que sale en χ_×.
- **Los cocientes:** F_tot(n) = P_A(n)/P_att(n) y F_reb(n) = P_A(n)/P_B(n).
  - F_tot incluye la huella del rebote y la de la era de radiación previa a la inflación.
  - F_reb aísla, en lo posible, la contracción y el rebote.

**Los modos:** por hijo, 16 n espaciados logarítmicamente entre n(ℓ = 2) y n(ℓ = 40), con n(ℓ) = ℓ/sen(D√|Ω_K,0|) como en U62, más un control n_c = 10·n_mem.

## 4. Puertas (si falla alguna, U62c es INCONCLUSO)

- **G0, la normalización de los modos de entrada:** ≤ 10⁻¹⁰.
- **G1, el wronskiano al final de la etapa 2:** Σ_j z²(ζ_jζ_j′* − c.c.) = i/2 con un error ≤ 10⁻³ (relativo) en todas las filas. Comprueba a la vez la normalización de la acción de la etapa 1 frente a la de la etapa 2.
- **G2, la tolerancia:** rtol 10⁻⁹ frente a 10⁻¹⁰ en n(ℓ = 2), n(ℓ = 10) y n(ℓ = 40) del primer hijo. P_A1 cambia ≤ 10⁻³.
- **G3, el control:** |F_tot(n_c) − 1| ≤ 0.03 con A2 en los cinco hijos.
- **G4, el cambio de etapa:** con ρ_r/ρ = 10⁻⁵ en lugar de 10⁻⁴, P_A1 cambia ≤ 1 % en los tres n de G2.
- **G5, el fondo:** el residuo de Friedmann es ≤ 3×10⁻⁵ en los tramos.

## 5. Veredicto

D_X(hijo) = máx_n |F_tot(n) − 1| sobre los 16 n de ℓ = 2–40, para cada estado X = A1, A2.
- **INCONCLUSO:** falla una puerta.
- **ESPECTRO_LISO:** D_A1 ≤ 0.1 y D_A2 ≤ 0.1 en los cinco hijos.
- **RASGOS_FUERTES:** mín(D_A1, D_A2) ≥ 0.5 en al menos 3 hijos.
- **RASGOS_MODERADOS:** cualquier otro caso.

**Se informan siempre:**
- F_tot(n) y F_reb(n) para A1, A2 y B;
- la dependencia del estado, |P_A1/P_A2 − 1|;
- los D restringidos a ℓ ≥ 10;
- χ_× y la adiabaticidad de A1.

## 6. Lo que no cubre

- **La proyección a C_ℓ,** la función de transferencia y la varianza cósmica. Los rasgos primordiales no son, por sí solos, una exclusión.
- **El estado de la madre antes del giro,** que está abierto, y otras prescripciones.
- **La curvatura y la radiación en la etapa 2** (errores ≲ 10⁻³).
- **Los tensores.**

## 7. Lo que se espera

- **G0–G5 pasan.** La menos segura es G3: la convención de normalización y el ajuste de la referencia.
- **F_tot** se acerca a 1 en n ≫ n_mem, y tiene rasgos de decenas de por ciento o más en n ≲ n_mem.
- **Veredicto:** RASGOS_FUERTES o RASGOS_MODERADOS. Confianza baja en cuál de los dos.
- **Habrá una dependencia del estado importante** en los n más bajos.

## 8. Salidas

- **Código:** `chain/u62c_espectro_rebote.py`.
- **Resultados:** `results/u62c_espectro_rebote.json` y `results/u62c_runs/`.
- **Informe:** sección U62c de `docs/informes/INFORME_U62.md`.
- Interno.
