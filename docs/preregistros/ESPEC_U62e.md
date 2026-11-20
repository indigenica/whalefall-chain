# U62e · Los hijos visibles del cuscutón frente a Planck: el C_ℓ completo con CLASS · Preregistro

*2026-11-16 22:24 (Madrid). Claude, en la nube, en la rama `cuscuton-explorer`, a petición del autor («adelante con el C_ℓ completo con CLASS»). Escrito antes del código de U62e. Interno: no se exporta.*

## 1. La pregunta

Con el C_ℓ completo, ¿son compatibles con el TT de baja ℓ de Planck 2018 los hijos visibles del cuscutón?
- U62d proyectó sus espectros primordiales solo con Sachs–Wolfe, y contra ΛCDM como dato.
- El C_ℓ completo añade el ISW temprano y tardío y el Doppler. El dato real es el de Commander, con su cuadrupolo bajo.

## 2. Comprobaciones previas (declaradas)

En `results/u62e_runs/previo/viabilidad_class.{py,out}`, sin ningún espectro con rasgos:
- **CLASS en plano** con l_max = 40 y todos los ℓ tarda 0.2 s.
- **Con q_linstep = 0.02** (frente a 0.45 por defecto), los C_ℓ cambian un 1.8 %: el muestreo en k por defecto no resuelve los ℓ bajos.
- **CLASS cerrado** con Ω_K = −2.7×10⁻⁵ no terminó en más de 20 min, porque las hiperesféricas con |Ω_K| tan pequeño no son viables aquí.

**Se sabe de U62d, y se declara:** el muestreo de F con el estado A1 en los dos hijos de A·H = 1 000 no pasó su prueba (G6 de U62d).

## 3. El montaje

- **Los espectros primordiales:** los F_tot(n) de U62d (`results/u62d_cl_rebote.json`), sin calcular modos nuevos, para los cuatro hijos de la hoja χ > 0 y los estados A1 y A2.
  - F se interpola linealmente en (ln n, ln F) y vale 1 fuera de la rejilla.
  - La correspondencia con k es k = √(n(n+2))·√|Ω_K,0|·H₀/c [Mpc⁻¹], con Ω_K,0 el de cada hijo.
  - P(k) = A_s·(k/0.05)^{n_s−1}·F(k).
- **CLASS** (classy):
  - **El fondo** es el de Planck 2018: h = 0.6736, ω_b = 0.02237, ω_cdm = 0.1200, τ = 0.0544, A_s = 2.100×10⁻⁹, n_s = 0.9649, N_ur = 3.046, sin neutrinos masivos.
  - **La geometría es plana.** La curvatura del hijo (|Ω_K| ≤ 5×10⁻⁵) entra solo a través de la correspondencia n → k. Su efecto sobre las funciones de transferencia es O(Ω_K) y no se comprueba (§2).
  - **Los ajustes:** salida TT sin lente, l_max = 40, todos los ℓ (l_logstep = 1.001, l_linstep = 1) y q_linstep = 0.02.
  - **El espectro** entra por `P_k_ini type = external_Pk`, con una tabla de k densa: cada n entero de la rejilla, más una malla logarítmica de 10⁻⁸ a 10 Mpc⁻¹.
- **La referencia:** el mismo CLASS con F ≡ 1 (tabla externa).
- **La verosimilitud:** en ℓ = 2–29, la aproximación de Wishart en cielo completo con f_sky = 0.8, −2 ln L = f_sky·Σ(2ℓ+1)[Ĉ_ℓ/C_ℓ + ln C_ℓ], con Ĉ_ℓ los D_ℓ de Commander.
  - **Δχ²_Planck** = −2 ln L(F) + 2 ln L(F ≡ 1), sobre el mismo fondo. Es negativo si los datos prefieren el espectro del hijo.
  - **Se informa también Δχ²_ΛCDM,** con Ĉ_ℓ = C_ℓ de la referencia, que es la comparación de U62d hecha con el C_ℓ completo.
- **Pares definidos.** Un par (hijo, estado) está definido si, con F interpolado sobre la mitad de la rejilla (los puntos pares), Δχ²_Planck cambia ≤ máx(2, 20 %).

## 4. Puertas (si falla alguna, U62e es INCONCLUSO)

- **G0, la tabla externa:** con F ≡ 1, el C_ℓ de la tabla externa coincide con el de la ley de potencias analítica de CLASS a ≤ 0.5 % en ℓ = 2–29.
- **G1, la precisión:** con q_linstep = 0.01, Δχ²_Planck cambia ≤ máx(0.5, 5 %) en todos los pares.
- **G2, la cobertura:** cada hijo tiene al menos un estado definido.

## 5. Veredicto

Sobre los pares definidos, con Δχ²_Planck. Un hijo queda:
- **excluido** si Δχ² ≥ 25 en todos sus estados definidos;
- **compatible** si Δχ² < 9 en todos;
- **preferido** si Δχ² ≤ −4 en todos.

| veredicto | condición |
|---|---|
| INCONCLUSO | falla una puerta |
| VISIBLES_EXCLUIDOS | los cuatro hijos, excluidos |
| VISIBLES_COMPATIBLES | los cuatro, compatibles |
| VISIBLES_PREFERIDOS | los cuatro, preferidos |
| MIXTO | cualquier otro caso |

**Se informan siempre:**
- C_ℓ/C_ℓ^ref en ℓ = 2–29;
- Δχ²_Planck y Δχ²_ΛCDM por par;
- el contraste con el R_ℓ de Sachs–Wolfe de U62d;
- el Δχ² de los ℓ = 2–5 por separado.

## 6. Lo que no cubre

- **La verosimilitud oficial de Commander** (no hay clik). La aproximación de Wishart con f_sky = 0.8 es algo optimista.
- **La polarización** (EE, TE) y la lente.
- **Un ajuste conjunto** de los parámetros: el fondo se fija en Planck 2018.
- **La geometría cerrada en CLASS.**
- **El hijo de la hoja χ < 0.**

## 7. Lo que se espera

- **G0–G2 pasan.**
- **Con A·H = 30 y 100, compatibles,** quizá con Δχ²_Planck algo positivo por el cuadrupolo alto frente a un dato bajo.
- **Con A·H = 1 000, excluidos** con A2 (Δχ² ≫ 25), porque el dato tiene el cuadrupolo bajo.
- **Veredicto: MIXTO.** Confianza media.

## 8. Salidas

- **Código:** `chain/u62e_cl_class.py`.
- **Resultados:** `results/u62e_cl_class.json` y `results/u62e_runs/`.
- **Informe:** sección U62e de `docs/informes/INFORME_U62.md`.
- Interno.
