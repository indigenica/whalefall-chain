# Preregistros

Cada prueba se escribió **antes** de ejecutar su código. Las enmiendas van al final de cada documento, con su fecha y hora y con el motivo, siempre antes de la producción a la que afectan. Las erratas posteriores a la ejecución se marcan como tales.

| Documento | Prueba | Resultado |
|---|---|---|
| `ESPEC_U1_CIERRE.md` | U1: cierre cosmológico condicional (barrido de la eficiencia ε de la red) | `docs/informes/INFORME_U1_CIERRE.md` |
| `ESPEC_U2.md` | U2: medida de ε en redes 3D RP³ frente a U(1) | `docs/informes/INFORME_U2.md` |
| `ESPEC_U3.md` | U3a (Lyα con Boltzmann) y U3b (redes en el atractor) | `docs/informes/INFORME_U3.md` |
| `ESPEC_U4.md` | U4a (Lyα con CLASS) y U4b (redes grandes) | U4a: CIERRE_SEGÚN_RAMA; U4b fase 1: GATES_PASSED; fase 2 en marcha (`docs/informes/INFORME_U4.md`) |
| `ESPEC_U5.md` | rama B4 con la pendiente de ξ medida, malla de forma ampliada, isocurvatura, más semillas | preregistrado el 8 oct, sin ejecutar |
| `ESPEC_U12.md` | cruce de la densidad de cuerdas creciente (B4a, B4b, B1) con q < 1, con la R_Γ de U11 | FRIA_SOBREVIVE_AL_CRECIMIENTO (`docs/informes/INFORME_U12.md`) |
| `ESPEC_U15b.md` | la cola tibia frente a Planck + BAO + Lyman-ABGD (MCMC, WSL) | **DETENIDO** (filas físicas superadas por U17–U21; la validación PASA: f₉₅ = 0.27 y 0.42 con 1 y 2 keV), `docs/informes/INFORME_U15b_U15c.md` |
| `ESPEC_U16.md` | resolución y volumen de la red: R_ξ y R_Γ hacia el continuo | RESOLUCION_ROBUSTA; volumen: R_Γ robusto, R_ξ indeterminado (`docs/informes/INFORME_U16.md`) |
| `ESPEC_U17.md` | K_sw frente a la energía de la red en H ≈ m_G | K_SW_FUERA_DEL_RANGO_TIBIA (`docs/informes/INFORME_U17.md`) |
| `ESPEC_U17b.md` | la población del colapso por CLASS (puertas de U12) | SOLO_XI_CONSTANTE_VIVE (`docs/informes/INFORME_U17b.md`) |
| `ESPEC_U18.md` | estimadores local e integral de R_Γ con instantáneas densas | INDETERMINADO; coinciden con muestreo denso (`docs/informes/INFORME_U18.md`) |
| `ESPEC_U19.md` | c/η del colapso medido en la red | INDETERMINADO (contador no válido) (`docs/informes/INFORME_U19.md`) |

| `ESPEC_U17c.md` | barrido en η del colapso para las ramas con ξ creciente (¿vive la rama cómoda?) | **RAMA_COMODA_FORZADA** (η_crit = 0.1; Lyα depende de ξℓη, no de c), `docs/informes/INFORME_U17c.md` |
| `ESPEC_U15c.md` | la rama viva (ξ constante con el colapso) frente a Planck + BAO + Lyman-ABGD (casos en los bordes F = 0.5 y 0.35·F₀, enmienda 3) | **RAMA_VIVA_PERMITIDA_LYA_COMPLETO** (p(s = 1)/p_max = 0.87 y 0.91; posterior plana, indistinguible de CDM), `docs/informes/INFORME_U15b_U15c.md` |
| `ESPEC_U15d.md` | el borde superior de U17f (F₀ con Y = 1.20 y 1.70; 0.85 y 0.7·F₀ con Y = 1.98) frente a la verosimilitud completa de Lyα, con la población medida en U19b | **F0_EXCLUIDO_LYA_COMPLETO; BORDE_0.7_SOSTENIDO**: p(s = 1)/p_max = 0.111, 0.047, 0.114 y 0.238; `docs/informes/INFORME_U15d.md` |
| `ESPEC_U17d.md` | filas decisivas del colapso con el balance corregido y el puente RP³ (auditoría de Codex) | **RAMA_VIVA_CON_COLAPSO_TIPO_U1** (solo por forma, sin margen), `docs/informes/INFORME_U17d.md` |
| `ESPEC_U19b.md` | η y época del colapso de RP³ en la red, con fase tardía de núcleo congelado (enmiendas 1–5) | **COLAPSO_BLANDO** en la frontera y frágil (Y = 1.70 ± 0.55), `docs/informes/INFORME_U19b.md`; con ℓ_G = 2.4, lectura (c) mixta (§7 y su corrección); control U(1): Y = 1.58 ± 0.50, RP³/U(1) ≈ 1.07 (§8) |
| `ESPEC_U19c.md` | X₀(ℓ₁) del colapso de RP³ en cuatro masas con el protocolo de ℓ_G = 3.2 relativo a ℓ₁ (enmienda 1: 576³; enmienda 2: la masa alta en una A100 con JAX) | **INDETERMINADO** (σ_p = 0.76 > 0.5; X₀ = 6.1, 6.0, 8.2 y 7.7 no siguen una potencia; 640³ y 768³ no validan), `docs/informes/INFORME_U19c.md` |
| `ESPEC_U19d.md` | las semillas 93–96 en las cuatro masas de U19c (32 redes: 256³ y 384³ en la nube, 512³ y 576³ en la A100), con el análisis y el veredicto de U19c | **X_UNIV_INDECISO (entre 244 y 1 367)**: p = 1.35 ± 0.31; F₀ queda por debajo del intervalo, `docs/informes/INFORME_U19d.md` |
| `ESPEC_U17e.md` | la población del colapso medida en U19b, pasada por CLASS en la rejilla de F de U21 (filas Y, Y+ y X) | **BANDA_RECORTADA** (0.35·F₀; Lyα solo, 0.7), `docs/informes/INFORME_U17e_U20b_U23.md` |
| `ESPEC_U17f.md` | el borde superior de la banda de F en todo el intervalo del colapso de U19b (X del universo ≈ 170–320), con F = 0.85·F₀ en la rejilla (enmienda 1: meseta de suma) | **BORDE_0.7_CONFIRMADO; F₀ PERMITIDO CON X ≤ 194** (solo por forma; con la cota de supermasivos, F₀ excluido), `docs/informes/INFORME_U17f.md` |
| `ESPEC_U17g.md` | el borde superior de la banda de F sobre el intervalo de X de U19d (244–1 367): seis familias Y = 2.5–8.44 con la tubería de U17f | **BORDE_BAJA_A_0.35**: 0.6·F₀ con X ≈ 400–490, 0.5·F₀ en el X central (556), 0.35·F₀ con X ≈ 900–1 370; `docs/informes/INFORME_U17g.md` |
| `ESPEC_U23.md` | superradiancia de agujeros negros con autointeracción sobre la banda de F (WM 2025 digitalizado; corrige la nota U23) | **SUPERRADIANCIA_PARTE_LA_BANDA** (estelares excluyen 0.045–0.0625·F₀), mismo informe |
| `ESPEC_U20b.md` | isocurvatura con las cotas observacionales de f_iso (Lyα GG 2026; UFD informativa), en vez del criterio aproximado de U20 | **ISOCURVATURA_NO_LIMITA**, mismo informe |
| `ESPEC_U21.md` | banda de F permitida por Lyα, isocurvatura (literatura), ΔN_eff y cuerdas del CMB | **BANDA_BAJO_EL_BENCHMARK** (F ≤ 0.7 F₀), `docs/informes/INFORME_U21.md` |
| `ESPEC_U20.md` | isocurvatura de ruido blanco: estimación propia sobre la banda de U21 | **ISOCURVATURA_RECORTA_BANDA** (F ≤ 0.5 F₀), `docs/informes/INFORME_U20.md` |
| `ESPEC_U22.md` | la meseta como inflación de polo de orden 4 en S1–S2 | **POLO_CAMBIA_NACIMIENTO**, `docs/informes/INFORME_U22.md` |
| `ESPEC_U22b.md` | estadística de nacimiento con el polo en toda la rejilla de C45h | **NACIMIENTO_ROBUSTO** (R = 0.94), `docs/informes/INFORME_U22b.md` |
| `ESPEC_U26.md` | la inflación de polo como potencial por defecto de la cadena (S1–S6, la acción común), con la suma como opción | **POLO_ADOPTADO**, `docs/informes/INFORME_U26.md` |
| `ESPEC_U24.md` | sensibilidad de R_Γ al drenaje (ℓ_d, Γ_d), 384³, 3 semillas (enmienda 1) | **DRENAJE_ROBUSTO** (max \|δR_Γ\| = 1.9 %; R_ξ −5.6 % con el doble de drenaje), `docs/informes/INFORME_U24.md` |
| `ESPEC_U25.md` | dτ/2 con el protocolo D0 de U24 (enmienda 1: sesgo b = 2·\|Δ\|) | **PASO_ROBUSTO** (sesgo < 0.8 % en R_Γ), `docs/informes/INFORME_U25.md` |

**Rutas.** Los preregistros conservan las rutas que existían cuando se escribieron, porque no se editan. En este repositorio corresponden a:

| Ruta en el preregistro | Ruta en el repositorio |
|---|---|
| `src/` | `chain/` |
| `vendor/` | `components/` |
| `inputs/` | `data/` |
| `u2_red/` | `network/` (resultados en `network/results/u2/`) |
| `u3/net/` | `network/` (resultados en `network/results/u3/`) |
| `u3/lya/` | `lya/` |
| `audit_u1/`, `u2_red/audit_u2/` | `audit/u1/`, `audit/u2/` |
| `docs/c45b_ballena5/` | `docs/modelo/`, `docs/literatura/`, `docs/evidencia/nacimiento_c45/` |
