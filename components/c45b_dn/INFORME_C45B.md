# C45b · Colapso con horizonte, rebote QT y lanzamiento (informe)

*3 de octubre de 2026 · claude_dual/c45b · banco preregistrado `results/c45b_dn_20261003T181915Z` · auditoría independiente en `audit/AUDIT_C45B.md`*

## 0. En una frase

Los núcleos de χ supercríticos forman horizonte, su centro llega a la curvatura límite de la QT y **rebota**. Pero, en todo lo que el código puede seguir con fiabilidad, **el centro rebotado sigue dominado por energía cinética y, donde llegamos a verlo, el manto que sigue cayendo lo vuelve a aplastar enseguida**.

No hay lanzamiento por rebote en el dominio fiable. **Veredicto formal: INCONCLUSIVE**, porque el dominio fiable acaba demasiado pronto tras el rebote para poder descartar un lanzamiento posterior.

## 1. Pregunta y herramienta

**Pregunta (ESPEC_C59A_C45B §3).** Un colapso empalmado con un exterior asintótico, que forma horizonte: ¿rebota dentro y lanza χ a la meseta como en C45?

**Variante elegida.** El núcleo es una estrella de χ (C59a), como acordamos.

**Código nuevo** (independiente de c58c): `dn/`, en doble nulo, ds² = −2e^{2σ}du dv + r²dΩ².
- **Ecuaciones.** NPQT con escalar real y potencial ballena; mismas ecuaciones y convenios que el evaluador auditado del equipo.
- **Esquema.** Diamante implícito de orden 2.
- **Centro regular.** e^{2σ} = 2r_v² y φ par en r.
- **Masa generalizada.** M = r³h(ψ)/2G se evoluciona con su ley de flujo, y ψ se obtiene de M. Esto garantiza |αψ| < 1 sin suelos ni recortes.
- **Enmascarado causal.** Un nodo que falla pasa a NaN y solo se pierde su futuro causal.

**Dominio fiable.** Se definió de antemano: el centro es fiable mientras N = 1600 y N = 3200 coinciden (|Δχ_c| < 0.02 y |Δρ_c| < 5 %). Los eventos solo se evalúan ahí.

**Validaciones:**
- **Espacio plano:** el esquema es exacto hasta el redondeo. El test original pedía un orden y fallaba por eso; se ha enmendado y está declarado.
- **Pulso débil en GR:** orden 2.0 en las restricciones y en la autoconvergencia (auditoría).
- **Estrella antes del horizonte:** orden 1.9–2.05.
- **Ecuaciones:** rederivadas a mano en el límite GR y comprobadas contra soluciones estáticas exactas de la NPQT, contra la propagación de restricciones y contra la ley de flujo de masa. Todo **CONFIRMADO** (auditoría §1–2).

## 2. Banco preregistrado (C45b-dn)

**Familias.** Las mismas en α = 0.1 y α = 0.03:
- estrellas de χ de C59a escaladas, s = 1.75–3 (M = 0.94–1.75);
- gaussianas;
- «top-hats», el análogo más cercano al parche homogéneo de C45.

**Malla.** L = 80, N = 800/1600/3200.

**Eventos:**
- **H:** horizonte, es decir, superficie atrapada.
- **B:** rebote QT del centro. H_c comóvil pasa de − a + con αψ_c ≥ 0.9.
- **L:** lanzamiento. Durante al menos 1/H_inf de tiempo propio se cumplen a la vez fracción cinética < 0.5, ρ_c < 10 V0, H_c > 0 y |χ_c| > χ_end.
- **R:** recompresión. ρ_c supera 10 veces su mínimo tras el rebote.
- **MI:** diagnóstico de inflación de masa.

**Reglas:**
- **PASS** si hay L en la familia estrella para cada α.
- **FAIL** si no hay ningún L con horizonte y todos los rebotes muestran R.
- **INCONCLUSIVE** en otro caso.

**Declaración.** El protocolo se escribió después de sondas exploratorias (gaussiana A = 0.25, w = 3; estrellas s = 1.5–3 con α = 0.1) y de una pasada rápida que solo sirvió para fijar el estimador de ρ_c y la tolerancia de fiabilidad.

| Config. | α | M | Horizonte (u, r) | Rebote | \|χ_c\| máx. tras rebote | Fracción cinética mín. tras rebote | Recompresión | L | u fiable | M_int/M_ADM (fiable) | αψ_c máx. |
|---|---|---|---|---|---|---|---|---|---|---|---|
| estrella s=1.75 | 0.1 | 0.94 | sí (39.3, 0.60) | u=42.7 | 1.79 | 0.98 | no | no | 44.1 | 0.75 | 0.978 |
| estrella s=2.0 | 0.1 | 1.12 | sí (16.7, 0.58) | u=19.0 | 1.52 | 0.04* | **u=19.8** | no | 22.7 | 1.05 | 0.998 |
| estrella s=2.5 | 0.1 | 1.45 | sí (7.8, 1.35) | u=12.1 | 5.19 | 0.99 | no | no | 14.2 | 2.89 | 0.995 |
| estrella s=3.0 | 0.1 | 1.75 | sí (4.1, 1.80) | u=9.1 | 7.04 | 0.99 | no | no | 11.1 | 5.64 | 1.000 |
| gauss 0.25/3 | 0.1 | 1.17 | sí (15.3, 0.88) | u=18.9 | 3.16 | 0.99 | no | no | 21.0 | 1.89 | 0.992 |
| gauss 0.15/5 | 0.1 | 2.13 | no (se dispersa) | — | — | — | — | no | 79.9 | 0.07 | 0.004 |
| gauss 0.1/8 | 0.1 | 4.29 | no (se dispersa) | — | — | — | — | no | 79.9 | 0.02 | 0.003 |
| top-hat R=3 | 0.1 | 1.64 | sí (3.8, 1.35) | u=6.4 | 3.42 | 0.97 | no | no | 7.2 | 1.77 | 0.999 |
| top-hat R=4 | 0.1 | 2.53 | sí (2.6, 4.13) | u=8.1 | 3.99 | 1.00 | no | no | 8.3 | 2.33 | 1.000 |
| top-hat R=6 | 0.1 | 4.35 | sí (1.3, 8.28) | no | — | — | — | **sí (inicial)** | 9.2 | 0.99 | 0.006 |
| estrella s=1.75 | 0.03 | 0.94 | sí (39.3, 0.56) | u=42.7 | 3.78 | 1.00 | no | no | 43.7 | 0.71 | 0.997 |
| estrella s=2.0 | 0.03 | 1.12 | sí (16.7, 0.60) | u=19.2 | 2.99 | 0.00* | **u=20.1** | no | 21.1 | 1.06 | 0.999 |
| estrella s=2.5 | 0.03 | 1.45 | sí (7.9, 1.29) | no (fiable acaba antes) | — | — | — | no | 12.2 | 4.78 | 1.000 |
| estrella s=3.0 | 0.03 | 1.75 | sí (4.1, 1.83) | no (fiable acaba antes) | — | — | — | no | 8.8 | 10.0 | 1.000 |
| gauss 0.25/3 | 0.03 | 1.17 | sí (15.3, 0.89) | u=19.1 | 4.44 | 1.00 | no | no | 19.4 | 4.64 | 1.000 |
| gauss 0.15/5 | 0.03 | 2.13 | no | — | — | — | — | no | 79.9 | 0.07 | 0.001 |
| gauss 0.1/8 | 0.03 | 4.29 | no | — | — | — | — | no | 79.9 | 0.02 | 0.001 |
| top-hat R=3 | 0.03 | 1.64 | sí (3.7, 1.37) | no | — | — | — | no | 6.3 | 4.62 | 1.000 |
| top-hat R=4 | 0.03 | 2.53 | sí (2.6, 4.13) | no | — | — | — | no | 7.4 | 4.54 | 0.336 |
| top-hat R=6 | 0.03 | 4.35 | sí (1.3, 8.28) | no | — | — | — | **sí (inicial)** | 9.2 | 0.99 | 0.002 |

\* Solo en el punto de retorno de χ_c (Δτ ≈ 0.02). Fuera de él, la fracción cinética es ≥ 0.97.

**Veredicto formal: INCONCLUSIVE.** La auditoría lo confirma: FAIL queda bloqueado dos veces, por los dos L del top-hat R = 6 y por 8 rebotes sin R dentro del dominio fiable.

![centro](figures/fig_c45b_centro.png)
![espaciotiempo](figures/fig_c45b_espaciotiempo.png)
![banco](figures/fig_c45b_banco.png)

## 3. Qué es robusto (convergente en N = 1600/3200/6400 y auditado)

1. **Horizonte.** Toda configuración que colapsa forma superficie atrapada. Las gaussianas anchas y poco densas no colapsan: se dispersan.

2. **Saturación QT y rebote central.**
   - El centro llega a αψ_c = 0.98–0.9999: es el «núcleo de de Sitter» de la QT, con ψ → 1/α.
   - Rebota exactamente donde ρ_c tiene su pico. Para s = 2.0 y α = 0.1: u = 19.000, ρ_c = 8.44 ≈ 970 V0, igual a las tres resoluciones.

3. **El rebote no deja χ en la meseta con V dominante.**
   - χ_c se dispara hasta 1.5–7 en unidades G = 1, es decir, 8–35 M̄: valores «de meseta».
   - Pero la energía es cinética (fracción ≥ 0.97) y ρ_c ~ 10²–10⁴ V0.

4. **Recompresión.**
   - En los dos casos donde el dominio fiable llega lo bastante lejos (s = 2.0, ambos α), el núcleo rebotado se vuelve a comprimir a Δu ≈ 0.8 del rebote (Δτ ≈ 0.2–0.3).
   - Para α = 0.1, ρ_c sube de ≈ V0 a 18.7 ≈ 2000 V0, por encima de la densidad del rebote.
   - Interpretación (no demostrada en detalle): es el manto que sigue cayendo.

5. **Horizonte interno y crecimiento de masa de tipo «inflación de masa» (mass-inflation-like).**
   - Se forma una capa en r ≈ √α (ψ ≈ 1/α, f ≈ 0): el borde del núcleo de de Sitter.
   - En ella la masa generalizada crece. En el dominio fiable llega a ×2–10 M_ADM en las configuraciones más masivas, y sin límite fuera de él.
   - Al mismo tiempo σ → −∞, el corrimiento al rojo exponencial de las coordenadas nulas.
   - Tiene la firma esperada de la inflación de masa en el horizonte interno de un agujero negro regular, pero **no está demostrado**: falta comprobar una ley ln M ∝ κ₋v convergente y seguir un invariante de curvatura (ver addendum). Ahí el código deja de converger, y eso fija el final del dominio fiable.

6. **El L del top-hat R = 6 no es un lanzamiento por rebote.**
   - Es el propio parche inicial: ya es mayor que su radio de Hubble y está dominado por V (ρ = 0.76 V0, χ = 0.3).
   - Infla unos 1.7 e-folds y rueda a cero, sin rebote, detrás de un horizonte que se forma fuera (r ≈ 8.3).
   - Es un parche inflacionario super-Hubble detrás de una superficie atrapada, **análogo cinemáticamente** a las construcciones de universo hijo (Blau–Guendelman–Guth usan paredes y una geometría de agujero de gusano que aquí no se demuestra). No lo produce la QT.

**Límite estructural (auditoría §6).** El tiempo propio fiable tras el rebote es solo de 0.08–1.09, y L exige ≥ 1/H_inf = 3.7. Un lanzamiento tardío **no es detectable por construcción**. Lo que sí está respaldado: no hay dominio de V antes de que termine el dominio fiable, y en s = 2.0 hay recompresión antes de cualquier lanzamiento.

## 4. Comparación homogénea (C45b-0, parche cerrado solo de χ)

**Configuración.**
- Port auditado de C45: reproduce la fila de referencia con 2 rebotes y N_inf = 59.6097, a 10⁻⁶.
- Enmienda declarada: la ventana temporal de C45 cortaba las inflaciones directas antes de los 60 e-folds; se amplió en 70/H_inf, y la referencia conserva su ventana original.
- Pasada de nube «rápida»: A·H_inf ∈ {0.96, 1.67}, 8 fases, al = 1e5, 120 s por run.

**Resultado:**
- **A·H = 0.96.** Ningún «fértil» según la definición de C45. En 4 de 8 fases el parche **infla directamente desde el turnaround, sin rebote**: 26 y 54 e-folds, con χ_0 ≈ 1.9–3 M̄. Hay 2 fases con 4 rebotes y poca inflación (1.7 e-folds) y 2 timeouts.
- **A·H = 1.67.** Las 8 fases agotan los 120 s.
- **Veredicto formal de la regla preregistrada: FAIL** (fertilidad 0 y 0). Pero con 10 de 16 timeouts no lo considero informativo. Lo decide el banco completo en WSL (J080, 600 s por run).

**Lectura.** En el caso homogéneo, un parche solo de χ con la densidad de C59a y el campo en el hombro tiende a inflar sin rebote, igual que el top-hat R = 6. El mecanismo de C45 (rebote → lanzamiento) necesita que **todo el parche** rebote a la vez. Una estrella que colapsa no lo hace: rebota el centro en una región de tamaño ~√α ≪ 1/H_inf, y el resto sigue cayendo.

## 5. Qué significa para ballena 4

1. **Lo que funciona.** La cadena «núcleo de χ → horizonte → rebote QT» ocurre de forma robusta y sin piezas a mano. La QT hace exactamente lo que prometía: curvatura acotada y rebote del centro.

2. **Lo que no se ve.** Que el rebote **lance** inflación.
   - Para lanzar, la región rebotada tendría que crecer un factor ~1/(H_inf√α): unas 12 veces con α = 0.1 y 21 con α = 0.03. Así la energía cinética caería ~10⁶ y V dominaría en una región mayor que 1/H_inf.
   - Lo observado es una recompresión en Δτ ≈ 0.2–0.3.
   - La hipótesis de homogeneidad de C45 (rebote coherente de todo el parche cerrado) **no la cumple** un núcleo estelar.

3. **Tensión nueva y concreta.** El lanzamiento por rebote requiere un colapso casi homólogo, en el que todas las capas rebotan a la vez, a una escala ≳ 1/H_inf. Una estrella de χ está concentrada en el centro y no colapsa así.

4. **El horizonte interno parece inestable** (crecimiento de masa de tipo inflación de masa, por confirmar), como en otros agujeros negros regulares. Cualquier versión «rebote hacia un agujero blanco o universo hijo» tiene que sobrevivir a eso.

## 6. Salidas posibles (para discutir; nada implementado)

- **(a) LTB-QT con polvo + χ.** En QT con polvo, cada capa evoluciona como un FRW cerrado. Es semianalítico, sin horizonte interno numérico. Permite medir directamente la **dispersión de tiempos de rebote** entre capas frente a 1/H_inf para perfiles de densidad realistas. Es la prueba limpia de cuándo el rebote es coherente, y la más barata.
- **(b) Código que cruce el horizonte interno.** Gauge adaptado o escisión de la capa r ≈ √α, para seguir la región anti-atrapada (rama de agujero blanco) que ya aparece en las simulaciones. Es caro, y el resultado puede ser simplemente que la inflación de masa la destruye.
- **(c) Cambiar el progenitor.** La inflación directa de parches super-Hubble (top-hat R = 6 y C45b-0) sugiere un camino sin rebote, análogo a los universos hijos de GR. Necesita condiciones iniciales que la cadena no fija hoy (el problema de Farhi–Guth), así que no lo recomiendo sin un argumento físico para ese parche.

**Mi recomendación:** (a) primero. Si incluso con colapso casi homólogo la dispersión de rebotes supera 1/H_inf, el mecanismo de lanzamiento por rebote queda descartado para núcleos estelares.

## 7. Enmiendas y límites declarados

- El protocolo C45b-dn se escribió tras las sondas exploratorias listadas. El criterio L cuenta también el estado inicial; por eso el top-hat R = 6 da L sin rebote. Se informa, sin reinterpretar el veredicto.
- `tests/test_dn_flat.py`: el test de orden se sustituyó por un test de exactitud, a raíz de la auditoría.
- C45b-0: se amplió la ventana temporal; la pasada de nube tiene muchos timeouts y el veredicto lo dará J080.
- Debilidades latentes señaladas por la auditoría, que no actúan en estos datos:
  - un H = NaN se trata como 0;
  - el H comóvil está mal condicionado cuando φ_t ≈ 0;
  - el criterio R es laxo tras un punto de retorno.

  En s = 2.0 la recompresión es real: supera la densidad del rebote.
- La pasada «fina» (N = 1600/3200/6400) va en J080. La auditoría ya la hizo para s = 2.0, α = 0.1: el dominio fiable pasa de 22.7 a 23.8 con los mismos eventos.

## 8. Addendum (3 de octubre, noche): revisión externa y pasos siguientes

La revisión externa mantiene **INCONCLUSIVE** y propone cambios que adopto para las rondas siguientes (C45b-R1). No reinterpreto el banco ya corrido:

1. **L solo tras B.** El lanzamiento QT se define como L posterior a un rebote (B ≺ L_B). El top-hat R = 6 queda como control de «inflación directa».
2. **Criterio regional, no central.** Hace falta un radio R_vac dentro del cual V > E_cin + E_grad, Θ > 0 y |χ| > χ_end, con **R_vac·H_inf ≳ 1**. Un punto en r = 0 no infla un universo.
3. **«Inflación de masa» pasa a ser «de tipo inflación de masa»** hasta demostrar ln M ∝ κ₋v convergente y una divergencia de un invariante. Si eso ocurre antes de que exista una región anti-atrapada útil, se convierte en **veto**.
4. **Detector de rebote geométrico.** Usar H_c = ȧ_c/a_c del desarrollo regular r(τ,ξ) = a_c(τ)ξ + O(ξ³), en lugar del H comóvil del campo, que está mal condicionado cuando φ_t ≈ 0. Mientras tanto, la proyección usa también el H de la ecuación de continuidad, H = −ρ̇/(3χ̇²), como control.
5. **Proyección a la cuenca de C45 y prueba LTB-QT.** Ya hechas: ver `INFORME_C45B_LTB.md`.

## 9. Reproducir

```
cd claude_dual/c45b
python -m tests.run_all
python -m dn_bank.run_c45b_dn          # ~7 min (nube)
python -m dn_bank.run_c45b_dn --fine   # N=1600/3200/6400 (WSL)
python -m hom.run_hom                  # banco homogéneo completo (WSL, horas)
python make_figures.py
```
