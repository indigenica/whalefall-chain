# U17f · El borde superior de la banda de F en todo el intervalo del colapso: 0.7·F₀ seguro; F₀ solo si el colapso extrapola por debajo de X ≈ 194

*19 de octubre de 2026, 22:43 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U17f.md`, con su enmienda 1. Los dos son anteriores a cualquier fila.
- *Código:* `lya/u17f_class.py`.
- *Resultado:* `lya/results/u17f/U17F_RESULT.json`; el registro de la ejecución, en `lya/results/u17f/u17f_run.log`.
- *Ejecución:* en la nube, 29 filas de CLASS con 4 procesos, ≈ 5 h en total.
  - Hubo un reinicio del contenedor con 19 filas hechas; esas salieron de la caché al relanzar.
  - La primera pasada completa falló al escribir el JSON (§5).

## 0. En una frase

**Veredicto preregistrado: BORDE_0.7_CONFIRMADO; F₀ PERMITIDO CON X ≤ 194.**
- **En todo el intervalo del colapso que deja U19b** (X del universo ≈ 170–320), la banda llega al menos a 0.7·F₀.
- **Si el colapso extrapola por debajo de X ≈ 194, pasa también F₀**, solo por forma. Ahí estarían el benchmark y LISA.
- **La extrapolación ∝ ℓ con el término de cuerdas da X ≈ 195:** cae justo en la frontera.

## 1. Los controles

- **La tubería:** la familia Y = 1.70 es la de U17e. En F = 1, 0.7 y 0.5·F₀ reproduce el δA de U17e exactamente: 0.4678, 0.2832 y 0.1540.
- **La dominancia:** la fila de control (Y = 1.98, F = 0.35·F₀) pasa la puerta primaria por CLASS, con δA = 0.090. Las filas F ≤ 0.35 que pasan Lyα por dominancia sobre las Y+ de U17e quedan validadas.
- **La regla de frío trivial:** se aplicó en F ≤ 0.125–0.25, según la familia. Su control es el de U17e, que pasó.

## 2. El borde superior

δA y la puerta primaria de U12 en las filas altas. La puerta usa la forma cuando el ajuste es representable, que lo fue en todas.

| Familia | X en F₀ | F₀ | 0.85·F₀ | 0.7·F₀ | 0.5·F₀ | F_max |
|---|---|---|---|---|---|---|
| Y = 1.05 | 168 | 0.347, **pasa por forma** | 0.262, pasa | 0.185, pasa | 0.098, pasa | **1** |
| Y = 1.20 | 193 | 0.378, **pasa por forma** | 0.295, pasa | 0.208, pasa | 0.112, pasa | **1** |
| Y = 1.36 | 219 | 0.409, falla | 0.328, pasa | 0.233, pasa | 0.126, pasa | **0.85** |
| Y = 1.70 | 275 | 0.468, falla | 0.383, pasa | 0.283, pasa | 0.154, pasa | **0.85** |
| Y = 1.98 | 321 | 0.510, falla | 0.421, falla | 0.323, pasa | 0.177, pasa | **0.7** |
| *(ii) Y = 1.0* | 160 | 0.335, pasa | 0.252, pasa | 0.177, pasa | 0.093, pasa | *1* |
| *(ii) Y = 0.5* | 78 | 0.183, pasa | 0.136, pasa | 0.094, pasa | 0.047, pasa | *1* |

**Referencias de δA:** la estricta es 0.234 (1 keV mezclado) y la ancha llega a 0.336 (térmico de 5.3 keV).

**Lo que hay que decir de F₀ en Y = 1.05 y 1.20:**
- **Pasa solo por el criterio de forma.** El ajuste a la familia térmica mezclada queda dentro de la cota estricta.
- **Su δA (0.35 y 0.38) supera incluso la referencia ancha.** La puerta preregistrada de U12 es esa, y U21 P pasó F₀ igual. Pero es un pase por un criterio y no por los dos.
- **0.85·F₀ pasa** en todas las familias hasta X ≈ 275, en las de X ≲ 220 con δA dentro de la referencia ancha.

**Las otras puertas no recortan arriba en ninguna fila:**
- la isocurvatura de U20b: f_iso ≤ 6.7×10⁻³ en el peor caso, frente al umbral de 8.4×10⁻³;
- ΔN_eff;
- las cuerdas del CMB.

## 3. El borde inferior

**Sin cambios: 0.09·F₀** en todas las familias primarias.
- La superradiancia estelar excluye 0.0625 y 0.045·F₀, y 0.03·F₀ vuelve a pasar. Es la misma banda partida que en U23.
- En Y = 1.98, 0.0625·F₀ cae en la «frontera» y no en «excluida», así que su tramo contiguo baja a 0.0625.
- En la familia informativa Y = 0.5 el borde sube a 0.125·F₀.

**La expectativa del ESPEC** («el borde inferior sube a ≈ 0.125·F₀ en el extremo favorable») **no se cumple en las primarias.** Solo pasa en Y = 0.5.

## 4. Con la cota secundaria de supermasivos (m_G < 1.7×10⁻¹⁷ eV excluida)

| X en F₀ | ≈ 320–275 | ≈ 220–170 | ≈ 160 (ii) | ≈ 78 (ii) |
|---|---|---|---|---|
| F_max | 0.5 | 0.7 | 0.85 | 1 |

**Con esta cota, F₀ queda excluido en todo el intervalo primario:** en Y ≤ 1.20, m_G = 6–8×10⁻¹⁸ eV. Solo lo salva la lectura (ii) más baja.

## 5. Notas de registro

- **El fallo al escribir.** La primera pasada completa calculó todo y falló al escribir `U17F_RESULT.json`. U23 (`f_cd`) devuelve `inf` cuando el modo 211 del supermasivo no es superradiante, y el JSON estricto no admite `inf`.
  - Se corrigió el código de salida, que escribe `inf` como `null` (`_finite`), y se volvió a ejecutar desde la caché de CLASS.
  - No cambia ninguna fila ni ningún criterio. Está anotado en el registro.
- **La meseta:** se usó la suma (enmienda 1). Con el polo, U26 mostró que m_G se mueve ≤ 0.4 % y que los veredictos de CLASS se arrastran.
- **La expectativa escrita antes** («BORDE_0.7_CONFIRMADO o BORDE_SUBE; si en el extremo favorable F₀ pasa por forma, con el calificador») se cumple en su segunda forma.

## 6. Qué significa

1. **La banda combinada llega al menos a 0.7·F₀** con cualquier extrapolación del colapso que crezca con ℓ. El 0.09–0.7·F₀ de U17e, U20b y U23 queda confirmado en todo el intervalo, no solo en X ≈ 275.
2. **El borde superior depende de X:**
   - 0.7·F₀ con X ≈ 320;
   - 0.85·F₀ con X ≈ 220–275;
   - F₀ con X ≲ 194.
3. **La pregunta del benchmark** (F₀, y con él la señal de cuerdas en LISA) **queda reducida a dos cosas:**
   - si el X del universo cae por debajo de ≈ 194. La extrapolación ∝ ℓ, la que mejor apoyan los datos de U19b, cae en la frontera;
   - si se da crédito a la cota secundaria de supermasivos, que excluye F₀ en todo el intervalo.

   Y aun pasando, F₀ pasa solo por forma.
4. **Una tercera masa de U19b** con el mismo protocolo relativo (ℓ_G ≈ 3.6 a 768³; `docs/NOTA_TERCERA_MASA_U19.md`, en la rama `red-jax-gpu`) **es ahora lo que decide F₀**. Es la condición que fijaba esa nota.
   - Si el exponente de X₀ resulta ≲ 1, X_univ quedaría por debajo de ≈ 194.
   - Si resulta ≳ 1.2, quedaría por encima.

## 7. Lo que U17f no resuelve

- **Cuál es la extrapolación correcta,** y con ella si F₀ vuelve.
- **La cota de supermasivos,** que sigue siendo secundaria y poco robusta (U23).
- **La normalización C de la isocurvatura,** y la forma de la población, que sale de ℓ₁ = 4.3.
