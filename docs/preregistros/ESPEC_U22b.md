# U22b · ¿Cuánto nacimiento queda con la inflación de polo? La rejilla completa de C45h · Preregistro

*16 de octubre de 2026 (Madrid; la hora que cuenta es la de su commit). Escrito antes de ejecutar U22b. Se ejecuta en la nube.*

## Motivo

U22 (POLO_CAMBIA_NACIMIENTO) solo rehízo las 90 ejecuciones de S1: los 20 hijos viables y 70 controles. Con el polo, los aterrizajes cambian (~3 veces más bajos) y el umbral eterno sube de 416 a 541. Falta la cifra que cuenta para el eslabón: **qué fracción de los nacimientos fértiles acaba eterna**, sobre toda la rejilla de C45h (1100 filas con frac_rad = 0, 900 fértiles).

## Qué se ejecuta (`chain/u22b_births.py`)

- **Las 1100 filas de C45h con frac_rad = 0,** con el puerto C de C45:
  - modo 0 (suma de 10 niveles), que debe reproducir lo guardado;
  - modo 4 (polo, C = 0.3516, el mismo que U22).
- **Por potencial** se calculan, con el V₀ común de S2 para cada potencial (el de la suma, 416; el del polo, 541):
  - la fracción fértil;
  - el umbral eterno χ_et;
  - **f_ES:** la fracción de fértiles en la ventana sep ∈ {10⁻⁴, 10⁻⁵, 10⁻⁶} con χ_land > χ_et y E_phys,max < 1, la definición de S1.
- Además, la desagregación por sep.

## Veredicto

Se mide con R = f_ES(polo)/f_ES(suma):
- **NACIMIENTO_ROBUSTO:** R ≥ 0.5;
- **NACIMIENTO_REDUCIDO:** 0.1 ≤ R < 0.5;
- **NACIMIENTO_FRAGIL:** R < 0.1;
- **INVALIDO:** el modo 0 no reproduce lo guardado (diferencia relativa en χ_land > 10⁻⁸).

## Lo que se espera, escrito antes de calcular

Por los 70 controles de U22 (7 eternos frente a ~35), R ≈ 0.2, es decir, **NACIMIENTO_REDUCIDO**.

## Lo que U22b no resuelve

- **El criterio de perturbaciones de C45j** (ζ₀,max), que da la «fracción viable» de 3.6 %, no se puede rehacer: el código de C45j no está en el repositorio, solo sus resultados. f_ES es la fracción antes de ese filtro.
