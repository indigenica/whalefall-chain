# U7c · ¿Sigue pasando Lyα la emisión calibrada con la R_Γ de las redes densas? · Preregistro

*9 de octubre de 2026, ≈22:00 (Madrid). Escrito después de ver U9 y antes de calcular ninguna fila de U7c.*

## Motivo

U9 (`docs/informes/INFORME_U9.md`) muestra que el cociente de emisión depende de la densidad inicial de las redes:
- R_Γ = 1.39 con ξ_i = 0.5 (U6);
- R_Γ = 1.69 con ξ_i = 2 (U9), en la misma ventana;
- moviendo el comienzo de la ventana, va de 1.50 a 2.09 con ξ_i = 2, y de 1.22 a 1.41 con ξ_i = 0.5.

Las filas de U7 se calcularon con R_Γ = 1.40 ± 0.13 de U4b, hasta 1.53. U7c repite las mismas filas con valores de R_Γ que cubren esa dispersión de protocolo.

## Filas

Las de U7b, sin cambios en nada más:
- familias B2κ (calibrada con ν_ax) y B2ζ (ξ* = 1.49);
- x₀ ∈ {10, 30} y K_sw ∈ {0, 1};
- la misma cadena S6;
- CLASS con el parche y 45 nodos;
- las mismas puertas de U4a y U5: la forma estricta como puerta primaria y δA ≤ 0.234;
- κ(x₀) de U7a sin recalcular.

R_Γ ∈ {1.39, 1.69, 2.09}: el central de U6, el de U9 con la ventana preregistrada y el máximo del barrido de ventana de U9. Son 24 filas.

`python -B lya/u7_class.py --workers 24 --u7c` → `lya/results/u7/U7C_RESULT.json`.

## Veredicto (sobre B2κ)

- **CALIBRADA_ROBUSTA:** con R_Γ = 2.09 hay al menos una fila permitida con x₀ = 10 y al menos una con x₀ = 30.
- **CALIBRADA_SENSIBLE_A_R:** no es robusta, pero hay alguna fila permitida con R_Γ = 1.69.
- **CALIBRADA_EXCLUIDA_CON_REDES_DENSAS:** ninguna fila de B2κ permitida con R_Γ = 1.69.

B2ζ se informa sin veredicto.

## Lo que U7c no resuelve

- **La pendiente espectral también depende del protocolo:** q ≈ 1.03 con ξ_i = 0.5 frente a ≈ 0.74 con ξ_i = 2. U7c fija el espectro en log-plano (q = 1), como U7.
- **Si q < 1 a ℓ grande**, la emisión deja de contar (U5f) y el sector oscuro pasa sin firma tibia, sea cual sea R_Γ. U7c solo explora la dirección peligrosa, la de más emisión.
