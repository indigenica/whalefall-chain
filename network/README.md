# network/ · redes 3D de cuerdas: RP³ frente a U(1) (U2, U3b)

> **Meseta (U26, 18 de octubre de 2026).** Las redes no dependen de la meseta de inflación: nada de esta carpeta cambia con `WHALEFALL_PLATEAU`.


- `kern.c`: kernel C con OpenMP, en era de radiación, con cuerdas físicas y leapfrog conforme.
  - Q RP³ (10 componentes) y control U(1) (2 componentes).
  - Preparación de núcleo opcional (`set_core`, U3b). Desactivada por defecto, en cuyo caso reproduce U2 exactamente (lo comprueba `tests/`).
- `lattice.py`: condiciones iniciales y medidas.
  - Densidad de cuerdas ξ por holonomía ℤ₂ o devanado.
  - Radiación Goldstone libre: 2 × energía cinética tangente fuera de los núcleos.
  - Espectro.
- `validate.py`: V1–V4 (ESPEC_U2 §2).
- `produce.py` (U2) y `produce_u3.py` (U3b): las redes de cada conjunto.
- `analyze.py`: R_ξ, R_Γ y sistemáticos (U2). `analyze_u3.py`: puertas G1 y G2 (U3b).
- `results/u2/runs/*.json`: las 20 redes de U2. `results/u3/runs/*.json`: el piloto y las 3 redes U(1) de U3b.
- Coste: una red de 256³ tarda ~4 min (U(1)) o ~20 min (RP³) con 2 núcleos y usa ~1.5 GB de RAM.

## U19b · η y la época del colapso de RP³ (ESPEC_U19b): producción en WSL

Antes de producir, la validación (V1–V5) y el piloto se hacen en la nube. La producción se niega a arrancar si `results/u19b/U19B_VALIDATION.json` no existe o no pasa.

```bash
git pull
cd network
export OMP_NUM_THREADS=12          # los núcleos que dejen libres las cadenas de U15b/U15c (44 en total)
nohup python -B u19b_produce.py produce 512 > results/u19b/produce.out 2>&1 &
tail -f results/u19b/runs.log      # una línea por instantánea: ell, xi, fracción de pared, sitios lejos del vacío, N_lin
python -B u19b_analyze.py          # cuando termine (o con parte de las redes: informa sin veredicto si faltan primarias)
```

- **Orden de las 12 redes,** con salto si ya existen (se puede parar y continuar):
  1. RP³ con ℓ_G = 3.2 (las 4 primarias);
  2. RP³ con ℓ_G = 2.4;
  3. U(1) con ℓ_G = 3.2.
- **Memoria:** ~30 GB de pico por red RP³ de 512³; ~8 GB por red U(1).
- **Tiempo:** el piloto de 192³ tardó 14 min con 4 núcleos. Escalado a 512³ hasta 3τ_A, sale:
  - con 12 hilos, ~6–8 h por red RP³ y ~2 h por red U(1);
  - con 40 hilos, ~2–3 h por red RP³.
  - En total, ~2–3 días con 12 hilos.
- **Resultados:** `results/u19b/runs/<TAG>.json` y `results/u19b/runs.log`. Hay que subirlos al repositorio, como los de U16.
- **El veredicto sale solo de las 4 redes primarias,** RP³ con ℓ_G = 3.2. Las demás se informan.
