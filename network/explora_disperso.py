"""Exploratory, no verdict (docs/informes/NOTA_U67C_LINEAL.md, follow-up written 2026-11-22 20:28 Madrid before running): radiation-era
RP3 networks at 192^3 with the U67 RD_SPARSE protocol but ten times sparser (xi_i = 0.05, the density at birth estimated in the linear
note), seeds 101-102, against the stored U67 RD_SPARSE networks (xi_i = 0.5). Reads how fast a very sparse network approaches the attractor.
    python3 -B network/explora_disperso.py  ->  network/results/explora_disperso/"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
import u67_produce as P

OUT = Path(__file__).resolve().parent / "results" / "explora_disperso"; OUT.mkdir(parents=True, exist_ok=True)
for sd in (101, 102):
    f = OUT / f"RP3_N192_XI005_s{sd}.json"
    if f.exists(): continue
    t0 = time.time()
    r = Lt.run(1, P.N, P.TAU_END, sd, sfac=1.0, xi_init=0.05, dtau=P.DTAU, nsnap=P.NSNAP, spectra_every=1000, damp_ell=P.ELL_D, log=lambda m: None)
    f.write_text(json.dumps(dict(seed=sd, xi_init=0.05, rows=[dict(ell=x["ell"], xi=x["xi"]) for x in r["rows"]] if isinstance(r, dict) else r,
                                 runtime_s=time.time() - t0), default=float))
    print("done", sd, time.time() - t0, flush=True)
