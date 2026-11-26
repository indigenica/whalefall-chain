"""U67c (docs/preregistros/ESPEC_U67c.md): RP3 networks at 256^3 with the U67 RD_SPARSE protocol, xi_i = 0.05 (row A) and
0.5 (row B), seeds 601-604.  Skips the networks already done.
    nohup python3 -B network/u67c_produce.py > network/results/u67c/produce.out 2>&1 &"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt

N, TAU_END, DTAU, NSNAP, ELL_D = 256, 128.0, 0.25, 30, 3.8
SEEDS = (601, 602, 603, 604)
ROWS = {"A": 0.05, "B": 0.5}
OUT = Path(__file__).resolve().parent / "results" / "u67c" / "runs"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for sd in SEEDS:
        for row, xi_i in ROWS.items():
            f = OUT / f"RP3_N{N}_{row}_s{sd}.json"
            if f.exists():
                continue
            t0 = time.time()
            r = Lt.run(1, N, TAU_END, sd, sfac=1.0, xi_init=xi_i, dtau=DTAU, nsnap=NSNAP, spectra_every=1000,
                       damp_ell=ELL_D, log=lambda m: None)
            rows = r["rows"] if isinstance(r, dict) else r
            f.write_text(json.dumps(dict(row=row, seed=sd, N=N, tau_end=TAU_END, xi_init=xi_i, damp_ell=ELL_D,
                                         rows=[dict(ell=x["ell"], xi=x["xi"]) for x in rows],
                                         runtime_s=time.time() - t0), default=float))
            print(f"done {row} s{sd} {time.time() - t0:.0f}s xi_end={rows[-1]['xi']:.3f}", flush=True)


if __name__ == "__main__":
    main()
