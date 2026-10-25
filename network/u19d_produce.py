"""U19d (docs/preregistros/ESPEC_U19d.md): seeds 93-96 in both protocols at the four U19c masses.  256^3, 384^3 and 576^3
use u19c_produce.network (the l_G = 3.2 protocol relative to l_1); 512^3 uses u19b_produce.network, the protocol of the four
primary U19b networks.  Validations are those of U19c and U19b (not repeated).  The order is by seed, then mass, then
protocol, so that an interrupted production loses as little as possible under the complete-seed rule of the analysis.

    python -B network/u19d_produce.py produce 256,384      (cloud, kern.c)
    WHALEFALL_BACKEND=jax python -B network/u19d_produce.py produce 512,576      (A100; U19D_OUT=<Drive folder> optional)
    -> network/results/u19d/runs/<TAG>.json and network/results/u19d/runs.log, skip-if-present"""
import json, os, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt                                                          # noqa: E402
import u19b_produce as B                                                      # noqa: E402
import u19c_produce as P                                                      # noqa: E402

HERE = Path(__file__).resolve().parent; OUT = Path(os.environ.get("U19D_OUT", HERE / "results" / "u19d"))
SEEDS = (93, 94, 95, 96)
LG = {256: 2.507, 384: 2.912, 512: 3.2, 576: 3.318}


def table(N):
    """the RP3 m_eff table of the mass, from the validation that already passed (U19c, or U19b for 512^3)"""
    if N == 512:
        v = json.loads((HERE / "results" / "u19b" / "U19B_VALIDATION.json").read_text()); assert v.get("PASS")
        return v["meff"]["RP3_3.2"]
    v = json.loads((HERE / "results" / "u19c" / "U19C_VALIDATION.json").read_text())["masses"][f"{LG[N]:g}"]; assert v.get("PASS")
    return v["meff"]["RP3"]


def produce(Ns):
    R = OUT / "runs"; R.mkdir(parents=True, exist_ok=True); log = open(OUT / "runs.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    tabs = {N: table(N) for N in Ns}
    lg(f"START N={Ns} {time.strftime('%Y-%m-%d %H:%M:%S')} threads={Lt.lib().nthreads()} backend={P.backend()}")
    for seed in SEEDS:
        for N in Ns:
            for prot in ("sparse", "dense"):
                tag = f"RP3_N{N}_lG{LG[N]:g}_{prot}_s{seed}"
                if (R / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
                if N == 512:
                    res = B.network(1, seed, prot, 3.2, 512, lg, tabs[N]); res.update(backend=P.backend(), spec="docs/preregistros/ESPEC_U19d.md")
                else:
                    res = P.network(1, seed, prot, N, lg, tabs[N]); res["spec"] = "docs/preregistros/ESPEC_U19d.md"
                (R / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    if sys.argv[1] == "produce": produce(tuple(int(x) for x in sys.argv[2].split(",")))
