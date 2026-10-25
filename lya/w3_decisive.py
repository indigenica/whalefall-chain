"""W3b (complement to the preregistered W3, declared after seeing that W3 had been evaluated on saturated rows, delta A ~ 0.99):
convergence of delta A on the DECISIVE rows (B2/B3, energy_logflat, K_sw = 0, central R): production settings (40 nodes, l_max 10,
phase 0.5) vs (64 nodes, l_max 16, phase 0.25) and a finer k grid.  Output: results/W3_decisive.json"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import run_u3a as R
import boltz as B


def main():
    d = json.loads((HERE / "results" / "U3A_RESULT.json").read_text())
    s6 = json.loads((R.RESULTS / "s6_measured_closure.json").read_text())
    U2 = json.loads((R.UNI / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
    Bm = R.benchmark(); F, Trh = Bm["F_GeV"], Bm["T_Gamma_GeV"]
    picks = [r for r in s6["rows"] if r.get("P1") and r["K_sw"] == 0 and r["R_label"] == "central" and r["spectrum"] == "energy_logflat"
             and r["branch"][:2] in ("B2", "B3")]
    out = []
    for kh, tag in ((R.KH, "k30"), (np.geomspace(0.05, 200.0, 44), "k44")):
        R.KH = kh
        D0, _ = B.solve(kh * B.h); R.D0 = D0; Pc = B.p_cdm_shape(kh * B.h, D0)
        for r in picks:
            pk = R.S6.source_eps(F, r["m_G_eV"], Trh, r["x0"], R.eps_fn_for(r, U2), spectrum=r["spectrum"])
            x = pk["q"] / pk["m_GeV"]
            for nq, lw, ph in ((40, 10, 0.5), (64, 16, 0.25)):
                sp, _ = B.from_packets(x, pk["nc"], 0.12 * r["f_em"] / B.h ** 2, nq=nq, lmax=lw)
                d1, _ = B.solve(kh * B.h, sp, phase_max=ph)
                dA = B.deltaA(kh, (d1 / D0) ** 2, Pc)
                out.append(dict(grid=tag, branch=r["branch"], x0=r["x0"], nq=nq, lmax=lw, phase=ph, deltaA=dA)); print(out[-1], flush=True)
    base = {(o["branch"], o["x0"]): o["deltaA"] for o in out if o["grid"] == "k30" and o["nq"] == 40}
    rel = max(abs(o["deltaA"] / base[(o["branch"], o["x0"])] - 1) for o in out)
    res = dict(note=__doc__, cases=out, max_rel_change=rel, passed=bool(rel < 0.02))
    (HERE / "results" / "W3_decisive.json").write_text(json.dumps(res, indent=1, default=float)); print("W3b", rel, res["passed"])


if __name__ == "__main__":
    main()
