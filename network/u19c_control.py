"""U19c enmienda 2 (b): the control network on the A100 (JAX) against the same network with kern.c in the cloud.
PASS if xi agrees within 2 % at every snapshot (where the kern.c network has no strings left, xi = 0, |xi_A100| <= 2 % of the
maximum xi) and X_lat (u19b_analyze.per_run) within 2 %.

    python -B network/u19c_control.py <A100 RP3_N256_lG2.507_sparse_s91.json>   ->  U19C_CONTROL.json in U19C_OUT (default network/results/u19c)"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import u19b_analyze as A                                                      # noqa: E402
import u19c_analyze as C                                                      # noqa: E402
import u19c_produce as P                                                      # noqa: E402

TAG, TOL = "RP3_N256_lG2.507_sparse_s91", 0.02
REPO = Path(__file__).resolve().parent / "results" / "u19c"                  # the kern.c reference, in the repository


def main(path):
    ref, gpu = json.loads((REPO / "runs" / f"{TAG}.json").read_text()), json.loads(Path(path).read_text())
    v = json.loads((REPO / "U19C_VALIDATION.json").read_text())["masses"]["2.507"]; val = C.val_u19b_format(v, 2.507)
    xr, xg = (np.array([r["xi"] for r in d["rows"] if "xi" in r]) for d in (ref, gpu))
    if len(xr) != len(xg): dxi = float("inf")
    else:                                                                     # where kern.c has no strings left, |xi_gpu| against max xi
        dxi = float(max(np.max(np.abs(xg[xr > 0] / xr[xr > 0] - 1)), np.max(np.abs(xg[xr == 0]), initial=0.0) / xr.max()))
    Xr, Xg = A.per_run(ref, val)["X_lat"], A.per_run(gpu, val)["X_lat"]; dX = abs(Xg / Xr - 1)
    out = dict(spec="docs/preregistros/ESPEC_U19c.md (enmienda 2)", tag=TAG, backend_ref=ref.get("backend", "kern.c"), backend_gpu=gpu.get("backend"),
               n_snap=[len(xr), len(xg)], max_rel_xi=dxi, X_lat_ref=Xr, X_lat_gpu=Xg, rel_X_lat=dX, identical=bool(dxi == 0 and dX == 0),
               PASS=bool(dxi <= TOL and dX <= TOL))
    (P.OUT / "U19C_CONTROL.json").write_text(json.dumps(out, indent=1, default=float))
    print("CONTROL", "PASS" if out["PASS"] else "FAIL", json.dumps(out, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
