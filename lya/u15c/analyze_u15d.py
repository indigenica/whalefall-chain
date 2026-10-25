"""U15d (docs/preregistros/ESPEC_U15d.md): the current upper edge of the F band against the full Lya likelihood, with the
collapse population measured in U19b (the U17f rows).  Per-case criterion of U15c unchanged (p(s = 1)/p_max >= e^-2 with 30
bins on the prior, stable with 20/45 bins and leaving one chain out, labels only after convergence and the U15b validation);
the helpers are those of analyze_u15c.py.  Verdict in two parts: F0 (D1, D2) and the worst-case edge at X ~ 320 (D3, D4).

    python -B lya/u15c/analyze_u15d.py ~/u15b/chains        ->  lya/results/u15d/U15D_RESULT.json"""
import json, os, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT / "lya" / "results" / "u15d"
sys.path.insert(0, str(ROOT / "lya" / "u15b")); sys.path.insert(0, str(ROOT / "lya" / "u15c"))
import analyze as AN                                                      # noqa: E402
import analyze_u15c as C                                                  # noqa: E402

D1, D2, D3, D4 = "U15d_F1_Y1.20", "U15d_F1_Y1.70", "U15d_F0.85_Y1.98", "U15d_F0.7_Y1.98"
CASES = [D1, D2, D3, D4]
S_EXTRA = ("p1_over_pmax", "p1_over_pmax_bins", "p1_over_pmax_jackknife", "posterior_prior_bins", "dchi2_min_s1", "n_points_near_1", "p0_over_pmax")


def case(base, c, val):
    f = os.path.join(base, c)
    if not os.path.isdir(f): return dict(error="missing")
    r = AN.case(f)
    if "error" in r: return r
    pr = C.prior_range(f)
    if pr is None or not (pr[0] <= 1.0 <= pr[1]): return dict(error=f"prior range of s_tail unknown or without s = 1: {pr}")
    names, ch = C.chains_of(f); j = 2 + names.index("s_tail")
    P = np.vstack(ch); s, w = P[:, j], P[:, 0]
    p1 = {b: C.ratio_at_1(s, w, *pr, b)[0] for b in (20, 30, 45)}
    _, edges, hist = C.ratio_at_1(s, w, *pr, 30)
    jack = [C.ratio_at_1(np.concatenate([x[:, j] for k, x in enumerate(ch) if k != i]),
                         np.concatenate([x[:, 0] for k, x in enumerate(ch) if k != i]), *pr, 30)[0] for i in range(len(ch))]
    near = np.abs(s - 1.0) < 0.05
    labels = {C.classify(x) for x in list(p1.values()) + jack}
    status = (None if not AN.converged(r) or val is not True else
              labels.pop() if len(labels) == 1 and None not in labels else "INDECISO (binning o muestreo)")
    r.update(prior=list(pr), p1_over_pmax=p1[30], p1_over_pmax_bins=p1, p1_over_pmax_jackknife=jack, p0_over_pmax=float(hist[0] / hist.max()),
             n_points_near_1=int(near.sum()), dchi2_min_s1=float(2 * (P[near, 1].min() - P[:, 1].min())) if near.any() else None,
             posterior_prior_bins=dict(edges=edges.tolist(), density=hist.tolist()), status=status)
    if status is None:                                   # nothing about s before convergence and validation (U15c enmienda 1)
        for k in AN.S_KEYS + S_EXTRA: r.pop(k, None)
    return r


def verdict(st):
    """ESPEC_U15d: '<F0 part>; <edge part>'"""
    und = lambda *xs: any(x is None or x.startswith("INDECISO") for x in xs)
    if und(st[D1], st[D2]): f0 = "F0_INCONCLUSO (p(1)/p_max en la frontera)"
    elif st[D1] == "EXCLUIDA" and st[D2] == "PERMITIDA": f0 = "F0_INCONCLUSO (no monótono en X)"
    elif st[D1] == "EXCLUIDA": f0 = "F0_EXCLUIDO_LYA_COMPLETO"
    elif st[D2] == "PERMITIDA": f0 = "F0_PERMITIDO_LYA_COMPLETO"
    else: f0 = "F0_PERMITIDO_SOLO_CON_X_195"
    if und(st[D3], st[D4]): edge = "BORDE_INCONCLUSO (p(1)/p_max en la frontera)"      # ESPEC_U15d: any INDECISO case makes its part inconclusive
    elif st[D3] == "PERMITIDA": edge = "BORDE_0.85_SOSTENIDO"
    elif st[D4] == "PERMITIDA": edge = "BORDE_0.7_SOSTENIDO"
    else: edge = "BORDE_BAJA_DE_0.7"
    return f"{f0}; {edge}"


def main(argv):
    base = argv[1] if len(argv) > 1 else str(Path.home() / "u15b" / "chains"); OUT.mkdir(parents=True, exist_ok=True)
    u15b = ROOT / "lya" / "results" / "u15b" / "U15B_RESULT.json"
    val = json.loads(u15b.read_text()).get("validation_passed") if u15b.exists() else None
    man = ROOT / "lya" / "u15b" / "tables" / "manifest_u15d.json"
    out = dict(spec="docs/preregistros/ESPEC_U15d.md", u15b_validation_passed=val,
               tables=json.loads(man.read_text()) if man.exists() else None, cases={c: case(base, c, val) for c in CASES})
    for c in CASES:
        r = out["cases"][c]; print(c, {k: r.get(k) for k in ("status", "p1_over_pmax", "R_minus_1_max", "R_minus_1_split_max", "n_chains", "error")}, flush=True)
    conv = all(AN.converged(out["cases"][c]) for c in CASES)
    if val is None: v = "PENDIENTE (validación de U15b todavía no disponible)"
    elif val is not True: v = "INCONCLUSO (validación de U15b fallida)"
    elif not conv: v = "PENDIENTE (cadenas sin converger: R-1 o R-1 por mitades >= 0.03, o un número de cadenas distinto de 4)"
    else: v = verdict({c: out["cases"][c].get("status") for c in CASES})
    out["converged"] = conv; out["verdict_U15d"] = v
    (OUT / "U15D_RESULT.json").write_text(json.dumps(AN.finite(out), indent=1, default=float, allow_nan=False))
    print("U15d:", v, flush=True)


if __name__ == "__main__":
    main(sys.argv)
