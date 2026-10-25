"""U15c (docs/preregistros/ESPEC_U15c.md): p(s = 1)/p_max >= e^-2 (Delta chi^2 <= 4) for the living-branch amplitude s?
Enmienda 1 (14 Oct 2026, codex audit, before any chain): histogram bins fixed on the prior range read from log.param, s = 1
must lie inside it, per-case labels only after convergence (U15b Nota 6 rule) and validation, and a label that changes with
20/30/45 bins or leaving one chain out is INDECISO.
Chains are read with the U15b rules (lya/u15b/analyze.py: restart copies skipped, points before the last proposal update
dropped, first 30 % of the rest dropped); the U15b validation result is required.

    python -B lya/u15c/analyze_u15c.py ~/u15b/chains        ->  lya/results/u15c/U15C_RESULT.json"""
import glob, json, math, os, re, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]; OUT = ROOT / "lya" / "results" / "u15c"
sys.path.insert(0, str(ROOT / "lya" / "u15b"))
import analyze as AN                                                      # noqa: E402
import chainread                                                          # noqa: E402

CASES = ["U15c_B2_F0.5_r1", "U15c_B2_F0.35_r3"]                          # enmienda 3


def chains_of(folder):
    """chains read with the U15b rules (ESPEC_U15b Nota 9: restart copies skipped, per-run non-Markovian cut; first 30 % dropped)"""
    names = AN.paramnames(folder); out = []
    rows_all, _ = chainread.read(folder, len(names) + 2)
    for rows in rows_all:
        a = np.array(rows)
        if len(a) >= 20: out.append(a[int(0.3 * len(a)):])
    return names, out


def prior_range(folder):
    """flat prior of s_tail from MontePython's log.param: data.parameters['s_tail'] = [mean, min, max, sigma, scale, role]"""
    f = os.path.join(folder, "log.param")
    if not os.path.exists(f): return None
    for line in open(f):
        m = re.match(r"\s*data\.parameters\[['\"]s_tail['\"]\]\s*=\s*\[([^\]]*)\]", line)
        if m:
            v = [x.strip() for x in m.group(1).split(",")]
            try: return float(v[1]), float(v[2])
            except ValueError: return None
    return None


def ratio_at_1(s, w, lo, hi, bins):
    """p(s = 1)/p_max with the bins fixed on the prior range [lo, hi] (enmienda 1: no clipping to the data range)"""
    hist, edges = np.histogram(s, bins=bins, range=(lo, hi), weights=w, density=True)
    j = min(int(np.searchsorted(edges, 1.0, side="right") - 1), bins - 1)
    return float(hist[j] / hist.max()) if hist.max() > 0 else float("nan"), edges, hist


def classify(x):
    return None if not np.isfinite(x) else ("PERMITIDA" if x >= math.exp(-2) else "EXCLUIDA")


def main(argv):
    base = argv[1] if len(argv) > 1 else str(Path.home() / "u15b" / "chains"); OUT.mkdir(parents=True, exist_ok=True)
    out = dict(spec="docs/preregistros/ESPEC_U15c.md", cases={})
    u15b = ROOT / "lya" / "results" / "u15b" / "U15B_RESULT.json"
    val = json.loads(u15b.read_text()).get("validation_passed") if u15b.exists() else None
    out["u15b_validation_passed"] = val
    for c in CASES:
        f = os.path.join(base, c)
        if not os.path.isdir(f): out["cases"][c] = dict(error="missing"); continue
        r = AN.case(f)
        if "error" in r: out["cases"][c] = r; continue
        pr = prior_range(f)
        if pr is None or not (pr[0] <= 1.0 <= pr[1]):
            out["cases"][c] = dict(error=f"prior range of s_tail unknown or without s = 1: {pr}"); continue
        names, ch = chains_of(f); j = 2 + names.index("s_tail")
        P = np.vstack(ch); s, w = P[:, j], P[:, 0]
        p1 = {b: ratio_at_1(s, w, *pr, b)[0] for b in (20, 30, 45)}
        _, edges, hist = ratio_at_1(s, w, *pr, 30)
        jack = [ratio_at_1(np.concatenate([x[:, j] for k, x in enumerate(ch) if k != i]),
                           np.concatenate([x[:, 0] for k, x in enumerate(ch) if k != i]), *pr, 30)[0] for i in range(len(ch))]
        near = np.abs(s - 1.0) < 0.05
        dchi2 = float(2 * (P[near, 1].min() - P[:, 1].min())) if near.any() else None
        labels = {classify(x) for x in list(p1.values()) + jack}
        status = (None if not AN.converged(r) or val is not True else
                  labels.pop() if len(labels) == 1 and None not in labels else "INDECISO (binning o muestreo)")
        r.update(prior=list(pr), p1_over_pmax=p1[30], p1_over_pmax_bins=p1, p1_over_pmax_jackknife=jack,
                 n_points_near_1=int(near.sum()), dchi2_min_s1=dchi2, posterior_prior_bins=dict(edges=edges.tolist(), density=hist.tolist()),
                 status=status)
        if status is None:                          # enmienda 1: nothing about s before convergence and validation
            for k in AN.S_KEYS + ("p1_over_pmax", "p1_over_pmax_bins", "p1_over_pmax_jackknife", "posterior_prior_bins", "dchi2_min_s1", "n_points_near_1"): r.pop(k, None)
        out["cases"][c] = r
        print(c, {k: r.get(k) for k in ("status", "p1_over_pmax", "R_minus_1_max", "R_minus_1_split_max", "n_chains")}, flush=True)
    conv = all(AN.converged(out["cases"][c]) for c in CASES)
    st = [out["cases"][c].get("status") for c in CASES]
    if val is None: v = "PENDIENTE (validación de U15b todavía no disponible)"
    elif val is not True: v = "INCONCLUSO (validación de U15b fallida)"
    elif not conv: v = "PENDIENTE (cadenas sin converger: R-1 o R-1 por mitades >= 0.03, o un número de cadenas distinto de 4)"
    elif any(x is None or x.startswith("INDECISO") for x in st): v = "INCONCLUSO (p(1)/p_max en la frontera)"
    elif all(x == "PERMITIDA" for x in st): v = "RAMA_VIVA_PERMITIDA_LYA_COMPLETO"
    elif all(x == "EXCLUIDA" for x in st): v = "RAMA_VIVA_EXCLUIDA_LYA_COMPLETO"
    else: v = "RAMA_VIVA_DEPENDE_DE_ETA"
    out["converged"] = conv; out["verdict_U15c"] = v
    (OUT / "U15C_RESULT.json").write_text(json.dumps(AN.finite(out), indent=1, default=float, allow_nan=False))
    print("U15c:", v, flush=True)


if __name__ == "__main__":
    main(sys.argv)
