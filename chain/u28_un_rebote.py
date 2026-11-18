"""Informative re-reading, no verdict (14 November 2026): how many of the viable children of U28 (dust parent with an
inhomogeneous chi) and of the audit of U36 (parent of pure radiation) bounce only once, and so do not rely on the
extrapolation of the coherence from one bounce to several (division of the tolerated amplitude by the number of
bounces). No new physics is computed: U28's own selection and per-child tolerances (chain/u28_c45k.py, unchanged) and
the stored result of the audit are counted by number of bounces.
    python3 -B chain/u28_un_rebote.py      ->  results/u28_un_rebote.json (~2 min)"""
import collections, json
from ucommon import ROOT, save
import u28_c45k as U28


def dust():
    ch, nfert, _ = U28.children()
    R = json.loads((ROOT / "results" / "u28_c45k.json").read_text())
    assert R["n_fertile"] == nfert
    out = {}
    for cls, g in (("E", "E-v"), ("D", "D-v")):
        s = R["class_LH1"][g]; by = collections.Counter(); n = [0, 0]
        for r in ch[cls]:                                   # the rule of u28_c45k.main, L.H = 1
            k = U28.tag(r) + "|LH1"
            own = R["per_child"].get(k, {})
            rl, rh = (own["rho_lo"], own["rho_hi"]) if "rho_lo" in own else (s["median_rho_lo"], s["median_rho_hi"])
            vlo, vhi = r["zeta0_LH1"] * rl >= U28.ZETA0, r["zeta0_LH1"] * rh >= U28.ZETA0
            n[0] += vlo; n[1] += vhi; by[(r["n_b"], 0)] += vlo; by[(r["n_b"], 1)] += vhi
        assert (n[0], n[1]) == (R["viable"][cls]["n_viable_lo"], R["viable"][cls]["n_viable_hi"])
        out[cls] = dict(n_viable_lo=n[0], n_viable_hi=n[1],
                        one_bounce_lo=by[(1, 0)], one_bounce_hi=by[(1, 1)],
                        f_one_bounce_lo=by[(1, 0)] / nfert, f_one_bounce_hi=by[(1, 1)] / nfert,
                        by_n_bounces={str(b): [by[(b, 0)], by[(b, 1)]] for b in sorted({k[0] for k in by})})
    return dict(n_fertile=nfert, classes=out)


def radiation():
    A = json.loads((ROOT / "audit" / "u36_coherencia" / "resultados.json").read_text())
    v = [c for c in A["per_child"].values() if c["cls"] == "D" and c.get("viable")]
    assert len(v) == A["branches"]["D"]["n_viable"]
    by = collections.Counter(c["n_b"] for c in v)
    one = [c for c in v if c["n_b"] == 1]
    return dict(n_fertile=A["n_fertile"], n_viable_D=len(v), by_n_bounces={str(k): by[k] for k in sorted(by)},
                one_bounce=len(one), f_one_bounce=len(one) / A["n_fertile"],
                one_bounce_margin_min=min(c["zeta0_max"] for c in one) / 4.6e-5)


if __name__ == "__main__":
    out = dict(note="informative re-reading, no verdict (chain/u28_un_rebote.py)", dust_U28=dust(), radiation_audit_U36=radiation())
    save("u28_un_rebote.json", out)
    print(json.dumps(out, indent=1))
