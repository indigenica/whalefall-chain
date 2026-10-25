"""Figures for the WHALEFALL (whale 5) Zenodo preprint.  Reads only result files already in the repository:

    network/results/u4/runs/*.json, network/results/u4/U6_RESULT.json
    lya/results/u4a/U4A_RESULT.json, lya/results/u5/U5_RESULT.json, lya/results/u7/U7_RESULT.json

Run from the repository root:  python -B paper/zenodo_whale/figures/make_figures.py
Output: paper/zenodo_whale/figures/fig_network.pdf, fig_lya.pdf"""
import glob, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "network"))
RUNS = ROOT / "network" / "results" / "u4" / "runs"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
C_RP3, C_U1, C_U1B = "#b2182b", "#2166ac", "#67a9cf"


def load(pattern):
    return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pattern)))]


def series(run, key="xi"):
    return np.array([(r["ell"], r[key]) for r in run["rows"] if key in r]).T


def fig_network():
    u6 = json.loads((ROOT / "network" / "results" / "u4" / "U6_RESULT.json").read_text())
    rp = load("RP3_N576_s6*_xi0.5.json"); u1 = load("U1_N576_s6*_xi0.5.json")
    u1b = load("U1_N768_s6*_xi0.5.json"); u1ic = load("U1_N768_s61_xi0.25.json") + load("U1_N768_s61_xi2.json")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
    a = ax[0]
    for i, r in enumerate(rp):
        a.plot(*series(r), color=C_RP3, lw=0.8, alpha=0.8, label=r"RP$^3$, $576^3$ (6 seeds)" if i == 0 else None)
    for i, r in enumerate(u1):
        a.plot(*series(r), color=C_U1, lw=0.8, alpha=0.8, label=r"U(1), $576^3$ (4 seeds)" if i == 0 else None)
    for i, r in enumerate(u1b):
        a.plot(*series(r), color=C_U1B, lw=0.8, ls="--", label=r"U(1), $768^3$ (4 seeds)" if i == 0 else None)
    for i, r in enumerate(u1ic):
        a.plot(*series(r), color="0.5", lw=0.7, ls=":", label=r"U(1), $768^3$, $\xi_i=0.25,\,2$" if i == 0 else None)
    for i, r in enumerate(load("RP3_N576_s6*_xi2.json")):
        a.plot(*series(r), color="#f4a582", lw=0.9, ls="-.", label=r"RP$^3$, $576^3$, $\xi_i=2$ (U9)" if i == 0 else None)
    for i, r in enumerate(load("U1_N576_s6*_xi2.json")):
        a.plot(*series(r), color="#92c5de", lw=0.9, ls="-.", label=r"U(1), $576^3$, $\xi_i=2$ (U9)" if i == 0 else None)
    w = u6["window"]; a.axvspan(*w, color="0.92", zorder=0)
    a.set_xlim(2.5, 6.05); a.set_ylim(0, 2.6); a.set_xlabel(r"$\ell=\ln(m_r/H)$"); a.set_ylabel(r"$\xi = L\tau^2/(4V)$")
    a.legend(loc="upper left", frameon=False, fontsize=5.8); a.set_title("(a) string density", fontsize=9)

    a = ax[1]
    g = np.array(u6["R_Gamma"]["grid"]); pe = np.array(u6["R_Gamma"]["per_ell"]); v, s = u6["R_Gamma"]["value"], u6["R_Gamma"]["boot_std"]
    a.plot(g, pe, "o-", color=C_RP3, ms=3, lw=0.8, label=r"$R_\Gamma(\ell)$, seed means")
    a.axhspan(v - s, v + s, color=C_RP3, alpha=0.15); a.axhline(v, color=C_RP3, lw=0.8)
    gx = np.array(u6["R_xi"]["grid"]); px = np.array(u6["R_xi"]["per_ell"]); vx, sx = u6["R_xi"]["value"], u6["R_xi"]["boot_std"]
    a.plot(gx, px, "s-", color="k", ms=3, lw=0.8, label=r"$R_\xi(\ell)$")
    a.axhspan(vx - sx, vx + sx, color="k", alpha=0.08); a.axhline(1.0, color="0.6", lw=0.6, ls=":")
    a.set_ylim(0.8, 2.4); a.set_xlabel(r"$\ell$"); a.set_ylabel("RP$^3$ / U(1)")
    a.text(0.03, 0.05, rf"$R_\Gamma={v:.2f}\pm{s:.2f}$" + "\n" + rf"$R_\xi={vx:.2f}\pm{sx:.2f}$", transform=a.transAxes, fontsize=7.5)
    a.legend(loc="upper left", frameon=False); a.set_title("(b) emission and density ratios", fontsize=9)

    a = ax[2]
    q = u6["q_RP3"]
    for i, ps in enumerate(q["per_seed"]):
        ps = np.array(ps); a.plot(ps[:, 0], ps[:, 1], "o", color=C_RP3, ms=2.5, alpha=0.6, label="RP$^3$ seeds" if i == 0 else None)
    lx = np.linspace(4.6, q["ell_end"], 20); fit = q["q0_at_ell5"] + q["q1_per_ell"] * (lx - 5)
    a.plot(lx, fit, color=C_RP3, lw=1.2, label=rf"fit: $q(5)={q['q0_at_ell5']:.2f}\pm{q['q0_std']:.2f}$")
    a.axhline(1.0, color="k", lw=0.7, ls="--", label=r"$q=1$ (Correia et al.)")
    a.axhline(1.025, color="0.4", lw=0.6, ls=":", label=r"$q_{\rm crit}=1.025$ (U5f)")
    a.set_xlabel(r"$\ell$"); a.set_ylabel(r"emission slope $q$"); a.set_ylim(0.4, 2.1)
    a.legend(loc="upper left", frameon=False); a.set_title("(c) RP$^3$ emission spectrum", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "fig_network.pdf"); plt.close(fig)


def fig_lya():
    u4 = json.loads((ROOT / "lya" / "results" / "u4a" / "U4A_RESULT.json").read_text())
    u5 = json.loads((ROOT / "lya" / "results" / "u5" / "U5_RESULT.json").read_text())
    u7 = json.loads((ROOT / "lya" / "results" / "u7" / "U7_RESULT.json").read_text())
    st, wd = u4["deltaA_strict"], u4["deltaA_wide"]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    a = ax[0]
    sty = {("B1_GHN", "energy_logflat"): ("v", "0.35", r"B1 ($\xi\propto\ell$)"),
           ("B1_GHN", "ir_delta"): ("v", "0.35", None),
           ("B2_standard_scaling", "energy_logflat"): ("o", C_RP3, r"B2 ($\xi_*$ const.), log-flat"),
           ("B3_trend_1_over_ell", "energy_logflat"): ("s", "#ef8a62", "B3 (measured trend), log-flat"),
           ("B2_standard_scaling", "ir_delta"): ("x", C_U1, "B2/B3, IR spectrum"),
           ("B3_trend_1_over_ell", "ir_delta"): ("x", C_U1, None)}
    for key, (mk, c, lab) in sty.items():
        rows = [r for r in u4["rows"] if (r["branch"], r["spectrum"]) == key and "deltaA" in r]
        a.scatter([r["f_em"] for r in rows], [r["deltaA"] for r in rows], marker=mk, s=14, color=c, label=lab, lw=0.8,
                  facecolors="none" if mk in "os" else None)
    b4 = [r for r in u5["rows"] if r.get("part") == "U5a" and "deltaA" in r]
    a.scatter([r["f_em"] for r in b4], [r["deltaA"] for r in b4], marker="^", s=12, color="#7b3294", label=r"B4 (growing $\xi$, U5a)",
              facecolors="none", lw=0.7)
    a7 = [r for r in u7["rows"] if "deltaA" in r]
    a.scatter([r["f_em"] for r in a7], [r["deltaA"] for r in a7], marker="*", s=22, color="#1b7837", label=r"calibrated to $\nu_{\rm ax}$ (U7)")
    a.axhline(st, color="k", lw=0.8, ls="--"); a.axhline(wd, color="k", lw=0.6, ls=":")
    a.text(0.012, st * 1.06, rf"strict $\delta A={st:.3f}$", fontsize=7); a.text(0.012, wd * 1.06, rf"wide $\delta A={wd:.3f}$", fontsize=7)
    a.set_xscale("log"); a.set_xlim(0.01, 1.2); a.set_ylim(0, 1.05)
    a.set_xlabel(r"emitted fraction $f_{\rm em}$"); a.set_ylabel(r"$\delta A$ (CLASS)")
    a.legend(loc="upper left", frameon=False, fontsize=6.3); a.set_title(r"(a) Ly$\alpha$ gate, all CLASS rows", fontsize=9)

    a = ax[1]
    for lab, mk, c in (("-1s", "v", "0.5"), ("central", "o", C_RP3), ("+1s", "^", "k")):
        rows = [r for r in u5["rows"] if r["branch"] == "B2" and r.get("part") in ("U5f", "ref") and r["R_label"] == lab
                and r["x0"] == 10 and r["K_sw"] == 0 and "deltaA" in r]
        rows.sort(key=lambda r: r["q"])
        a.plot([r["q"] for r in rows], [r["deltaA"] for r in rows], marker=mk, ms=4, lw=0.8, color=c,
               label={"-1s": r"$R_\Gamma-1\sigma$", "central": r"$R_\Gamma$ central", "+1s": r"$R_\Gamma+1\sigma$"}[lab])
    a.axhline(st, color="k", lw=0.8, ls="--"); a.axvline(u5["q_crit_B2"], color="0.4", lw=0.6, ls=":")
    q = json.loads((ROOT / "network" / "results" / "u4" / "U6_RESULT.json").read_text())["q_RP3"]
    a.axvspan(q["q0_at_ell5"] - q["q0_std"], q["q0_at_ell5"] + q["q0_std"], color=C_RP3, alpha=0.12, label=r"RP$^3$ lattice $q(\ell=5)$")
    a.set_xlabel(r"emission slope $q$ (B2, $x_0=10$, $K_{\rm sw}=0$)"); a.set_ylabel(r"$\delta A$")
    a.set_ylim(-0.02, 0.7); a.legend(loc="upper left", frameon=False); a.set_title(r"(b) sensitivity to $q$ (U5f)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "fig_lya.pdf"); plt.close(fig)


def fig_bracket():
    """U10b: the scaling attractor bracketed from below (U11, drained sparse starts) and from above (U10b, late dense starts)"""
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharex=True)
    for a, mod, title in ((ax[0], "RP3", r"(a) RP$^3$"), (ax[1], "U1", "(b) U(1)")):
        lo, hi = load(f"{mod}_N576_s6*_xi0.5_dr.json"), load(f"{mod}_N576_s7*_ab.json")
        for i, r in enumerate(lo):
            a.plot(*series(r), color=C_U1, lw=0.7, alpha=0.7, label=r"from below ($\xi_i=0.5$, U11)" if i == 0 else None)
        for i, r in enumerate(hi):
            a.plot(*series(r), color=C_RP3, lw=0.9, ls="-" if r["xi_init"] == 3 else "--",
                   label=(r"from above ($\xi_i=3$, $\ell_i=3.5$, U10b)" if r["xi_init"] == 3 else r"from above ($\xi_i=6$)")
                   if (i == 0 or r["xi_init"] == 6) else None)
        a.axvspan(4.36, 5.61, color="0.92", zorder=0)
        a.set_xlim(3.4, 5.8); a.set_ylim(0, 3.2 if mod == "RP3" else 2.4); a.set_xlabel(r"$\ell$"); a.set_ylabel(r"$\xi$")
        a.set_title(title, fontsize=9)
    ax[0].legend(loc="upper right", frameon=False, fontsize=6.5)
    fig.tight_layout(); fig.savefig(OUT / "fig_bracket.pdf"); plt.close(fig)


if __name__ == "__main__":
    fig_network(); fig_lya(); fig_bracket(); print("figures written to", OUT)
