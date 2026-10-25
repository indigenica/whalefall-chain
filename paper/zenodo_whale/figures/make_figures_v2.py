"""Figures for version 2 of the WHALEFALL (whale 5) Zenodo preprint (main-v2.tex).  Reads only result files already in
the repository:

    network/results/u19d/U19D_RESULT.json, network/results/u19c/U19C_RESULT.json
    lya/results/u17f/U17F_RESULT.json, lya/results/u17g/U17G_RESULT.json
    lya/results/u15c/U15C_RESULT.json, lya/results/u15d/U15D_RESULT.json

Run from the repository root:  python -B paper/zenodo_whale/figures/make_figures_v2.py
Output: paper/zenodo_whale/figures/fig_collapse.pdf, fig_band.pdf, fig_lyafull.pdf
(fig_bracket.pdf, shared with version 1, comes from make_figures.py.)"""
import json, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
C_RP3, C_U1, C_GREY, C_OK, C_NO = "#b2182b", "#2166ac", "#888888", "#1b7837", "#b2182b"
L_UNIV = 90.0


def J(*parts):
    return json.loads(ROOT.joinpath(*parts).read_text())


def wls(x, y, s, design):
    """weighted least squares, as network/u19c_analyze.wls: coefficients, covariance, chi2, dof"""
    A = design(np.asarray(x, float)); W = 1 / np.asarray(s, float) ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None])); b = cov @ (A.T @ (W * y))
    return b, cov, float(np.sum(W * (y - A @ b) ** 2)), len(y) - A.shape[1]


def fig_collapse():
    d = J("network", "results", "u19d", "U19D_RESULT.json"); c = J("network", "results", "u19c", "U19C_RESULT.json")
    ms = sorted(d["masses"], key=int)
    l1 = np.array([d["masses"][k]["ell_1"] for k in ms]); X0 = np.array([d["masses"][k]["affine"]["X0"] for k in ms])
    sX = np.array([d["masses"][k]["affine"]["X0_err"] for k in ms]); pl = d["power_law"]; ys = pl["strings_term"]
    lp = pl["pivot_ell_1"]
    beta, cov, _, _ = wls(np.log(l1 / lp), np.log(X0), sX / X0, lambda x: np.column_stack([np.ones_like(x), x]))
    assert abs(beta[1] - pl["p"]) < 1e-9
    bp, _, _, _ = wls(l1, X0, sX, lambda x: x[:, None])
    bl, _, _, _ = wls(l1, X0, sX, lambda x: np.column_stack([np.ones_like(x), x]))
    models = [(r"$X_0\propto\ell_1$", lambda l: bp[0] * l, "--"), (r"linear in $\ell_1$", lambda l: bl[0] + bl[1] * l, ":"),
              (rf"power, $p={pl['p']:.2f}\pm{pl['sigma_p']:.2f}$", lambda l, p=pl["p"]: np.exp(beta[0] + p * np.log(l / lp)), "-")]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
    a = ax[0]
    cm = [k for k in c["masses_in_fit"]]
    a.errorbar([c["masses"][str(k)]["ell_1"] - 0.012 for k in cm], [c["masses"][str(k)]["affine"]["X0"] for k in cm],
               [c["masses"][str(k)]["affine"]["X0_err"] for k in cm], fmt="o", ms=3, color=C_GREY, mfc="white", capsize=2,
               label="4 networks per mass")
    a.errorbar(l1, X0, sX, fmt="o", ms=4, color=C_RP3, capsize=2, label="12 networks per mass")
    ll = np.linspace(3.4, 4.6, 50)
    for lab, f, ls in models:
        a.plot(ll, f(ll), ls, color="k", lw=0.9, label=lab)
    a.set_xlabel(r"$\ell_1=\ln(m_r/H)$ at the onset of the collapse"); a.set_ylabel(r"soft momentum $X_0$ (string-independent part)")
    a.legend(frameon=False, loc="upper left"); a.set_title("(a) lattice, $256^3$--$576^3$", fontsize=9)
    a = ax[1]
    ll = np.geomspace(3.2, 110, 200)
    lo = np.exp(beta[0] + (pl["p"] - pl["sigma_p"]) * np.log(ll / lp)); hi = np.exp(beta[0] + (pl["p"] + pl["sigma_p"]) * np.log(ll / lp))
    a.fill_between(ll, lo, hi, color=C_RP3, alpha=0.15, lw=0)
    for lab, f, ls in models:
        a.plot(ll, f(ll), ls, color="k", lw=0.9)
    bc, _, _, _ = wls(l1, X0, sX, lambda x: np.ones((len(x), 1)))
    a.axhline(bc[0], color=C_GREY, lw=0.7, ls="-."); a.text(12, bc[0] * 0.8, "constant (disfavoured)", fontsize=7, color=C_GREY)
    a.errorbar(l1, X0, sX, fmt="o", ms=3, color=C_RP3, capsize=1.5)
    for k in ("X0_proporcional_a_l1", "X0_lineal_en_l1", "potencia_libre"):
        a.plot([L_UNIV], [d["models"][k]["X_univ"]], "s", ms=4, mfc="white", color="k")
    a.errorbar([L_UNIV * 1.08], [pl["X_univ"]], yerr=[[pl["X_univ"] - pl["X_lo"]], [pl["X_hi"] - pl["X_univ"]]], fmt="none",
               ecolor=C_RP3, capsize=2)
    a.text(L_UNIV * 0.93, 1700, rf"$X_{{\rm univ}}$ (+{ys:.0f})", fontsize=7, ha="right")
    a.axvline(L_UNIV, color=C_GREY, lw=0.6)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(3.2, 110); a.set_ylim(4, 2500)
    a.set_xlabel(r"$\ell$"); a.set_ylabel(r"$X_0(\ell)$; squares: $X_{\rm univ}$")
    a.set_title(r"(b) extrapolation to $\ell\simeq90$", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "fig_collapse.pdf"); plt.close(fig)


def fig_band():
    fam = []
    for f in (("lya", "results", "u17f", "U17F_RESULT.json"), ("lya", "results", "u17g", "U17G_RESULT.json")):
        r = J(*f)
        for k, b in r["bands"].items():
            Y = float(k[1:])
            if Y < 1.0: continue
            fam.append((162 * Y, b["F_max_over_F0"], b["F_min_contiguous_over_F0"], r["bands_with_SMBH"][k]["F_max_over_F0"]))
    fam.sort(); X = np.array([f[0] for f in fam]); up = np.array([f[1] for f in fam]); low = np.array([f[2] for f in fam])
    smbh = np.array([f[3] for f in fam])
    d = J("network", "results", "u19d", "U19D_RESULT.json"); pl = d["power_law"]
    u = J("lya", "results", "u15d", "U15D_RESULT.json")
    fig, a = plt.subplots(figsize=(4.6, 3.0))
    a.axvspan(pl["X_lo"], pl["X_hi"], color=C_RP3, alpha=0.08, lw=0)
    a.axvline(pl["X_univ"], color=C_RP3, lw=0.9)
    for k, ls in (("X0_proporcional_a_l1", "--"), ("X0_lineal_en_l1", ":")):
        a.axvline(d["models"][k]["X_univ"], color=C_RP3, lw=0.8, ls=ls)
    a.fill_between(X, low, up, step="mid", color=C_U1, alpha=0.18, lw=0, label="allowed (screening gates)")
    a.step(X, up, where="mid", color=C_U1, lw=1.2); a.step(X, low, where="mid", color=C_U1, lw=1.2)
    a.step(X, smbh, where="mid", color=C_U1, lw=0.9, ls="--", label="upper edge with the SMBH bound")
    for k, t in u["tables"].items():
        st = u["cases"][k]["status"]
        a.plot(t["X_univ"], t["F_over_F0"], "o" if st == "PERMITIDA" else "x", ms=6 if st != "PERMITIDA" else 5,
               mfc="none" if st == "PERMITIDA" else None, mew=1.3, color=C_OK if st == "PERMITIDA" else C_NO)
    a.plot([], [], "o", mfc="none", color=C_OK, label="full likelihood: allowed"); a.plot([], [], "x", color=C_NO, label="full likelihood: excluded")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(150, 1500); a.set_ylim(0.025, 1.25)
    a.set_yticks([0.03, 0.06, 0.1, 0.25, 0.5, 1]); a.set_yticklabels(["0.03", "0.06", "0.1", "0.25", "0.5", "1"])
    a.set_xticks([200, 300, 500, 1000]); a.set_xticklabels(["200", "300", "500", "1000"])
    a.set_xlabel(r"soft momentum of the collapse in the universe, $X_{\rm univ}$"); a.set_ylabel(r"$F/F_0$")
    a.minorticks_off(); a.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="lower left", fontsize=6.5)
    a.text(pl["X_univ"] * 1.03, 0.9, "central", fontsize=7, color=C_RP3)
    fig.tight_layout(); fig.savefig(OUT / "fig_band.pdf"); plt.close(fig)


def fig_lyafull():
    fig, a = plt.subplots(figsize=(4.6, 2.8))
    sets = [("u15c", "U15C_RESULT.json", {"U15c_B2_F0.5_r1": r"$0.5F_0$, earlier population", "U15c_B2_F0.35_r3": r"$0.35F_0$, earlier population"}, ":"),
            ("u15d", "U15D_RESULT.json", None, "-")]
    cols = {"U15d_F1_Y1.20": "#762a83", "U15d_F1_Y1.70": "#b2182b", "U15d_F0.85_Y1.98": "#ef8a62", "U15d_F0.7_Y1.98": "#1b7837"}
    for sub, f, labels, ls in sets:
        r = J("lya", "results", sub, f)
        for k, c in r["cases"].items():
            p = c["posterior_prior_bins"]; e = np.array(p["edges"]); dens = np.array(p["density"])
            s = 0.5 * (e[1:] + e[:-1]); y = dens / dens.max()
            if labels:
                a.plot(s, y, ls, color=C_GREY, lw=1.0, label=labels[k])
            else:
                t = r["tables"][k]
                a.plot(s, y, ls, color=cols[k], lw=1.2, label=(rf"${t['F_over_F0']:g}F_0$" if t["F_over_F0"] != 1 else "$F_0$") + rf", $X\approx{t['X_univ']:.0f}$")
    a.axhline(math.exp(-2), color="k", lw=0.7, ls="--"); a.text(0.9, math.exp(-2) * 1.1, r"$e^{-2}$", fontsize=7)
    a.axvline(1, color=C_GREY, lw=0.6)
    a.set_xlabel(r"amplitude $s$ of the non-cold component ($s=1$: the model, $s=0$: CDM)"); a.set_ylabel(r"$p(s)/p_{\max}$")
    a.set_xlim(0, 1.05); a.set_ylim(0, 1.1); a.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0.0, 0.2))
    fig.tight_layout(); fig.savefig(OUT / "fig_lyafull.pdf"); plt.close(fig)


if __name__ == "__main__":
    fig_collapse(); fig_band(); fig_lyafull()
    print("written:", ", ".join(str(OUT / n) for n in ("fig_collapse.pdf", "fig_band.pdf", "fig_lyafull.pdf")))
