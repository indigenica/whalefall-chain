"""Figure of the cuscuton birth for v4 of the WHALEFALL preprint (main-v4.tex).  Reads only stored results:

    results/u62e_cl_class.json, results/u64_cubico_unitario.json, data/planck2018/COM_PowerSpect_CMB-TT-full_R3.01.txt

Run from the repository root:  python3 -B paper/zenodo_whale/figures/make_figures_v4.py
Output: paper/zenodo_whale/figures/fig_cuscuton.pdf"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
COLS = {30.0: "#2166ac", 100.0: "#67a9cf", 1000.0: "#b2182b"}


def J(*parts):
    return json.loads(ROOT.joinpath(*parts).read_text())


def main():
    e = J("results", "u62e_cl_class.json"); c = J("results", "u64_cubico_unitario.json")
    pl = np.loadtxt(ROOT / "data" / "planck2018" / "COM_PowerSpect_CMB-TT-full_R3.01.txt")
    pl = pl[(pl[:, 0] >= 2) & (pl[:, 0] <= 29)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ell = np.array(e["ell"])
    a.errorbar(pl[:, 0], pl[:, 1], yerr=[pl[:, 2], pl[:, 3]], fmt="o", ms=2.5, color="k", lw=0.7, label="Planck 2018 (Commander)")
    ref = np.array(e["hosts"][0]["D_ref"])
    a.plot(ell, ref, color="#888888", lw=1.2, label="no imprint")
    seen = set()
    for h in e["hosts"]:
        AH = h["par"]["AH"]; st = h["states"]["A2"]
        lab = f"$AH={AH:.0f}$ ($\\Delta\\chi^2={st['dchi2_planck']:.1f}$)"
        a.plot(ell, ref * np.array(st["ratio"]), color=COLS[AH], lw=1.0, ls="-" if AH not in seen else "--", label=lab)
        seen.add(AH)
    a.set_xscale("log"); a.set_xticks([2, 3, 5, 10, 20, 30]); a.set_xticklabels(["2", "3", "5", "10", "20", "30"]); a.minorticks_off(); a.set_xlabel(r"$\ell$"); a.set_ylabel(r"$D_\ell\ [\mu{\rm K}^2]$"); a.set_xlim(1.8, 31)
    a.set_title("(a) visible cuscuton children against Planck", fontsize=8); a.legend(loc="upper right", frameon=False)
    for h in c["hosts"]:
        n = [p["n"] for p in h["per_n"]]; r3 = [p["rho3"] for p in h["per_n"]]
        AH = h["par"]["AH"]
        b.loglog(n, r3, "-o", ms=2, lw=0.8, color=COLS.get(AH, "#999999"))
    for AH, col in ((30.0, COLS[30.0]), (100.0, COLS[100.0]), (1000.0, COLS[1000.0]), (None, "#999999")):
        b.plot([], [], "-o", ms=2, lw=0.8, color=col, label=f"$AH={AH:.0f}$" if AH else "$AH=48$, $300$")
    b.axhline(1, color="k", lw=0.8, ls=":"); b.legend(loc="upper left", frameon=False)
    b.set_xlabel(r"harmonic $n$"); b.set_ylabel(r"cubic/quadratic, $\rho_3$")
    b.set_title("(b) cubic order in the unitary gauge (8 children)", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_cuscuton.pdf"); print("fig_cuscuton.pdf")


if __name__ == "__main__":
    main()
