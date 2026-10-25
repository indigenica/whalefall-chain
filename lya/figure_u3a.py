"""Figure of U3a: delta A of every closure row (S6) against the strict and wide Lyman-alpha thresholds."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def main():
    d = json.loads((HERE / "results" / "U3A_RESULT.json").read_text())
    rows = [r for r in d["rows"] if r["K_sw"] == 0]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    cols = {"B1": "C3", "B2": "C0", "B3": "C2"}; mk = {"ir_delta": "o", "energy_logflat": "s"}
    for r in rows:
        b = r["branch"][:2]; x = r["f_em"]
        ax.scatter(x, r["deltaA"], c=cols[b], marker=mk[r["spectrum"]], s=28 if r["R_label"] == "central" else 12,
                   alpha=1.0 if r["R_label"] == "central" else 0.45, edgecolors="k" if r["x0"] == 30 else "none")
    ax.axhline(d["deltaA_strict"], color="k", lw=1); ax.axhline(d["deltaA_wide"], color="k", lw=1, ls="--")
    ax.text(1e-3, d["deltaA_strict"] + 0.01, "umbral estricto (mezcla en el borde, 1 keV)", fontsize=7)
    ax.text(1e-3, d["deltaA_wide"] + 0.01, "umbral amplio (WDM térmico 5.3 keV)", fontsize=7)
    for b, c in cols.items():
        ax.scatter([], [], c=c, label={"B1": "B1 GHV (ξ ∝ ℓ)", "B2": "B2 Hindmarsh (ξ*=1.2)", "B3": "B3 tendencia medida"}[b])
    ax.scatter([], [], marker="o", c="gray", label="ir_delta"); ax.scatter([], [], marker="s", c="gray", label="energy_logflat")
    ax.set_xscale("log"); ax.set_xlabel("fracción de materia oscura emitida f_em"); ax.set_ylabel("δA (Murgia et al. 2017)")
    ax.set_title(f"U3a · Lyα con transferencia de Boltzmann (K_sw = 0; grande = R central, borde negro = x0 = 30)\n"
                 f"veredicto estricto: {d['verdict_strict']['verdict']} · amplio: {d['verdict_wide']['verdict']}", fontsize=9)
    ax.legend(fontsize=7, loc="lower right"); ax.grid(alpha=0.3, lw=0.3)
    fig.tight_layout(); fig.savefig(HERE / "results" / "u3a_lya.png", dpi=140); plt.close(fig)
    print("figure lya/results/u3a_lya.png")


if __name__ == "__main__":
    main()
