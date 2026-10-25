"""Figures of the unified closure + WDM-bound sensitivity of the allowed efficiency (post-processing, declared)."""
import json, math
from ucommon import RESULTS, FIGS, save
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from s4_dark import lambda_wdm_reference

STY = {("ir_delta", 10): ("C0", "-"), ("ir_delta", 30): ("C0", "--"), ("ir_delta", 100): ("C0", ":"),
       ("energy_logflat", 10): ("C3", "-"), ("energy_logflat", 30): ("C3", "--"), ("energy_logflat", 100): ("C3", ":")}


def main():
    d = json.loads((RESULTS / "s5_closure.json").read_text())
    rows = d["rows"]; g = d["gates"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for K, alpha in ((0, 1.0), (10, 0.35)):
        for (spec, x0), (c, ls) in STY.items():
            sub = sorted([r for r in rows if r["spectrum"] == spec and r["x0"] == x0 and r["K_sw"] == K and r["P1"]], key=lambda r: r["eps"])
            e = [r["eps"] for r in sub]
            lab = f"{spec}, x0={x0}" if K == 0 else None
            ax[0].loglog(e, [r["m_G_eV"] for r in sub], color=c, ls=ls, alpha=alpha, label=lab)
            ax[1].loglog(e, [r["lambda_em_over_WDM"] for r in sub], color=c, ls=ls, alpha=alpha)
            ax[2].loglog(e, [r["DeltaNeff_CMB"] for r in sub], color=c, ls=ls, alpha=alpha)
            ok = [r for r in sub if r["conservative"]]
            ax[0].scatter([r["eps"] for r in ok], [r["m_G_eV"] for r in ok], s=10, color=c, alpha=alpha)
    ax[0].axhspan(1e-23, g["M_LYA"], color="k", alpha=0.15, label="Lyα (m < 2e-20 eV)")
    for lo, hi in g["SR_WINDOWS"]:
        ax[0].axhspan(lo, hi, color="orange", alpha=0.12, hatch="//")
    ax[0].set_ylim(1e-22, 1e-13); ax[0].set_xlabel("eficiencia de la fuente ε (respecto a GHN U(1))"); ax[0].set_ylabel("m_G [eV] con Ωh² = 0.12")
    ax[0].set_title("Masa del pseudo-Goldstone\n(naranja: superradiancia, Stott 2020; la ventana SMBH exige f_a ≳ 2e16 GeV)", fontsize=9)
    ax[0].legend(fontsize=7, loc="lower left")
    ax[1].axhline(1, color="k", lw=1); ax[1].set_xlabel("ε"); ax[1].set_ylabel("λ_fs(emitida) / λ_fs(WDM 5.3 keV)")
    ax[1].set_title("Flujo libre de la componente emitida\n(puerta P5: por debajo de 1, o f_em ≤ 0.1)", fontsize=9)
    ax[2].axhline(g["DNEFF_CMB"], color="k", lw=1); ax[2].axhline(0.06, color="gray", lw=1, ls="--")
    ax[2].text(1.2e-4, 0.07, "CMB-S4 (2σ)", fontsize=8, color="gray"); ax[2].text(1.2e-4, 0.32, "Planck 95 %", fontsize=8)
    ax[2].set_xlabel("ε"); ax[2].set_ylabel("ΔN_eff (momento 3P a T = 0.26 eV)"); ax[2].set_title("Radiación oscura", fontsize=9)
    for a in ax: a.grid(alpha=0.3, which="both", lw=0.3)
    fig.suptitle("U1 · cierre condicional: líneas K_sw = 0 (opacas) y 10 (claras); puntos = nodos permitidos (P1–P5)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIGS / "u1_cierre.png", dpi=140); plt.close(fig)
    # sensitivity of eps_max to the WDM bound (post-processing; interpolates lambda_em(eps) of the scan)
    sens = {}
    for mk in (3.5, 5.3, 6.0):
        lw, _ = lambda_wdm_reference(mk)
        for K in (0, 1):
            for (spec, x0) in STY:
                sub = sorted([r for r in rows if r["spectrum"] == spec and r["x0"] == x0 and r["K_sw"] == K and r["P1"]], key=lambda r: r["eps"])
                ok = [r["P1"] and r["P2"] and r["P3"] and r["P4"] and (r["f_em"] <= 0.1 or r["lambda_em_Mpc"] < lw) for r in sub]
                n = next((i for i, o in enumerate(ok) if not o), len(ok))
                emax = sub[n - 1]["eps"] if n > 0 else None
                sens[f"{mk}keV|{spec}|x0={x0}|K={K}"] = emax
    save("s5_sensitivity_wdm.json", dict(note="largest eps of the contiguous allowed range from eps=0 (P1-P5) when the WDM reference mass of P5 is changed; 0.25-dex grid of the scan; post-processing, not preregistered", eps_max=sens))
    print("figure figures/u1_cierre.png; eps_max sensitivity:")
    for k, v in sens.items(): print(" ", k, v)


if __name__ == "__main__":
    main()
