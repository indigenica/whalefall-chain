"""U13 (docs/preregistros/ESPEC_U13.md): can the pseudo-Goldstone mass m_G come from the common action without new parameters?

1. O(4) invariance of every Q term of the declared action (U(q), S = Tr q^2 in xi R S/2, Tr(dq dq)), checked numerically with
   random Q and random O(4) elements, reflections included.
2. Tree level: the vacuum manifold is the orbit of a rank-1 projector (RP^3) and the Hessian there has exactly 3 zero modes.
3. Quantum-gravity control: m_G from V ~ M_Pl^2 f^2 exp(-c M_Pl/f) for c in {0.5, 1, 2, 3}, and the c that would be needed.
4. Informative: order n of a Z_N-allowed operator Q^n / M_Pl^(n-4) giving m_G = 1e-15 - 1e-14 eV with an O(1) coefficient.

    python3 action/u13_mG.py  ->  action/results_u13/U13_RESULT.json"""
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parent / "results_u13"
LAM, LP, ZQ = 320.0, 192.0, 400.0                       # docs/modelo/MODELO_ACCION_COMUN.md
MPL_GEV, F_GEV = 2.435e18, 1.070076182963697e15        # reduced Planck mass; F of the benchmark (f = F/2)
F_A = F_GEV / 2
TARGET_EV = (1e-15, 1e-14)
rng = np.random.default_rng(13)


def X0():
    """Q = X0 I with zero gradient of U: d/dX [LAM (4X-1)^2/4 + LP*4 (X^2-X)^2/4] = 0  ->  2LAM(4X-1) + 2LP(X^2-X)(2X-1) = 0."""
    r = np.roots([4 * LP, -6 * LP, 2 * LP + 8 * LAM, -2 * LAM])
    return float(r[np.abs(r.imag) < 1e-10].real[0])


def U(Q):
    A = Q @ Q - Q
    return LAM * (np.trace(Q) - 1) ** 2 / 4 + LP * np.trace(A @ A) / 4


def S(Q, x0): D = Q - x0 * np.eye(4); return ZQ * np.trace(D @ D)       # up to the constant Munit^2


def rand_sym():
    A = rng.normal(size=(4, 4)); return (A + A.T) / 2


def rand_O4():
    Qm, Rm = np.linalg.qr(rng.normal(size=(4, 4))); Qm = Qm @ np.diag(np.sign(np.diag(Rm)))
    if rng.random() < 0.5: Qm = Qm @ np.diag([1, 1, 1, -1])               # include reflections
    return Qm


def symmetry_check(n=2000):
    x0 = X0(); worst = dict(U=0.0, S=0.0, kinetic=0.0)
    for _ in range(n):
        Q = rand_sym(); dQ = [rand_sym() for _ in range(4)]; O = rand_O4()
        Qr = O @ Q @ O.T; dQr = [O @ d @ O.T for d in dQ]
        worst["U"] = max(worst["U"], abs(U(Qr) - U(Q)) / max(1.0, abs(U(Q))))
        worst["S"] = max(worst["S"], abs(S(Qr, x0) - S(Q, x0)) / max(1.0, abs(S(Q, x0))))
        k, kr = sum(np.trace(d @ d) for d in dQ), sum(np.trace(d @ d) for d in dQr)
        worst["kinetic"] = max(worst["kinetic"], abs(kr - k) / max(1.0, abs(k)))
    return worst


def basis():
    B = []
    for i in range(4):
        for j in range(i, 4):
            M = np.zeros((4, 4)); M[i, j] = M[j, i] = 1.0
            B.append(M / np.linalg.norm(M))
    return B


def vacuum_and_hessian():
    B = basis(); mat = lambda v: sum(c * b for c, b in zip(v, B))
    best = min((minimize(lambda v: U(mat(v)), rng.normal(size=10), method="BFGS", options=dict(gtol=1e-12)) for _ in range(20)), key=lambda r: r.fun)
    Qv = mat(best.x); ev = np.sort(np.linalg.eigvalsh(Qv))
    h = 1e-5; Hs = np.zeros((10, 10))
    for a in range(10):
        for b in range(10):
            ea, eb = np.eye(10)[a] * h, np.eye(10)[b] * h
            Hs[a, b] = (U(mat(best.x + ea + eb)) - U(mat(best.x + ea - eb)) - U(mat(best.x - ea + eb)) + U(mat(best.x - ea - eb))) / (4 * h * h)
    w = np.sort(np.linalg.eigvalsh(Hs))
    return dict(U_min=float(best.fun), Q_eigenvalues=ev.tolist(), hessian_eigenvalues=w.tolist(),
                n_zero_modes=int(np.sum(np.abs(w) < 1e-4 * np.max(np.abs(w)))))


def m_from_V(V_GeV4): return math.sqrt(V_GeV4) / F_A * 1e9           # eV, from V ~ m^2 f^2


def main():
    OUT.mkdir(parents=True, exist_ok=True); out = dict(spec="docs/preregistros/ESPEC_U13.md", X0=X0())
    out["symmetry_worst_rel"] = symmetry_check(); print("O(4) invariance, worst relative deviation:", out["symmetry_worst_rel"])
    out["tree_level"] = vacuum_and_hessian()
    print("vacuum Q eigenvalues", np.round(out["tree_level"]["Q_eigenvalues"], 6), "Hessian zero modes:", out["tree_level"]["n_zero_modes"])
    Lam = [math.sqrt(m * 1e-9 * F_A) for m in TARGET_EV]
    out["target"] = dict(m_G_eV=TARGET_EV, Lambda_keV=[l * 1e6 for l in Lam], V_over_f4=[(l / F_A) ** 4 for l in Lam])
    print("target Lambda = (m f)^1/2 =", [f"{x:.0f} keV" for x in out["target"]["Lambda_keV"]], "V/f^4 =", [f"{x:.1e}" for x in out["target"]["V_over_f4"]])
    r = MPL_GEV / F_A; qg = {}
    for c in (0.5, 1, 2, 3):
        lg = 2 * math.log10(MPL_GEV) + 2 * math.log10(F_A) - c * r / math.log(10)          # log10 V
        qg[str(c)] = dict(log10_V_GeV4=lg, log10_m_G_eV=0.5 * lg - math.log10(F_A) + 9)
    need = {str(m): (2 * math.log(MPL_GEV) + 2 * math.log(F_A) - 2 * math.log(m * 1e-9 * F_A)) / r for m in TARGET_EV}
    out["quantum_gravity"] = dict(M_Pl_over_f=r, by_c=qg, c_needed=need)
    print(f"M_Pl/f = {r:.0f};", {c: f"m_G ~ 10^{v['log10_m_G_eV']:.0f} eV" for c, v in qg.items()}, "c needed:", {k: round(v, 3) for k, v in need.items()})
    zn = {str(m): 4 + math.log((m * 1e-9) ** 2 / F_A ** 2) / math.log(F_A / MPL_GEV) for m in TARGET_EV}
    out["discrete_ZN_operator_order"] = zn; print("Z_N operator order n for c = 1:", {k: round(v, 1) for k, v in zn.items()})
    out["anomaly"] = "la acción declarada no tiene fermiones ni campos gauge acoplados a Q: O(4) no es anómala"
    sym_ok = max(out["symmetry_worst_rel"].values()) < 1e-10 and out["tree_level"]["n_zero_modes"] == 3
    out["routes"] = {
        "1 radiativas de xi R S": "NO_SALE" if sym_ok else "REVISAR",
        "2 acoplamiento sigma de chi": "NO_SALE (no contiene Q)",
        "3 instantones del sector Q": "NO_SALE" if sym_ok else "REVISAR",
        "4 U(q) y cinético (control)": "NO_SALE (masa cero a nivel árbol)" if sym_ok else "REVISAR",
        "5 gravedad cuántica (control)": "SALE_CON_HIPOTESIS" if any(0.5 <= v <= 3 for v in need.values()) else "NO_SALE",
    }
    out["verdict_U13"] = "M_G_DE_LA_ACCION" if any(v.startswith("SALE") for k, v in out["routes"].items() if k[0] in "123") else "M_G_REQUIERE_RUPTURA_NUEVA"
    print(json.dumps(out["routes"], ensure_ascii=False, indent=1)); print("U13:", out["verdict_U13"])
    (OUT / "U13_RESULT.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
