"""U59b (docs/preregistros/ESPEC_U59b.md): U59 with a new mother (adopted set, A.H = 300), a roundoff-aware G1, the energy gain
E = G^2 and a growth rule by decades on a dense scan; informative C.3 (n = 0, 1 against the explorer's C4) and C.5 (z^2, c_s^2 of
Boruah-Kim-Geshnizjani in the first 60 e-folds of inflation).  u59_escalares_rebote and u58_cuscuton_inflacion are imported
unchanged.
    python3 -B chain/u59b_escalares_rebote.py   ->  results/u59b_escalares_rebote.json
Internal, not exported."""
import json, math, time
import numpy as np
from ucommon import ROOT, save
import explora_cuscuton as EC
import u59_escalares_rebote as U
import u58_cuscuton_inflacion as U58

DENSE = tuple(int(x) for x in np.unique(np.round(np.logspace(2, math.log10(3e5), 28))))
C4_PUB = {0: dict(aux=[-0.004999999999930562, 0.0, 0.000467878824074332], kin=[0.0, 0.7265592702634625, 8.580000000202485]),
          1: dict(aux=[-6.012135567541188, 0.000467878824074332, 0.0071356162436229794], kin=[0.0, 0.0, 0.7265592702634625]),
          2: dict(aux=[-16.03234472115117, 0.000467878824074332, 0.027344851024213274], kin=[0.0, 0.7265592702634625, 2.6095240811535536])}


def c1_bounce_bg():
    ref = json.loads((ROOT / "data" / "explorador" / "cuscuton_c1" / "result.json").read_text())["bounce"]
    cu = EC.Cusc(**U.C1); _, Uv0, dU0, ddU0 = cu.all(0, 0.0); chi = ref["chi"]; Cp, V0 = EC.CP, EC.V0
    return dict(a=ref["a"], H=0.0, chid=ref["chi_dot"], phid=ref["phidot"], q=(4 * ref["rho_r"] / 3) ** 0.25, mu2=cu.mu2, Vv=EC.Vp(chi),
                Vp=EC.dVp(chi), Vpp=V0 * Cp * (Cp * Cp - 2 * chi * chi) / (Cp * Cp + chi * chi) ** 2.5, Uv=Uv0, Up=dU0, Upp=ddU0)


def exceptional(A, bg, n):
    """aux spectrum and full kinetic spectrum (q = zeta, x, w) for n = 0 (beta absent) and n >= 1"""
    m = A.matrices(bg, n)
    keep = [0, 2] if n == 0 else [0, 1, 2]
    S = m["S"]; T = m["T"]
    Ss = S[np.ix_(keep, keep)]; Ts = T[:, keep]
    K = m["K0"] - Ts @ np.linalg.solve(Ss, Ts.T); K = 0.5 * (K + K.T)
    sc = np.max(np.abs(K))
    return dict(aux_eig=sorted(float(e) for e in np.linalg.eigvalsh(S)), aux_rank=int(np.linalg.matrix_rank(S, tol=1e-10 * np.max(np.abs(S)))),
                kin_eig=sorted(float(e) if abs(e) > 1e-10 * sc else 0.0 for e in np.linalg.eigvalsh(K)))


def main():
    t0 = time.time(); A = U.Action(); out = dict(spec="docs/preregistros/ESPEC_U59b.md", dense_n=DENSE)
    bg1 = c1_bounce_bg()
    g0 = {}
    for n, (lo, hi) in U.PUB_C2.items():
        ev = np.linalg.eigvalsh(U.red2(A, bg1, n)[0])
        g0[n] = max([abs(ev[1] / hi - 1)] + ([abs(ev[0] / lo - 1)] if lo else []))
    G0 = all(v <= 1e-8 for v in g0.values())
    G1v = dict(n37=A.quad_check(37), n3e5=A.quad_check(300000)); G1 = all(v <= 1e-8 for v in G1v.values())
    print("G0", g0, "G1", G1v, flush=True)
    # the new mother
    rows = json.loads((ROOT / "results" / "u57b_runs" / "rows_rtol11.json").read_text())["rows"]
    k = sorted([r for r in rows if r["set"] == "eps_pequeno" and r["par"]["AH"] == 300.0 and r["outcome"] == "HIJO"], key=lambda r: r["par"]["f"])[0]
    cu, Crad, yb, o = U.host_state(dict(U.ADOPTED, AH=300.0, f=k["par"]["f"]))
    tr = U.Trajectory(cu, Crad, yb)
    kmin = math.inf; dmin = math.inf
    for ph in np.linspace(-0.2, 0.2, 81):
        t = tr.t_of_phi(ph) if abs(ph) > 1e-12 else 0.0; bg = tr.bg(t)
        for n in U.NS:
            K, _, _, dS = U.red2(A, bg, n); ev = np.linalg.eigvalsh(K)
            kmin = min(kmin, ev[0] / ev[1]); dmin = min(dmin, dS)
    tf = []
    for n in DENSE:
        r_ = U.transfer(A, tr, n, 1e-10); r_["E"] = r_["G"] ** 2 if r_["G"] else None; tf.append(r_)
        print(n, "%.3f" % r_["k_phys_bounce"], "E %.3f" % r_["E"], flush=True)
    # G2 at the three control n of ESPEC_U59b; a control n absent from the dense grid (10000) gets its own rtol 1e-10 reference
    # (code fix of 2026-11-17, after the Codex review of PR #89: the original comprehension silently dropped it)
    ref10 = {x["n"]: x["G"] for x in tf}
    conv = {}
    for n in (100, 10000, 300000):
        g10 = ref10[n] if n in ref10 else U.transfer(A, tr, n, 1e-10)["G"]
        conv[n] = abs(U.transfer(A, tr, n, 1e-9)["G"] / g10 - 1)
    E = {x["n"]: x["E"] for x in tf}
    dec = lambda lo, hi: max(v for n, v in E.items() if lo <= n < hi)
    M = {3: dec(1e3, 1e4), 4: dec(1e4, 1e5), 5: dec(1e5, 3.01e5), 2: dec(1e2, 1e3)}
    G2 = all(v <= 1e-4 for v in conv.values())
    G3 = all(x["symp"] <= 1e-6 for x in tf); G4 = all(x["energy_positive"] for x in tf)
    gates = dict(G0=G0, G1=G1, G2=G2, G3=G3, G4=G4)
    sane = bool(kmin > 0 and dmin > 1e-12 and M[4] <= 2 * M[3] and M[5] <= 2 * M[4])
    unstable = bool(kmin < -1e-9 or (M[3] < M[4] < M[5] and M[5] > 10 * M[3]))
    verdict = "INCONCLUSO" if not all(gates.values()) else ("ESCALARES_SANOS" if sane else ("INESTABILIDAD" if unstable else "INCONCLUSO"))
    # C.3 (informative): n = 0, 1, 2, 3 at the C1 bounce against C4, and at the adopted bounce
    c3 = {"C1": {n: exceptional(A, bg1, n) for n in (0, 1, 2, 3)}, "adoptado": {n: exceptional(A, tr.bg(0.0), n) for n in (0, 1, 2, 3)}}
    c3cmp = {}
    for n, pub in C4_PUB.items():
        mine = c3["C1"][n]
        rel = lambda a_, b_: max(abs(x - y) / max(abs(y), 1e-300) if abs(y) > 0 else abs(x) for x, y in zip(sorted(a_), sorted(b_)))
        c3cmp[n] = dict(aux_rel=rel(mine["aux_eig"], pub["aux"]), kin_rel=rel(mine["kin_eig"], pub["kin"]))
    # C.5 (informative): BKG z^2, c_s^2 in the first 60 e-folds of inflation, flat background, adopted set
    model = U58.Model(U.ADOPTED); s, Nend = U58.background(model, 1e-10); bgf = U58.make_bg(model, s, Nend)
    zmin, cmin = math.inf, math.inf
    for N in np.linspace(0.5, min(60.0, Nend - 1), 200):
        b = bgf(N)
        for kr in np.logspace(-3, 3, 25):
            z2, c2 = U58.z2_cs2(b, kr * b["aH"], N, False)
            zmin = min(zmin, z2 / (math.exp(2 * N) * b["alpha"])); cmin = min(cmin, c2)
    out.update(gates=gates, G0_rel=g0, G1=G1v, mother=dict(AH=300.0, f=k["par"]["f"], N=k["N_inf"], a_b=math.exp(yb[0])),
               min_K_ratio=kmin, min_detS_rel=dmin, transfers=tf, E=E, decade_max=M, conv_rtol=conv,
               would_be=("ESCALARES_SANOS" if sane else ("INESTABILIDAD" if unstable else "INCONCLUSO")), verdict=verdict,
               C3=c3, C3_vs_C4=c3cmp, C5=dict(min_z2_over_a2alpha=zmin, min_cs2=cmin, N_range=[0.5, min(60.0, Nend - 1)]),
               seconds=time.time() - t0)
    save("u59b_escalares_rebote.json", out)
    print(json.dumps({k_: out[k_] for k_ in ("gates", "G0_rel", "G1", "mother", "min_K_ratio", "min_detS_rel", "decade_max", "conv_rtol",
                                              "would_be", "verdict", "C3_vs_C4", "C5")}, indent=1, default=float))
    print("C3 adoptado", json.dumps(c3["adoptado"], default=float))


if __name__ == "__main__":
    main()
