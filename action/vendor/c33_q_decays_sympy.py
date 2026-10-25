#!/usr/bin/env python3
"""C33 — Do ALL massive modes of Q decay (to Goldstone pairs or lighter heavy pairs) fast enough?  sympy, tree level.

Q = R(theta) (P + h) R(theta)^T,  P = e1 e1^T,  R = exp(sum_j theta_j A_j),  A_j = e1 e_j^T - e_j e1^T (j = 2,3,4: Goldstones).
Gauge: h_1j = 0 (j >= 2), so h = a e11 + B (B symmetric in the 2-3-4 block): 7 heavy fields.
L = (ZQ/2) tr(dQ dQ) - U(Q),  U = LAM/4 (tr Q - 1)^2 + LP/4 ||Q^2 - Q||^2  (parent, p4_hermite at s = 0).
 - Goldstone-heavy couplings come only from the kinetic metric at theta = 0, linear in h (derivative couplings);
 - heavy masses and heavy cubic couplings from U(P + h) to cubic order.
Canonical: theta = pi / F, F^2 = 2 ZQ;  a = x/sqrt(ZQ);  B_jj = y_j/sqrt(ZQ);  B_jk (j<k) = z_jk/sqrt(2 ZQ).
Rates (massless Goldstones): L3 = phi * dpi^T G dpi  ->  Gamma = m^3/(16 pi) [sum_{j<k} G_jk^2 + 1/2 sum_j G_jj^2].
All quantities in code units (energy unit E0); physical translation: m_Q = m_lightest_heavy * E0 (feeds C34: F = sqrt(2 ZQ) E0).
Frozen labels in C33_C34_PREREGISTRATION.md.
"""
import json, math, sys
from pathlib import Path
import numpy as np
import sympy as sp

ZQ, LP, LAM = 400, 192, 320

def build():
    th = sp.symbols('t2:5'); dth = sp.symbols('dt2:5')
    a = sp.Symbol('a'); da = sp.Symbol('da')
    bs = {}; dbs = {}
    for j in range(1, 4):
        for k in range(j, 4):
            bs[(j, k)] = sp.Symbol(f'b{j+1}{k+1}'); dbs[(j, k)] = sp.Symbol(f'db{j+1}{k+1}')
    def mat(sym_a, sym_b):
        M = sp.zeros(4, 4); M[0, 0] = sym_a
        for (j, k), s in sym_b.items(): M[j, k] = s; M[k, j] = s
        return M
    P = sp.zeros(4, 4); P[0, 0] = 1
    h = mat(a, bs); dh = mat(da, dbs)
    A = []
    for j in range(1, 4):
        Aj = sp.zeros(4, 4); Aj[0, j] = 1; Aj[j, 0] = -1; A.append(Aj)
    # kinetic metric at theta = 0: dQ = sum dth_j [A_j, P + h] + dh
    dQ = dh + sum((dth[i] * (A[i] * (P + h) - (P + h) * A[i]) for i in range(3)), sp.zeros(4, 4))
    kin = sp.expand((dQ * dQ).trace())
    # Goldstone-Goldstone-heavy part: coefficient of dth_i dth_j, linear in h
    Gsym = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            c = sp.Poly(kin, *dth).coeff_monomial(dth[i] * dth[j]) if i != j else sp.Poly(kin, *dth).coeff_monomial(dth[i] ** 2)
            Gsym[i, j] = sp.expand(c if i == j else c / 2)
    # potential to cubic order
    Q0 = P + h
    U = sp.Rational(LAM, 4) * (Q0.trace() - 1) ** 2 + sp.Rational(LP, 4) * sum(e ** 2 for e in (Q0 * Q0 - Q0))
    U = sp.expand(U)
    fields = [a] + [bs[k] for k in sorted(bs)]
    return dict(G=Gsym, U=U, fields=fields, bkeys=sorted(bs), kin=kin, dth=dth)


def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    S = build(); fields = S['fields']
    # canonical normalisation factors: field = c * canonical
    cfac = [1 / math.sqrt(ZQ)] + [(1 / math.sqrt(ZQ) if j == k else 1 / math.sqrt(2 * ZQ)) for (j, k) in S['bkeys']]
    zero = {f: 0 for f in fields}
    # quadratic potential -> mass matrix in canonical fields: U2 = 1/2 m2_ij phi_i phi_j
    Hs = np.array([[float(sp.diff(S['U'], fi, fj).subs(zero)) for fj in fields] for fi in fields])
    C = np.diag(cfac); M2 = C @ Hs @ C
    assert abs(float(sp.diff(S['U'], fields[0]).subs(zero))) < 1e-12   # vacuum
    w, V = np.linalg.eigh(M2)
    masses = np.sqrt(np.clip(w, 0, None))
    # Goldstone couplings: kinetic term (ZQ/2) tr(dQ^2) -> (ZQ/2) * G_ij(h) dth_i dth_j ; theta = pi/F, F^2 = 2 ZQ
    # L3 = (ZQ/2)/F^2 * sum_ij G_ij^(1)(h) dpi_i dpi_j = (1/4) G^(1)(h) ... linear part in h
    G1 = [[sp.expand(S['G'][i, j] - S['G'][i, j].subs(zero)) for j in range(3)] for i in range(3)]
    G0 = np.array([[float(S['G'][i, j].subs(zero)) for j in range(3)] for i in range(3)])
    F2 = 2 * ZQ
    kin_norm_check = (ZQ / 2) * G0 / F2   # should be 1/2 * identity (canonical Goldstones)
    rows = []
    for n in range(len(fields)):
        vec = V[:, n]                                    # canonical eigenmode: phi_i = vec_i * s
        subs = {f: cfac[i] * vec[i] for i, f in enumerate(fields)}
        Gm = np.array([[float(sp.expand(G1[i][j]).subs(subs)) if G1[i][j] != 0 else 0.0 for j in range(3)] for i in range(3)])
        Gm = (ZQ / 2) * Gm / F2                          # L3 = s * dpi^T Gm dpi
        m = masses[n]
        # L3 = s * sum_jk Gm_jk dpi_j.dpi_k ; |M(j,k)| = Gm_jk m^2 ; Gamma = m^3/(16 pi) [sum_{j<k} Gm_jk^2 + 1/2 sum_j Gm_jj^2]
        Gam_GG = m ** 3 / (16 * math.pi) * (sum(Gm[j, k] ** 2 for j in range(3) for k in range(j + 1, 3)) + 0.5 * sum(Gm[j, j] ** 2 for j in range(3)))
        rows.append(dict(mode=n, m=m, Gamma_GG=Gam_GG, Gamma_GG_over_m=Gam_GG / m if m > 0 else float('nan'),
                         vec=[float(x) for x in vec]))
    # heavy cubic couplings for modes that cannot go to GG: decays into pairs of lighter heavy modes
    cub = sp.Poly(S['U'], *fields)
    def U3_in_modes():
        # substitute fields = sum_n cfac_i V_in s_n and take the cubic part
        s = sp.symbols(f's0:{len(fields)}')
        sub = {f: sum(cfac[i] * V[i, n] * s[n] for n in range(len(fields))) for i, f in enumerate(fields)}
        U3 = sum(c * sp.prod([f ** e for f, e in zip(fields, mon)]) for mon, c in cub.terms() if sum(mon) == 3)
        return sp.Poly(sp.expand(U3.subs(sub)), *s), s
    P3, s = U3_in_modes()
    for r in rows:
        n = r['mode']; mT = r['m']; tot = 0.0; ch = []
        for i in range(len(fields)):
            for j in range(i, len(fields)):
                if masses[i] + masses[j] >= mT - 1e-9: continue
                mon = [0] * len(fields); mon[n] += 1; mon[i] += 1; mon[j] += 1
                kappa = float(P3.coeff_monomial(tuple(mon)))
                if abs(kappa) < 1e-14: continue
                # U3 = kappa * s_n s_i s_j -> vertex = kappa * (symmetry factor): count permutations
                from collections import Counter
                cnt = Counter([n, i, j]); Mamp = kappa * float(np.prod([math.factorial(v) for v in cnt.values()]))   # -d3U/ds_n ds_i ds_j
                pcm = math.sqrt(max((mT ** 2 - (masses[i] + masses[j]) ** 2) * (mT ** 2 - (masses[i] - masses[j]) ** 2), 0)) / (2 * mT)
                g = Mamp ** 2 * pcm / (8 * math.pi * mT ** 2) * (0.5 if i == j else 1.0)
                tot += g; ch.append((i, j, g))
        r['Gamma_heavy'] = tot; r['Gamma_heavy_over_m'] = tot / mT if mT > 0 else float('nan'); r['n_heavy_channels'] = len(ch)
        r['has_channel'] = bool(r['Gamma_GG'] > 1e-14 * mT or tot > 1e-14 * mT)
    m_light = float(min(masses))
    res = dict(masses_code=[float(x) for x in masses], kinetic_norm_check=kin_norm_check.tolist(), modes=rows,
               m_lightest=m_light, F_over_mlightest=math.sqrt(2 * ZQ) / m_light,
               analytic_expectation=dict(m_5plet=math.sqrt(LP / (2 * ZQ)), m_trace=math.sqrt(2 * (LAM + LP / 4) / ZQ),
                                         Gamma5_over_m=(LP / (2 * ZQ)) / (32 * math.pi * ZQ)))
    # physical comparison: m_Q (lightest heavy) in units of the scalaron mass M; H at the end of reheating from C37
    try:
        c37 = json.loads((out.parent / 'c37_reheating_branching' / 'C37_RESULTS.json').read_text())
        pick = min(c37['scan'], key=lambda r_: abs(r_['xi_H'] - 1.1)); M = pick['M_GeV']
        H_rh = pick['Gamma_over_Gamma0'] * M ** 3 / (192 * math.pi * 2.435e18 ** 2)
    except Exception:
        M = 3.17e13; H_rh = 103 * M ** 3 / (192 * math.pi * 2.435e18 ** 2)
    phys = []
    for mu in (0.5, 0.7, 1.0, 2.0, 3.0):
        mQ = mu * M; E0 = mQ / m_light
        slowest = min((max(r['Gamma_GG'], r['Gamma_heavy']) for r in rows), default=0.0) * E0
        phys.append(dict(mu=mu, m_Q_GeV=mQ, slowest_Gamma_GeV=slowest, H_rh_GeV=H_rh, ratio=slowest / H_rh))
    res['physical'] = phys
    allch = all(r['has_channel'] for r in rows)
    fast = all(p['ratio'] > 1 for p in phys)
    res['verdict'] = ('C33_ALL_HEAVY_DECAY' if (allch and fast) else ('C33_SLOW_DECAY' if allch else 'C33_STABLE_RELIC'))
    (out / 'C33_RESULTS.json').write_text(json.dumps(res, indent=2, default=float))
    print('masses (code units):', np.round(masses, 4))
    print('kinetic norm check (should be 0.5*I):', np.round(kin_norm_check, 6).tolist())
    for r in rows:
        print(f"mode {r['mode']}: m={r['m']:.4f} Gamma_GG/m={r['Gamma_GG_over_m']:.3e} Gamma_heavy/m={r['Gamma_heavy_over_m']:.3e} channels={r['n_heavy_channels']}")
    print('analytic expectation:', res['analytic_expectation'])
    print('F/m_lightest =', round(res['F_over_mlightest'], 3))
    for p in phys: print(p)
    print('VERDICT', res['verdict'])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'claude/results/c33_q_decays')
