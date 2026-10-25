#!/usr/bin/env python3
"""C35b — RP^3 weld of C35 (after the external review of 27-sep): three corrections, same method.
 (1) isotropic misalignment on RP^3 = S^3/Z2: <n1^2> = 1/4 -> <sin^2 theta> = 3/4 -> rho_i = (3/8) m_G^2 F^2 (C35 used 1/4, a U(1) average);
 (2) F window with the corrected CMB convention: F in [57.7 m_Q(M/2), 1.26e15 GeV];
 (3) smooth g*(T), g_s(T) (log-interpolated table, declared approximate) instead of a step at 1 MeV.
K_sw (string+wall contribution) is still NOT calibrated for RP^3: 0, 1, 10 kept as a bracket.  CPU, seconds.
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

MPL = 2.435e18; OMEGA_DM = 0.12
_T = np.array([1e3, 1e2, 3e1, 1e1, 1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.01, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 1e-5])
_G = np.array([106.75, 106.75, 95.0, 86.25, 75.75, 70.0, 61.0, 40.0, 20.0, 17.25, 12.0, 10.75, 10.75, 10.6, 9.0, 5.5, 3.8, 3.36])
_GS = np.array([106.75, 106.75, 95.0, 86.25, 75.75, 70.0, 61.0, 40.0, 20.0, 17.25, 12.0, 10.75, 10.75, 10.6, 9.3, 6.0, 4.1, 3.91])
def _interp(T, arr):
    lt = np.log(np.clip(T, _T[-1], _T[0])); return float(np.interp(lt, np.log(_T[::-1]), arr[::-1]))
gstar = lambda T: _interp(T, _G); gs = lambda T: _interp(T, _GS)
C_MIS = 3.0 / 8.0


def omega(mG_eV, F, K_sw):
    m = mG_eV * 1e-9
    Tosc = brentq(lambda T: 3 * math.sqrt(math.pi ** 2 * gstar(T) / 90) * T * T / MPL - m, 1e-12, 1e8)
    rho = C_MIS * m * m * F * F
    s = 2 * math.pi ** 2 / 45 * gs(Tosc) * Tosc ** 3
    return 2.744e8 * m * (rho / m) / s * (1 + K_sw), Tosc


def m_needed(F, K_sw): return 10 ** brentq(lambda lg: math.log(omega(10 ** lg, F, K_sw)[0] / OMEGA_DM), -30, 5)


def minicluster_Msun(m_eV):   # Fairbairn et al. 1707.03310 characteristic mass (literature scaling)
    return 2.3e-7 * (m_eV / 1e-10) ** (-1.5)


def run(F_list):
    rows = []
    for F in F_list:
        for K in (0, 1, 10):
            m = m_needed(F, K); _, T = omega(m, F, K)
            rows.append(dict(F_GeV=F, K_sw=K, m_G_eV=m, T_osc_MeV=T * 1e3, M_mc_Msun=minicluster_Msun(m)))
    return rows


def main(root):
    root = Path(root); out = root / 'c35b_pseudo_goldstone_dm_rp3'; out.mkdir(parents=True, exist_ok=True)
    c34 = json.loads((root / 'c34_scale_window' / 'C34_RESULTS.json').read_text())
    lo = c34['F_over_mQ_toy'] * c34['window_primary'][0]; hi = 2 * 6.3e14
    rows = run([lo, (lo * hi) ** 0.5, hi])
    ms = [r['m_G_eV'] for r in rows]; Ts = [r['T_osc_MeV'] for r in rows]
    res = dict(F_window=[lo, hi], C_mis=C_MIS, rows=rows, m_G_range_eV=[min(ms), max(ms)], T_osc_range_MeV=[min(Ts), max(Ts)])
    res['verdict'] = 'C35B_DM_POSSIBLE' if all(t > 0.01 for t in Ts) else 'C35B_DM_IMPOSSIBLE'
    (out / 'C35B_RESULTS.json').write_text(json.dumps(res, indent=2, default=float))
    for r in rows: print(f"F={r['F_GeV']:.2e} K_sw={r['K_sw']:>2} m_G={r['m_G_eV']:.2e} eV T_osc={r['T_osc_MeV']:.2f} MeV M_mc~{r['M_mc_Msun']:.1e} Msun")
    print('VERDICT', res['verdict'], res['m_G_range_eV'], res['T_osc_range_MeV'])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'claude/results')
