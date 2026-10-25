#!/usr/bin/env python3
"""C46 + C47 — Ballena 4 (chi hierarchy inflaton): can chi reheat the universe, and does its Z2 leave a lethal relic?  (CPU, < 1 min)

Model: V = V0/Z sum_{n<10} 2^-n tanh^2(chi/(0.5 2^n))  (C44E; C44i attractor).  M_Pl (reduced) = 1 unless stated.  No scalaron (the R^2
term is dropped in ballena 4: inflation is chi, the bounce is the QT tower of C45).
A. CMB normalization, self-consistent with N_* (formula of C37: matter-like reheating, w = 0 declared).
B. C47 — EXACT Z2 (only Z2-even couplings, (g^2/2) chi^2 |H|^2): the chi condensate/quanta can only annihilate (chi chi -> h h), never decay.
   Boltzmann for rho_chi (matter) and rho_r, <sigma v> = g^4/(16 pi m^2) (4 real Higgs dof, s-wave, declared), plus an arm with the
   s-wave unitarity cap <sigma v> = 16 pi / m^2 (most optimistic depletion, independent of the coupling and of preheating details:
   any 2 -> 2 remnant is bounded by it).  Omega h^2 = 2.744e8 Y (m / GeV).
   Label C47_EXACT_Z2_FATAL if even the unitarity arm gives Omega_chi h^2 > 1.2 (10 x the observed DM).
C. C46 — SOFT Z2 breaking (the minimal physically sensible replacement if B is fatal): a trilinear sigma chi |H|^2 (the ONLY
   renormalizable operator that lets chi decay to the visible sector; technically natural: sigma -> 0 restores the Z2).
   Gamma = sigma^2 / (8 pi m_chi).  Gates (declared before running):
     G_EW   : sigma chi_end <= m_chi^2  (no tachyonic Higgs resonance while chi oscillates: EW-vacuum fragility, Ema+ 2016-type; conservative)
     G_RAD  : sigma^2/(16 pi^2) <= m_chi^2  (chi mass radiatively stable)
     G_PERT : Gamma < H_end
     G_BBN  : T_rh > 4 MeV
     G_DNEFF: DeltaN_eff(reheating) <= 0.076 (ACT-inspired gate of C37), with
              (i) gravitational/xi channel chi chi -> Q Q through xi Q^2 R, R = (4V - chidot^2)/M_Pl^2 (amplitude 6 xi m^2/M_Pl^2 per real
                  component, N_Q = 10, declared +-factor 3; all of it counted as dark radiation = upper estimate),
              (ii) thermal Goldstones through a Q-Higgs portal: equilibrate only if T_rh > T_G = (F^4/M_Pl)^(1/3) (F = 1e15 GeV, C34),
                   then 3 x 0.027; otherwise freeze-in 0.081 (T_rh/T_G)^3.
   Reported (prediction, not a gate): n_s(T_rh) and T_rh needed for n_s >= ACT - 1 sigma (0.9709).
   Labels: C46_SOFT_Z2_VIABLE (some sigma passes all gates), C46_NO_WINDOW.  Secondary: F1_GONE (DeltaN_eff(reheating) < 0.01 over the
   whole viable window, i.e. the ballena-3 F1 problem and its xi_H patch disappear without a new patch), NS_BAND (n_s over the window).
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from c44i_hierarchy_universality import predict, model

MPL = 2.435e18; AS = 2.1e-9
H0_over_c = 0.674 * 100 / 299792.458; KPIV = 0.05; T0 = 2.7255 * 8.617e-14; GS0 = 3.91; GSTAR = 106.75
Q_, S_, F_, L_ = 0.5, 2.0, 0.5, 10
DNEFF_GATE = 0.076; F_Q = 1e15


def chi_numbers():
    w = Q_ ** np.arange(L_); Z = w.sum(); sc = F_ * S_ ** np.arange(L_)
    m2_over_V0 = 2 * float((w / sc ** 2).sum()) / Z
    return m2_over_V0


def normalize(Trh_fun, N0=52.0):
    """Self-consistent N_*: Trh_fun(V0, m, H_end, rho_end) -> (T_rh [GeV], rho_rh [M_Pl^4])."""
    V, dV, _ = model(Q_, S_, F_, L_)
    m2V0 = chi_numbers(); N = N0
    for _ in range(80):
        p = predict(Q_, S_, F_, L_, N)
        V0 = 1.5 * math.pi ** 2 * p['r'] * AS / p['V_star_over_V0']
        m = math.sqrt(m2V0 * V0)
        rho_end = 1.5 * V0 * float(V(np.array([p['chi_end']]))[0]); H_end = math.sqrt(rho_end / 3)
        Hstar = math.sqrt(V0 * p['V_star_over_V0'] / 3)
        Trh, rho_rh = Trh_fun(V0, m, H_end, rho_end)
        H0 = H0_over_c * 299792.458 / 3.0857e19 * 6.582e-25 / MPL
        Nn = (-math.log(KPIV / H0_over_c) + math.log(Hstar / H0) + (1 / 3) * math.log(rho_rh / rho_end)
              + (1 / 3) * math.log(GS0 / GSTAR) + math.log(T0 / Trh))
        if abs(Nn - N) < 1e-7: break
        N = Nn
    return dict(N_star=N, n_s=p['n_s'], r=p['r'], V0=V0, m_chi=m, m_chi_GeV=m * MPL, H_end=H_end, H_end_GeV=H_end * MPL,
                H_star_GeV=Hstar * MPL, chi_end=p['chi_end'], rho_end=rho_end, T_rh_GeV=Trh)


def T_of_rho(rho):  # rho in M_Pl^4 -> T in GeV
    return (30 * rho / (math.pi ** 2 * GSTAR)) ** 0.25 * MPL


def relic_exact_z2(m, H_end, sv):
    """Boltzmann in N = ln a from the end of inflation; all energy starts in chi (matter-like).  Returns Omega_chi h^2 and Y."""
    rho0 = 3 * H_end ** 2
    def rhs(N, y):
        lc, lr = y; rc, rr = math.exp(lc), math.exp(lr); H = math.sqrt((rc + rr) / 3)
        ann = sv * rc * rc / m / H                       # d rho / dN transferred
        return [-3 - ann / rc, -4 + ann / rr]
    s = solve_ivp(rhs, (0, 60), [math.log(rho0), math.log(rho0 * 1e-30)], method='Radau', rtol=1e-8, atol=1e-10)
    rc, rr = math.exp(s.y[0, -1]), math.exp(s.y[1, -1])
    T = T_of_rho(rr) / MPL
    sden = 2 * math.pi ** 2 / 45 * GSTAR * T ** 3
    Y = (rc / m) / sden
    frac = rc / (rc + rr)
    om = 2.744e8 * Y * m * MPL
    return dict(Omega_h2=om, Y=Y, rho_chi_over_total_end=frac, N_integrated=float(s.t[-1]))


def main(root):
    out = Path(root) / 'c46_chi_reheating_relic'; out.mkdir(parents=True, exist_ok=True)
    res = {}
    # ---------- B. exact Z2 ----------
    base = normalize(lambda V0, m, He, re: (T_of_rho(re), re))        # instant-reheating reference (only for m, H_end)
    m, He = base['m_chi'], base['H_end']
    res['chi_sector'] = base
    ex = []
    for g in (1e-3, 1e-2, 0.1, 0.3, 1.0, math.sqrt(4 * math.pi)):
        sv = g ** 4 / (16 * math.pi * m * m); r = relic_exact_z2(m, He, sv); r.update(g=g, sigma_v=sv); ex.append(r)
    uni = relic_exact_z2(m, He, 16 * math.pi / (m * m)); uni.update(arm='unitarity_cap')
    res['exact_z2'] = dict(scan=ex, unitarity=uni)
    res['C47_EXACT_Z2_FATAL'] = bool(uni['Omega_h2'] > 1.2 and all(r['Omega_h2'] > 1.2 for r in ex))
    # ---------- C. soft Z2 ----------
    TG = (F_Q ** 4 / MPL) ** (1 / 3)
    rows = []
    for lsig in np.arange(-24, -6.9, 0.25):
        sig = 10 ** lsig
        Gam = sig ** 2 / (8 * math.pi * m)
        def Trf(V0, mm, H_end, rho_end, Gam=Gam):
            rho_rh = min(3 * Gam ** 2, rho_end); return T_of_rho(rho_rh), rho_rh
        nr = normalize(Trf)
        mm, He_, ce = nr['m_chi'], nr['H_end'], nr['chi_end']
        # xi channel: fraction of chi energy converted into Q quanta before reheating (upper estimate, all -> dark radiation)
        xi = 0.25; NQ = 10
        Mamp = 6 * xi * mm ** 2
        sv_xi = NQ * Mamp ** 2 / (32 * math.pi * mm ** 2) / 2
        # integrate d(f)/dN with rho_chi ~ matter until H = Gamma: rho_DR/rho_chi grows as a (redshift) -> f_DR(rh) = int sv n / H * (a/a_rh)^(-1)...
        Ns = np.linspace(0, max(1e-3, (2 / 3) * math.log(He_ / max(Gam, 1e-300))), 2000)
        H_N = He_ * np.exp(-1.5 * Ns); rho_chi = 3 * H_N ** 2
        dfrac = sv_xi * (rho_chi / mm) / H_N                      # d rho_DR / rho_chi per e-fold at injection
        # rho_DR/rho_chi at reheating: injected at N scales as exp(-(N_rh - N)) relative to matter
        Nrh = Ns[-1]; f_xi = float(np.trapezoid(dfrac * np.exp(-(Nrh - Ns)), Ns)) if Nrh > 0 else 0.0
        dn_xi = 2.857 * f_xi
        Trh = nr['T_rh_GeV']
        dn_th = 0.081 if Trh > TG else 0.081 * (Trh / TG) ** 3
        dn = dn_xi + dn_th
        gates = dict(G_EW=bool(sig * ce <= mm ** 2), G_RAD=bool(sig ** 2 / (16 * math.pi ** 2) <= mm ** 2), G_PERT=bool(Gam < He_),
                     G_BBN=bool(Trh > 4e-3), G_DNEFF=bool(dn <= DNEFF_GATE))
        row = dict(sigma=sig, sigma_GeV=sig * MPL, sigma_over_m=sig / mm, Gamma=Gam, T_rh_GeV=Trh, N_star=nr['N_star'], n_s=nr['n_s'],
                   r=nr['r'], DeltaNeff_xi_upper=dn_xi, DeltaNeff_thermalG=dn_th, DeltaNeff=dn, **gates, pass_all=bool(all(gates.values())))
        rows.append(row)
    res['soft_z2'] = rows
    ok = [r for r in rows if r['pass_all']]
    res['T_G_GeV'] = TG
    res['C46_verdict'] = 'C46_SOFT_Z2_VIABLE' if ok else 'C46_NO_WINDOW'
    if ok:
        res['window'] = dict(sigma_GeV=[ok[0]['sigma_GeV'], ok[-1]['sigma_GeV']], sigma_over_m=[ok[0]['sigma_over_m'], ok[-1]['sigma_over_m']],
                             T_rh_GeV=[ok[0]['T_rh_GeV'], ok[-1]['T_rh_GeV']], N_star=[ok[0]['N_star'], ok[-1]['N_star']],
                             n_s=[min(r['n_s'] for r in ok), max(r['n_s'] for r in ok)], r=[min(r['r'] for r in ok), max(r['r'] for r in ok)],
                             DeltaNeff_max=max(r['DeltaNeff'] for r in ok))
        ns_act_1sig = 0.9743 - 0.0034
        hit = [r for r in ok if r['n_s'] >= ns_act_1sig]
        res['T_rh_min_for_ACT_1sigma_GeV'] = (hit[0]['T_rh_GeV'] if hit else None)
        res['F1_GONE'] = bool(max(r['DeltaNeff'] for r in ok) < 0.01)
        res['binding_gate_upper'] = next((k for k in ('G_EW', 'G_RAD', 'G_PERT', 'G_DNEFF') if not rows[rows.index(ok[-1]) + 1][k]), None) if rows.index(ok[-1]) + 1 < len(rows) else None
    (out / 'C46_RESULTS.json').write_text(json.dumps(res, indent=2, default=float))
    b = base; print(f"chi: V0={b['V0']:.3e} m={b['m_chi_GeV']:.3e} GeV H*={b['H_star_GeV']:.3e} H_end={b['H_end_GeV']:.3e} GeV chi_end={b['chi_end']:.3f}")
    for r in ex: print(f"exactZ2 g={r['g']:.3g}  Omega h2={r['Omega_h2']:.3e}")
    print(f"exactZ2 unitarity Omega h2={uni['Omega_h2']:.3e}")
    for r in rows[::4]: print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k in ('sigma_GeV', 'T_rh_GeV', 'N_star', 'n_s', 'DeltaNeff', 'pass_all', 'G_EW', 'G_DNEFF')})
    print({k: res.get(k) for k in ('C47_EXACT_Z2_FATAL', 'C46_verdict', 'window', 'F1_GONE', 'binding_gate_upper', 'T_rh_min_for_ACT_1sigma_GeV')})


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'claude/results')
