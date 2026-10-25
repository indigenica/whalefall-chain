"""S4 · dark ledger: misalignment of the pseudo-Goldstone (C35b, vendored unchanged) + collisionless packets emitted by the
string network with the asymptotic GHN normalisation times an efficiency eps (cadena_cierre/dark.py, re-written here in
vectorised form and checked packet-by-packet against the vendored original).  New: comoving free-streaming length of the
emitted component, and the same integrator applied to a 5.3 keV thermal WDM reference (Lyman-alpha, Irsic et al. 2017)."""
import math
from ucommon import load, save, benchmark, provenance, PLATEAU
import numpy as np
from scipy.optimize import brentq
import c35b_pseudo_goldstone_dm_rp3 as C35

MPL = C35.MPL; T0 = 2.7255 * 8.617333262e-14; H100 = 100 / 3.085677581e19 * 6.582119569e-25
RHOC_H2 = 3 * MPL ** 2 * H100 ** 2
_lT = np.log(C35._T[::-1]); _G = C35._G[::-1]; _GS = C35._GS[::-1]
# late-time background for free streaming (Planck 2018 base LCDM)
h = 0.674; OM = 0.315; OR = 4.18e-5 / h ** 2; OL = 1 - OM - OR; C_OVER_H0_MPC = 2997.92458 / h
T_NU0_EV = (4 / 11) ** (1 / 3) * 2.7255 * 8.617333262e-5


def gstar(T): return np.interp(np.log(np.clip(T, C35._T[-1], C35._T[0])), _lT, _G)
def gs(T): return np.interp(np.log(np.clip(T, C35._T[-1], C35._T[0])), _lT, _GS)
def hubble(T): return np.sqrt(np.pi ** 2 * gstar(T) / 90) * T * T / MPL
def scale(T): return T0 / T * (3.91 / gs(T)) ** (1 / 3)


def source(F, m_eV, Trh, x0=10, efficiency=1.0, step=0.03, spectrum="ir_delta", spectral_bins=64):
    m = m_eV * 1e-9; fa = F / 2; mr = F / 57.735
    Tstop = math.exp(brentq(lambda l: math.log(hubble(math.exp(l)) / (m / 3)), math.log(1e-12), math.log(Trh)))
    edge = np.linspace(math.log(Trh), math.log(Tstop), int(math.ceil(math.log(Trh / Tstop) / step)) + 1)
    hi, lo = edge[:-1], edge[1:]
    T = np.exp((hi + lo) / 2); H = hubble(T); a = scale(T)
    dN = hi - lo + np.log(gs(np.exp(hi)) / gs(np.exp(lo))) / 3
    ell = np.log(mr / H); xi = 0.24 * ell; rho_s = 4 * xi * np.pi * fa ** 2 * ell * H ** 2
    dE = efficiency * 2 * rho_s * dN
    if spectrum == "ir_delta":
        X = np.full((len(T), 1), float(x0)); W = np.ones((len(T), 1))
    elif spectrum == "energy_logflat":
        qd, wt = np.polynomial.legendre.leggauss(spectral_bins)
        lo_p = math.log(x0); hi_p = np.log(mr / H)[:, None]
        X = np.exp((lo_p + hi_p) / 2 + qd[None, :] * (hi_p - lo_p) / 2); W = np.broadcast_to(wt / 2, X.shape)
    else:
        raise ValueError(spectrum)
    p = X * H[:, None]; E = np.sqrt(p * p + m * m)
    nc = dE[:, None] * W / E * a[:, None] ** 3
    A = np.broadcast_to(a[:, None], X.shape)
    return dict(q=(p * A).ravel(), nc=nc.ravel(), emitted=(dE[:, None] * W).ravel(), a=A.ravel().copy(),
                T_stop_GeV=Tstop, m_GeV=m, F_GeV=F, efficiency=efficiency, x0=x0, spectrum=spectrum)


def abundance(pack): return pack["m_GeV"] * float(np.sum(pack["nc"])) / RHOC_H2


def at_epoch(pack, T):
    a = float(scale(T)); act = pack["a"] <= a; n = pack["nc"][act] / a ** 3
    p = pack["q"][act] / a; m = pack["m_GeV"]; E = np.sqrt(m * m + p * p)
    rho = float(np.sum(n * E)); P = float(np.sum(n * p * p / (3 * E)))
    rho_nu_one = 7 * np.pi ** 2 / 120 * T ** 4 * ((4 / 11) ** (4 / 3) if T < 1e-4 else 1.0)
    return dict(T_GeV=T, rho=rho, DeltaNeff_radiation_moment=3 * P / rho_nu_one)


def _E(a): return np.sqrt(OR / a ** 4 + OM / a ** 3 + OL)


def lambda_fs(q_eV, m_eV, a_e, n_grid=1500):
    """comoving free-streaming length [Mpc] from emission a_e to today: (c/H0) int v/(a E) dln a, v = q/sqrt(q^2+m^2 a^2).
    q is today's physical momentum (a0 = 1).  Vectorised over packets with a common log grid and per-packet lower cut."""
    q = np.atleast_1d(q_eV); ae = np.atleast_1d(a_e)
    lna = np.linspace(math.log(ae.min()), 0.0, n_grid); a = np.exp(lna)
    out = np.empty(len(q))
    for s in range(0, len(q), 2000):
        qq = q[s:s + 2000, None]; v = qq / np.sqrt(qq ** 2 + (m_eV * a[None, :]) ** 2)
        f = v / (a * _E(a))[None, :]
        cum = np.concatenate([np.zeros((len(qq), 1)), np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * np.diff(lna)[None, :], axis=1)], axis=1)
        tot = cum[:, -1]
        at = np.array([np.interp(math.log(x), lna, c) for x, c in zip(ae[s:s + 2000], cum)])
        out[s:s + 2000] = (tot - at) * C_OVER_H0_MPC
    return out


def lambda_wdm_reference(m_keV=5.3):
    m = m_keV * 1e3; Tx = (0.12 * 94.0 / m) ** (1 / 3) * T_NU0_EV
    return float(lambda_fs(3.15 * Tx, m, 1e-12, n_grid=4000)[0]), dict(m_eV=m, Tx_over_Tnu=Tx / T_NU0_EV, mean_p_eV=3.15 * Tx)


def emitted_lambda(pack):
    """mass-weighted (= number-weighted) mean free-streaming length of the emitted packets [Mpc]"""
    if pack["efficiency"] == 0 or not np.sum(pack["nc"]):
        return 0.0
    keep = pack["nc"] > 1e-12 * pack["nc"].max()
    lam = lambda_fs(pack["q"][keep] * 1e9, pack["m_GeV"] * 1e9, pack["a"][keep])
    return float(np.sum(pack["nc"][keep] * lam) / np.sum(pack["nc"][keep]))


def dm_root(F, Trh, eff, spectrum, x0, K_sw, step=0.03):
    def budget(lm):
        m = 10 ** lm; mis = C35.omega(m, F, K_sw)[0]
        return mis + (abundance(source(F, m, Trh, x0, eff, step, spectrum)) if eff else 0.0)
    return 10 ** brentq(lambda lm: math.log(budget(lm) / 0.12), -27, -10, xtol=1e-7)


def main():
    """reproduction gates: (1) vectorised source == vendored dark.source packet by packet; (2) DM roots == dark.json;
    (3) free-streaming integrator: relativistic limit = comoving horizon, cold limit -> 0.
    U26: (1) and (2) are code gates against archives of the sum, so they use the sum's benchmark; with the pole the same roots are
    recomputed with its F and T_Gamma and their shift reported."""
    B = benchmark("sum"); F = B["F_GeV"]; Trh = B["T_Gamma_GeV"]
    import dark as D  # vendored original
    out = dict(stage="S4_dark_ledger", provenance=provenance())
    worst = 0.0
    for spec in ("ir_delta", "energy_logflat"):
        for x0 in (10, 100):
            m = 3e-17
            a = D.source(F, m, Trh, x0, 0.1, spectrum=spec); b = source(F, m, Trh, x0, 0.1, spectrum=spec)
            for k in ("q", "nc", "a"):
                worst = max(worst, float(np.max(np.abs(np.asarray(a[k]) / b[k] - 1))))
    out["vectorised_vs_vendored_max_rel"] = worst
    arch = load("cadena_cierre/dark.json")["conditional_DM_roots"]
    dl = []
    for r in arch:
        m = dm_root(F, Trh, r["efficiency"], r["spectrum"], 10, 0)
        dl.append(abs(math.log10(m / r["m_G_eV"])))
    out["roots_vs_dark_json_max_dlog10"] = max(dl)
    if PLATEAU == "pole":
        Bp = benchmark("pole")
        dp = [math.log10(dm_root(Bp["F_GeV"], Bp["T_Gamma_GeV"], r["efficiency"], r["spectrum"], 10, 0) / r["m_G_eV"]) for r in arch]
        out["pole_roots_vs_sum_dlog10"] = dict(max_abs=max(abs(x) for x in dp), min=min(dp), max=max(dp), n=len(dp))
    # free-streaming checks
    lam_rel = float(lambda_fs(1.0, 0.0, 1e-12)[0])
    from scipy.integrate import quad   # independent integrator for the comoving horizon
    eta0 = C_OVER_H0_MPC * quad(lambda la: 1 / (math.exp(la) * _E(math.exp(la))), math.log(1e-12), 0, limit=400, epsrel=1e-10)[0]
    out["fs_relativistic_vs_horizon_rel"] = abs(lam_rel / eta0 - 1)
    out["fs_cold_Mpc"] = float(lambda_fs(1e-20, 1.0, 1e-3)[0])
    lw, info = lambda_wdm_reference(); out["lambda_WDM_5p3keV_Mpc"] = lw; out["WDM_reference"] = info
    out["comoving_horizon_today_Mpc"] = eta0
    out["status"] = "PASS" if (worst < 1e-10 and max(dl) < 1e-4 and out["fs_relativistic_vs_horizon_rel"] < 1e-3) else "FAIL"
    save("s4_dark_checks.json", out)
    print(f"S4 {out['status']}: vectorised/vendored {worst:.1e}, roots dlog10 {max(dl):.1e}, horizon {eta0:.0f} Mpc "
          f"(rel {out['fs_relativistic_vs_horizon_rel']:.1e}), lambda_WDM(5.3 keV) = {lw:.4f} Mpc")
    return out


if __name__ == "__main__":
    main()
