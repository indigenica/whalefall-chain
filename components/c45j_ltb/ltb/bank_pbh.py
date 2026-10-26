"""C45b-PBH (decisive test): can a PHYSICALLY MOTIVATED overdensity of a dust-dominated parent
(primordial-black-hole-type curvature peak, growing mode) rebound coherently over >= 1/H_inf in QT?

python -m ltb.bank_pbh [--quick]
"""
import sys, json, time, math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from common.provenance import new_run_dir, write_json, write_manifest, source_hashes, environment
from common.potential import H_INF
from ltb.growing import GrowingShell
from ltb.shells import psi_of_y

WINDOW = (0.42, 9.6)
EPS = 1e-4
PROTOCOL = {
    "question": "Decisive test of the surviving niche of ballena 4 (dust parent + chi seed, C45 architecture): does any "
                "physically motivated PBH-type overdensity, in the growing mode of an Einstein-de Sitter dust parent, "
                "rebound coherently (no shell crossing) over a region of areal radius >= 1/H_inf, for a central "
                "turnaround scale A_c H_inf inside the C54 fertility window [0.42, 9.6]?",
    "model": "growing-mode LTB-QT (ltb/growing.py): M = mu x^3, -2E(x) = K(x) x^2 = 1 - (1 + r zeta')^2 with x = r e^zeta; "
             "exact QT shell equation; shells synchronized on their GR cycloids at t_i (alpha y_b = 0.1, capped at 0.3 of the "
             "earliest GR turnaround; the pre-t_i QT correction is <= 0.5 %; changed from 1e-3 after the smoke run because "
             "1e-3 put t_i after turnaround for A H ~ 1: a bug fix, declared)",
    "profiles": {
        "musco": "K(x) = K0 exp(-(x/x_m)^(2 a)/a) (Musco 2019 family), shape a in {1, 2.5, 6.33 (= monochromatic-equivalent), 15 (flatter than any Gaussian-peak profile, control)}",
        "peak": "zeta(r) = mu psi(r), psi = mean high-peak profile = normalized correlation function of zeta "
                "(BBKS; Yoo et al. 2018; Germani & Musco 2019): monochromatic (sinc) and lognormal spectra with width "
                "sigma in {0.5, 1, 2}",
        "control": "flat curvature: K(x) = K0 / (1 + exp((x - 1)/0.01)) (homogeneous closed core = the C45 patch)",
    },
    "amplitude": "max over x of -2E = K x^2 in {0.1, 0.2, ..., 0.9, 0.95, 0.99, 0.995} (type I; 1 is the type-II limit). "
                 "Refined from {0.3, 0.6, 0.9, 0.99} after the smoke run showed that the nonlinear map x = r e^zeta "
                 "inverts the bounce ordering near type II (outer shells first), so a flattest amplitude may exist in "
                 "between: a change that can only help PASS (declared)",
    "scale": "A_c H_inf in {0.42, 1, 2, 5, 9.6}, A_c = central turnaround scale (QT exact, reported)",
    "alpha": [0.1, 0.03],
    "shells": "600 shells uniform in x in [0, 1.2 x_m], x_m = argmax(-2E) (peak of the compaction); only bound shells (-2E > 0). "
              "Changed from uniform-in-mass after the convergence check: mass-uniform sampling puts only 1-2 shells inside "
              "x < 0.15 x_m and could not resolve small coherent cores (declared). Percentiles are mass weighted (x^2 dx).",
    "measures": "bounce-time spread over x <= x_m (P90-P10, times H_inf); first crossing; R_coh H_inf exactly as in C45b-LTB "
                "(largest inner region all bounced, all still expanding, no crossing among shells 0..k+1); ceiling R_ta(x_m) H_inf",
    "AMENDMENT_04oct": "After the independent audit found that outer shells fall through the re-expanding core "
                       "(not checked by the first criterion), the measure is recomputed with ltb/coherence.py: strict "
                       "(no intrusion) and tolerant (intruding outer mass <= f M_k, f = 0.01, 0.1, 0.3). Strict, 1 % and "
                       "10 % also kill the flat-curvature control (its soft edge falls onto the core within the very short "
                       "QT bounce time), so the verdict uses the LOOSEST tolerance f = 0.3, at which the control passes: a "
                       "FAIL there is the strongest possible FAIL. Thresholds 1 and 1/3 unchanged.",
    "verdict": {
        "PASS": "some PBH-type profile (musco a <= 6.33 or peak) reaches R_coh H_inf >= 1 for BOTH alphas at some amplitude and "
                "some A_c H_inf inside the window -> the niche survives, with the minimal amplitude as a condition",
        "FAIL": "every PBH-type configuration has R_coh H_inf < 1/3 -> no physically motivated progenitor rebounds coherently",
        "INCONCLUSIVE": "otherwise",
        "controls": "musco a = 15 and the flat-curvature control are reported but do not enter the verdict",
        "note_region": "R_coh requires a region that contains the centre; an annulus that collapses first is not a candidate patch",
    },
    "crossing_criterion": "AMENDED after the independent audit (audit_pbh/AUDIT_PBH.md): a region [0..k] ends when ANY outer "
                          "shell j > k enters it (R_j <= R_k), not only the neighbour k+1; the first version missed outer "
                          "shells falling through the re-expanding core (a definitional bug: the question asks for re-expansion "
                          "without shell crossing). Also: a local caustic R'(x_j, t) <= 0 measured with a twin "
                          "shell at x_j (1 + 1e-4) (continuum shell crossing; added after a convergence check showed that "
                          "pair crossing alone under-detects caustics at 600 shells: R_coh fell with shell number). Declared.",
    "time_grid": "uniform dt = sqrt(alpha)/100 from the first turnaround to max(t_ta + 2 tau_B) over x <= x_m, with an exact "
                 "early termination once no larger coherent region can qualify (added for the broad lognormal spectra, "
                 "checked to reproduce 11 earlier rows bit-for-bit)",
}


# ---------------- profiles: return (x, K x^2) on a fine grid, x_m ----------------
def musco_profile(a, Kmax, npts=20001):
    x = np.linspace(0, 6.0, npts)
    shape = x ** 2 * np.exp(-x ** (2 * a) / a)               # x_m = 1 is the maximum of x^2 exp(-x^{2a}/a)
    K2 = Kmax * shape / shape.max()
    return x, K2


def peak_psi(kind, sigma=None, rmax=12.0, npts=6001):
    r = np.linspace(0, rmax, npts)
    if kind == "mono":
        psi = np.sinc(r / np.pi)                              # sin(r)/r, k* = 1
        return r, psi
    lnk = np.linspace(-8 * sigma, 8 * sigma, 2001)
    k = np.exp(lnk)
    w = np.exp(-lnk ** 2 / (2 * sigma ** 2))
    psi = np.array([np.trapezoid(w * np.sinc(k * rr / np.pi), lnk) for rr in r])
    return r, psi / psi[0]


def peak_profile(kind, Kmax, sigma=None):
    r, psi = peak_psi(kind, sigma)
    g = r * np.gradient(psi, r)                              # r psi'
    # bound region: first stretch where g < 0
    gmin = g[: np.argmax((g > 0) & (r > 0.5)) or len(g)].min()
    mu = (-1 + math.sqrt(1 - Kmax)) / gmin                   # max of -mu g (2 + mu g) = Kmax
    zeta = mu * psi
    rz = mu * g
    x = r * np.exp(zeta)
    K2 = -rz * (2 + rz)
    return x, K2, mu


def control_profile(Kmax, npts=20001):
    x = np.linspace(0, 3, npts)
    K2 = Kmax * x ** 2 / (1 + np.exp((x - 1) / 0.01))
    return x, K2 / K2.max() * Kmax


def build(x, K2, AH, alpha, n_sh=600):
    """Shells uniform in mass in [0, 1.2 x_m]; mu set so that the central turnaround scale gives A_c H = AH."""
    km = int(np.argmax(K2)); x_m = x[km]
    # bound region contiguous from the centre
    bound_end = km + int(np.argmax(K2[km:] <= 0)) if (K2[km:] <= 0).any() else len(x) - 1
    x_hi = min(1.2 * x_m, x[bound_end - 1])
    # central curvature k0 = lim K2/x^2
    i0 = max(2, int(0.01 * km))
    k0 = K2[i0] / x[i0] ** 2
    # GR central closed FRW: A_c = 2 mu / k0^{3/2}; exact QT A_c from turnaround psi = 1/A^2:
    def A_exact(mu):
        M = mu * x[i0] ** 3; Kc = K2[i0]
        from ltb.growing import turnaround_radius
        Rta = turnaround_radius(M, Kc, alpha)
        psi = psi_of_y(2 * M / Rta ** 3, alpha)
        return 1 / math.sqrt(psi)
    target = AH / H_INF
    mu_gr = target * k0 ** 1.5 / 2                         # GR: A_c = 2 mu / k0^{3/2}
    lo, hi = 0.7 * mu_gr, 1.5 * mu_gr
    while A_exact(hi) < target:
        hi *= 1.5
    mu = brentq(lambda m_: A_exact(m_) - target, lo, hi, xtol=1e-14 * mu_gr)
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(x, K2)
    xs = x_hi * (np.arange(n_sh) + 0.5) / n_sh                 # uniform in x: resolves the central core
    K2s = cs(xs)
    # synchronization slice: background alpha*y_b = 0.1 (QT corrections to the earlier GR growing mode
    # ~ (alpha y)^2/2 <= 0.5 %), and in any case at most 0.3 of the earliest GR turnaround time
    y_i = 0.1 / alpha
    t_i = 2 / (3 * math.sqrt(y_i))
    t_ta_gr_min = min(math.pi * mu * xx ** 3 / kk ** 1.5 for xx, kk in zip(xs, K2s))
    t_i = min(t_i, 0.3 * t_ta_gr_min)
    shells = [GrowingShell(mu * xx ** 3, float(kk), alpha, t_i) for xx, kk in zip(xs, K2s)]
    # twin shells at x (1 + EPS): local caustic criterion R'(x, t) <= 0 (continuum shell crossing)
    xt = xs * (1 + EPS)
    build.twins = [GrowingShell(mu * xx ** 3, float(cs(xx)), alpha, t_i) for xx in xt]
    return shells, xs, x_m, mu, A_exact(mu) * H_INF


def analyse(shells, xs, x_m, alpha, dt=None, twins=None):
    if twins is None:
        twins = getattr(build, "twins", None)
    n = len(shells)
    cl = np.full(n, np.inf)                                 # local caustic time at x_j
    tB = np.array([s.t_B for s in shells]); tta = np.array([s.t_ta for s in shells])
    tauB = np.array([s.tau_B for s in shells])
    core = xs <= x_m
    w = xs[core] ** 2; o = np.argsort(tB[core]); cw = np.cumsum(w[o]) / w.sum()   # mass-weighted percentiles
    tq = lambda q: float(tB[core][o][np.searchsorted(cw, q)])
    d_tB = (tq(0.9) - tq(0.1)) * H_INF
    T0 = tta.min() * 0.999
    T_end = float((tta[core] + 2 * tauB[core]).max())
    if dt is None:
        dt = math.sqrt(alpha) / 100
    n_total = int(math.ceil((T_end - T0) / dt)) + 1
    c = np.full(n - 1, np.inf)
    chunk = max(1, int(1.5e7 // n))
    prevR = prevt = None
    tend_k = tta + 2 * tauB                       # each shell's own 'still expanding' limit
    cum_tend = np.minimum.accumulate(tend_k[:-1])
    cumB = np.maximum.accumulate(tB[:-1])
    n_done = 0
    for a in range(0, n_total, chunk):
        # early termination (exact): once a crossing among shells 0..k* has been found, where k* is the first
        # region that has not bounced yet, no larger region can ever qualify; the scan only has to continue
        # until the T_k of every smaller region is known.
        if a > 0:
            t_now = T0 + dt * (a - 1)
            cmin = np.minimum.accumulate(np.minimum(c, np.minimum(cl[:-1], cl[1:])))
            blocked = np.nonzero(cmin < cumB)[0]
            if len(blocked):
                kb = int(blocked[0])
                if t_now >= float(np.max(np.minimum(cmin[:kb + 1], cum_tend[:kb + 1]))):
                    break
        tt = T0 + dt * np.arange(a, min(a + chunk, n_total))
        n_done = a + len(tt)
        Rk = np.vstack([s.R_of_t(tt) for s in shells])
        if prevR is not None:
            Rk = np.hstack([prevR[:, None], Rk]); tt = np.concatenate([[prevt], tt])
        # intrusion (strict, after the independent audit): region [0..k] is spoiled as soon as ANY shell j > k
        # has R_j <= R_k (not only the neighbour k+1): suffix minimum over the outer shells
        suf = np.minimum.accumulate(Rk[::-1], axis=0)[::-1]
        neg = suf[1:] <= Rk[:-1]
        for j in np.nonzero(np.isinf(c) & neg.any(axis=1))[0]:
            q = int(np.argmax(neg[j]))
            if q == 0:
                c[j] = tt[0]; continue
            g = lambda t, j=j: float(min(s_.R_of_t(t) for s_ in shells[j + 1:]) - shells[j].R_of_t(t))
            try:
                c[j] = brentq(g, tt[q - 1], tt[q], xtol=1e-10)
            except ValueError:
                c[j] = tt[q]
        if twins is not None:
            Tw = np.vstack([s_.R_of_t(tt) for s_ in twins])
            negl = (Tw - Rk) <= 0
            for j in np.nonzero(np.isinf(cl) & negl.any(axis=1))[0]:
                q = int(np.argmax(negl[j]))
                if q == 0:
                    cl[j] = tt[0]; continue
                g = lambda t: float(twins[j].R_of_t(t) - shells[j].R_of_t(t))
                try:
                    cl[j] = brentq(g, tt[q - 1], tt[q], xtol=1e-10)
                except ValueError:
                    cl[j] = tt[q]
        prevR = Rk[:, -1]; prevt = tt[-1]
    c = np.minimum(c, np.minimum(cl[:-1], cl[1:]))         # pair crossing OR a local caustic at either end
    best = dict(R=0.0, k=-1, t=None)
    Bk = -np.inf; Tk = np.inf
    for k in range(n - 1):
        Bk = max(Bk, tB[k]); Tk = min(Tk, c[k], tta[k] + 2 * tauB[k])
        if Tk <= Bk:
            break
        Rk_ = float(shells[k].R_of_t(Tk - 1e-12))
        if Rk_ > best["R"]:
            best = dict(R=Rk_, k=k, t=float(Tk))
    jm = int(np.argmin(np.abs(xs - x_m)))
    jf = int(np.argmin(c)) if np.isfinite(c).any() else None
    return dict(dtB_core_H=d_tB, R_coh_H=best["R"] * H_INF, x_coh_over_xm=float(xs[best["k"]] / x_m) if best["k"] >= 0 else 0.0,
                M_coh_over_Mm=float((xs[best["k"]] / x_m) ** 3) if best["k"] >= 0 else 0.0, n_coh_shells=int(best["k"] + 1),
                ceiling_H=float(shells[jm].R_ta * H_INF), first_bounce_tH=float(tB.min() * H_INF),
                first_crossing_tH=float(c.min() * H_INF) if jf is not None else None,
                first_crossing_x_over_xm=float(xs[jf] / x_m) if jf is not None else None, n_t=n_done, n_t_full=n_total)


def shapes(quick):
    out = [("musco", dict(a=a)) for a in ([1, 6.33] if quick else [1, 2.5, 6.33, 15])]
    out += [("peak", dict(kind="mono"))] + ([] if quick else [("peak", dict(kind="logn", sigma=s)) for s in (0.5, 1, 2)])
    out += [("control", {})]
    return out


def profile_of(fam, par, Kmax):
    if fam == "musco":
        x, K2 = musco_profile(par["a"], Kmax); return x, K2, None
    if fam == "peak":
        return peak_profile(par["kind"], Kmax, par.get("sigma"))
    x, K2 = control_profile(Kmax); return x, K2, None


def main(quick=False):
    d = new_run_dir("c45b_pbh_amended" + ("_quick" if quick else "") + "_" + "_".join(a for a in sys.argv[1:] if a[0].isdigit()))
    write_json(d / "PROTOCOL.json", dict(PROTOCOL, sources=source_hashes(), env=environment()))
    # validation: GR limit of the growing-mode turnaround and crunch times
    from ltb.growing import GrowingShell as GS
    gs = GS(1.0, 0.5, 0.0, 1e-3)
    val = dict(GR_t_ta_rel=abs(gs.t_ta / (math.pi * 1.0 / 0.5 ** 1.5) - 1),
               GR_t_crunch_rel=abs((gs.t_ta + gs.tau_B) / (2 * math.pi / 0.5 ** 1.5) - 1))
    print("validation", val, flush=True)
    rows = []
    als = [float(a) for a in sys.argv[1:] if a.replace('.', '', 1).isdigit()] or ([0.1] if quick else [0.1, 0.03])
    for al in als:
        for fam, par in shapes(quick):
            for Kmax in ([0.6, 0.99] if quick else [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995]):
                for AH in ([1.0, 9.6] if quick else [0.42, 1.0, 2.0, 5.0, 9.6]):
                    t0 = time.time()
                    try:
                        x, K2, mu_z = profile_of(fam, par, Kmax)
                        shells, xs, x_m, mu, AH_exact = build(x, K2, AH, al)
                        from ltb.coherence import analyse_intrusion
                        res = analyse_intrusion(shells, build.twins, xs, x_m, al)
                        res["R_coh_H"] = res["tol0.3"]["R_coh_H"]          # amended verdict measure (loosest)
                        res.update(AH_exact=AH_exact, zeta_amplitude=mu_z)
                    except Exception as exc:
                        res = dict(error=repr(exc)[:200])
                    res.update(alpha=al, family=fam, par=par, Kmax=Kmax, AH=AH, wall=time.time() - t0)
                    rows.append(res)
                    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()
                                      if k in ("alpha", "family", "par", "Kmax", "AH", "R_coh_H", "ceiling_H", "error", "wall")}
                                     | {k: round(res[k]["R_coh_H"], 4) for k in ("strict", "tol0.1") if k in res}), flush=True)
            write_json(d / "results_partial.json", dict(rows=rows))
    pbh = [r for r in rows if "error" not in r and (r["family"] == "peak" or (r["family"] == "musco" and r["par"]["a"] <= 6.33))]
    alphas = sorted(set(r["alpha"] for r in rows))
    ok = {a: [r for r in pbh if r["alpha"] == a and WINDOW[0] <= r["AH"] <= WINDOW[1] and r["R_coh_H"] >= 1] for a in alphas}
    n_err = sum("error" in r and "no turnaround" not in r["error"] for r in rows)   # physical 'no turnaround' rows are excluded and reported
    if not pbh or n_err > 0.05 * len(rows):
        verdict = "INVALID (too many errors)"
    elif all(ok[a] for a in alphas):
        verdict = "PASS"
    elif all(r["R_coh_H"] < 1 / 3 for r in pbh):
        verdict = "FAIL"
    else:
        verdict = "INCONCLUSIVE"
    best = {str(a): max((r for r in pbh if r["alpha"] == a), key=lambda r: r["R_coh_H"], default=None) for a in alphas}
    write_json(d / "results.json", dict(rows=rows, verdict=verdict, best_pbh=best, validation=val, H_INF=H_INF))
    write_manifest(d)
    print("VERDICT", verdict)
    print("BEST", json.dumps({a: (b and {k: b[k] for k in ("family", "par", "Kmax", "AH", "R_coh_H", "ceiling_H")}) for a, b in best.items()}))
    return d


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
