"""U3a · minimal linear Einstein-Boltzmann solver (synchronous gauge, Ma & Bertschinger 1995) for the matter transfer of
a mixed cold + non-cold dark matter with an ARBITRARY momentum distribution f0(q) (ESPEC_U3.md, U3a).

Species: CDM; baryons that do not cluster before decoupling (a < 1/1090; Silk damping erases them on the scales of
interest, and the error at k < 1 h/Mpc is common to all runs) and are pressureless afterwards; massless radiation (photons + neutrinos, one free-streaming
hierarchy, N_eff = 3.044) with the radiation-streaming approximation after k tau > 45 (Blas, Lesgourgues & Tram 2011,
synchronous form: delta_r = 4 H h'/k^2 - 4 eta, theta_r = -h'/2), and one non-cold species with the full momentum-dependent
hierarchy Psi_l(q).  All k are integrated together with a fixed-order RK4 in ln(tau), with an adaptive step bounded by the
fastest active oscillation.  Only RATIOS to the CDM run with the same code are used (validated against Viel et al. 2005).

Units: Mpc, conformal time tau [Mpc], c = 1.  The non-cold species is described by x = q/m (today) only (scale invariance)."""
import math
import numpy as np

# ---- cosmology (Planck 2018-like; the same for every run, so it cancels in the ratios) ----
h = 0.674; OMB = 0.0224 / h ** 2; OMDM = 0.120 / h ** 2
OG = 2.473e-5 / h ** 2; ONU = 3.044 * 7 / 8 * (4 / 11) ** (4 / 3) * OG; OR = OG + ONU; RNU = ONU / OR
H0 = h / 2997.92458                     # Mpc^-1
NS = 0.9649
LMAX_R = 8
RSA_KTAU = 45.0
A_DEC = 1 / 1090.0


class Species:
    """non-cold component: nodes x_j = q_j/m (today), weights W_j (number per node, arbitrary norm), dlnf/dlnq at nodes,
    density today Omega (NR limit)."""
    def __init__(self, x, W, dlnf, Omega, lmax=12, label=""):
        self.x = np.asarray(x, float); self.W = np.asarray(W, float); self.dlnf = np.asarray(dlnf, float)
        self.Omega = float(Omega); self.lmax = int(lmax); self.label = label
        self.norm = float(np.sum(self.W * np.sqrt(self.x ** 2 + 1.0)))        # sum W eps(a=1)  (eps in units of m)

    def rho(self, a):        # in units of rho_crit0
        if self.Omega == 0: return 0.0, 0.0
        eps = np.sqrt(self.x ** 2 + a * a)
        r = self.Omega * np.sum(self.W * eps) / self.norm / a ** 4
        p = self.Omega * np.sum(self.W * self.x ** 2 / (3 * eps)) / self.norm / a ** 4
        return r, p


def thermal_fd(m_eV, omega_wdm, nq=48, lmax=12):
    """thermal fermion WDM (Viel et al. convention Omega h^2 = (T_x/T_nu)^3 m/94 eV); returns Species."""
    Tnu = (4 / 11) ** (1 / 3) * 2.7255 * 8.617333262e-5
    Tx = (omega_wdm * 94.0 / m_eV) ** (1 / 3) * Tnu
    y = np.geomspace(0.03, 25.0, nq); dl = np.gradient(np.log(y))
    f = 1 / (np.exp(y) + 1); W = y ** 3 * f * dl                                # number per node  ~ q^3 f dln q
    dlnf = -y * np.exp(y) / (np.exp(y) + 1)
    return Species(y * Tx / m_eV, W, dlnf, omega_wdm / h ** 2, lmax, f"thermal {m_eV/1e3:g} keV")


def from_packets(x_pk, n_pk, Omega, nq=48, bw=0.15, qtail=0.999, lmax=12):
    """KDE in ln x of number-weighted packets -> Species; drops the number tail above the qtail quantile (reported)."""
    lx = np.log(np.asarray(x_pk, float)); w = np.asarray(n_pk, float); keep = w > 0
    lx, w = lx[keep], w[keep]; o = np.argsort(lx); lx, w = lx[o], w[o]; c = np.cumsum(w) / w.sum()
    lo = lx[np.searchsorted(c, 0.0005)] - 3 * bw; hi = lx[min(np.searchsorted(c, qtail), len(lx) - 1)] + 3 * bw
    sel = lx <= lx[min(np.searchsorted(c, qtail), len(lx) - 1)]
    tail_number = float(1 - w[sel].sum() / w.sum())
    grid = np.linspace(lo, hi, nq); dl = grid[1] - grid[0]
    G = np.zeros(nq); dG = np.zeros(nq)
    for s in range(0, len(lx), 4000):
        d = grid[:, None] - lx[None, s:s + 4000]
        ww = (w[s:s + 4000] * sel[s:s + 4000])[None, :]
        k = np.exp(-0.5 * (d / bw) ** 2) / (bw * math.sqrt(2 * math.pi))
        G += np.sum(ww * k, 1); dG += np.sum(ww * k * (-d / bw ** 2), 1)
    G = np.maximum(G, 1e-300 * G.max())
    dlnf = dG / G - 3.0                                                        # f ~ G / q^3
    return Species(np.exp(grid), G * dl, dlnf, Omega, lmax, "packets"), tail_number


class Background:
    def __init__(self, sp=None, Omega_c=None):
        self.sp = sp
        self.Oc = (OMB + OMDM - (sp.Omega if sp else 0.0)) if Omega_c is None else Omega_c     # all pressureless (CDM + baryons)
        self.Ob = OMB
        r1, _ = sp.rho(1.0) if sp else (0.0, 0.0)
        self.OL = 1 - self.Oc - OR - r1
        la = np.linspace(math.log(1e-13), 0.0, 6000); a = np.exp(la)
        Hc = np.array([self.conf_H(v) for v in a])                             # conformal Hubble aH
        dtau = 1 / Hc                                                          # dtau/dln a = 1/(aH)
        tau = np.concatenate([[a[0] / (H0 * math.sqrt(OR))], np.zeros(len(a) - 1)])
        tau[1:] = tau[0] + np.cumsum(0.5 * (dtau[1:] + dtau[:-1]) * np.diff(la))
        self.lt, self.la = np.log(tau), la

    def rho_parts(self, a):
        rw, pw = self.sp.rho(a) if self.sp else (0.0, 0.0)
        return self.Oc / a ** 3, OR / a ** 4, rw, pw, self.OL

    def conf_H(self, a):
        c, r, w, _, L = self.rho_parts(a)
        return H0 * a * math.sqrt(c + r + w + L)

    def a_of(self, lt):
        return float(np.exp(np.interp(lt, self.lt, self.la)))


def solve(kvals, sp=None, a_end=0.25, dlt_max=0.004, phase_max=0.25):
    """returns delta_m(k) at a_end (cold + non-cold, density weighted) for the adiabatic growing mode with C = 1."""
    k = np.asarray(kvals, float); nk = len(k)
    bg = Background(sp)
    nq = len(sp.x) if sp else 0; lw = sp.lmax if sp else 0
    tau_i = 0.005 / k.max(); a_i = bg.a_of(math.log(tau_i)); tau_e = math.exp(np.interp(math.log(a_end), bg.la, bg.lt))
    # state: eta, dc, F[0..LMAX_R], Psi[q, 0..lw]
    nF = LMAX_R + 1
    eta = np.zeros(nk); dc = np.zeros(nk); db = np.zeros(nk); F = np.zeros((nk, nF)); Psi = np.zeros((nk, nq, lw + 1)) if sp else None
    kt = k * tau_i; C = 1.0
    R = RNU
    eta[:] = 2 * C - (5 + 4 * R) / (6 * (15 + 4 * R)) * C * kt ** 2
    dr = -(2 / 3) * C * kt ** 2; dc[:] = 0.75 * dr
    th_nu = -(1 / 18) * C * k ** 4 * tau_i ** 3 * (23 + 4 * R) / (15 + 4 * R); th_g = -(1 / 18) * C * k ** 4 * tau_i ** 3
    th_r = R * th_nu + (1 - R) * th_g; sig_r = R * 4 * C * kt ** 2 / (3 * (15 + 4 * R))
    F[:, 0] = dr; F[:, 1] = 4 * th_r / (3 * k); F[:, 2] = 2 * sig_r
    if sp:
        eps = np.sqrt(sp.x ** 2 + a_i ** 2)
        Psi[:, :, 0] = -0.25 * dr[:, None] * sp.dlnf[None, :]
        Psi[:, :, 1] = -(eps / (3 * sp.x))[None, :] / k[:, None] * th_nu[:, None] * sp.dlnf[None, :]
        Psi[:, :, 2] = -0.5 * (R * 4 * C * kt ** 2 / (3 * (15 + 4 * R)) / R if R else 0)[:, None] * sp.dlnf[None, :]
    l_r = np.arange(nF); l_w = np.arange(lw + 1)

    def rhs(lt, eta, dc, F, Psi, db):
        tau = math.exp(lt); a = bg.a_of(lt)
        c, r, w, pw, L = bg.rho_parts(a)
        Hc = H0 * a * math.sqrt(c + r + w + L)
        G = 1.5 * H0 ** 2 * a * a                                              # 4 pi G a^2 rho = G * (rho / rho_crit0)
        rsa = k * tau > RSA_KTAU
        drho_w = np.zeros(nk); thw = np.zeros(nk); sigw = np.zeros(nk)
        if sp:
            eps = np.sqrt(sp.x ** 2 + a * a); nrm = sp.Omega / sp.norm / a ** 4
            drho_w = nrm * Psi[:, :, 0] @ (sp.W * eps)
            thw = nrm * k * (Psi[:, :, 1] @ (sp.W * sp.x))
            sigw = nrm * (2 / 3) * (Psi[:, :, 2] @ (sp.W * sp.x ** 2 / eps))
        cb = bg.Ob / a ** 3; cc = c - cb
        drho_m = cc * dc + cb * db + drho_w
        hp = np.where(rsa,
                      (2 / Hc) * (k * k * eta + G * (drho_m - 4 * r * eta)) / (1 - 8 * G * r / (k * k)),
                      (2 / Hc) * (k * k * eta + G * (drho_m + r * F[:, 0])))
        th_r = np.where(rsa, -0.5 * hp, 0.75 * k * F[:, 1])
        etap = G * ((4 / 3) * r * th_r + thw) / (k * k)
        dF = np.zeros_like(F)
        full = ~rsa
        if full.any():
            kk = k[:, None]
            Fp1 = np.concatenate([F[:, 1:], ((2 * LMAX_R + 1) / (k * tau) * F[:, -1] - F[:, -2])[:, None]], 1)
            Fm1 = np.concatenate([np.zeros((nk, 1)), F[:, :-1]], 1)
            dF = kk / (2 * l_r + 1) * (l_r * Fm1 - (l_r + 1) * Fp1)
            dF[:, 0] = -k * F[:, 1] - (2 / 3) * hp
            dF[:, 1] = (k / 3) * (F[:, 0] - 2 * F[:, 2])
            dF[:, 2] = (k / 5) * (2 * F[:, 1] - 3 * F[:, 3]) + (4 / 15) * hp + (8 / 5) * etap
            dF[rsa] = 0.0
        dPsi = None
        if sp:
            eps = np.sqrt(sp.x ** 2 + a * a); om = (sp.x / eps)[None, :, None] * k[:, None, None]   # q k / eps
            Pp1 = np.concatenate([Psi[:, :, 1:], ((2 * lw + 1) * (eps / sp.x)[None, :] / (k[:, None] * tau) * Psi[:, :, -1] - Psi[:, :, -2])[:, :, None]], 2)
            Pm1 = np.concatenate([np.zeros((nk, nq, 1)), Psi[:, :, :-1]], 2)
            dPsi = om / (2 * l_w + 1) * (l_w * Pm1 - (l_w + 1) * Pp1)
            dl = sp.dlnf[None, :]
            dPsi[:, :, 0] = -om[:, :, 0] * Psi[:, :, 1] + (1 / 6) * hp[:, None] * dl
            dPsi[:, :, 1] = om[:, :, 0] / 3 * (Psi[:, :, 0] - 2 * Psi[:, :, 2])
            dPsi[:, :, 2] = om[:, :, 0] / 5 * (2 * Psi[:, :, 1] - 3 * Psi[:, :, 3]) - ((1 / 15) * hp + (2 / 5) * etap)[:, None] * dl
        # d/dln tau = tau d/dtau
        dbp = (-0.5 * hp) if a >= A_DEC else np.zeros(nk)
        return tau * etap, tau * (-0.5 * hp), tau * dF, (tau * dPsi if sp else None), tau * dbp, dict(c=c, w=w, drho_w=drho_w, cc=cc, cb=cb)

    lt, lt_e = math.log(tau_i), math.log(tau_e)
    state = [eta, dc, F, Psi, db]
    nsteps = 0
    while lt < lt_e - 1e-12:
        tau = math.exp(lt); a = bg.a_of(lt)
        fr = np.where(k * tau > RSA_KTAU, 0.0, k * tau).max()
        if sp:
            vmax = float(np.max(sp.x / np.sqrt(sp.x ** 2 + a * a)))
            fr = max(fr, k.max() * tau * vmax)
        dlt = min(dlt_max, phase_max / max(fr, 1e-30), lt_e - lt)
        def add(s, d, c):
            return [s[0] + c * d[0], s[1] + c * d[1], s[2] + c * d[2], (s[3] + c * d[3]) if sp else None, s[4] + c * d[4]]
        k1 = rhs(lt, *state)[:5]
        k2 = rhs(lt + dlt / 2, *add(state, k1, dlt / 2))[:5]
        k3 = rhs(lt + dlt / 2, *add(state, k2, dlt / 2))[:5]
        k4 = rhs(lt + dlt, *add(state, k3, dlt))[:5]
        state = [state[i] + dlt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) if (i != 3 or sp) else None for i in range(5)]
        lt += dlt; nsteps += 1
    aux = rhs(lt, *state)[5]
    dm = (aux["cc"] * state[1] + aux["cb"] * state[4] + aux["drho_w"]) / (aux["c"] + aux["w"])
    return dm, dict(nsteps=nsteps, a_end=a_end)


# ---------------------------------------------------------------- Lyman-alpha estimator (Murgia et al. 2017)
def p_cdm_shape(k, dm):
    """P(k) ~ k^(ns-4) delta_m^2 (curvature normalisation C = 1 per mode); arbitrary overall amplitude."""
    return k ** (NS - 4) * dm ** 2


def p1d(kgrid, P, kout, slope_tail=None):
    """P_1D(k) = (1/2 pi) int_k^inf k' P(k') dk' on a log grid; power-law tail beyond kgrid[-1]."""
    lk = np.log(kgrid); y = kgrid * kgrid * P                                  # integrand in dln k': k'^2 P
    if slope_tail is None:
        slope_tail = np.polyfit(lk[-6:], np.log(np.maximum(P[-6:], 1e-300)), 1)[0]
    tail = kgrid[-1] ** 2 * P[-1] / (-(slope_tail + 2)) if slope_tail < -2 else np.inf
    out = []
    for kk in kout:
        m = kgrid >= kk
        xs = np.concatenate([[math.log(kk)], lk[m]]); ys = np.concatenate([[np.interp(math.log(kk), lk, y)], y[m]])
        out.append((np.sum(0.5 * (ys[1:] + ys[:-1]) * np.diff(xs)) + tail) / (2 * math.pi))
    return np.array(out), slope_tail


def deltaA(kgrid_h, T2, Pcdm, kmin=0.5, kmax=20.0, n=400):
    """delta A = 1 - A/A_CDM, A = int_{kmin}^{kmax} P1D/P1D_CDM dk  (k in h/Mpc).  T2 beyond kgrid kept constant."""
    ko = np.linspace(kmin, kmax, n)
    p0, sl = p1d(kgrid_h, Pcdm, ko)
    # model: P = Pcdm T2 with T2 -> T2[-1] beyond the grid (same tail slope)
    p1, _ = p1d(kgrid_h, Pcdm * T2, ko, slope_tail=sl)
    # (the tail beyond the grid is the CDM power law times T2[-1], already included by p1d)
    r = p1 / p0
    A = np.sum(0.5 * (r[1:] + r[:-1]) * np.diff(ko))
    return float(1 - A / (kmax - kmin))


def viel_T(k_h, m_keV, Omega_x):
    al = 0.049 * m_keV ** -1.11 * (Omega_x / 0.25) ** 0.11 * (h / 0.7) ** 1.22
    nu = 1.12
    return (1 + (al * k_h) ** (2 * nu)) ** (-5 / nu)
