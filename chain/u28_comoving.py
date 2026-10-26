"""U28 (docs/preregistros/ESPEC_U28.md, amendments 3-4): spherical NPQT with dust + a real scalar with potential, in the
dust-comoving synchronous gauge over a whole closed universe,

    ds^2 = -dt^2 + e^{2 Lam} dchi^2 + r^2 dOmega^2,   chi in [0, pi]  (cell centred, parity at both poles).

Field equations: the covariant 2D form of the audited double-null equations of c45b/dn (equations.py; Bueno, Cano,
Hennigar & Murcia eqs. (33)-(34)), kappa = 8 pi G, Y = r^3 h(psi), h(psi) = psi/sqrt(1 - alpha^2 psi^2):
    (E1) grad_a grad_b r = -(kappa/2) r m^{3/2} (T_ab - g_ab T) + g_ab r psi (3m - 2)/2
    (E2) R2/2 = box r / r + 3 alpha^2 psi^3 - 6 alpha^2 psi m J - kappa m^{3/2} p_perp
    m = 1/(1 + alpha^2 y^2), y = Y/r^3, psi = y sqrt(m)   (never 1 - alpha^2 psi^2: 1 - alpha psi reaches 1e-10)
    J = b^2 - kappa^2 m [(rho + p_r)^2/4 - (Pi q)^2]/4,  b = kappa (rho - p_r) sqrt(m)/4 - 3 psi/2
Evolution (amendment 5, Misner-Sharp form): r_t = U, U_t from the tt part of E1, Gamma = e^{-Lam} r' with
Gamma_t = -(kappa/2) r m^{3/2} Pi q (tR part of E1, q = e^{-Lam} phi'), e^Lam = r'/Gamma, Klein-Gordon as a finite volume
with the exact spherical cell volume Vol = int d(r^3)/(3 Gamma) (amendment 7): P = Pi Vol/h,
P_t = [e^{-Lam} r^2 phi']_faces/h - (Vol/h) V'(phi), Y_t = -kappa r^2 (p_r U - Pi q Gamma) (generalized-mass flux, as in dn).
Dust: rho_d = D(chi)/(e^Lam r^2), exact.  The momentum constraint holds identically; the Hamiltonian one,
r^2 psi = 1 + U^2 - Gamma^2, is monitored.  E2 is not evolved (amendment 5: as an evolution equation for Lam it is
ill-posed).  Domain chi in [0, pi/2], cell centred, parity at the pole and mirror symmetry at the equator.
Initial data on a time-symmetric slice (U = 0, phi uniform), with the same discrete derivatives as the evolution."""
import math
import numpy as np
from scipy.optimize import brentq

C_POLE = 0.3516


# ---------------------------------------------------------------- potential (pole, frw_core mode 4) and homogeneous ODE
def V_pole(x, V0, C=C_POLE):
    return V0 * (1 - C / np.sqrt(x * x + C * C))


def dV_pole(x, V0, C=C_POLE):
    q = x * x + C * C
    return V0 * C * x / (q * np.sqrt(q))


def hom_rhs(y, al, V0, Aeff):
    """frw_core.c equations (8 pi G = 1), without its constraint relaxation"""
    a, H, x, xd = y
    rd = 3 * Aeff / a ** 3
    E = rd + 0.5 * xd * xd + V_pole(x, V0); yv = E / 3
    Hd = 1 / a ** 2 - (rd + xd * xd) * (1 + al * al * yv * yv) ** -1.5 / 2
    return np.array([a * H, Hd, xd, -3 * H * xd - dV_pole(x, V0)])


def hom_aeff(st):
    a0, al = st["a0"], st["al"]
    h0 = 1 / a0 ** 2; hh = h0 / math.sqrt(1 - al * al * h0 * h0)
    Em = 3 * hh - (0.5 * st["chid0"] ** 2 + V_pole(st["chi0"], st["V0"]))
    return max(Em, 0.0) * a0 ** 3 / 3


def hom_to_bounce(st, direction=+1, rtol=1e-12):
    """integrate the homogeneous patch from its turnaround (H = 0) forwards (+1) or backwards (-1) to the next bounce
    (H = 0 with a at a minimum); returns the state there and the elapsed time"""
    from scipy.integrate import solve_ivp
    Aeff = hom_aeff(st); al, V0 = st["al"], st["V0"]
    f = lambda t, y: direction * hom_rhs(y, al, V0, Aeff)
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = 1 if direction > 0 else -1
    y0 = [st["a0"], 0.0, st["chi0"], st["chid0"]]
    # leave the turnaround first (H = 0 there), then stop at the next H = 0
    T = 50 * st["a0"]
    s = solve_ivp(f, (0, T), y0, method="DOP853", rtol=rtol, atol=1e-30, events=ev, first_step=1e-6 * st["a0"])
    if not len(s.t_events[0]):
        raise RuntimeError("no bounce found")
    a, H, x, xd = s.y_events[0][0]
    return dict(a=a, H=H, chi=x, chid=xd, t=s.t_events[0][0], Aeff=Aeff)   # physical state (time runs backwards if direction < 0)


# ---------------------------------------------------------------- grid helpers
# parities (pole, equator) of r, U, Gamma, phi, P, Y
PAR = ((-1, 1), (-1, 1), (1, -1), (1, 1), (1, 1), (-1, 1))


def taper(x, x0=0.95):
    """C-infinity step: 1 for x <= x0, 0 at x = 1"""
    t = np.clip((1 - x) / (1 - x0), 0, 1)
    f = lambda s: np.where(s > 0, np.exp(-1 / np.where(s > 0, s, 1)), 0.0)
    return f(t) / (f(t) + f(1 - t))


def ext(f, par):
    """one ghost cell at the pole and one at the equator, with parities par = (pole, equator)"""
    return np.concatenate(([par[0] * f[0]], f, [par[1] * f[-1]]))


def d1(f, par, h):
    e = ext(f, par)
    return (e[2:] - e[:-2]) / (2 * h)


# ---------------------------------------------------------------- the model
class Model:
    def __init__(self, N, alpha, kappa, V0, D, with_field=True):
        self.N, self.h = N, (math.pi / 2) / N
        self.chi = (np.arange(N) + 0.5) * self.h
        self.x = np.sin(self.chi)
        self.alpha, self.kappa, self.V0, self.D, self.field = alpha, kappa, V0, D, with_field

    def V(self, p):
        return V_pole(p, self.V0) if self.field else 0.0 * p

    def dV(self, p):
        return dV_pole(p, self.V0) if self.field else 0.0 * p

    def geom(self, r, Y):
        y = Y / r ** 3
        m = 1 / (1 + (self.alpha * y) ** 2)
        return y, m, y * np.sqrt(m)

    def volume(self, r, G):
        """exact comoving cell volume / (4 pi): int e^Lam r^2 dchi = int d(r^3)/(3 Gamma) over the cell (amendment 7)"""
        re = ext(r, PAR[0]); rf = 0.5 * (re[1:] + re[:-1])
        rf[0] = 0.0
        return (rf[1:] ** 3 - rf[:-1] ** 3) / (3 * G)

    def derived(self, S):
        r, U, G, phi, P, Y = S
        h = self.h
        rp = d1(r, PAR[0], h)
        eL = rp / G
        php = d1(phi, PAR[3], h)
        q = php / eL
        Pi = P * h / self.volume(r, G)
        rho_d = self.D / (eL * r * r)
        V = self.V(phi)
        rho = rho_d + 0.5 * Pi * Pi + 0.5 * q * q + V
        pr = 0.5 * Pi * Pi + 0.5 * q * q - V
        return rp, eL, php, q, Pi, rho_d, rho, pr

    def rhs(self, S):
        r, U, G, phi, P, Y = S
        h, k = self.h, self.kappa
        y, m, psi = self.geom(r, Y)
        rp, eL, php, q, Pi, rho_d, rho, pr = self.derived(S)
        m15 = m * np.sqrt(m)
        Ud = -(k / 2) * r * m15 * pr - r * psi * (3 * m - 2) / 2
        Gd = -(k / 2) * r * m15 * Pi * q
        re, Ee, pe = ext(r, PAR[0]), ext(1 / eL, (1, 1)), ext(phi, PAR[3])
        flux = 0.5 * (Ee[1:] + Ee[:-1]) * (0.5 * (re[1:] + re[:-1])) ** 2 * (pe[1:] - pe[:-1]) / h
        Pd = (flux[1:] - flux[:-1]) / h - self.volume(r, G) / h * self.dV(phi)
        Yd = -k * r * r * (pr * U - Pi * q * G)
        return np.array([U, Ud, Gd, Pi, Pd, Yd])

    def hamiltonian(self, S):
        r, U, G, phi, P, Y = S
        y, m, psi = self.geom(r, Y)
        return float(np.max(np.abs(r * r * psi - (1 + U * U - G * G)) / (1 + U * U)))


def belly(K_of_x):
    """x_e with K(x_e) x_e^2 = 1: the maximum of the areal radius on the time-symmetric slice (amendment 6)"""
    f = lambda xx: K_of_x(np.array([xx]))[0] * xx * xx - 1
    return 1.0 if f(1.0) == 0 else brentq(f, 0.5, 1.5, xtol=1e-16, rtol=1e-15)


def time_symmetric_data(N, alpha, kappa, V0, Ybar, K_of_x, phi0, split, branch, with_field=True, Y_of_x=None, rho_of_x=None):
    """U = 0 slice.  Y(x) = Ybar x^3 (or Y_of_x), K(x) given; r from the Hamiltonian constraint r^2 psi(Y/r^3) = K x^2
    ('small' = bounce root, 'large' = turnaround root); Gamma = sqrt(1 - K x^2); e^Lam = r'/Gamma and
    rho = Y'/(kappa r^2 r') with the discrete derivatives of the evolution; split(rho) -> (Pi, rho_d)."""
    M = Model(N, alpha, kappa, V0, np.zeros(N), with_field)
    M.x = belly(K_of_x) * np.sin(M.chi)                   # shell label: mirror symmetric, closes at the belly
    x = M.x
    Yf = (lambda xx: Ybar * xx ** 3) if Y_of_x is None else Y_of_x
    psi = lambda yy: yy / np.sqrt(1 + (alpha * yy) ** 2)

    def r_of(xx):
        K = K_of_x(np.array([xx]))[0]; Yx = Yf(xx)
        if alpha == 0:
            return Yx / (K * xx * xx)
        F = lambda rr: rr * rr * psi(Yx / rr ** 3) - K * xx * xx
        rpk = (alpha * Yx) ** (1 / 3)                      # the peak of r^2 psi(Y/r^3) lies near (alpha Y)^{1/3}
        rs = np.geomspace(rpk * 1e-6, rpk * 1e6, 6000); v = F(rs)
        k = int(np.argmax(v))
        if v[k] <= 0:
            raise ValueError("no time-symmetric radius")
        if branch == "small":
            lo = np.nonzero(v[:k] < 0)[0][-1]
            return brentq(F, rs[lo], rs[lo + 1], xtol=1e-300, rtol=1e-15, maxiter=500)
        hi = k + int(np.argmax(v[k:] < 0))
        return brentq(F, rs[hi - 1], rs[hi], xtol=1e-300, rtol=1e-15, maxiter=500)
    r = np.array([r_of(xx) for xx in x])
    Y = np.array([Yf(xx) for xx in x])
    G = np.sqrt(np.maximum(1 - K_of_x(x) * x * x, 0.0))
    if np.any(1 - K_of_x(x) * x * x < 0):
        raise ValueError("1 - K x^2 < 0 inside the belly")
    rp = d1(r, PAR[0], M.h); eL = rp / G
    if rho_of_x is not None:                              # analytic density (family b: the one Y was built from)
        rho = rho_of_x(x)
    else:                                                  # analytic Y_x / (kappa r^2 r_x), relative step 1e-6
        dx = 1e-6
        rx = np.array([(r_of(xx * (1 + dx)) - r_of(xx * (1 - dx))) / (2 * dx * xx) for xx in x])
        Yx = np.array([(Yf(xx * (1 + dx)) - Yf(xx * (1 - dx))) / (2 * dx * xx) for xx in x])
        rho = Yx / (kappa * r * r * rx)
    Pi, rho_d = split(rho)
    M.D = rho_d * eL * r * r
    S = np.array([r, np.zeros(N), G, np.full(N, float(phi0)), Pi * M.volume(r, G) / M.h, Y])
    return M, S, dict(rho_d_min=float(np.min(rho_d)), eL_min=float(np.min(eL)), rho=rho)


def rk4_step(M, S, dt):
    k1 = M.rhs(S); k2 = M.rhs(S + dt / 2 * k1); k3 = M.rhs(S + dt / 2 * k2); k4 = M.rhs(S + dt * k3)
    return S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def evolve(M, S, cfl=0.25, t_max=None, mode="bounce", monitor_every=20, r_stop=None, track_rho=False):
    """'bounce': until every cell has bounced (U: - -> +) after a contraction; 'collapse' (GR): until r at the centre
    cell < r_stop.  Returns per-cell bounce times, the largest Hamiltonian residual, the final state."""
    N = M.N; t = 0.0
    tB = np.full(N, np.nan); was_neg = np.zeros(N, bool)
    cmax = M.hamiltonian(S); steps = 0
    dt0 = cfl * M.h * float(np.min(d1(S[0], PAR[0], M.h) / S[2]))
    rho_max = M.derived(S)[6].copy() if track_rho else None
    while True:
        eL = d1(S[0], PAR[0], M.h) / S[2]
        dt = cfl * M.h * float(np.min(eL))
        if mode == "collapse":
            dt = min(dt, 0.02 * float(S[0][0] / max(abs(S[1][0]), 1e-300)))
        if mode == "bounce" and (dt < 1e-7 * dt0 or steps > 2_000_000):      # collapsing step: a caustic forming
            return dict(status="STALLED", t=t, steps=steps, tB=tB, cmax=cmax, S=S, rho_max=rho_max)
        Snew = rk4_step(M, S, dt); steps += 1
        bad = (not np.all(np.isfinite(Snew))) or np.any(Snew[0] <= 0) or np.any(d1(Snew[0], PAR[0], M.h) / Snew[2] <= 0)
        if bad:
            return dict(status="FAILED", t=t, steps=steps, tB=tB, cmax=cmax, S=S)
        U0, U1 = S[1], Snew[1]
        was_neg |= U1 < 0
        cross = was_neg & (U0 < 0) & (U1 >= 0) & np.isnan(tB)
        tB[cross] = t + dt * (-U0[cross]) / (U1[cross] - U0[cross])
        S = Snew; t += dt
        if track_rho:
            rho_max = np.maximum(rho_max, M.derived(S)[6])
        if steps % monitor_every == 0:
            cmax = max(cmax, M.hamiltonian(S))
        if mode == "bounce" and not np.any(np.isnan(tB)):
            return dict(status="COMPLETE", t=t, steps=steps, tB=tB, cmax=max(cmax, M.hamiltonian(S)), S=S, rho_max=rho_max)
        if mode == "collapse" and S[0][0] < r_stop:
            return dict(status="COMPLETE", t=t, steps=steps, tB=tB, cmax=max(cmax, M.hamiltonian(S)), S=S)
        if t_max is not None and t > t_max:
            return dict(status="TIMEOUT", t=t, steps=steps, tB=tB, cmax=cmax, S=S)


# ---------------------------------------------------------------- C45 parents (8 pi G = 1 units of frw_core)
def parent_data(st, N, family, eps, LH_over_AH, branch="small"):
    """family 'a' (curvature, adiabatic split; amendments 3-4), 'b' (dPi/Pi = eps g, dust unperturbed), or None
    (homogeneous).  Background: the homogeneous C45 patch at its previous bounce (time symmetric with the next one)."""
    bw = hom_to_bounce(st, -1)
    a_b, phi0, Pib = bw["a"], bw["chi"], bw["chid"]
    al, V0 = st["al"], st["V0"]
    rho_db = 3 * bw["Aeff"] / a_b ** 3
    Vb = float(V_pole(phi0, V0))
    rhob = rho_db + 0.5 * Pib * Pib + Vb
    Ybar = a_b ** 3 * rhob / 3                              # kappa = 1: Y = r^3 h(psi) = kappa rho r^3 / 3
    fk = 0.5 * Pib * Pib / (rhob - Vb)
    g = lambda xx: np.sinc(xx / LH_over_AH)                # no taper (amendment 6): closure from the belly x_e
    sgn = 1.0 if Pib >= 0 else -1.0
    if family in (None, "a"):
        e = eps if family == "a" else 0.0
        Kf = lambda xx: 1 + e * g(xx)
        split = lambda rho: (sgn * np.sqrt(2 * fk * (rho - Vb)), (1 - fk) * (rho - Vb))
        M, S, info = time_symmetric_data(N, al, 1.0, V0, Ybar, Kf, phi0, split, branch)
    elif family == "b":
        M, S, info = field_family_data(N, al, V0, a_b, phi0, Pib, rho_db, eps, g)
    else:
        raise ValueError(family)
    info.update(a_b=a_b, phi0=phi0, Pib=Pib, rho_db=rho_db, fk=fk, t_half=bw["t"])
    return M, S, info


def field_family_data(N, al, V0, a_b, phi0, Pib, rho_db, eps, g):
    """family (b): K = 1, Pi = Pib (1 + eps g(x)), dust density rho_db; Y(x) from Y_x = kappa rho r^2 r_x with r from
    r^2 psi(Y/r^3) = x^2 (bounce root), integrated from the pole; then the same discrete construction as above."""
    from scipy.integrate import solve_ivp
    Vb = float(V_pole(phi0, V0))
    psi = lambda yy: yy / np.sqrt(1 + (al * yy) ** 2)
    dpsi = lambda yy: (1 + (al * yy) ** 2) ** -1.5
    rho_of = lambda xx: rho_db + 0.5 * (Pib * (1 + eps * g(xx))) ** 2 + Vb

    def r_of(xx, Y):
        F = lambda rr: rr * rr * psi(Y / rr ** 3) - xx * xx
        rpk = (al * Y) ** (1 / 3)
        rs = np.geomspace(rpk * 1e-6, rpk * 1e6, 6000); v = F(rs); k = int(np.argmax(v))
        lo = np.nonzero(v[:k] < 0)[0][-1]
        return brentq(F, rs[lo], rs[lo + 1], xtol=1e-300, rtol=1e-15, maxiter=500)

    def dY(xx, Yv):
        Y = Yv[0]; r = r_of(xx, Y); yy = Y / r ** 3
        Gr = 2 * r * psi(yy) - 3 * dpsi(yy) * yy / r * r     # d/dr [r^2 psi(Y/r^3)]
        GY = dpsi(yy) / r                                     # d/dY
        rho = rho_of(np.array([xx]))[0]
        rx = 2 * xx / (Gr + rho * r * r * GY)                 # kappa = 1
        return [rho * r * r * rx]
    x0 = 1e-4
    rho0 = rho_of(np.array([0.0]))[0]; r0 = x0 / math.sqrt(psi(rho0 / 3))
    Mtmp = Model(N, al, 1.0, V0, np.zeros(N))
    s = solve_ivp(dY, (x0, float(Mtmp.x[-1])), [rho0 * r0 ** 3 / 3], t_eval=Mtmp.x, rtol=1e-12, atol=1e-300, method="DOP853")
    from scipy.interpolate import CubicSpline
    xs = np.concatenate(([0.0], Mtmp.x)); Ys = np.concatenate(([0.0], s.y[0]))
    cs = CubicSpline(xs, Ys)
    Yfun = lambda xx: float(cs(xx))
    # the bounce root again on the grid, then the shared discrete construction (rho from the discrete Y')
    split = lambda rho: (Pib * (1 + eps * g(Mtmp.x)), rho - 0.5 * (Pib * (1 + eps * g(Mtmp.x))) ** 2 - Vb)
    return time_symmetric_data(N, al, 1.0, V0, None, lambda xx: np.ones_like(xx), phi0, split, "small", Y_of_x=Yfun,
                               rho_of_x=rho_of)


# ---------------------------------------------------------------- separate-universe continuation (amendment 3, point 5)
def hom_continue(a, H, chi, chid, rho_d, st, efolds=60.0, t_max_fac=240.0):
    """continue one cell as a homogeneous patch (frw_core equations, H != 0 allowed) until ln(a/a0) = efolds (launch),
    a < 1e-3 sqrt(al) (crunch) or t > t_max_fac a0; returns outcome, chi_end, chi_abs_max, E_max"""
    from scipy.integrate import solve_ivp
    al, V0, a0 = st["al"], st["V0"], st["a0"]
    Aeff = rho_d * a ** 3 / 3
    f = lambda t, y: hom_rhs(y, al, V0, Aeff)
    e1 = lambda t, y: math.log(y[0] / a0) - efolds; e1.terminal = True; e1.direction = 1
    e2 = lambda t, y: y[0] - 1e-3 * math.sqrt(al); e2.terminal = True
    s = solve_ivp(f, (0, t_max_fac * a0), [a, H, chi, chid], method="DOP853", rtol=1e-10, atol=1e-30, events=[e1, e2])
    A, Hs, X, Xd = s.y
    E = 3 * Aeff / A ** 3 + 0.5 * Xd * Xd + V_pole(X, V0)
    out = "launched60" if len(s.t_events[0]) else ("crunch" if len(s.t_events[1]) else "cycling")
    return dict(outcome=out, chi_end=float(X[-1]), chi_abs_max=float(np.max(np.abs(X))), E_max=float(E.max()))
