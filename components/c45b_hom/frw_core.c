/* C45g: C port of hom/frw.py run() (audited C45 equations, 8 pi G = 1), Dormand-Prince 5(4) adaptive.
   State y = (a, H, chi, chidot).  Equations (identical to the Python port, including its constraint relaxation):
     z = al (H^2 + 1/a^2) (clipped at 1 - 1e-14)
     Hdot = 1/a^2 - (rho_d + chidot^2) (1 - z^2)^{3/2} / 2,   rho_d = 3 Aeff / a^3
     E = rho_d + chidot^2/2 + V0 Vn(chi); psic = (E/3)/sqrt(1 + al^2 (E/3)^2); H2c = psic - 1/a^2
     if H2c > 0.05 psic:  Hc = sign(H) sqrt(H2c);  Hdot -= 10 |Hc| (H - Hc);  H -> Hc
     adot = a H, chiddot = -3 H chidot - V0 Vn'(chi)
   Online diagnostics as in run(): n_bounces (H: - -> +), chi_abs_max, N_inf (inflating steps), constraint max.
   Declared differences from the Python port: (1 - z^2)^{3/2} on shell, and the constraint measured as |psi - psi(E/3)|/psi(E/3)
   (the Python measure |3h(psi) - E|/E amplifies round-off near the throat by 1/(1 - z^2)).
   Compile: gcc -O3 -shared -fPIC -o frw_core.so frw_core.c -lm */
#include <math.h>

#define LH 10
static double W[LH], SC[LH], Z;
static double AL, V0, AEFF, BRAD;   /* dust rho_d = 3 AEFF/a^3; radiation rho_r = 3 BRAD/a^4 (C45h) */

static void init_pot(void) {
    Z = 0; for (int n = 0; n < LH; n++) { W[n] = pow(0.5, n); SC[n] = 0.5 * pow(2.0, n); Z += W[n]; }
}
static double Vh(double x) { double s = 0; for (int n = 0; n < LH; n++) { double t = tanh(x / SC[n]); s += W[n] * t * t; } return s / Z; }
static double dVh(double x) { double s = 0; for (int n = 0; n < LH; n++) { double t = tanh(x / SC[n]); s += W[n] * 2 * t * (1 - t * t) / SC[n]; } return s / Z; }
/* C45i: large-field deformation V = V0 Vh(chi) g(chi); mode 0 none, 1 rise g = 1 + k (chi/chic)^2, 2 fall g = 1/(1 + (chi/chic)^2),
   3 tilt g = 1 + k |chi|/chic.  Default 0 (identical to the undeformed potential). */
static int DMODE = 0; static double DCHIC = 1e30, DK = 0;
void set_deform(int mode, double chic, double k) { DMODE = mode; DCHIC = chic; DK = k; }
static double gdef(double x) { double u = x / DCHIC;
    switch (DMODE) { case 1: return 1 + DK * u * u; case 2: return 1 / (1 + u * u); case 3: return 1 + DK * fabs(u); default: return 1; } }
static double dgdef(double x) { double u = x / DCHIC;
    switch (DMODE) { case 1: return 2 * DK * u / DCHIC; case 2: return -2 * u / DCHIC / ((1 + u * u) * (1 + u * u));
                     case 3: return DK * (x >= 0 ? 1 : -1) / DCHIC; default: return 0; } }
/* U22 (adapted for unification, ESPEC_U22): mode 4 replaces the plateau by the regularized pole-inflation potential
   V = 1 - C/sqrt(chi^2 + C^2) with C = chic (exact V0 (1 - C/chi) tail).  Modes 0-3 unchanged. */
static double Vn(double x) { if (DMODE == 4) return 1 - DCHIC / sqrt(x * x + DCHIC * DCHIC); return Vh(x) * gdef(x); }
static double dVn(double x) { if (DMODE == 4) { double q = x * x + DCHIC * DCHIC; return DCHIC * x / (q * sqrt(q)); }
                              return dVh(x) * gdef(x) + Vh(x) * dgdef(x); }

static void rhs(const double *y, double *f) {
    double a = y[0], H = y[1], x = y[2], xd = y[3];
    double rd = 3 * AEFF / (a * a * a), rr = 3 * BRAD / (a * a * a * a);
    double E = rd + rr + 0.5 * xd * xd + V0 * Vn(x), yv = E / 3;
    /* (1 - z^2)^{3/2} evaluated ON SHELL as (1 + al^2 y^2)^{-3/2}, y = E/3 (identical when the constraint holds; the
       Python form 1 - (al psi)^2 underflows at strong scale separation, where 1 - z ~ 1/(2 al^2 y^2) < 1e-16). */
    double Hd = 1 / (a * a) - (rd + 4.0 / 3.0 * rr + xd * xd) * pow(1 + AL * AL * yv * yv, -1.5) / 2;   /* rho + p */
    double psic = yv / sqrt(1 + AL * AL * yv * yv), H2c = psic - 1 / (a * a);
    if (H2c > 0.05 * psic) { double Hc = copysign(sqrt(H2c), H); Hd -= 10 * fabs(Hc) * (H - Hc); H = Hc; }
    f[0] = a * H; f[1] = Hd; f[2] = xd; f[3] = -3 * H * xd - V0 * dVn(x);
}

/* returns outcome: 1 launched60, 2 crunch, 0 cycling (window end), -1 step limit, -2 step underflow
   out[]: n_bounces, chi_abs_max, chi_end, N_inf, constraint_max, t_end, a_min, z_max, efolds_total, n_steps,
          E_max (max total energy density, 8 pi G = 1 units, i.e. reduced-Planck^4), K_max (max chi kinetic density) */
int run(double a0, double al, double v0, double chi0, double chid0, double T, double efolds_stop, double rtol,
        long max_steps, double chi_end_thr, double *out, double frac_rad) {
    init_pot(); AL = al; V0 = v0;
    double h0 = 1 / (a0 * a0); double hh = h0 / sqrt(1 - al * al * h0 * h0);
    double Echi = 0.5 * chid0 * chid0 + v0 * Vn(chi0);
    /* matter budget at turnaround from the exact constraint; frac_rad of it in radiation, the rest dust (C45h) */
    double Em = 3 * hh - Echi; if (Em < 0) Em = 0;
    AEFF = (1 - frac_rad) * Em * a0 * a0 * a0 / 3;
    BRAD = frac_rad * Em * a0 * a0 * a0 * a0 / 3;
    static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
        a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
        a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
        b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
        e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
    double y[4] = {a0, 0.0, chi0, chid0}, k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4], yt[4], yn[4];
    double t = 0, h = 1e-3 * sqrt(al) + 1e-6 * a0;
    double scale[4];
    int nb = 0, outcome = 0; double xmax = fabs(chi0), Ninf = 0, cmax = 0, amin = a0, zmax = 0, Emax = 0, Kmax = 0;
    long ns = 0;
    rhs(y, k1);
    while (t < T) {
        if (ns++ > max_steps) { outcome = -1; break; }
        if (t + h > T) h = T - t;
        for (int i = 0; i < 4; i++) yt[i] = y[i] + h * a21 * k1[i]; rhs(yt, k2);
        for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]); rhs(yt, k3);
        for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); rhs(yt, k4);
        for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); rhs(yt, k5);
        for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); rhs(yt, k6);
        for (int i = 0; i < 4; i++) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        rhs(yn, k7);
        /* error norm: per-component scale; H and chidot get absolute floors tied to their natural sizes */
        scale[0] = rtol * fmax(fabs(y[0]), fabs(yn[0]));
        scale[1] = rtol * fmax(fmax(fabs(y[1]), fabs(yn[1])), 1 / y[0]) ;
        scale[2] = rtol * fmax(fmax(fabs(y[2]), fabs(yn[2])), 1.0);
        scale[3] = rtol * fmax(fmax(fabs(y[3]), fabs(yn[3])), sqrt(v0));
        double err = 0;
        for (int i = 0; i < 4; i++) {
            double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]) / scale[i];
            err += ei * ei;
        }
        err = sqrt(err / 4);
        if (!(err <= 1.0)) {
            double fac = (err == err) ? fmax(0.2, 0.9 * pow(err, -0.2)) : 0.2;
            h *= fac; if (h < 1e-14 * (1 + t)) { outcome = -2; break; }
            continue;
        }
        /* accepted: diagnostics on the step start (as the Python port: inflating mask at the start of each step) */
        {
            double a = y[0], H = y[1], x = y[2], xd = y[3];
            double rd = 3 * AEFF / (a * a * a), rr = 3 * BRAD / (a * a * a * a), V = v0 * Vn(x);
            double yv = (rd + rr + 0.5 * xd * xd + V) / 3;
            double Hd = 1 / (a * a) - (rd + 4.0 / 3.0 * rr + xd * xd) * pow(1 + AL * AL * yv * yv, -1.5) / 2;
            if (H > 0 && Hd + H * H > 0 && V > xd * xd && fabs(x) > chi_end_thr) Ninf += log(yn[0] / y[0]);
        }
        if (y[1] < 0 && yn[1] >= 0) nb++;
        for (int i = 0; i < 4; i++) y[i] = yn[i];
        for (int i = 0; i < 4; i++) k1[i] = k7[i];
        t += h;
        {
            double a = y[0], H = y[1], x = y[2], xd = y[3];
            double psi = H * H + 1 / (a * a), z = AL * psi; if (z > 1 - 1e-14) z = 1 - 1e-14;
            double E = 3 * AEFF / (a * a * a) + 3 * BRAD / (a * a * a * a) + 0.5 * xd * xd + v0 * Vn(x);
            double yv = E / 3, psiy = yv / sqrt(1 + AL * AL * yv * yv);
            double c = fabs(psi - psiy) / psiy; if (c > cmax) cmax = c;   /* constraint: psi = psi(E/3) (relative) */
            if (fabs(x) > xmax) xmax = fabs(x); if (a < amin) amin = a; if (z > zmax) zmax = z;
            if (E > Emax) Emax = E; if (0.5 * xd * xd > Kmax) Kmax = 0.5 * xd * xd;
        }
        if (y[0] < 1e-3 * sqrt(al)) { outcome = 2; break; }
        if (log(y[0] / a0) >= efolds_stop) { outcome = 1; break; }
        double fac = (err > 0) ? fmin(5.0, 0.9 * pow(err, -0.2)) : 5.0;
        h *= fac;
    }
    out[0] = nb; out[1] = xmax; out[2] = y[2]; out[3] = Ninf; out[4] = cmax; out[5] = t; out[6] = amin; out[7] = zmax;
    out[8] = log(y[0] / a0); out[9] = (double)ns; out[10] = Emax; out[11] = Kmax;
    return outcome;
}
