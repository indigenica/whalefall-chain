"""U62, post-verdict diagnostics (informative, written after seeing the U62 verdict): (i) the delta-chi content of each normal
mode at the window ends (kinetic-weighted fraction in x), (ii) the adiabaticity |dw/dt|/w^2 and k/aH at the ends, (iii) how far
past phi = +0.2 the frozen energy E of the uniform-radiation gauge stays positive (the U62 Sec. 3 assumption).
    python3 -B chain/u62_diagnostico.py   ->  results/u62_runs/diagnostico.json
Internal, not exported."""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from ucommon import ROOT
import explora_cuscuton as EC
import u59_escalares_rebote as U
import u62_bogoliubov_rebote as M

A = U.Action()


def bg_mode(cu, Crad, y, mode):
    Q, _ = EC.deriv(cu, mode, Crad, 0.0); qq = Q(y); phi, Uv, dU, ddU = cu.all(mode, y[3]); chi = y[1]; Cp, V0 = EC.CP, EC.V0
    return dict(a=math.exp(y[0]), H=qq["H"], chid=y[2], phid=qq["phidot"], q=(4 * qq["rho"] / 3) ** 0.25, mu2=cu.mu2, Vv=EC.Vp(chi),
                Vp=EC.dVp(chi), Vpp=V0 * Cp * (Cp * Cp - 2 * chi * chi) / (Cp * Cp + chi * chi) ** 2.5, Uv=Uv, Up=dU, Upp=ddU)


def emat(bg, n):
    K, B, C, _ = U.red2(A, bg, n); E = U.energy_matrix(K, B, C); return K, 0.5 * (E + E.T)


def xfrac(bg, n):
    K, E = emat(bg, n); w, F = M.normal_modes(E)
    return [float(K[1, 1] * abs(F[1, j]) ** 2 / np.real(F[:2, j].conj() @ K @ F[:2, j])) for j in range(2)], w


def main():
    out = []
    for h in M.visible_hosts():
        par = h["par"]; cu, Crad, yb, _ = U.host_state(par); tr = U.Trajectory(cu, Crad, yb)
        s = math.sin(M.D_LSS * math.sqrt(abs(h["OmegaK0"]))); nobs = {l: int(round(l / s)) for l in M.ELLS}
        # continuation past the window in mode +1 (d = phi_inf - phi), 10 e-folds
        y0 = list(tr.y(tr.tp)); y0[3] = cu.pinf - y0[3]; _, rhs = EC.deriv(cu, 1, Crad, 0.0)
        ev = lambda t, y: y[0] - y0[0] - 10; ev.terminal = True
        sol = solve_ivp(rhs, (tr.tp, tr.tp + 1e13), y0, method="DOP853", rtol=1e-12, atol=[1e-14, 1e-15, 1e-24, 1e-24], dense_output=True, events=ev)
        rows = []
        for n in sorted(set(nobs.values()) | {1000, 10000, 100000}):
            fx, w = xfrac(tr.bg(tr.tp), n); bgf = tr.bg(tr.tp); kc = math.sqrt(n * (n + 2.0)) / bgf["a"]
            ad = []
            for t, sg in ((tr.tm, 1), (tr.tp, -1)):
                dt = 1e-4 * abs(t); w0 = M.normal_modes(emat(tr.bg(t), n)[1])[0]; w1 = M.normal_modes(emat(tr.bg(t + sg * dt), n)[1])[0]
                ad.append((np.abs(w1 - w0) / dt / w0 ** 2).tolist())
            lost = None
            for dl in np.geomspace(1e-3, 10, 41):
                t = brentq(lambda tt: sol.sol(tt)[0] - y0[0] - dl, sol.t[0], sol.t[-1]) if dl < 10 else sol.t[-1]
                if np.min(np.linalg.eigvalsh(emat(bg_mode(cu, Crad, list(sol.sol(t)), 1), n)[1])) <= 0: lost = float(dl); break
            rows.append(dict(n=n, ell=[l for l, v in nobs.items() if v == n], w_over_k_end=(w / kc).tolist(), x_fraction_end=fx,
                             k_over_aH_end=kc / abs(bgf["H"]), adiab_start=ad[0], adiab_end=ad[1], E_indefinite_from_dlna=lost))
            print(par["AH"], "%.3f" % par["f"], n, "w/k", np.round(w / kc, 3), "x-frac", np.round(fx, 3), "k/aH %.1f" % (kc / abs(bgf["H"])),
                  "adiab", np.round(ad[1], 3), "E<0 from dlna", lost, flush=True)
        out.append(dict(par=par, n_obs=nobs, rows=rows))
    p = ROOT / "results" / "u62_runs" / "diagnostico.json"; p.write_text(json.dumps(dict(note=__doc__, hosts=out), indent=1))


if __name__ == "__main__":
    main()
