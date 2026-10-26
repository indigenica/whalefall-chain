"""Post-processing of a double-null run: centre history (proper time, field,
density, comoving expansion), apparent horizons, inner-horizon mass growth,
extent of the valid (finite) domain."""
import numpy as np
from dn.equations import psi_of_y, local_density


def centre_history(grid, V, alpha, G=1.0, n_off=(2, 3, 4)):
    """Along the centre u = v: proper time dtau = sqrt(2 e2) h, phi_c, rho_c (local),
    kinetic fraction, and the comoving Hubble rate H = (dr/dtau)/r of the
    scalar-comoving congruence, averaged over nodes j - i in n_off (r -> 0 limit)."""
    R, S, P, h = grid.r, grid.s, grid.p, grid.h
    N = grid.N
    out = {k: [] for k in ("k", "u", "tau", "phi", "rho", "Vc", "kin_frac", "H", "Hn", "apsi")}
    sc = np.array([S[k, k] for k in range(N + 1)])
    lap = np.sqrt(2 * np.exp(2 * sc))                     # dtau/du along u = v
    taus = np.concatenate([[0.0], np.cumsum(0.5 * (lap[1:] + lap[:-1]) * h)])
    for k in range(3, N - 1):
        if not (np.isfinite(S[k, k]) and np.isfinite(P[k, k])):
            break
        tau = taus[k]
        Hs, Hns, rhos = [], [], []
        for d in n_off:
            # node (i, j) with j - i = d on the anti-diagonal of the centre: i = k - d/2 (d even) or use (k-1, k-1+d)
            i = k - (d // 2); j = i + d
            if j + 1 > N or i - 1 < 0:
                continue
            blk = [R[i, j], R[i - 1, j], R[i + 1, j], R[i, j - 1], R[i, j + 1]]
            if not np.all(np.isfinite(blk)):
                continue
            ru = (R[i + 1, j] - R[i - 1, j]) / (2 * h); rv = (R[i, j + 1] - R[i, j - 1]) / (2 * h)
            pu = (P[i + 1, j] - P[i - 1, j]) / (2 * h); pv = (P[i, j + 1] - P[i, j - 1]) / (2 * h)
            e2 = np.exp(2 * S[i, j])
            Hns.append((ru + rv) / (R[i, j] * np.sqrt(2 * e2)))
            if pu * pv > 0:                       # timelike gradient: comoving frame exists
                nrm = np.sqrt(2 * pu * pv * e2)
                Hs.append((pv * ru + pu * rv) / (R[i, j] * nrm) * np.sign(pu + pv))
            rhos.append(local_density(S[i, j], pu, pv, P[i, j], V))
        if not Hns:
            break
        # intrinsic centre density: phi_u = phi_v = (1/2) d phi_c / du along u = v (second order),
        # rho_c = (phi_u^2 + phi_v^2)/(2 e2) + V = phi_v^2/e2 + V
        if not (np.isfinite(P[k + 1, k + 1]) and np.isfinite(P[k - 1, k - 1])):
            break
        pv_c = (P[k + 1, k + 1] - P[k - 1, k - 1]) / (4 * h)
        Vc = float(V(P[k, k]))
        rho = float(pv_c * pv_c / np.exp(2 * S[k, k]) + Vc)
        out["k"].append(k); out["u"].append(k * h); out["tau"].append(tau)
        out["phi"].append(P[k, k]); out["rho"].append(rho); out["Vc"].append(Vc)
        out["kin_frac"].append(1 - Vc / rho if rho > 0 else np.nan)
        out["H"].append(float(np.mean(Hs)) if Hs else np.nan); out["Hn"].append(float(np.mean(Hns)))
        out["apsi"].append(alpha * psi_of_y(8 * np.pi * G * rho / 3, alpha))
    return {k: np.array(v) for k, v in out.items()}


def horizon_map(grid, alpha, G=1.0):
    """Node classes from one-sided differences: 0 normal (ru<0<rv), 1 trapped (ru<0, rv<0),
    2 anti-trapped (ru>0, rv>0), 3 'reversed' (ru>0>rv).  Also f = 1 - r^2 psi(M)."""
    R, M, h = grid.r, grid.M, grid.h
    ru = np.full_like(R, np.nan); rv = np.full_like(R, np.nan)
    ru[1:, :] = (R[1:, :] - R[:-1, :]) / h
    rv[:, 1:] = (R[:, 1:] - R[:, :-1]) / h
    cls = np.full(R.shape, -1)
    ok = np.isfinite(ru) & np.isfinite(rv) & (R > 0)
    cls[ok & (ru < 0) & (rv > 0)] = 0
    cls[ok & (ru < 0) & (rv <= 0)] = 1
    cls[ok & (ru >= 0) & (rv > 0)] = 2
    cls[ok & (ru >= 0) & (rv <= 0)] = 3
    with np.errstate(invalid="ignore", divide="ignore"):
        f = 1 - R ** 2 * psi_of_y(2 * G * M / R ** 3, alpha)
    return cls, f, ru, rv


def summary(grid, V, alpha, G=1.0, M_adm=None):
    cls, f, ru, rv = horizon_map(grid, alpha, G)
    R, M, S = grid.r, grid.M, grid.s
    fin = np.isfinite(R)
    trapped = cls == 1
    anti = cls == 2
    out = dict(N=grid.N, h=grid.h, ok_until_uv=grid.ok_until * grid.h if grid.ok_until else None,
               n_failed_seeds=int(getattr(grid, "n_failed", 0)),
               frac_valid=float(fin.sum() / np.triu(np.ones_like(R)).sum()))
    ii, jj = np.nonzero(trapped)
    if ii.size:
        k = np.argmin(ii + jj)
        out.update(trapped=True, first_trapped_u=float(ii[k] * grid.h), first_trapped_v=float(jj[k] * grid.h),
                   first_trapped_r=float(R[ii[k], jj[k]]), r_trapped_max=float(np.nanmax(R[trapped])),
                   r_trapped_min=float(np.nanmin(R[trapped])))
    else:
        out.update(trapped=False)
    out["anti_trapped_nodes"] = int(anti.sum())
    if anti.sum():
        ia, ja = np.nonzero(anti)
        sel = (ja - ia) > 4
        out["anti_trapped_offcentre_nodes"] = int(sel.sum())
        if sel.any():
            out["anti_trapped_r_range"] = [float(np.nanmin(R[ia[sel], ja[sel]])), float(np.nanmax(R[ia[sel], ja[sel]]))]
    out["M_max_interior"] = float(np.nanmax(np.where(R < 2, M, np.nan)))
    out["M_min_interior"] = float(np.nanmin(np.where((R < 2) & (R > 0), M, np.nan)))
    if M_adm is not None:
        out["mass_inflation_ratio"] = out["M_max_interior"] / M_adm
    out["sigma_min"] = float(np.nanmin(S))
    return out
