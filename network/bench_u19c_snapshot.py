"""Cost of one late U19c snapshot on the CPU (lattice.measure with spectra, site_stats, the counter, linear and kinetic number),
on a synthetic RP3 field at N^3; no physics output.  With the evolution on a GPU (network/test_jax_kernel.py --bench N), these
measurements set the time of a network: about NSNAP x (this) + steps x (s/step).  Uses the l_G = 3.318 m_eff table of
U19C_VALIDATION.json for any N (timing only).

    python -B network/bench_u19c_snapshot.py N"""
import json, math, resource, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt                                                          # noqa: E402
import u19_produce as U19                                                     # noqa: E402
import u19b_produce as B                                                      # noqa: E402
import u19c_produce as P                                                      # noqa: E402


def main(N):
    lg = round(3.2 + math.log(N / 512), 3); B.END_FACTOR.setdefault(lg, 3.0); MG = B.set_mass(lg); g = B.geometry(lg, N / 2)
    s, tA, tB = g["s"], g["tau_A"], g["tau_B"]
    table = json.loads((P.OUT / "U19C_VALIDATION.json").read_text())["masses"]["3.318"]["meff"]["RP3"]
    t = time.time(); f, p = Lt.initial(1, N, 6.0, 91); t_ic = time.time() - t
    tau = tA; kb = B.kbins(N); L = U19.L(); L.set_core_late(tA)
    t = time.time(); r = Lt.measure(1, N, f, p, tau, s, spectra=True, tau_f=tA); r.update(U19.site_stats(1, N, f)); t_meas = time.time() - t
    t = time.time(); kw = B.counter_kw(table, tau, tA, MG, s); U19.linear_number(1, N, f, p, s * tau, kb, **kw); t_lin = time.time() - t
    t = time.time(); U19.kinetic_number(r, s * tau); t_kin = time.time() - t
    L.set_core_late(0.0)
    steps = int(round((tB - math.sqrt(math.exp(1.0) * tA)) / P.DTAU))
    print(f"N={N} l_G={lg:g} threads={Lt.lib().nthreads()}: initial {t_ic:.0f}s, measure+spectra {t_meas:.0f}s, linear number {t_lin:.0f}s, "
          f"kinetic {t_kin:.1f}s -> late snapshot {t_meas + t_lin + t_kin:.0f}s; ~{steps} steps, {P.NSNAP} snapshots; "
          f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.1f} GB", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
