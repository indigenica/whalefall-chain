"""The JAX evolution kernel (network/jax_kernel.py) against kern.c.  Not part of `make test` (CI has no JAX).

    python -B network/test_jax_kernel.py            # agreement tests (CPU: bit-identical; GPU: to TOL_GPU)
    python -B network/test_jax_kernel.py --bench N [--lowmem]  # time per leapfrog step at N^3, JAX device vs kern.c (OMP_NUM_THREADS);
                                                               # --lowmem: JAX evolve_lowmem (no force buffer, U19c)
    python -B network/test_jax_kernel.py --memory N  # device memory XLA plans for evolve / evolve_lowmem at N^3 (compile only)

1. 40 steps at N = 32 in an expanding background for RP3 and U(1): plain, with the U19b tilt and late frozen core, and with
   the U10a drain, with evolve and with evolve_lowmem (U19c: no force buffer); once in one z-slab and once in slabs of 4 planes
   with the eigh of the drain in chunks of 1000 matrices (with a remainder: the chunked path of production slabs) (the slab path of force and drain, with its halos).  On CPU the result is bit-identical to kern.c; on a GPU the arithmetic may be contracted differently
   (FMA), and the test asks max |delta| / max |x| <= TOL_GPU for f and p.
2. The same through the backend switch (WHALEFALL_BACKEND=jax in lattice.lib()): a 48^3 U(1) and RP3 network of the U2 driver
   (lattice.run), whose string density at every snapshot must agree with kern.c (identical on CPU, to TOL_XI on a GPU);
   and lattice.lib().evolve_lowmem through the switch (the call of u19c_produce.py at the high mass) against kern.c."""
import math, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import lattice as Lt                                                          # noqa: E402
import u19_produce as U19                                                     # noqa: E402

TOL_GPU, TOL_XI = 1e-5, 0.02
CASES = {"plain": (0.0, 0.0, 0.0, 0.0), "tilt+late core": (0.002, 12.0, 0.0, 0.0), "drain": (0.0, 0.0, 14.0, None)}


def _gpu():
    import jax
    return jax.devices()[0].platform != "cpu"


def agreement(N=32, steps=40):
    import jax
    import jax_kernel as JK
    ok = True
    for slab, eb in ((N, 1 << 15), (4, 1000)):
        JK.SLAB_SITES, JK.EIGH_BATCH = slab * N * N, eb; jax.clear_caches()
        print(f"-- z-slabs of {JK._cz(N)} planes, eigh in chunks of {eb} matrices")
        ok &= _agreement(N, steps)
    JK.SLAB_SITES, JK.EIGH_BATCH = 1 << 22, 1 << 15; jax.clear_caches()
    return ok


def _agreement(N, steps):
    import jax.numpy as jnp
    import jax_kernel as JK
    tol = TOL_GPU if _gpu() else 0.0; ok = True
    for model in (1, 0):
        m = math.sqrt(Lt.M2[model]); s = 1.0 / (m * 16.0); C = Lt.NCOMP[model]
        for case, (eps, tf, dtd, dg) in CASES.items():
            dg = m if dg is None else dg
            f, p = Lt.initial(model, N, 4.0, 61); fc, pc = f.copy(), p.copy(); sc = np.empty_like(fc)
            L = U19.L(); L.set_tilt(eps); L.set_core_late(tf); L.set_damp(dtd, dg)
            L.evolve(model, N, fc, pc, sc, 10.0, 0.25, s, steps)
            L.set_tilt(0.0); L.set_core_late(0.0); L.set_damp(0.0, 0.0)
            for nm, ev in (("evolve", JK.evolve), ("evolve_lowmem", JK.evolve_lowmem)):
                fj, pj = ev(model, jnp.asarray(f.reshape(C, N, N, N)), jnp.asarray(p.reshape(C, N, N, N)), 10.0, 0.25, s, steps,
                            eps_g=eps, tau_f=tf, damp_tau=dtd, damp_g=dg)
                fj, pj = np.asarray(fj).reshape(C, -1), np.asarray(pj).reshape(C, -1)
                df, dp = np.abs(fj - fc).max() / np.abs(fc).max(), np.abs(pj - pc).max() / np.abs(pc).max()
                good = df <= tol and dp <= tol; ok &= good
                print(f"{'PASS' if good else 'FAIL'} {'RP3' if model else 'U(1)':<4} {case:<15} {nm:<13} max|df|/max|f|={df:.1e} max|dp|/max|p|={dp:.1e} (tol {tol:g})")
    return ok


def backend_network(N=48):
    tol = TOL_XI if _gpu() else 0.0; ok = True
    for model in (0, 1):
        xs = {}
        for be in ("c", "jax"):
            os.environ["WHALEFALL_BACKEND"] = be; Lt._L = None
            m = math.sqrt(Lt.M2[model]); r = Lt.run(model, N, N / 2, 61, nsnap=8, spectra_every=100, log=lambda *a: None)
            xs[be] = np.array([x["xi"] for x in r["rows"] if "xi" in x])
        os.environ.pop("WHALEFALL_BACKEND"); Lt._L = None
        d = float(np.max(np.abs(xs["jax"] / xs["c"] - 1)))
        good = d <= tol; ok &= good
        print(f"{'PASS' if good else 'FAIL'} backend switch, {'RP3' if model else 'U(1)'} {N}^3 network (lattice.run): max |xi_jax/xi_c - 1| = {d:.1e} (tol {tol:g})")
    return ok


def drain_batch(N=96):
    """one RP3 drain on a production-sized batch (one slab of N^3 = 884736 matrices, eigh in chunks of EIGH_BATCH) against kern.c
    damp_massive: the call that failed in cuSolver on the A100 before the chunking"""
    import jax.numpy as jnp
    import jax_kernel as JK
    tol = TOL_GPU if _gpu() else 0.0
    f, _ = Lt.initial(1, N, 4.0, 61); p = (0.01 * np.random.default_rng(5).standard_normal(f.shape)).astype(np.float32); pc = p.copy()
    Lt.lib().damp_massive(1, N, f, pc, 0.7)
    pj = np.asarray(JK.damp_massive(1, jnp.asarray(f.reshape(10, N, N, N)), jnp.asarray(p.reshape(10, N, N, N)), 0.7, jnp.float64)).reshape(10, -1)
    d = float(np.abs(pj - pc).max() / np.abs(pc).max()); good = d <= tol
    print(f"{'PASS' if good else 'FAIL'} RP3 drain on one slab of {N}^3 matrices (eigh in chunks of {JK.EIGH_BATCH}): max|dp|/max|p|={d:.1e} (tol {tol:g})")
    return good


def backend_lowmem(N=32, steps=40):
    """lattice.lib().evolve_lowmem through the backend switch (the U19c high-mass call) against kern.c evolve_lowmem"""
    tol = TOL_GPU if _gpu() else 0.0; ok = True
    for model in (1, 0):
        m = math.sqrt(Lt.M2[model]); s = 1.0 / (m * 16.0); out = {}
        for be in ("c", "jax"):
            os.environ["WHALEFALL_BACKEND"] = be; Lt._L = None
            f, p = Lt.initial(model, N, 4.0, 61); L = Lt.lib(); L.set_core_late(12.0)
            L.evolve_lowmem(model, N, f, p, 10.0, 0.25, s, steps); L.set_core_late(0.0); out[be] = (f, p)
        os.environ.pop("WHALEFALL_BACKEND"); Lt._L = None
        (fc, pc), (fj, pj) = out["c"], out["jax"]
        d = max(np.abs(fj - fc).max() / np.abs(fc).max(), np.abs(pj - pc).max() / np.abs(pc).max())
        good = d <= tol; ok &= good
        print(f"{'PASS' if good else 'FAIL'} backend switch, {'RP3' if model else 'U(1)'} evolve_lowmem (U19c): max rel |delta| = {d:.1e} (tol {tol:g})")
    return ok


def bench(N, reps=5, lowmem=False):
    import jax
    import jax.numpy as jnp
    import jax_kernel as JK
    for model in (1, 0):
        m = math.sqrt(Lt.M2[model]); s = 1.0 / (m * N / 2); C = Lt.NCOMP[model]
        f, p = Lt.initial(model, N, 6.0, 61); sc = np.empty_like(f)
        Lt.lib().evolve(model, N, f, p, sc, 20.0, 0.25, s, 1)
        t = time.time(); Lt.lib().evolve(model, N, f, p, sc, 20.25, 0.25, s, reps); tc = (time.time() - t) / reps
        fj, pj = jnp.asarray(f.reshape(C, N, N, N)), jnp.asarray(p.reshape(C, N, N, N))
        ev = JK.evolve_lowmem if lowmem else JK.evolve
        fj, pj = ev(model, fj, pj, 20.0, 0.25, s, 1); jax.block_until_ready(fj)                         # compile
        t = time.time(); fj, pj = ev(model, fj, pj, 20.25, 0.25, s, reps); jax.block_until_ready(fj); tj = (time.time() - t) / reps
        print(f"{'RP3' if model else 'U(1)'} N={N}: kern.c {tc:.3f} s/step ({Lt.lib().nthreads()} threads), JAX{' lowmem' if lowmem else ''} on {jax.devices()[0]} {tj:.3f} s/step "
              f"-> x{tc / tj:.1f}")


def memory(N):
    """device memory that XLA plans for evolve and evolve_lowmem at N^3 (compile only, nothing is allocated)"""
    import jax
    import jax.numpy as jnp
    import jax_kernel as JK
    for model in (1, 0):
        sh = jax.ShapeDtypeStruct((Lt.NCOMP[model], N, N, N), jnp.float32)
        for nm in ("evolve", "evolve_lowmem"):
            m = getattr(JK, nm).lower(model, sh, sh, 20.0, 0.25, 0.01, 5, tau_f=12.0, damp_tau=14.0, damp_g=0.1).compile().memory_analysis()
            print(f"{'RP3' if model else 'U(1)'} N={N} {nm:<13} on {jax.devices()[0]}: arguments {m.argument_size_in_bytes / 1e9:.1f} GB, "
                  f"temporaries {m.temp_size_in_bytes / 1e9:.1f} GB (aliased {m.alias_size_in_bytes / 1e9:.1f} GB)", flush=True)


if __name__ == "__main__":
    if "--memory" in sys.argv:
        memory(int(sys.argv[sys.argv.index("--memory") + 1])); sys.exit(0)
    if "--bench" in sys.argv:
        bench(int(sys.argv[sys.argv.index("--bench") + 1]), lowmem="--lowmem" in sys.argv); sys.exit(0)
    ok = drain_batch(); ok &= agreement(); ok &= backend_network(); ok &= backend_lowmem()
    print("ALL PASSED" if ok else "FAILED"); sys.exit(0 if ok else 1)
