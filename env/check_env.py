"""Environment check for WHALEFALL (WSL):  python env/check_env.py
numpy/scipy/matplotlib/sympy, gcc + OpenMP, classy (a Planck-like run and a file-based ncdm distribution), RAM/cores plan."""
import importlib, json, os, shutil, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ok = True


def say(name, good, info=""):
    global ok
    ok &= bool(good); print(f"[{'OK' if good else 'FALLA'}] {name} {info}", flush=True)


print(f"python {sys.version.split()[0]} en {sys.executable}  (conda env: {os.environ.get('CONDA_DEFAULT_ENV', '?')})", flush=True)
for m in ("numpy", "scipy", "matplotlib", "sympy", "Cython"):
    try:
        mod = importlib.import_module(m); say(m, True, getattr(mod, "__version__", ""))
    except Exception as e:  # noqa: BLE001
        say(m, False, str(e))

try:  # informativo: U4 no usa la GPU
    import jax; print(f"[info] jax {jax.__version__}: {jax.devices()}", flush=True)
except Exception as e:  # noqa: BLE001
    print(f"[info] jax no disponible ({type(e).__name__})", flush=True)
gcc = shutil.which("gcc"); say("gcc", gcc is not None, gcc or "")
if gcc:
    with tempfile.TemporaryDirectory() as d:
        c = Path(d) / "t.c"; c.write_text("#include <omp.h>\n#include <stdio.h>\nint main(){printf(\"%d\",omp_get_max_threads());return 0;}\n")
        r = subprocess.run([gcc, "-fopenmp", str(c), "-o", str(Path(d) / "t")], capture_output=True, text=True)
        th = subprocess.run([str(Path(d) / "t")], capture_output=True, text=True).stdout if r.returncode == 0 else ""
        say("gcc -fopenmp", r.returncode == 0, f"threads={th} {r.stderr[:200]}")

try:
    import classy
    say("classy importa", True, f"{getattr(classy, '__version__', '?')} ({classy.__file__})")
    c = classy.Class(); t = time.time()
    c.set({"h": 0.674, "omega_b": 0.0224, "omega_cdm": 0.12, "A_s": 2.1e-9, "n_s": 0.9649, "output": "mPk", "P_k_max_h/Mpc": 10.0})
    c.compute(); p = c.pk(0.1, 0.0); c.struct_cleanup(); c.empty()
    say("classy compute (LCDM)", p > 0, f"P(0.1/Mpc, z=0) = {p:.4g} Mpc^3 in {time.time()-t:.1f}s")
    import numpy as np
    with tempfile.TemporaryDirectory() as d:
        q = np.geomspace(1e-3, 40, 800); f = Path(d) / "fd.dat"; np.savetxt(f, np.column_stack([q, 1 / (np.exp(q) + 1)]))
        base = {"h": 0.674, "omega_b": 0.0224, "omega_cdm": 0.0, "A_s": 2.1e-9, "n_s": 0.9649, "N_ur": 3.044, "output": "mPk",
                "P_k_max_h/Mpc": 30.0, "N_ncdm": 1, "m_ncdm": 3500.0, "T_ncdm": 0.1477 * (4 / 11) ** (1 / 3), "omega_ncdm": 0.12,
                "ncdm_fluid_approximation": 3, "z_pk": 3.0}
        res = []
        for extra in ({}, {"use_ncdm_psd_files": 1, "ncdm_psd_filenames": str(f)}):
            c = classy.Class(); c.set(dict(base, **extra)); t = time.time(); c.compute()
            res.append([c.pk(k * 0.674, 3.0) for k in (1.0, 10.0, 30.0)]); c.struct_cleanup(); c.empty()
        rel = max(abs(a / b - 1) for a, b in zip(res[1], res[0]))
        say("classy ncdm con fichero de distribución (WDM 3.5 keV)", rel < 5e-3, f"max rel = {rel:.2e} ({time.time()-t:.0f}s por ejecución)")
        # U4a needs the WHALEFALL quadrature patch (env/class_qm_trapz_log.patch): ncdm_quadrature_strategy = 4.  An unpatched
        # CLASS does not complain but leaves the momentum weights unset, so compare the background against the default sampling.
        bg = []
        for extra in ({}, {"ncdm_quadrature_strategy": 4, "ncdm_N_momentum_bins": 45}):
            bgp = {k: v for k, v in base.items() if k not in ("output", "P_k_max_h/Mpc", "z_pk")}   # background only
            c = classy.Class(); c.set(dict(bgp, use_ncdm_psd_files=1, ncdm_psd_filenames=str(f), **extra)); c.compute()
            b = c.get_background(); bg.append((np.array(b["z"]), np.array(b["(.)rho_ncdm[0]"]), np.array(b["(.)p_ncdm[0]"])))
            c.struct_cleanup(); c.empty()
        zz = np.geomspace(1e-2, 1e9, 60)
        rr = [np.interp(np.log(zz), np.log(z[::-1] + 1e-30), r[::-1]) for z, r, _ in bg]
        pp = [np.interp(np.log(zz), np.log(z[::-1] + 1e-30), p[::-1]) for z, _, p in bg]
        relb = float(max(np.max(np.abs(rr[1] / rr[0] - 1)), np.max(np.abs(pp[1] / pp[0] - 1))))
        say("CLASS con el parche de cuadratura logarítmica (strategy 4)", relb < 2e-3,
            f"max rel rho/p ncdm = {relb:.1e}" + ("" if relb < 2e-3 else "  -> aplicad env/class_qm_trapz_log.patch y reinstalad (INSTALAR_WSL.md §3)"))
except Exception as e:  # noqa: BLE001
    say("classy", False, f"{type(e).__name__}: {e}")

sys.path.insert(0, str(ROOT / "network"))
from u4b_plan import plan  # noqa: E402
print(json.dumps(plan(), indent=1))
print("\nENTORNO LISTO" if ok else "\nHAY FALLOS: ver env/INSTALAR_WSL.md (Problemas frecuentes)")
sys.exit(0 if ok else 1)
