"""ctypes wrapper of hom/frw_core.c (C port of hom.frw.run, C45 equations)."""
import ctypes, os, subprocess, math
import numpy as np
from frw import CHI_END  # adapted for unification (was: from hom.frw import CHI_END)

HERE = os.path.dirname(os.path.abspath(__file__))
_L = None


def lib():
    global _L
    if _L is None:
        src = os.path.join(HERE, "frw_core.c"); bdir = os.path.join(HERE, "_build"); os.makedirs(bdir, exist_ok=True)
        so = os.path.join(bdir, "frw_core.so")
        if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
            subprocess.check_call(["gcc", "-O3", "-shared", "-fPIC", "-o", so, src, "-lm"])
        L = ctypes.CDLL(so)
        d = ctypes.c_double
        L.run.restype = ctypes.c_int
        L.run.argtypes = [d, d, d, d, d, d, d, d, ctypes.c_long, d, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), d]
        L.set_deform.restype = None
        L.set_deform.argtypes = [ctypes.c_int, d, d]
        _L = L
    return _L


def set_deform(mode=0, chic=1e30, k=0.0):
    """C45i large-field deformation of the plateau (see frw_core.c); mode 0 restores the undeformed potential"""
    lib().set_deform(int(mode), float(chic), float(k))


OUT = {1: "launched60", 2: "crunch", 0: "cycling", -1: "step_limit", -2: "step_underflow"}


def run_c(a0, al, V0, chi0, chid0, max_cycles=60, efolds_stop=60.0, rtol=1e-11, max_steps=400_000_000, frac_rad=0.0):
    o = np.zeros(12)
    T = max_cycles * 4.0 * a0
    code = lib().run(a0, al, V0, chi0, chid0, T, efolds_stop, rtol, max_steps, CHI_END, o, frac_rad)
    r = dict(outcome=OUT[code], n_bounces=int(o[0]), chi_abs_max=o[1], chi_end=o[2], N_inf=o[3], constraint_max=o[4],
             t_end=o[5], a_min=o[6], z_max=o[7], efolds_total=o[8], n_steps=int(o[9]), E_max=o[10], K_max=o[11])
    r["fertile"] = bool(r["outcome"] == "launched60" and r["chi_abs_max"] > 5 and r["constraint_max"] < 1e-3)
    return r


def c45_state(A, lam, sigma, sign, sep):
    """C45/C54 initial state; sep = al/A^2 (scale-separation parameter; C54 SEP = 1e-3, TIGHT = 0.1)"""
    V0 = lam * 4.0 / (9.0 * A * A)
    return dict(a0=A, al=sep * A * A, V0=V0, chi0=0.0, chid0=sign * sigma * math.sqrt(3 * A / A ** 3))
