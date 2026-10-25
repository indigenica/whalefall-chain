"""Build the Python wrapper of clik (clik.lkl, clik.lkl_lensing) against an already installed libclik, with Cython 3.

The waf install builds libclik but skips the Cython extensions with Cython 3 (its version check reads 3.0 as older than 0.28) ("Cannot use clik wrapper (cause = 'No module
named clik.lkl')").  This compiles them from the clik sources and puts them next to the installed clik package.

    python lya/u15b/build_clik_wrapper.py /path/to/clik16"""
import glob, os, sys
from pathlib import Path

clik = Path(sys.argv[1]).resolve()
import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

# where waf installed the clik package (the folder that clik_profile.sh puts on PYTHONPATH): the one with clik/__init__.py
cands = [Path(d) for d in glob.glob(str(clik / "lib" / "python*" / "site-packages")) if (Path(d) / "clik" / "__init__.py").exists()]
if not cands:
    sys.exit(f"no installed clik package under {clik / 'lib'}: run the waf install first")
site = cands[0]
src = clik / "src"
# every source folder that holds headers (clik.h includes io.h and errorlist.h from src/minipmc, etc.)
hdirs = sorted({str(h.parent) for h in src.rglob("*.h")})
opts = dict(include_dirs=[np.get_include(), str(clik / "include"), str(src), str(src / "python" / "clik")] + hdirs,
            define_macros=[("HAS_LAPACK", None), ("LAPACK_CLIK", None), ("NOHEALPIX", None), ("CLIK_LENSING", None)],
            libraries=["clik"], library_dirs=[str(clik / "lib")], runtime_library_dirs=[str(clik / "lib")])
exts = [Extension("clik.lkl", [str(src / "python" / "clik" / "lkl.pyx")], **opts),
        Extension("clik.lkl_lensing", [str(src / "python" / "clik" / "lkl_lensing.pyx")], **opts)]
build_tmp = clik / "build" / "wrapper"
os.chdir(clik)
setup(name="clik_wrapper", ext_modules=cythonize(exts, compiler_directives={"language_level": "3"}, build_dir=str(build_tmp)),
      script_args=["build_ext", "--build-lib", str(site), "--build-temp", str(build_tmp)])
print("clik wrapper built into", site)
