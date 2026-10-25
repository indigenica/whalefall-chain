#!/bin/bash
# Setup script for Claude Code cloud sessions (claude.ai/code -> environment -> Setup script).  Self-contained: paste the
# whole file into the field.  Ubuntu 24.04, root, ~3-4 min (must stay under ~5 min so the environment gets cached).
# Installs the Python deps of the repository and CLASS v3.4.1 (commit 814ce59) with the WHALEFALL quadrature patch.
set -x
pip install --break-system-packages -q numpy scipy matplotlib sympy "cython>=3" || true
if ! python3 -c "import classy" 2>/dev/null; then
  mkdir -p /opt/src && cd /opt/src
  git clone -q https://github.com/lesgourg/class_public.git && cd class_public && git checkout -q 814ce59c380a5cccd0bed1747a384ff7bfc1f9cf || exit 0
  python3 - <<'PY'
from pathlib import Path
h = Path("include/quadrature.h"); s = h.read_text()
s = s.replace("enum ncdm_quadrature_method {qm_auto, qm_Laguerre, qm_trapz_indefinite, qm_trapz};",
              "enum ncdm_quadrature_method {qm_auto, qm_Laguerre, qm_trapz_indefinite, qm_trapz, qm_trapz_log};"); h.write_text(s)
c = Path("tools/quadrature.c"); s = c.read_text()
add = """  case (qm_trapz_log) :
    {
      double lqa, lqb, du;
      if ((qvec != NULL) && (qsiz > 1)) { lqa = log(qvec[0]); lqb = log(qvec[qsiz-1]); }
      else { lqa = log(1e-3); lqb = log(qmax); }
      class_test(N < 2, errmsg, "qm_trapz_log needs at least 2 momentum bins");
      du = (lqb - lqa)/(N - 1.0);
      for (i=0; i<N; i++){
        x[i] = exp(lqa + i*du);
        (*function)(params_for_function,x[i],&y);
        w[i] = y*x[i]*du;
        if ((i==0) || (i==N-1)) w[i] *= 0.5;
      }
    }
    return _SUCCESS_;
  case (qm_trapz_indefinite) :"""
s = s.replace("  case (qm_trapz_indefinite) :", add, 1); c.write_text(s)
PY
  # with build isolation: the Debian setuptools of the image fails with --no-build-isolation (install_layout); 12 Oct 2026
  make -j4 libclass.a > /tmp/class_build.log 2>&1 && CC=gcc pip install --break-system-packages --no-deps . >> /tmp/class_build.log 2>&1 || true
fi
python3 -c "import classy; print('classy', classy.__version__)" || true
exit 0
