#!/bin/bash
# U15b (docs/preregistros/ESPEC_U15b.md): Planck 2018 lite + BAO + Lyman-ABGD chains for the WHALEFALL warm tail, in WSL.
#
# Usage (from the repository, in the lighthouse-sage environment with the patched CLASS of env/INSTALAR_WSL.md):
#   bash lya/u15b/run_wsl.sh install     # 1) MontePython + Planck likelihood code and data + this likelihood   (~20-40 min)
#   bash lya/u15b/run_wsl.sh tables      # 2) T^2 tables with CLASS (51 runs; ~1 h with 40 workers)
#   bash lya/u15b/run_wsl.sh test        # 3) 30 steps of the 1 keV validation case, to check everything runs   (~10 min)
#   bash lya/u15b/run_wsl.sh chains      # 4) the active cases (CASES, default: the four U15d cases) x 4 MPI chains, in the background
#   bash lya/u15b/run_wsl.sh status [x]  #    state and convergence (R-1) of every case folder found (U15b, U15c, ...)
#   bash lya/u15b/run_wsl.sh diagnose    #    per chain file: points, update marks, restart segments, copies and forks (Nota 9)
#   bash lya/u15b/run_wsl.sh pause [x]   #    freeze every running chain (SIGSTOP: no CPU, memory kept); 'resume [x]' continues them
#   bash lya/u15b/run_wsl.sh stop [x]    #    end every running chain; 'chains [x]' later continues each one from its last point (-r)
#   bash lya/u15b/run_wsl.sh boost [x]   #    stop, refresh the proposal covmat from the chains, oversample the Lya nuisance block x10
#                                        #    (cases: CASES plus every case running, e.g. 'boost U15d')
#                                        #    (data.over_sampling, tested first) and continue from the last points (ESPEC_U15b Nota 8)
#   bash lya/u15b/run_wsl.sh fixrestart  #    move aside old chain files that a restart already copied (also done automatically)
#                                        #    [x]: optional part of the case names, e.g. x30 or val_thermal_1keV (default: all)
#   bash lya/u15b/run_wsl.sh analyze     # 5) posteriors, verdict -> lya/results/u15b/U15B_RESULT.json
#
# Paths can be changed with environment variables (defaults below).
set -e
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="${U15B_WORK:-$HOME/u15b}"                         # MontePython, Planck, chains (outside the repository: several GB)
CLASS_DIR="${CLASS_DIR:-$HOME/src/class_public}"       # the patched CLASS used by U4-U12
NPROC="${NPROC:-40}"
CHAINS_PER_CASE="${CHAINS_PER_CASE:-4}"
STEPS="${STEPS:-100000}"
MP="$WORK/montepython_public"
PLANCK="${PLANCK:-$WORK/planck/baseline/plc_3.0}"            # Planck 2018 baseline data
CLIK_DIR="${CLIK_DIR:-$WORK/planck/code/plc_3.0/plc-3.1}"     # Planck likelihood code (clik)
[ -f "$WORK/paths.env" ] && source "$WORK/paths.env"             # written by install: where the tarballs really unpacked
# Cases launched by 'chains' / 'boost' when CASES is not given: the ones still active.  U15b (validation) and U15c are finished
# and analyzed (INFORME_U15b_U15c.md), and the four physics rows U15b_B2_* are stopped (ESPEC_U15b Nota 7): relaunch any of them
# only with an explicit CASES=...  status, pause, resume and stop do not use this list: they act on every case folder found in
# $WORK/chains, or on every chain actually running (Nota 9).
CASES="${CASES:-U15d_F1_Y1.20 U15d_F1_Y1.70 U15d_F0.85_Y1.98 U15d_F0.7_Y1.98}"

install() {
  mkdir -p "$WORK" && cd "$WORK"
  # MontePython, the commit whose Lya_abgd was used in U15
  [ -d "$MP" ] || git clone https://github.com/brinckmann/montepython_public.git
  cd "$MP" && git checkout -q 5b365b093982099b926a0c0c169697c9dc6cbbc4 && cd "$WORK"
  python "$REPO/lya/u15b/patch_montepython_py313.py" "$MP"      # Python 3.13: likelihood import (PEP 667)
  python "$REPO/lya/u15b/patch_montepython_restart.py" "$MP"    # restarts: lock before truncating the new chain file (Nota 10)
  mkdir -p "$CLASS_DIR/python/build"                              # MontePython looks here; classy is the pip-installed one
  python -c "import mpi4py" 2>/dev/null || { echo "installing mpi4py (needs MPI: sudo apt install -y openmpi-bin libopenmpi-dev)"; conda install -y -c conda-forge mpi4py || pip install mpi4py; }
  python -c "import cython" 2>/dev/null || pip install cython      # the Planck wrapper (clik) is built with cython
  # Planck 2018 likelihood code (PLC 3.1) and baseline data, from the Planck Legacy Archive
  mkdir -p planck && cd planck
  [ -f data.tar.gz ] || wget -O data.tar.gz "https://pla.esac.esa.int/pla-sl/data-action?COSMOLOGY.COSMOLOGY_OID=151902"
  # Likelihood code: clik 16 (benabed/clik, the maintained successor of PLC 3.x, same .clik files, Python-3 wrapper)
  [ -d clik16 ] || git clone https://github.com/benabed/clik.git clik16
  ( cd clik16 && git checkout -q 9c6a194aece4175a3c499e4e8a2118b8338bc395 )
  [ -f .data_extracted ] || { tar xzf data.tar.gz && touch .data_extracted; }
  # locate what the tarballs contain (their folder names change between PLA releases)
  local waf="$WORK/planck/clik16/waf"
  local lite; lite=$(find "$WORK/planck" -maxdepth 8 -type d -name "plik_lite_v22_TTTEEE.clik" | head -1)
  if [ -z "$waf" ] || [ -z "$lite" ]; then
    echo "Could not find the Planck code (waf) or the plik_lite data. Contents of the tarballs:"
    ls "$WORK/planck/clik16" | head; echo ...; tar tzf data.tar.gz | grep -E "\.clik/?$" | head -20
    exit 1
  fi
  CLIK="$(dirname "$waf")"
  local planck_data; planck_data="$(cd "$lite/../../.." && pwd)"         # .../baseline/plc_3.0
  for f in low_l/commander/commander_dx12_v3_2_29.clik low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik; do
    [ -e "$planck_data/$f" ] || { echo "missing Planck data file $planck_data/$f"; exit 1; }
  done
  # The PLC 3.01 build system (waf 2.0.17) imports the module 'imp', removed in Python 3.12, and uses distutils:
  # a 10-line shim for 'imp' and setuptools' distutils are enough for Python 3.13.
  mkdir -p "$WORK/pyshim"
  cat > "$WORK/pyshim/imp.py" <<'SHIM'
"""Minimal stand-in for the 'imp' module removed in Python 3.12 (only what waf 2.0.x uses)."""
import importlib, importlib.util, types
def new_module(name): return types.ModuleType(name)
def reload(module): return importlib.reload(module)
def load_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
SHIM
  python -c "import setuptools" 2>/dev/null || pip install setuptools
  # Build with the system gcc/gfortran: the conda compiler wrappers of the environment (CC=x86_64-conda-linux-gnu-cc and
  # their CFLAGS/LDFLAGS) fail waf's test compile.
  [ -x /usr/bin/gcc ] || { echo "system gcc missing: sudo apt install -y build-essential"; exit 1; }
  [ -x /usr/bin/gfortran ] || { echo "system gfortran missing (clik needs it): sudo apt install -y gfortran"; exit 1; }
  # clik links against system lapack/blas and cfitsio
  for lib in liblapack.so libblas.so libcfitsio.so; do
    ldconfig -p | grep -q "$lib" || { echo "system $lib missing: sudo apt install -y liblapack-dev libblas-dev libcfitsio-dev"; exit 1; }
  done
  local BUILD_ENV=(env -u CFLAGS -u CPPFLAGS -u CXXFLAGS -u LDFLAGS -u FFLAGS -u FCFLAGS -u LDSHARED -u CC -u CXX -u FC -u F77 -u F90
                   CC=/usr/bin/gcc FC=/usr/bin/gfortran PYTHONPATH="$WORK/pyshim${PYTHONPATH:+:$PYTHONPATH}" SETUPTOOLS_USE_DISTUTILS=local)
  cd "$CLIK"
  if [ ! -f bin/clik_profile.sh ]; then
    rm -rf build
    if ! { "${BUILD_ENV[@]}" python ./waf configure --install_all_deps && "${BUILD_ENV[@]}" python ./waf install; }; then
      echo; echo "--- last lines of build/config.log ---"; tail -40 build/config.log 2>/dev/null
      echo; echo "The Planck likelihood did not build.  Paste the last lines above, or use the Cobaya route of lya/u15b/README.md."
      exit 1
    fi
  fi
  [ -f bin/clik_profile.sh ] || { echo "waf finished but bin/clik_profile.sh is missing"; exit 1; }
  python -c "import Cython" 2>/dev/null || pip install cython
  if ! ( source bin/clik_profile.sh && python -c "import clik.lkl" 2>/dev/null ); then
    echo "waf did not build the Python wrapper; building it with Cython"
    "${BUILD_ENV[@]}" LDSHARED="/usr/bin/gcc -shared" python "$REPO/lya/u15b/build_clik_wrapper.py" "$CLIK" || { echo "wrapper build failed"; exit 1; }
  fi
  ( source bin/clik_profile.sh && python -c "import clik; c = clik.clik('$planck_data/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik'); print('clik OK: plik_lite loads, lmax', c.get_lmax())" ) \
    || { echo "clik does not load plik_lite"; exit 1; }
  cd "$WORK"
  printf 'CLIK_DIR="%s"\nPLANCK="%s"\n' "$CLIK" "$planck_data" > "$WORK/paths.env"
  echo "Planck code: $CLIK"; echo "Planck data: $planck_data"
  # MontePython configuration
  cat > "$MP/default.conf" <<CONF
path['cosmo'] = '$CLASS_DIR'
path['clik'] = '$CLIK'
CONF
  # this likelihood
  rm -rf "$MP/montepython/likelihoods/Lya_abgd_whale"
  cp -r "$REPO/lya/u15b/Lya_abgd_whale" "$MP/montepython/likelihoods/"
  echo; echo "Done.  Before every session:  source $CLIK/bin/clik_profile.sh"
}

clik_env() { source "$CLIK_DIR/bin/clik_profile.sh"; }

make_param() {   # $1 = case name
  local tab="$REPO/lya/u15b/tables/$1.npz"
  [ -f "$tab" ] || { echo "missing table $tab: run 'tables' first"; exit 1; }
  local smax; smax=$(python -c "import numpy as np; print(float(np.load('$tab')['s'].max()))")
  mkdir -p "$WORK/params"
  sed -e "s#@TABLE@#$tab#" -e "s#@SMAX@#$smax#" -e "s#@PLANCK@#$PLANCK#" "$REPO/lya/u15b/params/template.param" > "$WORK/params/$1.param"
  echo "$WORK/params/$1.param"
}

tables() { cd "$REPO" && python -B lya/u15b/prepare_tables.py --workers "$NPROC"; }

test_run() {
  clik_env; cd "$MP"
  local p; p=$(make_param U15b_val_thermal_1keV)
  rm -rf "$WORK/chains/test"
  python montepython/MontePython.py run -p "$p" -o "$WORK/chains/test" -N 30 -c covmat/base2018TTTEEE_lite.covmat -j fast
  echo "If this printed accepted steps without errors, launch: bash lya/u15b/run_wsl.sh chains"
}

sel() { local c; for c in $CASES; do [[ "$c" == *"${1:-}"* ]] && echo "$c"; done; }   # cases matching the optional filter
sel_or_running() {   # boost: the CASES list plus every case with chains running, matching the optional filter (stopped rows stay out)
  { sel "${1:-}"; running_cases "${1:-}"; } | sort -u; }
found() {   # every case folder in $WORK/chains (with a log.param), optionally filtered
  local d c; for d in "$WORK"/chains/*/; do c=$(basename "$d"); [ "$c" = test ] && continue
    [ -f "$d/log.param" ] && [[ "$c" == *"${1:-}"* ]] && echo "$c"; done
}
running_cases() {   # every case with a MontePython chain process alive (running or paused), optionally filtered
  pgrep -af "^[^ ]*(python[0-9.]*|mpirun|mpiexec(.hydra)?) .*MontePython.py run -o $WORK/chains/" | sed -n "s#.*MontePython.py run -o $WORK/chains/\([^ ]*\).*#\1#p" | sort -u |
    while read -r c; do [[ "$c" == *"${1:-}"* ]] && echo "$c"; done
}
running() { pgrep -f "^[^ ]*(python[0-9.]*|mpirun|mpiexec(.hydra)?) .*MontePython.py run -o $WORK/chains/$1( |$)" > /dev/null; }

restart_case() {   # continue the chains of case $1 from their last points (MontePython -r copies each old chain into a new file)
  python "$REPO/lya/u15b/patch_montepython_restart.py" "$MP" > /dev/null   # Nota 10: no truncation race between the new chains
  local d="$WORK/chains/$1" last base files n first
  last=$(ls -t "$d"/*__*.txt | head -1); base="${last%%__*}"
  files=$(ls "$base"__*.txt | sort -V); n=$(echo "$files" | wc -l); first=$(echo "$files" | head -1)
  nohup mpirun -np "$n" python montepython/MontePython.py run -o "$d" -N "$STEPS" -r "$first" \
        -c covmat/base2018TTTEEE_lite.covmat -j fast --update 50 --superupdate 20 >> "$WORK/logs/$1.log" 2>&1 &
  # in the background: once the new chains have copied the old ones and stepped, the old files go to restarted_from/
  nohup bash "$REPO/lya/u15b/move_restarted.sh" "$d" "$n" >> "$WORK/logs/$1.log" 2>&1 &
  echo "continuing $1 from $(basename "$base") ($n chains; the old files will be moved to $d/restarted_from/)"
}

chains() {
  clik_env; cd "$MP"; mkdir -p "$WORK/chains" "$WORK/logs"
  local n; n=$(echo $CASES | wc -w)
  # one CLASS thread per chain (CLASS is built with OpenMP: without this, 24 chains x 44 threads oversubscribe the machine)
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(( NPROC / (n * CHAINS_PER_CASE) > 0 ? NPROC / (n * CHAINS_PER_CASE) : 1 ))}"
  export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"          # Intel MPI (conda) on a single machine; ignored by Open MPI
  echo "MPI check:"; mpirun -np 2 python -c "from mpi4py import MPI; c = MPI.COMM_WORLD; print('  rank', c.Get_rank(), 'of', c.Get_size())" \
    || { echo "mpirun + mpi4py do not work: see lya/u15b/README.md"; exit 1; }
  for c in $(sel "${1:-}"); do
    if running "$c"; then echo "$c already running (if paused: resume)"; continue; fi
    if ls "$WORK/chains/$c"/*__*.txt > /dev/null 2>&1; then restart_case "$c"; continue; fi
    local p; p=$(make_param "$c")
    # first call creates the folder and log.param; then the MPI chains
    [ -f "$WORK/chains/$c/log.param" ] || python montepython/MontePython.py run -p "$p" -o "$WORK/chains/$c" -f 0 > /dev/null
    python "$REPO/lya/u15b/patch_montepython_restart.py" "$MP" > /dev/null
    nohup mpirun -np "$CHAINS_PER_CASE" python montepython/MontePython.py run -o "$WORK/chains/$c" -N "$STEPS" \
          -c covmat/base2018TTTEEE_lite.covmat -j fast --update 50 --superupdate 20 >> "$WORK/logs/$c.log" 2>&1 &
    echo "launched $c (log: $WORK/logs/$c.log)"
  done
  echo "OMP_NUM_THREADS=$OMP_NUM_THREADS.  Keep this WSL window open (closing every WSL window can stop the distribution)."
}

status() {
  # convergence from analyze.py --check: the same R-1 (every parameter, after the same cuts) that decides the verdict.
  # (MontePython's info labels only the first parameter's R-1 line, so reading it showed omega_b alone: ESPEC_U15b, Nota 5.)
  local c
  for c in $(found "${1:-}"); do
    local d="$WORK/chains/$c" state="stopped" np=0
    if running "$c"; then
      np=$(pgrep -fc "^[^ ]*python[0-9.]* .*MontePython.py run -o $d( |$)" || true)
      local pids; pids=$(pgrep -f "^[^ ]*(python[0-9.]*|mpirun|mpiexec(.hydra)?) .*MontePython.py run -o $d( |$)" || true)
      state="running ($np processes)"
      if [ -n "$pids" ] && [ -z "$(ps -o stat= -p "$(echo $pids | tr ' ' ',')" | grep -v T)" ]; then state="paused ($np processes)"; fi
    fi
    local dup; dup=$(python "$REPO/lya/u15b/restarted.py" "$d" | head -1 | awk '{print $2}')
    echo "$c: $state, $(ls "$d"/*__*.txt 2>/dev/null | wc -l) chain files"
    if [ "${dup:-0}" -gt 0 ]; then echo "   $dup old file(s) already copied into restarted chains (not counted; 'fixrestart' moves them)"; fi
  done
  cd "$REPO" && python -B lya/u15b/analyze.py "$WORK/chains" --check
}

signal_cases() {   # $1 = signal, $2 = filter: every chain actually running (any case, any batch), not just the CASES list (Nota 9)
  local c k=0; for c in $(running_cases "${2:-}"); do
    pkill "-$1" -f "^[^ ]*(python[0-9.]*|mpirun|mpiexec(.hydra)?) .*MontePython.py run -o $WORK/chains/$c( |$)" && { echo "$1 -> $c"; k=1; }; done
  [ $k = 1 ] || echo "no matching chains running"
}
pause() { signal_cases STOP "${1:-}"; echo "Paused: no CPU used, memory kept.  'resume' continues them; a WSL restart loses them (then use 'chains')."; }
resume() { signal_cases CONT "${1:-}"; }
stop() { signal_cases CONT "${1:-}" > /dev/null; signal_cases TERM "${1:-}"; }

fixrestart() { python "$REPO/lya/u15b/restarted.py" "$WORK/chains" --move; }

OVERSAMPLE="${OVERSAMPLE:-10}"                          # Lya nuisance block oversampling (ESPEC_U15b Nota 8)
set_oversampling() {   # $1 = log.param, $2 = list: replace or append data.over_sampling
  if grep -q "^data.over_sampling" "$1"; then sed -i "s#^data.over_sampling.*#data.over_sampling=$2#" "$1"; else echo "data.over_sampling=$2" >> "$1"; fi
}
boost() {   # ESPEC_U15b Nota 8: same posterior, better sampler (proposal covmat from the chains + Lya nuisances oversampled)
  clik_env; cd "$MP"; mkdir -p "$WORK/tests" "$WORK/logs"
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"           # CLASS threads per chain: 16 chains x 2 leave room for U19b (Nota 8)
  export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"
  local c f="$OVERSAMPLE" stamp; stamp=$(date +%Y%m%d_%H%M%S)
  local todo; todo=$(sel_or_running "${1:-}")
  [ -n "$todo" ] || { echo "boost: no case matches '${1:-}' (neither in CASES nor running)"; return 1; }
  for c in $todo; do
    local d="$WORK/chains/$c"; [ -f "$d/log.param" ] || { echo "$c: no chains yet"; continue; }
    if running "$c"; then stop "$c" > /dev/null; local w=0; while running "$c" && [ $w -lt 120 ]; do sleep 5; w=$((w + 5)); done; fi
    running "$c" && { echo "$c: still running, not touched"; continue; }
    # 1) the block layout: Planck lite + low-l TT + low-l EE share A_planck, so MontePython makes 5 blocks [cosmo, A_planck, -, -, Lya];
    #    a 1-step run in a scratch folder checks the length before the real chains are touched
    local os="" cand t
    for cand in "[1, 1, 1, 1, $f]" "[1, 1, $f]" "[1, 1, 1, $f]"; do
      t="$WORK/tests/boost_${c}_$stamp"; mkdir -p "$t"; cp "$d/log.param" "$t/log.param"; set_oversampling "$t/log.param" "$cand"
      if timeout 900 python montepython/MontePython.py run -o "$t" -N 1 -f 0 -c covmat/base2018TTTEEE_lite.covmat > "$t/out.txt" 2>&1 \
         && ! grep -q "over_sampling field" "$t/out.txt"; then os="$cand"; break; fi
      mv "$t" "${t}_rejected_$(echo "$cand" | tr -dc '0-9_')"
    done
    [ -n "$os" ] || { echo "$c: no over_sampling layout accepted (see $WORK/tests/boost_${c}_$stamp*): chains left stopped; 'chains' restarts them unchanged"; continue; }
    # 2) proposal covmat from the chains themselves (the old one is kept next to it)
    local base; base=$(basename "$d"); [ -f "$d/$base.covmat" ] && cp "$d/$base.covmat" "$d/$base.covmat.before_nota8_$stamp"
    cp "$d/log.param" "$d/log.param.before_nota8_$stamp"
    python montepython/MontePython.py info "$d" --want-covmat --minimal --noplot --keep-fraction 0.5 > "$WORK/logs/$c.boost_info.log" 2>&1 \
      && echo "$c: proposal covmat refreshed from the chains" || echo "$c: covmat refresh failed (see $WORK/logs/$c.boost_info.log); keeping the old one"
    # 3) oversampling, then continue every chain from its last point
    set_oversampling "$d/log.param" "$os"; echo "$c: data.over_sampling=$os"
    restart_case "$c"
  done
}

analyze() { cd "$REPO" && python -B lya/u15b/analyze.py "$WORK/chains"; }

case "${1:-}" in
  install) install ;;
  tables) tables ;;
  test) test_run ;;
  chains) chains "${2:-}" ;;
  status) status "${2:-}" ;;
  diagnose) python "$REPO/lya/u15b/chainread.py" "$WORK/chains" ;;
  pause) pause "${2:-}" ;;
  resume) resume "${2:-}" ;;
  stop) stop "${2:-}" ;;
  boost) boost "${2:-}" ;;
  fixrestart) fixrestart ;;
  analyze) analyze ;;
  *) sed -n 2,16p "$0" ;;
esac
