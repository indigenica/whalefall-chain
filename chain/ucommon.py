"""Unification common layer: paths, vendored code, the single parameter set, JSON I/O, provenance.

Nothing outside the repository is read at run time: frozen inputs live in data/ and reused modules in components/
(see data/MANIFEST.json and components/MANIFEST.json for origin and SHA-256)."""
from pathlib import Path
import hashlib, json, os, platform, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
for p in (ROOT / "components" / "cadena_cierre", ROOT / "components", ROOT / "components" / "c45b_hom", ROOT / "chain"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "logs" / "mpl_config"))
import numpy as np

# the plateau of the chain (ESPEC_U26): 'pole' (order-4 pole inflation, the default since U26) or 'sum' (the 10-level sum of
# ballena 5 up to U25; reproduces every result stored before U26).  Read once, at import, from the environment.
PLATEAU = os.environ.get("WHALEFALL_PLATEAU", "pole")
if PLATEAU not in ("pole", "sum"):
    raise ValueError(f"WHALEFALL_PLATEAU must be 'pole' or 'sum', not {PLATEAU!r}")
INPUTS = ROOT / "data"
RESULTS = ROOT / "results" if PLATEAU == "pole" else ROOT / "results" / "plateau_sum"
FIGS = ROOT / "figures" if PLATEAU == "pole" else ROOT / "figures" / "plateau_sum"
for d in (RESULTS, RESULTS / "cadena_rerun", FIGS, ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)


def load(rel):
    return json.loads((INPUTS / rel).read_text(encoding="utf-8"))


def _conv(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, Path): return str(x)
    raise TypeError(type(x).__name__)


def save(name, data):
    path = RESULTS / name
    path.write_text(json.dumps(data, indent=1, default=_conv), encoding="utf-8")
    return path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def provenance():
    files = sorted(list((ROOT / "chain").glob("*.py")) + [p for p in (ROOT / "components").rglob("*") if p.suffix in (".py", ".c")])
    import scipy
    return dict(utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), plateau=PLATEAU, python=sys.version.split()[0], numpy=np.__version__,
                scipy=scipy.__version__, platform=platform.platform(),
                sources={str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in files})


# ---- the single parameter set of the unified chain (all taken from the closure campaigns, not refitted) ----
THERMAL = dict(sum="cadena_cierre/thermal.json", pole="cadena_cierre/thermal_pole.json")


def benchmark(plateau=None):
    """the parameter set of the chain for a plateau (default: PLATEAU); 'sum' reads thermal.json, 'pole' thermal_pole.json (U26)"""
    plateau = plateau or PLATEAU
    th = load(THERMAL[plateau])["continuous"]
    return dict(plateau=plateau,
        V0=th["V0"],                      # plateau height (reduced Planck units), from the pivot matching
        N_star=th["N_star"], n_s=th["n_s"], r=th["r"], chi_star=th["chi_star"], chi_end=th["chi_end"],
        T_rh_GeV=th["T_rh_GeV"], T_Gamma_GeV=th["T_Gamma_convention_GeV"], Gamma_GeV=th["Gamma_GeV"],
        H_star_GeV=th["H_star_GeV"], F_GeV=th["F_GeV"], f_a_GeV=th["F_GeV"] / 2, mu_GeV=th["mu_GeV"], xi=th["xi"],
        V0_BENCH_C45=1.67e-10,            # the plateau used in the C45/C54 birth scans (ballena 4 benchmark)
        A_C45=1e4,                        # C45 patch size in the C45h/C45j scans
    )
