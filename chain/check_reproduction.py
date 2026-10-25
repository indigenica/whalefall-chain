"""Compare freshly generated chain results with the committed ones (U26 audit): every number must agree to rtol = 1e-6 (or
|diff| <= 1e-12 for round-off diagnostics), every other value exactly; provenance and run times are ignored.
    python -B chain/check_reproduction.py results            # after `make chain` (pole)
    python -B chain/check_reproduction.py results/plateau_sum    # after `WHALEFALL_PLATEAU=sum make chain`"""
import json, math, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = ("s1_birth.json", "s2_pocket.json", "s3_sectors.json", "s4_dark_checks.json", "s5_closure.json", "s6_measured_closure.json",
         "u17_ksw.json", "CONTRATO_CIERRE.json")
SKIP = {"provenance", "runtime_s", "pivot_runtime_s", "grid_runtime_s", "seconds", "utc"}


def diff(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in SKIP: continue
            if k not in a or k not in b: out.append(f"{path}/{k}: only in {'committed' if k in b else 'new'}")
            else: diff(a[k], b[k], f"{path}/{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append(f"{path}: length {len(a)} != {len(b)}"); return
        for i, (x, y) in enumerate(zip(a, b)): diff(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if not (math.isfinite(a) and math.isfinite(b)):                              # nan/inf: must be the same non-finite value
            if not ((math.isnan(a) and math.isnan(b)) or a == b): out.append(f"{path}: {a!r} != {b!r}")
            return
        # relative-error diagnostics (ODE / quadrature tolerance level) may move with library versions: absolute floor 1e-7
        floor = 1e-7 if path.rsplit("/", 1)[-1].endswith(("rel", "rel_diff", "_max_rel")) or "rel_diff" in path or path.endswith("/rel") else 1e-12
        if a != b and abs(a - b) > max(1e-6 * max(abs(a), abs(b)), floor): out.append(f"{path}: {a!r} != {b!r}")
    elif a != b: out.append(f"{path}: {str(a)[:80]!r} != {str(b)[:80]!r}")


def main(folder):
    d = (ROOT / folder).resolve(); bad = 0
    for f in FILES:
        rel = (d / f).relative_to(ROOT).as_posix()
        old = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, text=True)
        if old.returncode: print(f"{rel}: not committed"); bad += 1; continue
        out = []; diff(json.loads((d / f).read_text()), json.loads(old.stdout), "", out)
        print(f"{rel}: {'OK' if not out else f'{len(out)} differences'}"); bad += bool(out)
        for line in out[:10]: print("   ", line)
    print("REPRODUCED" if not bad else f"{bad} FILE(S) DIFFER"); return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "results"))
