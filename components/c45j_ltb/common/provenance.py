"""Small provenance helpers: run directories, source hashes, JSON output."""
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_hashes():
    out = {}
    for p in sorted(ROOT.rglob("*.py")):
        if "results" in p.parts or "__pycache__" in p.parts:
            continue
        out[str(p.relative_to(ROOT))] = sha256(p)
    return out


def new_run_dir(tag, base=None):
    base = Path(base) if base else ROOT / "results"
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    d = base / f"{tag}_{stamp}"
    d.mkdir(parents=True, exist_ok=False)
    return d


def environment():
    import numpy, scipy
    return dict(python=sys.version.split()[0], numpy=numpy.__version__,
                scipy=scipy.__version__, platform=platform.platform(),
                host=platform.node(), cpu_count=os.cpu_count(),
                utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


def _default(o):
    import numpy as np
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_default, ensure_ascii=False)


def write_manifest(run_dir):
    files = {}
    for p in sorted(Path(run_dir).rglob("*")):
        if p.is_file() and p.name != "MANIFEST.json":
            files[str(p.relative_to(run_dir))] = sha256(p)
    write_json(Path(run_dir) / "MANIFEST.json", files)
