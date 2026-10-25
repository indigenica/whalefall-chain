"""Reading MontePython chain folders with restarts (ESPEC_U15b, Nota 9): which files count, and the non-Markovian cut per run.

A restarted chain file (-r) starts with a verbatim copy of the old file (points and comment lines) and then continues.  The
MontePython rule 'drop every point before the last update proposal' is applied per run: inside each copied segment (one
earlier process) and inside the new segment separately, never across runs, because a proposal update in a later process cannot
make the points of an earlier process (generated under its own fixed final proposal) non-Markovian.

Files at the top of the folder that a newer file copied entirely are skipped (Nota 3).  A 'fork' (an old file that kept growing
after a newer file copied its first points, i.e. the old process was not stopped) is reported; the shared prefix is counted once
(in the old file) and the newer file only from the end of the shared prefix.  Old files in restarted_from/ are never counted, but
they mark segment boundaries.  A 'continuation' (Nota 10: a restart whose copy was lost or truncated by a MontePython race) is read
together with the old file, as one chain.  Revised after the audit of 18 Oct (Nota 11): a continuation is recognised by content only
(after any leading NUL hole and partial line, its first points repeat the END of the old file; it counts from the end of that
overlap); a match by name alone makes it a separate chain, with a warning; the lineage follows copies (old -> continuation -> a
restart that copied the continuation is one chain); forks are ordered by (date, N, mtime).

    python lya/u15b/chainread.py ~/u15b/chains          # diagnosis of every case folder (files, segments, cuts, forks)"""
import glob, os, re, sys

MIN_PREFIX = 20


def parse(f):
    """(point lines, positions of 'update proposal' marks as the number of points before each mark)"""
    pts, marks = [], []
    with open(f) as fh:
        for line in fh:
            line = line.lstrip("\x00")                     # a NUL hole left by the truncation race (Nota 11)
            if line.startswith("#"):
                if "update proposal" in line: marks.append(len(pts))
            elif line.strip():
                pts.append(line)
    return pts, marks


def prefix_len(b, a):
    """number of leading points of a that b repeats (a's last line may be cut mid-write and merged into b's next line)"""
    n = min(len(a), len(b)); k = 0
    while k < n and b[k] == a[k]: k += 1
    if k == len(a) - 1 and k < len(b) and b[k].startswith(a[k].rstrip("\n")): k += 1
    return k


def steps_N(f):
    m = re.search(r"_(\d+)__\d+\.txt$", os.path.basename(f)); return int(m.group(1)) if m else 0


def file_key(f):
    """chronological order of chain files: (date in the name, N, mtime)"""
    return (os.path.basename(f)[:10], steps_N(f), os.path.getmtime(f))


def pkey(line):
    """a point without its multiplicity (first column), for matching copies"""
    return " ".join(line.split()[1:])


def overlap(old, new, max_skip=3):
    """index in `new` where its own points start, if `new` begins (after at most max_skip garbled lines) by repeating the last
    points of `old` up to its very end; 0 if it does not"""
    if not old or not new: return 0
    ko = [pkey(x) for x in old[-(len(new) + max_skip):]]; off = len(old) - len(ko)
    for j in range(min(max_skip, len(new))):
        kj = pkey(new[j])
        for i in (i for i, k in enumerate(ko) if k == kj):
            n = len(ko) - i
            if j + n <= len(new) and all(pkey(new[j + t]) == ko[i + t] for t in range(n)):
                return j + n
    return 0


def layout(folder):
    """{file: dict(points, marks, boundaries, start, status)} for the chain files at the top of the folder"""
    top = sorted(glob.glob(os.path.join(folder, "*__*.txt")))
    old = sorted(glob.glob(os.path.join(folder, "restarted_from", "*__*.txt")))
    P = {f: parse(f) for f in top + old}
    info = {f: dict(points=len(P[f][0]), marks=P[f][1], boundaries=[], start=0, status="cuenta", source=[]) for f in top}
    for b in top:
        pb = P[b][0]
        for a in top + old:
            if a == b or len(P[a][0]) < MIN_PREFIX: continue
            k = prefix_len(pb, P[a][0])
            if k < MIN_PREFIX: continue
            if k == len(P[a][0]) and len(pb) > k:                 # b copied all of a: a segment boundary of b
                info[b]["boundaries"].append((k, len(P[a][1]))); info[b]["source"].append(os.path.basename(a))
                if a in info: info[a]["status"] = "copiado (no cuenta)"
            elif a in info and k < len(P[a][0]) and k < len(pb):  # fork: shared prefix, both grew afterwards
                newer, older = (b, a) if file_key(b) > file_key(a) else (a, b)
                if newer == b:
                    info[b]["start"] = max(info[b]["start"], k); info[b]["status"] = f"bifurcación de {os.path.basename(a)} (cuenta desde el punto {k})"
    for f in info: info[f]["boundaries"] = sorted(set(info[f]["boundaries"]))
    # continuations (Nota 10, revised in Nota 11): by content only; a copied continuation passes its lineage to its copier
    copier = {}
    for b in top:
        for src in info[b]["source"]: copier[src] = b
    byname = {os.path.basename(f): f for f in top}
    for b in top:
        if info[b]["boundaries"] or info[b]["start"] or not P[b][0] or info[b].get("continues_from"): continue
        best = None
        for a in top:
            if a == b or file_key(a) >= file_key(b) or not P[a][0] or info[a].get("continued_by"): continue
            k = overlap(P[a][0], P[b][0])
            if k and (best is None or k > best[1]): best = (a, k)
        if not best:
            cand = [a for a in top if a != b and suffix(a) == suffix(b) and file_key(a) < file_key(b) and len(P[a][0]) >= MIN_PREFIX]
            if cand and not info[b]["status"].startswith("copiado"):
                info[b]["status"] += f" (aviso: mismo sufijo que {os.path.basename(max(cand, key=file_key))}, pero sin repetir su final: cadena aparte)"
            continue
        a, k = best; target = b
        while os.path.basename(target) in copier and info[target]["status"].startswith("copiado"):
            target = copier[os.path.basename(target)]                          # the restart that copied the continuation
        info[target]["start"] = max(info[target]["start"], k)
        info[target]["continues"] = a; info[a]["continued_by"] = target; info[b]["continues_from"] = a
        note = f"continuación de {os.path.basename(a)} (cuenta desde el punto {k})"
        if target == b: info[b]["status"] = note
        else: info[target]["status"] = f"{info[target]['status']}; {note}, a través de la copia de {os.path.basename(b)}"
    return P, info


def suffix(f):
    m = re.search(r"__(\d+)\.txt$", os.path.basename(f)); return int(m.group(1)) if m else -1


def same_point(l1, l2, rtol=1e-6):
    """same -lnL and parameters (the multiplicity, first column, may differ)"""
    a, b = l1.split()[1:], l2.split()[1:]
    if len(a) != len(b) or not a: return False
    try: return all(abs(float(x) - float(y)) <= rtol * max(abs(float(x)), abs(float(y)), 1e-300) for x, y in zip(a, b))
    except ValueError: return False


def kept_rows(pts, marks, boundaries, start):
    """indices of the points kept.  boundaries = [(L, M)]: an earlier run's file of L points with M update marks, copied at the
    start of this one (nested restarts give nested prefixes).  Per run, drop the points before that run's last update mark."""
    segs, lo, m0 = [], 0, 0
    for L, M in sorted(boundaries):
        if L <= lo or L >= len(pts): continue
        segs.append((lo, L, marks[m0:M])); lo, m0 = L, M
    segs.append((lo, len(pts), marks[m0:]))
    keep = []
    for lo, hi, mk in segs:
        cut = max(mk) if mk else lo
        keep += range(max(cut, lo, start), hi)
    return keep


def read(folder, ncols):
    """chains of a case folder: list of row lists (floats) after the per-run non-Markovian cut (before any burn-in), and a report"""
    P, info = layout(folder); chains, report = [], []
    def rows_of(f):
        pts, marks = P[f]; d = info[f]; rows = []
        for i in kept_rows(pts, marks, d["boundaries"], d["start"]):
            v = pts[i].split()
            if len(v) == ncols:
                try: rows.append([float(x) for x in v])
                except ValueError: pass
        report.append(dict(file=os.path.basename(f), **{k: (os.path.basename(v) if k in ("continues", "continued_by", "continues_from") else v)
                                                        for k, v in d.items() if k != "marks"}, n_marks=len(marks), kept=len(rows)))
        return rows
    for f, d in info.items():
        if d["status"].startswith("copiado"):
            report.append(dict(file=os.path.basename(f), **{k: (os.path.basename(v) if k in ("continues", "continued_by", "continues_from") else v)
                                                            for k, v in d.items() if k != "marks"}, n_marks=len(P[f][1]), kept=0)); continue
        if d.get("continues"): continue                                   # read with the file it continues
        rows, g = rows_of(f), f
        while info[g].get("continued_by"):                                # one chain: old run + its continuations, in order
            g = info[g]["continued_by"]; rows += rows_of(g)
        chains.append(rows)
    return chains, report


def main(argv):
    base = argv[1] if len(argv) > 1 else os.path.expanduser("~/u15b/chains")
    for d in sorted(x for x in glob.glob(os.path.join(base, "*")) if os.path.isdir(x) and os.path.exists(os.path.join(x, "log.param"))):
        _, info = layout(d); P = {f: parse(f) for f in info}
        print(f"{os.path.basename(d)}:  {len(info)} ficheros arriba, {len(glob.glob(os.path.join(d, 'restarted_from', '*__*.txt')))} en restarted_from/")
        for f, x in info.items():
            pts, marks = P[f]; idx = kept_rows(pts, marks, x["boundaries"], x["start"])
            if x.get("continued_by"): x = dict(x, status=x["status"] + f", sigue en {os.path.basename(x['continued_by'])}")
            print(f"   {os.path.basename(f):34s} {x['points']:6d} puntos, marcas 'update proposal' en {marks[-4:]}{'...' if len(marks) > 4 else ''}; "
                  f"tramos {x['boundaries']}; {x['status']}; se usan {len(idx) if not x['status'].startswith('copiado') else 0} (antes del 30 %)")


if __name__ == "__main__":
    main(sys.argv)
