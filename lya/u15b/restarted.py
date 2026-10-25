"""Find (and move aside) chain files that a MontePython restart (-r) has already copied into a newer file.

MontePython -r copies each old chain at the start of a new file and only deletes the old one when the new chain reaches
its N steps; until then both files are in the folder and every old point would be counted twice.  A file is 'contained'
when its points are exactly the first points of another, longer chain file of the same folder.

    python lya/u15b/restarted.py ~/u15b/chains            # report, for every case folder
    python lya/u15b/restarted.py ~/u15b/chains --move     # move the contained files to <case>/restarted_from/ (not deleted)

analyze.py uses contained() to leave those files out even if they were not moved."""
import glob, os, shutil, sys


def points(f):
    with open(f) as fh:
        return [l for l in fh if not l.startswith("#") and l.strip()]


def name_key(f):
    """(N, date) from a MontePython chain name YYYY-MM-DD_N__k.txt: a restart always writes a larger N than the file it copies"""
    b = os.path.basename(f)
    try:
        return (int(b.split("__")[0].split("_")[1]), b[:10])
    except (IndexError, ValueError):
        return (-1, b[:10])


def contained(folder):
    """[(old file, newer file that starts with all of its points)], for the chain files at the top of `folder`.
    Nota 12: the newer file must have a larger N in its name.  Without that check, a new file still being written (its copy of
    the old chain only half done) is a prefix of the OLD file and was taken as 'contained' in it: move_restarted.sh then moved
    the NEW file away, and MontePython, which reopens its output by name (io_mp.refresh_file), went on in a fresh file without
    the copy (U15c_B2_F0.5_r1, 19 Oct 2026)."""
    files = sorted(glob.glob(os.path.join(folder, "*__*.txt")))
    pts = {f: points(f) for f in files}
    out = []
    for a in files:
        pa = pts[a]
        if len(pa) < 5: continue                          # (fresh chains share their first point, the starting one)
        for b in files:
            if name_key(b)[0] <= name_key(a)[0]: continue  # only a later restart (larger N) can contain an older file
            pb = pts[b]
            # the last old line may be cut (chain killed mid-write): its copy then runs into the next point of the new chain
            if b != a and len(pb) > len(pa) and pb[:len(pa) - 1] == pa[:-1] and pb[len(pa) - 1].startswith(pa[-1].rstrip("\n")):
                out.append((a, b)); break
    return out


def main(argv):
    base = argv[1]; move = "--move" in argv
    folders = [base] if glob.glob(os.path.join(base, "*__*.txt")) else sorted(d for d in glob.glob(os.path.join(base, "*")) if os.path.isdir(d))
    for d in folders:
        c = contained(d)
        print(f"{os.path.basename(d)}: {len(c)} file(s) already copied into a restarted chain")
        for a, b in c:
            print(f"   {os.path.basename(a)}  ->  in {os.path.basename(b)}")
            if move:
                os.makedirs(os.path.join(d, "restarted_from"), exist_ok=True)
                shutil.move(a, os.path.join(d, "restarted_from", os.path.basename(a)))
        if move and c: print(f"   moved to {os.path.join(d, 'restarted_from')}/")


if __name__ == "__main__":
    main(sys.argv)
