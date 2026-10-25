"""Rejoin chain files that were split by a premature move (ESPEC_U15b, Nota 12).  Nothing is deleted.

What happened (U15c_B2_F0.5_r1, 19 Oct 2026): just after a restart, move_restarted.sh / restarted.py took three NEW chain
files, still copying their old chains, for files 'contained' in the old ones and moved them to restarted_from/.  MontePython
keeps writing to the moved file until it reopens its output by name (io_mp.refresh_file); from then on it writes into a new
file with the same name at the top of the folder.  The output of ONE process is therefore split in two files with the same name:
    restarted_from/<name>   = copy of the old chain + the first new points
    <name> (top)            = the later new points
This script, for a case whose chains are NOT running, writes <name> = restarted_from/<name> + <name> at the top and keeps
both originals in restarted_from/_rejoined_<stamp>/ (a subfolder: chainread and restarted.py do not read it), with a log.
Then restarted.py --move moves the old files that the rejoined ones now contain.

    python lya/u15b/rejoin_moved.py ~/u15b/chains/U15c_B2_F0.5_r1            # report only
    python lya/u15b/rejoin_moved.py ~/u15b/chains/U15c_B2_F0.5_r1 --apply    # rejoin (chains of the case must be stopped)"""
import glob, hashlib, json, os, shutil, subprocess, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import restarted  # noqa: E402


def sha(f):
    return hashlib.sha256(open(f, "rb").read()).hexdigest()


def running(folder):
    r = subprocess.run(["pgrep", "-f", f"MontePython.py run -o {folder.rstrip('/')}( |$)"], capture_output=True, text=True)
    return bool(r.stdout.strip())


def pairs(folder):
    """top files whose exact name is also in restarted_from/ (the two halves of one process's output)"""
    out = []
    for top in sorted(glob.glob(os.path.join(folder, "*__*.txt"))):
        moved = os.path.join(folder, "restarted_from", os.path.basename(top))
        if os.path.exists(moved):
            out.append((moved, top))
    return out


def main(argv):
    folder = os.path.abspath(os.path.expanduser(argv[1])); apply = "--apply" in argv
    P = pairs(folder)
    if not P:
        print(f"{os.path.basename(folder)}: nothing split (no top file has a namesake in restarted_from/)"); return 0
    rows = []
    for moved, top in P:
        a, b = restarted.points(moved), restarted.points(top)
        rows.append(dict(name=os.path.basename(top), first_part=len(a), second_part=len(b),
                         first_last=a[-1].split("\t")[:4] if a else None, second_first=b[0].split("\t")[:4] if b else None))
        print(f"   {os.path.basename(top)}: restarted_from/ {len(a)} points + top {len(b)} points")
    if not apply:
        print("report only; with --apply (chains of this case stopped) the halves are joined"); return 0
    if running(folder):
        print("the chains of this case are running: stop them first (bash lya/u15b/run_wsl.sh stop <case>)"); return 1
    stamp = time.strftime("%Y%m%d_%H%M%S"); keep = os.path.join(folder, "restarted_from", f"_rejoined_{stamp}")
    os.makedirs(os.path.join(keep, "top")); os.makedirs(os.path.join(keep, "restarted_from"))
    for (moved, top), row in zip(P, rows):
        row.update(sha_first=sha(moved), sha_second=sha(top))
        tmp = top + ".rejoin_tmp"
        with open(tmp, "w") as out:
            for f in (moved, top):
                with open(f) as fh:
                    shutil.copyfileobj(fh, out)
        shutil.move(top, os.path.join(keep, "top", os.path.basename(top)))
        shutil.move(moved, os.path.join(keep, "restarted_from", os.path.basename(moved)))
        os.replace(tmp, top); row["sha_joined"] = sha(top)
    json.dump(dict(folder=folder, when=stamp, note="ESPEC_U15b Nota 12", files=rows), open(os.path.join(keep, "REJOIN_LOG.json"), "w"), indent=1)
    print(f"joined {len(P)} file(s); originals and log in {keep}")
    for a, b in restarted.contained(folder):
        os.makedirs(os.path.join(folder, "restarted_from"), exist_ok=True)
        shutil.move(a, os.path.join(folder, "restarted_from", os.path.basename(a)))
        print(f"   {os.path.basename(a)} is contained in {os.path.basename(b)}: moved to restarted_from/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
