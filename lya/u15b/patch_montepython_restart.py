"""Fix a race in MontePython (commit 5b365b0) io_mp.create_output_files when several MPI chains start together (ESPEC_U15b Nota 10).

Each process opens its candidate chain file with mode 'w' (which truncates it) and only then tries to lock it.  When two processes
pick the same name, the loser truncates the winner's file, which may already hold the copy of the old chain written by a
restart (-r); the winner then keeps appending its new points after an empty start.  Seen on 17 Oct 2026: after 'boost', the
restarted chains of both U15b validation cases and two of U15c started without their copies.  This opens with mode 'a' (no
truncation), takes the lock, and only the owner truncates.  Idempotent; touches nothing else.

    python lya/u15b/patch_montepython_restart.py /path/to/montepython_public"""
import sys
from pathlib import Path

f = Path(sys.argv[1]) / "montepython" / "io_mp.py"
s = f.read_text()
if "# whalefall: lock before truncating" in s:
    print("already patched"); sys.exit(0)
old = """        while trying:
            data.out = open(os.path.join(
                command_line.folder, outname_base)+str(suffix)+'.txt', 'w')
            try:
                lock(data.out, fcntl.LOCK_EX | fcntl.LOCK_NB)
                trying = False
            except LockError:
                suffix += 1"""
new = """        while trying:   # whalefall: lock before truncating (opening with 'w' truncated other chains' files)
            data.out = open(os.path.join(
                command_line.folder, outname_base)+str(suffix)+'.txt', 'a')
            try:
                lock(data.out, fcntl.LOCK_EX | fcntl.LOCK_NB)
                data.out.seek(0); data.out.truncate()
                trying = False
            except LockError:
                data.out.close()
                suffix += 1"""
assert old in s, "unexpected io_mp.py: check the MontePython commit"
f.write_text(s.replace(old, new))
print("patched", f)
