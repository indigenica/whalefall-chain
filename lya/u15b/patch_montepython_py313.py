"""Make MontePython (commit 5b365b0) run on Python >= 3.13 (PEP 667).

montepython/data.py imports each likelihood class with exec("from likelihoods.X import X") and instantiates it in a second
exec; since Python 3.13 the name bound by the first exec is not visible to the second (NameError).  This replaces the two
exec calls by importlib + getattr.  Idempotent; touches nothing else.

    python lya/u15b/patch_montepython_py313.py /path/to/montepython_public"""
import sys
from pathlib import Path

f = Path(sys.argv[1]) / "montepython" / "data.py"
s = f.read_text()
if "importlib.import_module(\"likelihoods." in s:
    print("already patched"); sys.exit(0)
old1 = """                exec("from likelihoods.%s import %s" % (
                    elem, elem))"""
new1 = """                import importlib
                _lkl_class = getattr(importlib.import_module("likelihoods.%s" % elem), elem)"""
old2 = """                exec("self.lkl['%s'] = %s('%s/%s.data',\\
                    self, self.command_line)" % (
                    elem, elem, folder, elem))"""
new2 = """                self.lkl[elem] = _lkl_class('%s/%s.data' % (folder, elem), self, self.command_line)"""
assert old1 in s and old2 in s, "unexpected data.py: check the MontePython commit"
f.write_text(s.replace(old1, new1).replace(old2, new2))
print("patched", f)
