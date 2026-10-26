"""Fast checks: flat-space exactness of the dn scheme and the C45 reference of the homogeneous port."""
from tests.test_dn_flat import test_flat_exact
from hom.frw import c45_reference

test_flat_exact(); print("flat exact PASS")
r = c45_reference(1e4, 1e5, 1.1, 1e-4, 1)
ok = r["outcome"] == "launched60" and r["n_bounces"] == 2 and abs(r["N_inf"] - 59.6096740225035) < 1e-6
print("C45 reference", "PASS" if ok else "FAIL", r["N_inf"])
