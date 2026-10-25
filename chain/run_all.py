"""Run the whole unified chain: S1 birth, S2 pocket, S3 sectors, S4 dark-ledger checks, S5 closure, figures, S6, U17 (K_sw vs the network energy at H = m_G/3), contract.
Usage (from the repository root):  python -B chain/run_all.py      (~6 min on one CPU core; S1 needs gcc, else SKIPPED; S6 reads network/results/u2/U2_RESULT.json, the lattice runs are in network/)"""
import time
import ucommon  # noqa: F401
import s1_birth, s2_pocket, s3_sectors, s4_dark, s5_closure, figures, s6_measured, u17_ksw, u_contract

if __name__ == "__main__":
    t0 = time.time()
    for mod in (s1_birth, s2_pocket, s3_sectors, s4_dark, s5_closure, figures, s6_measured, u17_ksw, u_contract):   # U17 before the contract that reads it
        mod.main()
    print(f"run_all done in {time.time() - t0:.0f}s")
