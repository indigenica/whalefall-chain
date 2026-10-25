"""Fast reproducibility tests (~1-2 min on one core).  Usage:  python tests/run_tests.py   (also collectable by pytest).

Each test re-computes a small piece of the chain and compares it with the frozen numbers in results/ or data/."""
import json, math, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "chain")); sys.path.insert(0, str(ROOT / "lya")); sys.path.insert(0, str(ROOT / "network"))
import ucommon  # noqa: E402,F401  (sets the component paths)
import numpy as np  # noqa: E402


def test_potential_eternal_threshold():
    """same plateau: eternal threshold P_zeta = 1 at chi ~ 415.6 with the common V0 of the 10-level sum (541.0 with the pole, U26);
    flatness bound ~1.2e-6 per M_Pl"""
    from vendor_potential import SumPotential, PolePotential
    B = ucommon.benchmark("sum"); pot = SumPotential(B["V0"])
    assert abs(pot.eternal() - 415.575) < 0.01
    assert abs(math.sqrt(B["V0"] / 3) / (2 * math.pi) - 1.2038e-6) < 1e-9
    Bp = ucommon.benchmark("pole")
    assert abs(PolePotential(Bp["V0"]).eternal() - 540.953) < 0.01


def test_plateau_switch():
    """U26: the plateau switch picks the potential, the frozen pivot and the results folder; the chain results in results/ are the pole's"""
    import vendor_potential as VP
    assert ucommon.PLATEAU in ("pole", "sum") and VP.Potential is (VP.PolePotential if ucommon.PLATEAU == "pole" else VP.SumPotential)
    assert ucommon.benchmark("pole")["plateau"] == "pole" and abs(ucommon.benchmark("pole")["n_s"] - 0.97277136) < 1e-7
    s1 = json.loads((ROOT / "results" / "s1_birth.json").read_text()); u26 = json.loads((ROOT / "results" / "u26_pole_default.json").read_text())
    assert s1["provenance"]["plateau"] == "pole" and abs(s1["R_f_ES_pole_over_sum"] - 0.9356) < 1e-3 and u26["verdict_U26"] == "POLO_ADOPTADO"


def test_pivot_matching():
    """inflation end + reheating + pivot matching reproduces n_s, r, T_rh of data/cadena_cierre/thermal.json (the sum) and of
    thermal_pole.json (the pole, U26)"""
    import thermal, u22_pole as U22
    for nr, f in ((thermal.normalize_continuous(2e-8), "thermal.json"), (U22.s2_pole(U22.C_FIT, 2e-8), "thermal_pole.json")):
        arch = ucommon.load("cadena_cierre/" + f)["continuous"]
        for k in ("n_s", "r", "N_star", "T_rh_GeV", "F_GeV"):
            assert abs(nr[k] / arch[k] - 1) < 1e-6, (f, k)


def test_birth_children():
    """three C45j-viable children re-run with the C port: identical landing field and bounce count (needs gcc)"""
    if shutil.which("gcc") is None:
        raise Skip("no gcc")
    from frw_c import run_c, c45_state
    H = ucommon.load("c45b/c45h_results.json"); J = ucommon.load("c45b/c45j_results.json")
    viable = [c for c in J["children"] if c["zeta0_max_LH1"] >= 4.6e-5][:3]
    for c in viable:
        r = next(x for x in H["rows"] if x["frac_rad"] == 0.0 and x["sep"] == c["sep"] and x["lam"] == c["lam"] and x["sigma"] == c["sigma"])
        out = run_c(**c45_state(1e4, r["lam"], r["sigma"], 1, r["sep"]))
        assert out["n_bounces"] == r["n_bounces"] and abs(abs(out["chi_end"]) / r["chi_land"] - 1) < 1e-10


def test_dark_ledger():
    """vectorised emission ledger == original dark.py packet by packet; a DM root of data/cadena_cierre/dark.json"""
    import dark as D
    from s4_dark import source, dm_root
    B = ucommon.benchmark("sum"); F, Trh = B["F_GeV"], B["T_Gamma_GeV"]              # dark.json is an archive of the sum
    a = D.source(F, 3e-17, Trh, 10, 0.1); b = source(F, 3e-17, Trh, 10, 0.1)
    assert float(np.max(np.abs(np.asarray(a["nc"]) / b["nc"] - 1))) < 1e-12
    r = next(x for x in ucommon.load("cadena_cierre/dark.json")["conditional_DM_roots"] if x["efficiency"] == 0.01 and x["spectrum"] == "ir_delta")
    assert abs(math.log10(dm_root(F, Trh, 0.01, "ir_delta", 10, 0) / r["m_G_eV"])) < 1e-5


def test_closure_rows_frozen():
    """U1 closure verdict and the U2 network ratio stored in results/ are the ones documented"""
    s5 = json.loads((ROOT / "results" / "s5_closure.json").read_text())
    u2 = json.loads((ROOT / "network" / "results" / "u2" / "U2_RESULT.json").read_text())
    assert s5["verdict"] == "CIERRE_ESTRECHO" and s5["n_nodes"] == 309
    assert abs(u2["R_Gamma"]["value"] - 1.6945) < 1e-3 and abs(u2["R_Gamma"]["total"] - 0.287) < 5e-3


def test_boltzmann_cold_limit_and_wdm_anchor():
    """linear Boltzmann solver: very cold species == CDM; regression of the thermal 5.3 keV delta A (0.354 with k_max = 120;
    for reference the Viel-based fit of Dienes et al. 2112.09105 gives 0.362; that fit is NOT an independent validation)"""
    import boltz as Bz
    kh = np.geomspace(0.05, 120.0, 30)
    d0, _ = Bz.solve(kh * Bz.h)
    cold = Bz.Species(np.full(8, 1e-14), np.ones(8), -3 * np.ones(8), 0.12 / Bz.h ** 2)
    dc, _ = Bz.solve(kh * Bz.h, cold)
    assert float(np.max(np.abs(dc / d0 - 1))) < 1e-6
    sp = Bz.thermal_fd(5300.0, 0.12, nq=40, lmax=10)
    d1, _ = Bz.solve(kh * Bz.h, sp, phase_max=0.5)
    dA = Bz.deltaA(kh, (d1 / d0) ** 2, Bz.p_cdm_shape(kh * Bz.h, d0))
    assert abs(dA - 0.3541) < 2e-3


# Las tres pruebas del kernel de red viven en network/test_network.py, que se exporta tal cual al repositorio del
# artículo de la red (paper/prd_network/export/); aquí solo se importan para que `make test` las siga ejecutando.
from test_network import (test_network_kernel_small, test_network_reproduces_u2_run, test_distance_slab_exact,  # noqa: E402,F401
                          test_drain_rate_and_tau_grid_identity, test_evolve_lowmem_bit_identical, run_all, Skip)


def test_u17_collapse_solver_frozen():
    """U17/U17d collapse balance: m_G and K_eff of the B2 benchmark row (c_E = 3.199, eta = 1, R_rho = 1.56) reproduce results/lya u17d"""
    import u17_ksw as U17, u17b_class as U17B
    B = ucommon.benchmark("sum"); cE = U17B.c_energy(3.0, 0.5)
    assert abs(cE - 3.199) < 2e-3
    r = U17.solve(B["F_GeV"], B["T_Gamma_GeV"], "B2", 0.87, 10, c_eta=cE, R_rho=1.56)
    d = json.loads((ROOT / "lya" / "results" / "u17d" / "U17D_RESULT.json").read_text())
    ref = next(x for x in d["rows"] if x["tag"] == "U17d_B2_eta1_R1.56_r1")
    assert abs(r["m_G_eV"] / ref["m_G_eV"] - 1) < 1e-6 and abs(r["K"] / ref["K_eff"] - 1) < 1e-6


def test_u20b_u23_gates_frozen():
    """U20b f_iso and U23 stellar classification of two U17e rows reproduce results/u20b_iso.json and results/u23_superradiance.json"""
    import u20b_iso as I, u23_superradiance as S
    e = {r["tag"]: r for r in json.loads((ROOT / "lya" / "results" / "u17e" / "U17E_RESULT.json").read_text())["rows"]}
    fi = {(r["kind"], r["F_over_F0"]): r for r in json.loads((ROOT / "results" / "u20b_iso.json").read_text())["rows"] if r["source"] == "U17e"}
    for tag, key in (("U17e_Y_F1", ("Y", 1.0)), ("U17e_Y_F0.5", ("Y", 0.5))):
        assert abs(I.fiso(e[tag]["T_osc_MeV"], I.R_N) / fi[key]["f_iso"]["C1"] - 1) < 1e-9
    dark, light = S.boundary("dark"), S.boundary("light")
    assert S.classify(e["U17e_Y_F0.045"]["m_G_eV"], e["U17e_Y_F0.045"]["F_GeV"], dark, light)["stellar"] == "excluida"
    assert S.classify(e["U17e_Y_F0.35"]["m_G_eV"], e["U17e_Y_F0.35"]["F_GeV"], dark, light)["stellar"] == "libre"


def test_u19b_analysis_frozen():
    """U19b per-run analysis of the stored l_G = 3.2 networks reproduces the Y block of network/results/u19b/U19B_RESULT.json"""
    import u19b_analyze as A
    val = json.loads((A.D / "U19B_VALIDATION.json").read_text()); res = json.loads((A.D / "U19B_RESULT.json").read_text())
    rows = [A.per_run(json.loads(p.read_text()), val) for p in sorted((A.D / "runs").glob("RP3_N512_lG3.2_*.json"))]
    b = A.block(rows, val, "RP3", 3.2)
    assert len(rows) == 4 and abs(b["Y"] - res["blocks"]["RP3_3.2"]["Y"]) < 1e-9 and A.verdict(b, True) == res["verdict_U19b"]


def test_chainread_restart_segments():
    """U15b Nota 9-10: per-run cut, restart copies, forks and continuations on a synthetic chain folder"""
    import tempfile
    sys.path.insert(0, str(ROOT / "lya" / "u15b")); import chainread as C
    pt = lambda i, c: f"1 {i * 0.1 + c:.6e} {(i * 7 % 11) / 11:.6e} {(i * 3 % 13) / 13:.6e}\n"
    with tempfile.TemporaryDirectory() as d:
        d = Path(d); (d / "restarted_from").mkdir(); (d / "log.param").write_text("")
        g0 = [pt(i, 1) for i in range(100)]; g0.insert(30, "# After 30 accepted steps: update proposal with max(R-1) = 1\n")
        (d / "restarted_from" / "2026-10-12_1000__1.txt").write_text("".join(g0))
        g1 = g0 + [pt(100 + i, 1) for i in range(50)] + ["# x update proposal\n"] + [pt(150 + i, 1) for i in range(150)]
        (d / "2026-10-16_2000__1.txt").write_text("".join(g1))
        g2 = g1 + [pt(300 + i, 1) for i in range(20)] + ["# y update proposal\n"] + [pt(320 + i, 1) for i in range(80)]
        (d / "2026-10-17_3000__1.txt").write_text("".join(g2))
        old = [pt(i, 2) for i in range(300)]; (d / "2026-10-13_200000__2.txt").write_text("".join(old))
        # Nota 11: a continuation repeats the end of the old file (here after a NUL hole and a partial line) and counts after it;
        # a file that only shares the suffix is a separate chain
        hole = "\x00" * 40 + old[-3][5:] + old[-2] + old[-1]
        (d / "2026-10-17_300000__2.txt").write_text(hole + "".join([pt(1000 + i, 2) for i in range(50)]))
        (d / "2026-10-13_200000__3.txt").write_text("".join(pt(i, 3) for i in range(100)))
        (d / "2026-10-17_300000__3.txt").write_text("".join(pt(2000 + i, 3) for i in range(30)))
        ch, rep = C.read(d, 4)
        assert sorted(len(c) for c in ch) == [30, 100, 300, 350]
        st = {r["file"]: r["status"] for r in rep}
        assert st["2026-10-16_2000__1.txt"].startswith("copiado") and st["2026-10-17_300000__2.txt"].startswith("continuación")
        assert "cadena aparte" in st["2026-10-17_300000__3.txt"]



def test_restart_move_and_rejoin():
    """U15b Nota 12: a new file still copying its old chain is NOT taken as contained in the old one; the two halves of a
    file split by a premature move are rejoined, and then the old file is contained in the joined one"""
    import tempfile
    sys.path.insert(0, str(ROOT / "lya" / "u15b")); import restarted as R, rejoin_moved as J, chainread as C
    pt = lambda i, c: f"1\t{i * 0.1 + c:.6e}\t{(i * 7 % 11) / 11:.6e}\t{(i * 3 % 13) / 13:.6e}\n"
    with tempfile.TemporaryDirectory() as d:
        d = Path(d); (d / "restarted_from").mkdir(); (d / "log.param").write_text("")
        old = [pt(i, 1) for i in range(200)]; (d / "2026-10-17_200000__1.txt").write_text("".join(old))
        (d / "2026-10-19_300000__1.txt").write_text("".join(old[:80]))              # the copy, half written
        assert R.contained(d) == []                                                # before Nota 12: the NEW file was 'contained'
        # what happened: the half-copied file was moved; MontePython finished the copy in it, wrote 10 points, then reopened by name
        (d / "2026-10-19_300000__1.txt").rename(d / "restarted_from" / "2026-10-19_300000__1.txt")
        (d / "restarted_from" / "2026-10-19_300000__1.txt").write_text("".join(old + [pt(200 + i, 1) for i in range(10)]))
        (d / "2026-10-19_300000__1.txt").write_text("".join(pt(210 + i, 1) for i in range(40)))
        assert J.main(["x", str(d), "--apply"]) == 0
        joined = (d / "2026-10-19_300000__1.txt").read_text().splitlines(keepends=True)
        assert joined == old + [pt(200 + i, 1) for i in range(50)]
        assert not (d / "2026-10-17_200000__1.txt").exists() and (d / "restarted_from" / "2026-10-17_200000__1.txt").exists()
        ch, rep = C.read(d, 4)
        assert [len(c) for c in ch] == [250]


TESTS = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]

if __name__ == "__main__":
    sys.exit(run_all(TESTS))
