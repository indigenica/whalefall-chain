"""U26 (docs/preregistros/ESPEC_U26.md): the pole plateau as the chain default.  Reads the chain run with the pole (results/), the
same chain with WHALEFALL_PLATEAU=sum (results/plateau_sum/) and the common action with the pole (data/cierre_accion_pole/), and
evaluates the rows C1-C5 and the verdict fixed in the ESPEC.  Run `make u26` (both chains, the action and this script).
    python -B chain/u26_pole_default.py      ->  results/u26_pole_default.json"""
import json, re
from ucommon import ROOT, INPUTS

P, S = ROOT / "results", ROOT / "results" / "plateau_sum"
J = lambda d, n: json.loads((d / n).read_text())
CLASS_VERDICTS = ["U3a", "U4a", "U5a", "U5f", "U7", "U7c", "U12", "U15", "U17b", "U17c", "U17d", "U17e", "U21", "U20b", "U23"]


def numbers(d, p="", out=None):
    out = {} if out is None else out
    if isinstance(d, dict):
        for k, v in d.items(): numbers(v, f"{p}/{k}", out)
    elif isinstance(d, list):
        for i, v in enumerate(d): numbers(v, f"{p}[{i}]", out)
    elif isinstance(d, (int, float)) and not isinstance(d, bool): out[p] = d
    return out


def m_G_shift(name):
    """max |m_G(pole)/m_G(sum) - 1| over every m_G field of a stage result, paired by path"""
    a, b = numbers(J(P, name)), numbers(J(S, name))
    ka, kb = ({k for k in d if re.search(r"m_G", k.rsplit("/", 1)[-1])} for d in (a, b))
    keys = sorted(k for k in ka & kb if b[k])
    # a mass present in one run and missing (or zero in the sum) in the other is not comparable: C5 fails (Codex, PR #11)
    missing = sorted(ka ^ kb) + sorted(k for k in ka & kb if not b[k])
    return dict(n=len(keys), n_sum=len(kb), n_pole=len(ka), not_comparable=missing[:20], same_keys=not missing,
                max_abs_rel=max(abs(a[k] / b[k] - 1) for k in keys) if keys else float("inf"))


def archive_agrees(fresh, frozen, rtol=1e-6):
    """the C4 quantities of the frozen pole action (what S3 reads) equal those of the fresh run"""
    q = lambda d: dict(checks=J(d, "execution.json")["checks"], proof=J(d, "portal.json")["nonnegative_total_toy_potential_proof"],
                       lower=J(d, "portal.json")["first_level_lower_bound_coefficient"], neg=J(d, "portal.json")["negative_correction_bound"],
                       ratio=J(d, "decays.json")["min_Q_width_over_Gamma_chi"],
                       T=[c["final_visible_T_GeV"] for c in J(d, "late_fluid.json")["matching_scan"]])
    a, b = q(fresh), q(frozen)
    close = lambda x, y: abs(x - y) <= rtol * max(abs(x), abs(y))
    ok = (a["checks"] == b["checks"] and a["proof"] == b["proof"] and all(close(a[k], b[k]) for k in ("lower", "neg", "ratio"))
          and len(a["T"]) == len(b["T"]) and all(close(x, y) for x, y in zip(a["T"], b["T"])))
    return dict(agrees=bool(ok), rtol=rtol)


def main():
    s1, s2, s3, s4, s5, s6, u17 = (J(P, n) for n in ("s1_birth.json", "s2_pocket.json", "s3_sectors.json", "s4_dark_checks.json",
                                                    "s5_closure.json", "s6_measured_closure.json", "u17_ksw.json"))
    sums = {n: J(S, n) for n in ("s1_birth.json", "s2_pocket.json", "s3_sectors.json", "s4_dark_checks.json", "s5_closure.json", "s6_measured_closure.json", "u17_ksw.json")}
    for n, d in [("s1", s1), ("s2", s2), ("s3", s3), ("s4", s4), ("s5", s5), ("s6", s6)]:
        assert d["provenance"]["plateau"] == "pole", f"{n} in results/ was not run with the pole"
    assert u17.get("plateau") == "pole", "u17_ksw.json in results/ was not run with the pole"
    assert all((d.get("provenance") or {}).get("plateau", d.get("plateau")) == "sum" for d in sums.values()), "results/plateau_sum/ was not run with the sum"
    # C4 is evaluated on the fresh pole action of this `make u26` (action/results_pole/), and the frozen copy that S3 reads
    # (data/cierre_accion_pole/) must agree with it (Codex, PR #11)
    A, Z = ROOT / "action" / "results_pole", INPUTS / "cierre_accion_pole"
    if not (A / "execution.json").exists(): raise SystemExit("no action/results_pole/: run `make u26` (it runs the pole action)")
    ex, po, de = J(A, "execution.json"), J(A, "portal.json"), J(A, "decays.json")
    assert ex["plateau"] == "pole"
    archive = archive_agrees(A, Z)
    if "f_ES_pole" not in s1: raise SystemExit("S1 was skipped (no gcc): U26 needs the C port of C45 for row C2")
    rows = {}
    rows["C1_pivot"] = dict(pivot_status=s2["pivot_status"], pivot_max_rel=s2["pivot_max_rel"], **s2["pole_vs_sum"],
                            sum_control=sums["s2_pocket.json"]["pivot_status"],
                            passed=bool(s2["pivot_status"] == "PASS" and s2["pole_vs_sum"]["cmb_equivalent"] and sums["s2_pocket.json"]["pivot_status"] == "PASS"))
    rows["C2_births"] = dict(control_sum_mode0=s1["control_status_sum_mode0"], f_ES_pole=s1["f_ES_pole"], f_ES_sum=s1["f_ES_sum"],
                             R=s1["R_f_ES_pole_over_sum"], chi_et_pole=s1["chi_et_pole"], chi_et_sum=s1["chi_et_sum"],
                             viable_C45j_eternal_pole=s1["viable_eternal_pole"], n_viable_C45j=s1["n_viable_C45j"],
                             grid_mode0_max_rel_diff=s1["grid_mode0_max_rel_diff"],
                             passed=bool(s1["control_status_sum_mode0"] == "PASS" and s1["grid_mode0_max_rel_diff"] <= 1e-8 and s1["R_f_ES_pole_over_sum"] >= 0.5))
    rows["C3_eternal"] = dict(chi_et=s2["eternal"]["chi_Pzeta_1"], chi_vol=s2["eternal"]["chi_volume_Pzeta_1_over_6"],
                              flatness_bound_per_MPl=s2["eternal"]["flatness_bound_per_MPl"], max_ratio_above_threshold=s2["eternal"]["max_ratio_above_threshold"],
                              C45i_carried=s2["D3"]["C45i_verdict"], passed=bool(s2["eternal"]["undeformed_plateau_meets_bound"]))
    rows["C4_action"] = dict(checks=ex["checks"], min_Q_width_over_Gamma_chi=de["min_Q_width_over_Gamma_chi"],
                             portal_negative_correction=po["negative_correction_bound"], portal_lower_bound=po["first_level_lower_bound_coefficient"],
                             portal_proof=po["nonnegative_total_toy_potential_proof"], S3_status=s3["status"], T_late_over_T_rh_minus_1=s3["late_vs_S2_rel"],
                             F_identity_rel=s3["F_identity_rel"],
                             archive_agrees_with_fresh_run=archive,
                             passed=bool(ex["pass_execution"] and de["min_Q_width_over_Gamma_chi"] > 1 and po["nonnegative_total_toy_potential_proof"]
                                         and s3["status"] == "PASS" and archive["agrees"]))
    labels = dict(S4=(s4["status"], sums["s4_dark_checks.json"]["status"]), S5=(s5["verdict"], sums["s5_closure.json"]["verdict"]),
                  S5_nodes=(s5["n_nodes"], sums["s5_closure.json"]["n_nodes"]), S6=(s6["verdict"], sums["s6_measured_closure.json"]["verdict"]),
                  U17=(u17["summary"]["verdict_U17"], sums["u17_ksw.json"]["summary"]["verdict_U17"]))
    expected = dict(S4="PASS", S5="CIERRE_ESTRECHO", S5_nodes=309, S6="FALSADA_EN_FORMA_ACTUAL", U17="K_SW_FUERA_DEL_RANGO_TIBIA")
    shift = {n: m_G_shift(n) for n in ("s5_closure.json", "s6_measured_closure.json", "u17_ksw.json")}
    worst = max(x["max_abs_rel"] for x in shift.values())
    rows["C5_dark"] = dict(labels_pole_vs_sum=labels, labels_kept=all(labels[k][0] == expected[k] for k in expected),
                           pole_roots_vs_sum_dlog10=s4.get("pole_roots_vs_sum_dlog10"), m_G_shift=shift, m_G_max_abs_rel=worst,
                           class_verdicts_carried=bool(worst <= 0.02 and all(x["same_keys"] for x in shift.values())), class_verdicts=CLASS_VERDICTS)
    rows["C5_dark"]["same_m_G_keys"] = all(x["same_keys"] for x in shift.values())
    rows["C5_dark"]["passed"] = bool(rows["C5_dark"]["labels_kept"] and rows["C5_dark"]["class_verdicts_carried"] and rows["C5_dark"]["same_m_G_keys"])
    if not (rows["C1_pivot"]["passed"] and rows["C2_births"]["passed"]): v = "POLO_NO_ADOPTADO"
    elif all(rows[k]["passed"] for k in ("C3_eternal", "C4_action", "C5_dark")): v = "POLO_ADOPTADO"
    else: v = "POLO_ADOPTADO_CON_CAMBIOS"
    out = dict(spec="docs/preregistros/ESPEC_U26.md", rows=rows, verdict_U26=v,
               failed_rows=[k for k, x in rows.items() if not x["passed"]])
    (P / "u26_pole_default.json").write_text(json.dumps(out, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
    for k, x in rows.items(): print(k, "PASS" if x["passed"] else "FAIL")
    print(f"m_G max shift {worst:.2%}; f_ES R = {s1['R_f_ES_pole_over_sum']:.3f}; dn_s = {s2['pole_vs_sum']['dn_s']:+.1e}")
    print("U26:", v)
    return out


if __name__ == "__main__":
    main()
