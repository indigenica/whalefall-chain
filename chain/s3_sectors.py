"""S3 · after the pocket: Q ordering on RP3, Z2 strings, free species of the common action, vacuum decays and the late fluid.

Nothing here is re-simulated (the 3D network and octree runs are hours of GPU/CPU; see docs/cierre_*): S3 reads the
archived results, checks that they use the same parameters as S1-S2 and collects their status for the contract.

U26 (docs/preregistros/ESPEC_U26.md), row C4: with the pole plateau the common action is the one re-run with the pole
(data/cierre_accion_pole/, from `WHALEFALL_PLATEAU=pole make action`), and F is checked against thermal_pole.json; the F of the
cadena_cierre contract (the sum) stays as the control with WHALEFALL_PLATEAU=sum."""
from ucommon import ROOT, load, save, benchmark, provenance, PLATEAU, THERMAL

ACTION = dict(sum="cierre_accion", pole="cierre_accion_pole")


def main():
    B = benchmark(); A = ACTION[PLATEAU]
    sp = load("cierre_cosmologia/species.json"); de = load(f"{A}/decays.json"); lf = load(f"{A}/late_fluid.json")
    cs = load("cierre_cosmologia/closure_status.json"); sc = load("cierre_red/source_contract.json")
    e2e = load("cierre_cosmologia/end_to_end_control.json"); ct = load("cadena_cierre/contract.json")
    T_late = [c["final_visible_T_GeV"] for c in lf["matching_scan"]]
    dneff_late = max(c["Delta_Neff_massless_homogeneous_Goldstones"] for c in lf["matching_scan"])
    # with the pole, F is checked against the copy of the matching that the common action actually read (action/vendor/)
    F_ref = (ct["shared_parameters"]["F_GeV"] if PLATEAU == "sum" else
             __import__("json").loads((ROOT / "action" / "vendor" / "thermal_pole.json").read_text())["continuous"]["F_GeV"])
    same_F = abs(F_ref / B["F_GeV"] - 1)
    out = dict(stage="S3_sectors", plateau=PLATEAU, action_results=f"data/{A}/", provenance=provenance(),
               species=[dict(name=s["name"], multiplicity=s["multiplicity"], mass_code=s["mass_code"], width_code=s["width_code"]) for s in sp["species"]],
               stable_massive_species_in_action=sp["stable_massive_particle_species"],
               goldstone_mass_in_action=0.0,
               note_goldstone_mass="m_G of the pseudo-Goldstone DM is an added explicit-breaking term (C35/C35b), not part of the common action",
               Q_mass_GeV=de["masses_GeV"][0], Q_width_GeV=de["widths_GeV"][0], min_Q_width_over_Gamma_chi=de["min_Q_width_over_Gamma_chi"],
               late_fluid_T_GeV=[min(T_late), max(T_late)], T_rh_S2_GeV=B["T_rh_GeV"],
               late_vs_S2_rel=max(abs(t / B["T_rh_GeV"] - 1) for t in T_late),
               DeltaNeff_late_homogeneous_max=dneff_late, F_identity_rel=same_F,
               ordering=[c for c in ct["connections"] if c["connection"].startswith("χ -> curvatura")][0],
               network=dict(local_certified_perimeters=cs["local_certified_perimeters"],
                            representative_core_resolved_network=cs["representative_core_resolved_network"],
                            actual_network_source_calibrated=cs["actual_network_source_calibrated"],
                            source_contract_status=sc["status"], source_blockers=sc["blockers"]),
               downstream_controls_passed=bool(e2e["numerical_free_wave_controls_passed"] and e2e["energy_and_number_matching"]))
    out["status"] = "PASS" if (out["late_vs_S2_rel"] < 1e-3 and same_F < 1e-9) else "FAIL"
    save("s3_sectors.json", out)
    print(f"S3 {out['status']}: T_late/T_rh-1 = {out['late_vs_S2_rel']:.1e}, F identity {same_F:.1e}, dNeff_late {dneff_late:.1e}, "
          f"source {sc['status']}")
    return out


if __name__ == "__main__":
    main()
