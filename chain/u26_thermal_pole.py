"""U26 (docs/preregistros/ESPEC_U26.md), point 3: the frozen pivot matching of the pole plateau, data/cadena_cierre/thermal_pole.json.
Same continuous matching as data/cadena_cierre/thermal.json (rtol = 2e-10, convergence against 2e-8), with the U22 swap of potential,
predict and m_chi (chain/u22_pole.py, s2_pole).  Run once; data/MANIFEST.json records it, and action/vendor/ gets a copy.
    python -B chain/u26_thermal_pole.py      ->  data/cadena_cierre/thermal_pole.json, action/vendor/thermal_pole.json"""
import json, shutil
from ucommon import ROOT, INPUTS, sha, _conv
import u22_pole as U22

OUT = INPUTS / "cadena_cierre" / "thermal_pole.json"


def main():
    a, b = U22.s2_pole(U22.C_FIT, 2e-8), U22.s2_pole(U22.C_FIT, 2e-10)
    conv = {k: abs(a[k] / b[k] - 1) for k in ("N_star", "V0")}
    conv["T_final_GeV"] = abs(a["thermal_history"]["T_final_GeV"] / b["thermal_history"]["T_final_GeV"] - 1)
    out = dict(potential=f"order-4 pole, V = V0 [1 - C/sqrt(chi^2 + C^2)], C = {U22.C_FIT} (ESPEC_U22, ESPEC_U26)", continuous=b, convergence=conv,
               pass_energy_and_convergence=bool(b["thermal_history"]["energy_balance_relative"] < 1e-6 and max(conv.values()) < 1e-3),
               assumption="as data/cadena_cierre/thermal.json: exact homogeneous GR scalar endpoint epsilon_H=1; slow-roll CMB spectra; "
                          "subsequent averaged chi fluid w=0, g*=gs=106.75; no lattice preheating or Q backreaction; m_chi^2 = V0/C^2")
    OUT.write_text(json.dumps(out, indent=1, default=_conv), encoding="utf-8")
    shutil.copyfile(OUT, ROOT / "action" / "vendor" / "thermal_pole.json")
    man = INPUTS / "MANIFEST.json"; m = json.loads(man.read_text())
    m["files"]["cadena_cierre/thermal_pole.json"] = dict(origin="generated in this repository by chain/u26_thermal_pole.py (U26)", sha256=sha(OUT))
    man.write_text(json.dumps(m, indent=1), encoding="utf-8")
    print("thermal_pole:", {k: b[k] for k in ("N_star", "n_s", "r", "V0", "T_rh_GeV", "F_GeV", "m_chi_GeV")}, "conv", conv, "pass", out["pass_energy_and_convergence"])


if __name__ == "__main__":
    main()
