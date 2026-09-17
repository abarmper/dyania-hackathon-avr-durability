"""Command line: calibrate | validate | generate.

  python cli.py calibrate [--n-sim 20000] [--quick]
  python cli.py validate  [--params output/params_calibrated.json] [--n 10000]
  python cli.py generate  [--params output/params_calibrated.json] [--n 10000] [--seed 1] [--out ../../data/synthetic]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from params import SimParams, default_params, provenance_table  # noqa: E402

DATA_DICTIONARY = {
    "patients_valves": {
        "valve_id": "synthetic key (one row per implanted valve; D7: the valve is the unit of analysis)",
        "approach": "SAVR / TAVR (TVT DCF procedure; STS ACSD `VSAVPr`)",
        "valve_model": "named prosthesis (STS `AorticImplant` / UDI `VSAoImUDI`; TVT device fields)",
        "valve_class": "calibration stratum (Perimount-class, Trifecta-class, porcine/Mitroflow, balloon-expandable, self-expanding, new)",
        "is_new_model": "1 if the model has little published history (partial pooling target)",
        "implant_year / implant_date": "index procedure (STS/TVT procedure date; shifted dates in the provided notes)",
        "age_at_implant, sex, height_cm, weight_kg, bsa, bmi": "TVT Height 6000 / Weight 6005; STS `CalculatedBMI`; BSA by Mosteller",
        "diabetes, ckd, dialysis, smoking, af": "TVT history fields (Dialysis 13880, Tobacco 4625, AFib 13179); STS `Diabetes`, `Dialysis`, `TobaccoUse`, `AFib`",
        "anticoag_discharge, anticoag_type, anticoag_year1, anticoag_longterm": "TVT discharge medications block 10205; STS `DCDirOralAnticoag`",
        "lvef": "TVT LVEF 13305/13690; STS `PreLVEF`/`PostLVEF`; MIMIC-IV-ECHO measurement",
        "label_size_mm": "labelled prosthesis size (STS implant record; TVT device size)",
        "reference_source": "30d_echo (TVT 13675/13676; baseline echo, ESC/EACTS 2025 window within 3 months) or model_size_table (D1 missing-reference arm)",
        "ref_mean_gradient_mmHg, ref_av_area_cm2, ref_dvi, ref_ar_grade, ref_eoai_cm2_m2, ppm_class_ref": "reference (30-day) echo values; TVT AV Mean Gradient 13675, AV Area 13495, Central AR 14499; PPM per ASE 2024 Table 7",
        "table_eoa, table_mpg": "ASE 2024 normal values for the model x size (reference when the 30-day echo is missing)",
        "followup_end_years": "end of valve follow-up (death, reintervention, lost, administrative censoring)",
    },
    "echo_visits": {
        "valve_id, visit_time_years, visit_date": "one row per measurement per echo (long format as MIMIC-IV-ECHO `structured_measurement`)",
        "visit_type": "reference | scheduled | symptom_triggered | other_indication",
        "measurement": "av_mean_gradient_mmHg (TVT 13674-13676), av_peak_velocity_m_s (13703), av_area_cm2 (13481/13495/13669), doppler_velocity_index, stroke_volume_index_ml_m2 (Low Flow SVI 13700), lvef_pct, intraprosthetic_ar_grade (14499/14500), paravalvular_leak_grade (14503/14504)",
        "result": "value as reported (rounded as in clinical reports; EOA/DVI missing in a share of echoes as in PERIGON)",
    },
    "events": {
        "event_type": "death | endocarditis | reintervention | lost_to_follow_up | administrative_censoring (TVT follow-up form; STS `ValExpUDI`)",
        "detail": "reintervention: type|indication|urgency (ViV / redo; severe SVD, moderate SVD with symptoms, endocarditis); death: cause class",
    },
    "labels": {
        "<definition>_single2 / _single3": "first echo meeting the definition (Capodanno, Dvir, VARC-3 full, VARC-3 haemodynamic, oracle) at stage >=2 / stage 3, years",
        "<definition>_conf2 / _conf3": "first echo of a pair of consecutive echoes both meeting the stage (persistence-confirmed)",
        "primary_event, primary_time, primary_competing": "PRIMARY ENDPOINT: confirmed VARC-3 stage >=2 (time = first echo of the confirmed pair); competing 2 = death, 3 = reintervention",
    },
    "latent_truth": {"*": "NOT FOR TRAINING. Latent onset T_deg, mode, progression tau, thrombosis episodes, noise-free crossing times, symptom onset; used only for the oracle and the surveillance-delay metrics"},
}


def write_dataset(out, out_dir: Path, p):
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("patients_valves", "echo_visits", "events", "labels", "latent_truth"):
        df = out[name]
        df.to_csv(out_dir / f"{name}.csv", index=False)
        try:
            df.to_parquet(out_dir / f"{name}.parquet", index=False)
        except Exception:
            pass
    p.to_json(out_dir / "params_used.json")
    lines = ["# Synthetic cohort data dictionary\n", "Every column exists in a named real-world source (docs/04 section 3); values are simulated, the schema is not.\n"]
    for table, cols in DATA_DICTIONARY.items():
        lines += [f"\n## `{table}`\n", "| column(s) | meaning / real-world field |", "|---|---|"]
        lines += [f"| `{c}` | {m} |" for c, m in cols.items()]
    lines += ["\n## Parameter provenance\n", provenance_table()]
    (out_dir / "data_dictionary.md").write_text("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("calibrate"); c.add_argument("--n-sim", type=int, default=20000); c.add_argument("--seed", type=int, default=0); c.add_argument("--quick", action="store_true")
    v = sub.add_parser("validate"); v.add_argument("--params", default=str(HERE / "output" / "params_calibrated.json")); v.add_argument("--n", type=int, default=10000); v.add_argument("--seed", type=int, default=1); v.add_argument("--n-matched", type=int, default=20000)
    g = sub.add_parser("generate"); g.add_argument("--params", default=str(HERE / "output" / "params_calibrated.json")); g.add_argument("--n", type=int, default=10000); g.add_argument("--seed", type=int, default=1); g.add_argument("--out", default=str(HERE.parent.parent / "data" / "synthetic"))
    a = ap.parse_args(argv)
    if a.cmd == "calibrate":
        from calibrate import run_all; run_all(a.n_sim, a.seed, a.quick)
    elif a.cmd == "validate":
        from validate import run_all
        p = SimParams.from_json(a.params) if Path(a.params).exists() else default_params()
        run_all(p, a.n, a.seed, a.n_matched)
    else:
        from physics import fit_physics_constants
        from simulate import simulate_cohort
        p = SimParams.from_json(a.params) if Path(a.params).exists() else default_params()
        if not p.physics.c_model: fit_physics_constants(p)
        out = simulate_cohort(p, n=a.n, seed=a.seed)
        write_dataset(out, Path(a.out), p)
        print(f"wrote {a.out}: " + ", ".join(f"{k} {len(out[k])} rows" for k in ("patients_valves", "echo_visits", "events", "labels", "latent_truth")))


if __name__ == "__main__":
    main()
