"""Configuration for the modelling (Idea 1) and surveillance-policy (Idea 2) package.

Design: docs/05_modelling_and_surveillance_plan.md (recommended options). Data: data/synthetic/.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data" / "synthetic"
SIM = ROOT / "notebooks" / "simulator"
KM = ROOT / "notebooks" / "km_reconstruct"
OUT = HERE / "output"
for pth in (str(HERE), str(SIM), str(KM)):
    if pth not in sys.path:
        sys.path.insert(0, pth)

# ---- parallelism (255 CPUs, 8 GPUs on the development machine)
N_CPU = os.cpu_count() or 1
N_JOBS = int(os.environ.get("ANALYSIS_N_JOBS", min(N_CPU, 48)))        # loky worker processes for independent runs
OMP_THREADS = int(os.environ.get("ANALYSIS_OMP_THREADS", min(N_CPU, 64)))  # OpenMP threads for main-process boosted fits
os.environ.setdefault("OMP_NUM_THREADS", str(OMP_THREADS))
BACKEND = os.environ.get("ANALYSIS_BACKEND", "hgb")                      # "hgb" (CPU, default) | "xgb_gpu"

# ---- horizons and discretisation
PERIOD_Y = 0.5
HORIZONS_Y = (5.0, 8.0, 10.0)
LANDMARK_HORIZON_Y = 5.0                    # landmark model: periods 1..10; metrics at 3 and 5 y
LANDMARK_TAUS_Y = (3.0, 5.0)
IMPLANT_HORIZON_Y = 10.0                    # implant-time model: periods 1..20
INTERVAL_CANDIDATES_Y = (0.5, 1.0, 1.5, 2.0, 3.0)   # surveillance policy candidates
SEED = 20260916
TEST_FRAC, VAL_FRAC = 0.15, 0.15
BOOT_B = 500

# ---- label definitions: name -> (event-time column in labels.csv, description)
LABELS = {
    "varc3_conf": ("varc3_conf2_confirming", "VARC-3 full, stage >=2 confirmed on two consecutive echoes (PRIMARY); known at the confirming echo"),
    "varc3_single": ("varc3_single2", "VARC-3 full, single echo"),
    "varc3_haemo_conf": ("varc3_haemo_conf2_confirming", "NOTION haemodynamic variant (gradient/AR only), confirmed"),
    "capodanno_conf": ("capodanno_conf2_confirming", "Capodanno 2017, confirmed"),
    "dvir_conf": ("dvir_conf2_confirming", "Dvir 2018, confirmed"),
    "oracle": ("oracle_varc3_single2", "oracle: VARC-3 rule on the noise-free channels at the same visits"),
    "varc3_conf3": ("varc3_conf3", "VARC-3 stage 3 confirmed"),
    "varc3_single3": ("varc3_single3", "VARC-3 stage 3 single echo (triggers reintervention in the simulator)"),
}
PRIMARY_LABEL = "varc3_conf"

# ---- feature sets (column names produced by features.py)
BASELINE = ["age_at_implant", "female", "bsa", "bmi", "diabetes", "ckd", "dialysis", "smoking", "af", "anticoag_type",
            "approach", "valve_model", "valve_class", "is_new_model", "label_size_mm", "implant_year", "reference_from_table"]
REFERENCE_ECHO = ["ref_mpg", "ref_eoa", "ref_eoai", "ref_dvi", "ref_ar", "ppm_class_ref", "expected_eoa", "eoa_residual", "ref_mpg_flow_corrected"]
LONGITUDINAL = ["landmark_t", "n_echoes", "time_since_last", "mpg_now", "eoa_now", "dvi_now", "vmax_now", "svi_now", "lvef_now", "ar_now",
                "delta_mpg", "rel_eoa", "rel_dvi", "ewma_log_mpg", "p_true_mpg_ge20", "slope_since_ref", "slope_last2", "max_mpg_sofar",
                "stage_now_full", "stage_now_haemo", "unconfirmed_flag", "prior_positive_count", "consecutive_worsening",
                "transient_spike", "symptom_now", "any_symptom_before", "eoa_missing_now", "dvi_missing_now", "stage_now_unevaluable"]
FULL = BASELINE + REFERENCE_ECHO + LONGITUDINAL
# what TVT / STS + the current echo report hold, without engineered trajectory features
REGISTRY_ONLY = BASELINE + ["ref_mpg", "ref_eoa", "ref_eoai", "ref_dvi", "ref_ar", "ppm_class_ref"] + \
                ["landmark_t", "n_echoes", "time_since_last", "mpg_now", "eoa_now", "dvi_now", "vmax_now", "lvef_now", "ar_now", "delta_mpg", "symptom_now"]
# features whose meaning does not depend on the echo spacing (used by the surveillance policy model)
SPACING_ROBUST = BASELINE + REFERENCE_ECHO + ["landmark_t", "mpg_now", "eoa_now", "dvi_now", "ar_now", "delta_mpg", "rel_eoa", "rel_dvi",
                                              "ewma_log_mpg", "p_true_mpg_ge20", "slope_since_ref", "max_mpg_sofar", "stage_now_full",
                                              "stage_now_haemo", "unconfirmed_flag", "prior_positive_count", "symptom_now", "any_symptom_before",
                                              "eoa_missing_now", "dvi_missing_now", "stage_now_unevaluable"]
CATEGORICAL = ["anticoag_type", "approach", "valve_model", "valve_class"]
# fixed ordinal codes (trees split on them; SHAP then explains exactly the splits the model uses)
CODES = {
    "anticoag_type": {"none": 0, "VKA": 1, "DOAC": 2},
    "approach": {"SAVR": 0, "TAVR": 1},
    "valve_model": {m: i for i, m in enumerate(["Perimount Magna Ease", "Inspiris Resilia", "Trifecta", "Mitroflow", "Epic", "Mosaic",
                                                 "New SAVR valve", "Sapien XT", "Sapien 3", "CoreValve", "Evolut R", "New THV"])},
    "valve_class": {c: i for i, c in enumerate(["perimount", "trifecta", "porcine", "new_savr", "balloon", "self_expanding", "new_thv"])},
}
PARAMS_CALIBRATED = SIM / "output" / "params_calibrated.json"   # the ONLY parameters that reproduce data/synthetic (seed 1)
MPG_CAP = 100.0            # reported gradients are winsorised at 100 mmHg (physically implausible beyond)
# the B4 Table 2 predictor set (published risk-factor Cox comparator)
B4_SET = ["age_at_implant", "female", "bmi", "diabetes", "smoking", "ckd", "ppm_class_ref", "ref_mpg", "label_size_mm", "approach", "anticoag_type"]

LOG_MPG_SD = 0.14          # per-echo multiplicative Doppler error (A3-calibrated simulator noise)
EWMA_TAU_Y = 1.5           # time constant of the exponentially weighted log-gradient
