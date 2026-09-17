"""All simulator parameters, each with a source tag or an ASSUMPTION flag (see PROVENANCE).

Reference codes are the files in references/ (README there): A2 ESC 2021 (superseded by ESC/EACTS 2025, see references/README; still the source of the HALT figures), A3 Velders 2024, A4 Kermen 2022,
A5 Zoghbi 2024 (ASE), B1 NOTION 10-y, B4 Trimaille 2025 (JAHA review), B5 Johnston 2015, B7 Geisler 2026,
C1 Wakami 2022. Values marked CALIBRATED are overwritten by calibrate.py (params_calibrated.json).
"""
from __future__ import annotations

import dataclasses
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

CLASSES = ["perimount", "trifecta", "porcine", "new_savr", "balloon", "self_expanding", "new_thv"]


@dataclass
class Cohort:
    n: int = 10_000
    savr_frac: float = 0.55
    implant_year_min: float = 2010.0
    implant_year_max: float = 2021.0          # uniform on [2010-01, 2021-01)
    censor_year: float = 2026.71              # 2026-09-15
    seed: int = 1


@dataclass
class Baseline:
    age_savr: tuple = (70.0, 9.0, 40.0, 90.0)     # mean, sd, lo, hi
    age_tavr: tuple = (78.0, 6.0, 60.0, 95.0)
    female_frac: float = 0.45
    height_cm: dict = field(default_factory=lambda: {"M": (175.0, 7.0), "F": (162.0, 7.0)})
    bmi: tuple = (29.0, 5.0, 16.0, 55.0)
    diabetes: float = 0.27
    ckd: float = 0.11
    dialysis_given_ckd: float = 0.12
    smoking: float = 0.15
    af: dict = field(default_factory=lambda: {"SAVR": 0.35, "TAVR": 0.25})
    anticoag_if_af: float = 0.85
    anticoag_if_no_af: float = 0.10
    doac_from_year: float = 2015.0
    postop_vka_savr: float = 0.20             # 3-month VKA course after SAVR in some centres
    lvef: tuple = (59.0, 10.0, 20.0, 80.0)
    svi: tuple = (42.0, 4.0, 30.0, 58.0)      # mL/m2 (true between-patient spread; visit-to-visit flow variation is in Noise)
    svi_per_lvef_point: float = 0.20          # mild dependence of SVi on LVEF
    ejection_time_s: tuple = (0.30, 0.015)
    q_clip: tuple = (150.0, 450.0)            # mL/s
    reference_missing_frac: float = 0.30      # D1(a): no 30-day echo -> reference from the model x size table
    ar_ref_probs: tuple = (0.85, 0.14, 0.01)  # none / trace-mild / >= moderate intraprosthetic AR at reference
    pvl_ref_probs: dict = field(default_factory=lambda: {"SAVR": (0.95, 0.05, 0.0), "TAVR_early": (0.50, 0.40, 0.10), "TAVR_late": (0.65, 0.32, 0.03)})
    lvot_bias_sd: float = 0.06                # per-valve persistent LVOT measurement bias (lognormal sd)


@dataclass
class Physics:
    gorlin_c: float = 1.0 / 44.3 ** 2
    c_model: dict = field(default_factory=dict)      # CALIBRATED (fit to ASE tables) in physics.fit_physics_constants
    c_size: dict = field(default_factory=dict)       # per "model|size" constant (exact to the ASE row); c_model is the fallback
    peak_over_mean: tuple = (1.75, 0.10)
    eoa_floor: float = 0.4
    q_ref_ml_s: float = 280.0                        # flow at which EOA equals the tabulated EOA
    eoa_flow_exponent: float = 0.35                  # bioprosthetic EOA opens with flow: EOA_eff = EOA (Q/q_ref)^gamma (ASSUMPTION)


@dataclass
class Onset:
    # per-class Weibull baseline of LATENT onset; calibrate.py fits (scale, shape) per class, shift stays 0
    shape: dict = field(default_factory=lambda: {"perimount": 4.0, "trifecta": 1.6, "porcine": 2.0, "new_savr": 4.0,
                                                 "balloon": 2.0, "self_expanding": 2.5, "new_thv": 2.5})
    scale: dict = field(default_factory=lambda: {"perimount": 10.5, "trifecta": 22.0, "porcine": 12.0, "new_savr": 10.5,
                                                 "balloon": 30.0, "self_expanding": 16.0, "new_thv": 16.0})
    shift: dict = field(default_factory=lambda: {"perimount": 0.0, "trifecta": 0.0, "porcine": 0.0, "new_savr": 0.0,
                                                 "balloon": 0.0, "self_expanding": 0.0, "new_thv": 0.0})   # location shift (years), kept at 0
    calibration_age: dict = field(default_factory=lambda: {"perimount": 70.6, "trifecta": 78.0, "porcine": 79.0, "new_savr": 70.6,
                                                           "balloon": 78.0, "self_expanding": 79.2, "new_thv": 79.2})
    hr_age_per_year_savr: float = 0.951       # A4 Table 3 (stage 2/3)
    hr_age_per_year_thv: float = 0.97         # ASSUMPTION (B4: inconclusive)
    hr_smoking: float = 3.5                   # latent input chosen so that the EMERGENT Cox HR on observed confirmed labels reproduces B4's 2.58 (label noise dilutes effects)
    hr_diabetes: float = 1.33                 # B4 Table 2
    hr_bmi_per_unit: float = 1.08             # B4 Table 2, centred at bmi_centre
    bmi_centre: float = 29.0
    hr_ckd: float = 1.30                      # A4 1.29-1.34 / B4 1.10
    hr_no_anticoag_year1_tavr: float = 1.6    # biological residual of B4 3.35 (rest via thrombosis module)
    hr_vka_longterm_savr: float = 1.2         # ASSUMPTION (Salaun late factor, B4)
    hr_severe_ppm_residual: float = 1.3       # residual; total B4 1.85 is an emergent check
    hr_ref_mpg_ge15_residual: float = 1.6     # biological residual for a true reference gradient >= 15 mmHg (B5: higher gradient -> faster explant; B4 total 1.30 is the emergent check; the observed-reference regression to the mean works against it)
    hr_thv_le23_residual: float = 2.5         # residual; total B4 2.07 is an emergent check (small THVs also have high reference gradients, which regression to the mean turns protective on observed labels)
    hr_prior_thrombosis: float = 1.5          # ASSUMPTION (B4: HALT associated with SVD)
    hr_female: float = 1.0                    # fairness null
    p_regurgitant: dict = field(default_factory=lambda: {"perimount": 0.25, "trifecta": 0.60, "porcine": 0.30, "new_savr": 0.25,
                                                         "balloon": 0.05, "self_expanding": 0.05, "new_thv": 0.05})   # THV: NOTION severe SVD 1.5% at 10 y, intraprosthetic AR 5.8% mild/moderate


@dataclass
class Trajectory:
    drift_per_year: float = 0.011             # EOA loss/yr before onset -> MPG 14 -> 20 by year 15 (B5)
    # progression time constant after onset (EOA* decays as exp(-u/tau)); per class, ASSUMPTION checked on the
    # stage-3 curves (Kermen stage 3 12.8% vs stage 2/3 39% at 10 y -> tau ~ 7.5 y for Perimount-class). The
    # moderate/severe split of supra-annular THVs (NOTION: 15.4% vs 1.5%) comes from physics: from an 8 mmHg
    # baseline, stage 3 (>= 30 mmHg, +20) needs an almost fourfold gradient rise
    tau_median_y: dict = field(default_factory=lambda: {"perimount": 7.5, "trifecta": 6.0, "porcine": 7.5, "new_savr": 7.5,
                                                        "balloon": 10.0, "self_expanding": 10.0, "new_thv": 10.0})
    tau_sigma: float = 0.4
    # fraction of the orifice area eventually lost after onset (EOA* -> EOA_ref (1 - L)); surgical calcific
    # degeneration runs to completion, transcatheter deterioration mostly plateaus at a partial loss
    # (B1: NOTION TAVI 15.4% >= moderate but 1.5% severe at 10 y; ASSUMPTION on the spread)
    loss_fraction: dict = field(default_factory=lambda: {"perimount": (0.92, 0.06), "trifecta": (0.92, 0.06), "porcine": (0.92, 0.06), "new_savr": (0.92, 0.06),
                                                         "balloon": (0.40, 0.06), "self_expanding": (0.40, 0.06), "new_thv": (0.40, 0.06)})
    ar_severe_after_median_y: float = 2.0     # regurgitant mode: moderate at onset, severe after LN(median, sigma)
    ar_severe_after_sigma: float = 0.5
    regurgitant_eoa_factor: float = 0.95      # mild EOA change in regurgitant mode
    flow_slope_sd_per_year: float = 0.03      # per-valve slow drift of flow state (log Q per year); ASSUMPTION targeting A3's 5-y change PI


@dataclass
class Thrombosis:
    on: bool = True
    p_year1: dict = field(default_factory=lambda: {"SAVR": 0.03, "TAVR": 0.08})   # A2 / B4
    anticoag_factor: float = 0.5
    onset_window_y: tuple = (0.05, 1.0)
    duration_y: tuple = (0.25, 0.75)          # B4: resolves with 3-6 months of anticoagulation
    eoa_factor: float = 0.75                  # ASSUMPTION (+~10 mmHg at 13 mmHg)


@dataclass
class Noise:
    sigma_flow: float = 0.08                  # CALIBRATED (A3 decile table); per-visit flow multiplier, MPG x f^2
    sigma_lvot_vti: float = 0.08              # CALIBRATED; EOA and DVI x e_L
    sigma_cw: float = 0.06                    # CALIBRATED; MPG x e_A^2, EOA and DVI / e_A
    sigma_lvot_diam: float = 0.05             # CALIBRATED; EOA x e_D^2
    sigma_flow_ref_extra: float = 0.08        # extra flow variability at the 30-day reference echo (post-operative state; A3 Table E4: 30 d -> 3-6 mo interval has the largest regression to the mean)
    flow_ref_bias: float = 0.05               # hyperdynamic post-operative state: log-flow bias at the reference echo (ASSUMPTION)
    ar_misclass: float = 0.03                 # +/-1 grade; A3: >= moderate regurgitation 0.2% at 5 y, so misgrading must be rare
    eoa_missing: float = 0.08                 # A3 Table E1 (EOA missing more often than MPG)
    dvi_missing: float = 0.05
    lvef_missing: float = 0.05


@dataclass
class Visits:
    reference_time_y: float = 0.08            # 30-day echo (A2, D1a)
    max_years: int = 16
    jitter_sd_y: float = 0.125
    miss_by_year: tuple = (0.093, 0.145, 0.224, 0.350)   # A3 Table E1, years 1-4
    miss_saturation: float = 0.45             # ASSUMPTION beyond year 4
    miss_mult_well: float = 1.3               # ASSUMPTION
    miss_mult_symptomatic: float = 0.6        # ASSUMPTION
    symptom_h0: float = 0.25                  # ASSUMPTION: h = h0 * (MPG*/20)^2 per year when latent stage >= 2
    symptom_power: float = 2.0
    background_echo_rate: float = 0.02        # non-valve reasons for an echo, per year
    symptom_echo_delay_y: float = 0.04        # ~2 weeks
    max_symptom_echoes: int = 3


@dataclass
class Death:
    a: float = 0.020                          # CALIBRATED to Kermen overall survival (A4 Fig 2A)
    b_per_year: float = math.log(1.085)       # population Gompertz slope in attained age (ASSUMPTION; Kermen adjusted 1.033)
    age_ref: float = 75.0
    hr_lvef_per_point: float = 0.974          # A4 Table E2, centred at 60
    hr_dialysis: float = 2.0                  # ASSUMPTION
    hr_tavr_era: float = 1.8                  # ASSUMPTION: TAVR implanted <= 2018 = higher-risk population (NOTION STS 3.0, 63% dead at 10 y)
    tavr_era_until: float = 2018.0
    frailty_var: float = 0.6                  # gamma frailty (mean 1) on the death hazard; flattens late hazard (ASSUMPTION)


@dataclass
class Events:
    endocarditis_rate: float = 0.0075         # B1: 7.2% at 10 y
    endocarditis_explant_frac: float = 0.60
    endocarditis_death_frac: float = 0.25
    reint_h_stage2_symptomatic: float = 0.05  # per year, CALIBRATED
    reint_h_stage3: float = 1.50              # per year at age 70 (Kermen: 9 of 11 stage-3 patients reintervened, explant CIF 8.4% at 10 y)
    reint_age_slope: float = -0.30            # x exp(slope * (age - 70)), clipped to [0.05, 3]/y; NOTION at 79: 2-4% reintervention despite 10% severe SVD
    viv_logistic: tuple = (75.0, 0.25)        # P(ViV) = logistic((age - 75) * 0.25); ViV 79.0 +/- 5.3 vs redo 70.6 +/- 8.8 (B7)
    viv_min_year: float = 2012.0
    ltfu_rate: float = 0.01                   # per year


@dataclass
class SimParams:
    cohort: Cohort = field(default_factory=Cohort)
    baseline: Baseline = field(default_factory=Baseline)
    physics: Physics = field(default_factory=Physics)
    onset: Onset = field(default_factory=Onset)
    trajectory: Trajectory = field(default_factory=Trajectory)
    thrombosis: Thrombosis = field(default_factory=Thrombosis)
    noise: Noise = field(default_factory=Noise)
    visits: Visits = field(default_factory=Visits)
    death: Death = field(default_factory=Death)
    events: Events = field(default_factory=Events)

    def to_json(self, path):
        Path(path).write_text(json.dumps(dataclasses.asdict(self), indent=1))

    @classmethod
    def from_json(cls, path):
        d = json.loads(Path(path).read_text())
        p = cls()
        for group, vals in d.items():
            g = getattr(p, group)
            for k, v in vals.items():
                cur = getattr(g, k)
                setattr(g, k, tuple(v) if isinstance(cur, tuple) else v)
        return p

    def set(self, dotted: str, value):
        """p.set('noise.sigma_flow', 0.12) or p.set('onset.scale.perimount', 11.0)"""
        parts = dotted.split(".")
        obj = getattr(self, parts[0])
        if len(parts) == 2:
            setattr(obj, parts[1], value)
        else:
            getattr(obj, parts[1])[parts[2]] = value

    def get(self, dotted: str):
        parts = dotted.split(".")
        obj = getattr(getattr(self, parts[0]), parts[1])
        return obj[parts[2]] if len(parts) == 3 else obj


def default_params() -> SimParams:
    return SimParams()


# (parameter, source, assumption?) -- rendered into the data dictionary / data plan
PROVENANCE = [
    ("cohort.*", "docs/04 D5: 10,000 valves, SAVR:TAVR 55:45, implants 2010-2020, censor Sept 2026", False),
    ("baseline.age_*", "A3 Table 1 (70.2 +/- 9.0), B1 (79.1 +/- 4.8), B5", False),
    ("baseline.height/bmi", "A3 Table 1: BSA 2.0 +/- 0.2, BMI 29.4 +/- 5.4", False),
    ("baseline.diabetes/ckd/smoking", "A3 Table 1: 27% / 11% renal / smoking prevalence ASSUMED 15%", True),
    ("baseline.af, anticoag", "B1 new-onset AF 52-74%; A4 5.6% at implant; anticoagulation via AF ASSUMED shares", True),
    ("baseline.lvef", "A3 Table 1: 59 +/- 10", False),
    ("baseline.svi, ejection_time", "physiological ranges; low-flow < 35 mL/m2 ~16% (ASSUMPTION)", True),
    ("baseline.reference_missing_frac", "docs/04 D1(a) missing-reference arm 30% (ASSUMPTION; ~5/117 numeric TEE gradients in notes)", True),
    ("baseline.ar_ref/pvl_ref", "B1 (PVL 18-53% TAVI, 5.2% SAVR; intraprosthetic 5.8 / 2.2%), A3 (0.2% >= moderate at 5 y)", False),
    ("physics.c_model, A_eff", "A5 Appendix Tables A1/A2/A4 (fitted); Gorlin 44.3 constant", False),
    ("physics.peak_over_mean", "A5 Table A4 (Perimount 19.9/11.5, Mosaic 23.8/13.7)", False),
    ("onset.shape/scale", "CALIBRATED to km_reconstruct/output/targets.json (Kermen, Wakami, NOTION) and B4 point targets", False),
    ("onset.hr_age_per_year_savr", "A4 Table 3: 0.951 (stage 2/3)", False),
    ("onset.hr_age_per_year_thv", "ASSUMPTION (B4 Table 2: inconclusive)", True),
    ("onset.hr_smoking/diabetes/bmi/ckd", "B4 Table 2: 2.58 / 1.33 / 1.08 per unit / 1.10 (A4 1.29-1.34)", False),
    ("onset.hr_no_anticoag_year1_tavr", "biological residual of B4 3.35 (rest through thrombosis module)", True),
    ("onset.hr_vka_longterm_savr", "ASSUMPTION from Salaun late-HVD factor (B4)", True),
    ("onset.hr_severe_ppm_residual / hr_thv_le23_residual", "residuals; totals 1.85 / 2.07 (B4) are emergent checks", True),
    ("onset.p_regurgitant", "C1: 5/7 Trifecta failures were cusp tears; others ASSUMED", True),
    ("trajectory.drift_per_year", "B5: mean gradient 14 -> 20 mmHg by year 15", False),
    ("trajectory.tau_*", "ASSUMPTION; checked against Kermen stage 3 and NOTION severe curves", True),
    ("thrombosis.p_year1, anticoag_factor, duration", "A2 (HALT 12.4% vs 32.4%), B4 (SLT 6-15%, resolves in 3-6 months)", False),
    ("thrombosis.eoa_factor", "ASSUMPTION (no published gradient effect)", True),
    ("noise.sigma_*", "CALIBRATED to A3 Table E4 (decile regression to the mean) and 5-y change PI (-9.6, +7.5)", False),
    ("noise.eoa_missing/dvi_missing", "A3 Table E1", False),
    ("visits.miss_by_year", "A3 Table E1 (9.3 / 14.5 / 22.4 / 35.0 %)", False),
    ("visits.miss_saturation, miss_mult_*", "ASSUMPTION", True),
    ("visits.symptom_*", "ASSUMPTION (most stage-2 disease silent)", True),
    ("death.a", "CALIBRATED to A4 Fig 2A overall survival", False),
    ("death.b_per_year", "population Gompertz slope ~9%/yr (ASSUMPTION; A4 adjusted within-cohort HR 1.033)", True),
    ("death.hr_lvef_per_point", "A4 Table E2: 0.974", False),
    ("death.hr_dialysis, hr_tavr_era", "ASSUMPTION (NOTION 63% dead at 10 y is the check)", True),
    ("events.endocarditis_rate", "B1: 7.2-7.4% at 10 y", False),
    ("events.reint_*", "CALIBRATED to A4 explant CIF 8.4% at 10 y, 9/11 stage-3 reintervened; B1 2.2-4.3%", False),
    ("events.viv_logistic", "B7: ViV 79.0 +/- 5.3 vs redo 70.6 +/- 8.8 y", False),
    ("events.ltfu_rate", "A4: 2% lost over 6.6 y", False),
]


def provenance_table() -> str:
    lines = ["| parameter | source | assumption |", "|---|---|---|"]
    for name, src, asm in PROVENANCE:
        lines.append(f"| `{name}` | {src} | {'yes' if asm else ''} |")
    return "\n".join(lines) + "\n"
