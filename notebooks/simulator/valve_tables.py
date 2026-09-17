"""Reference haemodynamics by valve model and labelled size, valve classes, and paper cohort profiles.

Normal-value rows are transcribed from Zoghbi et al., ASE 2024 guideline (references/A5), Appendix
Tables A1 (SAPIEN family), A2 (CoreValve / Evolut R) and A4 (surgical valves): mean +/- SD of EOA (cm2),
mean gradient (mmHg) and, for THVs, DVI. Where a size is not in the table the value is extrapolated
from the neighbouring sizes of the same model (flagged EXTRAPOLATED). Mitroflow has only the 19 mm row
in the guideline; sizes 21-25 follow the Perimount size gradient (flagged).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Row:
    eoa: float
    eoa_sd: float
    mpg: float
    mpg_sd: float
    dvi: float | None = None
    dvi_sd: float | None = None
    note: str = ""


# model -> class. Classes are the calibration strata (one published curve or point target per class).
MODEL_CLASS = {
    "Perimount Magna Ease": "perimount",
    "Perimount (classic)": "perimount",     # reference rows only (A4); not drawn in the cohort
    "Inspiris Resilia": "perimount",        # new-generation Perimount-class; flagged "new" (little history)
    "Trifecta": "trifecta",
    "Mitroflow": "porcine",                 # class name = "porcine/early-failure" stented valves (Mitroflow is pericardial)
    "Epic": "porcine",
    "Mosaic": "porcine",
    "New SAVR valve": "new_savr",
    "Sapien XT": "balloon",
    "Sapien 3": "balloon",
    "CoreValve": "self_expanding",
    "Evolut R": "self_expanding",
    "New THV": "new_thv",
}
CLASS_PARENT = {"new_savr": "perimount", "new_thv": "self_expanding"}
SAVR_MODELS = ["Perimount Magna Ease", "Inspiris Resilia", "Trifecta", "Mitroflow", "Epic", "Mosaic", "New SAVR valve"]
TAVR_MODELS = ["Sapien XT", "Sapien 3", "CoreValve", "Evolut R", "New THV"]
SUPRA_ANNULAR = {"Trifecta", "CoreValve", "Evolut R"}

# ASE 2024 Appendix Table A4 (surgical). Perimount rows = "Baxter Perimount, stented bovine pericardial".
# Kermen's Magna Ease cohort (discharge MPG 12.6 +/- 3.0, EOA 1.6 +/- 0.3, sizes mostly 21/23) is
# reproduced by these rows (acceptance test in validate.py).
EOA_TABLE: dict[str, dict[int, Row]] = {
    # Magna Ease is not in the ASE table. Rows = the Perimount (Baxter) A4 rows scaled so that Kermen's Magna Ease
    # cohort (A4: EOA 1.6 +/- 0.3, EOAi 0.87 +/- 0.17, MPG 12.6 +/- 3.0 with sizes 19/21/23/25 = 19/37/32/13 %) is
    # reproduced; the classic Perimount rows (19: 1.3, 21: 1.3, 23: 1.6, 25: 1.6 cm2) give EOAi 0.74 for that mix.
    "Perimount Magna Ease": {
        19: Row(1.35, 0.25, 16.0, 4.5, note="Perimount A4 row scaled to Kermen"), 21: Row(1.50, 0.28, 13.5, 4.0, note="scaled to Kermen"),
        23: Row(1.70, 0.30, 11.5, 3.5, note="scaled to Kermen"), 25: Row(1.90, 0.35, 10.0, 3.5, note="scaled to Kermen"),
        27: Row(2.10, 0.40, 9.0, 3.5, note="EXTRAPOLATED"),
    },
    "Perimount (classic)": {
        19: Row(1.3, 0.2, 19.5, 5.5), 21: Row(1.3, 0.3, 13.8, 4.0), 23: Row(1.6, 0.3, 11.5, 3.9),
        25: Row(1.6, 0.4, 10.7, 3.8), 27: Row(2.0, 0.4, 4.8, 2.2, note="A4 row; 27 mm MPG unusually low"),
    },
    "Inspiris Resilia": {
        19: Row(1.1, 0.2, 17.6, 7.8), 21: Row(1.3, 0.3, 12.6, 4.7), 23: Row(1.6, 0.4, 10.1, 3.8),
        25: Row(1.8, 0.5, 9.6, 5.2), 27: Row(2.2, 0.6, 8.2, 3.5),
    },
    "Trifecta": {
        19: Row(1.41, 0.24, 10.7, 4.6), 21: Row(1.63, 0.29, 8.1, 3.5), 23: Row(1.81, 0.30, 7.2, 2.8),
        25: Row(2.02, 0.32, 6.2, 2.7), 27: Row(2.20, 0.20, 4.8, 2.0), 29: Row(2.35, 0.22, 4.7, 1.6),
    },
    "Mitroflow": {
        19: Row(1.1, 0.2, 13.1, 3.3),
        21: Row(1.2, 0.25, 11.0, 3.5, note="EXTRAPOLATED from 19 mm with Perimount size gradient"),
        23: Row(1.4, 0.3, 9.5, 3.5, note="EXTRAPOLATED"), 25: Row(1.6, 0.35, 8.5, 3.5, note="EXTRAPOLATED"),
        27: Row(1.8, 0.4, 7.5, 3.5, note="EXTRAPOLATED"),
    },
    "Epic": {
        19: Row(0.9, 0.3, 22.0, 8.0, note="EXTRAPOLATED"), 21: Row(1.0, 0.3, 19.1, 8.2), 23: Row(1.4, 0.5, 13.9, 6.0),
        25: Row(1.5, 0.5, 12.1, 5.1), 27: Row(1.6, 0.4, 11.4, 4.1), 29: Row(2.4, 1.1, 7.5, 3.3),
    },
    "Mosaic": {
        19: Row(1.2, 0.4, 16.0, 5.5, note="EXTRAPOLATED"), 21: Row(1.4, 0.4, 13.3, 5.3), 23: Row(1.6, 0.5, 11.8, 4.9),
        25: Row(1.8, 0.5, 10.6, 4.4), 27: Row(2.0, 0.5, 9.1, 4.0), 29: Row(2.3, 0.6, 8.6, 2.9),
    },
    # ASE Table A1 (SAPIEN family) and A2 (CoreValve / Evolut R), with DVI
    "Sapien XT": {
        23: Row(1.41, 0.30, 10.41, 3.74, 0.52, 0.10), 26: Row(1.74, 0.42, 9.24, 3.57, 0.54, 0.11),
        29: Row(2.06, 0.52, 8.36, 3.14, 0.53, 0.11),
    },
    "Sapien 3": {
        20: Row(1.22, 0.22, 16.23, 5.01, 0.42, 0.07), 23: Row(1.45, 0.26, 12.79, 4.65, 0.43, 0.08),
        26: Row(1.74, 0.35, 10.59, 3.88, 0.43, 0.09), 29: Row(1.89, 0.37, 9.28, 3.16, 0.40, 0.09),
    },
    "CoreValve": {
        23: Row(1.12, 0.36, 14.43, 5.72, 0.44, 0.09), 26: Row(1.74, 0.49, 8.27, 3.82, 0.59, 0.15),
        29: Row(1.97, 0.53, 8.85, 4.17, 0.54, 0.12), 31: Row(2.15, 0.72, 9.55, 3.44, 0.49, 0.12),
    },
    "Evolut R": {
        23: Row(1.09, 0.26, 14.97, 7.15, 0.42, 0.04), 26: Row(1.69, 0.40, 7.53, 2.65, 0.61, 0.13),
        29: Row(1.97, 0.54, 7.85, 3.08, 0.59, 0.14), 34: Row(2.60, 0.75, 6.30, 3.23, 0.58, 0.15),
    },
}
# "new" models have no published table: they borrow their parent's rows (flagged in provenance)
EOA_TABLE["New SAVR valve"] = EOA_TABLE["Inspiris Resilia"]
EOA_TABLE["New THV"] = EOA_TABLE["Evolut R"]

# effective LVOT area (cm2) used for DVI = EOA / A_eff, from the ASE all-size EOA / DVI ratios
# (Sapien 3 1.66/0.43, XT 1.67/0.53, CoreValve 1.88/0.55, Evolut R 2.01/0.59); surgical valves from
# Velders (A3) EOA 1.54 / DVI 0.49 = 3.14 cm2, i.e. an LVOT diameter of 2.0 cm
A_EFF = {"Sapien XT": 3.15, "Sapien 3": 3.86, "CoreValve": 3.42, "Evolut R": 3.41, "New THV": 3.41}
A_EFF_SAVR_DIAM_CM = (2.0, 0.15)     # per-valve LVOT diameter N(mean, sd) for surgical valves

# size distributions (spec 2.1), used when no paper profile is given
SIZE_PROBS = {
    "SAVR": {19: 0.15, 21: 0.35, 23: 0.30, 25: 0.15, 27: 0.05},
    "TAVR_balloon": {23: 0.30, 26: 0.45, 29: 0.25},
    "TAVR_self_expanding": {23: 0.10, 26: 0.40, 29: 0.40, 34: 0.10},
}
MODEL_PROBS_SAVR = {"Perimount Magna Ease": 0.45, "Inspiris Resilia": 0.10, "Trifecta": 0.25,
                    "Mitroflow": 0.06, "Epic": 0.05, "Mosaic": 0.04, "New SAVR valve": 0.05}
# TAVR: balloon 55% (Sapien XT if implanted <= 2015 else Sapien 3), self-expanding 40% (CoreValve <= 2014 else Evolut R), new 5%
MODEL_PROBS_TAVR = {"balloon": 0.55, "self_expanding": 0.40, "New THV": 0.05}

# published point targets for classes without a full curve (5-y observed-label SVD, %)
POINT_TARGETS = {
    "Sapien 3": {"t": 5, "svd_pct": 0.7, "source": "B4 (PARTNER 3 / Pibarot 2020, VARC-3)"},
    "Sapien XT": {"t": 5, "svd_pct": 1.6, "source": "B4 (PARTNER 2A, VARC-3)"},
    "Mitroflow": {"t": 5, "svd_pct": 8.0, "source": "A1 Table 1 (Senage 2014: freedom from SVD 92% at 5 y)"},
}


@dataclass(frozen=True)
class CohortProfile:
    """Baseline mix of a published cohort, used to simulate a matched cohort for calibration."""
    name: str
    approach: str
    models: dict            # model -> prob
    age_mean: float
    age_sd: float
    sizes: dict             # size -> prob
    reference_time_y: float # timing of the reference echo used by the paper's definition
    schedule_y: tuple       # scheduled echo times after the reference
    female_frac: float = 0.45
    diabetes: float = 0.27
    ckd: float = 0.11
    smoking: float = 0.15
    af: float = 0.30
    n_published: int = 0
    death_hr: float = 1.0   # population risk multiplier on the death hazard (NOTION: 63-64% dead at 10 y in BOTH arms, STS 3.0)
    notes: str = ""


# scheduled-visit missingness multiplier for trial / single-centre cohorts (PERIGON completed 97-99% of scheduled
# visits among survivors; Kermen follow-up 98% complete; NOTION 10-y echo missing in 14-23%); real-world cohort = 1.0
MISS_SCALE = {"perigon": 0.10, "kermen": 0.20, "wakami": 0.40, "notion_tavi": 0.50, "notion_savr": 0.50, "partner_xt": 0.30, "partner_s3": 0.30}

COHORT_PROFILES = {
    # Kermen 2022 (A4): Magna Ease, n=338, age 70.6 +/- 11.5, sizes 19/21/23/25 = 18.6/36.7/32.0/12.7 %,
    # diabetes 26%, AF 5.6%, LVEF 61.9; annual echo; VARC-3 staged vs discharge echo (single echo)
    "kermen": CohortProfile("kermen", "SAVR", {"Perimount Magna Ease": 1.0}, 70.6, 11.5,
                            {19: 0.186, 21: 0.367, 23: 0.320, 25: 0.127}, 0.08, tuple(range(1, 13)),
                            female_frac=0.453, diabetes=0.26, af=0.056, n_published=338,
                            notes="A4 Table 1; KM with death censored; annual echo"),
    # Wakami 2022 (C1): Trifecta, n=110, age ~78, BSA 1.54 (Japanese cohort), sizes 19/21/23/25 = 43/40/11/6 %
    "wakami": CohortProfile("wakami", "SAVR", {"Trifecta": 1.0}, 78.0, 5.5,
                            {19: 0.43, 21: 0.40, 23: 0.11, 25: 0.06}, 0.08, tuple(range(1, 9)),
                            female_frac=0.54, n_published=110, notes="C1 Table 1; BSA 1.54 +/- 0.18; CIF with death competing"),
    # NOTION (B1): age 79.1 +/- 4.8, 47% female, STS 3.0; 3-month reference; annual echo to 10 y
    "notion_tavi": CohortProfile("notion_tavi", "TAVR", {"CoreValve": 1.0}, 79.2, 4.9,
                                 {26: 0.45, 29: 0.45, 31: 0.10}, 0.25, tuple(range(1, 11)),
                                 female_frac=0.47, n_published=139, death_hr=1.8, notes="B1; haemodynamic VARC-3 vs 3-month echo; AJ"),
    # NOTION SAVR valve mix is NOT reported; assumed 50% Perimount-class, 25% Trifecta, 25% porcine/Mitroflow (flagged)
    "notion_savr": CohortProfile("notion_savr", "SAVR", {"Perimount Magna Ease": 0.50, "Trifecta": 0.25, "Mitroflow": 0.15, "Epic": 0.10},
                                 79.0, 4.7, {19: 0.15, 21: 0.40, 23: 0.35, 25: 0.10}, 0.25, tuple(range(1, 11)),
                                 female_frac=0.47, n_published=135, death_hr=1.8, notes="B1; valve mix ASSUMED; same mortality as the TAVI arm (64% at 10 y)"),
    # PERIGON / Velders (A3): Avalus (Perimount-class physics), n=1118, age 70.2 +/- 9.0, 25% female, BSA 2.0,
    # discharge reference, visits 3-6 mo, 1-5 y
    "perigon": CohortProfile("perigon", "SAVR", {"Perimount Magna Ease": 1.0}, 70.2, 9.0,
                             {19: 0.12, 21: 0.33, 23: 0.35, 25: 0.15, 27: 0.05}, 0.05, (0.4, 1, 2, 3, 4, 5),
                             female_frac=0.25, diabetes=0.27, ckd=0.11, n_published=1118,
                             notes="A3 Table 1; label persistence / regression to the mean targets"),
    # PARTNER 2A / 3 style cohorts for the balloon-expandable point targets (B4): age ~82 (XT, intermediate risk) and
    # ~73 (Sapien 3, low risk); 30-day reference; annual echo; VARC-3 full definition
    "partner_xt": CohortProfile("partner_xt", "TAVR", {"Sapien XT": 1.0}, 81.8, 6.5, {23: 0.40, 26: 0.40, 29: 0.20}, 0.08, tuple(range(1, 6)),
                                female_frac=0.45, n_published=1665, notes="B4: SVD 1.6% at 5 y (VARC-3)"),
    "partner_s3": CohortProfile("partner_s3", "TAVR", {"Sapien 3": 1.0}, 73.3, 6.0, {23: 0.35, 26: 0.40, 29: 0.25}, 0.08, tuple(range(1, 6)),
                                female_frac=0.30, n_published=495, notes="B4: SVD 0.7% at 5 y (VARC-3)"),
}
