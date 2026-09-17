"""Figure registry for the Kaplan-Meier / cumulative-incidence reconstruction.

Every published number used for validation is recorded here with its source so the
reconstruction can be audited. Manual inputs (typed from figures) are marked MANUAL.

Curve `kind`:
  "KM"  -- survival-type curve (non-increasing), death censored (Kermen Fig 2)
  "CIF" -- cumulative incidence with death competing (Aalen-Johansen); stored and
           reconstructed as 1 - CIF ("AJ-as-KM", see plan §5.5)

Anchor values are in percent on the curve's own scale (survival % for KM, CIF % for CIF).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
REF_DIR = REPO / "references"
OUT_DIR = HERE / "output"

PDF_KERMEN = REF_DIR / "A4_velders2022_varc3_durability_bovine_pericardial.pdf"
PDF_NOTION = REF_DIR / "B1_thyregod2024_notion_10year_tavi_vs_savr.pdf"
PDF_WAKAMI = REF_DIR / "C1_ppm_trifecta_early_svd_2022.pdf"


@dataclass
class RiskTable:
    times: list[float]
    n: list[int]


@dataclass
class Anchor:
    t: float                # years; use None for "end of follow-up" anchors
    value: float            # percent
    se: float | None = None # published standard error in percent, if given
    note: str = ""


@dataclass
class CurveSpec:
    id: str
    pdf: Path
    page: int               # 0-based page index
    figure: str
    panel: str
    name: str
    arm: str
    kind: str               # "KM" | "CIF"
    source: str             # "vector" | "raster_rgb" | "raster_gray"
    n0: int
    stratum: str
    estimator: str          # "KM_death_censored" | "AJ_death_competing"
    anchors: list[Anchor] = field(default_factory=list)
    risk: RiskTable | None = None
    tot_events: int | None = None
    legend_label: str | None = None     # vector: legend text to match
    colour_hint: str | None = None      # raster: "red" | "blue" | "black" | "gray"
    xobject: str | None = None          # raster: image XObject name
    panel_index: int | None = None      # raster: 0 = left panel, 1 = right panel
    y_tick_step: float | None = None    # raster: percent per y tick (MANUAL)
    x_max_years: float | None = None    # raster: last x tick value
    use_risk_table: bool = True         # False: printed row is not this curve's risk set -> Guyot 'no numbers at risk' case
    use_tot_events: bool = True         # False: published total inconsistent with the curve -> curve shape wins, total reported only
    censoring_from: str | None = None   # id of a curve from the SAME panel/population whose reconstructed censoring times are reused
    notes: str = ""


# ----------------------------------------------------------------------------------
# Kermen et al. 2022, JTCVS Open 11:72-80 (Magna Ease, n=338, mean FU 6.6 y, 2238 valve-y)
# Numbers below are quoted from the PDF text (Table 2, Figure 2, Figure 4, Results).
# Risk table (Fig 2A/2B): "Patients at risk, n 338 ... 247 ... 34 20" at 0, 5, 9.5, 10 y.
# Fig 4 risk table: "Patients at risk, n 248 338 20 34" -> 338, 248, 34, 20.
# Events: mortality 4 early + 84 late = 88; valve-related mortality 4 + 21 = 25;
# SVD stage 2/3 = 49; stage 3 = 11; explantation due to SVD = 9; non-valve mortality 63.
# ----------------------------------------------------------------------------------
KERMEN_RISK = RiskTable(times=[0, 5, 9.5, 10], n=[338, 247, 34, 20])
KERMEN_RISK_FIG4 = RiskTable(times=[0, 5, 9.5, 10], n=[338, 248, 34, 20])

KERMEN = [
    CurveSpec(
        id="kermen_fig2a_overall_survival", pdf=PDF_KERMEN, page=3, figure="Fig 2", panel="A",
        name="Overall survival", arm="all", kind="KM", source="vector", n0=338,
        stratum="SAVR_perimount_class_age70", estimator="KM_death_censored",
        anchors=[Anchor(5, 80.9, 2.2), Anchor(10, 66.7, 4.4)],
        risk=KERMEN_RISK, tot_events=88, legend_label="Overall survival",
        notes="Table 2: mortality 4 early + 84 late (3.75%/vy); KM 80.9 +/- 2.2 at 5 y, 66.7 +/- 4.4 at 10 y"),
    CurveSpec(
        id="kermen_fig2a_valve_related_survival", pdf=PDF_KERMEN, page=3, figure="Fig 2", panel="A",
        name="Valve-related survival", arm="all", kind="KM", source="vector", n0=338,
        stratum="SAVR_perimount_class_age70", estimator="KM_death_censored",
        anchors=[Anchor(5, 92.5, 1.5), Anchor(10, 86.0, 6.1)],
        risk=KERMEN_RISK, tot_events=25, legend_label="Valve-related survival",
        notes="Table 2: valve-related mortality 4 early + 21 late; KM 92.5 +/- 1.5, 86.0 +/- 6.1"),
    CurveSpec(
        id="kermen_fig2b_freedom_stage23", pdf=PDF_KERMEN, page=3, figure="Fig 2", panel="B",
        name="Freedom from moderate/severe SVD (VARC-3 stage 2/3)", arm="all", kind="KM",
        source="vector", n0=338, stratum="SAVR_perimount_class_age70", estimator="KM_death_censored",
        anchors=[Anchor(5, 98.5, 0.7), Anchor(10, 60.9, 7.0), Anchor(None, 26.7, None, "end of follow-up")],
        risk=KERMEN_RISK, tot_events=49, legend_label="Moderate/Severe SVD",
        notes="Table 2: SVD stage 2-3 49 events (2.19%/vy); KM 98.5 +/- 0.7 at 5 y, 60.9 +/- 7.0 at 10 y; "
              "text says 98.6 and 61.1. THIS IS THE PRIMARY-LABEL CURVE."),
    CurveSpec(
        id="kermen_fig2b_freedom_stage3", pdf=PDF_KERMEN, page=3, figure="Fig 2", panel="B",
        name="Freedom from severe SVD (VARC-3 stage 3)", arm="all", kind="KM",
        source="vector", n0=338, stratum="SAVR_perimount_class_age70", estimator="KM_death_censored",
        anchors=[Anchor(5, 99.6, 0.4), Anchor(10, 88.3, 5.0), Anchor(None, 58.9, None, "end of follow-up")],
        risk=KERMEN_RISK, tot_events=11, legend_label="Severe SVD", use_risk_table=False, use_tot_events=False,
        censoring_from="kermen_fig2b_freedom_stage23",
        notes="Table 2: SVD stage 3 11 events (0.49%/vy); KM 99.6 +/- 0.4, 88.3 +/- 5.0. "
              "The printed row (338/247/34/20) is the all-cause survival risk set: with 247 at risk the first "
              "drop (0.9 pp at 5.0 y) would be 2-3 events and the curve implies >= 14 events, not 11; the exact "
              "step inversion implies ~115 echo-followed patients at risk at 5 y. Forcing 11 events flattens the tail "
              "after 10 y (the curve falls 87 -> 58 there), so the curve shape wins: reconstructed with n0 and the "
              "censoring times of the stage 2/3 curve from the same panel (everyone at risk for stage 2/3 is at risk "
              "for stage 3); the published total is reported as informational."),
    CurveSpec(
        id="kermen_fig4a_explant_svd", pdf=PDF_KERMEN, page=6, figure="Fig 4", panel="A",
        name="Explantation due to SVD", arm="all", kind="CIF", source="vector", n0=338,
        stratum="SAVR_perimount_class_age70", estimator="AJ_death_competing",
        anchors=[Anchor(5, 0.0, None), Anchor(10, 8.0, 3.4)],
        risk=KERMEN_RISK_FIG4, tot_events=9, legend_label="Explantation due to SVD",
        notes="Fig 4A: cumulative incidence 0 at 5 y, 8.0 +/- 3.4 at 10 y; 9 explants for SVD"),
    CurveSpec(
        id="kermen_fig4a_valve_related_death", pdf=PDF_KERMEN, page=6, figure="Fig 4", panel="A",
        name="Valve-related death", arm="all", kind="CIF", source="vector", n0=338,
        stratum="SAVR_perimount_class_age70", estimator="AJ_death_competing",
        anchors=[Anchor(5, 7.1, 1.4), Anchor(10, 7.4, 1.5)],
        risk=KERMEN_RISK_FIG4, tot_events=25, legend_label="Valve-related death",
        notes="Fig 4A: 7.1 +/- 1.4 at 5 y, 7.4 +/- 1.5 at 10 y"),
    CurveSpec(
        id="kermen_fig4a_nonvalve_death", pdf=PDF_KERMEN, page=6, figure="Fig 4", panel="A",
        name="Non valve-related death", arm="all", kind="CIF", source="vector", n0=338,
        stratum="SAVR_perimount_class_age70", estimator="AJ_death_competing",
        anchors=[Anchor(5, 12.0, 1.8), Anchor(10, 25.4, 4.1)],
        risk=KERMEN_RISK_FIG4, tot_events=63, legend_label="Non valve-related death",
        notes="Fig 4A: 12.0 +/- 1.8 at 5 y, 25.4 +/- 4.1 at 10 y; Table 2 not valve-related mortality 63 late"),
]

# ----------------------------------------------------------------------------------
# NOTION 10-year, Thyregod et al. EHJ 2024;45:1116-1124. As-implanted n = 139 TAVI / 135 SAVR;
# echo cohort at risk at time 0 (Fig 3 risk rows) = 134 TAVI / 123 SAVR, used as n0.
# Fig 3 (page 6, /Im0 1266x656): >=moderate SVD (left) and severe SVD (right), VARC-3
# haemodynamic definition, cumulative incidence with death (and reintervention) competing.
# MANUAL: numbers at risk are inside the JPEG; type them from output/qa_notion_fig3_risk_table.png
# ----------------------------------------------------------------------------------
# Typed by hand from output/qa_notion_fig3_risk_table.png (Fig 3 "Patients at risk" rows, yearly 0..10).
_Y = list(range(11))
NOTION_RISK_MODSVD = {"TAVI": RiskTable(_Y, [134, 131, 128, 117, 109, 96, 82, 71, 56, 44, 30]),
                      "SAVR": RiskTable(_Y, [123, 122, 116, 107, 96, 84, 69, 61, 48, 41, 32])}
NOTION_RISK_SEVSVD = {"TAVI": RiskTable(_Y, [134, 132, 129, 118, 109, 96, 82, 73, 62, 51, 40]),
                      "SAVR": RiskTable(_Y, [123, 122, 119, 110, 100, 91, 79, 70, 58, 50, 39])}

NOTION = [
    CurveSpec(
        id="notion_fig3_modsvd_tavi", pdf=PDF_NOTION, page=5, figure="Fig 3", panel="left",
        name=">=moderate SVD", arm="TAVI", kind="CIF", source="raster_rgb", n0=134,
        stratum="TAVR_self_expanding_age79", estimator="AJ_death_competing",
        anchors=[Anchor(10, 15.4, None)], risk=NOTION_RISK_MODSVD["TAVI"], tot_events=None,
        colour_hint="red", xobject="/Im0", panel_index=0, y_tick_step=10.0, x_max_years=10,
        notes="Results: moderate or severe SVD TAVI 15.4% vs SAVR 20.8% at 10 y (HR 0.7)"),
    CurveSpec(
        id="notion_fig3_modsvd_savr", pdf=PDF_NOTION, page=5, figure="Fig 3", panel="left",
        name=">=moderate SVD", arm="SAVR", kind="CIF", source="raster_rgb", n0=123,
        stratum="SAVR_mixed_age79", estimator="AJ_death_competing",
        anchors=[Anchor(10, 20.8, None)], risk=NOTION_RISK_MODSVD["SAVR"], tot_events=None,
        colour_hint="blue", xobject="/Im0", panel_index=0, y_tick_step=10.0, x_max_years=10),
    CurveSpec(
        id="notion_fig3_sevsvd_tavi", pdf=PDF_NOTION, page=5, figure="Fig 3", panel="right",
        name="severe SVD", arm="TAVI", kind="CIF", source="raster_rgb", n0=134,
        stratum="TAVR_self_expanding_age79", estimator="AJ_death_competing",
        anchors=[Anchor(10, 1.5, None)], risk=NOTION_RISK_SEVSVD["TAVI"], tot_events=None,
        colour_hint="red", xobject="/Im0", panel_index=1, y_tick_step=10.0, x_max_years=10,
        notes="severe SVD 1.5% vs 10.0% (HR 0.2)"),
    CurveSpec(
        id="notion_fig3_sevsvd_savr", pdf=PDF_NOTION, page=5, figure="Fig 3", panel="right",
        name="severe SVD", arm="SAVR", kind="CIF", source="raster_rgb", n0=123,
        stratum="SAVR_mixed_age79", estimator="AJ_death_competing",
        anchors=[Anchor(10, 10.0, None)], risk=NOTION_RISK_SEVSVD["SAVR"], tot_events=None,
        colour_hint="blue", xobject="/Im0", panel_index=1, y_tick_step=10.0, x_max_years=10),
]

# ----------------------------------------------------------------------------------
# Wakami et al. 2022, J Cardiothorac Surg 17:174 (Trifecta, n=110, 7 SVD events, mean FU 66 mo).
# Fig 1 (page 4, /Im0 1772x844, grayscale): cumulative incidence of SVD, death competing.
# Published: SVD 4.8% at 5 y and 6.6% at 7 y. Fig 1 has a 'Number at risk' row (typed from
# output/qa_wakami_fig1_source.png): 110, 100, 84, 58, 24 at 0, 20, 40, 60, 80 months. x axis in MONTHS.
# ----------------------------------------------------------------------------------
WAKAMI = [
    CurveSpec(
        id="wakami_fig1_svd", pdf=PDF_WAKAMI, page=3, figure="Fig 1", panel="single",
        name="SVD (Trifecta)", arm="all", kind="CIF", source="raster_gray", n0=110,
        stratum="SAVR_trifecta_class_age78", estimator="AJ_death_competing",
        anchors=[Anchor(5, 4.8, None), Anchor(7, 6.6, None)],
        risk=RiskTable(times=[0, 20 / 12, 40 / 12, 60 / 12, 80 / 12], n=[110, 100, 84, 58, 24]), tot_events=7,
        colour_hint="black", xobject="/Im0", panel_index=0, y_tick_step=20.0, x_max_years=80 / 12,
        notes="Results: SVD rate 4.8% at 5 y, 6.6% at 7 y; 7 events; all-cause death 14.7%/23.8% at 5/7 y"),
]

ALL_CURVES = KERMEN + NOTION + WAKAMI
CURVES_BY_ID = {c.id: c for c in ALL_CURVES}

# Stratum mapping used by the simulator's onset module (docs/04 §2.1)
STRATA = {
    "SAVR_perimount_class_age70": "Perimount / Magna Ease class, SAVR, mean age 70 (Kermen)",
    "SAVR_mixed_age79": "Mixed stented pericardial/porcine, SAVR, mean age 79 (NOTION SAVR arm)",
    "TAVR_self_expanding_age79": "CoreValve self-expanding, TAVI, mean age 79 (NOTION TAVI arm)",
    "SAVR_trifecta_class_age78": "Trifecta class, SAVR, mean age 78 (Wakami)",
}
