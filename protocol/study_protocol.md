# Study Protocol

Supporting detail: `docs/03_ideas_brainstorm.md` §0 (endpoints), `docs/04_synthetic_cohort_spec.md` (prototype cohort), `docs/05_modelling_and_surveillance_plan.md` (analysis design), `notebooks/analysis/output/results_summary.md` (prototype results). Reporting follows TRIPOD+AI (`references/D8`).

---

## 1. Study Title and Objectives

**Title:**
Dynamic prediction of structural valve deterioration after bioprosthetic SAVR and TAVR to personalise echocardiographic surveillance: a retrospective multicentre development and external validation study.

**Primary Objective:**
To determine whether routinely collected post-implant echocardiographic and EHR data can predict an individual's time to confirmed moderate-or-greater haemodynamic structural valve deterioration (SVD) accurately enough to personalise the echo surveillance interval.

**Secondary Objectives:**
1. Differentiate SAVR vs TAVR, and valve-model-specific, deterioration trajectories.
2. Compare a risk-adapted echo schedule with the guideline schedule (ESC/EACTS 2025: within 3 months of implantation, at 1 year, then annually) on echoes used and detection delay.
3. Predict the reintervention window: time from first confirmed moderate SVD to severe SVD or reintervention, and the share of reinterventions that are urgent.
4. Quantify the reliability of the SVD label: agreement and persistence of the VARC-3, Capodanno 2017 and Dvir 2018 definitions.

---

## 2. Target Population

### Inclusion Criteria
- Adults (≥18 y) with a first bioprosthetic aortic valve implanted surgically (stented or stentless) or transcatheter, implanted 2010–2020 (development) or 2021–2023 (temporal validation).
- Valve model and labelled size recorded (STS `AorticImplant`/UDI or TVT device fields).
- A reference echo within 3 months of implant (the 30-day/discharge echo is accepted) **or** a model × size reference value, with the source flagged.
- At least one follow-up transthoracic echo ≥6 months after implant.

### Exclusion Criteria
- Mechanical, homograft or autograft (Ross) valves; valve-in-valve as the index procedure.
- Concomitant mitral or tricuspid valve replacement (multi-valve prosthesis confounds gradients).
- Active endocarditis as the indication for the index operation.
- Death or reintervention within 30 days, or no follow-up echo.
- Endocarditis, thrombosis and paravalvular leak **during follow-up are not exclusions**: they are competing events or covariates (§4 ground truth).

### Cohort Size Estimate
≥5,000 valves from ≥3 centres for development and ≥2,000 valves from a separate health system or registry for external validation, with median echo follow-up ≥5 y. Rationale: confirmed moderate SVD accrues slowly (freedom 98.5% at 5 y and 60.9% at 10 y in a Perimount-class cohort, `A4`), so 5,000 valves yield roughly 150–300 events by 5 y and >500 by 10 y. In the prototype, the 5-y AUC 95% CI width fell from 0.095 (500 valves) to 0.050 (1,000), 0.034 (2,000) and 0.022 (5,000); see §5.

---

## 3. Endpoints

### Primary Endpoint
Time from the reference echo to **VARC-3 moderate-or-greater (stage 2 or 3) haemodynamic SVD confirmed on two consecutive echoes**, with death, endocarditis and non-SVD reintervention as competing risks. Stage 2: mean gradient ≥20 mmHg and increase ≥10 mmHg from reference, with EOA decrease ≥0.3 cm² or ≥25% or DVI decrease ≥0.1 or ≥20%, or new ≥moderate intraprosthetic regurgitation. Reintervention for SVD before the confirming echo counts as an event. Reported as cumulative incidence at 5, 8 and 10 y. The event is dated at the first echo of the confirmed pair for incidence reporting and at the confirming echo for dynamic prediction (the time the label becomes known).

**Clinical usefulness threshold.** Statistical: time-dependent AUC ≥0.75 and C-index ≥0.70 at 5 y, calibration slope 0.8–1.2 at 5/8/10 y. Decision level: a risk-adapted schedule uses ≥25% fewer echoes per patient-decade with mean and 90th-percentile detection delay no more than 3 months longer than guideline surveillance.

### Secondary Endpoints

| Endpoint | Measurement | Timeframe |
|---|---|---|
| Bioprosthetic valve failure (VARC-3) | confirmed severe SVD, reintervention after a dysfunction diagnosis, or valve-related death; Aalen–Johansen cumulative incidence | 5 / 8 / 10 y |
| Time to aortic valve reintervention | ViV-TAVR or redo SAVR from procedure records; elective vs urgent | 10 y |
| Mean gradient progression | annualised change from reference, linear mixed model with random slopes | whole follow-up |
| Surveillance performance | echoes per patient-decade; delay from the first qualifying echo to confirmed detection; share of severe SVD first seen while moderate | 10 y |
| Non-structural dysfunction | ≥moderate paravalvular leak, clinical valve thrombosis, endocarditis; reported descriptively as competing events | 10 y |

---

## 4. Proposed Data Sources

| Source | Data Type | Access Pathway | Key Variables |
|---|---|---|---|
| Echocardiography reporting database (structured measurement tables) | serial TTE measurements | institutional data warehouse under IRB | mean and peak gradient, peak velocity, EOA, DVI, LVEF, stroke-volume index, intraprosthetic and paravalvular regurgitation grade, exam date |
| STS Adult Cardiac Surgery Database v4.20 | SAVR procedure, device, comorbidities | participant data file / STS research access | `AorticImplant`, UDI, size, `Diabetes`, `Dialysis`, `TobaccoUse`, `CalculatedBMI`, `PreLVEF`, discharge anticoagulant |
| STS/ACC TVT Registry v3.0 | TAVR procedure, device, 30-day and 1-y echo | TVT research proposal (analysis at registry centre) | device model and size, AV mean gradient (13674–13676), AV area, PVL, valve thrombosis, reintervention |
| EHR (labs, medications, notes) | comorbidities, anticoagulation, symptoms, missing implant details | FHIR/warehouse extract; NLP on operative and echo notes (`E1`) | AF, anticoagulant type, creatinine, valve model when uncoded, reintervention indication |
| CMS claims / National Death Index linkage | reintervention and death after loss to follow-up | CMS DUA; NDI application | ViV/redo procedure codes, date and cause of death |

**Ground Truth Definition:**
Hierarchy, highest first: (1) explant or autopsy pathology showing leaflet calcification, tear or pannus; (2) reintervention record with a documented SVD indication; (3) VARC-3 echo criteria (above) met on two consecutive echoes against the reference echo, with the flow-independent criterion (EOA or DVI) required when available; (4) structural heart team adjudication of discordant or unevaluable cases, blinded to model output. Endocarditis, clinical valve thrombosis (gradient normalising with anticoagulation), paravalvular leak and prosthesis–patient mismatch are adjudicated out of the SVD label at every level. Capodanno 2017, Dvir 2018 and single-echo VARC-3 labels are sensitivity analyses.

---

## 5. Statistical Analysis Plan

### Sample Size
Two criteria, the larger governs. (a) Events per candidate parameter ≥10 for the interpretable cause-specific Cox model with ~20 predictors → ≥200 primary events; Riley et al. minimum-sample-size criteria for time-to-event models are checked with `pmsampsize`. (b) Precision: simulation on the calibrated synthetic cohort gives a 5-y AUC 95% CI width of 0.034 with 2,000 valves and 0.022 with 5,000. Target ≥5,000 development valves (≥300 events) and ≥2,000 external-validation valves.

### Train / Validation / Test Split
Development cohort split by valve (70% train, 15% early-stopping validation, 15% internal test); every instance of a valve stays in one partition, and bootstrap resampling is by valve. Temporal validation on implants 2021–2023; external validation in a separate health system or registry. Features use only echoes at or before the prediction time; the event time for dynamic models is the confirming echo, so no future echo leaks into a prediction. Low event rates are handled by modelling hazards per 6-month period (probabilities remain calibrated) instead of resampling or class weights.

### Evaluation Metrics
Cause-specific time-dependent AUC (competing events counted as controls) at 3/5 y from each echo and 5/8/10 y from implant; truncated cause-specific C-index; Brier score; calibration intercept and slope with decile plots against Aalen–Johansen incidence; all with valve-level bootstrap 95% CIs. Metrics are also reported on echoes that are not already abnormal, because that is where a scheduling decision is made. Decision-level metrics are the surveillance endpoints above. Missed deterioration is weighted more than an extra echo: a missed severe SVD risks urgent redo surgery (22.6% vs 1.4% elective mortality, `B3`), while an extra echo is cheap and harmless.

### Subgroup Analyses
Implant approach (SAVR/TAVR), valve class and model (including new models with little history), age at implant (<60, 60–80, >80), labelled size and indexed EOA, sex, reference-echo source (30-day echo vs model × size table). These are the axes along which deterioration rates, echo noise and data completeness differ, so a pooled metric can hide failure.

### Comparator / Baseline
Time since implant alone; a published risk-factor Cox model with the predictor set of the 2025 JAHA review (`B4` Table 2: age, sex, BMI, diabetes, smoking, CKD, PPM, baseline gradient, valve size, anticoagulation); a registry-only feature set; for surveillance, the guideline schedule, guideline plus a confirmation echo, and rule-based de-escalation. STS-PROM enters only as a covariate for competing mortality. Cardiologist gestalt is assessed in a prospective silent-mode phase.

**Prototype evidence (synthetic, literature-calibrated cohort of 10,000 valves).** Dynamic model on echoes not already abnormal: 5-y AUC 0.875 (0.851–0.899), C 0.855, calibration slope 1.04; time-since-implant comparator 0.651; implant-time model 0.865/0.850/0.823 at 5/8/10 y vs 0.756/0.756/0.771 for the published risk-factor Cox. A risk-adapted schedule (1% threshold) used 44% fewer echoes than guideline-plus-confirmation surveillance with +0.7 months mean and +1.8 months P90 detection delay.

---

## 6. Ethical Considerations and Data Privacy

### IRB / Ethics Review
Retrospective development and validation: IRB review at each site with waiver of informed consent (minimal risk, impracticable, 45 CFR 46.116(f)); registry analyses under the STS/ACC research approval process. A prospective silent-mode phase (model runs, outputs hidden) needs IRB approval without consent; any interventional surveillance trial needs full review, consent and trial registration.

### Data Privacy
Analysis on a HIPAA limited data set under data use agreements, because time-to-event analysis needs exact intervals (Safe Harbor year-only dates break them). Direct identifiers never leave the covered entity; linkage uses hashed tokens; models train inside the health-system environment or in a federated set-up for multi-site work. Model artefacts contain no patient-level data; inference runs in the EHR with access logging. The public prototype uses synthetic data only.

### Algorithmic Fairness
Discrimination and calibration are reported within sex, age band, race and ethnicity (real cohort only), approach, valve manufacturer and model, and size, with pre-specified floors (AUC ≥0.70, slope 0.8–1.2). Women and small-annulus patients are a known risk group through prosthesis–patient mismatch, so sex is analysed as an effect modifier, not only a covariate. Race is not a predictor. Subgroups below the floor get subgroup recalibration or are flagged as out of scope in the output. Prototype: 5-y AUC 0.932 in women and 0.919 in men, calibration slopes 0.99 and 0.94.

### Clinical Transparency
Each inference shows the 5/8/10-y cumulative incidence, a risk tier, the top-3 contributing features, the recommended next echo interval (bounded 6–36 months), and flags for an unconfirmed abnormal echo, a missing reference echo or a valve model with little history. The system only schedules surveillance. It never recommends reintervention; a flagged valve goes to structural heart team review with the images. Any single abnormal echo triggers a confirmation echo at 3–6 months regardless of predicted risk, symptoms always override the schedule, and clinician overrides are logged for monitoring.
