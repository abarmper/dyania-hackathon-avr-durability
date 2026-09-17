# Model Approach

Design rationale and every design choice: `docs/05_modelling_and_surveillance_plan.md`. Code: `notebooks/analysis/`. Numbers below are from the synthetic, literature-calibrated cohort (`notebooks/analysis/output/results_summary.md`).

---

## 1. Problem Formulation

**Dynamic competing-risks prediction at every echo (landmark design).** Each post-implant echo is a decision point. From everything known at that echo, the model predicts the cumulative incidence of confirmed moderate-or-greater SVD over the next 0.5–5 years, with death, endocarditis and non-SVD reintervention competing. An implant-time model gives the 5/8/10-y incidence from the reference echo.

Why this framing:
- **Right censoring and competing death.** Patients are elderly and follow-up is incomplete; about 40% of the cohort dies during follow-up. A fixed-horizon classifier would discard or mislabel them; a cause-specific hazard keeps every period of follow-up.
- **The label is known late.** A confirmed label exists only at the second abnormal echo, a median 0.94 y after the first. The outcome time is therefore the confirming echo, so a prediction never uses information from the future.
- **Surveillance needs interval risk.** "What is the risk before the next echo if we wait 6, 12, 24 or 36 months?" is read directly from the same model.

---

## 2. Chosen Model(s)

| Model | Role | Justification |
|---|---|---|
| Discrete-time cause-specific hazard model: gradient-boosted multinomial classifier (scikit-learn `HistGradientBoostingClassifier`) over 6-month periods after each echo | Primary | handles censoring by person-period expansion; softmax over {none, SVD, competing} gives calibrated cause-specific hazards and a monotone cumulative incidence; nonlinear effects of trajectories; native missing-value handling; SHAP-explainable; trains in seconds |
| Implant-time version of the same model (10-y horizon) | Primary for 5/8/10-y incidence | baseline risk at implant for protocol reporting |
| Cause-specific Cox + Aalen–Johansen (lifelines) | Interpretable model | hazard ratios with CIs for clinical review (`D5`) |
| Parsimonious 8-feature model | Deployment fallback | same accuracy with features any registry holds |
| Time since implant only | Baseline / comparator | what guideline surveillance implicitly assumes |
| Published risk-factor Cox (`B4` Table 2 predictors) | Baseline / comparator | the current nomogram-style evidence |
| Registry-only feature set | Baseline / comparator | isolates the value of engineered trajectory features |

Gradient boosting on person-period rows scales linearly to registry volumes (millions of rows, OpenMP-threaded) and is transparent at the level clinicians need (per-patient contributions). Joint longitudinal–survival models (`D1`) and deep survival models (`D3`) are the next step at registry scale; they are not needed to show the effect at 10,000 valves.

**Prototype performance (test set, 1,401 valves).**

| Model | Instances | AUC at 5 y (95% CI) | C at 5 y | Calibration slope |
|---|---|---|---|---|
| Primary, all echoes | 6,251 | 0.925 (0.911–0.939) | 0.919 | 0.96 |
| Primary, echoes not already abnormal | 5,555 | 0.875 (0.851–0.899) | 0.855 | 1.04 |
| Registry-only features, not already abnormal | 5,555 | 0.872 (0.848–0.896) | 0.854 | 1.05 |
| Time since implant only, not already abnormal | 5,555 | 0.651 (0.614–0.686) | 0.638 | 0.87 |
| Implant-time boosted model at 5 / 8 / 10 y | 1,489 | 0.865 / 0.850 / 0.823 | 0.844 / 0.817 / 0.783 | 0.94 / 0.95 / 0.96 |
| Cause-specific Cox at 5 / 8 / 10 y | 1,489 | 0.814 / 0.809 / 0.810 | 0.784 / 0.762 / 0.750 | 1.21 / 1.17 / 1.08 |
| Published risk-factor Cox at 5 / 8 / 10 y | 1,489 | 0.756 / 0.756 / 0.771 | 0.723 / 0.705 / 0.699 | 1.05 / 0.99 / 0.98 |

---

## 3. Feature Engineering

### Input Variables
Age, sex, height/weight (BSA, BMI), diabetes, CKD/dialysis, smoking, atrial fibrillation, anticoagulation type; approach, valve model and class, labelled size, implant year; reference echo (mean gradient, EOA, indexed EOA, DVI, intraprosthetic regurgitation, PPM class, and whether it came from an echo or the model × size table); every follow-up echo (mean and peak gradient, peak velocity, EOA, DVI, stroke-volume index, LVEF, regurgitation grades, and whether the echo was scheduled or symptom-triggered). Selection followed a variable audit: a variable enters only with a quantified effect or a defined role and a field on the TVT/STS forms (`docs/04` §7).

### Engineered Features

| Feature | Derivation | Clinical Rationale |
|---|---|---|
| Change in mean gradient from reference | last − reference | VARC-3 defines SVD as change from the post-implant baseline |
| Current single-echo stage (VARC-3 full and haemodynamic) and count of prior abnormal echoes | VARC-3 rules on each echo so far | the "unconfirmed abnormal echo" flag; a positive echo is the strongest short-term signal, but 59–65% do not persist (`A3`) |
| Smoothed log gradient and probability that the true gradient ≥20 mmHg | exponentially weighted mean of log gradient (1.5-y time constant); normal posterior with 14% per-echo multiplicative error | Doppler error is multiplicative; smoothing separates flicker from progression |
| Gradient slope since reference and over the last two echoes | least squares on all echoes; last-step difference | progression rate; explant risk rises with gradient progression (`B5`) |
| Maximum gradient so far, consecutive-worsening count, transient-spike flag | running maximum; run length of increases; positive then negative | distinguishes sustained deterioration from a transient rise (thrombosis, high flow) |
| Relative change in EOA and DVI | current / reference − 1 | flow-independent confirmation required by VARC-3 |
| Expected EOA for model × size and residual | ASE 2024 normal-value tables (`A5`) | a small orifice for the model and size means PPM or early dysfunction |
| Indexed EOA and PPM class | EOA / BSA; ASE 2024 Table 7 cut-offs (BMI-adjusted) | PPM raises SVD hazard (`B4`, `C1`) |
| Time since last echo, number of echoes, symptom-triggered flag | from echo dates and indication | informative missingness and clinical concern |

**Selection and explainability results.** SHAP top features: approach, valve model, prediction horizon, change in gradient from reference, time since implant, age, peak velocity, current haemodynamic stage. Backward elimination along permutation importance kept 8 features (valve model, time since implant, approach, current VARC-3 stage, change in gradient from reference, age, current haemodynamic stage, smoking) with AUC 0.924 vs 0.925 for the full model. The engineered trajectory features added only 0.003 AUC over the registry-only set in this cohort; the gain comes from the dynamic framing and the label handling, not from feature engineering.

### Handling Missing Data
- Channel-level missingness (EOA missing in ~8% and DVI in ~5% of echoes, `A3`) is handled natively by the boosted trees plus explicit missing indicators; VARC-3 staging falls back to the haemodynamic variant when both EOA and DVI are missing and is flagged "unevaluable".
- Missed visits are not imputed: time since last echo and number of echoes carry the information, and features use whatever echoes exist.
- A missing reference echo (31% of valves) is replaced by the model × size normal value with a flag; AUC in that group was 0.911 vs 0.933 with an echo reference, so the output carries a "reference from table" warning.
- Outcomes are never imputed; censoring is handled by the hazard model.

---

## 4. Validation Strategy

- **Train / Validation / Test split:** 70/15/15 by valve; the validation valves drive early stopping; the test valves are untouched until evaluation; 500-resample valve-level bootstrap CIs, paired across models.
- **Cross-validation approach:** no hyperparameter search in the prototype (fixed, conservative settings with early stopping) to avoid optimistic selection; nested grouped cross-validation is planned for the real cohort.
- **Temporal validation:** planned on implants 2021–2023 in the real study (device generation and the VKA-to-DOAC shift change over time). Not run in the prototype because simulated implant year carries no drift of its own.
- **External validation:** planned in a separate health system or registry. A face-validity check on the provided de-identified notes (`notebooks/real_data_check/`; 26 patients, 11 with documented structural dysfunction) found: the VARC-3 rule with a model × size reference had 86% sensitivity and 89% specificity against clinician documentation, with both false positives caused by the table reference; the landmark model scored at the follow-up echo reached AUC 0.97 (0.88–1.00), partly circular because dysfunction is diagnosed from the same echo; the implant-time model was uninformative (C 0.59, 0.36–0.82) because ages are redacted and the sample is tiny. In the prototype, robustness was tested by regenerating the cohort under changed mechanisms: echo noise ×0.5 and ×1.5 (5-y AUC 0.938 and 0.927 vs 0.934), progression ×0.7 (0.931), missed visits ×2 (0.969, fewer and later instances). Predicting the true, noise-free deterioration time with the same features is harder (AUC 0.88) than predicting the observed label (0.925), because the observed label shares measurement noise with the latest echo; AUCs on real echo labels therefore overstate how well true deterioration is anticipated.

---

## 5. Expected Model Outputs

A single inference runs when an echo report is finalised and returns:
- cumulative incidence of confirmed moderate-or-greater SVD at 1–5 y from this echo, and at 5/8/10 y from implant;
- risk tier on 5-y incidence: low <5%, moderate 5–15%, high ≥15%;
- top-3 contributing features with direction (SHAP);
- recommended next echo interval: the longest of 6, 12, 18, 24 or 36 months with predicted confirmed-SVD risk ≤1% and severe-SVD risk ≤2% over the interval;
- flags: unconfirmed abnormal echo (confirmation echo due in 3–6 months), reference from table, valve model with little history, unevaluable stage.

**Output format:**
"Echo 9.9 y after TAVR: 5-y cumulative incidence of confirmed moderate SVD 5.1% (moderate); top-3: TAVR (↓), change in gradient from reference +4.9 mmHg (↑), time since implant 9.9 y (↑); next echo interval from the 1% rule; no unconfirmed abnormal echo." (a real prototype record; `notebooks/analysis/output/tables/example_inference_records.json`)

---

## 6. Clinical Integration

The model sits in echo surveillance scheduling. It is triggered when a prosthetic-valve echo is signed in the reporting system, writes a result card to the EHR, and proposes the recall date to the valve clinic's scheduling list.
- **Low risk:** interval extended up to 36 months.
- **Moderate risk:** annual echo, as today.
- **High risk or rising trajectory:** echo in 6 months and structural heart team review of the images (CT if leaflet thrombosis is suspected).
- **Any single abnormal echo:** confirmation echo in 3–6 months regardless of risk.
- **Symptoms:** always trigger an echo, independent of the schedule.

In the prototype surveillance simulation, which includes feedback from detection to reintervention and identical patients under every policy, the risk-adapted schedule used 4.7 echoes per patient-decade versus 8.4 for guideline-plus-confirmation surveillance (44% fewer) with +0.7 months mean and +1.8 months P90 delay to confirmed detection, and detected 69% versus 73% of true deteriorations during follow-up. Guideline surveillance without a confirmation echo detected deterioration a mean 1.5 y (P90 2.9 y) after it began; adding the confirmation echo alone cut this to 1.1 y (P90 2.2 y).

---

## 7. Limitations and Failure Modes

- **Synthetic development data.** Results come from a literature-calibrated simulator (`docs/04`); real-world discrimination will be lower. They show that the design works and where it breaks, not the expected real AUC.
- **Valve classes with poor fit:** porcine/Mitroflow-class calibration slope 0.53 and Trifecta-class AUC 0.84, both with early, often regurgitant failure. The output flags these classes until recalibrated locally.
- **New valve models** borrow from their class; performance is unknown until local events accrue (flagged).
- **Missing reference echo** lowers AUC (0.911 vs 0.933) and biases change-from-baseline features (flagged).
- **The label itself.** Even confirmed VARC-3 detects only 60% of true moderate deteriorations during follow-up, with a median lag of 1.5 y and 6% false positives; the model inherits this.
- **Echo variability** beyond the simulated inter-visit noise (sonographer, vendor, flow state) and pressure recovery in small supra-annular valves can create false trajectories; the smoothed gradient and confirmation rule mitigate but do not remove this.
- **Simulator assumptions** that favour the model: symptoms occur only with true deterioration, and death is independent of SVD. The symptom flag is reported with and without it (AUC 0.874 without vs 0.875 with).
- Clinicians see these limits as flags on the result card and in a model-facts sheet listing validated populations.
