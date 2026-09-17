# Literature, datasets and where the novelty is

Companion to `01_problem_and_data_understanding.md`. Everything cited here was verified against Europe PMC or a live publisher/dataset page. PDFs are in `references/` with an index in `references/README.md`.

---

## 1. What the evidence actually says

### 1.1 The label is the hard part, and this is our biggest finding

The field has **three competing definitions** of structural valve deterioration (SVD): EAPCI/ESC/EACTS 2017 (Capodanno), Dvir 2018, and VARC-3 2021. They are not interchangeable.

Velders and colleagues (*JTCVS Open* 2024, `A3`) took 1,118 patients who all received **the same stented bioprosthesis**, with echocardiograms read by an **independent core laboratory**, and applied all three definitions over five years:

| Definition | Patients ever labelled with SVD |
|---|---|
| Capodanno | 51 (4.6%) |
| Dvir | 32 (2.9%) |
| VARC-3 | 34 (3.0%) |

95% of patients were never labelled by any definition. Critically, **after a patient was first classified as having SVD, 59–65% had no SVD classification at the next visit.** The label flickers on and off.

Their conclusion, verbatim: the definitions are *"strong negative predictors but inconsistent positive discriminators."*

This is the single most consequential fact for our study design. It means:

- A naive label built by thresholding one echo is mostly measurement noise. Any model trained on it will learn the noise.
- The correct target is the **underlying trajectory** of valve function, not a single threshold crossing.
- We should require **persistence** (confirmation on a subsequent echo) before calling an event, and we must report how results change across all three definitions.

Most hackathon teams will not read this paper, will define SVD as "mean gradient ≥20 mmHg", and will report an impressive area under the curve on a noisy label. Building the label properly is both scientifically right and a visible differentiator.

The VARC-3 thresholds themselves (verified verbatim from the NOTION paper, `B1`):

- **Moderate SVD**: mean gradient ≥20 mmHg **and** increase ≥10 mmHg from the 3-month echo, **or** new ≥moderate intraprosthetic regurgitation.
- **Severe SVD**: mean gradient ≥30 mmHg **and** increase ≥20 mmHg from the 3-month echo, **or** new severe intraprosthetic regurgitation.
- **Valve failure**: valve-related death, or severe haemodynamic SVD, or reintervention after a dysfunction diagnosis.

Note the reference point is the **3-month** echo, not discharge. Gradients fall in the first weeks; using the discharge echo as baseline inflates apparent deterioration.

### 1.2 Event rates: rare, late, and competing with death

From NOTION at 10 years (`B1`), the longest randomised follow-up available:

| Endpoint at 10 years | TAVI | SAVR |
|---|---|---|
| ≥moderate SVD | 15.4% | 20.8% |
| Severe SVD | 1.5% | 10.0% |
| Valve failure | ~10% | 13.8% |
| Reintervention | 4.3% | 2.2% |
| Severe non-structural dysfunction (mostly mismatch) | 20.5% | 43.0% |
| Endocarditis | 7.2% | 7.4% |

PARTNER 2 at 5 years: SVD 2.20% for TAVR vs 4.38% for SAVR.

Three design consequences. First, **events are rare and late**, so the cohort must be large and followed long, and a fixed 5-year horizon will capture very little. Second, **death competes with SVD** in a population averaging 79 years, so a Fine-Gray or cause-specific competing-risk model is mandatory, not optional (`D5`). Third, **non-structural dysfunction is far more common than structural**, so a model that cannot separate mismatch from true degeneration will mostly predict mismatch.

### 1.3 Risk factors are known but have never been assembled into a model

Consistent predictors across the literature: younger age at implant, patient–prosthesis mismatch and small indexed orifice area (`C1`), specific valve models with known early failure such as Mitroflow and Trifecta, renal failure, diabetes and metabolic syndrome, dyslipidaemia, calcium-phosphate disturbance, elevated early post-implant gradient, and absence of early anticoagulation. Mechanisms are reviewed in `B2` (calcification, lipid infiltration, immune response to alpha-Gal, mechanical fatigue, thrombosis).

Salaun and colleagues' 2018 *Circulation* study of 1,387 patients found degeneration in 31%, **most of it subclinical**, with a dysmetabolic profile and leaflet calcification on computed tomography as independent predictors. That last point matters: deterioration is usually silent before it is symptomatic, which is exactly the window a surveillance model is trying to exploit.

### 1.4 What current practice does

The current 2025 ESC/EACTS guideline (Praz et al., §14.2) is explicit, verbatim:

> Serial TTE measurements of transprosthetic gradients, calculation of the effective valve area, and evaluation of leaflet motion and morphology should be performed in patients receiving a BHV within 3 months after valve implantation, again at 1 year, and annually thereafter, or sooner if new cardiovascular symptoms occur.

The 2021 edition (`A2`) set the baseline "within 30 days after valve implantation"; the rest of the schedule is unchanged.

So the standard of care is **a fixed annual interval, identical for every patient**, regardless of whether they have a 19 mm Mitroflow in a 55-year-old diabetic on dialysis or a 29 mm Evolut in an 85-year-old. That uniformity is the inefficiency our proposal targets, and having the guideline sentence quoted verbatim on a slide makes the argument concrete.

---

## 2. Datasets

### 2.1 The one that changes everything: MIMIC-IV-ECHO v1.0.1

**Verified live on the PhysioNet page.** Version 1.0.1, published **25 August 2026** — three weeks before this hackathon.

It now ships `structured_measurement.csv`: structured echocardiographic measurements for **206,488 studies from 91,372 unique patients, spanning 2008–2022**. Measurement categories include, verbatim, "chamber sizes and volumes, systolic and diastolic function, **valvular morphology and function**, aortic root and great vessel dimensions, and **Doppler-derived hemodynamics**."

This is the first time serial, structured, Doppler echo measurements have been publicly available and linkable to a full electronic health record. Joining on `subject_id` to MIMIC-IV's `procedures_icd` gives a bioprosthetic AVR cohort with repeated post-implant echoes and death outcomes.

Why this matters for novelty: **the resource is three weeks old.** Nobody has published an SVD durability analysis on it. Proposing it as the validation cohort is both credible and genuinely first-mover.

Caveats, stated honestly:
- The measurement dictionary sits behind the credentialed login, so the exact variable names for mean gradient, effective orifice area and dimensionless index are unconfirmed. The category description makes gradients very likely; **valve manufacturer and model are almost certainly absent** from structured fields.
- Echo **reports** are still not released; the page says they will come to the Note module later.
- Access requires PhysioNet credentialing: a CITI "Data or Specimens Only Research" course, a reference, and a signed data use agreement. Published turnaround is several days to two weeks.
- Single centre, intensive-care-enriched, and follow-up ends when the patient stops attending that hospital.

Cohort identification uses ICD-9 35.21 (tissue AVR) and 35.05 (transcatheter, from 2011), plus ICD-10-PCS 02RF08Z / 02RF0KZ (open, tissue) and 02RF3\*Z (transcatheter). The centre switched coding systems in late 2015, so both are needed. Published counts in MIMIC-IV v2.0 suggest roughly 400–600 bioprosthetic surgical cases and 800–1,500 transcatheter cases; v3.1 adds 2020–2022 and should be larger.

### 2.2 Everything else, ranked by usefulness

| Resource | Contains serial echo? | Valve model? | Access | Verdict |
|---|---|---|---|---|
| **MIMIC-IV + MIMIC-IV-ECHO** | Yes, structured, 2008–2022 | No | Credentialed, days to 2 weeks | Best public option by a wide margin |
| **MIMIC-IV-Note** | No echo reports | Possibly, in discharge summaries | Credentialed | Only public route to valve model/size, via text extraction |
| **MIMIC-III + ECHO-NOTE2NUM** | Free-text echo reports, 43,472 reports / 21,572 patients, 2001–2012 | No | Credentialed | Older, but longer death follow-up; ECHO-NOTE2NUM gives severity fields plus raw text |
| **STS/ACC TVT Registry** | Yes, with mandated follow-up | **Yes** | Research proposal; analysis done by their analytic centre, data never leaves | Ideal for the real study, impossible for a hackathon |
| **Trial data (PARTNER, CoreValve, NOTION) via Vivli / YODA / BioLINCC** | Yes, protocol-driven | Yes | Formal application, months | Real study only |
| **UK Biobank** | Cardiac MRI subset, not serial post-AVR echo | No | Application plus fee | Weak for this question |
| **EchoNext** (PhysioNet) | ECG with echo-derived labels | No | Registered plus agreement, same day | **Explicitly excluded prosthetic valves.** Not usable |
| **EchoNet-Dynamic / LVH / CAMUS / TMED** | Images only, no serial, no prostheses | No | Agreement, days | Pipeline demos only |
| **FDA MAUDE / post-approval studies** | No | Yes, by model | Open | Useful for model-specific failure signals and denominators |

**No public dataset anywhere contains prosthetic-valve or transcatheter echo images with durability labels.** Every echo image dataset we checked either excludes prostheses or does not annotate them. That is a gap worth stating out loud in the protocol.

### 2.3 The prototype route that needs no data access

Given a three-day window and credentialing that takes one to two weeks, the honest prototype path is a **literature-calibrated synthetic cohort**:

1. Take published Kaplan-Meier durability curves stratified by valve model, age band, and mismatch status (NOTION, PARTNER, Perimount, Trifecta, Mitroflow series).
2. Reconstruct patient-level time-to-event data from those curves using the Guyot algorithm (`D9`), which is an accepted, citable method rather than invented data.
3. Simulate serial gradient trajectories on top, with measurement error calibrated to the flicker rate observed by Velders (`A3`) — we know roughly 60% of first SVD classifications do not persist, which gives a defensible noise level.
4. Train and validate the pipeline on this, then present MIMIC-IV-ECHO as the pre-registered external validation cohort.

This is defensible precisely because it is labelled as synthetic and calibrated to real published curves, and because the reconstruction step has a methods citation behind it.

---

## 3. The novelty gap

I searched Europe PMC and the open web for existing models. Findings, stated plainly.

**Searches run against Europe PMC, with result counts:**

| Query | Hits |
|---|---|
| `TITLE:"machine learning" AND (TITLE:"structural valve deterioration" OR TITLE:"bioprosthetic valve")` | **0** |
| `TITLE:"risk score" AND (TITLE:"bioprosthetic" OR TITLE:"valve deterioration")` | **0** |
| `(TITLE:"surveillance" OR TITLE:"follow-up interval") AND TITLE:"prosthetic valve"` | **0** |
| `(TITLE:"artificial intelligence" OR TITLE:"machine learning") AND TITLE:"valve durability"` | **0** |
| `(ABSTRACT:"machine learning" OR ABSTRACT:"deep learning") AND ABSTRACT:"structural valve deterioration"` | 1, and it is a transapical TAVI survival/haemodynamics paper, not a durability model |
| `ABSTRACT:"predict" AND ABSTRACT:"structural valve deterioration" AND (ABSTRACT:"model" OR ABSTRACT:"algorithm")` | 3, all unrelated: a fluid-structure biomechanics simulation, a hypomagnesaemia association study, and a 20-year stentless valve series |

These are reproducible queries, not an impression. Anyone on the panel can re-run them.

Adjacent work exists, and none of it does what we propose:

- Prediction models after aortic valve replacement exist for **mortality**, heart-failure hospitalisation, technical failure, stroke and pneumonia. None for durability.
- Prediction models in this population have been kept current by **dynamic re-calibration over calendar time** for postoperative mortality (`D2`, 44,546 surgical patients), which shows static models drift here. Nobody has applied dynamic prediction from *serial echo measurements* to durability.
- Machine learning on echo exists for predicting **post-procedure transvalvular gradient waveform** from pre-procedure Doppler, and for choosing between surgical and transcatheter approach. Neither addresses long-term deterioration.
- The closest paper in spirit, Cartlidge and Dweck's *Detection and prediction of bioprosthetic aortic valve degeneration* (`B3`), predicts deterioration from a **PET-CT scan**, in 71 patients, with a single imaging biomarker. It is not a prediction model on routinely collected data, and PET is not a surveillance tool.
- Risk factors for SVD have been identified repeatedly with Cox models, but **never assembled into a validated, reported prediction model** with discrimination and calibration.
- Personalised imaging surveillance has been done in an adjacent field, endovascular aortic repair (`D10`), but not for prosthetic valves.
- The statistical machinery for scheduling the next measurement from a joint model already exists (`D1`, Rizopoulos) and has **never been applied to valve surveillance**.

### Where that leaves us — five defensible novelty claims

1. **First prediction model for bioprosthetic aortic valve durability.** Not a risk-factor analysis, a model, reported to TRIPOD+AI (`D8`).
2. **Label-noise-aware endpoint design.** Explicitly handle the Velders finding that 59–65% of SVD classifications do not persist: require confirmation, model the latent trajectory, and run all three definitions as a sensitivity analysis. To our knowledge nobody has done this.
3. **Surveillance interval as the model output, not a risk score.** Use the joint-model information-gain approach (`D1`) to answer "when is this patient's next echo due", replacing the guideline's fixed annual interval (`A2`) with a personalised one. This is the clinically actionable framing the rubric asks for.
4. **Valve-model-aware hierarchical structure.** Partial pooling across valve models so that a new prosthesis with little history borrows strength from its class rather than being unmodellable. This addresses the real deployment problem that new valves have no durability data, and it directly serves the challenge's stated interest in different valve types.
5. **First durability analysis on MIMIC-IV-ECHO**, a resource released three weeks ago, with free-text extraction of valve model and size from notes (`E1`, and the JAMIA 2025 paper reporting >96% accuracy for prosthetic valve detection) to recover the one variable the structured data lacks.

Claims 2 and 3 are the strongest, because they come from specific papers most teams will not have read, and both change what you actually build rather than just what you say.

---

## 4. What this implies for the protocol

- **Endpoint**: VARC-3 as primary, with persistence confirmation; Capodanno and Dvir as pre-specified sensitivity analyses. Baseline is the 3-month echo.
- **Model**: competing-risk survival with death as the competing event; landmark or joint model for dynamic updating with serial echoes; hierarchical valve-model effects; random survival forest and "time since implant alone" as baselines.
- **Metrics**: time-dependent area under the curve and C-index for discrimination, calibration plots at 5/8/10 years, and a decision-curve or surveillance-burden analysis showing echoes saved versus events missed.
- **Validation**: literature-calibrated synthetic cohort for the prototype, MIMIC-IV-ECHO as pre-registered external validation, registry data as the eventual target.
- **Honesty**: state that events are rare, that the label is noisy, and that no public dataset has valve model with serial echo. The panel includes clinicians who know this; naming the limitations first is worth more than hiding them.
