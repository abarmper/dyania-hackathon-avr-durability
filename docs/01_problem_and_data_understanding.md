# Problem and data understanding

Working notes for the Dyania Health Hackathon 2026 challenge (Sept 15–17). Purpose: understand what is being asked, what the clinical problem is, and what the three provided `.xlsx` files actually contain, so we can decide what to salvage and what to look for elsewhere.

---

## 1. What the challenge actually asks

**One-liner from the README:** "Build a study — using ML, a statistical model, or whatever approach you prefer — proposing a protocol to predict aortic valve durability in patients with a bioprosthetic aortic valve replacement."

Key points from reading the README, rubric, protocol/model/data templates and slide outline:

- **The deliverable is a study protocol, not a trained model.** `protocol/study_protocol.md` is the "main deliverable". A proof-of-concept notebook is optional ("evaluated positively if present").
- **It is deliberately vague.** No dataset is mandated. The data plan template explicitly invites "public datasets or registries (e.g. PARTNER trial public data extracts, STS/ACC TVT Registry public reports, synthetic echo follow-up datasets)". Section 5 of the data plan is literally "Synthetic or Proxy Data". So: any public dataset, synthetic data, or literature-derived simulation is fair game, as long as we are honest about it.
- **Scope is bioprosthetic AVR (both SAVR and TAVR).** Mechanical valves are out (they are listed as an exclusion example). The SAVR-vs-TAVR durability debate is mentioned repeatedly and is expected to be part of the framing. Differentiating failure trajectories by **valve type / manufacturer / model** is explicitly listed as a secondary objective and as a subgroup analysis, so taking the different valve generations into account is not just "a good idea", it is scored.
- **Scoring (rubric):** Clinical validity 25%, ML approach 25%, study design 20%, GitHub/docs 15%, presentation 15%. Half the score is "is it clinically grounded and is the ML appropriate". Study-design quality wants "publication-ready protocol": endpoints, cohort, sample size (events-per-variable), censoring, calibration, subgroups, comparator baseline.
- **Repeated hints of what they expect to see** (they appear in 3+ template files, so treat them as a checklist):
  - Time-to-event framing with right censoring (survival model / competing risks), or fixed-horizon classification at 5/8/10 years, or sequence model over serial echos.
  - Metrics: time-dependent AUROC, C-index, calibration at fixed horizons.
  - Ground truth per **VARC-3 / EAPCI-ESC-EACTS consensus** criteria for SVD, plus reintervention and explant pathology as a hierarchy.
  - Engineered features: mean gradient progression rate, EOA indexed to BSA, patient-prosthesis mismatch (PPM) flag, rate of change of peak velocity.
  - Informative missingness (a missed echo visit is itself a signal).
  - Output mapped to a clinical action: risk tier → recommended echo surveillance interval / structural heart referral.
  - Comparator: STS-PROM, published SVD nomograms, "time since implant alone", cardiologist gestalt.
  - SHAP / feature importance.
  - Fairness across demographics and across valve manufacturers.
- **Deadline:** PR open before the presentation session on Sept 17, branch `team/<name>`, do not merge. No real patient data in the repo.

---

## 2. The clinical problem, in plain terms

### 2.1 Background

- Aortic stenosis is the most common valve disease needing intervention. Treatment is aortic valve replacement (AVR), either **surgical (SAVR)** or **transcatheter (TAVR/TAVI)**.
- Valves are **mechanical** (last for life, but need lifelong warfarin) or **bioprosthetic** (bovine/porcine pericardium or porcine valve; no warfarin, but they wear out). TAVR valves are all bioprosthetic. Guidelines have shifted toward bioprosthetic valves in ever-younger patients, and TAVR is now used in low-risk 65–75 year olds, so **a lot of people will outlive their valve**.
- **Structural valve deterioration (SVD)** is the intrinsic wear-out of the bioprosthetic leaflets: calcification, fibrosis, pannus, leaflet tears. It causes prosthetic stenosis (rising gradient), prosthetic regurgitation (intraprosthetic leak), or both. It is distinct from **non-structural** dysfunction (paravalvular leak, patient-prosthesis mismatch, malposition), **valve thrombosis** (including subclinical leaflet thrombosis/HALT after TAVR) and **endocarditis**. The consensus umbrella term is **bioprosthetic valve dysfunction (BVD)**; when it produces symptoms or requires reintervention it becomes **bioprosthetic valve failure (BVF)**.
- Typical numbers to anchor the framing (verify in the literature step): surgical pericardial valves have roughly 10–15 year median durability; SVD incidence is very low in the first 5 years, then accelerates (roughly 5–10% at 10 years, 20–30%+ at 15 years, with wide variation by model and by patient age — younger patients calcify valves faster). Some models (Mitroflow, Trifecta, some Epic series) have notoriously worse durability. TAVR durability data beyond 5–8 years is still thin, which is the core of the "SAVR vs TAVR durability debate".

### 2.2 How SVD is defined (this is the ground truth the panel wants to see)

Two consensus documents are the reference:

**EAPCI/ESC/EACTS 2017 (Capodanno et al.)** — staged SVD by echo:
- *Morphological SVD*: leaflet thickening/calcification/tear without hemodynamic change.
- *Moderate hemodynamic SVD*: mean gradient ≥20 and <40 mmHg **and/or** ≥10 and <20 mmHg increase from baseline **and/or** new or worsening (≥1 grade) moderate intraprosthetic regurgitation.
- *Severe hemodynamic SVD*: mean gradient ≥40 mmHg **and/or** ≥20 mmHg increase from baseline **and/or** new or worsening (≥2 grades) severe intraprosthetic regurgitation.

**VARC-3 2021 (Généreux et al.)** — hemodynamic valve deterioration (HVD) and BVF:
- HVD Stage 1: morphological only. Stage 2 (moderate): mean gradient ≥20 mmHg with increase ≥10 mmHg from the early post-implant baseline (1–3 months), and/or new/worsening ≥moderate intraprosthetic AR. Stage 3 (severe): mean gradient ≥30 mmHg with increase ≥20 mmHg from baseline, and/or new/worsening severe AR.
- BVF Stage 1: any BVD with clinically expressive criteria (new symptoms, LV dilation/dysfunction, pulmonary hypertension) or irreversible Stage 3 HVD. Stage 2: aortic valve reintervention (redo SAVR or valve-in-valve TAVR). Stage 3: valve-related death.

The important modelling consequence: **the label needs a baseline echo (early post-implant) and at least one follow-up echo**, or a hard event (reintervention / explant / death). Change-from-baseline is part of the definition, so "mean gradient" alone is not a label.

### 2.3 Why this is a prediction/surveillance problem

- Current practice: echo at discharge/30 days, at 1 year, then every 1–3 years by fixed interval regardless of individual risk (guideline recommendations differ, adherence is patchy). Failing valves are often identified late, when the patient is symptomatic, with LV damage or in heart failure, and reintervention (redo surgery or valve-in-valve TAVR) carries higher risk when done urgently.
- Opportunity: **risk-stratify at implant and update at each follow-up** using routinely collected data (patient factors, valve model/size, PPM, early gradients, renal function, metabolic factors, anticoagulation), so that high-risk patients get closer surveillance and early elective reintervention planning, while low-risk patients get fewer echos.
- Known/suspected risk factors to look for in the literature step: younger age at implant, small valve size / PPM (indexed EOA ≤0.85 cm²/m²), specific valve models and generations (and anticalcification treatments), renal failure/dialysis, diabetes, metabolic syndrome, dyslipidemia/Lp(a), smoking, hyperparathyroidism and calcium-phosphate metabolism, female sex (mixed), absence of anticoagulation early after implant (thrombosis → SVD hypothesis), subclinical leaflet thrombosis (HALT) after TAVR, immune response (anti-Gal antibodies), valve-in-valve/small annulus, and elevated early post-implant gradients.

### 2.4 Modelling framings that fit the problem

- **Time-to-event with right censoring** (Cox, random survival forest, gradient-boosted survival, DeepSurv) predicting time to SVD (stage ≥2 HVD) or BVF, with **competing risk of death** (elderly cohort; Fine-Gray or cause-specific hazards).
- **Landmark / dynamic prediction** using serial echos (joint model of longitudinal gradient trajectory + survival; or landmark models at 1, 3, 5 years). This maps directly to "how do I update the surveillance interval at each visit".
- **Fixed-horizon classification** (SVD by 5/8/10 years) as a simpler comparator.
- Baselines the panel will accept: time since implant alone, STS-PROM, published SVD nomograms, valve-model-specific published SVD curves.

---

## 3. What is in the three data files

All three are single-sheet Excel exports ("Query result") from an Epic-style EHR (Caboodle keys, RxNorm, LOINC columns). Dates are **shifted and reduced to a year** (an integer like 2013; some rows say 2027, which is an artifact of date shifting). Free text is de-identified with tokens like `[NAME]`, `[DATE]`, `[AGE]`, `<PERSON>`. Two overlapping patient sets exist:

| Patients | Notes | Labs | Meds |
|---|---|---|---|
| `Patient_001`–`Patient_100` | yes (op report ± progress note) | no | no |
| `Patient_101`–`Patient_117` | yes (1 note each) | yes | yes |

So there is **no patient with the full trio of rich notes + labs + meds**. The 100 note-only patients are the clinically interesting ones; the 17 labs/meds patients are a peri-operative EHR sample.

### 3.1 `notes_deidentified.xlsx` — 215 notes, 117 patients (the valuable file)

| Field | Content |
|---|---|
| `Profile Key` | patient id |
| `Type` | Operative Report (112), Progress Notes (83), Procedures (15, mostly TAVR cath-lab notes and intra-op TEE reports), Discharge Summary (5) |
| `Service` | Cardiac Surgery 92, Cardiovascular Surgery 14, Cardiovascular Disease/Medicine 14, Thoracic Surgery 6, Anesthesiology 5, blank 84 |
| `Service Date` | year only, 2006–2026 (bulk 2010–2022) |
| `Notes` | free text, median ~4,000 chars, max ~25,000 |
| Signed Status / provider type / specialty | mostly Signed, Physician; 13 notes flagged Deleted |

Structure per patient: 100 patients have an index **operative report** (or TAVR procedure note); **71 of them also have a later progress note**. Follow-up span from index operation to last progress note:

| Follow-up (years) | Patients |
|---|---|
| 0–1 (peri-op / early) | 37 |
| 2–4 | 8 |
| 5–9 | 15 |
| 10–15 | 11 |

Index procedure mix (rough, from the procedure line of each op report; a few op reports are unreadable templates or non-valve surgery):

- **TAVR ~40 patients**: Edwards Sapien XT / S3 / S3 Ultra (sizes 20–29 mm), Medtronic CoreValve (26–31 mm), Evolut. Several from the PARTNER II era (2013–2015).
- **Surgical bioprosthetic AVR ~50 patients**: Carpentier-Edwards Perimount / Magna / Magna Ease (19–27 mm), St Jude Trifecta and Trifecta GT (19–29 mm), St Jude Biocor, St Jude Epic, Medtronic Freestyle, Inspiris, homograft, one Bio-Bentall (Konect). Many with concomitant CABG, myectomy, maze, aortoplasty, mini-sternotomy approach.
- **Mechanical or non-AVR surgery ~10 patients**: CABG-only, VSD closure, chest re-exploration, mechanical valves. These would be excluded in a study.

What the text contains that is extractable:
- **Index-procedure features**: valve manufacturer/model, labelled size, approach (full sternotomy / mini / transfemoral / axillary), native valve morphology (bicuspid vs tricuspid, ~17 patients bicuspid), annular enlargement, concomitant procedures, cross-clamp details, intra-op TEE result (gradient, "well seated, no paravalvular leak"), pre-op echo (AVA, indexed AVA, peak/mean gradient, DVI, EF), pre-op cath, PMH list (diabetes ~43, AF ~46, CKD/dialysis ~18, smoking ~56, statin ~70 patients mention it), BMI (~56 patients), weight (~77).
- **Follow-up echo values in progress notes**: 48 progress notes carry "mean gradient", 56 carry a dimensionless valve index (DVI), 38 carry EF, 68 mention a TTE; many quote the echo conclusion verbatim (e.g. "Trifecta prosthetic aortic valve (size #27). Signs of prosthetic valve stenosis ... peak gradient 46 mmHg, mean gradient 22 mmHg").
- **Outcome events (SVD / BVF) — around 15 patients.** From reading the contexts, the following look like true durability events (to be hand-adjudicated later):

| Patient | Index valve (year) | Event described |
|---|---|---|
| 017 | Trifecta 23 mm (2014) | severe prosthetic AS + 2+ AR from prosthetic dysfunction → redo sternotomy, Bio-Bentall (2026) |
| 030 | Trifecta 27 mm (2011) | prosthetic stenosis, mean gradient 22 mmHg (2024) |
| 032 | CE 23 mm (2013) | prosthetic AS + moderate AI, possible PPM (2024) |
| 035 | Biocor 27 mm (2010) | prosthetic stenosis → valve-in-valve TAVR; later 2+ AI incl. paravalvular |
| 036 | Biocor 25 mm (2009) | → ViV TAVR with 26 mm Evolut FX (2024) |
| 038 | Perimount 25 mm (2018) | progressive prosthetic stenosis → ViV TAVR 12/2023 (only ~5 years!), then query HALT |
| 042 | bioprosthetic AVR (2016) | moderate prosthetic AS, managed medically (2023) |
| 044 | Magna 25 mm (2009/10) | → ViV TAVR (2024) |
| 048 | CE 23 mm redo AVR (2010) | prosthetic valve dysfunction, evaluated for redo surgery (2024) |
| 058 | Trifecta GT 27 mm (2016) | prosthetic AS → ViV TAVR 08/2024 |
| 061 | bio-AVR (2015) | prosthetic AI → homograft root replacement (2024) |
| 068 | CE 21 mm (2012) | "prosthetic aortic valve failure" → ViV TAVR S3 23 mm |
| 071 | Freestyle 23 mm (2013) | "severe AS with structural valve deterioration" → redo (2018) |
| 089 | 21 mm bioprosthesis + CABG (2016) | severe prosthetic AI → planned redo AVR/CABG (2025) |
| 104 | (labs cohort) | "prosthetic aortic valve failure" in problem list |
| 103 | TAVR (2019) | endocarditis → redo AVR with Epic 23 mm (2023) — **non-SVD** failure, useful as a competing/excluded event |

Note the cohort is clearly **enriched for failures** (a random 100-patient AVR sample would not have ~15 reinterventions). It is a curated teaching sample, not an incidence estimate. Also note that valve-brand names are sometimes eaten by the de-identifier (`29mm [name] aortic valve`, `[address] jude trifecta`), so brand extraction needs fuzzy rules.

### 3.2 `labs_deidentified.xlsx` — 43,550 rows, 17 patients (`Patient_101`–`117`)

Long-format results: `Lab Component Name`, LOINC (only ~6% coded), `Result Date` (year), `Numeric Value` (62% numeric), `String Value`, `Flag`, `Unit`, reference range. **~40,700 rows have no `Lab Type`** and are actually respiratory-therapy flowsheet rows, ECG measurements and cardiac-device (pacemaker/CRT/ICD) interrogation values, not chemistry. Highlights:

- `Transcription` (8,677 rows): fragments of respiratory-therapy notes ("Work of breathing: ...", "Albuterol 2.5 mg unit dose"). Noise for our purpose.
- Peri-op ICU chemistry: arterial/venous blood gases, whole-blood Na/K/glucose/ionized Ca, lactate — hundreds of rows in the surgery year only.
- Standard chemistry/CBC/coag: creatinine + eGFR, BUN, electrolytes, calcium, magnesium, phosphorus (37 rows), LFTs, albumin, CBC, PT/INR, aPTT.
- Cardiac/metabolic markers: NT-proBNP (35), BNP (7), troponin T (47), HbA1c (27), lipid panel (~19 each), Lp(a) (1), CRP (1), ferritin (3).
- **`LV Ejection Fraction` (85 rows across 16 patients)** — the only echo-derived quantity in this file. No gradients, no AVA, no DVI.
- Device data: implant date, model, serial number, leads (Medtronic / St Jude / Boston Scientific), pacing thresholds, AT/AF burden — for pacemaker/CRT/ICD, **not** for the valve.
- Anatomic pathology (54 rows): specimen = "AORTIC VALVE" with clinical history ("severe AS; 1+ AI", "calcific AS") for ~8 patients, i.e. the explanted **native** valve, not an explanted prosthesis.
- Temporal coverage: 11 of 17 patients have labs in only 1–2 years (the surgical admission). Six patients (101, 102, 103, 105, 106, 117) have 6–23 distinct years, which is the only real longitudinal lab signal.

### 3.3 `medications_deidentified.xlsx` — 5,807 rows, same 17 patients

One row per order: proper name, generic, therapeutic class / pharmaceutical class / subclass, strength, form, route, start/end/discontinued year, `Mode` (Inpatient 5,003 / Outpatient 804), frequency, dose, discontinue reason (mostly "Auto DC at discharge"), RxNorm codes, prescriber specialty (26% filled; Thoracic Surgery, Cardiology, Anesthesiology), associated diagnoses (almost always empty; the few present include I35.0, Z95.2 "presence of prosthetic heart valve", T82.01XA "breakdown of heart valve prosthesis").

Mostly the **inpatient peri-operative MAR** (insulin sliding scale, metoprolol, furosemide, heparin, opioids, laxatives, PPIs, amiodarone). Chronic outpatient orders exist for the 6 longitudinal patients. SVD-relevant classes are extractable per patient: statins (all 17), warfarin (6), DOAC (3), antiplatelet (16), ACEi/ARB/ARNI, beta-blockers, loop diuretics, MRA, SGLT2 (3), insulin, calcium/vitamin D (all), one bisphosphonate. No phosphate binders.

---

## 4. What can be salvaged, and for what

Honest verdict: **none of the three files can train or validate an SVD model.** 17 patients with labs/meds, 71 patients with any follow-up, ~15 events, year-level dates, no structured echo table. But they are very useful as a *reference for what real-world data looks like* and for a few concrete prototype pieces:

1. **Feature schema and extraction demo (notes, 100 patients).** Build (regex + LLM) extraction of index-procedure fields from op reports: valve type, model, size, approach, native morphology, concomitant CABG, annular enlargement, intra-op TEE gradient/PVL, pre-op AVA/gradient/EF, PMH flags, BMI/weight. This directly answers the panel's "how realistic is access to these data" and "how do you handle free-text echo reports", and gives a small real table to show in the slides.
2. **Label-construction demo (progress notes, ~48 with gradients).** Parse follow-up echo conclusions into mean gradient, DVI, AR grade, EF; apply VARC-3 / EAPCI staging rules against the intra-op or pre-discharge baseline; show which patients would be Stage 2/3 HVD or BVF. The 15 event patients make a credible **case series / failure-mode illustration** (Trifecta and Biocor failures at 8–14 years, Perimount ViV at 5 years, 21–23 mm valves failing with PPM, endocarditis as non-SVD competitor).
3. **Valve-model coverage argument.** The sample already spans ~10 models and 3 TAVR generations. This supports a design decision to model valve generation as a feature (or stratify) and to require manufacturer-coded implant data in deployment.
4. **Ancillary EHR variables (labs/meds, 17 patients).** Use as the template for the "routinely collected" feature list: creatinine/eGFR trajectory, calcium/phosphorus, HbA1c, lipids, NT-proBNP, INR/anticoagulation exposure, statin exposure. Also a good example of *missingness*: only 6 of 17 patients have multi-year labs, phosphorus/PTH/Lp(a) are nearly absent, which argues for missingness indicators rather than heavy imputation.
5. **Data-quality lessons for the protocol**: year-only shifted dates (need month-level for gradient rate), brand names lost to de-identification, "Deleted" notes, flowsheet noise mislabelled as labs, no BSA (needed for indexed EOA; BMI/weight are in notes), device-registry data present for pacemakers but not for valves (the argument for a valve implant registry link).

What we still need from elsewhere (next step: literature and public-data search):
- A **cohort with serial post-implant echo values and time-to-SVD/reintervention**, or published Kaplan–Meier / hazard information by valve model, age and PPM that we can use to build a **literature-calibrated synthetic cohort** for the proof-of-concept (clearly labelled as such, with a plan to replace it with a registry).
- Candidate public sources to check: trial data-sharing platforms (YODA / Vivli / BioLINCC for PARTNER, CoreValve US Pivotal, NOTION, SURTAVI extracts), MIMIC-IV (notes + echo module; prosthetic AV patients identifiable, follow-up short), UK Biobank, FDA post-approval study summaries and MAUDE adverse-event reports by valve model, published SVD nomograms and long-term single-centre series (Perimount, Mitroflow, Trifecta, Epic, Hancock II, Sapien/Evolut 5–10 year data).

---

## 5. Quick reference: patient-level view of the 117 patients

- 001–100: index op/procedure note; 71 also have a later progress note; 26 have ≥5 years of follow-up; ~15 have a documented BVD/BVF event.
- 101–117: one note each (post-op progress note or discharge summary) + labs + meds; 6 of them have multi-year labs/meds; 104 has "prosthetic aortic valve failure" in the problem list; 103 had TAVR → endocarditis → redo SAVR.

Scripts used for this profiling were run ad hoc with pandas + openpyxl (system Python lacks both; a venv is needed). A reusable EDA notebook can be built from the same queries once we decide what to extract.
